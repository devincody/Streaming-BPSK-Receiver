from dataclasses import dataclass
from collections import deque
import numpy as np

class FIR():
    def __init__(self, taps):
        self.taps = taps
        self.N = len(self.taps)
        self.shift_reg = np.zeros(self.N, complex)
        
        # Instrumentation
        self.data_out = deque()
        
    def update(self, sample):
        for i in reversed(range(self.N-1)):
            self.shift_reg[i+1] = self.shift_reg[i]
        self.shift_reg[0] = sample
        out = np.dot(self.shift_reg, self.taps)
        self.data_out.append(out)
        return out
    
class decimator():
    def __init__(self, M, taps):
        self.M = M
        
        # Zero pad taps
        self.num_taps = int(np.ceil(len(taps)/self.M)*self.M)
        self.taps = np.zeros(self.num_taps)
        self.taps[:len(taps)] = taps
        
        # Create filter structure
        self.filters = []
        for i in range(self.M):
            self.filters.append(FIR(self.taps[i::self.M]))
            
        # Memory for rolling sum
        self.nsamples = 0
        self.sum_out = 0
            
        # Instrumentation
        self.data_out = deque()
    
    def update(self, sample):
        self.sum_out += self.filters[self.nsamples].update(sample)
        if self.nsamples == self.M - 1:
            tmp = self.sum_out
            self.sum_out = 0
            self.nsamples = 0
            self.data_out.append(tmp)
            return tmp
        else:
            self.nsamples += 1
            
class interpolator:
    def __init__(self, L, taps):
        self.L = L
        
        # Zero pad taps
        self.num_taps = int(np.ceil(len(taps)/self.L)*self.L)
        self.taps = np.zeros(self.num_taps)
        self.taps[:len(self.taps)] = taps

        # Create filter structure
        self.filters = []
        for i in range(self.L):
            self.filters.append(FIR(self.taps[i::self.L]))
            
        # Instrumentation
        self.data_out = deque()      
        
    def update(self, sample):
        out = [filt.update(sample) for filt in self.filters]
        self.data_out.append(out)
        return out

@dataclass
class mandm:
    """ Mueller and Muller method for clock recovery"""
    mu_gain: float = 0.008
    iL: int = 8 # input oversampling rate
    L: int = 16 # Intermediary Interpolation Rate
    Ltaps: np.array = np.ones(1)
    
    def __post_init__(self):
        self.past_samples = np.zeros(3, complex)
        self.interp = interpolator(self.L, self.Ltaps)
        self.input_count = 0
        self.index = 0
        self.in_rate_mu = 0
        self.mu: int = 0
        self.res_mu = 0
        
        # Instrumentation
        self.log_out = []
        self.log_mu = []
        self.log_idx = []
        self.log_in_rate_mu = []
        self.log_interp_rate_mu = []
        self.log_res_mu = []
        
    def update(self, sample):
        fine_samples = self.interp.update(sample)
        
        if self.input_count == self.in_rate_mu: # coarse synchronization
            # Update shift register
            for i in range(2):
                self.past_samples[i] = self.past_samples[i+1]
            self.past_samples[2] = fine_samples[self.mu] # fine synchronization
            
            # Timing Error Detection
            out_rail = [np.sign(x.real) + 1j*np.sign(x.imag) for x in self.past_samples]
            x = (out_rail[2] - out_rail[0]) * np.conj(self.past_samples[1])
            y = (self.past_samples[2] - self.past_samples[0]) * np.conj(out_rail[1])
            mm_val = np.real(y - x)
            
            self.res_mu += self.iL*self.L + self.mu_gain*mm_val # offset to next sample
            self.index += np.floor(self.res_mu) # index of next sample
            
            self.log_mu.append(self.res_mu)
            
            self.res_mu -= (self.L - self.mu) # account for samples aready produced by interp filter
            self.in_rate_mu = self.res_mu//self.L
            self.mu = int(self.res_mu % self.L)
            
            #running error
            self.res_mu -= self.in_rate_mu*self.L + self.mu
            self.input_count = 0
            
            # Logging
            self.log_out.append(self.past_samples[2])
            self.log_idx.append(self.index)
            self.log_in_rate_mu.append(self.in_rate_mu)
            self.log_interp_rate_mu.append(self.mu)
            self.log_res_mu.append(self.res_mu)
            
            return self.past_samples[2]
        else:
            self.input_count += 1
    
@dataclass
class costas_loop:
    alpha: float = 0.132
    beta: float = 0.00932
    vco_phase: float = 0
    vco_freq: float = 0
        
    def __post_init__(self):
        self.log_out = deque()
        self.log_vco_out = deque()
        self.log_error = deque()
        self.log_vco_freq = deque()
        self.log_vco_phase = deque()

    def update(self, sample: complex) -> complex:
        vco_out = np.exp(-1j*self.vco_phase)
        out = sample*vco_out
        error = out.real * out.imag
        
        self.vco_freq += self.beta * error
        self.vco_phase += self.alpha * error + self.vco_freq
        
        # Record results
        self.log_out.append(out)
        self.log_vco_out.append(vco_out)
        self.log_error.append(error)
        self.log_vco_freq.append(self.vco_freq)
        self.log_vco_phase.append(self.vco_phase)
        
        return out # demodulated data