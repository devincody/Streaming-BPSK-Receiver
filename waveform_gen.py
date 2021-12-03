import numpy as np
import commpy

def gen_waveform(bits_enc, L, f_offset, EbN0dB):
    N = len(bits_enc)
    EbN0lin = 10**(EbN0dB/10)
    N0 = 1/EbN0lin*L

    bits_pulse_train = np.pad(bits_enc.reshape((N, 1)), (0, L-1)).ravel()

    # Apply Root Raised Cosine Filter
    alpha = 0.35
    num_taps_rrcos = 101
    _, rrcos_taps = commpy.filters.rrcosfilter(num_taps_rrcos, alpha = 0.35, Ts = 1, Fs = L)
    bits_filt = np.convolve(bits_pulse_train, rrcos_taps)

    # Impairments
    bits_clk_offset = bits_filt * np.exp(2j*np.pi*np.arange(len(bits_filt))*f_offset)
    noise = np.sqrt(N0/2)*(np.random.randn(len(bits_clk_offset)) + 1j*np.random.randn(len(bits_clk_offset)))
    return bits_clk_offset + noise