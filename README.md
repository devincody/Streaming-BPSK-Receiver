# Streaming-BPSK-Receiver
Python &amp; Jupyter Notebook-based simulation of a streaming BPSK receiver.

![image](https://user-images.githubusercontent.com/5545285/144552554-49c4e1ac-bcc0-4dd6-9502-ff2b22e86f23.png)

**Figure 0:** A celebratory Elon Musk, obviously excited about the correct transmission, and reception of a simulated BPSK signal with a Eb/N0 of 4dB. This image required the transmission of 320Kb of data, which when simulated on my laptop, took ~5 real minutes meaning a link rate of ~1kb/sec. One cool thing from this image is that the first ~10 lines are completely jumbled, this is because the feedback loops haven't properly converged. The rest of the "speckles" are bit errors which occur at a rate of ~1% at a Eb/N0 ratio of 4dB.

The two main files in this repository are "BPSK Demodulator Single Run.ipynb" and "BPSK Receiver Performance.ipynb". The "Single Run" file uses the receiver to demoduled an image which has been encoded as a BPSK waveform (see Figure 0). The "Performance" notebook evaluates the bit-error rate of the receiver as a function of the Eb/N0 ratio (see Figure 2). Figure 1 below shows a rough outline of the receiver. We first apply a root-raised cosine filter to the data to reduce intersymbol interference. Next, we synchronize the sampling of the received data to the transmitting clock using a Mueller and Muller feedback loop. Lastly, we apply fine-frequency synchronization to the samples to remove any residual phase. 

![image](https://user-images.githubusercontent.com/5545285/144552830-7f345aca-4ce1-4856-867d-11888927b8aa.png)

The performance of the demodulator was better than I had hoped, but was only made possible by careful tuning of the feedback parameters by hand. If I had more time, I would certainly try to tune the parameters further in simulink or maybe PySim. As you can see in Figure 2, the performance of the receiver is very close to the theoretical predictions for the bit error rate of a coherently-demodulated BPSK receiver.

![image](https://user-images.githubusercontent.com/5545285/144552845-1488bc95-ff06-47a8-a818-de914f99d780.png)
