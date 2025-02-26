import numpy as np
from scipy import stats
from privacy_amplification import fft_toeplitz as th
from privacy_amplification.FFT.FFT import *
import struct

class Eve:
    def __init__(self, target_trace, hash_len, key_len, Toeplitz_spec, sigma):
        self.target_trace = target_trace
        self.hash_len = hash_len
        self.key_len = key_len
        self.Toeplitz_spec = Toeplitz_spec 
        self.sigma = sigma


    def binary(num):
        return ''.join('{:0>8b}'.format(c) for c in struct.pack('!f', num))

    ##compute HW model for Butterfly FFT
    def HW_Model(self, X):
        return DIT_FFT_leaky(X)
        
    def guess_string(self, num, hyp_size):
        b_string = binary(num)
        guess = np.zeros(hyp_size)
        for i in range(len(b_string)):
            guess[i] = int(b_string[i])
        return guess

    def run_attack(self, hyp_size):
        best_guesses = []
        for i in range(2**hyp_size):
            guess = self.guess_string(i, hyp_size)
            x, HW = DIT_FFT_leaky(guess)
            ##calculate correlation
            
            ##if correlation is max,,




















    #def isolate_FFT():

    #def isolate_DFT():

    def DFT_SCA(self, N, trace_start, trace_end):
        DFT_trace = self.Traces[trace_start: trace_end]
        X = np.zeros(N)
         
        # Step 1: For each of the N secret bits, the FFT operates over 
        # these bits N times. That is, there is N samples corresponding to
        # the ith bit. 

        # The first N samples correspond to bit 0. the Nth to 2Nth to bit 1...

        # For the first bit, if HW is non-zero, then bit is 1. If 0, then bit is
        # 0 also.

        # For bit 1, if the HW trace is the same across the N samples, then this
        # implies the bit value is 0. If it changes, this implies 1.

        # Repeat this check on all neighbouring bit ranges.

        samples_left = DFT_trace[0:N]
        if (samples_left != 0).any():
            X[0] = 1
            #Default is 0
        
        for check in range(1,N):
            samples_left = DFT_trace[(check-1)*N:check*N]
            samples_right = DFT_trace[check*N:(check+1)*(N)]

            ##Most basic version, will need CPA in realisitc model
            if (samples_left == samples_right).all():
                X[check] = 0
            else:
                X[check] = 1

        return X