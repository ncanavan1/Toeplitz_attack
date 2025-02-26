import numpy as np
import struct
import scipy
import matplotlib.pyplot as plt
import scipy.linalg
from privacy_amplification.FFT.FFT import FFT, DIT_FFT, DIT_FFT_leaky


class FFT_Toeplitz_hashing:
    def __init__(self, hash_len, key_len, ec_key, TM_first_row_coeff=[], TM_first_col_coeff=[], sigma=0.1):

        self.hash_len = hash_len ##output hash length
        self.key_len = key_len ##input EC key length
        self.hashed_key = np.zeros(self.hash_len)
        self.ec_key = ec_key
        
        if len(TM_first_row_coeff) == 0 and len(TM_first_col_coeff) == 0:
            self.TM_first_row_coeff = np.random.randint(0,2,self.key_len)
            self.TM_first_col_coeff = np.random.randint(0,2,self.hash_len)
            self.TM_first_col_coeff[0] = self.TM_first_row_coeff[0]

        else:
            self.TM_first_row_coeff = TM_first_row_coeff
            self.TM_first_col_coeff = TM_first_col_coeff
        self.sigma = sigma
        self.simulated_power = []
        self.circulant_coeff = np.zeros(self.hash_len + self.key_len - 1)
        self.padded_key = np.zeros(self.hash_len + self.key_len - 1)

    def embed_toeplitz_on_circulant(self):
        ## This is how scipy does it https://github.com/scipy/scipy/blob/v1.15.1/scipy/linalg/_basic.py#L1960
        self.circulant_coeff = np.concatenate((self.TM_first_col_coeff, self.TM_first_row_coeff[-1:0:-1]))
        #self.circulant_coeff = np.array([0,1,1,0,1,0,1])

        self.padded_key[:self.key_len] = self.ec_key



    def print_matrices(self):

        print("\nFirst Row Coefficents: \n")
        print(self.TM_first_row_coeff)

        print("\nFirst Col Coefficients: \n")
        print(self.TM_first_col_coeff)

        print("\nToeplitz Matrix: \n")

        TM = np.zeros([self.hash_len, self.key_len])
        for j in range(self.hash_len):
            for i in range(self.key_len):
                if i >= j:
                    TM[j,i] = self.TM_first_row_coeff[i-j] 
                else:
                    TM[j,i] = self.TM_first_col_coeff[j-i]
        print(TM)

        print("\nCirculant Matrix:\n")
        CM = np.zeros([self.hash_len + self.key_len -1, self.hash_len + self.key_len -1])
        for j in range(self.hash_len + self.key_len - 1):
            CM[j,:] = np.roll(self.circulant_coeff[:],j)

        print(CM)

    
    def FFT_hash_recursive(self):

        N = len(self.circulant_coeff)
        v = FFT(np.asarray(self.circulant_coeff,dtype=complex))#, record=False)
        y = FFT(np.asarray(self.padded_key,dtype=complex))#, record=False)
        u = np.multiply(v,y)
        cx_p = np.conj(FFT(np.conj(u)))/N
        cx = np.round(cx_p).real%2
        return cx[:self.hash_len], v, y, u, cx_p


    def FFT_hash_itterative(self):

        N = len(self.circulant_coeff)
        v, HW = DIT_FFT_leaky(self.circulant_coeff)
        y = DIT_FFT(self.padded_key)
        u = np.multiply(v,y)
        cx_p = np.conj(DIT_FFT(np.conj(u)))/N
        cx = np.round(cx_p).real%2
        return cx[:self.hash_len], v, y, u, cx_p, HW

