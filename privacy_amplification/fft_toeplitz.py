import numpy as np
import struct
import scipy
import matplotlib.pyplot as plt
import scipy.linalg


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

        
        
    def DFT(self, x):
        N = len(x)
        n = np.arange(N)
        k = n.reshape((N,1))
        e = np.exp(-2j * np.pi * k * n / N)
        X = np.dot(e,x)
        return X
    
    def binary(self,num):
        return ''.join('{:0>8b}'.format(c) for c in struct.pack('!f', num))
    
    def HW_calc(self, x):
        x_re = x.real
        x_im = x.imag
        b_re = self.binary(x_re)
        b_im = self.binary(x_im)
        HW = 0
        for b in b_re:
            if b == "1":
                HW += 1
        for b in b_im:
            if b == "1":
                HW += 1
        return HW



    def DFT_explicit(self, x):
        N = len(x)
        n = np.arange(N)
        k = n.reshape((N,1))
        X = np.zeros(N,dtype=complex)
        HW_seq = []
        bit_val_seq = []

        ##for each k
        for n_f in n:
            for k_f in k:
                e = np.exp(-2j * np.pi * k_f * n_f / N)
                X[k_f] = X[k_f] + x[n_f] * e
                HW_seq.append(self.HW_calc(X[k_f]))
                if x[n_f] == 0:
                    bit_val_seq.append(0)
                else:
                    bit_val_seq.append(1)
        bit_val_seq = np.asarray(bit_val_seq) * max(HW_seq)
        HW_seq = np.asarray(HW_seq)

        HW_diff = []
        HW_diff.append(HW_seq[0])
        for i in range(len(HW_seq)-1):
            HW_diff.append(np.abs(HW_seq[i]-HW_seq[i+1]))


        plt.plot(HW_seq, label="HW")
        plt.plot(HW_diff, label="HW Diff")
        plt.plot(bit_val_seq, label="Bit Value")
        plt.xlabel("Trace Sample Point")
        plt.ylabel("Hamming Weight")
        plt.legend()
        plt.tight_layout()
        plt.show()
        return X
    
    def IDFT(self, X):
        N = X.size
        n = np.arange(N)
        k = n.reshape((N,1))
        e = np.exp(2j * np.pi * k * n / N)
        x = np.dot(e, X)/N
        return x
    
    def DFT_hash(self):
        v = self.DFT(self.circulant_coeff)
        y = self.DFT_explicit(self.padded_key)
        u = np.multiply(v,y)
        cx = self.IDFT(u)
        cx = np.round(cx).real%2
        return cx[:self.hash_len]
    
    ##Cooly-Turkey 1-D Method
    def FFT(self, X):
        N = X.size

        if N % 2 > 0:
            p = 0
            while(2**p < N):
                p = p+1
            N_new = 2**p
            X = np.concatenate((X,np.zeros(N_new-N)))
            N = N_new

        if N<=16:
            return self.DFT(X)
        else:
            X_even = self.FFT(X[::2])
            X_odd = self.FFT(X[1::2])
            factor = np.exp(-2j * np.pi * np.arange(N) / N)
            return np.concatenate([X_even + factor[:int(N / 2)] * X_odd, 
                                   X_even + factor[int(N / 2):] * X_odd])

    ###Cooley-Turky 1-D method
    def IFFT(self, X):
        N = X.size

        if N % 2 > 0:
            p = 0
            while(2**p < N):
                p = p+1
            N_new = 2**p
            X = np.concatenate((X,np.zeros(N_new-N)))
            N = N_new
            
        if N<=16:
            return self.IDFT(X)
        else:
            X_even = self.IFFT(X[::2])
            X_odd = self.IFFT(X[1::2])
            factor = np.exp(-2j * np.pi * np.arange(N) / N)
            return np.concatenate([X_even + factor[:int(N / 2)] * X_odd, 
                                   X_even + factor[int(N / 2):] * X_odd])
    
    def FFT_hash(self):
        v = self.FFT(self.circulant_coeff)
        y = self.FFT(self.padded_key)
        u = np.multiply(v,y)
        cx = self.IFFT(u)
        cx = np.round(cx).real%2
        return cx[:self.hash_len]