import numpy as np
import struct
import scipy
import matplotlib.pyplot as plt
import scipy.linalg
from adversary import FFT_HW_attack


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

    
    def binary(self,num):
        return ''.join('{:0>8b}'.format(c) for c in struct.pack('!f', num))
    
    def float_to_bin(self, num):
        return format(struct.unpack('!I', struct.pack('!f', num))[0], '032b')

    def HW_calc(self, x):
        x_re = x.real
        x_im = x.imag
        b_re = self.float_to_bin(x_re)
        b_im = self.float_to_bin(x_im)
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

        Eve = FFT_HW_attack.Eve(HW_seq, 0,0,0,0,0)
        Eve_key = Eve.DFT_SCA(N,0,N**2)

        return X
    
    def DFT(self, x):
        """Compute the discrete Fourier Transform of the 1D array x"""
        N = x.size
        n = np.arange(N)
        k = n.reshape((N,1))
        M = np.exp(-2j * np.pi * k * n / N)
        return np.dot(M, x)
    
    
    def FFT(self, x):
        #"""A recursive implementation of the 1D Cooley-Tukey FFT"""
        
        N = x.shape[0]
        
        #if N % 2 > 0:
         #   raise ValueError("size of x must be a power of 2")
        if N <= 16:  # this cutoff should be optimized
            return self.DFT(x)
        else:
            X_even = self.FFT(x[::2])
            X_odd = self.FFT(x[1::2])
            x_out = np.zeros(N,dtype=complex)
            w_n = np.exp(-2j * np.pi / N)
            w = 1
            for i in range(int(N/2)):
                x_out[i] = X_even[i] + w*X_odd[i]
                x_out[int(N/2) + i] = X_even[i] - w*X_odd[i]
                w = w_n*w
            #return np.concatenate([X_even + factor[:int(N / 2)] * X_odd,
            #                    X_even + factor[int(N / 2):] * X_odd])
            return x_out

    
    def reverse_bits(self, n, bitSize):
        result = 0
        for i in range(bitSize):
            if n & (1 << i):
                result |= 1 << (bitSize - 1 - i)
        return result
    
    def bit_reverse(self, x):
        N = x.shape[0]
        s = int(np.log2(N))
        for i in range(N):
            rb = self.reverse_bits(i,s)
            if(i < rb):
                tmp = x[i]
                x[i] = x[rb]
                x[rb] = tmp
        return x

    

    def DIT_FFT(self, x):
        N = x.shape[0]
        s = int(np.log2(N))
        x = np.asarray(self.bit_reverse(x),dtype=complex)
        for stage in range(1,s+1):
            p=0
            q=0 + 2**(stage-1)
            n=0
            while(n <= 2**(stage-1) and q <=N):
                w = np.exp(-2j*np.pi*n/(2**stage))
                y = x[p] + w*x[q]
                z = x[p] - w*x[q]
                x[p] = y
                x[q] = z
                p=p+1
                q=q+1
                n=n+1
                if(q%2**stage == 0):
                    p = p + 2**(stage-1)
                    q = q + 2**(stage-1)
                    n=0
        return x


    
    def FFT_hash_recursive(self):

        N = len(self.circulant_coeff)
        v = self.FFT(np.asarray(self.circulant_coeff,dtype=complex))#, record=False)
        y = self.FFT(np.asarray(self.padded_key,dtype=complex))#, record=False)
        u = np.multiply(v,y)
        cx_p = np.conj(self.FFT(np.conj(u)))/N
        cx = np.round(cx_p).real%2
        return cx[:self.hash_len], v, y, u, cx_p


    def FFT_hash_itterative(self):

        N = len(self.circulant_coeff)
        v = self.DIT_FFT(self.circulant_coeff)
        y = self.DIT_FFT(self.padded_key)
        u = np.multiply(v,y)
        cx_p = np.conj(self.FFT(np.conj(u)))/N
        cx = np.round(cx_p).real%2
        return cx[:self.hash_len], v, y, u, cx_p

