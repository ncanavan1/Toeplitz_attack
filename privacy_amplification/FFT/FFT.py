import numpy as np
import matplotlib.pyplot as plt
import struct


    
def binary(num):
    return ''.join('{:0>8b}'.format(c) for c in struct.pack('!f', num))

def float_to_bin(num):
    return format(struct.unpack('!I', struct.pack('!f', num))[0], '032b')

def HW_calc(x):
    x_re = x.real
    x_im = x.imag
    b_re = float_to_bin(x_re)
    b_im = float_to_bin(x_im)
    HW = 0
    for b in b_re:
        if b == "1":
            HW += 1
    for b in b_im:
        if b == "1":
            HW += 1
    return HW



def DFT_explicit(x):
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
            HW_seq.append(HW_calc(X[k_f]))
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

   # Eve = FFT_HW_attack.Eve(HW_seq, 0,0,0,0,0)
   # Eve_key = Eve.DFT_SCA(N,0,N**2)

    return X

def DFT(x):
    """Compute the discrete Fourier Transform of the 1D array x"""
    N = x.size
    n = np.arange(N)
    k = n.reshape((N,1))
    M = np.exp(-2j * np.pi * k * n / N)
    return np.dot(M, x)


def FFT(x):
    #"""A recursive implementation of the 1D Cooley-Tukey FFT"""
    
    N = x.shape[0]
    
    #if N % 2 > 0:
        #   raise ValueError("size of x must be a power of 2")
    if N <= 16:  # this cutoff should be optimized
        return DFT(x)
    else:
        X_even = FFT(x[::2])
        X_odd = FFT(x[1::2])
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


def reverse_bits(n, bitSize):
    result = 0
    for i in range(bitSize):
        if n & (1 << i):
            result |= 1 << (bitSize - 1 - i)
    return result

def bit_reverse(x):
    N = x.shape[0]
    s = int(np.log2(N))
    for i in range(N):
        rb = reverse_bits(i,s)
        if(i < rb):
            tmp = x[i]
            x[i] = x[rb]
            x[rb] = tmp
    return x



def DIT_FFT(x):
    N = x.shape[0]
    s = int(np.log2(N))
    x = np.asarray(bit_reverse(x),dtype=complex)
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


def DIT_FFT_leaky(x):
    N = x.shape[0]
    s = int(np.log2(N))
    x = np.asarray(bit_reverse(x),dtype=complex)
    HW_trace = []
    for stage in range(1,s+1):
        p=0
        q=0 + 2**(stage-1)
        n=0
        while(n <= 2**(stage-1) and q <=N):
            w = np.exp(-2j*np.pi*n/(2**stage))
            y = x[p]
            z = x[q]
            
            HW_trace.append(HW_calc(x[p]) + HW_calc())

            z = z*w
            ###HW before and after this assigment
            x[p] = y + z
            x[q] = y - z
            #####################################
            HW_trace.append(HW_calc(x[p]) + HW_calc(x[q]))


            p=p+1
            q=q+1
            n=n+1
            if(q%2**stage == 0):
                p = p + 2**(stage-1)
                q = q + 2**(stage-1)
                n=0
    return x, HW_trace


if __name__=="__main__":
    arr = np.zeros(32)
    for i in range(32):
        arr[i] = i
    arr_r = bit_reverse(arr)
    k=7
    x, HW = DIT_FFT_leaky(np.asarray([0,1,0,1,0,1,0,0]))
    #DIT_FFT(np.asarray([1,0,1,1,0,1,0,0,1,0,0,1,0,1,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0]))