import numpy as np
import matplotlib.pyplot as plt

####Key and mask are both in binary format
N = 16
x = np.random.randint(0,2,N) ##key bit  values
m = np.random.randint(0,2,N) ##mask values
c = np.random.randint(0,2,N) ##this is binary, hence the key is binary. Can we byte format this?

def binary_format_test(x,m,c,N):
    ###compute correct, unmasked version for refernce
    v = np.fft.fft(c)
    y = np.fft.fft(x)
    print(y)
    u = np.multiply(v,y)
    K_pa = np.array(np.round(np.fft.ifft(u)).real%2,dtype=int)

    #print("y: {0}".format(y))
    #print("u: {0}".format(u))
    print("K_pa_x:    {0}".format(K_pa))

    ###compute mask only version
    y_m = np.fft.fft(m)
    u_m = np.multiply(v,y_m)
    K_pa_m = np.array(np.round(np.fft.ifft(u_m)).real%2,dtype=int)
    print("K_pa_m:    {0}".format(K_pa_m))

    ###compute masked key version
    xm = x^m
    y_xm = np.fft.fft(xm)
    u_xm = np.multiply(v,y_xm)
    K_pa_xm = np.array(np.round(np.fft.ifft(u_xm)).real%2,dtype=int)
    print("K_pa_xm:   {0}".format(K_pa_xm))

    K_pa_xm_m = K_pa_xm ^ K_pa_m
    print("K_pa_xm_m: {0}".format(K_pa_xm_m))

    if (K_pa_xm_m == K_pa).all():
        return True
    else:
        return False



result = binary_format_test(x,m,c,N)
print("Masking correct: {0}".format(result))

####Key and mask are both in byte format

def BitsToBytes(y):
    """
    Converts a bit string into a byte string using little endian order
    
    Input: A but string y of length alpha

    Output: A byte string z of length ceil(alpha/8)
    """
    alpha = len(y)
    z = np.zeros(int(np.ceil(alpha/8)), dtype=int)
    for i in range(alpha): ##alpha not alpha-1
        z[int(np.floor(i/8))] = z[int(np.floor(i/8))] + y[i]*2**(i%8)
    return z    


def BytesToBits(z):
    """
    Converts a byte string into a bit string using little endian order

    Input: A byte string z of length alpha

    Output: A bit string y of length alpha*8
    """
    alpha = len(z)
    y = np.zeros(8*alpha, dtype=int)
    z_dash = z.copy()
    for i in range(alpha): ##alpha not alpha-1
        for j in range(8):
            y[(8*i) + j] = z_dash[i]%2
            z_dash[i] = np.floor(z_dash[i]/2)
    return y

x_bytes = BitsToBytes(list(x))
m_bytes = BitsToBytes(list(m))
c_bytes = BitsToBytes(list(c))
N_bytes = np.zeros(int(np.ceil(len(list(x))/8)), dtype=int)

def bytes_format_test(x,m,c,N):
    v = np.fft.fft(c)
    y = np.fft.fft(x)
    u = np.multiply(v,y)
    K_pa = np.fft.ifft(u).real
    return K_pa

bytes_res = bytes_format_test(x_bytes, m_bytes, c_bytes, N_bytes)
in_bits = BytesToBits(bytes_res)

print("Bytes version: {0}".format(bytes_res))



def y_distribution_experiment(repeats,N):

    real_vals = []
    imag_vals = []

    for _ in range(repeats):
        x = np.random.randint(0,2,N)
        y = np.fft.fft(x)

        real_vals.extend(y.real)
        imag_vals.extend(y.imag)

    num_distinct_real = len(set(real_vals))
    print("Number of Distinct Real values: {0}".format(num_distinct_real))

    num_distinct_imag = len(set(imag_vals))
    print("Number of Distinct Imag values: {0}".format(num_distinct_imag))

    plt.figure(figsize=(8,5))
    plt.hist(real_vals, bins = 500, density=True, alpha=0.4, label="Real")
    plt.hist(imag_vals, bins = 500, density=True, alpha=0.4, label="Imaginary")
    plt.xlabel("FFT Magnitude")
    plt.ylabel("Probability Density")
    plt.title(f"Distribution of FFT Magnitudes (N={N}, trials={repeats})")
    plt.grid(True)
    plt.legend()
    plt.yscale("log")
    plt.tight_layout()
    plt.show()

import struct

##a is a string of binary float rep
def float_bin_array(a):
    N = len(a)
    a_arr = np.zeros(N, dtype=int)
    for i in range(N):
        a_arr[i] = a[i]
    return a_arr


def convert_bin_array_to_float(arr):
    bit_string = ''.join(str(b) for b in arr)

    # Convert bit string to integer
    int_val = int(bit_string, 2)

    # Pack as 32-bit unsigned int, then unpack as float
    float_val = struct.unpack('!f', int_val.to_bytes(4, byteorder='big'))[0]
    return float_val


def fp_masked_multiply(a,b):

    c = a*b
    a_bin = np.binary_repr(np.float32(a).view(np.int32),width=32)
    b_bin = np.binary_repr(np.float32(b).view(np.int32),width=32)


    a_bin_arr = float_bin_array(a_bin)
    b_bin_arr = float_bin_array(b_bin)


    mask_a = np.random.random()
    mask_a_bin = float_bin_array(np.binary_repr(np.float32(mask_a).view(np.int32),width=32))
    mask_b = np.random.random()
    mask_b_bin = float_bin_array(np.binary_repr(np.float32(mask_b).view(np.int32),width=32))

    a_m = np.bitwise_xor(a_bin_arr, mask_a_bin)
    b_m = np.bitwise_xor(b_bin_arr, mask_b_bin)

    a_m_float = convert_bin_array_to_float(a_m)
    b_m_float = convert_bin_array_to_float(b_m)

    c_m = a_m_float * b_m_float
    mask_c = mask_a * mask_b

    c_m_bin = float_bin_array(np.binary_repr(np.float32(c_m).view(np.int32),width=32))
    mask_c_bin = float_bin_array(np.binary_repr(np.float32(mask_c).view(np.int32),width=32))



    c_unmasked = np.bitwise_xor(c_m_bin, mask_c_bin)

    c_unmasked_float = convert_bin_array_to_float(c_unmasked)

    return [c_unmasked == c]


a = np.random.random()
b = np.random.random()
print(fp_masked_multiply(a,b))
            
#y_distribution_experiment(20000,1024)

