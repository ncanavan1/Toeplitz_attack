import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from privacy_amplification import fft_toeplitz as th
from privacy_amplification import toeplitz_hashing as th_basic


def gen_ec_key(N):
    return np.random.randint(0,2,N,dtype=int)



##the fft solutions and basic matrix multiplications should be equal
def verify_correct_outputs(key_len, hash_len):
    ec_key = gen_ec_key(key_len)


    TH_Bas = th_basic.Toeplitz_hashing(hash_len, ec_key, key_len)
    TH_Bas.calc_pa_key()
    basic_hash = TH_Bas.pa_key
    print("\n{0}\n".format(basic_hash))

        
    TH_FFT = th.FFT_Toeplitz_hashing(hash_len, key_len, ec_key,
                                     TM_first_row_coeff=TH_Bas.TM_coeff_row, TM_first_col_coeff=TH_Bas.TM_coeff_col)
    TH_FFT.embed_toeplitz_on_circulant()
    fft_hash = TH_FFT.fft_hash()

    print("\n{0}\n".format(fft_hash))





def main():
    key_len = 10
    hash_len = 8
    ec_key = gen_ec_key(key_len)

    TH = th.FFT_Toeplitz_hashing(hash_len, key_len, ec_key)
    TH.embed_toeplitz_on_circulant()
    TH.print_matrices()

    hashed_key = TH.fft_hash()
    
    print("\nHashed Key FFT:\n")
    print(hashed_key)

    TH_B = th_basic.Toeplitz_hashing(hash_len, ec_key, key_len, 
                                     TM_coeff_row=TH.TM_first_row_coeff, TM_coeff_col=TH.TM_first_col_coeff)
    
    TH_B.calc_pa_explicit()
    basic_key = TH_B.pa_key

    print("\nBasic Key:\n")
    print(basic_key)


    short_fft = TH.DFT(ec_key)
    explicit_fft = TH.DFT_explicit(ec_key)

    print(short_fft)
    print("\n\n")
    print(explicit_fft)



if __name__=="__main__":
    #main()
    verify_correct_outputs(10,8)