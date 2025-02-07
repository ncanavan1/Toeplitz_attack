import numpy as np
from privacy_amplification import toeplitz_hashing as th
from adversary import power_analysis
import matplotlib.pyplot as plt
from graphing import plotting

def gen_ec_key(N):
    return np.random.randint(0,2,N,dtype=int)


def main():
    N = 100
    ec_key = gen_ec_key(N) ##should remain unknown to Eve
    pa_len = 70
    l_key = N
    sample_number = 10
    sigma_range = np.arange(0,2,0.2)
    sigma_acc = []


    TM_test = th.Toeplitz_hashing(pa_len,ec_key,l_key)
    row_spec = TM_test.TM_coeff_row
    col_spec = TM_test.TM_coeff_col
    TM_test.calc_pa_explicit()
    pa_test = TM_test.pa_key

    print("ROW SPEC: {0}".format(row_spec))
    print("COL SPEC: {0}".format(col_spec))
    print("EC KEY: {0}".format(ec_key))
    print("PA KEY: {0}".format(pa_test))


    for sigma in sigma_range:
        ##Alice's genuine computation
        TM_target = th.Toeplitz_hashing(pa_len,ec_key,l_key,sigma=sigma)
        TM_target.calc_pa_explicit()
        target_trace = TM_target.simulated_power

        ##Alice send's Bob toeplitz matrix specification, Eve listens
        TM_coeff_row_pub = TM_target.TM_coeff_row
        TM_coeff_col_pub = TM_target.TM_coeff_col

        ##Eve attacks
        Eve = power_analysis.Eve(target_trace,pa_len,l_key,TM_coeff_row_pub, 
                                TM_coeff_col_pub, sample_number, sigma)
        Eve.run_attack()
        ec_guess = Eve.ec_key_guess

        accuracy = 0
        for i in range(N):
            if ec_key[i] == ec_guess[i]:
                accuracy += 1
        accuracy = accuracy/N
        sigma_acc.append(accuracy)
    print("Key Recovery Accuracy for sigma = {0}: {1}".format(sigma,accuracy))

    plt.figure()
    plt.plot(sigma_range,sigma_acc)
    plt.show()






if __name__=="__main__":
    main()