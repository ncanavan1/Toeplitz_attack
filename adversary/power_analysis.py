import numpy as np
import matplotlib.pyplot as plt
from privacy_amplification import toeplitz_hashing as th
from graphing import plotting
from scipy.stats import pearsonr

class Eve:
    def __init__(self, target_trace, pa_len, l_key, TM_coeff_row, 
                 TM_coeff_col, sample_number, sigma): ##for now work with the case N=l_key
        self.target_trace = target_trace
        self.pa_len = pa_len
        self.l_key = l_key
        self.TM_coeff_row = TM_coeff_row
        self.TM_coeff_col = TM_coeff_col
        self.reference_traces = []
        self.sample_number = sample_number
        self.sigma = sigma
        self.ec_key_guess = np.zeros(l_key)

    def generate_reference_traces(self, bit_pos):
        ##We only need to observe traces where all possible key bits are 0
        ##excepts for target bit, which we try 0 and 1. Repeat for N times

        for bit_value in range(0,2): ##don't really need to generate traces for 0 hypothesis as will average out to 0
            ec_key_hypothesis = np.zeros(self.l_key)
            ec_key_hypothesis[bit_pos] = bit_value
            for sample in range(self.sample_number):
                TM = th.Toeplitz_hashing(self.pa_len, ec_key_hypothesis, self.l_key, sigma=self.sigma, 
                                            TM_coeff_row=self.TM_coeff_row, TM_coeff_col=self.TM_coeff_col)
                TM.calc_pa_explicit()
                self.reference_traces.append([bit_pos, bit_value, TM.simulated_power])

    ##finds the sample points that are related to bit_pos
    def find_sample_points(self, bit_pos):
        sample_points = []
        for i in range(self.pa_len): ##number of matrix entries
            sample_points.append(bit_pos + i*self.l_key)
        return sample_points
    
    def find_relevant_power(self, bit_pos, trace):
        power_sample = []
        for i in range(self.pa_len):
            power_sample.append(trace[bit_pos + i*self.l_key])
        return power_sample
    
    def determine_ec_key_bit_value_pearson(self, hyp0_avg, hyp1_avg, target_relevant_power):
        corr0, _ = pearsonr(hyp0_avg, target_relevant_power)
        corr1, _ = pearsonr(hyp1_avg, target_relevant_power)
        if corr0 > corr1:
            return 0, corr0
        else:
            return 1, corr1

    def determine_ec_key_bit_value_MSE(self, hyp0_avg, hyp1_avg, target_relevant_power):
        corr0 = np.square(np.subtract(hyp0_avg, target_relevant_power)).mean()
        corr1 = np.square(np.subtract(hyp1_avg, target_relevant_power)).mean()
        if corr0 < corr1:
            return 0, corr0
        else:
            return 1, corr1


    def run_attack(self):

        for bit_pos in range(self.l_key): ##uncomment later
            self.generate_reference_traces(bit_pos)
            hyp0_avg = []
            hyp1_avg = []
            for trace in self.reference_traces:
                if trace[0] == bit_pos:
                    if trace[1] == 0:
                        hyp0_avg.append(self.find_relevant_power(bit_pos,trace[2]))
                    else:
                        hyp1_avg.append(self.find_relevant_power(bit_pos,trace[2]))
            hyp0_avg = np.mean(hyp0_avg,axis=0)
            hyp1_avg = np.mean(hyp1_avg,axis=0)
            target_relevant_power = self.find_relevant_power(bit_pos,self.target_trace)

            guess, corr = self.determine_ec_key_bit_value_MSE(hyp0_avg, hyp1_avg, target_relevant_power)
            print("Bit {0}: Guess {1}. Correlation Score: {2}".format(bit_pos, guess, corr))
            self.ec_key_guess[bit_pos] = guess

            #plt.figure()
            #plt.plot(hyp0_avg, label="hyp0")
            #plt.plot(hyp1_avg, label="hyp1")
            #plt.plot(target_relevant_power, label="target")
            #plt.title("Bit {0} relevant power traces".format(bit_pos))
            #plt.legend()
            #plt.savefig("/home/40265864@ecit.qub.ac.uk/Toeplitz_attack/results/bit_{0}_relevant_pwr.png".format(bit_pos))
        #plotting.visualise_vary_hypothisis(self.reference_traces, self.target_trace, 0, 7000)

