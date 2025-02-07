import numpy as np

class Toeplitz_hashing:
    def __init__(self, pa_len, ec_key, l_key, sigma=0.1, TM_coeff_row=[], TM_coeff_col=[]):
        self.pa_len = pa_len ##output len (number of rows)
        self.pa_key = np.zeros(self.pa_len)
        self.ec_key = ec_key
        self.l_key = l_key ##input len (number of coloumns)
        if len(TM_coeff_row) == 0 and len(TM_coeff_col) == 0:
            self.TM_coeff_row = np.zeros(self.l_key,dtype=int)
            self.TM_coeff_col = np.zeros(self.pa_len,dtype=int)
            self.gen_TM_matrix_coeff()
        else:
            self.TM_coeff_row = TM_coeff_row
            self.TM_coeff_col = TM_coeff_col
        self.sigma = sigma
        self.simulated_power = []


    def gen_TM_matrix_coeff(self):
        ##To generate toplitz matrix, only need to know values fromfirst row and first coloumn
        for i in range(self.l_key):
            self.TM_coeff_row[i] = np.random.randint(0,2,1)

        self.TM_coeff_col[0] = self.TM_coeff_row[0]

        for i in range(1,self.pa_len):
            self.TM_coeff_col[i] = np.random.randint(0,2,1)


    
    def gen_explicit_TM(self):

        TM = np.ones([self.pa_len,self.l_key]) ##place holder matrix
        for j in range(self.pa_len):
            for i in range(self.l_key):
                if i >= j:
                    TM[j,i] = self.TM_coeff_row[i-j] 
                else:
                    TM[j,i] = self.TM_coeff_col[j-i]
        
        return TM

    ##Need a more explicit version that will simulate power
    def calc_pa_key(self):
        TM = self.gen_explicit_TM()
        self.pa_key = np.matmul(TM,self.ec_key)%2

    def calc_pa_explicit(self):
        ##must make sure that calc_pa_key has bden called first such that explicit version can be verified.
        result = np.zeros(self.pa_len, dtype=int)

        for j in range(self.pa_len):
            for i in range(self.l_key):
                intermediate_val = 0
                if i >= j:
                    intermediate_val = int(self.TM_coeff_row[i-j]) & int(self.ec_key[i])
                    result[j] ^= intermediate_val
                else:
                    intermediate_val = int(self.TM_coeff_col[j-i]) & int(self.ec_key[i])
                    result[j] ^= intermediate_val
                
                ##power model
                self.simulated_power.append(intermediate_val + np.random.normal(0,self.sigma))

        self.pa_key = result
        
    

