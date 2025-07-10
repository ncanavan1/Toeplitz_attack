##########Contains the tools and algorithms for the template attacks####
import chipwhisperer as cw
import numpy as np 
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error 
from scipy import stats
import sys
import time
import scipy.signal

class Tools:
    def __init__(self, scope, target):
        self.scope = scope
        self.target = target


    def plot_overlay(self, traces):

        fig = cw.plot()
        for trace in traces:
            fig *= cw.plot(trace)
        return fig
    
    def plot_difference(self, trace1, trace2):
        return cw.plot(trace1 - trace2) 


    def align_traces(self, traces):
        ref_trace = traces[0]  # Use the first trace as reference
        aligned_traces = []

        for trace in traces:
            correlation = np.correlate(trace, ref_trace)#, mode="full")  # Compute cross-correlation
            shift = np.argmax(correlation) - (len(trace) - 1)  # Find best alignment
            aligned_trace = np.roll(trace, -shift)  # Shift the trace
            aligned_traces.append(aligned_trace)
        
        return np.array(aligned_traces)

    def get_trace(self, x):
        #num_char = target.in_waiting()
        #while num_char > 0:
            #target.read(num_char, 10)
            #time.sleep(0.01)
            #num_char = target.in_waiting()
        time.sleep(0.01)
        #target.flush()
        self.scope.arm()
        self.target.write(x)
        if self.scope.capture():
            raise RuntimeError("Capture failed")
        trace_segments = self.scope.get_last_trace_segmented()
        alligned_traces = self.align_traces(trace_segments)
        return alligned_traces
    
    ##################################
    ### USED IN FORMAL MGD METHOD ####
    ##################################

    def calculate_mean_sample(self, samples):
        return np.mean(samples,axis=0)
    
    
    def compute_noise_vector(self, sample, mean, points):
        n_vec = []
        for p in points:
            n_vec.append(np.abs(sample[p] - mean[p]))
        return np.asarray(n_vec)
    
    def compute_noise_covarience_matrix(self, samples, mean, points):
        mat_size = points.shape[0]
        cov_mat = np.zeros([mat_size, mat_size])
        p = 0
        for u in points:
            q = 0
            for v in points:
                n1 = self.compute_noise_vector(u, mean, points)
                n2 = self.compute_noise_vector(v, mean, points)
                cov =  np.cov(n1,n2)[0][1]
                cov_mat[p,q] = cov
                q=q+1
            p=p+1
        return cov_mat
    
 
    def compress_samples(self, sample_matrix_in, points):
        sample_size = sample_matrix_in.shape[0]
        sample_matrix_out = np.zeros([sample_size, points.shape[0]])

        for s in range(sample_size):
            q = 0
            for p in points:
                sample_matrix_out[s,q] = sample_matrix_in[s,p]
                q = q+1
        return sample_matrix_out
    
    #assumes that input matrix is already reduced to dimension of use
    def compute_noise_cov_mat(self, sample_matrix, mean):
        sample_size = sample_matrix.shape[0]
        sample_points = sample_matrix.shape[1]
        cov_mat = np.zeros([sample_points, sample_points])

        for i in range(sample_size):
            n_vec = sample_matrix[i] - mean
            mat = np.matmul(n_vec.T,n_vec)
            cov_mat = cov_mat + mat

        return (1/(sample_size-1))*cov_mat
    
    def calc_p_vec(self, n, det_sig, inv_sig):
        N = n.shape[0]
        e1 = np.matmul(n.T,inv_sig)
        e = np.exp(-0.5*np.matmul(e1,n))
        s = 1/(np.sqrt(2*np.pi)**N * det_sig)
        return s*e
    
    def compute_p_vec(self, sample, means, points, dets, invs):
        test_vec0 = self.compute_noise_vector(sample, means[0], points)
        test_vec1 = self.compute_noise_vector(sample, means[1], points)
        test_vec2 = self.compute_noise_vector(sample, means[2], points)
        test_vec3 = self.compute_noise_vector(sample, means[3], points)

        p_0 = self.calc_p_vec(test_vec0, dets[0], invs[0])
        p_1 = self.calc_p_vec(test_vec1, dets[1], invs[1])
        p_2 = self.calc_p_vec(test_vec2, dets[2], invs[2])
        p_3 = self.calc_p_vec(test_vec3, dets[3], invs[3])

        p_sum = p_0 + p_1 + p_2 + p_3
        p_0 = p_0/p_sum
        p_1 = p_1/p_sum
        p_2 = p_2/p_sum
        p_3 = p_3/p_sum


        return np.asarray([p_0, p_1, p_2, p_3])
    
    
    def compute_pearson_vec(self, sample, means):
        p0, d0 = stats.pearsonr(sample, means[0])
        p1, d1 = stats.pearsonr(sample, means[1])
        p2, d2 = stats.pearsonr(sample, means[2])
        p3, d3 = stats.pearsonr(sample, means[3])

        d0 = 1 - d0
        d1 = 1 - d1
        d2 = 1 - d2
        d3 = 1 - d3

        sum = d0 + d1 + d2 + d3
        d0 = d0/sum
        d1 = d1/sum
        d2 = d2/sum
        d3 = d3/sum

        #sum = p0 + p1 + p2 + p3
        #p0 = p0/sum
        #p1 = p1/sum
        #p2 = p2/sum
        #p3 = p3/sum




        return np.asarray([d0,d1,d2,d3])
    
    def class_accuracy_exp(self, samples, points, means, dets, invs, L):
        classes = 4
        Confusion_matrix = np.zeros([classes,classes])
        confidence_matrix = np.zeros([classes,classes])
        ##4 classes
        for k in range(classes):
            for l in range(L):
                p_vec = self.compute_p_vec(samples[k*L + l], means, points, dets, invs)
                #p_vec = compute_pearson_vec(samples[k*L + l], means)
                winner = np.argmax(p_vec)
                Confusion_matrix[k,winner] = Confusion_matrix[k,winner]+1
                confidence_matrix[k] = confidence_matrix[k] + p_vec
        confidence_matrix = confidence_matrix/L
        return Confusion_matrix, confidence_matrix
    

    ################################
    ### USED IN INFORMAL METHOD ####
    ################################


    def array_to_bin_input(self, arr):
        op = ""
        for x in arr:
            op = op + str(int(x))
        op = op + "\n"
        return op


    def infer_from_traces(self, target_trace, template_traces, segment):
        ##built for 2 bit guesser

        MSE = []
        PC = []
        print("\n\nSegment {0}".format(segment))
        for i in range(4):
            MSE.append(mean_squared_error(target_trace[segment], template_traces[i][segment]))
            PC.append(stats.pearsonr(target_trace[segment],template_traces[i][segment]))
            print("Guess {0}. MSE: {1}, PC: {2}".format(i,MSE[i], PC[i]))
            
        best = np.argsort(MSE)[0]
        result = []
        match best:
            case 0:
                result = [0,0]
            case 1:
                result = [0,1]
            case 2:
                result = [1,0]
            case 3:
                result = [1,1]
        
        return result


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

    def gen_pairs(self, N):
        pairs = []
        order_list = np.arange(N)
        reverse_order = self.bit_reverse(order_list)
        for i in range(int(N/2)):
            pairs.append([reverse_order[int(2*i)], reverse_order[int(2*i + 1)]])
        return pairs


    def guess_sequentially_by_2(self, target_trace, L, N, sigma):
        ##modify to generate this on function input N (lenght of FFT)
        pairs = self.gen_pairs(N)
        key_hyp = np.zeros(N)
        segment = 0

        for pair in pairs:

            self.reset_target()
            for _ in range(1):
                warmup_str = self.array_to_bin_input(np.ones(N))
                warmup = self.get_trace(warmup_str)

            traces_avg = []
            guesses = [[0,0], [0,1], [1,0], [1,1]]
            for guess in guesses:
                key_hyp[pair[0]] = guess[0]
                key_hyp[pair[1]] = guess[1]
                hyp_string = self.array_to_bin_input(key_hyp)
                traces = []
                for _ in range(L):
                    ## changed this line in 6th may to include noise in Eves device
                    traces.append(self.get_trace(hyp_string) + np.random.normal(0,sigma,self.scope.adc.samples-1))
                traces_avg.append(np.mean(traces,axis=0))
            
            result = self.infer_from_traces(target_trace, traces_avg, segment)
            key_hyp[pair[0]] = result[0]
            key_hyp[pair[1]] = result[1]
            segment += 1
        return key_hyp


    def run_attack(self, L, sigma, N, seg_count):
        self.reset_target()
        correct_key = np.random.randint(2,size=N)
        correct_key_str = self.array_to_bin_input(correct_key)
        print(correct_key)

        for _ in range(5):
            warmup_str = self.array_to_bin_input(np.ones(N))
            warmup = self.get_trace(warmup_str)

        target_trace = self.get_trace(correct_key_str)
        trace_pre = target_trace.copy()
        ##add noise to target trace
        for seg in range(seg_count):
                target_trace[seg] = target_trace[seg] + np.random.normal(0,sigma,self.scope.adc.samples-1)


        start = time.time()
        sigma_Eve = 0
        guess = self.guess_sequentially_by_2(target_trace, L, N, sigma_Eve)
        end = time.time()

        time_taken = end - start

        guess = np.asarray(guess,dtype=int)
        correct_key = np.asarray(correct_key,dtype=int)
        print("\n{0} : Eve's Guess".format(guess))
        print("{0} : Correct Key".format(correct_key))
        cor_count = 0
        for i in range(N):
            if guess[i] == correct_key[i]:
                cor_count = cor_count+1
            
        if(guess == correct_key).all():
            print("Successfull Key Recovery")
        else:
            print("Incorrect Key Recovery, {0}% Correct".format(cor_count/N *100))
        print("Time Taken: {0}s".format(time_taken))

        return cor_count, time_taken
    
    def reset_target(self):
        self.scope.io.nrst = False
        time.sleep(0.001)
        self.scope.io.nrst = True
        time.sleep(0.001)
