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
            correlation = np.correlate(trace, ref_trace, mode="full")  # Compute cross-correlation
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
        time.sleep(0.1)
        #target.flush()
        self.scope.arm()
        self.target.write(x)
        if self.scope.capture():
            raise RuntimeError("Capture failed")
        trace_segments = self.scope.get_last_trace_segmented()
        alligned_traces = self.align_traces(trace_segments)
        return alligned_traces

""""
def array_to_bin_input(arr):
    op = ""
    for x in arr:
        op = op + str(int(x))
    op = op + "\n"
    return op


def infer_from_traces(target_trace, template_traces, segment):
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

def gen_pairs(N):
    pairs = []
    order_list = np.arange(N)
    reverse_order = bit_reverse(order_list)
    for i in range(int(N/2)):
        pairs.append([reverse_order[int(2*i)], reverse_order[int(2*i + 1)]])
    return pairs


def guess_sequentially_by_2(target_trace, repeats, N):
    ##modify to generate this on function input N (lenght of FFT)
    pairs = gen_pairs(N)
    key_hyp = np.zeros(N)
    segment = 0

    for pair in pairs:

        reset_target(scope)
        for _ in range(5):
            warmup_str = array_to_bin_input(np.ones(N))
            warmup = get_trace(warmup_str)

        traces_avg = []
        guesses = [[0,0], [0,1], [1,0], [1,1]]
        for guess in guesses:
            key_hyp[pair[0]] = guess[0]
            key_hyp[pair[1]] = guess[1]
            hyp_string = array_to_bin_input(key_hyp)
            traces = []
            for r in range(repeats):
                traces.append(get_trace(hyp_string))
            traces_avg.append(np.mean(traces,axis=0))
        
        result = infer_from_traces(target_trace, traces_avg, segment)
        key_hyp[pair[0]] = result[0]
        key_hyp[pair[1]] = result[1]
        segment += 1
    return key_hyp


def run_attack(repeats, sigma):
    reset_target(scope)
    correct_key = np.random.randint(2,size=N)
    correct_key_str = array_to_bin_input(correct_key)
    print(correct_key)

    for _ in range(5):
        warmup_str = array_to_bin_input(np.ones(N))
        warmup = get_trace(warmup_str)

    target_trace = get_trace(correct_key_str)
    trace_pre = target_trace.copy()
    ##add noise to target trace
    for seg in range(seg_count):
            target_trace[seg] = target_trace[seg] + np.random.normal(0,sigma,scope.adc.samples-1)


    start = time.time()
    guess = guess_sequentially_by_2(target_trace, 1, N)
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

    return cor_count
"""