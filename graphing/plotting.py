import numpy as np
import matplotlib.pyplot as plt

def visualise_power_trace(traces, start=0,end=100):
    """Displays power traces overlayed
    
    :param trace: list of power traces
    start: sample point to begin from
    end: sample point to end at
    """
    for trace in traces:
        plt.plot(trace[start:end])
    plt.show()

##traces should come with bit pos and value
def visualise_vary_hypothisis(ref_traces, target_trace, start=0, end=100):
    
    plt.plot(target_trace[start:end],color="k")
    for trace in ref_traces:
        if trace[1] == 0:
            plt.plot(trace[2][start:end],color="r")
        else:
            plt.plot(trace[2][start:end],color="b")
    plt.show()