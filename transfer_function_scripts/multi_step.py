# Installed Libs
import math
import numpy as np

# t: time vector 
# tc: time changes vector [minutes] 
# ts: sample time
# a: amplitude vector for time changes

def step_multiple(t,tc,a,ts):
    var_multi_step = np.ones(len(t))
    tc = tc*(60/ts)
    
    for i in range (len(tc)-1):
        var_multi_step[int(tc[i]):int(tc[i+1])] = var_multi_step[int(tc[i]):int(tc[i+1])]*a[i]
        
    return var_multi_step
