# Installed Libs

import math
import numpy as np
# import matplotlib.pyplot as plt

# t: time vector 
# t: time changes vector [minutes] 
# ts: sample time
# a: amplitude vector for time changes

def step_multiple(t,tc,a,ts):
    var_multi_step = np.ones(len(t))
    tc = tc*(60/ts)
    print(tc)
    for i in range (len(tc)-1):
        # print("i: ",i)
        var_multi_step[int(tc[i]):int(tc[i+1])] = var_multi_step[int(tc[i]):int(tc[i+1])]*a[i]
        # print(var_multi_step)
    print("im working")
    return var_multi_step



# Time definition
ts = 0.1		#sampling time [s]
tf = 35	#simulation time [min]
tf = 60.0*tf
t = np.arange(0.0, tf, ts)

tc = np.array([10, 20, 30, 32])
a = np.array([2, 3, 4])

var_multi_step = np.ones(len(t))

# print(tc[0])
# print(tc[1])
# print(a[1])
# print(var_multi_step[tc[0]:tc[1]]*a[1])
# var_multi_step[tc[0]:tc[1]] = var_multi_step[tc[0]:tc[1]]*a[1]
# print(var_multi_step)
print(step_multiple(t,tc,a,ts))


# plt.figure(1)
# # plt.plot(t1,y1,'b--',linewidth=3,label='Transfer Fcn')
# # plt.plot(t2,y2,'g:',linewidth=2,label='State Space')
# plt.plot(t,step_multiple(t,tc,a,ts),'r-',linewidth=1,label='ODE Integrator')
# plt.xlabel('Time [s]')
# plt.ylabel('Response (y)')
# plt.grid()
# plt.legend(loc='best')
# plt.show()