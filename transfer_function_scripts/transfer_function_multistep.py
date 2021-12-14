##### FIRST ORDER TRANSFER FUNCTION CONTROL SIMULATION WITH MULTISTEP #####

import numpy as np
from sklearn.metrics import mean_squared_error
import pandas as pd
from scipy import signal
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from control import *
from multi_step import *

# Simulate taup * dy/dt = -y + K*u
Kp = 1.0
taup = 1.0

# (1) Transfer Function
num = [Kp]
den = [taup,1]
sys1 = signal.TransferFunction(num,den)
t1,y1 = signal.step(sys1)

# Time definition
ts = 0.1	#sampling time [s]
tf = 20	#simulation time [min]
tf = 60.0*tf      #-----> samples per minute = 60.0/ts
t = np.arange(0.0, tf, ts)
tt = t/60 # set minutes
# tt = t # set seconds


# ODE Integrator
def model3(y,t,u):
    return (-y + Kp * u)/taup

t3 = np.linspace(0,100,100)

tc = np.array([2, 3, 7])
a = np.array([2, 3, 5])

# ctrl = pid(ts, 0.08863, 0.5455, 0.0036)   # control 1_1
ctrl = pid(ts, 1.415, 1.415, 0.0)   # control 1_2c
sp = step_multiple(t, tc, a, ts)
# sp = 1
u = 1
y3_0 = 0

# print("LEN P: ", len(sp))
Y = np.zeros(t.shape)

for k in range(t.size):

    y3 = odeint(model3,y3_0,np.array([0.0, ts]),args=(u,))
    # print("Y3: ",y3)
    
    yk = y3[1]
    # print("yk: ",yk)  
    # print("k: ",k)  
    Y[k] = yk
    # print("U BEFORE es: ",u)
    # u = ctrl.solve(sp[k]+0.5*np.random.rand(1),yk) # with noise
    u = ctrl.solve(sp[k],yk)
    # u = ctrl.solve(sp,yk)
    # print("sp: ",sp)
    y3_0 = yk



plt.figure(1)
# plt.plot(t1,y1,'b--',linewidth=3,label='Transfer Fcn')
# plt.plot(t2,y2,'g:',linewidth=2,label='State Space')
plt.plot(tt,Y,'r-',linewidth=1,label='ODE Integrator')
plt.xlabel('Time [s]')
plt.ylabel('Response (y)')
plt.grid()
plt.legend(loc='best')
plt.show()


# # MATLAB response plotting

# df = pd.read_csv (r'control_data\control1_1.csv', header = None)
# y_matlab = df.to_numpy()
# # print(data[0][1]+1)
# print("LEN MATLAB: ",len(y_matlab[0]))
# t_plot_matlab = np.arange(0.0, len(y_matlab[0]), 1)
# # print(t_plot_matlab)

# plt.figure(1)
# # plt.plot(t1,y1,'b--',linewidth=3,label='Transfer Fcn')
# # plt.plot(t2,y2,'g:',linewidth=2,label='State Space')
# plt.plot(tt,Y,'r-',linewidth=1,label='ODE Integrator')
# plt.xlabel('Time [s]')
# plt.ylabel('Response (y)')
# plt.grid()
# plt.legend(loc='best')
# plt.show()
# plt.plot(t_plot_matlab,y_matlab[0],'b-',linewidth=1,label='MATLAB data')
# plt.xlabel('Samples')
# plt.ylabel('Response (y)')
# plt.grid()
# plt.legend(loc='best')
# plt.show()

# RMS_ERROR = np.sqrt(((y_matlab[0]-Y)**2).mean())
# RMS_ERROR_SKLEARN = mean_squared_error(y_matlab[0],Y, squared=False)
# print("My error [%]: ",RMS_ERROR*100)
# print("SKLEARN error [%]: ",RMS_ERROR_SKLEARN*100) 