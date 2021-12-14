import numpy as np
from sklearn.metrics import mean_squared_error
import pandas as pd
from scipy import signal
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from control_2 import *
from multi_step import *


# Simulate SECOND ORDER TRANSFER FUNCTION
K = 1.0
wn = 1.0
ro = 1.0

# Time definition
tf = 60   #simulation time [s]
ts = 0.01	#sampling time [s]
# tf = 60.0*tf      #-----> samples per minute = 60.0/ts
t = np.linspace(0.0, tf, num=int(tf/ts)+1, endpoint=True) # include last sample
# t = np.arange(0.0, tf, ts)
# tt = t/60 # set minutes
tt = t # set seconds


# ODE Integrator
def model_2nd_order(x,t,u):
    y, theta = x
    dxdt = [theta, K*wn*wn*u-2*ro*wn*theta-wn*wn*y]    # TF: k*wn^2/(s^2+2*ro*wn*s+wn^2) 
    return dxdt


tc = np.array([2, 3, 7])
a = np.array([2, 3])

ctrl = pid(ts, 8.864, 0.4351, 0.0)   # control 2_1
# ctrl = pid(ts, 12.8941, 8.8302, 4.707)   # control 2_2
# sp = step_multiple(t, tc, a, ts)
sp = 1
u = 1
x_0 = 0, 0

# print("LEN P: ", len(sp))
Y = np.zeros(t.shape)

for k in range(t.size):
    Y[k] = x_0[0]
    sol = odeint(model_2nd_order,x_0,np.array([0.0, ts]),args=(u,))
    
    yk = sol[1,0]
    thetak = sol[1,1]
 
    # Y[k] = yk
    
    # u = ctrl.solve(sp[k],yk)  # for multi step
    u = ctrl.solve(sp,yk)
    
    x_0 = yk, thetak  

# print("Y is: ", Y)
print("Y len: ", len(Y))


# # MATLAB response plotting

df = pd.read_csv (r'control_data\control2_1.csv', header = None)
y_matlab = (df.to_numpy()).reshape(Y.shape)
no_samples = len(y_matlab)
print("LEN MATLAB: ",no_samples)

### RMS CALCULATION

RMS_ERROR = np.sqrt(((y_matlab-Y)**2).mean())
RMS_ERROR_SKLEARN = mean_squared_error(y_matlab,Y, squared=False)
print("My error [%]: ",RMS_ERROR*100)
print("SKLEARN error [%]: ",RMS_ERROR_SKLEARN*100) 


### PLOTTING 

fig = plt.figure(1)
plt.plot(tt,Y,'r-',linewidth=1,label='Control Python')
plt.plot(tt,y_matlab,'b-',linewidth=1,label='Control MATLAB')
plt.xlabel('Time [s]')
plt.ylabel('Response (y)')
plt.grid()
plt.legend(loc='best')
plt.title("Python vs Matlab RMS error: " + str(format(RMS_ERROR_SKLEARN*100, '.2f'))+"%")
fig.canvas.set_window_title('SECOND ORDER: PID control response')
plt.show()


