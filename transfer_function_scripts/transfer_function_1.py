import numpy as np
from sklearn.metrics import mean_squared_error
import pandas as pd
from scipy import signal
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from control import *

# Simulate FIRST ORDER TRANSFER FUNCTION taup * dy/dt = -y + K*u
Kp = 1.0
taup = 1.0

# (1) Transfer Function
num = [Kp]
den = [taup,1]
sys1 = signal.TransferFunction(num,den)
t1,y1 = signal.step(sys1)

# Time definition
ts = 0.01	#sampling time [s]
tf = 7	   #simulation time [s]
# tf = 60.0*tf      #-----> samples per minute = 60.0/ts
t = np.linspace(0.0, tf, num=int(tf/ts)+1, endpoint=True) # include last sample
# tt = t/60 # set minutes
tt = t # set seconds


# ODE Integrator
def ODE_model(y,t,u):
    return (-y + Kp * u)/taup

# ctrl = pid(ts, 0.08863, 0.5455, 0.0036)   # control 1_1
ctrl = pid(ts, 1.415, 1.415, 0.0)   # control 1_2
sp = 1
u = 1
y_0 = 0

# print("LEN P: ", len(sp))
Y = np.zeros(t.shape)

for k in range(t.size):
    Y[k] = y_0
    sol = odeint(ODE_model,y_0,np.array([0.0, ts]),args=(u,))
    
    yk = sol[1] 
    # Y[k] = yk
    u = ctrl.solve(sp,yk)
    y_0 = yk

print("Y shape: ", Y.shape)


#### MATLAB response plotting

df = pd.read_csv (r'control_data\control1_2.csv', header = None)
y_matlab = (df.to_numpy()).reshape(Y.shape)

print("LEN MATLAB: ",len(y_matlab))


### CALC RMS ERROR

RMS_ERROR = np.sqrt(((y_matlab-Y)**2).mean())
RMS_ERROR_SKLEARN = mean_squared_error(y_matlab,Y, squared=False)
print("My error [%]: ",RMS_ERROR*100)
print("SKLEARN error [%]: ",RMS_ERROR_SKLEARN*100) 


### PLOTTING

fig = plt.figure(1)

plt.plot(tt,Y,'r-',linewidth=1,label='ODE Integrator')
plt.plot(tt,y_matlab,'b-',linewidth=1,label='MATLAB data')
plt.xlabel('Time [s]')
plt.ylabel('Response (y)')
plt.grid()
plt.legend(loc='best')
plt.title("Python vs Matlab RMS error: " + str(format(RMS_ERROR_SKLEARN*100, '.2f'))+"%")
fig.canvas.set_window_title('FIRST ORDER: PID control response')
plt.show()