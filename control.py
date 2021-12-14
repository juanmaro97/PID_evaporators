#-*- coding: utf-8 -*-
#! python

# Installed Libs
import math
import numpy as np

class pid:
    '''
    Parameters:
    sp, Setpoint
    pv, Process Variable
    mv, Manipulated Variable
    '''

    def __init__(self, ts, kp, ki, kd):
        # Create pid properties
        self.ts = ts
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.error0 = 0.0
        self.i0 = 0.0

    def solve(self, sp, pv):
        # Mohammad Saadeh (2021). MATLAB based PID controller 
        # (https://www.mathworks.com/matlabcentral/fileexchange/35163-matlab-based-pid-controller), 
        # MATLAB Central File Exchange. Retrieved November 16, 2021.
        
        self.sp = sp
        self.pv = pv
        self.error = self.sp - self.pv

        prop = self.error
        i = 0.5*(self.error + self.error0)*self.ts
        integral = i + (self.i0)
        derivative = (self.error - self.error0)/self.ts

        self.mv = self.kp*prop + self.ki*integral + self.kd*derivative
        self.error0 = self.error # saving error for next loop
        self.i0 = integral
        # print("ERROR: ",self.error)
        return self.mv


        