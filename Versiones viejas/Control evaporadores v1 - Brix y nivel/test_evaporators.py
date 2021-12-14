#-*- coding: utf-8 -*-
#! python

# Installed Libs
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Local Libs
from physicochemical_properties import *
from evaporators import *
from control import *

liquor_prpty=liquor_properties()
vapor_prpty=vapor_properties()


# ============================================================================================
#                                       Model Initialization
# ============================================================================================
# Time definition
# ts = 1.0		#sampling time
# tf = 100.0		#simulation time
# t = np.arange(0.0, tf, ts)

# Parameters
A = 800.0
h = 2.0
Ne = 1.0
Op = 0.0
Ls = 0.0
At = 10.0 # Transversal Area of juice

# Juice Properties 
Mjin = 100.0 #t/h
Mjin = Mjin/3.6 #kg/s
Tjin = 115.0
Bjin = 0.15
Zjin = 0.85
pjin = liquor_prpty.density(Tjin,Bjin,Zjin)
Fjin = (Mjin)/pjin

# Vapor Properties
Pvin = 253037.6 # 36.7 psia -> 22 psig
Tvin = vapor_prpty.temperature(Pvin)

Pvv = 184090 # 26.7 psia -> 12 psig

# Initial juice temperature (out)
Bjout0 = 0.185
Lje0 = 0.36
Tje = 118.07
Tjout = Tje

# Juice out
Mjout = 0.8*Mjin*Bjin/Bjout0
pjout = liquor_prpty.density(Tjout,Bjout0,Zjin)
Fjout = (Mjout)/pjout

# Model parmeters
u = (A, h, At, Ls, Ne, Op, Fjin, Tjin, Bjin, Zjin, Tvin, Fjout, Tjout, Pvv)

# ============================================================================================
#                                      Test HTC Function
# ============================================================================================
Bx_out = 0.22 # Temporal variable for HTC
#evap = evaporator_roberts(A, h, Ls, Ne, Op)
htc = evaporator_roberts.htc_calc(Op, Ne, Bx_out, Tje, Lje0)

print('\n')
print("=================================================================================")
print("                                  HTC, Roberts                                   ")
print("=================================================================================")
print("Evaporator Design Properties")
print("A[m2]: " + str(A))
print("h[m]: " + str(h))
print("Ne[]: " + str(Ne))
print("Op[d]: " + str(Op))
print("Ls[%]: " + str(Ls))
print("Bjout [%]: " + str(Bx_out*100.0))
print("Tje [%]: " + str(Tje))
print("Lje [%]: " + str(Lje0*100.0))
print("HTC [W/m2.K]: " + str(htc))
print("Fjout [m3/s]: " + str(Fjout))
print("Fjin [m3/s]: " + str(Fjin))

# # ============================================================================================
# #                                       HTC Calc
# # ============================================================================================
# L = np.arange(0.2, 1.0, 0.01)
# htc = evaporator_roberts.htc_calc(Op, Ne, Bx_out, Tje, L)

# # ============================================================================================
# #                                       Plotting
# # ============================================================================================
# plt.plot(L, htc, 'b', label='HTC')
# plt.legend(loc='best')
# plt.title('Heat Transfer Coefficient')
# plt.ylabel('HTC [W/m2.K]')
# plt.xlabel('Level [%]')
# plt.grid()
# plt.show()

# ============================================================================================
#                                      Evaporator Model
# ============================================================================================
print('\n')
print("=================================================================================")
print("                            Evaporator Model, Roberts                            ")
print("=================================================================================")
print("Evaporator Design Properties")
print("A[m2]: " + str(A))
print("h[m]: " + str(h))
print("Ne[]: " + str(Ne))
print("Op[d]: " + str(Op))
print("Ls[%]: " + str(Ls))


print("\nJuice Properties")
print("Mjin[kg/s]: " + str(Mjin)) #100 t/h
print("Tjin[oC]: " + str(Tjin))
print("Bjin[kg/kg]: " + str(Bjin))
print("Zjin[kg/kg]: " + str(Zjin))


print("\nVapor in Properties")
print("Pvin[Pa]: " + str(Pvin))
print("Tvin[oC]: " + str(Tvin))

print("\nVapor out Properties")
print("Pvv[Pa]: " + str(Pvv))
print("Tvv[oC]: " + str(vapor_prpty.temperature(Pvv)))

print("\nInitial Conditions")
print("Bjout [%]: " + str(Bjout0*100.0))
print("Tje [%]: " + str(Tje))
print("Lje [%]: " + str(Lje0*100.0))

# ============================================================================================
#                                       Simulation
# ============================================================================================

# Time definition
ts = 1		#sampling time [s]
tf = 60		#simulation time [min]
tf = 60*tf  #simulation time [s]
t = np.arange(0, tf, ts)
tt = t/60

# Run Model
x0 = Bjout0, Lje0
sol = odeint(evaporator_roberts.model, x0, t, args=u)
Bx = sol[:,0]
Lv = sol[:,1]

# ============================================================================================
#                                       Plotting
# ============================================================================================

f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(tt, Bx, 'b', label='Brix')
axarr[0].legend(loc='best')
axarr[0].grid()
axarr[0].set_title('Evaporator')
axarr[1].plot(tt, Lv, 'r', label='Level')
axarr[1].legend(loc='best')
axarr[1].grid()

axarr[0].set(xlabel = 'Time [min]', ylabel = 'Brix [%]')
axarr[1].set(xlabel = 'Time [min]', ylabel = 'Level [%]')
f.canvas.set_window_title('No control simulation') 
plt.show()

# ============================================================================================
#                                 Evaporator Level Control
# ============================================================================================
print('\n')
print("=================================================================================")
print("                            Evaporator Level Control                             ")
print("=================================================================================")

# Control setpoints
sp_lvl = 0.4
sp_brx = 0.25

# PID parameters
Kp_brx = -5
Ki_brx = 0.3
Kd_brx = 0 
Kp_lvl = 20
Ki_lvl = 3 
Kd_lvl = 1

print("PID Controllers Parameters")
print("Brix control Kp: " + str(Kp_brx))
print("Brix control Ki: " + str(Ki_brx))
print("Brix control Kd: " + str(Kd_brx))
print("Level control Kp: " + str(Kp_lvl))
print("Level control Ki: " + str(Ki_lvl))
print("Level control Kd: " + str(Kd_lvl))

print("\nInitial Conditions")
print("Bjout [%]: " + str(Bjout0*100.0))
print("Lje [%]: " + str(Lje0*100.0))
print("Ts [s]: " + str(ts))

print("\nSetpoint values")
print("Brix sp [%]: " + str(sp_brx))
print("Level sp [%]: " + str(sp_lvl))

# ============================================================================================
#                                       Simulation
# ============================================================================================

# Time definition
# ts = 1.0		#sampling time [s]
# tf = 60.0		#simulation time [min]
# tf = 60.0*tf
# t = np.arange(0.0, tf, ts)
# tt = t/60.0

# Control
ctrl = pid(ts, -2.5, 0.3, 0.0)
ctrl_brx = pid(ts, Kp_brx, Ki_brx, Kd_brx)
ctrl_lvl = pid(ts, Kp_lvl, Ki_lvl, Kd_lvl)

# Run Model
x0 = Bjout0, Lje0

BxBx = np.zeros(t.shape)
LvLv = np.zeros(t.shape)
Fjout_t = np.zeros(t.shape)
Fjin_t = np.zeros(t.shape)

for k in range(t.size):
	sol = odeint(evaporator_roberts.model, x0, np.array([0.0, ts]), args=u)
	Bx = sol[1,0]
	Lv = sol[1,1]

	if Lv >= 1:      # level limiter
		Lv = 1
	if Lv <= 0:
		Lv = 0

	x0 = Bx, Lv

	# if math.isnan(Bx):    
	# 	break

	BxBx[k] = Bx
	LvLv[k] = Lv

	Fjout = ctrl_brx.solve(sp_brx, Bx)    # Fjout is the manipulated variable for Brix control
	Fjin = ctrl_lvl.solve(sp_lvl, Lv)	  # Fjin is the manipulated variable for level control

	# Fjout = ctrl.solve(sp_lvl, Lv)       # Fjout is the manipulated variable for level control
 
	if Fjout <= 0:      # Fjout and Fjin limiter to avoid negative values
		Fjout = 0
	
	if Fjin <= 0:
		Fjin = 0

	Fjout_t[k] = Fjout    
	Fjin_t[k] = Fjin
	
	u = (A, h, At, Ls, Ne, Op, Fjin, Tjin, Bjin, Zjin, Tvin, Fjout, Tjout, Pvv)


# ============================================================================================
#                                       Plotting
# ============================================================================================

f, axarr = plt.subplots(4, sharex=True)
axarr[0].plot(tt, BxBx, 'b', label='Brix')
axarr[0].legend(loc='best')
axarr[0].grid()
axarr[0].set_title('Evaporator')
axarr[1].plot(tt, LvLv, 'r', label='Level')
axarr[1].legend(loc='best')
axarr[1].grid()
axarr[2].plot(tt, Fjout_t, 'g', label='Fjout')
axarr[2].legend(loc='best')
axarr[2].grid()
axarr[3].plot(tt, Fjin_t, 'm', label='Fjin')
axarr[3].legend(loc='best')
axarr[3].grid()

axarr[0].set(xlabel = 'Time [min]', ylabel = 'Brix [%]')
axarr[1].set(xlabel = 'Time [min]', ylabel = 'Level [%]')
axarr[2].set(xlabel = 'Time [min]', ylabel = 'Fjout [m3/s]')
axarr[3].set(xlabel = 'Time [min]', ylabel = 'Fjin [m3/s]')
f.canvas.set_window_title('PID control simulation') 
plt.show()