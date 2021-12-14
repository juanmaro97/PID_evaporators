#-*- coding: utf-8 -*-
#! python

# Installed Libs
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Local Libs
from physicochemical_properties import *
from evaporators_v2 import *
from control_2 import *  # using version 2 of the PID controller
from multi_step import *

liquor_prpty=liquor_properties()
vapor_prpty=vapor_properties()


# ============================================================================================
#                                       Model Initialization
# ============================================================================================
# Time definition
ts = 1		#sampling time [s]
tf = 60		#simulation time [min]
tf = 60*tf  #simulation time [s]
t = np.arange(0, tf, ts)
tt = t/60

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

# variate Bjin
Bjin = 0.15
Zjin = 0.85
pjin = liquor_prpty.density(Tjin,Bjin,Zjin)
Fjin = (Mjin)/pjin

# Initial juice temperature (out)
Bjout0 = 0.185
Lje0 = 0.36
Tje = 118.07
Tjout = Tje

# Juice out
Mjout = 0.8*Mjin*Bjin/Bjout0
pjout = liquor_prpty.density(Tjout,Bjout0,Zjin)
Fjout = (Mjout)/pjout

# ============================================================================================
#                                      Test HTC Function
# ============================================================================================
Bx_out = 0.22 # Temporal variable for HTC
#evap = evaporator_roberts(A, h, Ls, Ne, Op)
htc = evaporator_roberts.htc_calc(Op, Ne, Bx_out, Tje, Lje0)

# Vapor Properties
Pvin = 253037.6 # 36.7 psia -> 22 psig
Tvin = vapor_prpty.temperature(Pvin)
Hvin = vapor_prpty.enthalpy(Tvin,Pvin) - water_prpty.enthalpy(Tvin)
Mvin = htc*A*(Tvin - Tjout)/Hvin
p_vin = (vapor_prpty.density(Pvin))
Fvin = Mvin/p_vin	

Pvv = 184090 # 26.7 psia -> 12 psig

# Model parmeters
u = (A, h, At, Ls, Ne, Op, Fjin, Tjin, Bjin, Zjin, Fvin, Tvin, Pvin, Fjout, Tjout, Pvv)


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
print("Fvin[m3/s]: " + str(Fvin))

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

# # Control setpoints
sp_lvl = 0.40
# sp_brx = 0.25

# sp_lvl = step_multiple(tt, tc_lvl, a_lvl, ts)
sp_brx = 0.25

# PID parameters
Kp_brx = -20
Ki_brx = 0.002
Kd_brx = 0
Kp_lvl = 10
Ki_lvl = 1
Kd_lvl = 0.5

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
print("Brix sp [%]: " + str(sp_brx*100.0))
print("Level sp [%]: " + str(sp_lvl*100.0))

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
# ctrl = pid(ts, -2.5, 0.3, 0.0)
ctrl_brx = pid(ts, Kp_brx, Ki_brx, Kd_brx)
ctrl_lvl = pid(ts, Kp_lvl, Ki_lvl, Kd_lvl)

# Run Model
x0 = Bjout0, Lje0

BxBx = np.zeros(t.shape)
LvLv = np.zeros(t.shape)
Fjout_t = np.zeros(t.shape)
Fjin_t = np.zeros(t.shape)
Fvin_t = np.zeros(t.shape)

tc_Bjin = np.array([0, 30, 40, 50, 60])
a_Bjin = np.array([0.15, 0.17, 0.13, 0.15])
Bjin_t = step_multiple(tt, tc_Bjin, a_Bjin, ts)

for k in range(t.size):
	# print("k is: ",k)
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

	# Accumulate the samples of the variables
	BxBx[k] = Bx
	LvLv[k] = Lv
	Fjout_t[k] = Fjout    
	Fjin_t[k] = Fjin
	Fvin_t[k] = Fvin
	Bjin = Bjin_t[k]

	Fjin = ctrl_lvl.solve(sp_lvl, Lv)	  # Fjin is the manipulated variable for level control
	Fjout = ctrl_brx.solve(sp_brx, Bx)    # Fjout is the manipulated variable for Brix control
	# Fvin = ctrl_brx.solve(sp_brx, Bx)
	# Fjin = ctrl_lvl.solve(sp_lvl, Lv)	  # Fjin is the manipulated variable for level control

	# For multi step
	# Fvin = ctrl_brx.solve(sp_brx, Bx)
	

	# Fjout = ctrl.solve(sp_lvl, Lv)       # Fjout is the manipulated variable for level control
 
	if Fjout <= 0:      # Fjout and Fjin limiter to avoid negative values
		Fjout = 0
	
	if Fjin <= 0:
		Fjin = 0

	if Fvin <= 0:
		Fvin = 0

	if Fvin >= 20:
		Fvin = 20
	
	u = (A, h, At, Ls, Ne, Op, Fjin, Tjin, Bjin, Zjin, Fvin, Tvin, Pvin, Fjout, Tjout, Pvv)


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
axarr[3].plot(tt, Bjin_t, 'm', label='Bjin')
axarr[3].legend(loc='best')
axarr[3].grid()
# axarr[4].plot(tt, Bjin_t, 'b', label='Bjin')
# axarr[4].legend(loc='best')
# axarr[4].grid()

axarr[0].set(xlabel = 'Time [min]', ylabel = 'Brix [%]')
axarr[1].set(xlabel = 'Time [min]', ylabel = 'Level [%]')
axarr[2].set(xlabel = 'Time [min]', ylabel = 'Fjout [m3/s]')
axarr[3].set(xlabel = 'Time [min]', ylabel = 'Bjin [%]')
# axarr[4].set(xlabel = 'Time [min]', ylabel = 'Bjin [%]')
f.canvas.set_window_title('PID control simulation') 
plt.show()