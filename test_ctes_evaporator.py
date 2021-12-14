#-*- coding: utf-8 -*-
#! python

# Installed Libs
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import auc

# Local Libs
from physicochemical_properties import *
from evaporators import *
from control import * 
from multi_step import *
from IDEE_CEM import *

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

# Juice Properties (Tjin from csv file)
Mjin = 100.0 #t/h
Mjin = Mjin/3.6 #kg/s
Tjin = 115.0
Bjin = 0.15
Zjin = 0.85
pjin = liquor_prpty.density(Tjin,Bjin,Zjin)
Fjin = (Mjin)/pjin

# Initial juice temperature (out)
Bjout0 = 0.185
Lje0 = 0.30
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
htc = evaporator_roberts.htc_calc(Op, Ne, Bx_out, Tje, Lje0)

# Vapor Properties (Pvin from csv file)
Pvin =  253037.6 # 36.7 psia -> 22 psig
Tvin = vapor_prpty.temperature(Pvin)
# Hvin = vapor_prpty.enthalpy(Tvin,Pvin) - water_prpty.enthalpy(Tvin)
# Mvin = htc*A*(Tvin - Tjout)/Hvin
Mvin = evaporator_roberts.mass_vapor_in(Tjout,Pvin,Tvin,htc,A)
p_vin = (vapor_prpty.density(Pvin))
Fvin = Mvin/p_vin	

Pvv = 184090 # 26.7 psia -> 12 psig
Tvv = vapor_prpty.temperature(Pvv)
Hvv = vapor_prpty.enthalpy(Tvv,Pvv)

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
print("Tjin[°C]: " + str(Tjin))
print("Bjin[kg/kg]: " + str(Bjin))
print("Zjin[kg/kg]: " + str(Zjin))


print("\nVapor in Properties")
print("Pvin[Pa]: " + str(Pvin))
print("Tvin[°C]: " + str(Tvin))
print("Fvin[m3/s]: " + str(Fvin))
print("p_vin[kg/m3]: " + str(p_vin))

print("\nVapor out Properties")
print("Pvv[Pa]: " + str(Pvv))
print("Tvv[°C]: " + str(Tvv))

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

# ============================================================================================
#                                 Evaporator Level Control
# ============================================================================================
print('\n')
print("=================================================================================")
print("                            Evaporator Level Control                             ")
print("=================================================================================")

# # Control setpoints
sp_lvl = 0.33
sp_brx = 0.25

# PID parameters
Kp_brx = -20
Ki_brx = 0.002
Kd_brx = 0
Kp_lvl = 30
Ki_lvl = 0.05
Kd_lvl = 0.05

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

# Control
ctrl_brx = pid(ts, Kp_brx, Ki_brx, Kd_brx)
ctrl_lvl = pid(ts, Kp_lvl, Ki_lvl, Kd_lvl)

# Run Model
x0 = Bjout0, Lje0

# Variable arrays for accumulate and plot
BxBx = np.zeros(t.shape)
LvLv = np.zeros(t.shape)
Fjout_t = np.zeros(t.shape)
Mjout_t = np.zeros(t.shape)
Fjin_t = np.zeros(t.shape)
Mjin_t = np.zeros(t.shape)
Fvin_t = np.zeros(t.shape)
Mvin_t = np.zeros(t.shape)
Mvv_t = np.zeros(t.shape)


for k in range(t.size):

	sol = odeint(evaporator_roberts.model, x0, np.array([0.0, ts]), args=u)
	Bx = sol[1,0]
	Lv = sol[1,1]

	if Lv >= 1:      # level limiter
		Lv = 1
	if Lv <= 0:
		Lv = 0

	x0 = Bx, Lv

	# Accumulate the samples of the variables
	BxBx[k] = Bx
	LvLv[k] = Lv
	Fjout_t[k] = Fjout    
	Fjin_t[k] = Fjin

	# Control loop
	Fjin = ctrl_lvl.solve(sp_lvl, Lv)	  # Fjin is the manipulated variable for level control
	Fjout = ctrl_brx.solve(sp_brx, Bx)    # Fjout is the manipulated variable for Brix control
	# Fvin = ctrl_brx.solve(sp_brx, Bx)     # Fvin is the manipulated variable for Brix control


	if Fjout <= 0:      # Fjout and Fjin limiter to avoid negative values
		Fjout = 0
	
	if Fjin <= 0:
		Fjin = 0

	if Fjin >= Fjin_t[0]*1.4:
		Fjin = Fjin_t[0]*1.4

	if Fvin <= 0:
		Fvin = 0

	if Fvin >= 20:
		Fvin = 20
	
	# # # # Calculate instant mass flow rate in kg/h
	# Fjout [m3/s] -> Mjout [kg/min]
	pjout_k = liquor_prpty.density(Tjout,BxBx[k],Zjin)
	Mjout_t[k] = 60*Fjout_t[k]*pjout_k
	Mjin_t[k] = 60*Fjin_t[k]*pjin

	# # # # Calculate instant htc, Mvv and Mvin
	Uk = evaporator_roberts.htc_calc(Op, Ne, BxBx[k], Tje, LvLv[k])
	Mvin_t[k] = 60*evaporator_roberts.mass_vapor_in(Tjout, Pvin, Tvin, Uk, A)
	Fvin_t[k] = evaporator_roberts.mass_vapor_in(Tjout, Pvin, Tvin, Uk, A)/p_vin
	Mvv_t[k] = 60*evaporator_roberts.mass_vapor_out(Tjin,Zjin,Tjout, (Mjin_t[k]/60), Bjin, BxBx[k], Tvin, Pvv, Uk, A)

	Fvin = Fvin_t[k]

	u = (A, h, At, Ls, Ne, Op, Fjin, Tjin, Bjin, Zjin, Fvin, Tvin, Pvin, Fjout, Tjout, Pvv)

# ============================================================================================
# #   Calculate area under curve to know the amount of processed juice and incoming vapor
# ============================================================================================
print('\n')
print("=================================================================================")
print("                               Simulation Results                                ")
print("=================================================================================")
incoming_juice = auc(tt,Mjin_t)
print("\nINCOMING JUICE IS: "+str(format(incoming_juice/1000, '.3f'))+" [ton] in "+str(tf/60)+ " mins.")
processed_juice = auc(tt,Mjout_t)
print("PROCESSED JUICE IS: "+str(format(processed_juice/1000,'.3f'))+" [ton] in "+str(tf/60)+ " mins.")
incoming_vapor = auc(tt,Mvin_t)
print("INCOMING VAPOR IS: "+str(format(incoming_vapor/1000,'.3f'))+" [ton] in "+str(tf/60)+ " mins.")
vegetal_vapor = auc(tt,Mvv_t)
print("VEGETAL VAPOR IS: "+str(format(vegetal_vapor/1000,'.3f'))+" [ton] "+str(tf/60)+ " mins.\n")

# Calculate IDEE and CEM values
print("IDEE, CEM: ",IDEE_CEM_calc(BxBx,A,Pvin,processed_juice,tf))

# Calculate steady state error
sse_brix = sp_brx - BxBx[t.size-1]
sse_lvl = sp_lvl - LvLv[t.size-1]

print("\nSTEADY STATE ERROR BRIX: ",sse_brix)
print("STEADY STATE ERROR LEVEL: ",sse_lvl,"\n")

# ============================================================================================
#                                       Plotting
# ============================================================================================

f, axarr = plt.subplots(2, 2)
f.suptitle('Roberts Evaporator')
axarr[0,0].plot(tt, BxBx, 'b', label='Brix')
axarr[0,0].plot(tt, sp_brx*np.ones(t.shape), 'r--', label='Setpoint')
axarr[0,0].legend(loc='best')
axarr[0,0].grid()
axarr[0,1].plot(tt, LvLv, 'r', label='Level')
axarr[0,1].plot(tt, sp_lvl*np.ones(t.shape), 'b--', label='Setpoint')
axarr[0,1].legend(loc='best')
axarr[0,1].grid()
axarr[1,0].plot(tt, Mjout_t, 'g', label='Mjout')
axarr[1,0].legend(loc='best')
axarr[1,0].grid()
axarr[1,1].plot(tt, Mjin_t, 'm', label='Mjin')
axarr[1,1].legend(loc='best')
axarr[1,1].grid()

axarr[0,0].set(xlabel = 'Time [min]', ylabel = 'Brix [%]')
axarr[0,1].set(xlabel = 'Time [min]', ylabel = 'Level [%]')
axarr[1,0].set(xlabel = 'Time [min]', ylabel = 'Mjout [kg/min]')
axarr[1,1].set(xlabel = 'Time [min]', ylabel = 'Mjin [kg/min]')
f.canvas.set_window_title('PID control simulation')

fig_csv = plt.figure(3)
plt.plot(tt, Mvv_t, 'g', label='Mvv')
plt.plot(tt, Mvin_t, 'y', label='Mvin')    
plt.legend(loc='best')
plt.grid()
plt.title('Mvin / Mvv from evaporator capture')
plt.xlabel('Time [min]')
plt.ylabel('Mvv [kg/min]')
fig_csv.canvas.set_window_title('Vapor variables') 

plt.show()

