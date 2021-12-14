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
from evaporators_v2 import *
from control_2 import *  # using version 2 of the PID controller
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

# Import CSV data (Pvin and Tjin)
Pvin_df = pd.read_csv (r'datos evaporador\pvin_evap.csv', header = None) # CSV data is psig
Pvin_t = (Pvin_df.to_numpy())
Pvin_t = Pvin_t.reshape(len(Pvin_t),)
Pvin_t = (Pvin_t.reshape(len(Pvin_t),)+14.69595)/0.0001450377 # psig to Pa

Tjin_df = pd.read_csv (r'datos evaporador\tjin_evap.csv', header = None) # CSV data is °C
Tjin_t = (Tjin_df.to_numpy())
Tjin_t = Tjin_t.reshape(len(Tjin_t),)

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
Tjin = Tjin_t[0]

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
htc = evaporator_roberts.htc_calc(Op, Ne, Bx_out, Tje, Lje0)

# Vapor Properties (Pvin from csv file)
Pvin = Pvin_t[0]  
Tvin = vapor_prpty.temperature(Pvin)
Hvin = vapor_prpty.enthalpy(Tvin,Pvin) - water_prpty.enthalpy(Tvin)
Mvin = htc*A*(Tvin - Tjout)/Hvin
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
sp_lvl = 0.40
sp_brx = 0.25

# # Noise percent in the multistep variable
# noise_percent = 0.05

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
Tvin_t = np.zeros(t.shape)

# Multi step variable definition
# tc_Bjin = np.array([0, 30, 40, 50, 60])
# a_Bjin = np.array([0.15, 0.16, 0.17, 0.18])
# Bjin_t = step_multiple(tt, tc_Bjin, a_Bjin, ts)
# tc_Tjin = np.array([0, 15, 30, 45, 60])
# a_Tjin = np.array([115, 110, 108, 122])
# Tjin_t = step_multiple(tt, tc_Tjin, a_Tjin, ts)
# tc_Pvin = np.array([0, 20, 30, 40, 60])
# a_Pvin = np.array([253037.6, 260000, 250000, 230000])
# Pvin_t = step_multiple(tt, tc_Pvin, a_Pvin, ts)
# tc_Fvin = np.array([0, 25, 35, 45, 60])
# a_Fvin = np.array([Fvin, Fvin*1.05, Fvin*0.95, Fvin*1.025])
# Fvin_t = step_multiple(tt, tc_Fvin, a_Fvin, ts)



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
	Fvin_t[k] = Fvin
    
    # Update Temperature and Pressure (csv variables)
	Pvin = Pvin_t[k]
	Tjin = Tjin_t[k]
	
	# Bjin = Bjin_t[k]+Bjin_t[0]*float(noise_percent*np.random.uniform(-1,1, (1,)))   # add noise 
	# Bjin_t[k] = Bjin

	# Tjin = Tjin_t[k]+Tjin_t[0]*float(noise_percent*np.random.uniform(-1,1, (1,)))   # add noise 
	
	 
	# Pvin = Pvin_t[k]+Pvin_t[0]*float(noise_percent*np.random.uniform(-1,1, (1,)))   # add noise 
	 
	# Tvin_t[k] = vapor_prpty.temperature(Pvin_t[k])

	# Fvin = Fvin_t[k]+Fvin_t[0]*float(noise_percent*np.random.uniform(-1,1, (1,)))   # add noise 
	# Fvin_t[k] = Fvin
	Fjin = ctrl_lvl.solve(sp_lvl, Lv)	  # Fjin is the manipulated variable for level control
	Fjout = ctrl_brx.solve(sp_brx, Bx)    # Fjout is the manipulated variable for Brix control
	# Fvin = ctrl_brx.solve(sp_brx, Bx)     # Fvin is the manipulated variable for Brix control

	# For multi step
	# Fvin = ctrl_brx.solve(sp_brx, Bx)

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
	Mvin_t[k] = 60*Fvin_t[k]*p_vin

	# # # # Calculate instant htc, Mvv and Mvin
	Uk = evaporator_roberts.htc_calc(Op, Ne, BxBx[k], Tje, LvLv[k])
	Hjin = Tjin*liquor_prpty.heat_capacity(Tjin,Bjin,Zjin)
	Hjout = Tjout*liquor_prpty.heat_capacity(Tjout,BxBx[k],Zjin) # say Zjin = Zjout
	Mvv_t[k] = 60*(Uk*A*(Tvin - Tjout) + (Mjin_t[k]/60)*Hjin - (Mjin_t[k]/60)*Hjout)/(Hvv - Hjout) # kg/min
	 
	u = (A, h, At, Ls, Ne, Op, Fjin, Tjin, Bjin, Zjin, Fvin, Tvin, Pvin, Fjout, Tjout, Pvv)

# ============================================================================================
# #   Calculate area under curve to know the amount of processed juice and incoming vapor
# ============================================================================================

incoming_juice = auc(tt,Mjin_t)
print("\nINCOMING JUICE IS: "+str(format(incoming_juice/1000, '.3f'))+" [ton] in "+str(tf/60)+ " mins.")
processed_juice = auc(tt,Mjout_t)
print("PROCESSED JUICE IS: "+str(format(processed_juice/1000,'.3f'))+" [ton] in "+str(tf/60)+ " mins.")
incoming_vapor = auc(tt,Mvin_t)
print("INCOMING VAPOR IS: "+str(format(incoming_vapor/1000,'.3f'))+" [ton] in "+str(tf/60)+ " mins.")
vegetal_vapor = auc(tt,Mvv_t)
print("VEGETAL VAPOR IS: "+str(format(vegetal_vapor/1000,'.3f'))+" [ton] "+str(tf/60)+ " mins.\n")

# Calculate IDEE and CEM values

# idee, cem = IDEE_CEM_calc(BxBx,A,Pvin,processed_juice,tf)

print("IDEE, CEM: ",IDEE_CEM_calc(BxBx,A,Pvin,processed_juice,tf))

# ============================================================================================
#                                       Plotting
# ============================================================================================

f, axarr = plt.subplots(2, 2)
f.suptitle('Roberts Evaporator')
axarr[0,0].plot(tt, BxBx, 'b', label='Brix')
axarr[0,0].legend(loc='best')
axarr[0,0].grid()
axarr[0,1].plot(tt, LvLv, 'r', label='Level')
axarr[0,1].legend(loc='best')
axarr[0,1].grid()
axarr[1,0].plot(tt, Mjout_t, 'g', label='Mjout')
axarr[1,0].legend(loc='best')
axarr[1,0].grid()
axarr[1,1].plot(tt, Fjin_t, 'm', label='Fjin')
axarr[1,1].legend(loc='best')
axarr[1,1].grid()

axarr[0,0].set(xlabel = 'Time [min]', ylabel = 'Brix [%]')
axarr[0,1].set(xlabel = 'Time [min]', ylabel = 'Level [%]')
axarr[1,0].set(xlabel = 'Time [min]', ylabel = 'Mjout [kg/min]')
axarr[1,1].set(xlabel = 'Time [min]', ylabel = 'Fjin [m3/s]')
f.canvas.set_window_title('PID control simulation') 

fig_csv, ax_csv = plt.subplots(2)
ax_csv[0].plot(tt, Tjin_t, 'g', label='Tjin') 
ax_csv[0].legend(loc='best')
ax_csv[0].grid()
ax_csv[0].set_title('Tjin / Pvin from evaporator capture')
ax_csv[0].set(xlabel = 'Time [min]', ylabel = 'Tjin [°C]')
ax_csv[1].plot(tt, (Pvin_t*0.0001450377)-14.69595, 'y', label='Pvin')    
ax_csv[1].legend(loc='best')
ax_csv[1].grid()
ax_csv[1].set(xlabel = 'Time [min]', ylabel = 'Pvin [psig]')
fig_csv.canvas.set_window_title('CSV variables') 

plt.show()

