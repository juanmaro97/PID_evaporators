#-*- coding: utf-8 -*-
#! python

# Installed Libs
import math
import numpy as np

# Local Libs
from physicochemical_properties import *

liquor_prpty = liquor_properties()
water_prpty = water_properties()
vapor_prpty = vapor_properties()

class evaporator_roberts:
	'''
	Parameters:
	A, Heat transfer area [m2]
	h, Calandria heigth [m]
	Np, Number of pipes []
	Di, Pipe internal diamenter [in]
	Dd, Downtake diameter [m]
	hc, Inferior cone heigth [m]
	Ne, Effect number []
	Op, Operation days [d]
	Ls, Heat losses [%]
	Lje, height of juice in the evaporator [] 0 to 1
	'''

	def __init__(self, A, h, N, Di, Dd, hc, Ne, Op, Ls, Lje):
		# Create evaporator desing properties
		self.A = A
		self.h = h
		self.Np = Np
		self.Di = Di
		self.Dd = Dd
		self.hc = hc
		self.Ne = Ne
		self.Op = Op
		self.Ls = Ls
		self.Lje = Lje

		self.At = math.pi*(Dd/2)**2 + Np*math.pi*(Di/2)**2

	def in_out(self, juice_in, vapor_in, juice_out, vapor_out, condensate_out):
		self.juice_in = juice_in
		self.vapor_in = vapor_in
		self.juice_out = juice_out
		self.vapor_out = vapor_out
		self.condensate_out = condensate_out

	def solve(self, time):

		return self.juice_out, self.vapor_out, self.condensate_out


	@staticmethod
	def htc_aus(Bjout, Tje): # Australian typical OHTC for roberts evaporators, Wright 2008
		# print("Bjout for HTC is: ",Bjout)
		# print("Tje for HTC is: ",Tje)
		htc = 16.9390445*(Tje**1.0174)*(Bjout/(0.86 - Bjout))**-0.2695352
		return htc

	@staticmethod
	def htc_scaling(Op, Ne): # Scale in evaporators based on Rosero 2008
		scale = {1: 0.0625,
					2: 0.0556,
					3: 0.0498,
					4: 0.0376,
					5: 0.0626,
					} #dictionary with scale values depending of effect
		htc = math.exp(-scale[Ne]*Op)
		return htc

	@staticmethod
	def htc_level(Lje): # HTC factor depending on level, regression data from Hugot 1977
		htc = -6.1936*Lje**4 + 18.414*Lje**3 - 19.869*Lje**2 + 8.6037*Lje - 0.277
		return htc

	@classmethod # class method to be used with or without instance
	def htc_calc(cls, Op, Ne, Bjout, Tjout, Lje):
		htc = cls.htc_aus(Bjout, Tjout)*cls.htc_scaling(Op, Ne)*cls.htc_level(Lje)
		return htc

	# Internal use only, call in_out(...) method first
	def htc(self):
		Bjout = self.juice_out.Bj
		Tjout = self.juice_out.Tj

		htc = self.htc_aus(Bjout, Tjout)*self.htc_scaling(self.Op, self.Ne)*self.htc_level(self.Lje)
		return htc

	def mass_vapor_in(self): 
		Tjout = self.juice_out.Tj
		Pvin = self.vapor_in.Pv
		Tvin = self.vapor_in.Tv

		Hvin = vapor_prpty.enthalpy(Tvin,Pvin) - water_prpty.enthalpy(Tvin)
		U = self.htc()

		Mvin = U*self.A*(Tvin - Tjout)/Hvin
		return Mvin

	def mass_vapor_out(self):
		U = self.htc()
		Mjin = self.juice_in.Mj
		Hjin = self.juice_in.Hj

		Hjout = self.juice_out.Hj

		Tvin = self.vapor_in.Tv
		Tjout = self.vapor_out.Tj
		Hvv = self.vapor_out.Hv

		Mvv = (U*self.A*(Tvin - Tjout) + Mjin*Hjin - Mjin*Hjout)/(Hvv - Hjout)
		return Mvv

	def hydrostatic_pressure(self):
		Tje = self.juice_out.Tj
		Bje = self.juice_out.Bj
		Zje = self.juice_out.Zj

		p = liquor_prpty.density(Tje,Bje,Zje)
		g = 9.7775 # Gravity in Cali, by latitude and height
		h = self.h + self.heat_capacity

		self.P_hyd = p*g*h
		return self.P_hyd

	def residence_time(self):
		tr = 0
		return 0

	def sucrose_losses(self):
		#liquor_prpty.sucrose_losses()  
		#sucrose_losses(self, time, Temperature, Brix, SolIn, Purity, pH)
		return 0

	@staticmethod # static method to be used in scipy.integrate.odeint
	def model(x, t, A, h, At, Ls, Ne, Op,
						 Fjin, Tjin, Bjin, Zjin, Tvin, Fjout, Tjout, Pvv):

		# Initial conditions
		Bjout, Lje = x

		# Vapor in
		# Tvin = vapor_prpty.temperature(Pvin)
		# Hvin = vapor_prpty.enthalpy(Tvin,Pvin) - water_prpty.enthalpy(Tvin)

		# Vapor out
		Tvv = vapor_prpty.temperature(Pvv)
		Hvv = vapor_prpty.enthalpy(Tvv,Pvv) 

		# Heat capacity
		Cpjin = liquor_prpty.heat_capacity(Tjin,Bjin,Zjin)
		Cpjout = liquor_prpty.heat_capacity(Tjout,Bjout,Zjin)

		# Density
		pjin = liquor_prpty.density(Tjin,Bjin,Zjin)
		pjout = liquor_prpty.density(Tjout,Bjout,Zjin)

		# Heat transfer coefficient
		U = evaporator_roberts.htc_calc(Op, Ne, Bjout, Tjout, Lje)

		# Brix differential equation
		dBjout_dt = (1/(pjout*At*h*Lje))*( pjin*Fjin*Bjin + \
					pjin*Fjin*Bjout*((Cpjin*Tjin - Cpjout*Tjout)/(Hvv - Cpjout*Tjout) - 1) + \
					Bjout*U*A*(Tvin - Tjout)/(Hvv - Cpjout*Tjout))

		# Level differential equation
		dLje_dt = (1/(pjout*At*h))*( pjin*Fjin*(Hvv - Cpjin*Tjin)/(Hvv - Cpjout*Tjout) - \
					pjout*Fjout - U*A*(Tvin - Tjout)/(Hvv - Cpjout*Tjout))

		dx_dt = [dBjout_dt, dLje_dt]
		return dx_dt