# -*- coding: utf-8 -*-
#! python

# Installed Libs
import math

class water_properties:

	def density(self, Temperature):
		'''
		Water density from temperature
		Based on Llano 2005. Paper: Modeling and Simulation of Vertical Continuous Cooling Crystallizers for the Sugar Industry

		Parameters:
		Water temperature [C]

		Result:
		Water density [kg/m3]
		'''
		self.Tmp = Temperature

		q0 = 999.83952 + 16.952577*self.Tmp - 7.9905127e-3*self.Tmp**2;
		q1 = -46.241757e-6*self.Tmp**3 + 105.84601e-9*self.Tmp**4;
		q2 = -281.03006e-12*self.Tmp**5;
		q3 = 1 + 16.887236e-3*self.Tmp;

		pw = (q0 + q1 + q2)/q3;
		return pw

	def enthalpy(self, Temperature):
		'''
		Specific water enthalpy from temperature
		Based on Damour et. al. 2011. Paper: Multivariable linearizing control of an industrial sugar crystallization process 

		Parameters:
		Water temperature [C]

		Result:
		Water Enthalpy [J/kg]
		'''
		self.Tmp = Temperature

		hw = 2323.3 + 4106.7*self.Tmp;
		return hw



class liquor_properties:

	def heat_capacity(self, Temperature, Brix, Purity):
		'''
		Specific heat capacity of mother liquor
		Based on Llano 2005. Paper: Modeling and Simulation of Vertical Continuous Cooling Crystallizers for the Sugar Industry

		Parameters:
		Mother liquor temperature [C]
		Mother liquor brix [kg/kg]
		Mother liquor purity [kg/kg]

		Result:
		Specific heat capacity of liquor [J/(kg.K)]
		'''
		self.Tmp = Temperature
		self.Brx = Brix
		self.Pty = Purity

		Cpl = 4184 - 2971*(self.Brx) + 460*(self.Brx)*(self.Pty) + 7.5*(self.Brx)*(self.Tmp)
		return Cpl

	def sucrose_solution_density(self, Brix, Temperature):
		'''
		Density of a pure sucrose solution. 
		Based on Llano 2005. Paper: Modeling and Simulation of Vertical Continuous Cooling Crystallizers for the Sugar Industry

		Parameters:
		Pure sucrose solution brix [kg/kg]
		Pure sucrose solution temperature [C]

		Result:
		Pure sucrose solution density [kg/m3]
		'''
		self.Brx = Brix
		self.Tmp = Temperature
		self.t = ((self.Tmp) - 20.0)/100.0

		c = [[385.1761,135.3705,40.9299,-3.9643,13.4853,-17.2890],
		[-46.2720,-7.1720, 1.1597, 5.1126,17.5254,0],
		[59.7712,7.2491, 12.3630,-35.4791,0,0],
		[-47.2207,-21.6977,27.6301,0,0,0],
		[18.3184,12.3081,0,0,0,0]]

		p1 = (c[0][0]*(self.Brx)) + (c[0][1]*((self.Brx)**2)) + (c[0][2]*((self.Brx)**3)) + (c[0][3]*((self.Brx)**4)) + (c[0][4]*((self.Brx)**5)) + (c[0][5]*((self.Brx)**6))
		p2 = (self.t)*((c[1][0]*(self.Brx)) + (c[1][1]*((self.Brx)**2)) + (c[1][2]*((self.Brx)**3)) + (c[1][3]*((self.Brx)**4)) + (c[1][4]*((self.Brx)**5)))
		p3 = self.t**2*(c[2][0]*(self.Brx) + c[2][1]*(self.Brx)**2 + c[2][2]*(self.Brx)**3 + c[2][3]*(self.Brx)**4)
		p4 = self.t**3*(c[3][0]*(self.Brx) + c[3][1]*(self.Brx)**2 + c[3][2]*(self.Brx)**3)
		p5 = self.t**4*(c[4][0]*(self.Brx) + c[4][1]*(self.Brx)**2)

		wt_dnst = water_properties().density(self.Tmp)
		ps = wt_dnst + p1 + p2 + p3 + p4 + p5;
		return ps

	def density(self, Temperature, Brix, Purity):
		'''
		Density of mother liquor
		Based on Llano 2005. Paper: Modeling and Simulation of Vertical Continuous Cooling Crystallizers for the Sugar Industry

		Parameters:
		Mother liquor temperature [C]
		Mother liquor brix [kg/kg]
		Mother liquor purity [kg/kg]

		Result:
		Mother liquor density [kg/m3]
		'''
		self.Tmp = Temperature
		self.Brx = Brix
		self.Pty = Purity

		scr_dst = self.sucrose_solution_density(self.Brx, self.Tmp)
		pl = scr_dst - 1 + math.exp((1 - (self.Pty))*(6.927*(self.Brx)**2 + 1.165*(self.Brx)))
		return pl

	def viscosity(self, Temperature, Brix, Purity):
		'''
		Viscosity of mother liquor
		Based on Sugar Technologist Manual - Bubnnik et all 1995

		Parameters:
		Mother liquor brix [kg/kg]
		Mother liquor purity [kg/kg]
		Mother liquor temperature [C]

		Result:
		Mother liquor viscosity [Pa s]
		'''
		self.Tmp = Temperature
		self.Brx = Brix
		self.Pty = Purity

		x2 = (30.0 - self.Tmp)/(91.0 + self.Tmp);
		x3 = 0.85 + 0.15*self.Pty;
		x4 = self.Brx*(0.962 + ((1 - 0.962)*self.Pty));
		x1 = x4/(1900.0 - (18.0*x4));

		vl = 0.001*(10.0**((22.46*x1) - (0.114) + (x2*(1.1 + (43.1*x3*(x1**1.25))))));
		return vl

	def thermal_conductivity(self, Temperature, Brix):
		'''
		Thermal conductivity of liquor
		Based on Llano 2005. Paper: Modeling and Simulation of Vertical Continuous Cooling Crystallizers for the Sugar Industry
		
		Parameters:
		Mother liquor temperature [C]
		Mother liquor brix [kg/kg]
		
		Result:
		Thermal conductivity of liquor [W/(m.K)]
		'''
		self.Tmp = Temperature
		self.Brx = Brix
		tcl = (self.Brx*((5.466*(10**-6)*(self.Tmp**2)) - (1.176*(10**-3)*self.Tmp) - 0.3024)) - (7.847*(10**-6)*(self.Tmp**2))+(1.976*(10**-3)*self.Tmp)+0.563
		return tcl

	def loss_saccharose(self, time, Temperature, Brix, SolIn, Purity, pH):
		'''
		Loss of saccharose in liquor
		Based on Rein 2012, Ingenieria de la caña de ázucar

		Parameters:
		Residence time [min]
		Mother liquor temperature [C]
		Mother liquor brix [kg/kg]
		Mother liquor insoluble solids [kg/kg]
		Mother liquor purity [kg/kg]
		Mother liquor pH []

		Result:
		Percentage of saccharose loss [%]
		'''

		cw = ((100.0 - (Brix*100.0) - (SolIn*100.0))*(self.density(Temperature,Brix,Purity)))*0.001
		k = (8.12831)*(10**16)*cw*math.exp((-13055.7/(Temperature + 273.15)) - 2.30259*pH)

		lss_sac = (1 - math.exp(-k*time))*100.0
	
		return lss_sac

class vapor_properties:

	def temperature(self, Pressure):
		'''
		Saturated vapor temperature from pressure
		
		Parameters:
		Vapor pressure [Pa]

		Result:
		Vapor temperature [C]
		'''
		self.Prs = Pressure
		#----Buck Equation----
		#Pressure=0.61121*math.exp((18.678*(Temperature/234.5))*(Temperature/(257.14+Temperature)));

		#Solving Buck Equation
		#For Temperature>0 C
		Ts = -117.25*(math.sqrt(((math.log(1/(((self.Prs)*0.001))))**2)+40.7576*math.log(1/((self.Prs)*0.001))+328.56)-math.log(1/((self.Prs)*0.001))-18.1857);
		return Ts

	def density(self, Pressure):
		'''
		Saturated vapor density from pressure
		Based on W. M. Haynes 2015, Handbook of Chemistry and physics %Ideal gas formula
		
		Parameters:
		Vapor pressure  [Pa]

		Result:
		Vapor density  [kg/m3]
		'''
		self.Prs = Pressure
		self.Tv = self.temperature(self.Prs) + 273.15; #vapor temperature in Kelvin

		#Rg = 8.31447; % universal gas constant kg m2 s?2 K?1 mol?1
		#m = 18.01528/1000; % molar mass of water kg/mole
		#pv = (Pv*m)/(Rg*Tv) % Ideal gas formula
		
		Pv = (0.0022*self.Prs)/self.Tv
		return Pv

	def enthalpy(self, Temperature, Pressure):
		'''
		Specific vapor enthalpy from temperature
		Based on Damour et. al. 2011. Paper: Multivariable linearizing control of an industrial sugar crystallization process 

		Parameters:
		Temperature [C]
		Pressure [Pa]

		Result:
		Enthalpy [J/kg]
		'''
		self.Tmp = Temperature
		self.Prs = Pressure
		hv = 2499980 - 24186*((self.Prs)*0.00001) + (1891.1+106.1*((self.Prs)*0.00001))*self.Tmp;
		return hv

	def thermal_conductivity(self, Temperature):
		'''
		Thermal conductivity of vapor 
		Based on Felipe Ospina - Cenicana 2017
		
		Parameters:
		Vapor temperature [C]

		Result:
		Thermal conductivity of vapor [W/(m.C)]
		'''
		self.Tmp = Temperature

		Tcv = (-0.0000087935*(self.Tmp**2)) + (0.0020729683*self.Tmp) + (0.5608807885)
		return Tcv

	def viscosity(self, Temperature):
		'''
		Viscosity of vapor
		Based on Felipe Ospina - Cenicana 2017

		Parameters:
		Vapor temperature [C]

		Result:
		Viscosity of vapor [Pa s]
		'''
		self.Tmp = Temperature

		vv = ((0.0000138297*self.Tmp**4) - (0.004990139*self.Tmp**3) + (0.6914895681*self.Tmp**2) - (47.4893894146*self.Tmp) + 1722.8465450434)*10**-6
		return vv