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
		self.sp = sp
		self.pv = pv
		self.error = self.sp - self.pv

		p = self.kp*self.error
		i = self.ki*0.5*(self.error + self.error0)*self.ts
		# i = self.i0+((self.error+self.error0)*0.5)*self.ts
		d = self.kd*(self.error - self.error0)/self.ts

		# self.mv = p + self.ki*i + d
		self.mv = p + i + d
		self.error0 = self.error # saving error for next loop
		self.i0 = i
		return self.mv


		