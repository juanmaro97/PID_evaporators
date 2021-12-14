#-*- coding: utf-8 -*-
#! python

# Installed Libs
import numpy as np

# Pvin[Pa]
# A [m2]
# mjont [kg]


def IDEE_CEM_calc(BxBx,A,Pvin,mjout,tf):
    
    final_brix = BxBx[len(BxBx)-1]*100
    init_brix = BxBx[0]*100
    Pvin = Pvin*0.0001450377 # Pa -> psi
    tf = tf/3600 # s -> h
    mjout = mjout/(1000*tf) # ton/h
    
    print("A: ",A)
    print("init brix: ",init_brix)
    print("final brix: ",final_brix)
    print("Pvin: ",Pvin)
    print("tf: ",tf)
    print("Mjout: ",mjout)

    IDEE = ((mjout/tf)*(final_brix-init_brix))/(A*Pvin)
    CEM = ((mjout/tf)*(final_brix-init_brix)/final_brix)/(A*Pvin)

    return IDEE, CEM
