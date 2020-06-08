# -*- coding: utf-8 -*-
"""
Simple heat pump model with CoolProp (high-level interface). 

Only work with subcritical cycles and single substances

Ther thermodynamic properties are called via CoolProp and the property plots are generated with the dedicated CoolProp functions

Cycles conditions are imposed and the secondary fluid temperature profiles are recalculated.

@author: Sylvain Quoilin


"""

# clear variables:
from IPython import get_ipython
get_ipython().magic('reset -sf')

"Imports"
import CoolProp.CoolProp as CP
import numpy as np

from misc.utils import NoStdStreams
from components.heat_exchangers import hx
        

p = np.zeros(8, dtype=float)
v = np.zeros(8, dtype=float)
T = np.zeros(8, dtype=float)
x = np.zeros(8, dtype=float)
h = np.zeros(8, dtype=float)
s = np.zeros(8, dtype=float)

T_hw = np.zeros(8, dtype=float)     # hot water
T_cw = np.zeros(8, dtype=float)     # cold water

# Cycle parameters:
fluid ='R245fa'
T_ev = 273.15 + 80
T_cd = 273.15 + 145

DELTAT_sh = 5
DELTAT_sc = 5
epsilon_s = 0.6

# Heat sink parameters:
T_w_su_cd =273.15 + 120        # cooling water inlet temperature /K
DELTAT_cd=10                    # pinch point temperature of cooling water /K
Q_dot_cd=10000                  # heat capacity flowrates of cooling water kW/K

# Heat source parameters:
T_w_su_ev =273.15 + 95        # cooling water inlet temperature /K
DELTAT_ev=10                    # pinch point temperature of cooling water /K

# get fluid properties:
#T_crit = CP.PropsSI("Tcrit",fluid)
#p_crit = CP.PropsSI("Pcrit",fluid)
p_low = CP.PropsSI('P','Q', 0.5, 'T', T_ev, fluid)
p_high = CP.PropsSI('P','Q', 0.5, 'T', T_cd, fluid)
#h_crit = CP.PropsSI('H','T', T_crit, 'P', p_crit, fluid)
#s_crit = CP.PropsSI('S','T', T_crit, 'P', p_crit, fluid)

#Evaporator outlet:
p[0] = p_low
s[0] = CP.PropsSI('S','Q', 1, 'P', p_low, fluid)
T[0] = CP.PropsSI('T','Q', 1, 'P', p_low, fluid)
h[0] = CP.PropsSI('H','Q', 1, 'P', p_low, fluid)

# Compressor inlet:
p_su_cp = p_low
T_su_cp = CP.PropsSI('T','Q', 1, 'P', p_low, fluid) + DELTAT_sh
h_su_cp = CP.PropsSI('H','T', T_su_cp, 'P', p_su_cp, fluid)
s_su_cp = CP.PropsSI('S','T', T_su_cp, 'P', p_su_cp, fluid)
s[1] = s_su_cp
T[1] = T_su_cp
h[1] = h_su_cp
p[1] = p_su_cp

#Compressor outlet:
p_ex_cp = p_high
h_ex_cp_s = CP.PropsSI('H','S', s_su_cp, 'P', p_ex_cp, fluid)
h_ex_cp = h_su_cp - (h_su_cp - h_ex_cp_s)/epsilon_s
T_ex_cp = CP.PropsSI('T','H', h_ex_cp, 'P', p_ex_cp, fluid)
s_ex_cp = CP.PropsSI('S','H', h_ex_cp, 'P', p_ex_cp, fluid)
p[2] = p_high
s[2] = s_ex_cp
T[2] = T_ex_cp
h[2] = h_ex_cp

#Saturated vapor in the condenser:
p[3] = p_high
s[3] = CP.PropsSI('S','Q', 1, 'P', p_high, fluid)
T[3] = CP.PropsSI('T','Q', 1, 'P', p_high, fluid)
h[3] = CP.PropsSI('H','Q', 1, 'P', p_high, fluid)

#Saturated liquid in the condenser:
p[4] = p_high
s[4] = CP.PropsSI('S','Q', 0, 'P', p_high, fluid)
T[4] = CP.PropsSI('T','Q', 0, 'P', p_high, fluid)
h[4] = CP.PropsSI('H','Q', 0, 'P', p_high, fluid) 
  
T[5] = T[4] - DELTAT_sc
p[5] = p_high
h[5] = CP.PropsSI('H','T', T[5], 'P', p[5], fluid) 
s[5] = CP.PropsSI('S','T', T[5], 'P', p[5], fluid)
  
#Inlet of the evaporator:
h[6] = h[5]
p[6] = p_low
T[6] = CP.PropsSI('T','H', h[5], 'P', p[6], fluid)
s[6] = CP.PropsSI('S','H', h[5], 'P', p[6], fluid)
    
h[7] = h[0]
p[7] = p[0]
T[7] = T[0]
s[7] = s[0]

print("The temperature of each state")
print(T)
print("The pressure of each state (Pa) ")
print(p)
print("The enthalpy of each state (J/kg)")
print(h)
print("The entropy of each state (J/kg/K)")
print(s)

# heat sink:
cp_w = 4800
fluid2 ='water'
M_dot = Q_dot_cd / (h[2] - h[5])
T_hw[5] = T_w_su_cd
T_hw[3] = T[3] - DELTAT_cd
M_dot_hw = M_dot * (h[3] - h[5])/(cp_w*(T_hw[3] - T_hw[5]))
T_hw[4] = T_hw[5] + M_dot * (h[4] - h[5]) / (cp_w * M_dot_hw)
T_hw[2] = T_hw[5] + M_dot * (h[2] - h[5]) / (cp_w * M_dot_hw)

print ("The mass flow rate of working fluid is %.2f kg/s" %(M_dot))
print("The mass flowrate of hot water is %.2f kg/s" %(M_dot_hw))

# Heat source:
T_cw[1] = T_w_su_ev
T_cw[6] = T[6] + DELTAT_ev
M_dot_cw = M_dot * (h[1] - h[6])/(cp_w*(T_cw[1] - T_cw[6]))
T_cw[0] = T_cw[6] + M_dot * (h[0] - h[6]) / (cp_w * M_dot_cw)

# Temperature profile in the condenser:
s_cd,T_cd,T_hf,pinch_cd = hx(fluid,10,M_dot,h[2],h[5],p_high,p_high,T_w_su_cd,M_dot_hw*cp_w)

# Temperature profile in the evaporator:
s_ev,T_ev,T_cf,pinch_ev = hx(fluid,10,M_dot,h[6],h[1],p_low,p_low,T_w_su_ev,M_dot_cw*cp_w)

#%%

from CoolProp.Plots import PropertyPlot
from CoolProp.Plots.SimpleCycles import StateContainer
import pickle
import os

cache_plot = False
if cache_plot:
    filename = fluid + '.p'
    if os.path.isfile(filename):   # Load previously saved plot
        pp = pickle.load(open( filename, "rb" ) )
    else:
        states = StateContainer()
        states_hf = StateContainer()
        states_cf = StateContainer()
        pp = PropertyPlot('HEOS::'+fluid, 'TS')
        with NoStdStreams():
            pp.calc_isolines()
        with open(filename, 'wb') as f: 
            pickle.dump(pp, f) 
else:
    states = StateContainer()
    states_hf = StateContainer()
    states_cf = StateContainer()
    pp = PropertyPlot('HEOS::'+fluid, 'TS')
    with NoStdStreams():
        pp.calc_isolines()
        
    
for i in range(3):
    states[i,'T'] = T[i]
    states[i,"S"] = s[i]
    
for i,Tx in enumerate(T_cd):
    states.append({'T':Tx,'S':s_cd[i]})

for i in range(4,len(T)):
    states.append({'T':T[i],'S':s[i]})  
states.append({'T':T[1],'S':s[1]})    # for some reasons, the second point needs to be repeated to close the cycle
    
for i,Tx in enumerate(T_hf):
    states_hf.append({'T':Tx,'S':s_cd[i]})
    
for i,Tx in enumerate(T_cf):
    states_cf.append({'T':Tx,'S':s_ev[i]})

with NoStdStreams():
    pp.draw_process(states,line_opts={'color':'green'})
    pp.draw_process(states_hf,line_opts={'color':'red', 'linestyle':'dashed'})
    pp.draw_process(states_cf,line_opts={'color':'blue', 'linestyle':'dashed'})

pp.show()




