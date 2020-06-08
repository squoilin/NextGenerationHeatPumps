# -*- coding: utf-8 -*-
"""

Heat pump model with CoolProp (low-level interface). 

Suitable both for subcritical or supercritical cycles, with single substances or zeotropic mixtures.

The thermodynamic properties are called via CoolProp and the property plots are generated with the dedicated CoolProp functions

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
from components.heat_exchangers import hx_pinch

p = np.zeros(8, dtype=float)
v = np.zeros(8, dtype=float)
T = np.zeros(8, dtype=float)
x = np.zeros(8, dtype=float)
h = np.zeros(8, dtype=float)
s = np.zeros(8, dtype=float)

# Cycle parameters:
fluid ='R245fa'
fluid = 'Butane'
#fluid = 'HEOS::Pentane[0.697615]&Hexane[0.302385]'
fluid1 = 'Butane'
fluid2 = 'Hexane'
x_molar = 0.5
fluid = fluid1 + '&' + fluid2
fluid = 'R1336MZZ(Z)'
#fluid = 'R1234ze(E)'

library = 'HEOS'


fluid = 'R1336MZZ(Z)'
library = 'REFPROP'

DELTAT_sh = 7
DELTAT_sc = 32
epsilon_s = 0.7

# Heat sink parameters:
T_w_su_cd =273.15 + 130        # cooling water inlet temperature /K
DELTAT_cd=7.5                    # pinch point temperature of cooling water /K
Q_dot_cd=10000                  # heat capacity flowrates of cooling water kW/K

# Heat source parameters:
T_w_su_ev =273.15+135        # cooling water inlet temperature /K
DELTAT_ev=7.5                    # pinch point temperature of cooling water /K
T_w_ex_ev = 0    # set to 0 if the cf T profile is calculated with the pinch

limits = [1.25,2,250,550]   # limits for plotting

cycle  = []
for i in range(8):
    cycle.append(CP.AbstractState(library, fluid))
    if '&' in fluid:
        cycle[i].set_mole_fractions([x_molar, 1 - x_molar])
#    c.build_phase_envelope("dummy")

T_crit = cycle[0].T_critical()
p_crit = cycle[0].p_critical()

impose_T = False  # Set to true to impose the saturation temperatures and not the pressure

if impose_T:
    T_ev = 273.15 + 100
    T_cd = 273.15 + 160
    # mid-evaporation and mid-condensation:
    ev_mid = CP.AbstractState(library, fluid)
    if '&' in fluid:
        ev_mid.set_mole_fractions([x_molar, 1 - x_molar])
    ev_mid.build_phase_envelope("dummy")
    ev_mid.update(CP.QT_INPUTS,0.5,T_ev)
    
    cd_mid = CP.AbstractState(library, fluid)
    if '&' in fluid:
        cd_mid.set_mole_fractions([x_molar, 1 - x_molar])
    cd_mid.build_phase_envelope("dummy")
    cd_mid.update(CP.QT_INPUTS,0.5,T_cd)

    # get fluid properties:
    p_low = ev_mid.keyed_output(CP.iP)
    p_high = cd_mid.keyed_output(CP.iP)
else:
    p_low = 10.5E5
    p_high = 32E5

#Evaporator outlet:
p[0] = p_low
cycle[0].update(CP.PQ_INPUTS,p[0],1)
s[0] = cycle[0].keyed_output(CP.iSmass)
T[0] = cycle[0].keyed_output(CP.iT)
h[0] = cycle[0].keyed_output(CP.iHmass)

# Compressor inlet:
p[1] = p_low
T[1] = T[0] + DELTAT_sh
cycle[1].update(CP.PT_INPUTS,p[1],T[1])
s[1] = cycle[1].keyed_output(CP.iSmass)
h[1] = cycle[1].keyed_output(CP.iHmass)

#Compressor outlet:
cp_ex = CP.AbstractState(library, fluid)
if '&' in fluid:
    cp_ex.set_mole_fractions([x_molar, 1 - x_molar])
cp_ex.build_phase_envelope("dummy")
cp_ex.specify_phase(CP.iphase_gas)
cp_ex.update(CP.PSmass_INPUTS,p_high,s[1])
h_ex_cp_s = cp_ex.keyed_output(CP.iHmass)

p[2] = p_high
h[2] = h[1] - (h[1] - h_ex_cp_s)/epsilon_s
cycle[2].update(CP.HmassP_INPUTS,h[2],p[2])
s[2] = cycle[2].keyed_output(CP.iSmass)
T[2] = cycle[2].keyed_output(CP.iT)

if p_high < p_crit:
    #Saturated vapor in the condenser:
    p[3] = p_high
    cycle[3].update(CP.PQ_INPUTS,p[3],1)
    s[3] = cycle[3].keyed_output(CP.iSmass)
    T[3] = cycle[3].keyed_output(CP.iT)
    h[3] = cycle[3].keyed_output(CP.iHmass)
    
    #Saturated liquid in the condenser:
    p[4] = p_high
    cycle[4].update(CP.PQ_INPUTS,p[4],0)
    s[4] = cycle[4].keyed_output(CP.iSmass)
    T[4] = cycle[4].keyed_output(CP.iT)
    h[4] = cycle[4].keyed_output(CP.iHmass)

else:
    p[3],s[3],T[3],h[3] = p[2],s[2],T[2],h[2]
    p[4],s[4],T[4],h[4] = p[2],s[2],T[2],h[2]
  
p[5] = p_high
if p_high < p_crit:
    T[5] = T[4] - DELTAT_sc
else:
    T[5] = T_crit - DELTAT_sc
cycle[5].update(CP.PT_INPUTS,p[5],T[5])
s[5] = cycle[5].keyed_output(CP.iSmass)
h[5] = cycle[5].keyed_output(CP.iHmass)
  
#Inlet of the evaporator:
h[6] = h[5]
p[6] = p_low
cycle[6].update(CP.HmassP_INPUTS,h[6],p[6])
s[6] = cycle[6].keyed_output(CP.iSmass)
T[6] = cycle[6].keyed_output(CP.iT)
    
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

# Mass flow rate:
M_dot = Q_dot_cd / (h[2] - h[5])

# COP:
COP = (h[2] - h[5])/(h[2] - h[1])

print("The mass flow rate of working fluid is %.3f kg/s" %(M_dot))
print("The COP is %.2f" %(COP))

# Temperature profile in the condenser:
s_cd,T_cd,T_hf,C_dot_cd,CPstates_cd = hx_pinch(library,fluid,40,M_dot,cycle[2],cycle[5],T_w_su_cd,DELTAT_cd)

# Temperature profile in the evaporator:
if T_w_ex_ev > 0:
    s_ev,T_ev,T_cf,pinch_ev,CPstates_ev = hx_pinch(library,fluid,20,M_dot,cycle[6],cycle[1],T_w_su_ev,DELTAT_ev,T_sf_ex=T_w_ex_ev)
else:
    s_ev,T_ev,T_cf,pinch_ev,CPstates_ev = hx_pinch(library,fluid,20,M_dot,cycle[6],cycle[1],T_w_su_ev,DELTAT_ev)
print("Heat sink: %.1f°C to %.1f°C" %(T_hf[-1]-273.15, T_hf[0] - 273.15))
print("Heat source: %.1f°C to %.1f°C" %(T_cf[-1]-273.15, T_cf[0] - 273.15))


#%%

from CoolProp.Plots import PropertyPlot
from CoolProp.Plots.SimpleCycles import StateContainer

# for the property plot, the composition needs to be embedded in the fluid name:
if '&' in fluid:
    fluid_full = library+'::'+fluid1 + '[' + str(x_molar) + ']&' + fluid2 + '[' + str(1-x_molar) + ']'
else:
    fluid_full =  library+'::'+fluid

states = StateContainer()
states_hf = StateContainer()
states_cf = StateContainer()
pp = PropertyPlot(fluid_full, 'TS')
if len(limits)>0:
    pp.set_axis_limits(limits)
with NoStdStreams():
    pass
    pp.calc_isolines()
        
    
for i in range(3):
    states[i,'T'] = T[i]
    states[i,"S"] = s[i]
    
for i,Tx in enumerate(T_cd):
    states.append({'T':Tx,'S':s_cd[i]})

for i in range(5,len(T)):
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

def draw_ts():
    '''
    simplified isolines for mixtures in which the calc_isolines() function does not work

    Returns
    -------
    None.

    '''
    states_sat = StateContainer()
    N = 100
    sat_l = CP.AbstractState(library, fluid)
    sat_v = CP.AbstractState(library, fluid)
    Tsat = np.linspace(limits[2],T_crit-0.00001,N)
    s_l = np.zeros(N)
    s_v = np.zeros(N)
    
    if '&' in fluid:
        sat_l.set_mole_fractions([x_molar, 1 - x_molar])
        sat_v.set_mole_fractions([x_molar, 1 - x_molar])
    for i,T in enumerate(Tsat):
        sat_l.update(CP.QT_INPUTS,0,T)
        s_l[i] = sat_l.keyed_output(CP.iSmass)
        sat_v.update(CP.QT_INPUTS,1,T)
        s_v[i] = sat_v.keyed_output(CP.iSmass)
    for i in range(N):
        states_sat.append({'T':Tsat[i],'S':s_l[i]})
    for i in range(N):
        states_sat.append({'T':Tsat[N-1-i],'S':s_v[N-1-i]})
    
    with NoStdStreams():
        pp.draw_process(states_sat,line_opts={'color':'black','linewidth':0.3})
    
    for x in range(1,10):
        states_x = StateContainer()
        for i in range(3,N):
             s = x/10 * s_l[i] + (1-x/10) * s_v[i]
             states_x.append({'T':Tsat[i],'S':s})
        with NoStdStreams():
            pp.draw_process(states_x,line_opts={'color':'black','linewidth':0.1})

#draw_ts()

pp.figure.set_figheight(8)
pp.figure.set_figwidth(10)
x1,x2,y1,y2 = pp.get_axis_limits()
pp.figure.text(0.2,0.85,"COP = %.2f" %(COP),{'size':16})
pp.title(fluid)

pp.show()




