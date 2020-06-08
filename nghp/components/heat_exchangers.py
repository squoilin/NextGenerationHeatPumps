# -*- coding: utf-8 -*-
"""
Heat exchanger models

@author: Sylvain Quoilin
"""
import CoolProp.CoolProp as CP
import numpy as np
from scipy.optimize import fsolve


def hx(fluid,N,M_dot,h_su,h_ex,p_su,p_ex,T_sf_su,C_dot_sf):
    '''
    Simple discretized heat exchanger model 
    Uses the CoolProp high-level interface
    
    Parameters
    ----------
    fluid : String
        Primary working fluid.
    N : Int
        Number of discretization elements.
    M_dot : Float
        Primary working fluid mass flow rate.
    h_su : Float
        Inlet enthalpy of the primary working fluid.
    h_ex : Float
        Outlet enthalpy of the primary working fluid.
    h_su : Float
        Inlet pressure of the primary working fluid.
    h_ex : Float
        Outlet pressure of the primary working fluid.
    T_sf_su : Float
        Inlet temperature of the secondary working fluid.
    C_dot_sf : Float
        Heat capacity flow rate for the secondary fluid.

    Returns
    -------
    Temperature profiles, pinch point temperature difference.

    '''
    T_sf = np.zeros(N)
    Ti = np.zeros(N)
    si = np.zeros(N)
    T_sf_ex = T_sf_su + M_dot/C_dot_sf * (h_su - h_ex)
    
    for i in range(N):
         hi= h_su + (h_ex-h_su) / (N - 1) * i 
         Pi = p_su + (p_su-p_ex) / (N - 1) * i 
         Ti[i] = CP.PropsSI('T','P',Pi,'H',hi,fluid) 
         si[i] = CP.PropsSI('S','P',Pi,'H',hi,fluid) 
         T_sf[i] = T_sf_ex + (T_sf_su - T_sf_ex) / (N - 1) * i 

    if T_sf[N-1] > Ti[0]:
        pinch = np.min(T_sf - Ti)
    else:
        pinch = np.min(Ti - T_sf)
    print("The pinch temperature difference is %.2f K" % (pinch))
    return si,Ti,T_sf,pinch




def hx_pinch(library,fluid,N,M_dot,CPstate_su,CPstate_ex,T_sf_su,pinch,T_sf_ex=None):
    '''
    Descretized heat exchanger model with imposed pinch point temperature 
    difference or with imposed secondary fluid outlet temperature
    Uses the CoolProp low-level interface
    
    Parameters
    ----------
    fluid : String
        Primary working fluid.
    N : Int
        Number of discretization elements.
    M_dot : Float
        Primary working fluid mass flow rate.
    CPstate_su : CoolProp.CoolProp.AbstractState
        Inlet state of the primary working fluid.
    CPstate_ex : CoolProp.CoolProp.AbstractState
        Outlet stete of the primary working fluid.
    T_sf_su : Float
        Inlet temperature of the secondary working fluid.
    pinch : Fload
        Imposed pinch point temperature difference
    T_sf_ex : Float
        Secondary fluid outlet temperature (if specified, the pinch point is not imposed anymore)

    Returns
    -------
    Temperature profiles

    '''
    T_sf = np.zeros(N)
    Ti = np.zeros(N)
    si = np.zeros(N)

    h_su = CPstate_su.keyed_output(CP.iHmass)
    h_ex = CPstate_ex.keyed_output(CP.iHmass)
    p_su = CPstate_su.keyed_output(CP.iP)
    p_ex = CPstate_ex.keyed_output(CP.iP)

    CPstates = []
    for i in range(N):
        CPstates.append(CP.AbstractState(library, fluid))
        if  '&' in fluid:
            CPstates[i].set_mole_fractions([x_molar, 1 - x_molar])
        hi= h_su + (h_ex-h_su) / (N - 1) * i 
        Pi = p_su + (p_su-p_ex) / (N - 1) * i 
        CPstates[i].update(CP.HmassP_INPUTS,hi,Pi)
        Ti[i] = CPstates[i].keyed_output(CP.iT)
        si[i] = CPstates[i].keyed_output(CP.iSmass)
    
    def pinch_calc(C_dot):
        T_sf_ex = T_sf_su + M_dot/C_dot[0] * (h_su - h_ex)
        T_sf = np.linspace(T_sf_ex,T_sf_su,N)
        if T_sf[N-1] > Ti[0]:
            DELTAT = np.min(T_sf - Ti)
        else:
            DELTAT = np.min(Ti - T_sf)
        return DELTAT - pinch
    if h_su > h_ex and Ti[-1] - T_sf_su < pinch or h_su < h_ex and T_sf_su - Ti[-1] < pinch:
        print('The pinch point requirement cannot be respected')
        sys.exit(1)
    if T_sf_ex is not None:
        C_dot_sf = (T_sf_ex - T_sf_su)/(M_dot * (h_su - h_ex))
    else:
        C_dot_sf = fsolve(pinch_calc,10)[0]
        T_sf_ex = T_sf_su + M_dot/C_dot_sf * (h_su - h_ex)
    T_sf = np.linspace(T_sf_ex,T_sf_su,N)    
    if T_sf[N-1] > Ti[0]:
        DELTAT = np.min(T_sf - Ti)
    else:
        DELTAT = np.min(Ti - T_sf)        
    print("The pinch temperature difference is %.2f K" % (DELTAT))
    return si,Ti,T_sf,C_dot_sf,CPstates