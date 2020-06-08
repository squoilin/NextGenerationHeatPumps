# -*- coding: utf-8 -*-
"""
Simple script to generate T-s property plots

@author: Sylvain Quoilin
"""

# clear variables:
from IPython import get_ipython
get_ipython().magic('reset -sf')

"Imports"
import CoolProp.CoolProp as CP
from CoolProp.Plots import PropertyPlot
from CoolProp.Plots.SimpleCycles import StateContainer
import pickle
from misc.utils import NoStdStreams



# Cycle parameters:
fluid ='R245fa'
#fluid = 'HEOS::Pentane[0.697615]&Hexane[0.302385]'
fluid1 = 'Pentane'
fluid2 = 'Hexane'
x_molar = 0.5
fluid = fluid1 + '&' + fluid2
#fluid = 'Pentane[0.697615]&Hexane[0.302385]'
#fluid = 'Butane[0.5]&Hexane[0.5]'
fluid = 'Butane'

library = 'HEOS'
library = 'REFPROP'


limits = [0,1.75,250,550]   # limits for plotting

cache_plot = False
if cache_plot:
    filename = fluid + '.p'
    if os.path.isfile(filename):   # Load previously saved plot
        pp = pickle.load(open( filename, "rb" ) )
    else:
        states = StateContainer()
        states_hf = StateContainer()
        states_cf = StateContainer()
        pp = PropertyPlot(library+'::'+fluid, 'TS')
        with NoStdStreams():
            pp.calc_isolines()
        with open(filename, 'wb') as f: 
            pickle.dump(pp, f) 
else:
    states = StateContainer()
    states_hf = StateContainer()
    states_cf = StateContainer()
    pp = PropertyPlot(library+'::'+fluid, 'TS')
#    pp.set_axis_limits(limits)
    with NoStdStreams():
        pp.calc_isolines()
        
pp.figure.set_figheight(8)
pp.figure.set_figwidth(10)
x1,x2,y1,y2 = pp.get_axis_limits()
pp.title(fluid)

pp.show()


