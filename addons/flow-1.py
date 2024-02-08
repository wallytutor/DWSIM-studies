# -*- coding: utf-8 -*-
from DWSIM.Thermodynamics import *
from DWSIM import *

mass_flow_in  = ims1
#mass_flow_out = oms1

# XXX: get Nm³/h and output kg/s.
mdot = mass_flow_in.GetProp("totalFlow" , "Overall", None, "", "mass")

# mass_flow_out.Clear()
#mass_flow_out.SetProp("totalFlow", "Overall", None, "", "mass", mdot)
#Flowsheet.WriteMessage("")