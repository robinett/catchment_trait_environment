import numpy as np
import sys

class harmonize:
    def __init__(self):
        print('bb_g1 should be unitless')
        print('temp need to be in K')
        print('rel humidity should be fractional')
        print('will return med g1 in kPa^0.5')
    def get_med_g1_from_bb_g1(self,bb_g1,temp,rel_hum):
        e_star = self.calculate_e_star(temp)
        vpd = self.calculate_vpd(e_star,rel_hum)
        med_g1 = self.harmonize_to_bb_g1(bb_g1,rel_hum,vpd)
        return med_g1
    def calculate_e_star(self,temp):
        e_0 = 0.6113 #kPa
        Lv = 2.5e6 #J/kg
        Rv = 461 #J/K/kg
        T0 = 273.15 #K
        e_star = e_0*np.exp((Lv/Rv)*((1/T0) - (1/temp)))
        return e_star
    def calculate_vpd(self,e_star,rel_hum):
        e_act = e_star*rel_hum
        vpd = e_star - e_act
        return vpd
    def harmonize_to_bb_g1(self,g1_bb,rel_hum,vpd):
        g1_med = (np.sqrt(vpd))*((g1_bb - rel_hum)/1.6 - 1)
        return g1_med

