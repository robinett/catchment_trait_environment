import netCDF4 as nc
import numpy as np
import sys

class inspect:
    def inspect_station_info(self,fname):
        s_info = nc.Dataset(fname)
        print('dataset information:')
        print(s_info)
        print(s_info.variables)
        # get the files that we are interested in
        var_names = ['Station_ID','N_cells_basin','SMAPID','SMAPFrac']
        for this_var_name in var_names:
            if this_var_name == 'Station_ID':
                this_var = np.array(s_info[this_var_name]).astype(np.float64)
                for stat in this_var:
                    stat_str = ''
                    for num in stat:
                        stat_str = stat_str + str(int(num))
                    print(stat_str)
            else:
                this_var = np.array(s_info[this_var_name])
            print(this_var_name)
            print(this_var)
