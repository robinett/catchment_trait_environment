import numpy as np
import pandas as pd
import netCDF4 as nc

class confirm:
    def get_tiles(self,run_tiles_fname,point_force_example_fname):
        # get the tiles that we are going to run as defined in step 1
        run_tiles = pd.read_csv(run_tiles_fname,header=None)
        # turn this into an np array
        run_tiles = np.array(run_tiles).astype(int)
        # make a nice np array
        run_tiles = run_tiles.T
        run_tiles = run_tiles[0]
        # get from the forcing data
        force_data = nc.Dataset(point_force_example_fname)
        force_tiles = np.array(force_data['tile'])
        force_tiles = force_tiles[0]
        return [run_tiles,force_tiles]
    def check_tiles(self,run_tiles,force_tiles):
        num_missing = 0
        for t,ti in enumerate(run_tiles):
            if ti not in force_tiles:
                num_missing += 1
                print(
                    'catchment tile {} is missing from forcing'.format(
                        ti
                    )
                )
        print(
            'the total number of tiles missing from the' +
            'forcing dataset is: {}'.format(
                num_missing
            )
        )
