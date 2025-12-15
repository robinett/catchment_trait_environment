import pandas as pd
import sys
import numpy as np

class process_streamflow:
    def __init__(self):
        pass
    def get_tiles(self,tile_fname):
        # get the tiles that we are going to run as defined in step 1
        tiles = pd.read_csv(tile_fname,header=None)
        # turn this into an np array
        tiles = np.array(tiles).astype(int)
        # make a nice np array
        tiles = tiles.T
        tiles = tiles[0]
        # save to self
        self.tiles = tiles
        # index for this tile would start at 0, so subtract 1 to get index
        tiles_idx = tiles - 1
        self.tiles_idx = tiles_idx
        return tiles
    def get_hucs(self,huc_dir):
        hucs = pd.read_csv(huc_dir)
        # follow same strategy as for formatting tiles above
        hucs = np.array(hucs).astype(int)
        hucs = hucs.T
        hucs = hucs[0]
        return hucs
    def get_streamflow(self,strm_data_dir):
        # dictionary to hold the streamflow data while we manually process it
        data_dict = {}
        # open file to manually process because formatted weird
        with open(strm_data_dir) as f:
            # iterate through all lines
            for l,line in enumerate(f):
                # splite the lines by tab
                line_spl = line.split()
                # line zeros lists the dates, while the others list streamflow
                # slightly different data types hence the if statements
                # for line 0 holding dates
                if l == 0:
                    data_dict[line_spl[0]] = [
                        int(l) for l in line_spl[1:]
                    ]
                # for all other lines holding streamflow information
                else:
                    data_dict[int(line_spl[0])] = [
                        float(l) for l in line_spl[1:]
                    ]
        # turn the temporary dictionary into a pd dataframe for later use
        strm_data = pd.DataFrame(data_dict)
        # let's rename this column as date for simplicity
        strm_data = strm_data.rename(columns={'huc_cd':'time'})
        # set the index as the data
        strm_data = strm_data.set_index('time')
        print(strm_data)
        return strm_data
    def trim_data(self,data,start,end,hucs):
        # format the start and end dates
        start_fmt = int(start.strftime('%Y%m'))
        end_fmt = int(end.strftime('%Y%m'))
        # select data within the start and end dates
        data_start = data[data.index >= start_fmt]
        data_start_end = data_start[data_start.index <= end_fmt]
        # select only for subset hucs and return
        data_start_end_hucs = data_start_end[hucs]
        return data_start_end_hucs
    def save_data(self,data,fname,start,end):
        start_str = start.strftime('%Y%m')
        end_str = end.strftime('%Y%m')
        final_fname = fname.format(
            start=start_str,
            end=end_str
        )
        data.to_csv(final_fname)
