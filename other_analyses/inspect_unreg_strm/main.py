import os
from inspect_stream import inspect

def main():
    # where are we now
    base_dir = '/shared/pso/other_analyses/inspect_unreg_strm'
    # where is the data located that we want to process
    data_dir = os.path.join(
        base_dir,
        'data'
    )
    # where is the station info located
    station_info_fname = os.path.join(
        data_dir,
        'Riverflow_Station_Information_NA-M36.nc'
    )

    # get an instance
    insp = inspect()
    insp.inspect_station_info(station_info_fname)

if __name__ == '__main__':
    main()
