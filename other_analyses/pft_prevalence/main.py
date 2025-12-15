import os
from get_pft_distribution import pft_dist

def main():
    # where is this script located?
    main_dir = '/shared/pso/other_analyses/pft_prevalence'
    # where is the data located?
    data_dir = os.path.join(main_dir,'data')
    # where should the plots be placed
    plots_dir = os.path.join(main_dir,'plots')
    # where is the pft data located
    pft_data_fname = (
        '/shared/pso/step_5_analyze_outputs/outputs/pft_distribution.csv'
    )
    # where is tile to watershed info located?
    intersection_info_fname = (
        '/shared/pso/step_1_choose_tiles/outputs/intersection_info.pkl'
    )

    # load the pft data
    p = pft_dist()
    pft_data = p.get_info(pft_data_fname)
    # load intersection info
    intersection_info = p.get_intersection_info(
        intersection_info_fname
    )
    # convert to the pfts that we want
    p.get_pso_pfts(pft_data,intersection_info)



if __name__ == '__main__':
    main()
