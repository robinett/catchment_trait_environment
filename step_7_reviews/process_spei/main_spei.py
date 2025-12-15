import os
import datetime
from process_spei import process_spei


def main():
    # --------- EDIT THESE PATHS ---------
    base_dir = (
        "/discover/nobackup/trobinet/from_aws/"
        "pso/step_7_reviews"
    )
    spei_nc_path = (
        "/discover/nobackup/trobinet/from_aws/pso/"
        "step_7_reviews/data/spei_raw/nclimgrid-spei-gamma-06.nc"
    )
    output_dir = os.path.join(base_dir, "outputs")
    os.makedirs(output_dir, exist_ok=True)

    # Tiles that MUST appear as columns in the output:
    intersecting_tiles_fname = (
        "/discover/nobackup/trobinet/from_aws/pso/"
        "step_6x_extrapolate/outputs/"
        "extrapolation_pix.csv"
    )

    # Catchment geometry for each tile:
    catch_tile_info_fname = (
        "/discover/nobackup/trobinet/from_aws/pso/"
        "step_1x_choose_tiles_large/data/tile_coord.csv"
    )

    # Monthly range for SPEI(06)
    start_date = datetime.date(1992, 1, 1)
    end_date = datetime.date(2015, 12, 31)

    # Weights cache
    weights_fname = os.path.join(
        output_dir, "spei_weights_selected_tiles.pkl"
    )

    # Output file
    truth_fname_selected_tiles = os.path.join(
        output_dir,
        "spei06_truth_nclimgrid_gamma_{start}_{end}_extrap.csv",
    )

    # Toggles
    load_weights = False
    save_weights = True
    save_selected_truth = True
    # ------------------------------------

    pro = process_spei(base_dir)

    # Only these tiles will be enforced as columns in outputs
    pro.get_tiles(intersecting_tiles_fname)

    # Build/Load weights ONLY for those tiles; tiles with no overlap
    # will still be present with empty intersections.
    pro.get_weights(
        catch_tile_info_fname=catch_tile_info_fname,
        load_weights=load_weights,
        save_weights=save_weights,
        weights_fname=weights_fname,
        example_spei_fname=spei_nc_path,
    )

    # Build monthly truth; every requested tile gets a column
    pro.create_truth_spei(
        start=start_date,
        end=end_date,
        spei_nc_path=spei_nc_path,
        save_selected_truth=save_selected_truth,
        truth_fname_selected_tiles=truth_fname_selected_tiles,
    )


if __name__ == "__main__":
    main()

