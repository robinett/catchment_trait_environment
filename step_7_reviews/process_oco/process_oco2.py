import sys
import netCDF4 as nc
import os
import pandas as pd
import numpy as np
import datetime
import copy
import pickle

class process_oco2:
    def __init__(self, base_dir):
        self.base_dir = base_dir

    def get_tiles(self, tiles_fname):
        tiles = pd.read_csv(tiles_fname, header=None)
        tiles = np.array(tiles).astype(int).T[0]
        self.tiles = tiles
        self.tiles_idx = tiles - 1

    def _latlon_from_example(self, example_oco2_fname):
        # open one year to get centers and bounds
        ds = nc.Dataset(example_oco2_fname)
        lon = np.array(ds['lon'])
        lat = np.array(ds['lat'])
        # spacing and half-spacing
        dlon = np.abs(lon[1] - lon[0])
        dlat = np.abs(lat[1] - lat[0])
        hdlon = dlon / 2.0
        hdlat = dlat / 2.0
        lon_left = lon - hdlon
        lon_right = lon + hdlon
        lat_lower = lat - hdlat
        lat_upper = lat + hdlat
        return {
            'center_lon': lon,
            'center_lat': lat,
            'lon_left': lon_left,
            'lon_right': lon_right,
            'lat_lower': lat_lower,
            'lat_upper': lat_upper,
        }

    def get_weights(
        self,
        catch_tile_info_fname,
        load_weights,
        save_weights,
        weights_fname,
        example_oco2_fname
    ):
        # get grid geometry from example OCO2 file
        g = self._latlon_from_example(example_oco2_fname)

        # load list of all catchment tiles
        catch_info = pd.read_csv(catch_tile_info_fname)
        catch_info = catch_info.set_index('tile_id')
        all_catch_tiles = np.array(catch_info.index)
        self.all_catch_tiles = all_catch_tiles

        if not load_weights:
            # build overlap weights mapping
            print('building OCO2→catchment weights ...')
            catch_pixel_weights = {}
            for ti in all_catch_tiles:
                if ti % 1000 == 0:
                    print(f'  working pixel {ti}')
                lat_box = [
                    catch_info.loc[ti, 'min_lat'],
                    catch_info.loc[ti, 'max_lat'],
                ]
                lon_box = [
                    catch_info.loc[ti, 'min_lon'],
                    catch_info.loc[ti, 'max_lon'],
                ]
                com_lat = catch_info.loc[ti, 'com_lat']
                com_lon = catch_info.loc[ti, 'com_lon']

                # find intersecting OCO2 rows/cols
                lat_inter = self._intersect_indices_1d(
                    g['lat_lower'], g['lat_upper'], lat_box
                )
                lon_inter = self._intersect_indices_1d(
                    g['lon_left'], g['lon_right'], lon_box
                )

                n_la = len(lat_inter)
                n_lo = len(lon_inter)
                total = n_la * n_lo
                areas = np.zeros(total)
                flux_pixels = np.zeros((total, 2), dtype=int)

                # accumulate overlap "areas" in lat/lon-space
                k = 0
                for j, jj in enumerate(lon_inter):
                    # width along lon
                    if j == 0:
                        w = np.abs(g['lon_right'][jj] - lon_box[0])
                    elif j + 1 == n_lo:
                        w = np.abs(g['lon_left'][jj] - lon_box[1])
                    else:
                        w = np.abs(
                            g['lon_right'][jj] -
                            g['lon_right'][lon_inter[j-1]]
                        )
                    for i, ii in enumerate(lat_inter):
                        # height along lat
                        if i == 0:
                            h = np.abs(g['lat_lower'][ii] - lat_box[1])
                        elif i + 1 == n_la:
                            h = np.abs(g['lat_upper'][ii] - lat_box[0])
                        else:
                            h = np.abs(
                                g['lat_lower'][ii] -
                                g['lat_lower'][lat_inter[i-1]]
                            )
                        areas[k] = w * h
                        flux_pixels[k, 0] = ii  # lat idx
                        flux_pixels[k, 1] = jj  # lon idx
                        k += 1

                weights = areas / np.sum(areas) if np.sum(areas) > 0 else areas

                # centers/bounds for those pixels
                num = flux_pixels.shape[0]
                flux_coords = np.zeros((num, 2))
                flux_lower = np.zeros(num)
                flux_upper = np.zeros(num)
                flux_left = np.zeros(num)
                flux_right = np.zeros(num)
                for k in range(num):
                    ii = flux_pixels[k, 0]
                    jj = flux_pixels[k, 1]
                    flux_coords[k, 0] = g['center_lat'][ii]
                    flux_coords[k, 1] = g['center_lon'][jj]
                    flux_lower[k] = g['lat_lower'][ii]
                    flux_upper[k] = g['lat_upper'][ii]
                    flux_left[k] = g['lon_left'][jj]
                    flux_right[k] = g['lon_right'][jj]

                this_catch = {
                    'fluxcom_all_lon': g['center_lon'],
                    'fluxcom_all_lat': g['center_lat'],
                    'fluxcom_pixels': flux_pixels,
                    'weights': weights,
                    'coords': flux_coords,
                    'coords_lower': flux_lower,
                    'coords_upper': flux_upper,
                    'coords_left': flux_left,
                    'coords_right': flux_right,
                    'catch_coords': np.array([com_lon, com_lat]),
                    'catch_lower': lat_box[0],
                    'catch_upper': lat_box[1],
                    'catch_left':  lon_box[0],
                    'catch_right': lon_box[1],
                }
                catch_pixel_weights[ti] = this_catch

            self.fluxcom_to_catch = catch_pixel_weights
            print('weights dictionary created.')
            if save_weights:
                with open(weights_fname, 'wb') as f:
                    pickle.dump(catch_pixel_weights, f)
                print(f'weights saved: {weights_fname}')
        else:
            with open(weights_fname, 'rb') as f:
                catch_pixel_weights = pickle.load(f)
            self.fluxcom_to_catch = catch_pixel_weights
            self.all_catch_tiles = np.array(list(catch_pixel_weights.keys()))
            print(f'weights loaded: {weights_fname}')

    @staticmethod
    def _intersect_indices_1d(low_edges, high_edges, box):
        # return indices whose cell intersects [box_min, box_max]
        in_right = np.where(high_edges > box[0])[0]
        in_left = np.where(low_edges < box[1])[0]
        inter = np.intersect1d(in_right, in_left)
        return inter.tolist()

    def create_truth_oco2(
        self,
        start,
        end,
        oco2_fname_fmt,
        save_selected_truth,
        truth_fname_selected_tiles
    ):
        """
        Build daily regridded OCO2 SIF at your catchment tiles.
        For each date in [start, end], if a DOY exists in that
        year's file we compute a weighted avg; else NaN.
        """
        weights = self.fluxcom_to_catch

        curr = copy.deepcopy(start)
        days_processed = 0
        truth_df_selected = None
        last_year = None
        ds = None
        sif_3d = None
        doys = None

        while curr <= end:
            print(f'start: {start}, end: {end}, curr: {curr}')
            y = curr.year
            if y != last_year:
                fpath = oco2_fname_fmt.format(y)
                if os.path.exists(fpath):
                    ds = nc.Dataset(fpath)
                    sif_3d = np.array(ds['all_daily_sif'])
                    # read DOY variable (assumed 1..N)
                    doys = np.array(ds['doy']).astype(int)
                else:
                    ds = None
                    sif_3d = None
                    doys = None
                last_year = y

            # prepare one row (selected tiles only)
            day_vals = np.repeat(np.nan, len(self.tiles))

            if (ds is not None) and (doys is not None):
                jan1 = datetime.date(y, 1, 1)
                doy = (curr - jan1).days + 1
                # is this DOY present in file?
                idx = np.where(doys == doy)[0]
                if idx.size > 0:
                    k = int(idx[0])
                    # extract SIF slice for that DOY
                    day_sif = sif_3d[k, :, :]
                    # replace missing with NaN
                    day_sif = np.where(day_sif == -999.9,
                                       np.nan, day_sif)

                    # weighted avg for each selected tile
                    for p, pi in enumerate(self.tiles):
                        pix_idx = weights[pi]['fluxcom_pixels']
                        wts = weights[pi]['weights'].copy()
                        vals = np.zeros(len(wts))
                        for t in range(len(wts)):
                            ii = int(pix_idx[t, 0])
                            jj = int(pix_idx[t, 1])
                            vals[t] = day_sif[ii, jj]
                        mask = ~np.isnan(vals)
                        if mask.sum() < (0.67 * len(vals)):
                            day_vals[p] = np.nan
                        else:
                            day_vals[p] = np.average(vals[mask],
                                                     weights=wts[mask])
                    # no unit conversion for SIF
                # else: leave NaNs if DOY absent

            # write into DF
            date_str = curr.strftime('%Y-%m-%d')
            if days_processed == 0:
                truth_df_selected = pd.DataFrame(
                    day_vals[np.newaxis, :],
                    index=[date_str],
                    columns=self.tiles
                )
            else:
                truth_df_selected.loc[date_str] = day_vals

            curr += datetime.timedelta(days=1)
            days_processed += 1

        truth_df_selected.index.name = 'time'
        if save_selected_truth:
            s = start.strftime('%Y-%m-%d')
            e = end.strftime('%Y-%m-%d')
            out = truth_fname_selected_tiles.format(start=s, end=e)
            os.makedirs(os.path.dirname(out), exist_ok=True)
            truth_df_selected.to_csv(out)
            print(f'OCO2 SIF truth (selected tiles) saved: {out}')

