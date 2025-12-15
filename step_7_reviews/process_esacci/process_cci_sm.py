import os
import re
import netCDF4 as nc
import pandas as pd
import numpy as np
import datetime
import copy
import pickle

# Compile once: filenames contain "...-YYYYMMDDhhmmss-..."
_FNAME_DT_RE = re.compile(r'(\d{8})\d{6}')  # capture YYYYMMDD, ignore hhmmss


class process_cci_sm:
    def __init__(self, base_dir):
        self.base_dir = base_dir

    # -------------------------
    # Tile / grid setup
    # -------------------------
    def get_tiles(self, tiles_fname):
        tiles = pd.read_csv(tiles_fname, header=None)
        tiles = np.array(tiles).astype(int).T[0]
        self.tiles = tiles
        self.tiles_idx = tiles - 1

    def _latlon_from_example(self, example_cci_sm_fname):
        # Open one file to get centers & bounds (regular 0.25° grid)
        ds = nc.Dataset(example_cci_sm_fname)
        lon = np.array(ds['lon'])
        lat = np.array(ds['lat'])
        dlon = float(np.abs(lon[1] - lon[0]))
        dlat = float(np.abs(lat[1] - lat[0]))
        hdlon = dlon / 2.0
        hdlat = dlat / 2.0
        return {
            'center_lon': lon,
            'center_lat': lat,
            'lon_left':  lon - hdlon,
            'lon_right': lon + hdlon,
            'lat_lower': lat - hdlat,
            'lat_upper': lat + hdlat,
        }

    def get_weights(
        self,
        catch_tile_info_fname,
        load_weights,
        save_weights,
        weights_fname,
        example_cci_sm_fname
    ):
        # Grid geometry from example CCI file
        g = self._latlon_from_example(example_cci_sm_fname)

        # Load list of all catchment tiles
        catch_info = pd.read_csv(catch_tile_info_fname).set_index('tile_id')
        all_catch_tiles = np.array(catch_info.index)
        self.all_catch_tiles = all_catch_tiles

        if not load_weights:
            print('building CCI-SM→catchment weights ...')
            catch_pixel_weights = {}

            for ti in all_catch_tiles:
                if ti % 1000 == 0:
                    print(f'  working pixel {ti}')

                lat_box = [catch_info.loc[ti, 'min_lat'],
                           catch_info.loc[ti, 'max_lat']]
                lon_box = [catch_info.loc[ti, 'min_lon'],
                           catch_info.loc[ti, 'max_lon']]
                com_lat = catch_info.loc[ti, 'com_lat']
                com_lon = catch_info.loc[ti, 'com_lon']

                # Intersecting rows/cols
                lat_inter = self._intersect_indices_1d(
                    g['lat_lower'], g['lat_upper'], lat_box
                )
                lon_inter = self._intersect_indices_1d(
                    g['lon_left'], g['lon_right'], lon_box
                )

                n_la, n_lo = len(lat_inter), len(lon_inter)
                total = n_la * n_lo
                areas = np.zeros(total)
                flux_pixels = np.zeros((total, 2), dtype=int)

                k = 0
                for j, jj in enumerate(lon_inter):
                    # width along lon
                    if j == 0:
                        w = abs(g['lon_right'][jj] - lon_box[0])
                    elif j + 1 == n_lo:
                        w = abs(g['lon_left'][jj] - lon_box[1])
                    else:
                        w = abs(g['lon_right'][jj] -
                                g['lon_right'][lon_inter[j-1]])
                    for i, ii in enumerate(lat_inter):
                        # height along lat
                        if i == 0:
                            h = abs(g['lat_lower'][ii] - lat_box[1])
                        elif i + 1 == n_la:
                            h = abs(g['lat_upper'][ii] - lat_box[0])
                        else:
                            h = abs(g['lat_lower'][ii] -
                                    g['lat_lower'][lat_inter[i-1]])

                        areas[k] = w * h
                        flux_pixels[k, 0] = ii  # lat idx
                        flux_pixels[k, 1] = jj  # lon idx
                        k += 1

                weights = areas / np.sum(areas) if np.sum(areas) > 0 else areas

                num = flux_pixels.shape[0]
                flux_coords = np.zeros((num, 2))
                flux_lower = np.zeros(num)
                flux_upper = np.zeros(num)
                flux_left  = np.zeros(num)
                flux_right = np.zeros(num)
                for kk in range(num):
                    ii = flux_pixels[kk, 0]
                    jj = flux_pixels[kk, 1]
                    flux_coords[kk, 0] = g['center_lat'][ii]
                    flux_coords[kk, 1] = g['center_lon'][jj]
                    flux_lower[kk] = g['lat_lower'][ii]
                    flux_upper[kk] = g['lat_upper'][ii]
                    flux_left[kk]  = g['lon_left'][jj]
                    flux_right[kk] = g['lon_right'][jj]

                catch_pixel_weights[ti] = {
                    'fluxcom_all_lon': g['center_lon'],
                    'fluxcom_all_lat': g['center_lat'],
                    'fluxcom_pixels': flux_pixels,
                    'weights': weights,
                    'coords': flux_coords,
                    'coords_lower': flux_lower,
                    'coords_upper': flux_upper,
                    'coords_left':  flux_left,
                    'coords_right': flux_right,
                    'catch_coords': np.array([com_lon, com_lat]),
                    'catch_lower': lat_box[0],
                    'catch_upper': lat_box[1],
                    'catch_left':  lon_box[0],
                    'catch_right': lon_box[1],
                }

            self.ccism_to_catch = catch_pixel_weights
            print('weights dictionary created.')
            if save_weights:
                with open(weights_fname, 'wb') as f:
                    pickle.dump(catch_pixel_weights, f)
                print(f'weights saved: {weights_fname}')
        else:
            with open(weights_fname, 'rb') as f:
                catch_pixel_weights = pickle.load(f)
            self.ccism_to_catch = catch_pixel_weights
            self.all_catch_tiles = np.array(list(catch_pixel_weights.keys()))
            print(f'weights loaded: {weights_fname}')

    @staticmethod
    def _intersect_indices_1d(low_edges, high_edges, box):
        # return indices whose cell intersects [box_min, box_max]
        in_right = np.where(high_edges > box[0])[0]
        in_left  = np.where(low_edges  < box[1])[0]
        return np.intersect1d(in_right, in_left).tolist()

    # -------------------------
    # Fast date indexing
    # -------------------------
    def _date_from_fname(self, fn):
        """
        Extract date from filename (YYYYMMDDhhmmss), returning datetime.date.
        Returns None if pattern not found.
        """
        m = _FNAME_DT_RE.search(fn)
        if not m:
            return None
        ymd = m.group(1)  # 'YYYYMMDD'
        return datetime.datetime.strptime(ymd, '%Y%m%d').date()

    def _decode_date_quick(self, path):
        """
        Fallback for files without a parsable name.
        Prefer scalar 'time'; else read a single element from 't0'
        to avoid loading the full (time,lat,lon) array.
        """
        try:
            with nc.Dataset(path) as ds:
                if 'time' in ds.variables:
                    t = ds['time']
                    return nc.num2date(
                        t[0], t.units,
                        calendar=getattr(t, 'calendar', 'standard')
                    ).date()

                if 't0' in ds.variables:
                    t0 = ds['t0']
                    fill = getattr(t0, '_FillValue', -9999.0)
                    # Read a single corner, avoiding a full plane load
                    v = t0[0, 0, 0]
                    if v == fill:
                        v = t0[0, -1, -1]
                        if v == fill:
                            return None
                    return nc.num2date(
                        v, t0.units,
                        calendar=getattr(t0, 'calendar', 'standard')
                    ).date()
        except Exception:
            return None
        return None

    def _index_files_by_date(self, root_dir, start, end, progress_every=2000):
        """
        Build {date -> [paths]} ONLY for dates within [start, end].
        Uses filename parsing; tiny fallback read if needed.
        """
        by_date = {}
        scanned = 0
        kept = 0

        for dirpath, _, filenames in os.walk(root_dir):
            for fn in sorted(filenames):  # deterministic order
                if not fn.endswith('.nc'):
                    continue
                scanned += 1
                if scanned % progress_every == 0:
                    print(f'... scanned {scanned:,} files, kept {kept:,} (within range)')

                d = self._date_from_fname(fn)
                path = os.path.join(dirpath, fn)
                if d is None:
                    d = self._decode_date_quick(path)
                    if d is None:
                        continue

                if d < start or d > end:
                    continue

                by_date.setdefault(d, []).append(path)
                kept += 1

        print(f'Indexed {kept:,} files across {len(by_date):,} unique dates.')
        return by_date

    # -------------------------
    # Data read / processing
    # -------------------------
    def _read_sm_slice(self, path):
        """Return 2D sm (lat,lon) as float with NaNs applied, else None."""
        try:
            with nc.Dataset(path) as ds:
                sm = np.array(ds['sm'])[0, :, :]  # (time=1, lat, lon)
                fill = getattr(ds['sm'], '_FillValue', -9999.0)
                sm = np.where(sm == fill, np.nan, sm)
                vmin, vmax = getattr(ds['sm'], 'valid_range', (None, None))
                if vmin is not None:
                    sm = np.where(sm < float(vmin), np.nan, sm)
                if vmax is not None:
                    sm = np.where(sm > float(vmax), np.nan, sm)
                return sm
        except Exception:
            return None

    # -------------------------
    # Main producer
    # -------------------------
    def create_truth_cci_sm(
        self,
        start,
        end,
        cci_root_dir,
        save_selected_truth,
        truth_fname_selected_tiles
    ):
        """
        Build DAILY regridded CCI ACTIVE SM (%) at your catchment tiles.
        For each calendar date in [start, end], if files exist for that
        day we average them; else NaNs. Mirrors your OCO2 daily loop with
        the same progress print per day.
        """
        weights = self.ccism_to_catch

        # Fast indexing with filename-first strategy
        print('Indexing CCI files by date ...')
        by_date = self._index_files_by_date(cci_root_dir, start, end)

        curr = copy.deepcopy(start)
        days_processed = 0
        truth_df_selected = None

        while curr <= end:
            # progress print to mirror your OCO2 create_truth_oco2
            print(f'start: {start}, end: {end}, curr: {curr}')

            # prepare one row (selected tiles only)
            day_vals = np.repeat(np.nan, len(self.tiles))

            if curr in by_date:
                # Multiple files for the same day → average across files
                per_file_tile = []

                for fpath in by_date[curr]:
                    sm2d = self._read_sm_slice(fpath)
                    if sm2d is None:
                        continue

                    vals_today = np.repeat(np.nan, len(self.tiles))
                    for p, pi in enumerate(self.tiles):
                        pix_idx = weights[pi]['fluxcom_pixels']
                        wts = weights[pi]['weights'].copy()
                        vals = np.zeros(len(wts))
                        for t in range(len(wts)):
                            ii = int(pix_idx[t, 0])
                            jj = int(pix_idx[t, 1])
                            vals[t] = sm2d[ii, jj]
                        mask = ~np.isnan(vals)
                        if mask.sum() < (0.67 * len(vals)):
                            vals_today[p] = np.nan
                        else:
                            vals_today[p] = np.average(vals[mask],
                                                       weights=wts[mask])
                    per_file_tile.append(vals_today)
                if len(per_file_tile) > 0:
                    arr = np.vstack(per_file_tile)
                    day_vals = np.nanmean(arr, axis=0)

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
            print(f'CCI SM truth (selected tiles) saved: {out}')

