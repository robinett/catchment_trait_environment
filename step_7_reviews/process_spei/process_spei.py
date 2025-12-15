import os
import netCDF4 as nc
import pandas as pd
import numpy as np
import datetime
import pickle


class process_spei:
    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.tiles = None            # columns to guarantee in output
        self.tiles_idx = None
        self.spei_to_catch = None    # weights dict
        self.grid_meta = None        # lat/lon centers and edges

    # -------------------------
    # Tile / grid setup
    # -------------------------
    def get_tiles(self, tiles_fname):
        """
        Load the intersecting tile list. These are the tiles that
        MUST appear as columns in the final output, even if there’s
        no overlap with the SPEI grid (in which case values are NaN).
        """
        tiles = pd.read_csv(tiles_fname, header=None)
        tiles = np.array(tiles).astype(int).T[0]
        self.tiles = tiles
        self.tiles_idx = tiles - 1

    def _latlon_from_example(self, example_spei_fname):
        """
        Read grid centers and cell edges from the example NetCDF.
        """
        with nc.Dataset(example_spei_fname) as ds:
            lon = np.array(ds["lon"])
            lat = np.array(ds["lat"])
        dlon = float(np.abs(lon[1] - lon[0]))
        dlat = float(np.abs(lat[1] - lat[0]))
        hdlon = dlon / 2.0
        hdlat = dlat / 2.0
        meta = {
            "center_lon": lon,
            "center_lat": lat,
            "lon_left": lon - hdlon,
            "lon_right": lon + hdlon,
            "lat_lower": lat - hdlat,
            "lat_upper": lat + hdlat,
        }
        self.grid_meta = meta
        return meta

    @staticmethod
    def _intersect_indices_1d(low_edges, high_edges, box):
        """
        Return indices of grid cells whose [low, high] edges
        intersect the requested [min, max] box.
        """
        in_right = np.where(high_edges > box[0])[0]
        in_left  = np.where(low_edges  < box[1])[0]
        return np.intersect1d(in_right, in_left).tolist()

    # -------------------------
    # Weights (SPEI -> tile boxes)
    # -------------------------
    def get_weights(
        self,
        catch_tile_info_fname,
        load_weights,
        save_weights,
        weights_fname,
        example_spei_fname,
    ):
        """
        Build a dictionary mapping each requested tile to the set of
        overlapping SPEI pixels and their area weights. Tiles with NO
        overlap are still included with empty intersections so they stay
        in the outputs as all-NaN.
        """
        g = self._latlon_from_example(example_spei_fname)

        catch_info = pd.read_csv(catch_tile_info_fname).set_index("tile_id")
        requested_tiles = set(self.tiles.tolist())

        if load_weights and os.path.exists(weights_fname):
            with open(weights_fname, "rb") as f:
                prev = pickle.load(f)
            # Filter to only the requested tiles; add empty stubs if missing
            weights = {}
            for ti in requested_tiles:
                if ti in prev:
                    weights[ti] = prev[ti]
                else:
                    weights[ti] = self._empty_weight_entry(ti, catch_info, g)
            self.spei_to_catch = weights
            return

        print("Building SPEI→tile weights ...")
        weights = {}
        for n, ti in enumerate(sorted(requested_tiles)):
            if (n % 1000) == 0:
                print(f"  working tile {ti}")

            if ti not in catch_info.index:
                weights[ti] = self._empty_weight_entry(ti, None, g)
                continue

            lat_box = [
                float(catch_info.loc[ti, "min_lat"]),
                float(catch_info.loc[ti, "max_lat"]),
            ]
            lon_box = [
                float(catch_info.loc[ti, "min_lon"]),
                float(catch_info.loc[ti, "max_lon"]),
            ]
            com_lat = float(catch_info.loc[ti, "com_lat"])
            com_lon = float(catch_info.loc[ti, "com_lon"])

            lat_inter = self._intersect_indices_1d(
                g["lat_lower"], g["lat_upper"], lat_box
            )
            lon_inter = self._intersect_indices_1d(
                g["lon_left"], g["lon_right"], lon_box
            )

            n_la, n_lo = len(lat_inter), len(lon_inter)
            if (n_la == 0) or (n_lo == 0):
                weights[ti] = {
                    "fluxcom_all_lon": g["center_lon"],
                    "fluxcom_all_lat": g["center_lat"],
                    "fluxcom_pixels": np.zeros((0, 2), dtype=int),
                    "weights": np.zeros(0, dtype=float),
                    "coords": np.zeros((0, 2), dtype=float),
                    "coords_lower": np.zeros(0, dtype=float),
                    "coords_upper": np.zeros(0, dtype=float),
                    "coords_left":  np.zeros(0, dtype=float),
                    "coords_right": np.zeros(0, dtype=float),
                    "catch_coords": np.array([com_lon, com_lat]),
                    "catch_lower": lat_box[0],
                    "catch_upper": lat_box[1],
                    "catch_left":  lon_box[0],
                    "catch_right": lon_box[1],
                }
                continue

            total = n_la * n_lo
            areas = np.zeros(total, dtype=float)
            flux_pixels = np.zeros((total, 2), dtype=int)
            k = 0
            for j, jj in enumerate(lon_inter):
                if j == 0:
                    w = abs(g["lon_right"][jj] - lon_box[0])
                elif (j + 1) == n_lo:
                    w = abs(g["lon_left"][jj] - lon_box[1])
                else:
                    w = abs(
                        g["lon_right"][jj] - g["lon_right"][lon_inter[j - 1]]
                    )
                for i, ii in enumerate(lat_inter):
                    if i == 0:
                        h = abs(g["lat_lower"][ii] - lat_box[1])
                    elif (i + 1) == n_la:
                        h = abs(g["lat_upper"][ii] - lat_box[0])
                    else:
                        h = abs(
                            g["lat_lower"][ii]
                            - g["lat_lower"][lat_inter[i - 1]]
                        )
                    areas[k] = w * h
                    flux_pixels[k, 0] = ii  # lat idx
                    flux_pixels[k, 1] = jj  # lon idx
                    k += 1

            if areas.sum() > 0.0:
                wts = areas / areas.sum()
            else:
                wts = np.zeros(0, dtype=float)
                flux_pixels = np.zeros((0, 2), dtype=int)

            num = flux_pixels.shape[0]
            flux_coords = np.zeros((num, 2))
            flux_lower = np.zeros(num)
            flux_upper = np.zeros(num)
            flux_left  = np.zeros(num)
            flux_right = np.zeros(num)
            for kk in range(num):
                ii = int(flux_pixels[kk, 0])
                jj = int(flux_pixels[kk, 1])
                flux_coords[kk, 0] = g["center_lat"][ii]
                flux_coords[kk, 1] = g["center_lon"][jj]
                flux_lower[kk] = g["lat_lower"][ii]
                flux_upper[kk] = g["lat_upper"][ii]
                flux_left[kk]  = g["lon_left"][jj]
                flux_right[kk] = g["lon_right"][jj]

            weights[ti] = {
                "fluxcom_all_lon": g["center_lon"],
                "fluxcom_all_lat": g["center_lat"],
                "fluxcom_pixels": flux_pixels,
                "weights": wts,
                "coords": flux_coords,
                "coords_lower": flux_lower,
                "coords_upper": flux_upper,
                "coords_left":  flux_left,
                "coords_right": flux_right,
                "catch_coords": np.array([com_lon, com_lat]),
                "catch_lower": lat_box[0],
                "catch_upper": lat_box[1],
                "catch_left":  lon_box[0],
                "catch_right": lon_box[1],
            }

        self.spei_to_catch = weights
        if save_weights:
            with open(weights_fname, "wb") as f:
                pickle.dump(weights, f)
            print(f"weights saved: {weights_fname}")

    def _empty_weight_entry(self, ti, catch_info, g):
        """
        Construct a well-formed empty weights entry for a tile.
        """
        if (catch_info is not None) and (ti in catch_info.index):
            com_lat = float(catch_info.loc[ti, "com_lat"])
            com_lon = float(catch_info.loc[ti, "com_lon"])
            lat_box = [
                float(catch_info.loc[ti, "min_lat"]),
                float(catch_info.loc[ti, "max_lat"]),
            ]
            lon_box = [
                float(catch_info.loc[ti, "min_lon"]),
                float(catch_info.loc[ti, "max_lon"]),
            ]
        else:
            com_lat = np.nan
            com_lon = np.nan
            lat_box = [np.nan, np.nan]
            lon_box = [np.nan, np.nan]

        return {
            "fluxcom_all_lon": g["center_lon"],
            "fluxcom_all_lat": g["center_lat"],
            "fluxcom_pixels": np.zeros((0, 2), dtype=int),
            "weights": np.zeros(0, dtype=float),
            "coords": np.zeros((0, 2), dtype=float),
            "coords_lower": np.zeros(0, dtype=float),
            "coords_upper": np.zeros(0, dtype=float),
            "coords_left":  np.zeros(0, dtype=float),
            "coords_right": np.zeros(0, dtype=float),
            "catch_coords": np.array([com_lon, com_lat]),
            "catch_lower": lat_box[0],
            "catch_upper": lat_box[1],
            "catch_left":  lon_box[0],
            "catch_right": lon_box[1],
        }

    # -------------------------
    # NetCDF helpers
    # -------------------------
    @staticmethod
    def _apply_scale_and_mask(var_obj, arr2d):
        """
        Apply scale_factor/add_offset, and honor valid_min/max or valid_range
        in addition to _FillValue. Returns a float array with NaNs for invalids.
        """
        arr = arr2d.astype(float)

        # _FillValue
        fill = getattr(var_obj, "_FillValue", None)
        if fill is not None:
            arr[arr == fill] = np.nan

        # scale/offset
        scale = getattr(var_obj, "scale_factor", 1.0)
        offset = getattr(var_obj, "add_offset", 0.0)
        if (scale is not None) and (scale != 1.0):
            arr = arr * float(scale)
        if (offset is not None) and (offset != 0.0):
            arr = arr + float(offset)

        # valid range
        vmin = getattr(var_obj, "valid_min", None)
        vmax = getattr(var_obj, "valid_max", None)
        vrng = getattr(var_obj, "valid_range", None)

        if vrng is not None and len(vrng) == 2:
            vmin, vmax = float(vrng[0]), float(vrng[1])

        if vmin is not None:
            arr[arr < float(vmin)] = np.nan
        if vmax is not None:
            arr[arr > float(vmax)] = np.nan

        return arr

    @staticmethod
    def _normalize_to_latlon(var_obj, arr2d):
        """
        Ensure the 2D slice is [lat, lon] for downstream [ii, jj] indexing.
        If var dims end with ('lon','lat'), transpose.
        """
        dims = var_obj.dimensions
        # Expecting dims like ('time','lat','lon') OR ('time','lon','lat')
        if len(dims) < 3:
            raise RuntimeError(f"Unexpected SPEI dims: {dims}")
        tail = dims[-2:]
        if tail == ("lat", "lon"):
            return arr2d  # already [lat, lon]
        if tail == ("lon", "lat"):
            return arr2d.T  # transpose to [lat, lon]
        raise RuntimeError(f"Unexpected 2D dims order: {tail}")

    # -------------------------
    # Producer (monthly)
    # -------------------------
    def create_truth_spei(
        self,
        start,
        end,
        spei_nc_path,
        save_selected_truth,
        truth_fname_selected_tiles,
        min_valid_fraction=0.67,
        var_name="spei_06",
        time_name="time",
    ):
        """
        Build MONTHLY regridded SPEI(6-month) at the requested tiles (self.tiles).
        Every requested tile becomes a column, even with no overlap (all-NaN).

        min_valid_fraction applies ONLY when a tile has overlapping pixels;
        otherwise the value is NaN by construction.
        """
        assert self.tiles is not None, "Tiles not loaded."
        assert self.spei_to_catch is not None, "Weights not built."

        with nc.Dataset(spei_nc_path) as ds:
            if var_name not in ds.variables:
                raise KeyError(f"Variable '{var_name}' not found in {spei_nc_path}")
            v = ds[var_name]
            if time_name not in ds.variables:
                raise KeyError(f"Time variable '{time_name}' not found in {spei_nc_path}")
            time_var = ds[time_name]
            times = nc.num2date(
                np.array(time_var),
                units=time_var.units,
                calendar=getattr(time_var, "calendar", "gregorian"),
            )

            # Map (year, month) -> index in NetCDF
            ym_to_k = {}
            for ti, tdt in enumerate(times):
                ym_to_k[(tdt.year, tdt.month)] = ti

            curr = datetime.date(start.year, start.month, 1)
            end_m = datetime.date(end.year, end.month, 1)

            idx = []
            rows = []
            tile_list = list(self.tiles)

            missing_time_count = 0

            while curr <= end_m:
                print(f"start: {start}, end: {end}, curr: {curr}")
                k = ym_to_k.get((curr.year, curr.month), None)

                # default row = NaNs for ALL tiles
                row_vals = np.full(len(tile_list), np.nan, dtype=float)

                if k is not None:
                    # read slice, normalize shape, scale/mask
                    slice2d = np.array(v[k, :, :])
                    slice2d = self._normalize_to_latlon(v, slice2d)
                    slice2d = self._apply_scale_and_mask(v, slice2d)

                    for p, ti in enumerate(tile_list):
                        entry = self.spei_to_catch.get(ti, None)
                        if entry is None:
                            # no weights entry (should not occur)
                            continue
                        pix_idx = entry["fluxcom_pixels"]
                        wts = entry["weights"]

                        if (wts is None) or (len(wts) == 0):
                            # No overlap → NaN
                            continue

                        vals = np.zeros(len(wts))
                        for t in range(len(wts)):
                            ii = int(pix_idx[t, 0])  # lat index
                            jj = int(pix_idx[t, 1])  # lon index
                            vals[t] = slice2d[ii, jj]

                        mask = ~np.isnan(vals)
                        if mask.sum() < (min_valid_fraction * len(vals)):
                            row_vals[p] = np.nan
                        else:
                            row_vals[p] = np.average(vals[mask], weights=wts[mask])
                else:
                    missing_time_count += 1

                rows.append(row_vals)
                idx.append(curr.strftime("%Y-%m-01"))

                # increment month
                if curr.month == 12:
                    curr = datetime.date(curr.year + 1, 1, 1)
                else:
                    curr = datetime.date(curr.year, curr.month + 1, 1)

        # DataFrame: all requested tiles as columns
        truth_df_selected = pd.DataFrame(
            np.vstack(rows), index=idx, columns=tile_list
        )
        truth_df_selected.index.name = "time"

        if missing_time_count > 0:
            print(f"NOTE: {missing_time_count} months not present in NetCDF; kept as NaN rows.")

        if save_selected_truth:
            s = f"{start.year}-{start.month:02d}-01"
            e = f"{end.year}-{end.month:02d}-01"
            out = truth_fname_selected_tiles.format(start=s, end=e)
            os.makedirs(os.path.dirname(out), exist_ok=True)
            truth_df_selected.to_csv(out)
            print(f"SPEI truth (selected tiles) saved: {out}")

        return truth_df_selected

