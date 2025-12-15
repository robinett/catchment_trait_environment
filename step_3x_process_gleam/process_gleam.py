import sys
import netCDF4 as nc
import os
import pandas as pd
import numpy as np
import datetime
import copy
import pickle
from dateutil.relativedelta import relativedelta
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches
import matplotlib as mpl
import matplotlib.cm as cm

class process:
    def __init__(self,base_dir):
        '''
        Initializer for the class.
        Inputs:
            base_dir: the base directory where data to be used by class is stored.
            catch_pixels: pixels numbers, in terms of the first set of pixels, to do this analysis at
        Outputs:
            Nothing
        '''
        self.base_dir = base_dir
    def get_tiles(self,tiles_fname):
        # get the tiles that we are going to run as defined in step 1
        tiles = pd.read_csv(tiles_fname,header=None)
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
    def inspect_fluxcom(self,example_fluxcom_fname):
        '''
        Function that prints the general information and varaibles for a fluxcom
        .nc file.
        Inputs:
            Nothing. The file to examine must be changed by editing the code itself.
        Outputs:
            Nothing, but prints general .nc information and variables stored in .nc.
        '''
        this_file = os.path.join(
            self.base_dir,
            'data',
            example_fluxcom_fname.format('2001')
        )
        fluxcom = nc.Dataset(this_file)
        print(fluxcom)
        print(fluxcom.variables)
    def get_weights(self,catch_tile_info_fname,load_weights,save_weights,weights_fname,
                    spatial_plot_weights,spatial_plot_name,fluxcom_plot_name,
                    example_fluxcom_fname):
        '''
        Function that, given lat/lon information about one set of pixels, will give which pixels of
        the second set of pixels intersect the first set of pixels and the percent area covered by
        this second set of pixels.
        Inputs:
            catch_tile_fname: information about the first set of pixels
            load_weights: logical variables about whether to do all the calculations laid out below
                          or to load a .pkl file of these calculations done previously.
            save_weights: save the results of these calculations as a .pkl?
            weights_fname: if loading weights, the fname of these weights; if calculating weights, the
                           fname at which to save these results.
            spatial_plot_weights: logical variable, should the code plot a map showing where the pixels
                                  overlap?
            spatial_plot_name: if plotting said map, what should it be called?
        Ouputs:
            Nothing. However, saves a dictionary to self (and optionally to a .pkl file). The keys of this
            dictionary are the catchment tiles indicated by catch_pixels. Each of these entries is itself
            a dictionary. Then entries of this nested dictionary are:
                fluxcom_pixels: the fluxcom pixels that intersect this catchment pixel
                weights: the percent of the catchment pixel covered by the corresponding fluxcom pixel
                coords: the coordinates of the center of each fluxcom pixel (this and all following
                        coordinates are in lat/lon)
                coords_lower: the lower bounds of each fluxcom pixel
                coords_upper: the upper bound of each fluxcom pixel
                coords_left: the left bound of each fluxcom pixel
                coords_right: the right bound of each fluxcom pixel
                catch_coords: the center of this catchment pixel
                catch_lower: the lower bound of the catchment pixel
                catch_upper: the upper bound of the catchment pixel
                catch_left: the left bound of the catchment pixel
                cath_right: the right bound of the catchment pixel
        '''
        # get catchment pixels out of self
        catch_tiles = self.tiles
        # get an example fluxcom lat/lon
        this_file = os.path.join(
            self.base_dir,
            'data',
            example_fluxcom_fname.format('2003')
        )
        fluxcom = nc.Dataset(this_file)
        fluxcom_le = np.array(fluxcom['E'])
        summer_le = np.array(fluxcom_le[180,:,:])
        # we're going to need the lat and lon bounds
        fluxcom_lon = np.array(fluxcom['lon'])
        # let's get the spacing between lons
        lon_spacing = np.abs(fluxcom_lon[1] - fluxcom_lon[0])
        half_lon_spacing = lon_spacing/2
        # now get the upper and lower lon bounds for each cell
        fluxcom_lon_left = fluxcom_lon - half_lon_spacing
        fluxcom_lon_right = fluxcom_lon + half_lon_spacing
        # let's do the same thing for lats
        fluxcom_lat = np.array(fluxcom['lat'])
        # lets get the spacing between lats
        lat_spacing = np.abs(fluxcom_lat[1] - fluxcom_lat[0])
        half_lat_spacing = lat_spacing/2
        # now get the upper and lower bounds for each cell
        fluxcom_lat_upper = fluxcom_lat + half_lat_spacing
        fluxcom_lat_lower = fluxcom_lat - half_lat_spacing
        # might as well get the centers to be consistant
        center_lon = fluxcom_lon
        center_lat = fluxcom_lat
        if load_weights == False:
            # get the lon and lats for the fluxcom tiles
            #fluxcom_lon = np.array(fluxcom['lon_bnds'])
            #fluxcom_lon_left = fluxcom_lon[:,1]
            #fluxcom_lon_right = fluxcom_lon[:,0]
            #fluxcom_lat = np.array(fluxcom['lat_bnds'])
            #fluxcom_lat_upper = fluxcom_lat[:,1]
            #fluxcom_lat_lower = fluxcom_lat[:,0]
            #center_lon = np.array(fluxcom['lon'])
            #center_lat = np.array(fluxcom['lat'])

            # get the tile coords from catchment
            catch_tile_info = pd.read_csv(
                catch_tile_info_fname
            )
            catch_tile_info = catch_tile_info.set_index('tile_id')
            all_catch_tiles = np.array(catch_tile_info.index)
            self.all_catch_tiles = all_catch_tiles
            # loop over all catchment tiles
            # get the area intersecting of each fluxcom tile
            catch_pixel_weights = {}
            for t,ti in enumerate(all_catch_tiles):
                if ti%1000 == 0:
                    print('working for pixel {}'.format(ti))
                # get the lat/lon box for each catchment pixel
                lat_box = [catch_tile_info['min_lat'].loc[ti], catch_tile_info['max_lat'].loc[ti]]
                lon_box = [catch_tile_info['min_lon'].loc[ti], catch_tile_info['max_lon'].loc[ti]]
                com_lat = catch_tile_info['com_lat'].loc[ti]
                com_lon = catch_tile_info['com_lon'].loc[ti]
                # find the fluxcom lats/lons that are in or below/above/left/right each catchment pixel
                # intersection of these lists tells you the fluxcom pixel boundaries that are inside of
                #   a catchment pixel
                lats_in_above = np.where(fluxcom_lat_upper>lat_box[0])[0]
                lats_in_below = np.where(fluxcom_lat_upper<lat_box[1])[0]
                lat_inter = set(lats_in_above).intersection(lats_in_below)
                lat_inter = list(lat_inter)
                lats_in_above = np.where(fluxcom_lat_lower>lat_box[0])[0]
                lats_in_below = np.where(fluxcom_lat_lower<lat_box[1])[0]
                lat_inter_2 = set(lats_in_above).intersection(lats_in_below)
                lat_inter_2 = list(lat_inter_2)
                lat_inter = np.unique(
                    np.concatenate((lat_inter,lat_inter_2),0)
                )
                lat_inter = sorted(lat_inter)
                lons_in_right = np.where(fluxcom_lon_left>lon_box[0])[0]
                lons_in_left = np.where(fluxcom_lon_left<lon_box[1])[0]
                lon_inter = set(lons_in_left).intersection(lons_in_right)
                lon_inter = list(lon_inter)
                lons_in_right = np.where(fluxcom_lon_right>lon_box[0])[0]
                lons_in_left = np.where(fluxcom_lon_right<lon_box[1])[0]
                lon_inter_2 = set(lons_in_left).intersection(lons_in_right)
                lon_inter_2 = list(lon_inter_2)
                lon_inter = np.unique(
                    np.concatenate((lon_inter,lon_inter_2),0)
                )
                lon_inter = sorted(lon_inter)
                num_lat_pixels = len(lat_inter)
                num_lon_pixels = len(lon_inter)
                total_pixels = num_lat_pixels*num_lon_pixels
                areas = np.zeros(total_pixels)
                flux_pixels = np.zeros((total_pixels,2))
                for lo in range(num_lon_pixels):
                    if lo == 0:
                        dist_lo = np.abs(
                            fluxcom_lon_right[lon_inter[lo]] - lon_box[0]
                        )
                    elif lo+1 == num_lon_pixels:
                        dist_lo = np.abs(
                            fluxcom_lon_left[lon_inter[lo]] - lon_box[1]
                        )
                    else:
                        dist_lo = np.abs(
                            fluxcom_lon_right[lon_inter[lo]] -
                            fluxcom_lon_right[lon_inter[lo-1]]
                        )
                    for la in range(num_lat_pixels):
                        if la == 0:
                            dist_la = np.abs(
                                fluxcom_lat_lower[lat_inter[la]] - lat_box[1]
                            )
                        elif la+1 == num_lat_pixels:
                            dist_la = np.abs(
                                fluxcom_lat_upper[lat_inter[la]] - lat_box[0]
                            )
                        else:
                            dist_la = np.abs(
                                fluxcom_lat_lower[lat_inter[la]] -
                                fluxcom_lat_lower[lat_inter[la-1]]
                            )
                        areas[lo*num_lat_pixels+la] = dist_lo*dist_la
                        total_area = np.sum(areas)
                        weights = areas/total_area
                        flux_pixels[lo*num_lat_pixels+la,0] = lat_inter[la]
                        flux_pixels[lo*num_lat_pixels+la,1] = lon_inter[lo]
                ## create a dictionary corresponding to this catchment pixel
                num_pixels_int = np.shape(flux_pixels)[0]
                flux_coords = np.zeros((num_pixels_int,2))
                flux_lower = np.zeros(num_pixels_int)
                flux_upper = np.zeros(num_pixels_int)
                flux_left = np.zeros(num_pixels_int)
                flux_right = np.zeros(num_pixels_int)
                for i in range(np.shape(flux_pixels)[0]):
                    flux_coords[i,0] = center_lat[int(flux_pixels[i,0])]
                    flux_coords[i,1] = center_lon[int(flux_pixels[i,1])]
                    flux_lower[i] = fluxcom_lat_lower[int(flux_pixels[i,0])]
                    flux_upper[i] = fluxcom_lat_upper[int(flux_pixels[i,0])]
                    flux_left[i] = fluxcom_lon_left[int(flux_pixels[i,1])]
                    flux_right[i] = fluxcom_lon_right[int(flux_pixels[i,1])]
                #else:
                #    for i in range(np.shape(flux_pixels)[0]):
                #        flux_coords[i,0] = np.nan
                #        flux_coords[i,1] = np.nan
                #        flux_lower[i] = np.nan
                #        flux_upper[i] = np.nan
                #        flux_left[i] = np.nan
                #        flux_right[i] = np.nan
                this_catch = {
                    'fluxcom_all_lon':center_lon,
                    'fluxcom_all_lat':center_lat,
                    'fluxcom_pixels':np.array(flux_pixels),
                    'weights':weights,
                    'coords':np.array(flux_coords),
                    'coords_lower':np.array(flux_lower),
                    'coords_upper':np.array(flux_upper),
                    'coords_left':np.array(flux_left),
                    'coords_right':np.array(flux_right),
                    'catch_coords':np.array([com_lon,com_lat]),
                    'catch_lower':lat_box[0],
                    'catch_upper':lat_box[1],
                    'catch_left':lon_box[0],
                    'catch_right':lon_box[1]
                }
                # creates a dictionary of dictionaries
                # each first-level entry is a catchment pixel
                # second level entry fluxcom_pixels is a list of the fluxcom pixels that are inside of said
                #   catchment pixel
                # second level entry weights is the percent of the catchment pixel that is covered by said
                #   fluxcom pixel
                catch_pixel_weights[ti] = this_catch
            # save this dictionary to self
            self.fluxcom_to_catch = catch_pixel_weights
            print('final dictionary mapping catchment pixels to fluxcom pixels created:')
            if save_weights == True:
                # save to the described .pkl file if directed
                with open(weights_fname,'wb') as f:
                    pickle.dump(catch_pixel_weights,f)
                print('weights saved as {}'.format(weights_fname))

        else:
            with open(weights_fname,'rb') as f:
                catch_pixel_weights = pickle.load(f)
                self.fluxcom_to_catch = catch_pixel_weights
                all_catch_tiles = np.array(list(catch_pixel_weights.keys()))
                self.all_catch_tiles = all_catch_tiles
            print('weights loaded from {}'.format(weights_fname))

        if spatial_plot_weights == True:
            # if we want to plot a map of each catchment pixel and the overlapping fluxcom pixels
            print('plotting pixel location maps')
            # get the pixels of interest for plotting purposes
            catch_pixel = self.tiles
            # we will use the plot_spatial_cells from the site_specifics class
            # turn the -9999 values to nans for this general map for asthetic purposes
            summer_nans = np.where(summer_le == -999)
            summer_le[summer_nans] = np.nan
            summer_le = summer_le*28.94
            for p,pi in enumerate(catch_tiles):
                # go through and plot a map for each pixel
                print('plotting for pixel {}'.format(pi))
                # add in all of the tiles for the fluxcom dataset
                left_list = catch_pixel_weights[pi]['coords_left']
                bottom_list = catch_pixel_weights[pi]['coords_lower']
                width_list = np.abs(left_list - catch_pixel_weights[pi]['coords_right'])
                height_list = np.abs(bottom_list - catch_pixel_weights[pi]['coords_upper'])
                type_list = np.repeat('fluxcom',len(left_list))
                weights_list = catch_pixel_weights[pi]['weights']
                fluxcom_centers = catch_pixel_weights[pi]['coords']
                center_lon = catch_pixel_weights[pi]['fluxcom_all_lon']
                center_lat = catch_pixel_weights[pi]['fluxcom_all_lat']
                # add in the catchment pixel
                left_list = np.append(left_list,catch_pixel_weights[pi]['catch_left'])
                bottom_list = np.append(bottom_list,catch_pixel_weights[pi]['catch_lower'])
                width_list = np.append(
                    width_list,np.abs(left_list[-1] - catch_pixel_weights[pi]['catch_right'])
                )
                height_list = np.append(
                    height_list,np.abs(bottom_list[-1] - catch_pixel_weights[pi]['catch_upper'])
                )
                type_list = np.append(type_list,'catch')
                # map them to the color of the squares to be plotted
                type_to_color = {
                    'fluxcom':'blue',
                    'catch':'red'
                }
                savename = spatial_plot_name.format(pi)
                extent = 'one_degree'
                self.plot_spatial_cells(
                    left_list,bottom_list,width_list,height_list,savename,
                    extent=extent,background=True,text=True,
                    back_data=summer_le,back_lon=center_lon,back_lat=center_lat,
                    tile_centers=fluxcom_centers,tile_text=weights_list,type_to_color=type_to_color,
                    type_list=type_list
                )
            # let's plot some tile numbers to make sure that we have that
            # correct as well
            #for p,pi in enumerate(catch_tiles):
            #    # go through and plot a map for each pixel
            #    print('plotting for pixel {}'.format(pi))
            #    # add in all of the tiles for the fluxcom dataset
            #    left_list = catch_pixel_weights[pi]['coords_left']
            #    bottom_list = catch_pixel_weights[pi]['coords_lower']
            #    width_list = np.abs(left_list - catch_pixel_weights[pi]['coords_right'])
            #    height_list = np.abs(bottom_list - catch_pixel_weights[pi]['coords_upper'])
            #    type_list = np.repeat('fluxcom',len(left_list))
            #    weights_list = catch_pixel_weights[pi]['fluxcom_pixels']
            #    fluxcom_centers = catch_pixel_weights[pi]['coords']
            #    center_lon = catch_pixel_weights[pi]['fluxcom_all_lon']
            #    center_lat = catch_pixel_weights[pi]['fluxcom_all_lat']
            #    # add in the catchment pixel
            #    left_list = np.append(left_list,catch_pixel_weights[pi]['catch_left'])
            #    bottom_list = np.append(bottom_list,catch_pixel_weights[pi]['catch_lower'])
            #    width_list = np.append(
            #        width_list,np.abs(left_list[-1] - catch_pixel_weights[pi]['catch_right'])
            #    )
            #    height_list = np.append(
            #        height_list,np.abs(bottom_list[-1] - catch_pixel_weights[pi]['catch_upper'])
            #    )
            #    type_list = np.append(type_list,'catch')
            #    # map them to the color of the squares to be plotted
            #    type_to_color = {
            #        'fluxcom':'blue',
            #        'catch':'red'
            #    }
            #    savename = 'fluxcom_pixel_number_{}.png'.format(pi)
            #    extent = 'one_degree'
            #    self.plot_spatial_cells(
            #        left_list,bottom_list,width_list,height_list,savename,
            #        extent=extent,background=True,text=True,
            #        back_data=summer_le,back_lon=center_lon,back_lat=center_lat,
            #        tile_centers=fluxcom_centers,tile_text=weights_list,type_to_color=type_to_color,
            #        type_list=type_list
            #    )
            # lets get a conus-wide snapshot of what these le values look like on this summer day
            #   for fluxcom
            extent = 'conus'
            savename = fluxcom_plot_name
            self.plot_spatial_cells(
                left_list,bottom_list,width_list,height_list,savename,
                extent=extent,background=True,text=True,
                back_data=summer_le,back_lon=center_lon,back_lat=center_lat,
                tile_centers=fluxcom_centers,tile_text=weights_list,type_to_color=type_to_color,
                type_list=type_list
            )

    def plot_spatial_cells(self,left_list,bottom_list,width_list,height_list,savename,
                           extent='global',background=False,text=False,facecolor=False,vals=np.nan,
                           back_data=np.nan,back_lon=np.nan,back_lat=np.nan,tile_centers=np.nan,
                           tile_text=np.nan,type_to_color=np.nan,type_list=np.nan,bounds=[-200]):
        '''
        Function for plotting geospatial data where you have information about the tile
        but not a complete lat/lon grid listing
        Always need left_list,bottom_list,width_list,height_list,savename
        To plot tranparent tiles use facecolor=False
        To plot filled tiles with colorbar use facecolor=True, with vals as the val for each tile
        To plot a backgroun with transparent tiles use background=True,back_vals,back_lon,back_lat
        To add text to the center of each tranparent tile use text=True,tile_centers,tile_text
        '''
        # create the figure
        f,ax = plt.subplots(1, 1, subplot_kw=dict(projection=ccrs.PlateCarree()))
        ax.coastlines(resolution='110m')
        if background:
            # plot the background map
            minval = np.nanmin(back_data)
            maxval = np.nanmax(back_data)
            xx,yy = np.meshgrid(back_lon,back_lat)
            m = plt.pcolormesh(xx,yy,back_data,
                               vmin=0,vmax=200)
            cbar = plt.colorbar(m,cmap="PiYG",orientation = 'horizontal',aspect=35)
            cbar.set_label(label=r'',fontsize=20)
            cbar.ax.tick_params(labelsize=20)
        # plot the pixels that we want to highlight
        if not facecolor:
            for l,le in enumerate(left_list):
                ax.add_patch(mpatches.Rectangle(xy=[left_list[l],bottom_list[l]],
                                               width=width_list[l],height=height_list[l],
                                               facecolor='none',edgecolor=type_to_color[type_list[l]],
                                               alpha=0.8,transform=ccrs.PlateCarree()
                                              )
                            )
            if text:
                for c,ce in enumerate(tile_centers):
                    ax.text(ce[1],ce[0],("%.2f" % tile_text[c]),transform=ccrs.PlateCarree(),
                            fontsize=8,horizontalalignment='center',verticalalignment='center'
                           )
        elif facecolor:
            if bounds[0] == -200:
                minval = np.nanmin(vals)
                maxval = np.nanmax(vals)
                abs_max = np.nanmax([np.abs(minval),np.abs(maxval)])
                norm = mpl.colors.Normalize(vmin=-abs_max,vmax=abs_max)
            else:
                minval = bounds[0]
                maxval = bounds[1]
                norm = mpl.colors.Normalize(vmin=minval,vmax=maxval)
            cmap = cm.PiYG
            m = cm.ScalarMappable(norm=norm,cmap=cmap)
            for l,le in enumerate(left_list):
                face_val = m.to_rgba(vals[l])
                ax.add_patch(mpatches.Rectangle(xy=[left_list[l],bottom_list[l]],
                                                width=width_list[l],height=height_list[l],
                                                facecolor=face_val,transform=ccrs.PlateCarree()
                                               )
                            )
            divider = make_axes_locatable(ax)
            ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
            f.add_axes(ax_cb)
            cb = mpl.colorbar.ColorbarBase(ax_cb,cmap=cmap,norm=norm)
            cb.set_label('correlation')
        if extent == 'conus':
            # set extent causes error with this weird version of cartopy. Do
            # this manually instead
            #ax.set_extent([-127,-60,20,55],crs=ccrs.PlateCarree())
            ax.set_xlim(-127,-60)
            ax.set_ylim(20,55)
        elif extent == 'percent':
            min_lat = np.amin(lat_list)
            max_lat = np.amax(lat_list)
            min_lon = np.amin(lon_list)
            max_lon = np.amax(lon_list)
            lat_diff = max_lat - min_lat
            lon_diff = max_lon - min_lon
            lat_append = lat_diff*0.15
            lon_append = lon_diff*0.15
            # set extent causes error with this weird version of cartopy. Do
            # this manually instead
            #ax.set_extent(
            #    [min_lon-lon_append,
            #     max_lon+lon_append,
            #     min_lat-lat_append,
            #     max_lat+lat_append
            #    ],
            #    crs=ccrs.PlateCarree()
            #)
            ax.set_xlim(min_lon-lon_append,max_lon+lon_append)
            ax.set_ylim(min_lat-lat_append,max_lat+lat_append)
        elif extent == 'one_degree':
            min_lat = np.amin(bottom_list)
            max_lat = np.amax(bottom_list) + np.unique(
                np.amax(height_list[np.where(left_list == np.amax(left_list))])
            )
            min_lon = np.amin(left_list)
            max_lon = np.amax(left_list) + np.unique(
                np.amax(width_list[np.where(bottom_list == np.amax(bottom_list))])
            )
            max_lon = max_lon[0]
            max_lat = max_lat[0]
            set_min_lon = np.amax([min_lon-1,-180])
            set_max_lon = np.amin([max_lon+1,180])
            set_min_lat = np.amax([min_lat,-90])
            set_max_lat = np.amin([max_lat,90])
            # set extent causes error with this weird version of cartopy. Do
            # this manually instead
            #ax.set_extent(
            #    [set_min_lon,
            #     set_max_lon,
            #     set_min_lat,
            #     set_max_lat
            #    ],
            #    crs=ccrs.PlateCarree()
            #)
            ax.set_xlim(set_min_lon,set_max_lon)
            ax.set_ylim(set_min_lat,set_max_lat)
        elif extent == 'global':
            # set extent causes error with this weird version of cartopy. Do
            # this manually instead
            ax.set_global()
            ax.set_xlim(-180,180)
            ax.set_ylim(-90,90)
        plt.savefig(savename,bbox_inches='tight',dpi=400)
        plt.close()
    def create_truth(self,start,end,fname,save_all_truth,save_selected_truth,
                     truth_fname_all_tiles,truth_fname_selected_tiles,
                     gleam_start,gleam_end
                    ):
        '''
        Function that uses the information from get_weights() to create a truth dataset of LE in W/m2.
        Does this at a daily timescale, given daily fluxcom data
        Inputs:
            start: what day should the truth dataset start?
            end: what day should the truth dataset end (all days in between will be filled)?
            fname: what is the fname for fluxcom? see example in main() below for how this should be foramtted
            save_truth: do we save the truth dataset we have created?
            truth_fname_all_tiles: fname at which to save the truth dataset for all catchment tiles
            truth_fname_selected_tiles: fname at which to save truth dataset for selected catchment tiles
        Outpus:
            Nothing, however saves the final truth dataset if directed to do so.
        '''
        # create the current time
        curr = copy.deepcopy(start)
        # let's get a local reference variable for weights dict
        weights = self.fluxcom_to_catch
        # load the fluxcom file for the first year
        last_year = curr.year
        this_file = fname.format(last_year)
        curr_fluxcom = nc.Dataset(this_file)
        this_le = np.array(curr_fluxcom['E'])
        # loop over each day that we want fluxcom data for
        days_processed = 0
        while curr <= end:
            if curr < gleam_start or curr > gleam_end:
                #num_all_pixels = len(self.all_catch_tiles)
                num_sel_pixels = len(self.tiles)
                #day_avg_le_all = np.repeat(np.nan,num_all_pixels)
                day_avg_le_selected = np.repeat(np.nan,num_sel_pixels)
                if days_processed == 0:
                    # if this is the first time throught the loop, we need to create the dataframe
                    #truth_df_all = pd.DataFrame(
                    #    day_avg_le_all[:,np.newaxis].T,index=[curr.strftime('%Y-%m-%d')],
                    #    columns=self.all_catch_tiles
                    #)
                    truth_df_selected = pd.DataFrame(
                        day_avg_le_selected[:,np.newaxis].T,index=[curr.strftime('%Y-%m-%d')],
                        columns=self.tiles
                    )
                else:
                    # otherwise just append to the dataframe
                    #truth_df_all.loc[curr.strftime('%Y-%m-%d')] = day_avg_le_all
                    truth_df_selected.loc[curr.strftime('%Y-%m-%d')] = day_avg_le_selected
            else:
                if curr.month == 1 and curr.day == 1:
                    print('spatially averaging fluxcom to catchment for year {}'.format(curr.year))
                # if its a new year, load the fluxcom file for that year
                this_year = curr.year
                if this_year != last_year:
                    this_file = os.path.join(
                        self.base_dir,
                        'data',
                        fname.format(this_year)
                    )
                    curr_fluxcom = nc.Dataset(this_file)
                    this_le = np.array(curr_fluxcom['E'])
                # extract the LE data for that day
                # begin by getting the number of days since the beginning of the year
                delta = curr - datetime.date(this_year,1,1)
                day_num = delta.days
                day_le = this_le[day_num,:,:]
                day_avg_le_all = np.zeros(len(self.all_catch_tiles))
                day_avg_le_selected = np.zeros(len(self.tiles))
                # get the weighted average le value for each catchment pixel
                print('year {}, month {}, day {}'.format(curr.year,curr.month,curr.day))
                for p,pi in enumerate(self.tiles):
                    # extract the relevant data from the saved dictionary
                    fluxcom_pix = weights[pi]['fluxcom_pixels']
                    pix_weights = weights[pi]['weights']
                    pix_le_vals = np.zeros(len(pix_weights))
                    # loop over each fluxcom pixel in the catchment pixel
                    for w,we in enumerate(pix_weights):
                        # extract the spatial location of this fluxcom pixel
                        curr_pix_idx = fluxcom_pix[w]
                        curr_pix_idx_2 = [
                            int(p) for p in curr_pix_idx
                        ]
                        curr_pix_idx = curr_pix_idx_2
                        # get the le value from that fluxcom pixel
                        this_pix_le = day_le[curr_pix_idx[0],curr_pix_idx[1]]
                        if this_pix_le == -999:
                            this_pix_le = np.nan
                            #print(pi)
                            #print('missing')
                        elif this_pix_le < 0:
                            this_pix_le = 0
                        # add it to the array saving all fluxcom le values
                        pix_le_vals[w] = this_pix_le
                    # get the average le
                    tot_pix = len(pix_weights)
                    nan_pix_idx = np.where(
                        np.isnan(pix_le_vals) == True
                    )[0]
                    not_nan_pix_idx = np.where(
                        np.isnan(pix_le_vals) == False
                    )
                    num_nan_pix = len(nan_pix_idx)
                    perc_nan = num_nan_pix/tot_pix
                    if perc_nan > 0.33:
                        this_avg_le = np.nan
                    else:
                        pix_le_vals = pix_le_vals[not_nan_pix_idx]
                        pix_weights = pix_weights[not_nan_pix_idx]
                        this_avg_le = np.average(pix_le_vals,weights=pix_weights)
                    # convert from mm/day to W/m2
                    this_avg_le = this_avg_le*28.94
                    day_avg_le_selected[p] = this_avg_le
                #day_avg_le_selected = day_avg_le_all[self.tiles_idx]
                if days_processed == 0:
                    # if this is the first time throught the loop, we need to create the dataframe
                    #truth_df_all = pd.DataFrame(
                    #    day_avg_le_all[:,np.newaxis].T,index=[curr.strftime('%Y-%m-%d')],
                    #    columns=self.all_catch_tiles
                    #)
                    truth_df_selected = pd.DataFrame(
                        day_avg_le_selected[:,np.newaxis].T,index=[curr.strftime('%Y-%m-%d')],
                        columns=self.tiles
                    )
                else:
                    # otherwise just append to the dataframe
                    #truth_df_all.loc[curr.strftime('%Y-%m-%d')] = day_avg_le_all
                    truth_df_selected.loc[curr.strftime('%Y-%m-%d')] = day_avg_le_selected
            last_year = curr.year
            curr += datetime.timedelta(days=1)
            days_processed += 1
        #truth_df_all.index.name = 'time'
        truth_df_selected.index.name = 'time'
        if save_all_truth or save_selected_truth:
            # save the truth df as specified
            start_str = start.strftime('%Y-%m-%d')
            end_str = end.strftime('%Y-%m-%d')
            #this_truth_fname_all_tiles = truth_fname_all_tiles.format(start=start_str,end=end_str)
            this_truth_fname_selected_tiles = truth_fname_selected_tiles.format(start=start_str,end=end_str)
            #if save_all_truth:
            #    truth_df_all.to_csv(this_truth_fname_all_tiles)
            #    print(
            #        'fluxcom le truth dataset for all tiles saved to {}'.format(this_truth_fname_all_tiles)
            #    )
            if save_selected_truth:
                truth_df_selected.to_csv(this_truth_fname_selected_tiles)
                print(
                    'fluxcom le truth dataset for selected tiles saved to {}'.format(
                        this_truth_fname_selected_tiles
                    )
                )
                test = pd.read_csv(this_truth_fname_selected_tiles)
