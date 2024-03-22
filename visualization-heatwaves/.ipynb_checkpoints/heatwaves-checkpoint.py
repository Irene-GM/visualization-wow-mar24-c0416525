from collections import defaultdict
import xarray as xr
import pandas as pd
import numpy as np
import warnings
from datetime import datetime, timedelta
import itertools
warnings.filterwarnings("ignore")
import timeit
import matplotlib.pyplot as plt
import networkx as nx
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
from netCDF4 import Dataset,num2date,date2num


start = timeit.default_timer()

pd.options.display.max_columns = 500
pd.options.display.width = 0

dic_coef_names = {"DC": "Average Degree of Centrality", "BC": "Average Betweenness Coefficient",
                  "CC": "Average Clustering Coefficient", "ID": "Average Input Degree",
                  "OD": "Average Output Degree", "IO": "Average Input/Output Degree"}

dic_vmin_vmax_cmap = {"DC": (0, 1, plt.cm.YlOrBr) , "BC": (0, 1, plt.cm.YlOrBr),
                  "CC": (0, 1, plt.cm.YlOrBr), "ID": (0, 300, plt.cm.YlOrBr),
                  "OD": (0, 100, plt.cm.YlOrBr), "IO": (-300, 100, plt.cm.RdGy)}


#############
# Functions #
#############


# Utility functions
# ----------------------------------------------------------------------------------------------------------------------

def get_bounds_for_extent_name(extent_name):
    if extent_name == "Iberia":
        return [-10, 5, 45, 35]
    elif extent_name == "Mediterranean":
        return [-10, 19, 47, 35]
    elif extent_name == "Europe":
        return [-26, 41, 71, 35.8 ]


def generate_panel_titles_from_dates(date_range):
    n = int(np.divide(len(date_range), 4)) + 1
    date_template = "{0} - {1}"
    titles = []
    for i in range(0, len(date_range), n):
        cur_dates_title = date_range[i:i+n]
        cur_title = date_template.format(cur_dates_title[0].strftime("%Y-%m-%d"), cur_dates_title[-1].strftime("%Y-%m-%d"))
        titles.append(cur_title)
    return titles

def prepare_coefficients_for_plotting(dic_stack, sd, ed):

    dic_coefs_per_period = defaultdict(list)

    date_range = pd.date_range(start=sd, end=ed, freq='D')
    right_bound = int(np.divide(len(date_range), 4)) + 1

    for key in sorted(dic_stack.keys()):

        stack_p1 = dic_stack[key][:, :, 0:right_bound*1]
        stack_p2 = dic_stack[key][:, :, right_bound*1:right_bound*2]
        stack_p3 = dic_stack[key][:, :, right_bound*2:right_bound*3]
        stack_p4 = dic_stack[key][:, :, right_bound*3:right_bound*4]

        mean_p1 = np.mean(stack_p1, axis=2)
        mean_p2 = np.mean(stack_p2, axis=2)
        mean_p3 = np.mean(stack_p3, axis=2)
        mean_p4 = np.mean(stack_p4, axis=2)

        # Patch weird numbers
        mean_p1[mean_p1 < 0] = 0
        mean_p2[mean_p2 < 0] = 0
        mean_p3[mean_p3 < 0] = 0
        mean_p4[mean_p4 < 0] = 0

        dic_coefs_per_period[key] = [mean_p1, mean_p2, mean_p3, mean_p4]

    # Manually adding the Input/Output degree
    mean_id = dic_coefs_per_period["ID"]
    mean_od = [-1 * item for item in dic_coefs_per_period["OD"]]
    mean_io = [mean_id[i] + mean_od[i] for i in range(len(mean_id))]
    dic_coefs_per_period["IO"] = mean_io

    return dic_coefs_per_period

def save_as_raster_files(dic_stack, path_folder_out, extent_name, sd, ed, resolution, dimensions):

    # This receives a dictionary with one tri-dimensional array per 'key'.
    # The Z dimension is the number of days in the heatwave. These partial
    # daily values will be stored as a single raster file, to enable a deeper
    # exploration of the heatwave.

    date_range = pd.date_range(start=sd, end=ed, freq="D")

    file_out_template = "CN_{0}_{1}_{2}_{3}_{4}.nc"

    unout = 'days since 2003-01-01 00:00:00'

    ny, nx = dimensions

    for key in sorted(dic_stack.keys()):

        short_description = "Complex network analysis: {0} ({1})".format(dic_coef_names[key], key)
        file_out_cur = file_out_template.format(extent_name, resolution, key, sd.strftime("%Y-%m-%d"), ed.strftime("%Y-%m-%d"))
        path_out = path_folder_out.format(file_out_cur)

        # Necessary to match the dimensions. Do not touch!
        arr_reshaped = np.swapaxes(np.swapaxes(dic_stack[key], 0, 2), 1, 2)

        # For some reason I had to generate a np.linspace, not a masked xr.DataArray
        # Used example here: https://stackoverflow.com/questions/55956270/convert-a-numpy-dataset-to-netcdf
        left, right, north, south = get_bounds_for_extent_name(extent_name)
        lons = np.linspace(left, right, nx)
        lats = np.linspace(south, north, ny)

        ncout = Dataset(path_out, 'w', 'NETCDF4')
        ncout.createDimension('lon', len(lons))
        ncout.createDimension('lat', len(lats))
        ncout.createDimension('time', len(date_range))

        lonvar = ncout.createVariable('lon', 'float32', ('lon'))
        lonvar[:] = lons

        latvar = ncout.createVariable('lat', 'float32', ('lat'))
        latvar[:] = lats

        timevar = ncout.createVariable('time', 'float64', ('time'))
        timevar.setncattr('units', unout)
        timevar[:] = date2num(date_range.tolist(), unout)

        coef = ncout.createVariable('coefficient', 'float32', ('time', 'lat', 'lon'))
        coef.standard_name = short_description
        coef.setncattr('units', 'N/A')
        coef.setncattr('grid_mapping', 'spatial_ref')
        coef[:] = arr_reshaped

        crs = ncout.createVariable('spatial_ref', 'i4')
        crs.spatial_ref = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'

        ncout.close()



# Visualization functions
# ----------------------------------------------------------------------------------------------------------------------



def plot_coefficients_per_period(dic_coefs_per_period, all_coords_bbox, bounds, sd, ed, extent_name, resolution, path_out):

    rows = range(0, 2)
    cols = range(0, 2)
    positions = list(itertools.product(rows, cols))

    date_range = pd.date_range(start=sd, end=ed, freq="D")
    titles = generate_panel_titles_from_dates(date_range)
    min_lat, max_lat, min_lon, max_lon = [item for item in bounds]
    lats_bbox, lons_bbox = [item for item in all_coords_bbox]

    file_out_template = "runs/AN_{0}_{1}_{2}_{3}_{4}.png"

    for key in sorted(dic_coefs_per_period.keys()):

        file_out_cur = file_out_template.format(extent_name, resolution, key, sd.strftime("%Y-%m-%d"), ed.strftime("%Y-%m-%d"))

        fig, ax = plt.subplots(nrows=2, ncols=2, subplot_kw={'projection': ccrs.PlateCarree()})
        fig.set_size_inches(15, 12, forward=True)
        plt.suptitle(dic_coef_names[key], fontsize=32, fontweight="bold")
        ocean = cfeature.NaturalEarthFeature('physical', 'ocean', scale='50m', edgecolor='none', facecolor='#CCCCCC')
        land = cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor=cfeature.COLORS['land'])

        current_coefficient = dic_coefs_per_period[key]
        vmin, vmax, cmap = dic_vmin_vmax_cmap[key]

        for i in range(len(current_coefficient)):

            r, c = positions[i]

            cb = ax[r, c].pcolormesh(lons_bbox, lats_bbox, current_coefficient[i], transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, cmap=cmap)

            # if key in ["ID", "OD"]:
            #
            # elif key in ["DC", "BC", "CC"]:
            #     cb = ax[r, c].pcolormesh(lons_bbox, lats_bbox, current_coefficient[i], transform=ccrs.PlateCarree(), cmap=plt.cm.viridis, vmin=0, vmax=1)
            # else:
            #     cb = ax[r, c].pcolormesh(lons_bbox, lats_bbox, current_coefficient[i], transform=ccrs.PlateCarree(), cmap=plt.cm.RdGy)

            # Setting up some map features
            gl = ax[r, c].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='black', alpha=0.5, linestyle='-')
            gl.top_labels = False
            gl.right_labels = False
            gl.xlabel_style = {'size': 18, 'color': 'gray'}
            gl.ylabel_style = {'size': 18, 'color': 'gray'}

            ax[r, c].set_xlim([min_lon, max_lon])
            ax[r, c].set_ylim([min_lat, max_lat])
            ax[r, c].add_feature(cfeature.BORDERS, linestyle='-', linewidth=.9)
            ax[r, c].add_feature(cfeature.RIVERS, linestyle='-', linewidth=.7)
            ax[r, c].add_feature(ocean)
            ax[r, c].set_title(titles[i], fontsize=28, loc='right')

        cbar_ax = fig.add_axes([0.92, 0.15, 0.05, 0.7])  # how close to border of graph, how centered wrt graph, width of the bar, height of bar
        colorbar = fig.colorbar(cb, cax=cbar_ax, orientation='vertical')
        colorbar.set_label(dic_coef_names[key], size=18, labelpad=25)
        cbar_ax.tick_params(labelsize=16)

        maximize()
        fig.savefig(path_out.format(file_out_cur), dpi=300)
        # plt.show()
        plt.gcf()
        plt.clf()

def plot_multiple_minimaps(dataset, coef, extent_name, resolution, sd, ed, path_out):

    file_out = "minimaps/CN_minimaps_{0}_{1}_{2}_{3}_{4}.png".format(extent_name, resolution, coef, sd.strftime("%Y-%m-%d"), ed.strftime("%Y-%m-%d"))

    ocean = cfeature.NaturalEarthFeature('physical', 'ocean', scale='50m', edgecolor='none', facecolor='lightgray')
    left, right, north, south = get_bounds_for_extent_name(extent_name)
    xlocs = np.arange(left, right, 3)
    ylocs = np.arange(south, north, 3)

    vmin, vmax, cmap = dic_vmin_vmax_cmap[coef]

    if dataset.dims["time"] == 21:
        rows = range(0, 3)
        cols = range(0, 7)
        pairs = list(itertools.product(rows, cols))
        fig, ax = plt.subplots(nrows=3, ncols=7, subplot_kw={'projection': ccrs.PlateCarree()})
        fig.set_size_inches(24, 15, forward=True)
        plt.subplots_adjust(wspace=0.4, hspace=0.1)
    elif dataset.dims["time"] == 30:
        pass

    plt.suptitle("Time-series of coefficient: {0} ({1})".format(dic_coef_names[coef], coef), fontsize=32, fontweight="bold")
    for i in range(dataset.dims["time"]):
        r, c = pairs[i]
        arr = dataset.coefficient.isel(time=i)

        img = arr.plot(ax=ax[r, c], add_colorbar=False, cmap=cmap, vmin=vmin, vmax=vmax)

        # Setting up some map features
        gl = ax[r, c].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, xlocs=xlocs, ylocs=ylocs, linewidth=0.5, color='grey', alpha=0.5, linestyle='-')
        ax[r, c].set_yticks(ylocs, crs=ccrs.PlateCarree())
        ax[r, c].set_xticks(xlocs, crs=ccrs.PlateCarree())
        ax[r, c].xaxis.set_tick_params(labelsize=14)
        ax[r, c].yaxis.set_tick_params(labelsize=14)
        ax[r, c].set_title(ax[r, c].get_title().split(" ")[2], fontsize=16)
        ax[r, c].xaxis.get_label().set_fontsize(16)
        ax[r, c].yaxis.get_label().set_fontsize(16)
        ax[r, c].xaxis.labelpad = 15
        ax[r, c].yaxis.labelpad = 15

        gl.top_labels = False
        gl.right_labels = False
        gl.left_labels = False
        gl.bottom_labels = False

        # Adding some map features
        ax[r, c].set_xlim([left, right])
        ax[r, c].set_ylim([south, north])
        ax[r, c].add_feature(cfeature.BORDERS, linestyle='-', linewidth=.9)
        ax[r, c].add_feature(ocean)

    # Adding the colorbar
    cbar = fig.colorbar(img, ax=ax.ravel().tolist(), fraction=0.03, orientation="horizontal", anchor=(1.0, 8))
    cbar.ax.tick_params(labelsize=24)

    maximize()
    fig.savefig(path_out.format(file_out), dpi=300)
    # plt.show()
    plt.gcf()
    plt.clf()


def maximize(backend=None,fullscreen=False):
    """Maximize window independently on backend.
    Fullscreen sets fullscreen mode, that is same as maximized, but it doesn't have title bar (press key F to toggle full screen mode)."""
    if backend is None:
        import matplotlib
        backend=matplotlib.get_backend()
    mng = plt.get_current_fig_manager()
    if fullscreen:
        mng.full_screen_toggle()
    else:
        if backend == 'wxAgg':
            mng.frame.Maximize(True)
        elif backend == 'Qt4Agg' or backend == 'Qt5Agg':
            mng.window.showMaximized()
        elif backend == 'TkAgg':
            # mng.window.state('iconic') #works fine on Windows! # In fedora: normal, iconic, or withdrawn
            mng.resize(*mng.window.maxsize()) # This one works on Linux/Fedora
        else:
            print ("Unrecognized backend: ",backend) #not tested on different backends (only Qt)

    plt.show(block=False)
    plt.pause(2)
    plt.close()





# Is Heatwave? Calculation
# ----------------------------------------------------------------------------------------------------------------------


def gendates(date, years_past):
    past_dates = list(reversed([datetime(date.year - i, date.month, date.day) for i in range(1, years_past + 1)]))
    return past_dates

def context_around_date(date):
    # Here, for each date in the context (i.e. (-2, -1, 0, 1, 2), we select the frame of
    # that date in the previous N years, in this case, the climatology (I think).

    left_bound = date - timedelta(days=2)
    right_bound = date + timedelta(days=2)
    date_range = pd.date_range(left_bound, right_bound, freq="D")

    dic = defaultdict(list)
    for i in range(len(date_range)):
        cur_date = date_range[i]
        past_years = gendates(cur_date, 30)
        dic[i-2] = past_years
    return dic

def get_context_dates_dict(date_range):
    dic_past_years = {}
    for cur_date in date_range:
        dic_ctxt_clim = context_around_date(cur_date)
        dic_past_years[cur_date] = dic_ctxt_clim
    return dic_past_years

def get_context_arrays_from_climatology(cur_date, dic_past_years, dataset, bounds):

    min_lat, max_lat, min_lon, max_lon = [item for item in bounds]

    dic_ctxt_dates = dic_past_years[cur_date]
    dic_ctxt_data = defaultdict(list)
    for key in sorted(dic_ctxt_dates.keys()):
        # Keeping the context date to the cube, to make comparisons in the next step
        dates_to_fetch = dic_ctxt_dates[key] + [cur_date + timedelta(days=key)]
        # arr = dataset.sel(time=dates_to_fetch, latitude=slice(min_lat, max_lat), longitude=slice(min_lon, max_lon))
        arr = dataset.sel(time=dates_to_fetch, lat=slice(min_lat, max_lat), lon=slice(min_lon, max_lon))
        # df = pd.DataFrame(dates_month.sel(time=dates_month.time[i-2])["tx"][lat_smaller_filter,lon_smaller_filter][:,:])
        newkey = (cur_date, key)
        dic_ctxt_data[newkey].append(arr)
    return dic_ctxt_data

def calculate_statistics(dic):

    dic_stat = {}

    # This loop iterates over the real dates for the heatwave period, but the
    # program also accounts for the context in the dates around it.
    for key in sorted(dic.keys()):
        pos_in_ctxt = key[1]
        cur_date = key[0]
        ctxt_date = cur_date + timedelta(days=pos_in_ctxt)
        cube = dic[key][0]

        # Select current slice and drop it from the climatology
        slice_today_ctxt = cube.sel(time=ctxt_date).to_array()[0,:,:] # No clue why is prepending an extra dimension
        cube = cube.drop_sel(time=ctxt_date).to_array()[0,:,:]

        # We compress the temporal dimension (i.e. axis=0 in this case) and for
        # each of the pixels we find the 90-th percentile. Both outputs are
        # 10x10 (in the NL case, larger if extent is larger)
        slice_avg = np.nanmean(cube, axis=0)
        slice_90p = np.percentile(cube, 90, method="closest_observation", axis=0)

        # Now we comparethe current context day to the climatology, and pixel
        # values are masked and binarized.
        slice_ctxt_compared_clim = slice_today_ctxt - slice_90p
        masked_bgtz = slice_ctxt_compared_clim.where(slice_ctxt_compared_clim >= 0, 0)
        slice_bin = masked_bgtz.where(masked_bgtz<=0, 1)

        dic_stat[key] = slice_bin

    return dic_stat


def is_heatwave(dic, date, dimensions_board):

    # Loading the 5 binary arrays associated to the context dates
    keys = [-2, -1, 0, 1, 2]
    dic_ctxt_bin = {}
    for key in keys:
        tup = (date, key)
        arr = dic[tup]
        dic_ctxt_bin[key] = arr

    # Make decision on whether it's a heatwave or not
    board = np.zeros((dimensions_board[0], dimensions_board[1]))
    today, yday, tomo, minus2, plus2 = [dic_ctxt_bin[0], dic_ctxt_bin[-1], dic_ctxt_bin[1], dic_ctxt_bin[-2], dic_ctxt_bin[2]]

    # Operations to vectorize from the if-else structure
    today_minus2 = today + minus2
    today_yday = today + yday
    today_tomo = today + tomo
    today_plus2 = today + plus2

    # Not entirely sure about this either, but it is more compact!
    part1 = ((today_minus2==2) & (today_yday==2))
    part2 = ((today_yday==2) & (today_tomo==2))
    part3 = ((today_tomo==2) & (today_plus2==2))

    where_twos_in_arrays = np.where(part1 | part2 | part3)
    board[where_twos_in_arrays] = 1

    return board


def heatwave_calc(sd, ed, dataset, bounds, dimensions_board):

    dic_heatwave = {}

    # Precompute the bounds for each date
    date_range = pd.date_range(start=sd, end=ed, freq='D')
    dic_past_years = get_context_dates_dict(date_range)

    # For the selected date_range, get the context around each of the
    # dates, also looking back in time and fetching the climatology
    # for that particular date. Then calculate the statistics for each
    # of the prepared data cubes.
    for cur_date in date_range:
        dic_ctxt_arr_clim = get_context_arrays_from_climatology(cur_date, dic_past_years, dataset, bounds)
        dic_stat = calculate_statistics(dic_ctxt_arr_clim)
        arr_hw = is_heatwave(dic_stat, cur_date, dimensions_board)
        dic_heatwave[cur_date] = arr_hw

    return dic_heatwave



# Event Synchronization Calculation
# ----------------------------------------------------------------------------------------------------------------------

def forward_calc(today, tomo):
    # Oct22, about the 1st if-clause: the original version uses the indices
    # of the matrix, thus compares the indices and if they are the same, it
    # returns a zero. This part is now vectorized, so it seems that the clause
    # can be removed now, and the required zeros added outside. Same for next
    # function.

    # if tomo == tomo:
    #     return 0.0
    if tomo - today == 1:
        return 1.0
    elif tomo + today == 2:
        return 0.5
    else:
        return 0.0

def backward_calc(today, tomo):
    # if today == today:
    #     return 0.0
    if today - tomo == 1:
        return 1.0
    elif today + tomo == 2:
        return 0.5
    else:
        return 0.0

def big_q_calc(sqrt_today_tomo, J_forw_T, J_back):
    # Make operation
    if sqrt_today_tomo != 0:
        v = ((J_forw_T + J_back) / sqrt_today_tomo)
    else:
        v = 0.0

    # Apply threshold
    if v > 1:
        return 1
    else:
        return v

def small_q_calc(sqrt_today_tomo, J_forw, J_back):
    # Make operation
    if sqrt_today_tomo != 0:
        v = ((J_forw - J_back) / sqrt_today_tomo)
    else:
        v = 0.0

    # Apply thresholds
    if v > 1:
        return 1
    elif v < -1:
        return -1
    else:
        return v

def adjacency_big_q_calc(Q):
    if Q > 0:
        return 1
    else:
        return 0

def adjacency_small_q_calc(q):
    if q > 0:
        return 1
    else:
        return 0

def eventsync_calc(sd, ed, heatwaves):

    dic_A_Q, dic_A_q = [{} for i in range(2)]
    date_range = pd.date_range(start=sd, end=ed, freq='D')[:-1]

    func_forward = np.vectorize(forward_calc, otypes=[np.float32]) # This was np.float16, but plt.imshow can't visualize it
    func_backward = np.vectorize(backward_calc, otypes=[np.float32])
    func_bigq = np.vectorize(big_q_calc, otypes=[np.float32])
    func_smallq = np.vectorize(small_q_calc, otypes=[np.float32])
    func_adj_bigq = np.vectorize(adjacency_big_q_calc, otypes=[np.float32])
    func_adj_smallq = np.vectorize(adjacency_small_q_calc, otypes=[np.float32])

    for cur_date in date_range:
        tomo_date = cur_date + timedelta(days=1)
        today, tomo = np.meshgrid(heatwaves[cur_date], heatwaves[tomo_date])

        # Calculate forward/backward event synchronization and the
        # square root of both things. Requires filling diagonal
        # with zeros
        J_ij = func_forward(today, tomo)
        J_ji = func_backward(today, tomo)
        sqrt_si_sj = np.sqrt((today.T + tomo) * (today + tomo.T))

        np.fill_diagonal(J_ij, 0)
        np.fill_diagonal(J_ji, 0)
        np.fill_diagonal(sqrt_si_sj, 0)

        Q_ij = func_bigq(sqrt_si_sj, J_ij.T, J_ji)
        q_ij = func_smallq(sqrt_si_sj, J_ij, J_ji)

        A_Q = func_adj_bigq(Q_ij)
        A_q = func_adj_smallq(q_ij)

        dic_A_Q[cur_date] = A_Q
        dic_A_q[cur_date] = A_q

    return [dic_A_Q, dic_A_q]


# Modelling data with Complex Networks
# ----------------------------------------------------------------------------------------------------------------------

def complexnetwork_modelling(sd, ed, dic_A_Q, dic_A_q):

    dic_cn_coefs = {}
    date_range = pd.date_range(start=sd, end=ed, freq='D')[:-1]

    for cur_date in date_range:

        # Fetch adjacency matrices associated to a date
        A_Q = dic_A_Q[cur_date]
        A_q = dic_A_q[cur_date]

        # Built directed and undirected graphs
        G = nx.from_numpy_array(A_Q, create_using=nx.Graph)
        g = nx.from_numpy_array(A_q, create_using=nx.DiGraph)

        # Calculate some coefficients
        degree_centrality = nx.degree_centrality(G)
        betweenness_centrality = nx.betweenness_centrality(G)
        clustering_coefficient = nx.clustering(G)
        in_degree = [item[1] for item in g.in_degree()] # Both return list of tuples
        ou_degree = [item[1] for item in g.out_degree()]

        # Pack results and return
        pack = [degree_centrality, betweenness_centrality, clustering_coefficient, in_degree, ou_degree]
        dic_cn_coefs[cur_date] = pack

    return dic_cn_coefs

def unpack_and_stack_rasters(dic_cn_coefs, date_range, all_coords_bbox):

    dic_stacks = {}

    lats_bbox, lons_bbox = [item for item in all_coords_bbox]
    stack_dc, stack_bc, stack_cc, stack_id, stack_od, stack_io = [np.empty(shape=(len(lats_bbox), len(lons_bbox))) for i in range(6)]

    for cur_date in date_range:

        dic_dc = dic_cn_coefs[cur_date][0]
        dic_bc = dic_cn_coefs[cur_date][1]
        dic_cc = dic_cn_coefs[cur_date][2]
        lst_id = dic_cn_coefs[cur_date][3]
        lst_od = dic_cn_coefs[cur_date][4]

        arr_dc = np.array(list(dic_dc.values())).reshape((len(lats_bbox), len(lons_bbox)))
        arr_bc = np.array(list(dic_bc.values())).reshape((len(lats_bbox), len(lons_bbox)))
        arr_cc = np.array(list(dic_cc.values())).reshape((len(lats_bbox), len(lons_bbox)))
        arr_id = np.array(lst_id).reshape((len(lats_bbox), len(lons_bbox)))
        arr_od = np.array(lst_od).reshape((len(lats_bbox), len(lons_bbox)))

        stack_dc = np.dstack((stack_dc, arr_dc))
        stack_bc = np.dstack((stack_bc, arr_bc))
        stack_cc = np.dstack((stack_cc, arr_cc))
        stack_id = np.dstack((stack_id, arr_id))
        stack_od = np.dstack((stack_od, arr_od))

    dic_stacks["DC"] = stack_dc
    dic_stacks["BC"] = stack_bc
    dic_stacks["CC"] = stack_cc
    dic_stacks["ID"] = stack_id
    dic_stacks["OD"] = stack_od

    return dic_stacks

def complexnetwork_summary_coefs(sd, ed, dic_cn_coefs, all_coords_bbox):

    lats_bbox, lons_bbox = [item for item in all_coords_bbox]
    date_range = pd.date_range(start=sd, end=ed, freq='D')[:-1]
    dic_stacks = unpack_and_stack_rasters(dic_cn_coefs, date_range, all_coords_bbox)

    mean_dc = np.round(np.mean(dic_stacks["DC"], axis=2), 2)
    mean_bc = np.round(np.mean(dic_stacks["BC"], axis=2), 2)
    mean_cc = np.round(np.mean(dic_stacks["CC"], axis=2), 2)
    mean_id = np.round(np.mean(dic_stacks["ID"], axis=2), 2)
    mean_od = np.round(np.mean(dic_stacks["OD"], axis=2), 2)

    pack_summary_coefs = [mean_dc, mean_bc, mean_cc, mean_id, mean_od]

    return [dic_stacks, pack_summary_coefs]


