import os, re, copy, string, glob, subprocess as sp, numba
import numpy as np
import xarray as xr
import uxarray as ux
from stockwell import st
from scipy.signal import savgol_filter
# from scipy.signal import tukey
import scipy
import cmocean
import datetime, cftime
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import hapy
from pyproj import Geod
# from sknni import SkNNI # why did I add this before...?
geod_obj = Geod(ellps='clrk66') # Use Clarke 1866 ellipsoid.
deg_to_rad = np.pi/180. ; rad_to_deg = 180./np.pi
grid_root = '/global/cfs/cdirs/m4842/whannah/files_grid'
#---------------------------------------------------------------------------------------------------
xr.set_options(use_new_combine_kwarg_defaults=True)
#---------------------------------------------------------------------------------------------------
# define constants
g    = 9.81   # m/s^2
cp   = 1004.  # J/kg/K
#---------------------------------------------------------------------------------------------------
@numba.njit
def calc_great_circle_distance(lat1,lat2,lon1,lon2):
    '''
    input should be in degrees
    output will be in km
    '''
    dlon = lon2 - lon1
    cos_dist = np.sin(lat1*deg_to_rad)*np.sin(lat2*deg_to_rad) + \
                        np.cos(lat1*deg_to_rad)*np.cos(lat2*deg_to_rad)*np.cos(dlon*deg_to_rad)
    # print( str(cos_dist.min()) +"   "+ str(cos_dist.max()) )
    # cos_dist = np.clip(cos_dist, -1.0, 1.0) # this doesn't work with my current version of numba
    cos_dist = np.where(cos_dist> 1.0, 1.0,cos_dist)
    cos_dist = np.where(cos_dist<-1.0,-1.0,cos_dist)
    EARTH_RADIUS_KM = 6371.0  # mean radius in km
    dist = np.arccos( cos_dist ) * EARTH_RADIUS_KM
    return dist
#---------------------------------------------------------------------------------------------------
@numba.njit
def calc_great_circle_bearing(lat1_in,lat2_in,lon1_in,lon2_in):
    deg_to_rad = np.pi/180.
    rad_to_deg = 180./np.pi
    lat1 = lat1_in * deg_to_rad
    lat2 = lat2_in * deg_to_rad
    lon1 = lon1_in * deg_to_rad
    lon2 = lon2_in * deg_to_rad

    dlon = lon1 - lon2

    atan_tmp1 = np.sin(lon2-lon1)*np.cos(lat2)
    atan_tmp2 = np.cos(lat1)*np.sin(lat2)-np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1)
    bearing = np.arctan2( atan_tmp1, atan_tmp2 )

    bearing = bearing * rad_to_deg
    return bearing
#---------------------------------------------------------------------------------------------------
@numba.njit()
def calc_dist_array(lat,lon,center_lat,center_lon):
    ncol = len(center_lat)
    dist = np.zeros(ncol)
    for n in range(ncol):
        dist[n] = calc_great_circle_distance(lat,center_lat[n],lon,center_lon[n])
    return dist
#---------------------------------------------------------------------------------------------------
@numba.njit()
def calc_dist_bear_array(lat,lon,center_lat,center_lon):
    ncol = len(center_lat)
    dist = np.zeros(ncol)
    bear = np.zeros(ncol)
    for n in range(ncol):
        dist[n] = calc_great_circle_distance(lat,center_lat[n],lon,center_lon[n])
        bear[n] = calc_great_circle_bearing(lat,center_lat[n],lon,center_lon[n])
    return (dist,bear)
#---------------------------------------------------------------------------------------------------
@numba.njit()
def find_closest_cells(lat,lon,center_lat,center_lon,num_cells=1):
    dist = calc_dist_array(lat,lon,center_lat,center_lon)
    min_dist_ncol = np.argsort(dist)[0:num_cells]
    return min_dist_ncol
#---------------------------------------------------------------------------------------------------
@numba.njit()
def find_closest_cells_and_dist(lat,lon,center_lat,center_lon,num_cells=1):
    dist = calc_dist_array(lat,lon,center_lat,center_lon)
    min_dist_ncol = np.argsort(dist)[0:num_cells]
    return ( min_dist_ncol, dist[min_dist_ncol] )
#---------------------------------------------------------------------------------------------------
@numba.njit()
def find_path_ncol_wgt( npts, nlev, ncll, path_lat, path_lon, center_lat, center_lon, area ):
    ncol_idx = np.zeros((npts,ncll))
    dist_wgt = np.zeros((npts,ncll)) # inverse distance
    area_wgt = np.zeros((npts,ncll))
    for n in range(npts):
        ( min_dist_ncol, min_dist_val ) = find_closest_cells_and_dist( path_lat[n], path_lon[n],
                                                                       center_lat, center_lon,
                                                                       num_cells=ncll )
        ncol_idx[n,:] = min_dist_ncol[:ncll]
        dist_wgt[n,:] = 1./min_dist_val[:ncll]
        area_wgt[n,:] = area[min_dist_ncol[:ncll]]
    return (ncol_idx, dist_wgt, area_wgt)
#---------------------------------------------------------------------------------------------------
@numba.njit()
def interpolate_to_path_numba( npts, nlev, ncll, data, path_lat, path_lon, center_lat, center_lon ):
    data_interp = np.zeros((npts,nlev))
    for n in range(npts):
        (min_dist_ncol,min_dist_val) = find_closest_cells_and_dist(path_lat[n],path_lon[n],center_lat,center_lon,num_cells=ncll)
        ncol_idx = min_dist_ncol[:ncll]
        wgt = 1./min_dist_val[:ncll]
        for k in range(nlev):
            data_interp[n,k] = np.sum( data[ncol_idx,k] * wgt ) / np.sum(wgt)
    return data_interp
#---------------------------------------------------------------------------------------------------
# def interpolate_to_path(npts,nlev,ncll,data_in,path_lat,path_lon,path_coord,lev_coord,center_lat,center_lon):
#     data_interp = interpolate_to_path_numba(npts,nlev,ncll,data_in,path_lat,path_lon,center_lat,center_lon)
#     coords_out = {'path_coord':path_coord}
#     if 'lev' in data_in.dims: coords_out['lev'] = data_in['lev']
#     data_interp = xr.DataArray(data_interp,coords=coords_out)
#     return data_interp
#---------------------------------------------------------------------------------------------------
def interpolate_to_path(npts, nlev, ncll, data_in, path_lat, path_lon, path_coord, center_lat, center_lon):
    is_dataarray = isinstance(data_in, xr.DataArray)
    data_values = data_in.values if is_dataarray else data_in

    center_lat_values = center_lat.values if isinstance(center_lat, xr.DataArray) else center_lat
    center_lon_values = center_lon.values if isinstance(center_lon, xr.DataArray) else center_lon

    data_interp = interpolate_to_path_numba(npts, nlev, ncll, data_values,
                                            path_lat, path_lon,
                                            center_lat_values, center_lon_values)

    coords_out = {'path_coord': path_coord}
    if is_dataarray and 'lev' in data_in.dims:
        coords_out['lev'] = data_in['lev']
    elif lev_coord is not None and data_interp.ndim == 2:
        coords_out['lev'] = lev_coord

    return xr.DataArray(data_interp, coords=coords_out)
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
# define path outward from a given center location
''' usage notes
(npts,path_coord,path_lat,path_lon) = calculate_path(xlat,xlon,slat,slon)
'''
def calculate_path(xlat,xlon,slat,slon,path_len_km,path_spc_km):
    # Calculate bearing between tangent and instrument locations 
    bearing1 = calc_great_circle_bearing(lat1_in=xlat,lat2_in=slat,lon1_in=xlon,lon2_in=slon)
    bearing2 = bearing1+180
    if bearing2>360: bearing2 -= 360
    # define end points for each side of the path outward along the bearing
    elon1, elat1, back_azimuth = geod_obj.fwd(xlon, xlat, bearing1, path_len_km*1e3/2.)
    elon2, elat2, back_azimuth = geod_obj.fwd(xlon, xlat, bearing2, path_len_km*1e3/2.)
    # define points along each path segment
    npts_half = int((path_len_km/2.)/path_spc_km)
    path1 = geod_obj.npts(lon1=xlon, lat1=xlat, lon2=elon1, lat2=elat1, npts=npts_half, radians=False)
    path2 = geod_obj.npts(lon1=xlon, lat1=xlat, lon2=elon2, lat2=elat2, npts=npts_half, radians=False)
    # flip one side of the path around
    path1 = path1[::-1]
    # concatenate the path along with the mid-point and end points
    path = [(elon1,elat1)] + path1 + [(xlon,xlat)] + path2 + [(elon2,elat2)]
    # extract separate arrays for lat and lon values
    path_lon = [p[0] for p in path]
    path_lat = [p[1] for p in path]
    # calculate a coordinate for other variables using the distance from the mid-point
    # path_dist = calc_dist_array(xlat,xlon,path_lat,path_lon)
    (path_dist,path_bear) = calc_dist_bear_array(xlat,xlon,path_lat,path_lon)
    path_coord = xr.DataArray(path_dist,dims='path_coord')
    path_coord = xr.where(path_bear<0,path_coord*-1,path_coord)
    npts = len(path_lat)
    return (npts,path_coord,path_lat,path_lon)
#---------------------------------------------------------------------------------------------------
def reduce_file_list_to_target_time(file_list,target_time):
    '''
    file_list       list
    target_time     cftime      yyyy-mm-dd HH:MM
    '''
    #---------------------------------------------------------------------------
    # reduce file list to what's needed for slice to speed-up loading time
    # matching_files = []
    tsec_diff_list = []
    for f in file_list:
        file_suffix = r'\.nc$'
        if ".reduced." in f: file_suffix = r'\.reduced\.nc$'
        if 'IMERG' in f:
            # ex: 3B-HHR.MS.MRG.3IMERG.20230623-S163000-E165959.0990.V07B.HDF5
            match = re.search(r'3IMERG\.(\d{4})(\d{2})(\d{2})-S(\d{2})(\d{2})(\d{2}).*?\.(\d{4})\.V', f)
            year, month, day, hour, minute, second, granule = match.groups()
            seconds = int(granule)*60
        else:
            match = re.search(r'(\d{4})-(\d{2})-(\d{2})-(\d+)'+file_suffix, f)
            year, month, day, seconds = match.groups()
        # Create datetime for the date at midnight
        base_time = datetime.datetime(int(year), int(month), int(day))
        # Add the seconds
        file_time = base_time + datetime.timedelta(seconds=int(seconds))
        # Convert to cftime for comparison
        if 'IMERG' in f:
            file_cftime = cftime.DatetimeJulian(file_time.year, file_time.month, 
                                                file_time.day, file_time.hour, 
                                                file_time.minute, 0)
        else:
            file_cftime = cftime.DatetimeNoLeap(file_time.year, file_time.month, 
                                                file_time.day, file_time.hour, 
                                                file_time.minute, 0)
        # Check if this file contains our target time (assuming 10min output)
        tsec_diff = abs((file_cftime - target_time).total_seconds())
        tsec_diff_list.append(tsec_diff)
        # print(' '*4+f'file_cftime: {file_cftime}  target_time: {target_time}  tsec_diff: {tsec_diff}')
        # if abs((file_cftime - target_time).total_seconds()) <= 600:  # within 10 min
        #     matching_files.append(f)
        #     print(' '*4+f'Time matched file: {f} -> {file_cftime}')
    # if matching_files!=[]: file_list = matching_files
    first_tsec_diff_min_index = min(range(len(tsec_diff_list)), key=lambda i: tsec_diff_list[i])
    first_tsec_diff_min_value = tsec_diff_list[first_tsec_diff_min_index]
    file_list = [file_list[first_tsec_diff_min_index]]
    #---------------------------------------------------------------------------
    return file_list
#---------------------------------------------------------------------------------------------------
def get_grid_file(case_in):
    grid_file = None
    if '256x2-eq-ind-v1' in case_in: grid_file = f'{grid_root}/2025-sohip-256x2-eq-ind-v1-pg2_scrip.nc'
    if '256x2-eq-ind-v1' in case_in: grid_file = f'{grid_root}/2025-sohip-256x2-eq-ind-v1-pg2_scrip.nc'
    if '256x2-ptgnia-v1' in case_in: grid_file = f'{grid_root}/2025-sohip-256x2-ptgnia-v1-pg2_scrip.nc'
    if '256x2-sc-ind-v1' in case_in: grid_file = f'{grid_root}/2025-sohip-256x2-sc-ind-v1-pg2_scrip.nc'
    if '256x2-sc-pac-v1' in case_in: grid_file = f'{grid_root}/2025-sohip-256x2-sc-pac-v1-pg2_scrip.nc'
    if '256x2-se-pac-v1' in case_in: grid_file = f'{grid_root}/2025-sohip-256x2-se-pac-v1-pg2_scrip.nc'
    if '256x2-sw-ind-v1' in case_in: grid_file = f'{grid_root}/2025-sohip-256x2-sw-ind-v1-pg2_scrip.nc'
    if '256x3-eq-ind-v1' in case_in: grid_file = f'{grid_root}/2025-sohip-256x3-eq-ind-v1-pg2_scrip.nc'
    if '256x3-eq-ind-v1' in case_in: grid_file = f'{grid_root}/2025-sohip-256x3-eq-ind-v1-pg2_scrip.nc'
    if '256x3-ptgnia-v1' in case_in: grid_file = f'{grid_root}/2025-sohip-256x3-ptgnia-v1-pg2_scrip.nc'
    if '256x3-sc-ind-v1' in case_in: grid_file = f'{grid_root}/2025-sohip-256x3-sc-ind-v1-pg2_scrip.nc'
    if '256x3-sc-pac-v1' in case_in: grid_file = f'{grid_root}/2025-sohip-256x3-sc-pac-v1-pg2_scrip.nc'
    if '256x3-se-pac-v1' in case_in: grid_file = f'{grid_root}/2025-sohip-256x3-se-pac-v1-pg2_scrip.nc'
    if '256x3-sw-ind-v1' in case_in: grid_file = f'{grid_root}/2025-sohip-256x3-sw-ind-v1-pg2_scrip.nc'
    return grid_file
#---------------------------------------------------------------------------------------------------
def get_coi_mask(height,vwav):
    # cone of influence (COI) for S-transform:
    # require N_cycles full cycles to fit between level and edge
    N_cycles   = 2        # increase this to be more conservative
    coi_top    = (height.max() - height) / N_cycles
    coi_bottom = (height - height.min()) / N_cycles
    coi        = np.minimum(coi_top, coi_bottom)
    coi_mask   = vwav[:, np.newaxis] <= coi[np.newaxis, :]
    return coi_mask
#---------------------------------------------------------------------------------------------------
def calc_gw_ep(da, poly_order=7, min_lz=None, max_lz=5000., bkgd=None):
    """
    Calculate gravity wave potential energy Ep(z) from temperature profiles
    using the Stockwell transform.
    ----------------------------------------------------------------------------
    da : xr.DataArray
        Temperature data with dims (time, path_coord, height) or any subset.
        Height coordinate must be in meters, uniformly spaced.
    poly_order : int
        Polynomial order for background removal.
    min_lz : float or None
        Minimum vertical wavelength [m] to include. Default 2*dz.
    max_lz : float
        Maximum vertical wavelength [m] to include. Default 5000.
    ----------------------------------------------------------------------------
    Returns: ep : xr.DataArray
        GW potential energy, same dims as da except height is preserved.
        Units are dimensionless (T'/T_bar)^2, multiply by 0.5*(g/N)^2 for J/kg.
    """
    height = da.height.values
    dz     = height[1] - height[0]

    if min_lz is None: min_lz = dz * 2
    #---------------------------------------------------------------------------
    # determine which dims to loop over
    loop_dims = [d for d in da.dims if d != 'height']
    #---------------------------------------------------------------------------
    # stack all non-height dims for simple iteration
    if loop_dims:
        da_stacked = da.stack(sample=loop_dims)  # (height, sample)
    else:
        da_stacked = da.expand_dims('sample').stack(sample=['sample'])
    n_samples = da_stacked.sizes['sample']
    ep_vals   = np.zeros((len(height), n_samples))
    #---------------------------------------------------------------------------
    if bkgd is not None:
        if loop_dims:
            bkgd_stacked = bkgd.stack(sample=loop_dims)
        else:
            bkgd_stacked = bkgd.expand_dims('sample').stack(sample=['sample'])
    #---------------------------------------------------------------------------
    # reduce height dimension based on missing/NaN values
    min_hgt_pts = len(height)
    nan_mask = np.zeros((len(height), n_samples), dtype=bool)
    for i in range(n_samples):
        tmp_data = da_stacked.isel(sample=i).values
        nan_mask[:,i] = ~np.isnan(tmp_data) # Create a boolean mask for non-NaN values
        # tmp_data = tmp_data[nan_mask[:,1]]
        # min_hgt_pts = min(min_hgt_pts,len(tmp_data))
    # da_stacked = da_stacked.isel(height=slice(0,min_hgt_pts+1))
    # height     = height[0:min_hgt_pts+1]
    #---------------------------------------------------------------------------
    # vertical wavelength axis and masks
    vert_frq  = np.fft.rfftfreq(len(height), d=dz)
    vwav      = 1.0 / vert_frq[1:]
    vwav_mask = (vwav >= min_lz) & (vwav <= max_lz)
    vwav      = vwav[vwav_mask]
    #---------------------------------------------------------------------------
    # # cone of influence
    # dist_edge = np.minimum(height - height.min(), height.max() - height)
    # coi_mask  = vwav[:, np.newaxis] <= 2 * dist_edge[np.newaxis, :]  # (n_wav, n_hgt)
    #---------------------------------------------------------------------------
    # #---------------------------------------------------------------------------
    # # determine which dims to loop over
    # loop_dims = [d for d in da.dims if d != 'height']
    # #---------------------------------------------------------------------------
    # # stack all non-height dims for simple iteration
    # if loop_dims:
    #     da_stacked = da.stack(sample=loop_dims)  # (height, sample)
    # else:
    #     da_stacked = da.expand_dims('sample').stack(sample=['sample'])
    # n_samples = da_stacked.sizes['sample']
    # ep_vals   = np.zeros((len(height), n_samples))
    #---------------------------------------------------------------------------
    for i in range(n_samples):
        profile = da_stacked.isel(sample=i).values
        #--------------------------------------------------
        # restrict to valid values
        profile = profile[nan_mask[:,i]]
        tmp_height = height[nan_mask[:,i]]
        #--------------------------------------------------
        # vertical wavelength axis and masks
        vert_frq  = np.fft.rfftfreq(len(tmp_height), d=dz)
        vwav      = 1.0 / vert_frq[1:]
        vwav_mask = (vwav >= min_lz) & (vwav <= max_lz)
        vwav      = vwav[vwav_mask]
        coi_mask = get_coi_mask(tmp_height,vwav)
        #--------------------------------------------------
        # add a simple check here that dimensions are correct
        if i==0 and profile.shape != (len(tmp_height),): profile = profile.T
        # skip if all NaN
        if np.all(np.isnan(profile)):
            ep_vals[:,i] = np.nan
            continue
        #--------------------------------------------------
        # remove background via polynomial fit
        if bkgd is None:
            bkgd_loc = np.polyval(np.polyfit(tmp_height, profile, poly_order), tmp_height)
        else:
            bkgd_loc = bkgd_stacked.values[nan_mask[:,i],i]
        anom = (profile - bkgd_loc) / bkgd_loc  # normalized T'/T_bar
        #--------------------------------------------------
        # S-transform
        S_tx = st.st(anom)           # (n_hgt, n_hgt) complex
        S_tx = S_tx[1:,:][vwav_mask] # (n_wav, n_hgt), drop zero-freq row and mask wavelengths
        #--------------------------------------------------
        # calcualte power
        S_sq = np.abs(S_tx)**2
        #--------------------------------------------------
        # calculate N for Ep
        dTdz = np.gradient(bkgd_loc, tmp_height)       # K/m
        N2   = (g / bkgd_loc) * (dTdz + g/cp)      # rad^2/s^2
        N2   = np.maximum(N2, 1e-6)            # guard against negative values near tropopause
        N    = np.sqrt(N2)                     # rad/s
        #--------------------------------------------------
        # calculate Ep
        ep_vals[nan_mask[:,i],i] = 0.5 * np.sum((g/N)**2 * S_sq * coi_mask, axis=0)
        #--------------------------------------------------
        ep_avg = np.mean(ep_vals[nan_mask[:,i],i])
        if ep_avg==np.nan:
            # print(f'ep_avg: {ep_avg}   ep_max: {ep_max}')
            hapy.print_stat(ep_vals[:,i],name='ep_vals')
            hapy.print_stat(S_sq,name='S_sq')
            hapy.print_stat(N,name='N')
            hapy.print_stat(coi_mask,name='coi_mask')
            hapy.print_stat(anom,name='anom')
            hapy.print_stat(bkgd_loc,name='bkgd_loc')
            hapy.print_stat(S_tx,name='S_tx')
            hapy.print_stat(S_sq,name='S_sq')
            # hapy.print_stat()
            exit()
            
    #---------------------------------------------------------------------------
    # unstack back to original shape
    coords = {'height': height, 'sample': da_stacked.coords['sample']}
    ep_stacked = xr.DataArray( ep_vals, dims=('height', 'sample'), coords=coords )
    ep = ep_stacked.unstack('sample').transpose(*da.dims)
    ep.attrs['long_name'] = 'GW potential energy (dimensionless)'
    ep.attrs['note']      = f'poly_order={poly_order}, min_lz={min_lz}m, max_lz={max_lz}m'
    #---------------------------------------------------------------------------
    return ep
#---------------------------------------------------------------------------------------------------
