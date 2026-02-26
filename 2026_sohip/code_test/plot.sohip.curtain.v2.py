import os, ngl, subprocess as sp, numpy as np, xarray as xr, dask, copy, string, cmocean, numba
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
# import pg_checkerboard_utilities as pg
from sknni import SkNNI
from pyproj import Geod
geod_obj = Geod(ellps='clrk66') # Use Clarke 1866 ellipsoid.
deg_to_rad = np.pi/180. ; rad_to_deg = 180./np.pi
#---------------------------------------------------------------------------------------------------
'''
nohup time python -u  code_sohip/plot.sohip.curtain.v2.py > sohip.curtain.v2.out &
pspy ; echo ; tail sohip.curtain.v2.out
'''
#---------------------------------------------------------------------------------------------------

# data_case = f'helene.ne1024pg2_ne1024pg2.F2010-SCREAMv1.may14_helene_pert1_001'
# data_root = f'/pscratch/sd/e/ebercosh/e3sm_scratch/pm-gpu/may14_helene_pert1_001/{data_case}/run'
# data_file = f'{data_root}/scream_output.helene.3hourlyINST.ne30.h.INSTANT.nhours_x3.2024-09-27-00000.nc'
### /pscratch/sd/e/ebercosh/e3sm_scratch/pm-gpu/may14_helene_pert1_001/helene.ne1024pg2_ne1024pg2.F2010-SCREAMv1.may14_helene_pert1_001/scream_output*2024-09-27*

data_file = '/pscratch/sd/e/ebercosh/e3sm_scratch/pm-gpu/jul25_katrina/katrina.ne1024pg2_ne1024pg2.F2010-SCREAMv1.jul25_katrina.pert01.seed2/run/scream_output.katrina.3hourlyINST.h.INSTANT.nhours_x3.2005-08-28-00000.nc'


# clat1,clon1 = 40.7128, -74.0060 # new york
# clat2,clon2 = 51.5072, 0. # london

# # across italy
# clat1,clat2 = 40,45
# clon1,clon2 = 10,20

# Hurricane Katrina
slat,slon = 25.2, -86.0 # 8/28 00Z
dx = 2.5
clat1,clat2 = slat-dx,slat+dx
clon1,clon2 = slon-dx,slon+dx

path_spc_km = 1

# data_var,data_clev = 'LiqWaterPath',np.logspace(-3,0.8,num=40)
data_var,data_clev = 'wind_speed_10m',None

recalculate = False

# Define mask
mdx = 10
mlat1,mlat2 = clat1-mdx,clat2+mdx
mlon1,mlon2 = clon1-mdx,clon2+mdx
if mlon1<0: mlon1 = (mlon1+360)%360
if mlon2<0: mlon2 = (mlon2+360)%360

#---------------------------------------------------------------------------------------------------
# # get Katrina storm location from IBTrACS
# ibtracs_ds = xr.open_dataset('/global/homes/w/whannah/Research/E3SM/IBTrACS.NA.v04r01.nc')
# storm_idx = None
# for i,name in enumerate(ibtracs_ds['name'].values):
#   if name==b'KATRINA': storm_idx = i
# if storm_idx is None: raise ValueError('Storm not found!')
# print()
# # for var in ibtracs_ds.data_vars:
# #   if 'usa_' in var: print(f'{var:15}  {ibtracs_ds[var].long_name}')
# slat = ibtracs_ds['usa_lat'].isel(storm=storm_idx)
# slon = ibtracs_ds['usa_lon'].isel(storm=storm_idx)
# time = slat['time'].values
# for t in range(len(time)):
#   if np.isfinite(time[t]): print(f'{time[t]}  {slat[t].values}, {slon[t].values}')
# print()
# exit()
#---------------------------------------------------------------------------------------------------
fig_file,fig_type = os.getenv('HOME')+'/Research/E3SM/figs_sohip/sohip.curtain.v2','png'
tmp_file = '/global/homes/w/whannah/Research/E3SM/data_temp/sohip.curtain.v2.tmp_data.nc'
#---------------------------------------------------------------------------------------------------
@numba.njit
def calc_great_circle_distance(lat1,lat2,lon1,lon2):
  '''
  input should be in degrees
  '''
  dlon = lon2 - lon1
  cos_dist = np.sin(lat1*deg_to_rad)*np.sin(lat2*deg_to_rad) + \
            np.cos(lat1*deg_to_rad)*np.cos(lat2*deg_to_rad)*np.cos(dlon*deg_to_rad)
  # print( str(cos_dist.min()) +"   "+ str(cos_dist.max()) )
  cos_dist = np.where(cos_dist> 1.0, 1.0,cos_dist)
  cos_dist = np.where(cos_dist<-1.0,-1.0,cos_dist)
  dist = np.arccos( cos_dist )
  return dist
#---------------------------------------------------------------------------------------------------
@numba.njit()
def calc_dist_array(lat,lon,center_lat,center_lon):
  ncol = len(center_lat)
  dist = np.zeros(ncol)
  for n in range(ncol):
    dist[n] = calc_great_circle_distance(lat,center_lat[n],lon,center_lon[n])
  return dist
#---------------------------------------------------------------------------------------------------
def find_closest_cells(lat,lon,center_lat,center_lon,num_cells=1):
  dist = calc_dist_array(lat,lon,center_lat,center_lon)
  min_dist_ncol = np.argsort(dist)[0:num_cells]
  return min_dist_ncol
#---------------------------------------------------------------------------------------------------
def find_closest_cells_and_dist(lat,lon,center_lat,center_lon,num_cells=1):
  dist = calc_dist_array(lat,lon,center_lat,center_lon)
  min_dist_ncol = np.argsort(dist)[0:num_cells]
  return ( min_dist_ncol, dist[min_dist_ncol] )
#---------------------------------------------------------------------------------------------------
# @numba.njit()
# def interpolate_to_path(path_lat,path_lon,center_lat,center_lon,halo_size=4)
#   npts = len(path_lat)
#   data_interp = np.zeros(npts)
#   for n in range(npts):
#     (min_dist_ncol,min_dist_val) = find_closest_cells_and_dist(path_lat[n],path_lon[n],center_lat,center_lon,num_cells=4)
#     data_interp[n] = np.sum( data[min_dist_ncol] * min_dist_val ) / np.sum(min_dist_val)
#---------------------------------------------------------------------------------------------------
# calculate path mid-point for map view
Bx = np.cos(np.deg2rad(clat2)) * np.cos(np.deg2rad(clon2 - clon1))
By = np.cos(np.deg2rad(clat2)) * np.sin(np.deg2rad(clon2 - clon1))
# lat_mid = np.arctan2(np.sin(np.deg2rad(clat1)) + np.sin(np.deg2rad(clat2)), np.sqrt((np.cos(np.deg2rad(clat1)) + Bx)**2 + By**2))
lon_mid = np.deg2rad(clon1) + np.arctan2(By, np.cos(np.deg2rad(clat1)) + Bx)
#---------------------------------------------------------------------------------------------------
# define path - set point spacing for consistent distance
dist_rd = calc_great_circle_distance(clat1,clat2,clon1,clon2)
dist_km = dist_rd*180/np.pi*111
npts = int(dist_km/path_spc_km)
path = geod_obj.npts(lon1=clon1, lat1=clat1, lon2=clon2, lat2=clat2, npts=npts, radians=False)
# # Add the start and end points
# path = [(clon1, clat1)] + [(lon, lat) for lon, lat, _ in path] + [(clon2, clat2)]
path.insert(0,(clon1,clat1))
path.append(  (clon2,clat2))
path_lon = [p[0] for p in path]
path_lat = [p[1] for p in path]
npts = len(path_lat)

#---------------------------------------------------------------------------------------------------

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

plot = []

# map resources
map_res = hs.res_contour_fill_map()
map_res.tmYLLabelFontHeightF    = 0.008
map_res.tmXBLabelFontHeightF    = 0.008
map_res.lbLabelFontHeightF      = 0.01
map_res.lbLabelBarOn            = False
map_res.mpOutlineBoundarySets   = 'Geophysical'
# map_res.mpGeophysicalLineColor  = 'white'
map_res.mpDataBaseVersion       = 'MediumRes' # LowRes / MediumRes / HighRes

map_res.cnFillMode              = 'CellFill'
map_res.cnCellFillEdgeColor     = 'black'
map_res.cnFillOpacityF          = 0.0 # disable fill colors

map_res.mpCenterLonF            = lon_mid
map_res.mpLimitMode             = 'LatLon'
map_res.mpMinLatF               = np.min([clat1,clat2])#-3
map_res.mpMaxLatF               = np.max([clat1,clat2])#+3
map_res.mpMinLonF               = np.min([clon1,clon2])#-3
map_res.mpMaxLonF               = np.max([clon1,clon2])#+3

# map_res.cnFillPalette           = 'MPL_viridis'
map_res.cnFillPalette           = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
if data_clev is not None:
  map_res.cnLevelSelectionMode = "ExplicitLevels"
  map_res.cnLevels = data_clev

# # create black=>white color map for clouds
# nclr = 256
# ramp = np.linspace(0,1,nclr)
# clr_map_cld = np.zeros([nclr,4])
# clr_map_cld[:,0] = ramp
# clr_map_cld[:,1] = ramp
# clr_map_cld[:,2] = ramp
# clr_map_cld[:,3] = ramp
# tres.cnFillPalette = clr_map_cld

# marker resources
mrk_res = hs.res_xy()
mrk_res.xyMarkLineMode  = 'Markers'

lin_res = hs.res_xy()
lin_res.vpHeightF = 0.3
#---------------------------------------------------------------------------------------------------
if recalculate:
  # load SCREAM data
  data_ds = xr.open_dataset( data_file ).isel(time=1)
  data = data_ds[data_var]
  ncol_orig = len(data)
  ncol_coord = data_ds['ncol']
  #-----------------------------------------------------------------------------
  # load scrip grid data
  # scrip_ds = xr.open_dataset('/global/cfs/projectdirs/m3312/whannah/HICCUP/files_grid/scrip_ne30pg2.nc')
  # scrip_ds = xr.open_dataset('/global/cfs/projectdirs/m3312/whannah/HICCUP/files_grid/scrip_ne256pg2.nc')
  scrip_ds = xr.open_dataset('/pscratch/sd/w/whannah/e3sm_scratch/files_grid/ne1024pg2_scrip.nc')
  scrip_ds = scrip_ds.rename({'grid_size':'ncol'})
  center_lon = scrip_ds['grid_center_lon'][:]
  center_lat = scrip_ds['grid_center_lat'][:]
  corner_lon = scrip_ds['grid_corner_lon'][:,:]
  corner_lat = scrip_ds['grid_corner_lat'][:,:]
  #-----------------------------------------------------------------------------
  if data.shape!=center_lat.shape:
    print()
    print(f'  data file shape  : {data.shape}')
    print(f'  scrip data shape : {center_lat.shape}')
    raise ValueError('data shapes do not match')
  #-----------------------------------------------------------------------------
  # print(f'{hc.tclr.RED}WARNING - using dummy data!{hc.tclr.END}')
  # data = scrip_ds['grid_area'][:] # dummy plotting data
  #-----------------------------------------------------------------------------
  print()
  print(f'  mask mlat1: {mlat1}')
  print(f'  mask mlat2: {mlat2}')
  print(f'  mask mlon1: {mlon1}')
  print(f'  mask mlon2: {mlon2}')
  #-----------------------------------------------------------------------------
  # create the mask
  lat = data_ds['lat'].isel(time=0,missing_dims='ignore',drop=True)
  lon = data_ds['lon'].isel(time=0,missing_dims='ignore',drop=True)
  mask = xr.DataArray( np.ones([len(ncol_coord)],dtype=bool), coords=[('ncol', ncol_coord.values)], dims='ncol' )
  mask = mask & (lat>=mlat1) & (lat<=mlat2)
  mask = mask & (lon>=mlon1) & (lon<=mlon2)
  #-----------------------------------------------------------------------------
  # apply the mask
  data       = data.where(mask,drop=True)
  center_lon = center_lon.where(mask,drop=True).values
  center_lat = center_lat.where(mask,drop=True).values
  corner_lon = corner_lon.where(mask,drop=True).values
  corner_lat = corner_lat.where(mask,drop=True).values
  #-----------------------------------------------------------------------------
  ncol = len(data)
  #-----------------------------------------------------------------------------
  ratio = ncol/ncol_orig
  print()
  print(f'  ncol before mask : {ncol_orig}')
  print(f'  ncol after mask  : {ncol}')
  print(f'  ncol ratio       : {ratio}')
  print()
  #-----------------------------------------------------------------------------
  # update cell corner location for plotting outlines
  corner_lon_alt = np.empty([ncol,5])
  corner_lat_alt = np.empty([ncol,5])
  corner_lon_alt[:,:4] = corner_lon
  corner_lat_alt[:,:4] = corner_lat
  corner_lon_alt[:,4] = corner_lon_alt[:,0]
  corner_lat_alt[:,4] = corner_lat_alt[:,0]
  #-----------------------------------------------------------------------------
  # interpolate data to path
  print(); print('starting interpolation...')
  ncll = 20
  data_interp_dist1 = np.zeros([npts,ncll])
  data_interp_dist2 = np.zeros([npts,ncll])
  for n in range(npts):
    (min_dist_ncol,min_dist_val) = find_closest_cells_and_dist(path_lat[n],path_lon[n],center_lat,center_lon,num_cells=ncll)
    wgt1 = np.power(min_dist_val,1)
    wgt2 = np.power(min_dist_val,2)
    for nc in range(ncll):
      data_interp_dist1[n,nc] = np.sum( data[min_dist_ncol[:nc]] * wgt1[:nc] ) / np.sum(wgt1[:nc])
      data_interp_dist2[n,nc] = np.sum( data[min_dist_ncol[:nc]] * wgt2[:nc] ) / np.sum(wgt2[:nc])

  print(); print('calculating distance coordinate...')
  path_dist = calc_dist_array(path_lat[0],path_lon[0],path_lat,path_lon)


  # write to file
  ds_tmp = xr.Dataset()
  path_coord = xr.DataArray(path_dist,dims='path_coord')
  ncll_coord = xr.DataArray(np.arange(0,ncll),dims='ncll')
  coords = {'path_coord':path_coord,'ncll':ncll_coord}
  ds_tmp['data_interp_dist1'] = xr.DataArray(data_interp_dist1,coords=coords)
  ds_tmp['data_interp_dist2'] = xr.DataArray(data_interp_dist2,coords=coords)
  ds_tmp.to_netcdf(path=tmp_file,mode='w')
else:
  ds_tmp = xr.open_dataset(tmp_file)
  data_interp_dist1 = ds_tmp['data_interp_dist1'].values
  data_interp_dist2 = ds_tmp['data_interp_dist2'].values
  path_coord = ds_tmp['path_coord']

path_coord = path_coord.values * 1e3
#---------------------------------------------------------------------------------------------------

# #-------------------------------------------------------------------------------
# # create map plot
# map_res.cnFillMode    = 'CellFill'
# map_res.sfXArray      = center_lon
# map_res.sfYArray      = center_lat
# map_res.sfXCellBounds = corner_lon
# map_res.sfYCellBounds = corner_lat

# map_plot = ngl.contour_map(wks, data.values, map_res)

# plot.append(map_plot)

# #-------------------------------------------------------------------------------
# # overlay path end location markers
# mrk_res.xyMarker        = 14
# mrk_res.xyMarkerThicknessF = 2.
# mrk_res.xyMarkerSizeF   = 0.01
# mrk_res.xyMarkerColor   = 'red'
# mrk_plot = ngl.xy(wks, np.array([clon1,clon2]), np.array([clat1,clat2]), mrk_res)
# ngl.overlay( plot[0], mrk_plot )

# # overlay path markers
# mrk_res.xyMarker        = 16
# mrk_res.xyMarkerSizeF   = 0.0005 # 0.001
# mrk_res.xyMarkerColor   = 'red'
# pth_plot = ngl.xy(wks, path_lon, path_lat, mrk_res)
# ngl.overlay( plot[0], pth_plot )

# #-------------------------------------------------------------------------------
# # find nearest cells to path points
# idx_list = []
# for n in range(npts):
#   tmp_idx_list = find_closest_cells(path_lat[n],path_lon[n],center_lat,center_lon,num_cells=4)
#   for idx in tmp_idx_list:
#     if idx not in idx_list: idx_list.append(idx)
# # outline intersected cells
# lin_res.xyLineThicknessF = 6
# lin_res.xyLineColor = 'blue'
# for idx in idx_list:
#   ngl.overlay( plot[0], ngl.xy(wks, corner_lon_alt[idx,:], corner_lat_alt[idx,:], lin_res) )
# #-------------------------------------------------------------------------------

# create line plot of interpolated data
print(); print('creating plot...')
lin_res.xyLineThicknessF = 6
lin_res.tiXAxisString = 'path distance coordinate [km]'
lin_res.tiYAxisString = f'{data_var}'

lin_res.trXMinF = np.max(path_coord) * 0.2
lin_res.trXMaxF = np.max(path_coord) * 0.6

plot = [None]*4

lin_res.xyLineColor = 'red' ; tplot1 = ngl.xy(wks, path_coord, data_interp_dist1[:,4-1], lin_res)
lin_res.xyLineColor = 'blue'; tplot2 = ngl.xy(wks, path_coord, data_interp_dist2[:,4-1], lin_res)
pidx=0 ; plot[pidx]=tplot1 ; ngl.overlay( plot[pidx], tplot2 )
hs.set_subtitles(wks, plot[pidx], center_string='dist^1 (red) vs dist^2 (blue) weighting (N=4)', font_height=0.015)

lin_res.xyLineColor = 'red' ; tplot1 = ngl.xy(wks, path_coord, data_interp_dist1[:,10-1], lin_res)
lin_res.xyLineColor = 'blue'; tplot2 = ngl.xy(wks, path_coord, data_interp_dist2[:,10-1], lin_res)
pidx=1 ; plot[pidx]=tplot1 ; ngl.overlay( plot[pidx], tplot2 )
hs.set_subtitles(wks, plot[pidx], center_string='dist^1 (red) vs dist^2 (blue) weighting (N=10)', font_height=0.015)

lin_res.xyLineColor = 'red'   ; tplot1 = ngl.xy(wks, path_coord, data_interp_dist1[:, 2-1], lin_res)
lin_res.xyLineColor = 'green' ; tplot2 = ngl.xy(wks, path_coord, data_interp_dist1[:, 4-1], lin_res)
lin_res.xyLineColor = 'blue'  ; tplot3 = ngl.xy(wks, path_coord, data_interp_dist1[:, 8-1], lin_res)
lin_res.xyLineColor = 'purple'; tplot4 = ngl.xy(wks, path_coord, data_interp_dist1[:,16-1], lin_res)
pidx=2 ; plot[pidx]=tplot1
ngl.overlay( plot[pidx], tplot2 )
ngl.overlay( plot[pidx], tplot3 )
ngl.overlay( plot[pidx], tplot4 )
hs.set_subtitles(wks, plot[pidx], center_string='dist^1, N=2, 4, 8, 16 (red, green, blue, purple)', font_height=0.015)

lin_res.xyLineColor = 'red'   ; tplot1 = ngl.xy(wks, path_coord, data_interp_dist2[:, 2-1], lin_res)
lin_res.xyLineColor = 'green' ; tplot2 = ngl.xy(wks, path_coord, data_interp_dist2[:, 4-1], lin_res)
lin_res.xyLineColor = 'blue'  ; tplot3 = ngl.xy(wks, path_coord, data_interp_dist2[:, 8-1], lin_res)
lin_res.xyLineColor = 'purple'; tplot4 = ngl.xy(wks, path_coord, data_interp_dist2[:,16-1], lin_res)
pidx=3 ; plot[pidx]=tplot1
ngl.overlay( plot[pidx], tplot2 )
ngl.overlay( plot[pidx], tplot3 )
ngl.overlay( plot[pidx], tplot4 )
hs.set_subtitles(wks, plot[pidx], center_string='dist^2, N=2, 4, 8, 16 (red, green, blue, purple)', font_height=0.015)

#---------------------------------------------------------------------------------------------------
num_plot_col = 2#len(plot)
layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent = 5

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------