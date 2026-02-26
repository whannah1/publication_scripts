import os, ngl, subprocess as sp, numpy as np, xarray as xr, dask, copy, string, cmocean, numba
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
# import pg_checkerboard_utilities as pg
from sknni import SkNNI
from pyproj import Geod
geod_obj = Geod(ellps='clrk66') # Use Clarke 1866 ellipsoid.
deg_to_rad = np.pi/180. ; rad_to_deg = 180./np.pi
#---------------------------------------------------------------------------------------------------
'''
nohup time python -u  code_sohip/plot.sohip.curtain.v3.py > sohip.curtain.v3.ncll_3.out &
pspy ; echo ; tail sohip.curtain.v3.ncll_3.out
'''
#---------------------------------------------------------------------------------------------------

# data_case = f'helene.ne1024pg2_ne1024pg2.F2010-SCREAMv1.may14_helene_pert1_001'
# data_root = f'/pscratch/sd/e/ebercosh/e3sm_scratch/pm-gpu/may14_helene_pert1_001/{data_case}/run'
# data_file = f'{data_root}/scream_output.helene.3hourlyINST.ne30.h.INSTANT.nhours_x3.2024-09-27-00000.nc'
### /pscratch/sd/e/ebercosh/e3sm_scratch/pm-gpu/may14_helene_pert1_001/helene.ne1024pg2_ne1024pg2.F2010-SCREAMv1.may14_helene_pert1_001/scream_output*2024-09-27*

# data_file = '/pscratch/sd/e/ebercosh/e3sm_scratch/pm-gpu/jul25_katrina/katrina.ne1024pg2_ne1024pg2.F2010-SCREAMv1.jul25_katrina.pert01.seed2/run/scream_output.katrina.3hourlyINST.h.INSTANT.nhours_x3.2005-08-28-00000.nc'
data_file = '/pscratch/sd/e/ebercosh/SCREAM/Helene_ne1024/control/E0/scream_output.helene.1hourlyINST_3D.h.INSTANT.nhours_x1.2024-09-26-00000.nc'

# elat1,elon1 = 40.7128, -74.0060 # new york
# elat2,elon2 = 51.5072, 0. # london

# # across italy
# elat1,elat2 = 40,45
# elon1,elon2 = 10,20

### original method based on given end points for path
# slat,slon = 25.2, -86.0 # Hurricane Katrina - 8/28 00Z
# slat,slon = 22.96, 273.53 # Hurricane Helene - 9/26 00Z
# dx = 2.5
# elat1,elat2 = slat-dx,slat+dx
# elon1,elon2 = slon-dx,slon+dx

# new method to define path outward from center location
xlat,xlon   = 22.96, 273.53 # center point - Hurricane Helene - 9/26 00Z
slat,slon   = 27, 271       # faux ISS point used for path bearing
dist_km     = 800e3         # total path distance in meters
path_spc_km = 2             # spacing between interpolated path points

ncll = 1

# data_var,data_clev = 'LiqWaterPath',np.logspace(-3,0.8,num=40)
# data_var,data_clev = 'wind_speed_10m',None
# data_var,data_clev = 'omega',None

recalculate = False

var_3D = 'omega'
var_2D = 'wind_speed_10m'


#---------------------------------------------------------------------------------------------------
fig_file,fig_type = os.getenv('HOME')+'/Research/E3SM/figs_sohip/sohip.curtain.v3','png'
# tmp_file = '/global/homes/w/whannah/Research/E3SM/data_temp/sohip.curtain.v3.tmp_data.nc'
# tmp_file = '/pscratch/sd/w/whannah/tmp_data/sohip.curtain.v3.tmp_data.nc'
tmp_data_root = '/pscratch/sd/w/whannah/tmp_data'
tmp_file = f'{tmp_data_root}/sohip.curtain.v3.tmp_data.nc'
tmp_file = tmp_file.replace('.nc',f'.var3D_{var_3D}.var2D_{var_2D}.nc')
tmp_file = tmp_file.replace('.nc',f'.nc_{ncll}.nc')
tmp_file = tmp_file.replace('.nc',f'.x_{xlat}_{xlon}.nc')
tmp_file = tmp_file.replace('.nc',f'.s_{slat}_{slon}.nc')
tmp_file = tmp_file.replace('.nc',f'.dist_km_{int(dist_km)}.nc')
#---------------------------------------------------------------------------------------------------
print() ; print(f'tmp_file: {tmp_file}') ; print()
# exit()
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
#-----------------------------------------------------------------------------
@numba.njit()
def interpolate_to_path(npts,nlev,ncll,data,path_lat,path_lon,center_lat,center_lon):
  data_interp = np.zeros((npts,nlev))
  for n in range(npts):
    (min_dist_ncol,min_dist_val) = find_closest_cells_and_dist(path_lat[n],path_lon[n],center_lat,center_lon,num_cells=ncll)
    ncol_idx = min_dist_ncol[:ncll]
    wgt = min_dist_val[:ncll]
    for k in range(nlev):
      data_interp[n,k] = np.sum( data[ncol_idx,k] * wgt ) / np.sum(wgt)
  return data_interp
#---------------------------------------------------------------------------------------------------
# define path outward from a given center location

# Calculate bearing between tangent and instrument locations 
bearing1 = calc_great_circle_bearing(lat1_in=xlat,lat2_in=slat,lon1_in=xlon,lon2_in=slon)
bearing2 = bearing1+180
if bearing2>360: bearing2 -= 360
# define end points for each side of the path outward along the bearing
elon1, elat1, back_azimuth = geod_obj.fwd(xlon, xlat, bearing1, dist_km/2.)
elon2, elat2, back_azimuth = geod_obj.fwd(xlon, xlat, bearing2, dist_km/2.)
# define points along each path segment
npts_half = int((dist_km/2.)/path_spc_km)
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

#---------------------------------------------------------------------------------------------------
# calculate path mid-point for map view
Bx = np.cos(np.deg2rad(elat2)) * np.cos(np.deg2rad(elon2 - elon1))
By = np.cos(np.deg2rad(elat2)) * np.sin(np.deg2rad(elon2 - elon1))
lon_mid = np.deg2rad(elon1) + np.arctan2(By, np.cos(np.deg2rad(elat1)) + Bx)
#---------------------------------------------------------------------------------------------------
# Define mask limits
mdx = 10
mlat1,mlat2 = elat1-mdx,elat2+mdx
mlon1,mlon2 = elon1-mdx,elon2+mdx
if mlon1<0: mlon1 = (mlon1+360)%360
if mlon2<0: mlon2 = (mlon2+360)%360
#---------------------------------------------------------------------------------------------------

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

plot = [None]*2

# map resources
map_res = hs.res_contour_fill_map()
map_res.tmYLLabelFontHeightF    = 0.008
map_res.tmXBLabelFontHeightF    = 0.008
map_res.lbLabelFontHeightF      = 0.01
# map_res.lbLabelBarOn            = False
map_res.mpOutlineBoundarySets   = 'Geophysical'
# map_res.mpGeophysicalLineColor  = 'white'
map_res.mpDataBaseVersion       = 'MediumRes' # LowRes / MediumRes / HighRes

map_res.cnFillMode              = 'CellFill'
# map_res.cnCellFillEdgeColor     = 'black'
# map_res.cnFillOpacityF          = 0.0 # disable fill colors

map_res.mpCenterLonF            = lon_mid
map_res.mpLimitMode             = 'LatLon'
map_res.mpMinLatF               = np.min([elat1,elat2,slat])-2
map_res.mpMaxLatF               = np.max([elat1,elat2,slat])+2
map_res.mpMinLonF               = np.min([elon1,elon2,slon-360 if slon>180 else slon]) - 4
map_res.mpMaxLonF               = np.max([elon1,elon2,slon-360 if slon>180 else slon]) + 4

# map_res.cnFillPalette           = 'MPL_viridis'
map_res.cnFillPalette           = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )

#---------------------------------------------------------------------------------------------------
# print()
# print(f'  elat1: {elat1}')
# print(f'  elat2: {elat2}')
# print(f'  slat : {slat}')
# print(f'  elon1: {elon1}')
# print(f'  elon2: {elon2}')
# print(f'  slon : {slon}  ({(slon-360 if slon>180 else slon)})')
# print()
# print(f'  min lat: {np.min([elat1,elat2,slat])}')
# print(f'  max lat: {np.max([elat1,elat2,slat])}')
# print(f'  min lon: {np.min([elon1,elon2,slon-360 if slon>180 else slon])}')
# print(f'  max lon: {np.max([elon1,elon2,slon-360 if slon>180 else slon])}')
# print()
# exit()
#---------------------------------------------------------------------------------------------------

# if data_clev is not None:
#   map_res.cnLevelSelectionMode = "ExplicitLevels"
#   map_res.cnLevels = data_clev

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
  data    = data_ds[var_3D]
  data_2D = data_ds[var_2D]
  ncol_orig = len(data_ds.ncol)
  ncol_coord = data_ds['ncol']
  nlev = len(data_ds.lev)
  lev_coord = data_ds['lev']
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
  # if data.shape!=center_lat.shape:
  #   print()
  #   print(f'  data file shape  : {data.ncol.shape}')
  #   print(f'  scrip data shape : {center_lat.shape}')
  #   raise ValueError('data shapes do not match')
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
  print(); print('  creating mask...')
  lat = data_ds['lat'].isel(time=0,missing_dims='ignore',drop=True)
  lon = data_ds['lon'].isel(time=0,missing_dims='ignore',drop=True)
  mask = xr.DataArray( np.ones([len(ncol_coord)],dtype=bool), coords=[('ncol', ncol_coord.values)], dims='ncol' )
  mask = mask & (lat>=mlat1) & (lat<=mlat2)
  mask = mask & (lon>=mlon1) & (lon<=mlon2)
  #-----------------------------------------------------------------------------
  # apply the mask
  print(); print('  applying mask...')
  # data       = data   .where(mask,drop=True)
  data_2D    = data_2D.where(mask,drop=True)
  center_lon = center_lon.where(mask,drop=True).values
  center_lat = center_lat.where(mask,drop=True).values
  corner_lon = corner_lon.where(mask,drop=True).values
  corner_lat = corner_lat.where(mask,drop=True).values
  #-----------------------------------------------------------------------------
  ncol = len(corner_lat)
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
  print(); print('  starting interpolation...')
  data_interp = interpolate_to_path(npts,nlev,ncll,data.values,path_lat,path_lon,center_lat,center_lon)
  #-----------------------------------------------------------------------------
  # write to file
  print(); print(f'  writing file...  => {tmp_file}')
  ds_tmp = xr.Dataset()
  coords = {'path_coord':path_coord,'lev':lev_coord}
  ds_tmp['data_interp'] = xr.DataArray(data_interp,coords=coords)
  ds_tmp['data_2D']     = data_2D
  ds_tmp['center_lon']  = center_lon
  ds_tmp['center_lat']  = center_lat
  ds_tmp['corner_lon']  = (['ncol','ncorner'],corner_lon)
  ds_tmp['corner_lat']  = (['ncol','ncorner'],corner_lat)
  ds_tmp['path_lon']    = path_lon
  ds_tmp['path_lat']    = path_lat
  ds_tmp.to_netcdf(path=tmp_file,mode='w')
else:
  print(); print(f'  reading file...  => {tmp_file}')
  ds_tmp = xr.open_dataset(tmp_file)
  # print() ; print(ds_tmp) ; print() ; exit()
  # path_coord  = ds_tmp['path_coord'].values
  # path_lon    = ds_tmp['path_lon'].values
  # path_lat    = ds_tmp['path_lat'].values
  lev_coord   = ds_tmp['lev'].values
  data_interp = ds_tmp['data_interp'].values
  data_2D     = ds_tmp['data_2D']#.values
  center_lon  = ds_tmp['center_lon'].values
  center_lat  = ds_tmp['center_lat'].values
  corner_lon  = ds_tmp['corner_lon'].values
  corner_lat  = ds_tmp['corner_lat'].values

if isinstance(path_coord, xr.DataArray): path_coord = path_coord.values
if isinstance( lev_coord, xr.DataArray):  lev_coord = lev_coord.values

# stride = 6
# path_lat, path_lon = path_lat[::stride], path_lon[::stride]
# path_coord, data_interp = path_coord[::stride], data_interp[::stride,:]


buffer_frac = 0.3
n1,n2 = int(npts*buffer_frac), int(npts*(1-buffer_frac))
# path_lat, path_lon = path_lat[n1:n2], path_lon[n1:n2]
path_coord, data_interp = path_coord[n1:n2], data_interp[n1:n2,:]



#---------------------------------------------------------------------------------------------------  
path_coord = path_coord * 1e3

# print(); print(path_coord)
# print(); print(np.min(path_coord))
# print(); print(np.max(path_coord))
# print()
# exit()

#---------------------------------------------------------------------------------------------------
print(); print('  creating plot...')

# create map plot
map_res.cnFillMode    = 'CellFill'
map_res.sfXArray      = center_lon
map_res.sfYArray      = center_lat
map_res.sfXCellBounds = corner_lon
map_res.sfYCellBounds = corner_lat

plot[0] = ngl.contour_map(wks, data_2D.values, map_res)

# overlay path center markers
mrk_res.xyMarker        = 14
mrk_res.xyMarkerThicknessF = 2.
mrk_res.xyMarkerSizeF   = 0.02
mrk_res.xyMarkerColor   = 'red'
mrk_plot = ngl.xy(wks, np.array([xlon,xlon]), np.array([xlat,xlat]), mrk_res)
ngl.overlay( plot[0], mrk_plot )

# overlay path end markers
mrk_res.xyMarker        = 14
mrk_res.xyMarkerThicknessF = 2.
mrk_res.xyMarkerSizeF   = 0.02
mrk_res.xyMarkerColor   = 'red'
# mrk_plot = ngl.xy(wks, np.array([elon1,elon2]), np.array([elat1,elat2]), mrk_res)
mrk_res.xyMarker        = 15
mrk_plot = ngl.xy(wks, np.array([slon,slon]), np.array([slat,slat]), mrk_res)
ngl.overlay( plot[0], mrk_plot )

# overlay path markers
mrk_res.xyMarker        = 16
mrk_res.xyMarkerSizeF   = 0.0005 # 0.001
mrk_res.xyMarkerColor   = 'red'
pth_plot = ngl.xy(wks, path_lon, path_lat, mrk_res)
ngl.overlay( plot[0], pth_plot )

hs.set_subtitles(wks, plot[0], center_string=f'{var_2D}', font_height=0.015)


#-------------------------------------------------------------------------------
res = hs.res_contour_fill()
res.vpHeightF = 0.3
res.tmYLLabelFontHeightF    = 0.008
res.tmXBLabelFontHeightF    = 0.008
res.lbLabelFontHeightF      = 0.01
# res.lbLabelBarOn            = False # True
res.trYReverse              = True
res.cnFillMode              = 'CellFill'
res.sfXArray                = path_coord
res.sfYArray                = lev_coord
res.cnFillPalette           = np.array( cmocean.cm.diff(np.linspace(0,1,256)) )

# res.trYMinF                 = 200
# res.trXMinF                 = np.max(path_coord) * 0.4
# res.trXMaxF                 = np.max(path_coord) * 0.6

# (cmin,cmax,cint) = ngl.nice_cntr_levels(np.min(data_interp), np.max(data_interp), \
#                     cint=None, max_steps=21, returnLevels=False, aboutZero=True )
# res.cnLevels = np.linspace(cmin,cmax,num=21)
res.cnLevelSelectionMode = 'ExplicitLevels'
res.cnLevels = np.linspace(-3.8,3.8,num=31)

plot[1] = ngl.contour(wks, np.transpose(data_interp), res)

hs.set_subtitles(wks, plot[1], center_string=f'{var_3D} (N={ncll})', font_height=0.015)

#---------------------------------------------------------------------------------------------------
num_plot_col = 2#len(plot)
layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent = 5

# ngl.panel(wks,plot,layout,pnl_res)

ngl.panel(wks,[plot[1]],[1,1],pnl_res)

ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------