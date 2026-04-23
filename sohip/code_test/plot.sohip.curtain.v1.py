import os, ngl, subprocess as sp, numpy as np, xarray as xr, dask, copy, string, cmocean, numba
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
# import pg_checkerboard_utilities as pg
from sknni import SkNNI
from pyproj import Geod
geod_obj = Geod(ellps='clrk66') # Use Clarke 1866 ellipsoid.
deg_to_rad = np.pi/180. ; rad_to_deg = 180./np.pi
#---------------------------------------------------------------------------------------------------

data_case = f'helene.ne1024pg2_ne1024pg2.F2010-SCREAMv1.may14_helene_pert1_001'
data_root = f'/pscratch/sd/e/ebercosh/e3sm_scratch/pm-gpu/may14_helene_pert1_001/{data_case}/run'
data_file = f'{data_root}/scream_output.helene.3hourlyINST.ne30.h.INSTANT.nhours_x3.2024-09-27-00000.nc'

# clat1,clon1 = 40.7128, -74.0060 # new york
# clat2,clon2 = 51.5072, 0. # london

# # across italy
# clat1,clat2 = 40,45
# clon1,clon2 = 10,20

# Hurricane Katrina
clat1,clat2 = 15,25
clon1,clon2 = 360-80,360-70

path_spc_km = 10

#---------------------------------------------------------------------------------------------------
fig_file,fig_type = os.getenv('HOME')+'/Research/E3SM/figs_sohip/sohip.curtain.v1','png'
#---------------------------------------------------------------------------------------------------
# load scrip grid data
# scrip_ds = xr.open_dataset('/global/cfs/projectdirs/m3312/whannah/HICCUP/files_grid/scrip_ne30pg2.nc')
scrip_ds = xr.open_dataset('/global/cfs/projectdirs/m3312/whannah/HICCUP/files_grid/scrip_ne256pg2.nc')
# scrip_ds = xr.open_dataset('/pscratch/sd/w/whannah/e3sm_scratch/files_grid/ne1024pg2_scrip.nc')
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
def find_neighbors(center_lat, center_lon, corner_lat, corner_lon, neighbors, edge_flag):
   '''
   Use FV cell information from a SCRIP file to identify adjacent cells
   center_lat[num_col,num_corner]    - input latitude at cell centers
   center_lon[num_col,num_corner]    - input longitude at cell centers
   corner_lat[num_col,num_corner]    - input latitude at cell corners
   corner_lon[num_col,num_corner]    - input longitude at cell corners
   neighbors [num_col,num_neighbors] - output column indices of neighbors                       
   edge_flag [num_col,num_neighbors] - output flag indicating if neighbor shares edge or corner
   '''
   max_num_neighbors = 8
   ncol = len(center_lat)
   # Exact equality doesn't always work so we need to adjust the coordinates a little
   for n in range(ncol):
      for c in range(4):
         corner_lat[n,c] = np.round(corner_lat[n,c],max_num_neighbors)
         corner_lon[n,c] = np.round(corner_lon[n,c],max_num_neighbors)
         # also make sure all prime meridian points have lon=0 rather than 360
         if corner_lon[n,c]==360: corner_lon[n,c] = 0
   # loop through all points that will define the center of the neighborhood
   for n in range(ncol):
      neighbors[n,0] = n
      cnt = 1
      # loop over all 4 corners of the FV cell
      for c in range(4):
         # loop over all points again to find points that share each corner
         for nn in range(ncol):
            # make sure this potential neighbor has not been found already
            if nn not in neighbors[n,0:cnt] :
               # check if corner locations match 
               if  corner_lat[n,c] in corner_lat[nn,:] \
               and corner_lon[n,c] in corner_lon[nn,:] :
                  neighbors[n,cnt] = nn
                  # count how many cell corners are shared
                  common_corner_cnt = 0
                  for cc in range(4):
                     if (corner_lat[n,cc] in corner_lat[nn,:]) and \
                        (corner_lon[n,cc] in corner_lon[nn,:]):
                        common_corner_cnt += 1
                  # if 2 corners are shared then this is an edge neighbor
                  edge_flag[n,cnt] = 1 if common_corner_cnt==2 else 0
                  # increment the neighbor count
                  cnt += 1
         # stop after finding 8 neighbors
         if cnt==9: break
   return
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
# calculate mid-point for map view
Bx = np.cos(np.deg2rad(clat2)) * np.cos(np.deg2rad(clon2 - clon1))
By = np.cos(np.deg2rad(clat2)) * np.sin(np.deg2rad(clon2 - clon1))
lat_mid = np.arctan2(np.sin(np.deg2rad(clat1)) + np.sin(np.deg2rad(clat2)), np.sqrt((np.cos(np.deg2rad(clat1)) + Bx)**2 + By**2))
lon_mid = np.deg2rad(clon1) + np.arctan2(By, np.cos(np.deg2rad(clat1)) + Bx)
#---------------------------------------------------------------------------------------------------
# define path - set point spacing for consistent distance
dist_rd = calc_great_circle_distance(clat1,clat2,clon1,clon2)
dist_km = dist_rd*180/np.pi*111
npts = int(dist_km/path_spc_km)
path = geod_obj.npts(lon1=clon1, lat1=clat1, lon2=clon2, lat2=clat2, npts=npts, radians=False)
# # Add the start and end points
# path = [(start_lon, start_lat)] + [(lon, lat) for lon, lat, _ in points] + [(end_lon, end_lat)]
path_lon = [p[0] for p in path]
path_lat = [p[1] for p in path]
#---------------------------------------------------------------------------------------------------
# identify neighbors
center_lon = scrip_ds['grid_center_lon'][:].values
center_lat = scrip_ds['grid_center_lat'][:].values
corner_lon = scrip_ds['grid_corner_lon'][:,:].values
corner_lat = scrip_ds['grid_corner_lat'][:,:].values
ncol = len(center_lat)
neighbors = np.full([ncol,8+1],-1)
edge_flag = np.full([ncol,8+1],-1)
find_neighbors(center_lat, center_lon, corner_lat, corner_lon, neighbors, edge_flag)
#---------------------------------------------------------------------------------------------------
# update cell corner location for plotting outlines
corner_lon_alt = np.empty([ncol,5])
corner_lat_alt = np.empty([ncol,5])
corner_lon_alt[:,:4] = corner_lon
corner_lat_alt[:,:4] = corner_lat
corner_lon_alt[:,4] = corner_lon_alt[:,0]
corner_lat_alt[:,4] = corner_lat_alt[:,0]
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
map_res.mpDataBaseVersion       = 'MediumRes' # LowRes / MediumRes / HighRes
# map_res.mpProjection            = 'Orthographic' # Robinson / Mollweide / Orthographic
map_res.cnFillPalette           = 'MPL_viridis'
map_res.cnFillMode              = 'CellFill'
map_res.cnCellFillEdgeColor     = 'black'
map_res.cnFillOpacityF          = 0.0
# map_res.mpCenterLatF            = lat_mid
# map_res.mpCenterLonF            = lon_mid
map_res.mpCenterLonF            = lon_mid
map_res.mpLimitMode             = 'LatLon'
map_res.mpMinLatF               = np.min([clat1,clat2])-3
map_res.mpMaxLatF               = np.max([clat1,clat2])+3
map_res.mpMinLonF               = np.min([clon1,clon2])-3
map_res.mpMaxLonF               = np.max([clon1,clon2])+3


# marker resources
mrk_res = hs.res_xy()
mrk_res.xyMarkLineMode  = 'Markers'

lin_res = hs.res_xy()
#---------------------------------------------------------------------------------------------------

# load SCREAM data for plot of clouds
# data_ds = xr.open_dataset( data_file_path ).isel(time=0)

#-------------------------------------------------------------------------------
# create map plot
map_res.cnFillMode    = 'CellFill'
map_res.sfXArray      = center_lon
map_res.sfYArray      = center_lat
map_res.sfXCellBounds = corner_lon
map_res.sfYCellBounds = corner_lat

area = scrip_ds['grid_area'].values

map_plot = ngl.contour_map(wks, area, map_res)

plot.append(map_plot)

#-------------------------------------------------------------------------------
# overlay path end location markers
mrk_res.xyMarker        = 14
mrk_res.xyMarkerThicknessF = 2.
mrk_res.xyMarkerSizeF   = 0.01
mrk_res.xyMarkerColor   = 'red'
mrk_plot = ngl.xy(wks, np.array([clon1,clon2]), np.array([clat1,clat2]), mrk_res)
ngl.overlay( plot[0], mrk_plot )

# overlay path markers
mrk_res.xyMarker        = 16
mrk_res.xyMarkerSizeF   = 0.001
mrk_res.xyMarkerColor   = 'red'
pth_plot = ngl.xy(wks, path_lon, path_lat, mrk_res)
ngl.overlay( plot[0], pth_plot )

#-------------------------------------------------------------------------------
# # find all cells intersected by the path
# idx_list = []
# for n in range(npts):
#   idx = find_closest_cells(path_lat[n],path_lon[n],center_lat,center_lon)
#   if idx not in idx_list: idx_list.append(idx)
# # outline intersected cells
# lin_res.xyLineThicknessF = 8
# lin_res.xyLineColor = 'blue'
# for idx in idx_list:
#   ngl.overlay( plot[0], ngl.xy(wks, corner_lon_alt[idx,:], corner_lat_alt[idx,:], lin_res) )
#-------------------------------------------------------------------------------
# # define test point at mid-point of path
# pp = int(npts/2)
# xx = np.array([1,1])
# # add special marker for test point
# mrk_res.xyMarker        = 3
# mrk_res.xyMarkerSizeF   = 0.01
# mrk_res.xyMarkerColor   = 'blue'
# pth_plot = ngl.xy(wks, xx*path_lon[pp], xx*path_lat[pp], mrk_res)
# ngl.overlay( plot[0], pth_plot )
# # outline cell around test point
# idx = find_closest_cells(path_lat[pp],path_lon[pp],center_lat,center_lon)
# lin_res.xyLineThicknessF = 12
# lin_res.xyLineColor = 'blue'
# nnb_plot = ngl.xy(wks, corner_lon_alt[idx,:], corner_lat_alt[idx,:], lin_res)
# ngl.overlay( plot[0], nnb_plot )
# # outline adjacent neighborhood cells for test point
# lin_res.xyLineThicknessF = 12
# lin_res.xyLineColor = 'magenta'
# lin_res.xyDashPattern = 2
# for n in range(8):
#   nidx = neighbors[idx,1+n]
#   if nidx!=-1:
#     nx = corner_lon_alt[nidx,:]
#     ny = corner_lat_alt[nidx,:]
#     ngl.overlay( plot[0], ngl.xy(wks, nx, ny, lin_res) )
#-------------------------------------------------------------------------------
# find nearest cells to path points
idx_list = []
for n in range(npts):
  tmp_idx_list = find_closest_cells(path_lat[n],path_lon[n],center_lat,center_lon,num_cells=4)
  for idx in tmp_idx_list:
    if idx not in idx_list: idx_list.append(idx)
# outline intersected cells
lin_res.xyLineThicknessF = 6
lin_res.xyLineColor = 'blue'
for idx in idx_list:
  ngl.overlay( plot[0], ngl.xy(wks, corner_lon_alt[idx,:], corner_lat_alt[idx,:], lin_res) )
#-------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
num_plot_col = 1
layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent = 5

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------