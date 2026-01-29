import os, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
import pg_checkerboard_utilities as pg
import cmocean
data_dir,data_sub = None,None
case,name,clr,dsh = [],[],[],[]
def add_case(case_in,n='',c='black',d=0):
   global name,case
   case.append(case_in); name.append(n); clr.append(c); dsh.append(d)
#-------------------------------------------------------------------------------
add_case('E3SM.CHX.ne30pg2_ne30pg2.F-MMF1.00',            n='E3SM-MMF')
# add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.00',         n='E3SM')
# add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.no_dcape.00',n='E3SM (no DCAPE)')
#-------------------------------------------------------------------------------

var = 'TGCLDLWP'

fig_type = 'png'
fig_file = 'figs/F04-pattern-examples'

lat1,lat2 = -40,40
lon1,lon2 = 90,260

htype,years,months,first_file,num_files = 'h1',[],[],5*12,1

use_remap,remap_grid = False,'180x360' # 90x180 / 180x360

recalc_patterns = False

var_x_case = True

num_plot_col = 2

common_colorbar = True

set_list = []
set_list.append('[0 1 0 1 0 1 0 1]')
set_list.append('[0 0 0 0 1 1 1 1]')
set_list.append('[0 0 0 1 0 1 0 1]')
set_list.append('[0 1 0 1 1 1 1 1]')

#---------------------------------------------------------------------------------------------------
# generate sets of possible neighbor states
#---------------------------------------------------------------------------------------------------

### full set of patterns
rotate_sets = True; sets = pg.all_possible_sets

(num_set,set_len) = sets.shape
set_coord,nn_coord = np.arange(num_set),np.arange(set_len)
sets.assign_coords(coords={'set':set_coord,'neighbors':nn_coord})

set_labels = pg.get_set_labels(sets)

set_idx = []
for s1 in set_list:
   found = False
   for s,s2 in enumerate(set_labels):
      if s1==s2:
         found = True
         set_idx.append(s)
   if not found: exit(f'ERROR: set not found! {s1}')

# print(set_labels)
# print(set_list)
# print(set_idx)
# exit()
#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)

if 'lev' not in vars(): lev = np.array([0])

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix

wks = ngl.open_wks(fig_type,fig_file,wkres)
# plot = [None]*(num_var*num_case)
plot = [None]*len(set_list)
   
res = hs.res_contour_fill_map()

res.tmYLLabelFontHeightF         = 0.008
res.tmXBLabelFontHeightF         = 0.008
res.lbLabelFontHeightF           = 0.018
res.tmXBOn                       = False
res.tmYLOn                       = False
# res.mpGeophysicalLineColor       = 'white'

res.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
res.cnLevelSelectionMode = 'ExplicitLevels'
if var=='PRECT'    : res.cnLevels = np.arange(5,150+5,5)/1e1
# if var=='TGCLDLWP' : res.cnLevels = np.arange(0.01,0.25,0.015)
if var=='TGCLDLWP' : res.cnLevels = np.logspace( -2, 0.25, num=60).round(decimals=2)

if common_colorbar:
   res.lbLabelBarOn = False
else:
   res.lbLabelBarOn = True

var_str = var
if var=='PRECT':     var_str = 'Precipitation'
if var=='TGCLDLWP':  var_str = 'Liquid Water Path'

scrip_file_path = 'scrip_files/ne30pg2_scrip.nc'
scrip_ds = xr.open_dataset(scrip_file_path)
res.cnFillMode    = "CellFill"
res.sfXArray      = scrip_ds['grid_center_lon'].values
res.sfYArray      = scrip_ds['grid_center_lat'].values
res.sfXCellBounds = scrip_ds['grid_corner_lon'].values
res.sfYCellBounds = scrip_ds['grid_corner_lat'].values

### resources for box around pattern of interest
pgres = ngl.Resources()
pgres.nglDraw,pgres.nglFrame = False,False
pgres.gsLineColor = 'red'
pgres.gsLineThicknessF = 12

#---------------------------------------------------------------------------------------------------
# Load SCRIP data
#---------------------------------------------------------------------------------------------------
center_lat = scrip_ds['grid_center_lat'].values
center_lon = scrip_ds['grid_center_lon'].values
corner_lon = scrip_ds['grid_corner_lon'].values
corner_lat = scrip_ds['grid_corner_lat'].values

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_data_file_name(case,var):
   data_file = f'{tmp_file}.daily.{case}.{var}.sml_{subset_min_length}.nc'
   return data_file
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

hc.printline()
print('  var: '+hc.tcolor.MAGENTA+var+hc.tcolor.ENDC)
data_list,area_list,lat_list,lon_list = [],[],[],[]
std_list,cnt_list = [],[]

for c in range(num_case):

   print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)

   case_obj = he.Case( name=case[c], time_freq='daily')
   case_obj.set_coord_names(var)
   use_remap,remap_str = False,'remap_ne30pg2'
   
   #-------------------------------------------------------------------------
   # read the data
   #-------------------------------------------------------------------------
   data = case_obj.load_data(var, component='eam', htype=htype, ps_htype=htype,
                             first_file=first_file, num_files=num_files)
   lat = case_obj.load_data('lat', htype=htype,num_files=1)
   lon = case_obj.load_data('lon', htype=htype,num_files=1)

   # Convert to daily mean
   data = data.resample(time='D').mean(dim='time')

   #-------------------------------------------------------------------------
   # find neighbors and count the occurrence of each set
   #-------------------------------------------------------------------------
   scratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl/chx-detection-data'
   pattern_file=f'{scratch}/chx-examples-ocurrence.{var}.{case}.nc'
   if recalc_patterns:
      scripfile_path = 'scrip_files/ne30pg2_scrip.nc'
      cnt_ds = pg.find_neighbors_and_count_sets(data, scripfile_path, 
                                                sets, rotate_sets, 
                                                keep_time=True, 
                                                verbose=True)
      cnt_ds.to_netcdf(path=pattern_file,mode='w')
   else:
      cnt_ds = xr.open_dataset( pattern_file )

   #-------------------------------------------------------------------------
   # apply spatial mask to limit search area
   #-------------------------------------------------------------------------
   cnt_ds['lat'],cnt_ds['lon'] = lat,lon

   mask = xr.DataArray( np.ones([len(cnt_ds['ncol'])],dtype=bool), 
                        coords=[('ncol', cnt_ds['ncol'])], dims='ncol' )
   mask = mask & (cnt_ds['lat']>=lat1) & (cnt_ds['lat']<=lat2)
   mask = mask & (cnt_ds['lon']>=lon1) & (cnt_ds['lon']<=lon2)

   cnt_ds = cnt_ds.where( mask, drop=True)

   #-------------------------------------------------------------------------
   # Find single instance of each set
   def find_first_set_occurrence(flag,time,ncol):
      for t,tt in enumerate(time):
         for n,nn in enumerate(ncol):
            if flag[t,n]==1:
               return t,nn
      return None,None
               
   ncol_list, time_list, lat_list, lon_list = [],[],[],[]
   cnt = 0
   for idx in set_idx:
      t,n = find_first_set_occurrence(cnt_ds['cnt'].isel(set=idx).values,
                                      cnt_ds['time'].values,
                                      cnt_ds['ncol'].values)
      if n is not None and t is not None: 
         time_list.append(t)
         ncol_list.append(n)
      else:
         exit(f'ERROR: something bad happened! t: {t}  n: {n}')

   #-------------------------------------------------------------------------
   # get local neighbor indices
   #-------------------------------------------------------------------------
   def get_neighbors(ilat,ilon):
      dist = pg.calc_great_circle_distance( ilat, center_lat[:], ilon, center_lon[:] )
      mcol = np.argsort(dist,kind='mergesort')[0]
      num_neighbors = 8
      neighbors = np.full([len(center_lat),num_neighbors+1],-1)
      bearings  = np.full([len(center_lat),num_neighbors+1],-1)
      edge_flag = np.full([len(center_lat),num_neighbors+1],-1)
      pg.find_neighbors(center_lat, center_lon, corner_lat, corner_lon, neighbors, bearings, edge_flag)
      neighborhood,neighborbear = neighbors[mcol,:], bearings[mcol,:]
      neighborbear,neighborhood = np.array(neighborbear), np.array(neighborhood)
      tmp_neighborhood = neighborhood[1:]
      neighborhood[1:] = tmp_neighborhood[ np.argsort(neighborbear[1:]) ]
      return neighborhood
   #-------------------------------------------------------------------------
   # Create plot
   #-------------------------------------------------------------------------
   # ip = v*num_case+c if var_x_case else c*num_var+v
   tres = copy.deepcopy(res)
   
   for i,idx in enumerate(set_idx):
      t,n = time_list[i],ncol_list[i]
      ilat, ilon = lat[n].values, lon[n].values
      dx = 4
      tlat1,tlat2 = ilat-dx,ilat+dx
      tlon1,tlon2 = ilon-dx,ilon+dx
      map_data = data.isel(time=t).values
      tres.mpMinLatF,tres.mpMaxLatF = tlat1,tlat2
      tres.mpMinLonF,tres.mpMaxLonF = tlon1,tlon2

      ip = i

      plot[ip] = ngl.contour_map(wks,map_data,tres) 

      ### draw a box around adjacent neighbors
      neighborhood = get_neighbors(ilat,ilon)
      clat = corner_lat[ neighborhood[[1,3,5,7,1]], [0,3,2,1,0] ]
      clon = corner_lon[ neighborhood[[1,3,5,7,1]], [0,3,2,1,0] ]
      pdum = ngl.add_polyline(wks,plot[ip],clon,clat,pgres)

      ### add annotations that correspond to the anomaly flag
      ttres = ngl.Resources()
      ttres.txFontColor = 'red'
      ttres.txFontHeightF = 0.04
      txid = [None]*9
      for nn in range(1,9):
         nnid = neighborhood[nn]
         ty,tx = center_lat[ nnid ], center_lon[ nnid ]
         if map_data[nnid]<=map_data[neighborhood[0]]: txt_flag = '0'
         if map_data[nnid]> map_data[neighborhood[0]]: txt_flag = '1'
         txid[nn] = ngl.add_text(wks, plot[ip], txt_flag, tx, ty, ttres)

      # hs.set_subtitles(wks, plot[ip], name[c], '', var_str, font_height=0.015)
      hs.set_subtitles(wks, plot[ip], '', set_list[i], '', font_height=0.015)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
   
pnl_res = hs.setres_panel()

### add panel labels
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.015

if common_colorbar: 
   pnl_res.nglPanelLabelBar = True
   pnl_res.lbTitleFontHeightF = 0.015
   pnl_res.nglPanelLabelBarLabelFontHeightF = 0.015
   pnl_res.lbTitlePosition = 'Bottom'
   if var=='PRECT':    pnl_res.lbTitleString = 'mm/day'
   if var=='TGCLDLWP': pnl_res.lbTitleString = 'kg/m2'


pnl_res.nglPanelYWhiteSpacePercent = 5

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
