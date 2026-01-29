# v3 is modified to use data from v2 and add complimentary map view
#-------------------------------------------------------------------------------
import os, ngl, copy, xarray as xr, numpy as np, glob, dask, numba, cmocean
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
data_dir,data_sub = None,None
#-------------------------------------------------------------------------------
case_name,case,case_dir,case_sub = [],[],[],[]
clr,dsh,mrk = [],[],[]
obs_flag = []
def add_case(case_in,n='',p=None,s='',g=None,d=0,c='black',m=0,r=False,obs=False):
   global case_name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); case_name.append(n)
   case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
   obs_flag.append(obs)
#-------------------------------------------------------------------------------
var = []
file_type_list = []
obs_var_list = []
obs_file_type_list = []
def add_var(var_name,file_type,obs_var=None,obs_file_type=None): 
   var.append(var_name)
   file_type_list.append(file_type)
   obs_var_list.append(obs_var)
   obs_file_type_list.append(obs_file_type)
#-------------------------------------------------------------------------------
scrip_file_sim = '/pscratch/sd/w/whannah/e3sm_scratch/files_grid/ne1024pg2_scrip.nc'
scrip_file_obs = '/pscratch/sd/w/whannah/e3sm_scratch/files_grid/IMERG_1800x3600_scrip.nc'
data_root = '/global/cfs/cdirs/e3sm/gsharing/EAMxx'
obs_data_root = '/pscratch/sd/w/whannah/Obs'

# add_case('DYAMOND2_SCREAMv1' ,n='SCREAMv1 Jan 2020',p=data_root)
# add_case('Apr1_2013_SCREAMv1',n='SCREAMv1 Apr 2013',p=data_root)
# add_case('DYAMOND1_SCREAMv1' ,n='SCREAMv1 Aug 2016',p=data_root)
# add_case('Oct1_2013_SCREAMv1',n='SCREAMv1 Oct 2013',p=data_root)

# add_case('Apr1_2013_SCREAMv1',n='SCREAMv1 Apr 2013',p=data_root)
add_case('IMERG',n='IMERG Apr 2013',p=obs_data_root,s='hourly_nc',obs=True)

#-------------------------------------------------------------------------------

# add_var('LW_flux_up@tom','output.scream.TOMVars.INSTANT')

# add_var('VapWaterPath','output.scream.VertIntegrals.INSTANT')
# add_var('IceWaterPath','output.scream.VertIntegrals.INSTANT')

add_var('precip','output.scream.SurfVars.INSTANT','precipitation','3B-HHR.MS.MRG.3IMERG')

# add_var('precip_liq_surf_mass','output.scream.SurfVars.INSTANT','precipAvg','3B-DAY.MS.MRG.3IMERG.2016')
# add_var('precip_ice_surf_mass','output.scream.SurfVars.INSTANT')
# add_var('surf_evap','output.scream.SurfVars.INSTANT')
# add_var('horiz_winds@bot','output.scream.SurfVars.INSTANT')

#-------------------------------------------------------------------------------
# tmp_data_path = os.getenv('HOME')+'/Research/E3SM/pub_figs/2023_screamv1_4season/data'
tmp_data_path = 'data'

fig_file,fig_type = 'figs/hovmoller.v3.CONUS','png'


first_day,num_days =  0,40

day_center,day_span = 1,0.5

# lat1,lat2 = 30,48 # Carbone & Tuttle 2008
# lat1,lat2 = 32,45 # Surcel et al. 2010
# lat1,lat2 = 38,49 # Clark et al. 2007 (Fig 7)
lat1,lat2 = 35,45
lon1,lon2 = 360-120,360-78
# lon1,lon2 = 360-114,360-78

# dlon = 1; lon_bins = np.arange(50,180,dlon)
# dlon = 1; lon_bins = np.arange(1,360,dlon)
# dlon = 0.5; lon_bins = np.arange(0.5,359.5,dlon)
dlon = 0.1; lon_bins = np.arange(lon1+dlon/2,lon2-dlon/2,dlon)

var_x_case = False

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case = len(case)
num_var  = len(var)

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

plot = [None]*num_case*num_var*2
res = hs.res_contour_fill()
# res.vpHeightF = 0.6
# res.vpWidthF  = 0.3
# res.tmYLLabelFontHeightF         = 0.008
# res.tmXBLabelFontHeightF         = 0.008
res.tiYAxisString = 'Time [days]'
res.tiXAxisString = 'Longitude'
# res.cnFillPalette = 'MPL_viridis'
res.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
# res.cnFillPalette = np.array( cmocean.cm.amp(np.linspace(0,1,256)) )

map_res = hs.res_contour_fill()
# map_res.vpWidthF   = 0.2
# map_res.vpHeightF  = 0.2
map_res.mpLimitMode           = 'LatLon'
map_res.mpMinLatF             = lat1-10
map_res.mpMaxLatF             = lat2+10
map_res.mpMinLonF             = lon1-20
map_res.mpMaxLonF             = lon2+20
map_res.tmYLLabelFontHeightF  = 0.008
map_res.tmXBLabelFontHeightF  = 0.008
map_res.lbLabelFontHeightF    = 0.01
map_res.tmXBOn                = False
map_res.tmYLOn                = False
map_res.mpGridAndLimbOn       = False
map_res.cnFillPalette         = res.cnFillPalette
map_res.lbLeftMarginF         = 0.5
map_res.lbRightMarginF        = 0.5

lres = hs.res_xy()
lres.xyDashPattern    = 1
lres.xyLineThicknessF = 1
lres.xyLineColor      = 'black'

pgres = ngl.Resources()
pgres.nglDraw,pgres.nglFrame = False,False
pgres.gsLineColor = 'red'
pgres.gsLineThicknessF = 6

if 'lev' not in locals(): lev = np.array([0])

pdum = [None]*len(plot) # for map box

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

sim_grid_ds = xr.open_dataset( scrip_file_sim ).rename({'grid_size':'ncol'})
# sim_grid_ds['grid_center_lon'] = xr.where(sim_grid_ds['grid_center_lon']<0,sim_grid_ds['grid_center_lon']+360,sim_grid_ds['grid_center_lon'])
# sim_grid_ds['grid_corner_lon'] = xr.where(sim_grid_ds['grid_corner_lon']<0,sim_grid_ds['grid_corner_lon']+360,sim_grid_ds['grid_corner_lon'])

obs_grid_ds = xr.open_dataset( scrip_file_obs ).rename({'grid_size':'ncol'})
obs_grid_ds['grid_center_lon'] = xr.where(obs_grid_ds['grid_center_lon']<0,obs_grid_ds['grid_center_lon']+360,obs_grid_ds['grid_center_lon'])
obs_grid_ds['grid_corner_lon'] = xr.where(obs_grid_ds['grid_corner_lon']<0,obs_grid_ds['grid_corner_lon']+360,obs_grid_ds['grid_corner_lon'])

#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   time_list = []
   lon_list = []
   hov_list = []
   map_list = []
   for c in range(num_case):
      if obs_flag[c]:
         tfile_type = obs_file_type_list[v]
         tvar = obs_var_list[v]
      else:
         tfile_type = file_type_list[v]
         tvar = var[v]
      if c==0: print(' '*2+'var: '+hc.tcolor.GREEN+tvar+hc.tcolor.ENDC)
      print('\n'+' '*4+'case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)
      #-------------------------------------------------------------------------
      # Load precalculated hovmoller data from v2
      tmp_file = f'{tmp_data_path}/CONUS.hovmoller.horz.v1.tmp.{tvar}.{case[c]}'
      if 'lat1' in globals(): tmp_file += f'.lat1_{lat1}.lat2_{lat2}'
      if 'lon1' in globals(): tmp_file += f'.lon1_{lon1}.lon2_{lon2}'
      tmp_file += '.nc'
      print(' '*6+f'Reading hov data from file: {hc.tcolor.MAGENTA}{tmp_file}{hc.tcolor.ENDC}')
      tmp_ds = xr.open_dataset( tmp_file )
      hov    = tmp_ds[tvar]
      #-------------------------------------------------------------------------
      hov['time'] = hov['time'] - hov['time'][0]
      hov = hov.sel(time=slice(day_center-day_span, day_center+day_span),drop=True)

      lon_bins = hov['lon'].values
      time     = hov['time'].values

      time_list.append(time)
      lon_list.append(lon_bins)
      hov_list.append(hov.values)

      hc.print_stat(hov,name=tvar,compact=True,indent=' '*6)

      #-------------------------------------------------------------------------
      # idenfity the file to load for map plot
      file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/{tfile_type}*'
      file_list = sorted(glob.glob(file_path))
      f_ind = day_center if not obs_flag[c] else day_center*48
      #-------------------------------------------------------------------------
      # load data for map
      print(' '*6+f'Reading map data from file: {hc.tcolor.LIGHTRED}{file_list[f_ind]}{hc.tcolor.ENDC}')
      group_name = 'Grid' if case[c]=='IMERG' else None
      data_ds = xr.open_dataset( file_list[f_ind], group=group_name)
      if tvar=='precip':
         data = data_ds['precip_liq_surf_mass'] + data_ds['precip_ice_surf_mass']
      else:
         data = data_ds[tvar]
      data = data.isel(time=0)
      if obs_flag[c]:
         data['lon'] = xr.where(data['lon']<0,data['lon']+360,data['lon'])
         data = data.stack(ncol=('lat','lon'))
      #-------------------------------------------------------------------
      # convert units
      if var[v]=='precip':
         if obs_flag[c]:
            data = data*24. # mm/hr to mm/day
         else:
            data = (data/100)*86400 # kg/m2 to mm/day using dtime=100 (ne1024)
      #-------------------------------------------------------------------
      hc.print_stat(data,name=tvar,compact=True,indent=' '*6)

      map_list.append(data.values)
      #-------------------------------------------------------------------------

   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   # unit_str = ''
   # if var in ['PRECT','PRECC','PRECL']   : unit_str = '[mm/day]'
   # res.tiXAxisString = unit_str
   # var_str = var
   # if var=='PRECT' : var_str = "Precipitation"
   # res.trXMinF = bin_ds['bins'].min().values
   # res.trXMaxF = bin_ds['bins'].max().values

   #-------------------------------------------------------------------------
   # Create hovmoller plot
   #-------------------------------------------------------------------------
   # min_val = np.min([np.min(h) for h in hov_list]) #; print(f'min_val: {min_val}')
   # max_val = np.max([np.max(h) for h in hov_list]) #; print(f'max_val: {max_val}')
   # num_clev,aboutZero = 21, False
   # clev_tup = ngl.nice_cntr_levels(min_val, max_val, \
   #                                 cint=None, max_steps=num_clev, \
   #                                 returnLevels=False, aboutZero=aboutZero )
   # if clev_tup==None: 
   #    print('SWITCHING TO AUTOMATIC CONTOUR LEVELS!')
   #    res.cnLevelSelectionMode = "AutomaticLevels"   
   # else:
   #    cmin,cmax,cint = clev_tup
   #    res.cnLevels = np.linspace(cmin,cmax,num=num_clev)
   #    res.cnLevelSelectionMode = 'ExplicitLevels'
   if var[v]=='precip': 
      res.cnLevels = np.logspace( -1, 1.35, num=20).round(decimals=4)
      res.cnLevelSelectionMode = 'ExplicitLevels'
   #-------------------------------------------------------------------------
   for c in range(num_case):
      ip = v*num_case+c if var_x_case else c*num_var*2+v
      tres = copy.deepcopy(res)
      # tres.tmXBMode      = 'Explicit'
      # tres.tmXBValues    = np.arange( len(lon_bins) )
      # tres.tmXBLabels    = lon_bins
      if var[v]=='precip': tres.lbTitleString = 'mm/day'
      tres.sfYArray = time_list[c]
      tres.sfXArray = lon_list[c]
      plot[ip] = ngl.contour(wks, hov_list[c] ,tres) 
      #-------------------------------------------------------------------------
      # add line to show time of map
      ngl.overlay( plot[ip], ngl.xy(wks, np.array([-1e3,1e8]), day_center*np.array([1,1]), lres) )
      #-------------------------------------------------------------------------
      hs.set_subtitles(wks, plot[ip], case_name[c], '', tvar, font_height=0.008)

   #----------------------------------------------------------------------------
   # Create map plot
   #----------------------------------------------------------------------------
   # min_val = np.min([np.min(h) for h in map_list]) #; print(f'min_val: {min_val}')
   # max_val = np.max([np.max(h) for h in map_list]) #; print(f'max_val: {max_val}')
   # num_clev,aboutZero = 21, False
   # clev_tup = ngl.nice_cntr_levels(min_val, max_val, \
   #                                 cint=None, max_steps=num_clev, \
   #                                 returnLevels=False, aboutZero=aboutZero )
   # if clev_tup==None: 
   #    print('SWITCHING TO AUTOMATIC CONTOUR LEVELS!')
   #    map_res.cnLevelSelectionMode = 'AutomaticLevels'
   # else:
   #    cmin,cmax,cint = clev_tup
   #    map_res.cnLevels = np.linspace(cmin,cmax,num=num_clev)
   #    map_res.cnLevelSelectionMode = 'ExplicitLevels'
   if var[v]=='precip': 
      # map_res.cnLevels = np.logspace( -3, 0.5, num=30).round(decimals=4)
      map_res.cnLevels = np.logspace( -2, 2, num=20).round(decimals=4)
      # map_res.cnLevels = res.cnLevels
      map_res.cnLevelSelectionMode = 'ExplicitLevels'
   #----------------------------------------------------------------------------
   for c in range(num_case):
      tres = copy.deepcopy(map_res)
      grid_ds = sim_grid_ds if not obs_flag[c] else obs_grid_ds
      tres.cnFillMode    = 'CellFill'
      tres.sfXArray      = grid_ds['grid_center_lon'].values
      tres.sfYArray      = grid_ds['grid_center_lat'].values
      tres.sfXCellBounds = grid_ds['grid_corner_lon'].values
      tres.sfYCellBounds = grid_ds['grid_corner_lat'].values
      if var[v]=='precip': tres.lbTitleString = 'mm/day'
      #-------------------------------------------------------------------------
      ip = (v+1)*num_case+c if var_x_case else c*num_var*2+(v+1)
      plot[ip] = ngl.contour_map(wks,map_list[c],tres)
      #-------------------------------------------------------------------------
      # draw box around map region used for inset
      bx = np.array([lon1,lon2,lon2,lon1,lon1])
      by = np.array([lat1,lat1,lat2,lat2,lat1])
      pdum[ip] = ngl.add_polyline(wks,plot[ip],bx,by,pgres)
      #-------------------------------------------------------------------------
      hs.set_subtitles(wks, plot[ip], case_name[c], '', tvar, font_height=0.008)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
# layout = [1,len(plot)]
layout = [num_var*2,num_case] if var_x_case else [num_case,num_var*2]
ngl.panel(wks,plot[0:len(plot)],layout,hs.setres_panel())
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
