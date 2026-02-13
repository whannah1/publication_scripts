#-------------------------------------------------------------------------------
import os, ngl, copy, string, xarray as xr, numpy as np, glob, dask, numba, cmocean, subprocess as sp, pandas as pd
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
#-------------------------------------------------------------------------------

fig_file,fig_type = f'figs/fig-MJO-DYNAMO-hovmoller-v1','png'
tmp_file_head = 'data/MJO-DYNAMO-hovmoller-v1'

time_mean_opt    = '1D' # none / 1D / 5D / 10D

lat1,lat2 = -15,15
dlon = 2

recalculate = True

# num_files = 10 # comment to disable and load all files

print_stats = True

use_anomaly = False


#---------------------------------------------------------------------------------------------------

# tmp_file_root, tmp_file_prefix = os.getenv('HOME')+f'/Research/E3SM/data_temp','scream.hov.v1'
def get_tmp_file(case,var):
   return f'{tmp_file_head}.{case}.{var}.{time_mean_opt}.nc'
#---------------------------------------------------------------------------------------------------
case,case_name,case_dir,case_sub = [],[],[],[]
obs_flag = []
def add_case(case_in,n='',p=None,s='',g=None,d=0,c='black',m=0,r=False,obs=False):
   global case_name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); case_name.append(n); case_dir.append(p); case_sub.append(s)
   obs_flag.append(obs)
#---------------------------------------------------------------------------------------------------
var = []
var_str = []
var_unit = []
file_type_list = []
obs_var_list = []
obs_file_type_list = []
lev_list = []
def add_var(var_name, file_type, obs_var=None, obs_file_type=None, name='', unit='', lev=None):
   var.append(var_name)
   file_type_list.append(file_type)
   obs_var_list.append(obs_var)
   obs_file_type_list.append(obs_file_type)
   var_str.append(name)
   var_unit.append(unit)
   lev_list.append(lev)
#---------------------------------------------------------------------------------------------------
# scrip_file_sim = '/pscratch/sd/w/whannah/e3sm_scratch/files_grid/ne1024pg2_scrip.nc'
# scrip_file_obs = '/pscratch/sd/w/whannah/e3sm_scratch/files_grid/IMERG_1800x3600_scrip.nc'
# scrip_file_sim = '/lustre/orion/cli115/proj-shared/hannah6/HICCUP/data/scrip_ne1024pg2.nc'
# scrip_file_sim = '/lustre/orion/cli115/proj-shared/hannah6/HICCUP/data/scrip_ne30pg2.nc'
# scrip_file_obs = '~/HICCUP/data_scratch/scrip_ERA5_721x1440.nc'
hx_data_root = '/global/cfs/cdirs/e3smdata/simulations/scream-decadal/DYNAMO_MJO_hindcasts'
hx_data_sub  = 'data_remap_73x144'

add_case('Obs',obs=False)

### 2025 DYNAMO hindcasts - testing rain frac for decadal paper
mn_min,dy_min,mn_max,dy_max = 11,10,11,30
add_case('SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-10.rfrac_fix_0', n='SCREAMv1 no rain frac fix',       p=hx_data_root,s=hx_data_sub)
add_case('SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-10.rfrac_fix_1', n='SCREAMv1 w/ rain frac fix',       p=hx_data_root,s=hx_data_sub)

# mn_min,dy_min,mn_max,dy_max = 11,15,12, 9
# add_case('SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-15.rfrac_fix_0', n='SCREAMv1 3-km control',           p=hx_data_root,s=hx_data_sub)
# add_case('SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-15.rfrac_fix_1', n='SCREAMv1 3-km w/ cold pool fix',  p=hx_data_root,s=hx_data_sub)

# mn_min,dy_min,mn_max,dy_max = 11,20,12,20
# add_case('SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-20.rfrac_fix_0', n='SCREAMv1 no rain frac fix',       p=hx_data_root,s=hx_data_sub)
# add_case('SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-20.rfrac_fix_1', n='SCREAMv1 w/ rain frac fix',       p=hx_data_root,s=hx_data_sub)

# scrip_file = 'files_grid/ne30pg2_scrip.nc'
# file_type = 'output.scream.2D.1hr.ne30pg2.AVERAGE.nhours_x1'

scrip_file = 'files_grid/73x144_scrip.nc'
file_type = 'output.scream.2D.1hr.ne30pg2.AVERAGE.nhours_x1'

num_plot_col = len(case)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

add_var('precip_total_surf_mass_flux', file_type, name='precip', unit='mm/day',  obs_var='', obs_file_type='daily_QC_Jan_2020')
add_var('LW_flux_up_at_model_top',     file_type, name='OLR',    unit='W/m2',    obs_var='', obs_file_type='')
add_var('U_at_850hPa',                 file_type, name='U850',   unit='m/s',     obs_var='u', obs_file_type='ERA5.daily.atm')

# add_var('VapWaterPath',                file_type, name='VWP')
# add_var('LiqWaterPath',                file_type, name='LWP')

# add_var('surf_sens_flux',                 file_type, name='SHF')
# add_var('surface_upward_latent_heat_flux',file_type, name='LHF')

# add_var('U_at_model_bot',              file_type, name='UBOT')

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_var = len(case),len(var)

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
# wks = ngl.open_wks(fig_type,fig_file,wkres)
   
# plot = [None]*(num_var*num_case)
plot = [None]*(num_case)

res = hs.res_contour_fill()
res.vpHeightF = 0.4
res.tmYLLabelFontHeightF   = 0.015
res.tmXBLabelFontHeightF   = 0.015
res.tiXAxisFontHeightF     = 0.025
res.tiYAxisFontHeightF     = 0.025
res.tiYAxisString          = 'Date'
res.tiXAxisString          = 'Longitude'
res.lbOrientation          = 'Vertical'
res.lbLabelFontHeightF     = 0.015
res.lbLabelBarOn           = False

res.tmYLLabelAngleF = 45


lres = hs.res_xy()
lres.xyDashPattern    = 1
lres.xyLineThicknessF = 2
lres.xyLineColor      = 'black'
#---------------------------------------------------------------------------------------------------
def run_cmd(cmd,verbose=True,indent='  '):
   cmd_str = indent + hc.tcolor.GREEN + cmd + hc.tcolor.ENDC
   if verbose: print('\n'+cmd_str)
   proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, universal_newlines=True)
   (msg, err) = proc.communicate()
   if verbose and msg!='': print(f'  msg: {msg}')
   if err!='' and not verbose: print(cmd_str)
   if err!='': print(f'err: {err}'); exit()
   return msg
#---------------------------------------------------------------------------------------------------
def get_obs_name(case,var):
   obs_name, obs_var = None, None
   if case=='Obs' and var in ['PRECT']:                        obs_var = 'PRECT'; obs_name = 'IMERG'
   if case=='Obs' and var in ['OLR','FLNT','FLUT']:            obs_var = 'olr'  ; obs_name = 'NOAA'
   if case=='Obs' and var in ['U850']:                         obs_var = 'U850' ; obs_name = 'ERA5'
   if case=='Obs' and var in ['precip_total_surf_mass_flux']:  obs_var = 'PRECT'; obs_name = 'IMERG'
   if case=='Obs' and var in ['LW_flux_up_at_model_top']:      obs_var = 'olr'  ; obs_name = 'NOAA'
   if case=='Obs' and var in ['U_at_850hPa']:                  obs_var = 'U850' ; obs_name = 'ERA5'
   

   if obs_name is None: raise ValueError(f'get_obs_name(): obs_name cannot be None!  case: {case}  var: {var}')
   return obs_name, obs_var
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

sim_grid_ds = xr.open_dataset( scrip_file ).rename({'grid_size':'ncol'})

grid_rank = len(sim_grid_ds['grid_rank'])

if grid_rank==1:
   sim_area = sim_grid_ds['grid_area']
   sim_lat  = sim_grid_ds['grid_center_lat']
   sim_lon  = sim_grid_ds['grid_center_lon']

   tmp_data = np.ones(len(sim_area),dtype=bool)
   sim_mask = xr.DataArray( tmp_data, coords=sim_area.coords )
   if 'lat1' in locals(): sim_mask = sim_mask & (sim_lat>=lat1) & (sim_lat<=lat2)
   # if 'lon1' in locals(): sim_mask = sim_mask & (sim_lon>=lon1) & (sim_lon<=lon2)

   sim_area = sim_area.where(sim_mask,drop=True)
   sim_lat  = sim_lat .where(sim_mask,drop=True)
   sim_lon  = sim_lon .where(sim_mask,drop=True)

if grid_rank==2:
   sim_area = sim_grid_ds['grid_area']
   sim_lon  = sim_grid_ds['grid_center_lon']

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# def get_var_data(ds,var,obs_var,obs_flag):
#    tvar = var
#    if var=='horiz_winds_at_model_bot_u': tvar = 'horiz_winds_at_model_bot'
#    if var=='horiz_winds_at_model_bot_v': tvar = 'horiz_winds_at_model_bot'
#    #----------------------------------------------------------------------------
#    if obs_flag:
#       data = ds[obs_var]
#    else:
#       data = ds[var]
#    # if var=='precip':
#    #    data = ds['precip_liq_surf_mass'] + ds['precip_ice_surf_mass']
#    # elif var=='precip_total_surf_mass_flux':
#    #    data = ds['precip_liq_surf_mass_flux'] + ds['precip_ice_surf_mass_flux']
#    # else:
#    #    data = ds[var]
#    #----------------------------------------------------------------------------
#    if obs_flag:
#       if var=='U_at_850hPa': data = data.sel(pressure_level=850)
#       if var=='V_at_850hPa': data = data.sel(pressure_level=850)
#       if var=='U_at_500hPa': data = data.sel(pressure_level=500)
#       if var=='V_at_500hPa': data = data.sel(pressure_level=500)
#    else:
#       if var=='horiz_winds_at_model_bot_u': data = data.isel(dim2=0)
#       if var=='horiz_winds_at_model_bot_v': data = data.isel(dim2=1)
#    #----------------------------------------------------------------------------
#    # convert units
#    if var=='precip_total_surf_mass_flux':
#       data = data*86400*1e3 # m/s to mm/day
#    if var=='precip':
#       if     obs_flag: data = data*24.          # mm/hr to mm/day
#       if not obs_flag: data = (data/100)*86400  # kg/m2 to mm/day using dtime=100 (ne1024)
#    #----------------------------------------------------------------------------
#    if obs_flag:
#       # if 'latitude'  in data.dims: data = data.rename({'latitude':'lat'})
#       # if 'longitude' in data.dims: data = data.rename({'longitude':'lon'})
#       # data = data.stack(ncol=('lat','lon'))
#       data = data.rename({'valid_time':'time'})
#    #----------------------------------------------------------------------------
#    return data
#---------------------------------------------------------------------------------------------------
fig_file_list = []
for v in range(num_var):
   fig_file_var = fig_file+f'.{var[v]}'
   wks_var = ngl.open_wks(fig_type,fig_file_var,wkres)
   fig_file_list.append(fig_file_var)
   #----------------------------------------------------------------------------
   data_list = []
   time_list = []
   lon_list  = []
   for c in range(num_case):
      tcase = case[c]
      if obs_flag[c]:
         tcase,tvar = get_obs_name(case[c],var[v])
         # tfile_type = obs_file_type_list[v]
         # tvar = obs_var_list[v]
      else:
         tfile_type = file_type_list[v]
         tvar = var[v]
      if c==0: print(''+' '*2+'var: '+hc.tcolor.GREEN+tvar+hc.tcolor.ENDC)
      print(''+' '*4+'case: '+hc.tcolor.CYAN+tcase+hc.tcolor.ENDC)
      #-------------------------------------------------------------------------
      tmp_file = get_tmp_file(case[c],var[v])
      if recalculate :
         print(' '*6+'recalculating...')
         #----------------------------------------------------------------------
         # idenfity the files to load
         if case[c]=='Obs':
            tcase,tvar = get_obs_name(case[c],var[v])
            # dst_root = None
            if tcase=='NOAA':  file_list = ['/global/cfs/cdirs/m3312/whannah/obs_data/OLR/olr.day.mean.nc']
            if tcase=='IMERG': file_list = ['/global/cfs/cdirs/m3312/whannah/obs_data/IMERG/IMERG_Daily_PRECT_200101_202012.remap_73x144.nc']
            if tcase=='ERA5':  file_list = ['/global/cfs/cdirs/m3312/whannah/obs_data/ERA5/U850_198001_202212.remap_73x144.nc']
            # if dst_root is None: raise ValueError('dst_root cannot be None!')
         else:
            file_path = f'{case_dir[c]}/{tcase}/{case_sub[c]}/{tfile_type}*'
            file_list = sorted(glob.glob(file_path))
         #----------------------------------------------------------------------
         if 'num_files' in locals(): file_list = file_list[:num_files]
         #----------------------------------------------------------------------
         if file_list==[]: print('ERROR: Empty file list:'); print(); print(file_path); exit()
         #----------------------------------------------------------------------
         print(' '*6+f'Loading data ({tvar})...')
         #----------------------------------------------------------------------
         group_name = None # group_name = 'Grid' if case[c]=='IMERG' else None
         #----------------------------------------------------------------------
         # print()
         # for f in file_list[:2]: print(f)
         # print()
         # ds = xr.open_mfdataset( file_list[:2], group=group_name)
         # print(); print(ds); print()
         # exit()
         #----------------------------------------------------------------------
         ds = xr.open_mfdataset( file_list, group=group_name)
         #----------------------------------------------------------------------
         # Load the data
         # data = get_var_data(ds,tvar,obs_flag[c])
         # data = get_var_data(ds,var[v],obs_var_list[v],obs_flag[c])
         #----------------------------------------------------------------------
         data = ds[tvar]
         #----------------------------------------------------------------------
         if tvar=='precip_total_surf_mass_flux':   data = data*86400*1e3   # m/s to mm/day
         if tcase=='IMERG':                        data = data*86400*1e3   # m/s to mm/day
         # if tcase=='IMERG':                        data = data*24.         # mm/hr to mm/day
         #----------------------------------------------------------------------
         # apply mask
         if grid_rank==1: data = data.where(sim_mask,drop=True)
         if grid_rank==2: data = data.where((data.lat>=lat1) & (data.lat<=lat2),drop=True)
         # if not obs_flag[c]: data = data.where(sim_mask,drop=True)
         # if     obs_flag[c]: data = data.where(obs_mask,drop=True)
         #----------------------------------------------------------------------
         if time_mean_opt!='none': 
            data = data.resample(time=time_mean_opt).mean(dim='time')
         #-------------------------------------------------------------------
         if 'lev' in data.dims and lev_list[v] is not None: 
            if lev_list[v]<0: data = data.isel({'lev':np.absolute(lev_list[v])})
         #-------------------------------------------------------------------------
         # print(); print(data)
         # print(); print(sim_area)
         # print(); print(sim_lon)
         # print()
         #-------------------------------------------------------------------------
         # Calculate hovmoller
         print(' '*6+f'Calculating hovmoller: {hc.tclr.MAGENTA}{tmp_file}{hc.tclr.END}')
         # if     obs_flag[c]: area, lon = obs_area, obs_lon
         # if not obs_flag[c]: area, lon = sim_area, sim_lon
         if grid_rank==1:
            area, lon = sim_area, sim_lon
            hov_ds = hc.bin_YbyX(data, lon, bin_min=lon.min().values, bin_max=lon.max().values,
                                 bin_spc=dlon, wgt=area, keep_time=True ).drop_vars(['bin_std'])
            hov_ds['bin_val'] = hov_ds['bin_val'].transpose('time','bin')
            hov_ds['bin_cnt'] = hov_ds['bin_cnt'].transpose('time','bin')
            hov_ds['time'] = data['time']
            del area; del lon
         if grid_rank==2:
            hov_ds = xr.Dataset()
            hov_ds['data'] = data.mean(dim='lat')
         #----------------------------------------------------------------------
         # write to file
         print(' '*6+f'writing to file: {tmp_file}')
         hov_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         print(' '*6+f'writing from file: {hc.tclr.YELLOW}{tmp_file}{hc.tclr.END}')
         hov_ds = xr.open_dataset( tmp_file )
      #-------------------------------------------------------------------------
      # print some summary stats after spatial averaging
      # if print_stats: hc.print_stat(hov_ds['bin_val'],name=tvar,compact=True,indent=' '*6)
      #-------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------
      if grid_rank==1:
         if use_anomaly: hov_ds['bin_val'] = hov_ds['bin_val'] - hov_ds['bin_val'].mean(dim='time')
         data_list.append( hov_ds['bin_val'].values )
         time_list.append( hov_ds['time'] )
         lon_list.append( hov_ds['bin'].values )
      if grid_rank==2:
         if use_anomaly: hov_ds['data'] = hov_ds['data'] - hov_ds['data'].mean(dim='time')
         data_list.append( hov_ds['data'].values )
         time_list.append( hov_ds['time'] )
         lon_list.append(  hov_ds['lon'].values )
   #----------------------------------------------------------------------------
   # Create plot
   tres = copy.deepcopy(res)
   data_min = np.nanmin([np.nanmin(d) for d in data_list])
   data_max = np.nanmax([np.nanmax(d) for d in data_list])
   #----------------------------------------------------------------------------
   print(f'  data_min / data_max : {data_min} / {data_max}')
   #----------------------------------------------------------------------------
   # Set colors
   tres = copy.deepcopy(res)
   pcp_clr_map = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
   olr_clr_map = np.array( cmocean.cm.dense(np.linspace(0,1,256)) )
   wnd_clr_map = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
   # tres.cnFillPalette = 'MPL_viridis'
   tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
   # if var[v]=='U'                          : tres.cnFillPalette = bal_clr_map
   # if var[v]=='T_2m'                       : tres.cnFillPalette = amp_clr_map
   # if var[v]=='precip_total_surf_mass_flux': tres.cnFillPalette = rain_clr_map
   # if var[v]=='precip_ice_surf_mass'       : tres.cnFillPalette = rain_clr_map
   # if var[v]=='precip_liq_surf_mass'       : tres.cnFillPalette = rain_clr_map
   # if var[v]=='LiqWaterPath'               : tres.cnFillPalette = rain_clr_map
   # if var[v]=='IceWaterPath'               : tres.cnFillPalette = rain_clr_map
   #----------------------------------------------------------------------------
   # Set explicit contour levels
   tres.cnLevelSelectionMode = 'ExplicitLevels'
   # if var[v]=='precip_total_surf_mass_flux': tres.cnLevels = np.logspace( -2, 2, num=20).round(decimals=2)
   # if var[v]=='FLNT' : tres.cnLevels = np.logspace(  0.0, 2.0, num=30).round(decimals=4)
   # if var[v]=='precip_total_surf_mass_flux': tres.cnLevels = np.arange( 4, 24,4)
   # if var[v]=='precip_total_surf_mass_flux': tres.cnLevels = np.arange( -24, 24+4,4)
   # if var[v]=='LW_flux_up_at_model_top':     tres.cnLevels = np.arange( -60, 60+12,12)
   # if var[v]=='U_at_850hPa':                 tres.cnLevels = np.arange( -5, 5+1,1)
   # if var[v]=='LiqWaterPath'               : tres.cnLevels = np.arange( -0.04, 0.04+0.01,0.01)
   # if var[v]=='IceWaterPath'               : tres.cnLevels = np.arange( -8, 8+2,2)

   # if var[v]=='precip_total_surf_mass_flux': dx= 4; n=11; tres.cnLevels = np.linspace( -n*dx/2, n*dx/2, n+1)
   # if var[v]=='LW_flux_up_at_model_top':     dx=10; n=11; tres.cnLevels = np.linspace( -n*dx/2, n*dx/2, n+1)
   # if var[v]=='U_at_850hPa':                 dx= 1; n=11; tres.cnLevels = np.linspace( -n*dx/2, n*dx/2, n+1)

   if use_anomaly: 
      if var[v]=='surf_sens_flux':                 tres.cnFillPalette = wnd_clr_map; tres.cnLevels = np.linspace( -14, 14, 11)
      if var[v]=='surface_upward_latent_heat_flux':tres.cnFillPalette = wnd_clr_map; tres.cnLevels = np.linspace( -50, 50, 11)
   else:
      if var[v]=='precip_total_surf_mass_flux': tres.cnFillPalette = pcp_clr_map; tres.cnLevels = np.linspace( 2, 30, 11)
      if var[v]=='LW_flux_up_at_model_top':     tres.cnFillPalette = olr_clr_map; tres.cnLevels = np.linspace( 180, 300, 11)
      if var[v]=='U_at_850hPa':                 tres.cnFillPalette = wnd_clr_map; tres.cnLevels = np.linspace( -8, 8, 11)

      if var[v]=='surf_sens_flux':                 tres.cnFillPalette = wnd_clr_map; tres.cnLevels = np.linspace( 10, 80, 11)
      if var[v]=='surface_upward_latent_heat_flux':tres.cnFillPalette = wnd_clr_map; tres.cnLevels = np.linspace( 20, 220, 11)



   # tres.cnLevels = np.linspace( data_min, data_max, 11)
   
   #----------------------------------------------------------------------------
   for c in range(num_case):
      # ip = v*num_case+c
      ip = c
      #-------------------------------------------------------------------------
      time = time_list[c]
      time_coord = ( ( time - time[0] ).astype('float') / 86400e9 ).values
      #-------------------------------------------------------------------------
      # if obs_flag[c]:
      #    print(); print(time[0].values)
      #    # time = pd.Timestamp(time)
      #    # time = time.convert_calendar("noleap")
      #    print(); print(pd.Timestamp(time[0].values).year)
      #    # print(); print(time[0].values.convert_calendar("noleap"))
      #    # pd_timestamp = pd.Timestamp(time[0].values)
      #    # print(); print(pd_timestamp.year)
      #    # print(); print(time['month'])
      #    print()
      #    exit()
      #-------------------------------------------------------------------------
      ntime  = len(time.values)

      # if obs_flag[c]:
      #    date_labels = [ f'{pd.Timestamp(t).month}-{pd.Timestamp(t).day}' for t in time.values ]
      # else:
      #    # date_labels = [ t.strftime("%Y-%m-%d") for t in time.values ]
      #    date_labels = [ t.strftime("%m-%d") for t in time.values ]

      if case[c]=='Obs':
         date_labels = [ f'{pd.Timestamp(t).year}-{pd.Timestamp(t).month:02}-{pd.Timestamp(t).day:02}' for t in time.values ]
      else:
         date_labels = [ t.strftime("%Y-%m-%d") for t in time.values ]

      #-------------------------------------------------------------------------
      tmin,tmax = None,None
      for t,tt in enumerate(time.values):
         # if obs_flag[c]:
         if case[c]=='Obs':
            yr = pd.Timestamp(tt).year
            mn = pd.Timestamp(tt).month
            dy = pd.Timestamp(tt).day
         else:
            yr = tt.year
            mn = tt.month
            dy = tt.day
         if yr==2011:
            # print(f'mn-dy: {mn}-{dy}')
            if mn==mn_min and dy==dy_min: tmin = t
            if mn==mn_max and dy==dy_max: tmax = t
            
            # if mn==11 and dy==10: tmin = t
            # if mn==12 and dy== 5: tmax = t

            # if mn==11 and dy==15: tmin = t
            # if mn==12 and dy== 9: tmax = t

            # if mn==11 and dy==20: tmin = t
            # if mn==12 and dy==20: tmax = t



      if tmin is None: raise ValueError('ERROR: time values not found!')
      # tres.trYMinF = tmin
      # tres.trYMaxF = tmax
      tres.tmYLMode = 'Explicit'
      tres.tmYLValues = time_coord  [tmin:tmax:5]
      tres.tmYLLabels = date_labels [tmin:tmax:5]
      tres.sfYArray = time_coord    [tmin:tmax]
      tres.sfXArray = lon_list[c]
      #-------------------------------------------------------------------------
      tres.trXMinF = 0
      tres.trXMaxF = 180
      #-------------------------------------------------------------------------
      # tres.cnFillPalette = "MPL_viridis"
      if not hasattr(tres,'cnLevels') : tres.cnLevelSelectionMode = 'AutomaticLevels'
      # data_list[c] = data_list[c] - np.mean(data_list[c])
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      # if c==num_case-1: tres.lbLabelBarOn = True # only enable colorbar on RHS panels
      plot[ip] = ngl.contour(wks_var, np.ma.masked_invalid(  data_list[c][tmin:tmax,:] ), tres) 
      #-------------------------------------------------------------------------
      tname = case_name[c]
      if case[c]=='Obs': tname,tvar = get_obs_name(case[c],var[v])
      # hs.set_subtitles(wks, plot[ip], tname, '', var_str[v], font_height=0.01)
      # hs.set_subtitles(wks_var, plot[ip], var_str[v], '', tname, font_height=0.009)
      hs.set_subtitles(wks_var, plot[ip], '', tname, '', font_height=0.012)
   #----------------------------------------------------------------------------
   pnl_res = hs.setres_panel()
   pnl_res.nglPanelYWhiteSpacePercent       = 5
   pnl_res.nglPanelXWhiteSpacePercent       = 5
   pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase[v*3:v*3+3+1])
   pnl_res.nglPanelFigureStringsJust        = 'TopLeft'
   pnl_res.nglPanelFigureStringsFontHeightF = 0.01
   pnl_res.nglPanelLabelBar                 = True
   pnl_res.nglPanelLabelBarLabelFontHeightF = 0.008
   pnl_res.nglPanelLabelBarWidthF           = 0.4
   pnl_res.lbTitleString                    = f'{var_str[v]} [{var_unit[v]}]'
   pnl_res.lbTitleFontHeightF               = 0.012
   # pnl_res.lbTitleJust                      = 'BottomCenter'
   pnl_res.lbTitlePosition                  = 'Bottom'
   ngl.panel(wks_var,plot,[1,num_case],pnl_res)
   hc.trim_png(fig_file_var)
#---------------------------------------------------------------------------------------------------
for fig_file_var in fig_file_list:
   cmd = f"convert -trim  {fig_file_var}.{fig_type} -format '%wx%h\n' info:"
   (msg, err) = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, universal_newlines=True).communicate()
   dim_list = np.array( [ int(d) for d in msg.replace('\n','').split('x') ] )
   # max_dim = np.max(dim_list)
   # nx,ny = dim_list[0],dim_list[1]+100
   nx,ny = dim_list[0],dim_list[1]+50
   run_cmd(f'convert -trim -gravity center -extent {nx}x{ny} {fig_file_var}.{fig_type} {fig_file_var}.{fig_type}')
#-------------------------------------------------------------------------------
fig_file_list_str = ''
for fig_file_var in fig_file_list:
   fig_file_list_str += f' {fig_file_var}.{fig_type}'

run_cmd(f'montage {fig_file_list_str} -geometry +10 -tile 1x{len(fig_file_list)} {fig_file}.{fig_type}')
hc.trim_png(fig_file)
ngl.end()
# #---------------------------------------------------------------------------------------------------
# if 'num_plot_col' in locals():
#    layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
# else:
#    layout = [num_var,num_case]

# ngl.panel(wks,plot,layout,hs.setres_panel())
# ngl.end()

# hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

