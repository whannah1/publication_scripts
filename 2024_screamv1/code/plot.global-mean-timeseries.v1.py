# v4 is similar to v3 except panels are placed 'manually' using ImageMagick
#-------------------------------------------------------------------------------
import os, ngl, copy, string, xarray as xr, numpy as np, glob, dask, numba, cmocean, subprocess as sp
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
var_str = []
var_unit = []
file_type_list = []
obs_var_list = []
obs_file_type_list = []
def add_var(var_name,file_type,obs_var=None,obs_file_type=None,name='',unit=''): 
   var.append(var_name)
   file_type_list.append(file_type)
   obs_var_list.append(obs_var)
   obs_file_type_list.append(obs_file_type)
   var_str.append(name)
   var_unit.append(unit)
#-------------------------------------------------------------------------------
scrip_file_sim = '/pscratch/sd/w/whannah/e3sm_scratch/files_grid/ne1024pg2_scrip.nc'
scrip_file_obs = '/pscratch/sd/w/whannah/e3sm_scratch/files_grid/IMERG_1800x3600_scrip.nc'
scrip_file_era = '~/HICCUP/files_grid/scrip_ERA5_721x1440.nc'
sim_data_root  = '/global/cfs/cdirs/e3sm/gsharing/EAMxx'
obs_data_root  = '/pscratch/sd/w/whannah/Obs'

add_case('DYAMOND2_SCREAMv1' ,n='SCREAMv1 Jan 2020',p=sim_data_root,c='red')
add_case('Apr1_2013_SCREAMv1',n='SCREAMv1 Apr 2013',p=sim_data_root,c='green')
add_case('DYAMOND1_SCREAMv1' ,n='SCREAMv1 Aug 2016',p=sim_data_root,c='blue')
# add_case('Oct1_2013_SCREAMv1',n='SCREAMv1 Oct 2013',p=sim_data_root,c='purple')

# add_case('IMERG_Jan_2020',n='IMERG Jan 2020',p=obs_data_root,s='daily_QC_Jan_2020',obs=True,d=1,c='red')
# add_case('IMERG_Apr_2013',n='IMERG Apr 2013',p=obs_data_root,s='daily_QC_Apr_2013',obs=True,d=1,c='green')
# add_case('IMERG_Aug_2016',n='IMERG Aug 2016',p=obs_data_root,s='daily_QC_Aug_2016',obs=True,d=1,c='blue')
# add_case('IMERG_Oct_2013',n='IMERG Oct 2013',p=obs_data_root,s='daily_QC_Oct_2013',obs=True,d=1,c='purple')

#-------------------------------------------------------------------------------

# add_var('SW_flux_up@tom',  'output.scream.TOMVars.INSTANT',name='TOA SW Up')
# add_var('LW_flux_up@tom',  'output.scream.TOMVars.INSTANT',name='TOA LW Up')
# add_var('VapWaterPath',    'output.scream.VertIntegrals.INSTANT',name='water vapor path')
# add_var('IceWaterPath',    'output.scream.VertIntegrals.INSTANT',name='ice water path')
# add_var('T_2m',            'output.scream.SurfVars.INSTANT',name='2m Temperature')
# add_var('ps',              'output.scream.SurfVars.INSTANT',name='Psfc')
# add_var('wind_speed_10m',  'output.scream.SurfVars.INSTANT',name='') # inconsistent trends?
# add_var('qv_2m',           'output.scream.SurfVars.INSTANT',name='') # error with data?

add_var('precip','output.scream.SurfVars.INSTANT','precipitation','3B-HHR.MS.MRG.3IMERG',name='Precipitation')
# add_var('precip','output.scream.SurfVars.INSTANT','precipAvg','3B-DAY.MS.MRG.3IMERG',name='Precipitation')

# add_var('precip_liq_surf_mass','output.scream.SurfVars.INSTANT','precipAvg','3B-DAY.MS.MRG.3IMERG.2016')
# add_var('precip_ice_surf_mass','output.scream.SurfVars.INSTANT')
# add_var('surf_evap','output.scream.SurfVars.INSTANT')
# add_var('horiz_winds@bot','output.scream.SurfVars.INSTANT')

#-------------------------------------------------------------------------------
# tmp_data_path = os.getenv('HOME')+'/Research/E3SM/pub_figs/2023_screamv1_4season/data'
tmp_data_path = 'data'

fig_file,fig_type = f'figs/global-mean-timeseries-v1','png'

nday_data = 10
nday_plot = 10

recalculate = False

# var_x_case = False
num_plot_col = 2

subtitle_font_height = 0.025

print_stats = False
#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_var = len(case),len(var)

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*(num_var)

res = hs.res_xy()
res.vpHeightF = 0.5
res.vpHeightF = 0.4
res.tmYLLabelFontHeightF         = 0.015
res.tmXBLabelFontHeightF         = 0.015
res.tiXAxisFontHeightF           = 0.015
res.tiYAxisFontHeightF           = 0.015
res.xyLineThicknessF = 5

lres = hs.res_xy()
lres.xyDashPattern    = 1
lres.xyLineThicknessF = 2
lres.xyLineColor      = 'black'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

sim_grid_ds = xr.open_dataset( scrip_file_sim ).rename({'grid_size':'ncol'})
# sim_grid_ds['grid_center_lon'] = xr.where(sim_grid_ds['grid_center_lon']<0,sim_grid_ds['grid_center_lon']+360,sim_grid_ds['grid_center_lon'])
# sim_grid_ds['grid_corner_lon'] = xr.where(sim_grid_ds['grid_corner_lon']<0,sim_grid_ds['grid_corner_lon']+360,sim_grid_ds['grid_corner_lon'])
sim_area = sim_grid_ds['grid_area'].values

obs_grid_ds = xr.open_dataset( scrip_file_obs ).rename({'grid_size':'ncol'})
obs_grid_ds['grid_center_lon'] = xr.where(obs_grid_ds['grid_center_lon']<0,obs_grid_ds['grid_center_lon']+360,obs_grid_ds['grid_center_lon'])
obs_grid_ds['grid_corner_lon'] = xr.where(obs_grid_ds['grid_corner_lon']<0,obs_grid_ds['grid_corner_lon']+360,obs_grid_ds['grid_corner_lon'])
obs_area = obs_grid_ds['grid_area'].values

era_grid_ds = xr.open_dataset( scrip_file_era ).rename({'grid_size':'ncol'})
era_grid_ds['grid_center_lon'] = xr.where(era_grid_ds['grid_center_lon']<0,era_grid_ds['grid_center_lon']+360,era_grid_ds['grid_center_lon'])
era_grid_ds['grid_corner_lon'] = xr.where(era_grid_ds['grid_corner_lon']<0,era_grid_ds['grid_corner_lon']+360,era_grid_ds['grid_corner_lon'])

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def calculate_obs_area(lon,lat,lon_bnds,lat_bnds):
   re = 6.37122e06  # radius of earth
   nlat,nlon = len(lat),len(lon)
   area = np.empty((nlat,nlon),np.float64)
   for j in range(nlat):
      for i in range(nlon):
         dlon = np.absolute( lon_bnds[j,1] - lon_bnds[j,0] )
         dlat = np.absolute( lat_bnds[j,1] - lat_bnds[j,0] )
         dx = re*dlon*np.pi/180.
         dy = re*dlat*np.pi/180.
         area[j,i] = dx*dy
   return area
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   time_list = []
   data_list = []
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
      tmp_file = f'{tmp_data_path}/global-mean-timeseries.v1.tmp.nday_{nday_data}.{tvar}.{case[c]}.nc'

      if recalculate :
         print(' '*6+'recalculating...')
         #----------------------------------------------------------------------
         # idenfity the files to load
         tcase = case[c]
         if 'IMERG' in tcase: tcase = 'IMERG'
         file_path = f'{case_dir[c]}/{tcase}/{case_sub[c]}/{tfile_type}*'
         file_list = sorted(glob.glob(file_path))
         # trim down file list
         nf = nday_data
         # if obs_flag[c]: nf = nday_data*48 # IMERG frequency is 30min
         file_list = file_list[:nf] # use initial files
         #----------------------------------------------------------------------
         print(' '*6+f'Loading data ({tvar})...')
         for f in range(len(file_list)):
            print(f' '*8+f'f: {f:03d}  file: {hc.tcolor.YELLOW}{file_list[f]}{hc.tcolor.ENDC}')
            #-------------------------------------------------------------------
            group_name = None # group_name = 'Grid' if case[c]=='IMERG' else None
            ds = xr.open_dataset( file_list[f], group=group_name)
            #-------------------------------------------------------------------
            # Load the data
            if tvar=='precip':
               data = ds['precip_liq_surf_mass'] + ds['precip_ice_surf_mass']
            else:
               data = ds[tvar]
            if obs_flag[c]: data = data.stack(ncol=('lat','lon'))
            #-------------------------------------------------------------------
            # build time coodinate
            time_tmp = data['time.year'].values*365    \
                      +data['time.dayofyear'].values   \
                      +data['time.hour'].values/24     \
                      +data['time.minute'].values/24/60
            #-------------------------------------------------------------------
            # calculate global mean
            if f==0:
               if obs_flag[c]:
                  area = obs_area
                  # xy_dims = ('longitude','latitude')
                  # xy_dims = ('lat','lon')
                  # area = calculate_obs_area(ds['lon'].values,ds['lat'].values,ds['lon_bnds'].values,ds['lat_bnds'].values)
               else:
                  area = sim_area
               xy_dims = ('ncol')
            #-------------------------------------------------------------------
            # gbl_mean = ( (data*area).sum(dim=xy_dims) / area.sum(dim=xy_dims) )
            gbl_mean_tmp = (data*area).sum(dim=xy_dims) / np.sum(area)
            # np.ma.masked_invalid(var_out.values)
            gbl_mean_tmp['time'] = time_tmp
            # gbl_mean_tmp = xr.DataArray( np.zeros( (ntime,nbins) ), \
            #                         coords=[('time',time_tmp),('lon',lon_bins)], \
            #                         dims=['time','lon'] )
            #-------------------------------------------------------------------
            if f==0:
               gbl_mean = gbl_mean_tmp
               time     = time_tmp
            else:
               gbl_mean = xr.concat([gbl_mean,gbl_mean_tmp], dim='time')
               time     = np.append(time,time_tmp)
         #-------------------------------------------------------------------
         # convert units
         if var[v]=='precip':
            if obs_flag[c]:
               gbl_mean = gbl_mean*24. # mm/hr to mm/day
            else:
               gbl_mean = (gbl_mean/100)*86400 # kg/m2 to mm/day using dtime=100 (ne1024)
         #----------------------------------------------------------------------
         # print some summary stats after averaging globally
         hc.print_stat(gbl_mean,name=tvar,compact=True,indent=' '*6)
         #----------------------------------------------------------------------
         # Write to file 
         print(' '*6+f'Writing data to file: {tmp_file}')
         gbl_mean.name = tvar
         tmp_ds = xr.Dataset()
         # tmp_ds['time'] = time
         tmp_ds[tvar]   = gbl_mean
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         print(' '*6+f'Reading pre-calculated data from file: {hc.tcolor.MAGENTA}{tmp_file}{hc.tcolor.ENDC}')
         tmp_ds = xr.open_dataset( tmp_file )
         # time_tmp = tmp_ds['time'].values
         gbl_mean = tmp_ds[tvar]
      #-------------------------------------------------------------------------
      time_tmp = gbl_mean['time'].values
      time_tmp = time_tmp - time_tmp[0]
      time_list.append( time_tmp )
      data_list.append( gbl_mean.values )
   #----------------------------------------------------------------------------
   # # skip days 0-4 to focus on later adjustment
   # sample_per_day = int(24*60/15) # 15 minute sampling fro SCREAM
   # for c in range(num_case):
   #    time_list[c] = time_list[c][5*sample_per_day:]
   #    data_list[c] = data_list[c][5*sample_per_day:]
   #----------------------------------------------------------------------------
   # Create plot
   data_min = np.min([np.min(d) for d in data_list])
   data_max = np.max([np.max(d) for d in data_list])
   time_min = np.min([np.min(d) for d in time_list])
   time_max = np.max([np.max(d) for d in time_list])
   res.trYMinF = data_min
   res.trYMaxF = data_max
   res.trXMinF = time_min
   if 'nday_plot' in locals():
      res.trXMaxF = nday_plot
   else:
      res.trXMaxF = time_max
   res.tiXAxisString = 'Time [days]'
   # res.tiYAxisString = f'{var_str[v]}'
   # res.tiYAxisString = f'[{var_unit[v]}]'
   for c in range(num_case):
      #-------------------------------------------------------------------------
      # calculate ref value as mean of last 2 days
      sample_per_day = int(24*60/15) # 15 minute sampling fro SCREAM
      nsample = 2*sample_per_day 
      ref_value = np.mean(data_list[c][-1*nsample:])
      #-------------------------------------------------------------------------
      # subtract reference value to collapse curves
      data_list[c] = data_list[c] - ref_value
   data_min = np.min([np.min(d) for d in data_list])
   data_max = np.max([np.max(d) for d in data_list])
   res.trYMinF = data_min
   res.trYMaxF = data_max
   for c in range(num_case):
      #-------------------------------------------------------------------------
      tres = copy.deepcopy(res)
      tres.xyLineColor   = clr[c]
      tres.xyDashPattern = dsh[c]
      tplot = ngl.xy(wks, time_list[c], data_list[c], tres) 
      if c==0:
         plot[v] = tplot
      else:
         ngl.overlay(plot[v],tplot)
      #-------------------------------------------------------------------------
      # add line to indicate reference value
      ref_value = np.mean(data_list[c][-1*nsample:])
      ref_array = np.ones(time_list[c].shape) * ref_value
      lres.xyLineColor   = clr[c]
      lres.xyDashPattern = 1
      lplot = ngl.xy(wks, time_list[c], ref_array, lres)
      ngl.overlay(plot[v],lplot)
   #----------------------------------------------------------------------------
   hs.set_subtitles(wks, plot[v], '', '', tvar, font_height=subtitle_font_height)
#---------------------------------------------------------------------------------------------------
if 'num_plot_col' in locals():
   layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
else:
   layout = [1,num_var]


ngl.panel(wks,plot,layout,hs.setres_panel())
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
