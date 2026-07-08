#-------------------------------------------------------------------------------
import os, copy, string, xarray as xr, numpy as np, glob, dask, numba, cmocean, subprocess as sp
import hapy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
fig_file,fig_type = f'figs/global-mean-timeseries-v1','png'
#-------------------------------------------------------------------------------
case_name,case,case_dir,case_sub = [],[],[],[]
clr,dsh,mrk = [],[],[]
obs_flag = []
htype_alt_list = []
def add_case(case_in,n='',p=None,s='',g=None,d=0,c=None,m=0,r=False,obs=False,htype_alt=False):
   global case_name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); case_name.append(n)
   case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
   obs_flag.append(obs)
   htype_alt_list.append(htype_alt)
#-------------------------------------------------------------------------------
var = []
var_str = []
var_unit = []
file_type_list = []
obs_var_list = []
obs_file_type_list = []
lev_list = []
unit_fac_list = []
def add_var(var_name,file_type,obs_var=None,obs_file_type=None,name=None,unit='',lev=None,unit_fac=None):
   var.append(var_name)
   file_type_list.append(file_type)
   obs_var_list.append(obs_var)
   obs_file_type_list.append(obs_file_type)
   # var_str.append(name)
   var_str.append(var_name if name is None else name)
   var_unit.append(unit)
   lev_list.append(lev)
   unit_fac_list.append(unit_fac)
#-------------------------------------------------------------------------------
# scrip_file_sim = '/pscratch/sd/w/whannah/e3sm_scratch/files_grid/ne1024pg2_scrip.nc'
scrip_file_obs = '/pscratch/sd/w/whannah/e3sm_scratch/files_grid/IMERG_1800x3600_scrip.nc'
# scrip_file_era = '~/HICCUP/files_grid/scrip_ERA5_721x1440.nc'
# sim_data_root  = '/global/cfs/cdirs/e3sm/gsharing/EAMxx'
# obs_data_root  = '/pscratch/sd/w/whannah/Obs'

# scrip_file_sim = '/lustre/orion/cli115/proj-shared/hannah6/HICCUP/data/scrip_ne1024pg2.nc'
# scrip_file_sim = '/lustre/orion/cli115/proj-shared/hannah6/HICCUP/data/scrip_ne30pg2.nc'
# scrip_file_obs = '~/HICCUP/data_scratch/scrip_ERA5_721x1440.nc'
# obs_data_root  = '/lustre/orion/cli115/proj-shared/hannah6/obs_data'


# tmp_scratch='/lustre/orion/cli115/proj-shared/noel/e3sm_scratch/s10-feb7/dd1024'
tmp_scratch='/global/cfs/cdirs/e3sm/jpaige3/dy1ne1024'
tmp_sub = 'SCREAM.2024-autocal-00.ne1024pg2/run'

# for n in range( 1,12+1): add_case(f'm{n:04d}', p=tmp_scratch,s=tmp_sub)
# for n in range(13,18+1): add_case(f'm{n:04d}', p=tmp_scratch,s=tmp_sub)
# for n in range(20,24+1): add_case(f'm{n:04d}', p=tmp_scratch,s=tmp_sub)

add_case('seed',c='black', p=tmp_scratch,s=tmp_sub)
add_case('m0001',p=tmp_scratch,s=tmp_sub)
add_case('m0003',p=tmp_scratch,s=tmp_sub)
add_case('m0004',p=tmp_scratch,s=tmp_sub)
add_case('m0005',p=tmp_scratch,s=tmp_sub)
add_case('m0006',p=tmp_scratch,s=tmp_sub)
add_case('m0009',p=tmp_scratch,s=tmp_sub)
add_case('m0015',p=tmp_scratch,s=tmp_sub)
add_case('m0017',p=tmp_scratch,s=tmp_sub)
add_case('m0018',p=tmp_scratch,s=tmp_sub)
add_case('m0020',p=tmp_scratch,s=tmp_sub)
add_case('m0023',p=tmp_scratch,s=tmp_sub)
add_case('m0026',p=tmp_scratch,s=tmp_sub)
add_case('m0027',p=tmp_scratch,s=tmp_sub)
add_case('m0028',p=tmp_scratch,s=tmp_sub)
add_case('m0032',p=tmp_scratch,s=tmp_sub)
add_case('m0033',p=tmp_scratch,s=tmp_sub)
add_case('m0034',p=tmp_scratch,s=tmp_sub)
add_case('m0035',p=tmp_scratch,s=tmp_sub)
add_case('m0037',p=tmp_scratch,s=tmp_sub)
add_case('m0038',p=tmp_scratch,s=tmp_sub)
add_case('m0039',p=tmp_scratch,s=tmp_sub)
add_case('m0040',p=tmp_scratch,s=tmp_sub)
add_case('m0042',p=tmp_scratch,s=tmp_sub)
add_case('m0049',p=tmp_scratch,s=tmp_sub)
add_case('m0050',p=tmp_scratch,s=tmp_sub)
add_case('m0051',p=tmp_scratch,s=tmp_sub)
add_case('m0053',p=tmp_scratch,s=tmp_sub)
add_case('m0054',p=tmp_scratch,s=tmp_sub)
add_case('m0057',p=tmp_scratch,s=tmp_sub)
add_case('m0058',p=tmp_scratch,s=tmp_sub)
add_case('m0060',p=tmp_scratch,s=tmp_sub)
add_case('m0062',p=tmp_scratch,s=tmp_sub)
add_case('m0063',p=tmp_scratch,s=tmp_sub)
add_case('m0065',p=tmp_scratch,s=tmp_sub)
add_case('m0066',p=tmp_scratch,s=tmp_sub)
add_case('m0067',p=tmp_scratch,s=tmp_sub)
add_case('m0068',p=tmp_scratch,s=tmp_sub)
add_case('m0069',p=tmp_scratch,s=tmp_sub)
add_case('m0073',p=tmp_scratch,s=tmp_sub)
add_case('m0074',p=tmp_scratch,s=tmp_sub)
add_case('m0079',p=tmp_scratch,s=tmp_sub)
add_case('m0080',p=tmp_scratch,s=tmp_sub)
add_case('m0085',p=tmp_scratch,s=tmp_sub)
add_case('m0087',p=tmp_scratch,s=tmp_sub)
add_case('m0088',p=tmp_scratch,s=tmp_sub)
add_case('m0089',p=tmp_scratch,s=tmp_sub)
add_case('m0090',p=tmp_scratch,s=tmp_sub)
add_case('m0091',p=tmp_scratch,s=tmp_sub)
add_case('m0092',p=tmp_scratch,s=tmp_sub)
add_case('m0094',p=tmp_scratch,s=tmp_sub)
add_case('m0095',p=tmp_scratch,s=tmp_sub)
add_case('m0096',p=tmp_scratch,s=tmp_sub)
add_case('m0099',p=tmp_scratch,s=tmp_sub)
add_case('m0101',p=tmp_scratch,s=tmp_sub)
add_case('m0106',p=tmp_scratch,s=tmp_sub)
add_case('m0109',p=tmp_scratch,s=tmp_sub)
add_case('m0110',p=tmp_scratch,s=tmp_sub)
add_case('m0111',p=tmp_scratch,s=tmp_sub)
add_case('m0114',p=tmp_scratch,s=tmp_sub)
add_case('m0115',p=tmp_scratch,s=tmp_sub)
add_case('m0117',p=tmp_scratch,s=tmp_sub)
add_case('m0119',p=tmp_scratch,s=tmp_sub)
add_case('m0121',p=tmp_scratch,s=tmp_sub)
add_case('m0124',p=tmp_scratch,s=tmp_sub)
add_case('m0125',p=tmp_scratch,s=tmp_sub)
add_case('m0128',p=tmp_scratch,s=tmp_sub)
add_case('m0131',p=tmp_scratch,s=tmp_sub)
add_case('m0132',p=tmp_scratch,s=tmp_sub)
add_case('m0133',p=tmp_scratch,s=tmp_sub)
add_case('m0134',p=tmp_scratch,s=tmp_sub)
add_case('m0135',p=tmp_scratch,s=tmp_sub)
add_case('m0136',p=tmp_scratch,s=tmp_sub)
add_case('m0137',p=tmp_scratch,s=tmp_sub)
add_case('m0140',p=tmp_scratch,s=tmp_sub)
add_case('m0147',p=tmp_scratch,s=tmp_sub)
add_case('m0148',p=tmp_scratch,s=tmp_sub)
add_case('m0153',p=tmp_scratch,s=tmp_sub)
add_case('m0154',p=tmp_scratch,s=tmp_sub)
add_case('m0159',p=tmp_scratch,s=tmp_sub)
add_case('m0160',p=tmp_scratch,s=tmp_sub)
add_case('m0161',p=tmp_scratch,s=tmp_sub)
add_case('m0162',p=tmp_scratch,s=tmp_sub)
add_case('m0163',p=tmp_scratch,s=tmp_sub)
add_case('m0164',p=tmp_scratch,s=tmp_sub)
add_case('m0165',p=tmp_scratch,s=tmp_sub)
add_case('m0166',p=tmp_scratch,s=tmp_sub)
add_case('m0167',p=tmp_scratch,s=tmp_sub)
add_case('m0168',p=tmp_scratch,s=tmp_sub)
add_case('m0170',p=tmp_scratch,s=tmp_sub)
add_case('m0174',p=tmp_scratch,s=tmp_sub)
add_case('m0176',p=tmp_scratch,s=tmp_sub)
add_case('m0180',p=tmp_scratch,s=tmp_sub)
add_case('m0181',p=tmp_scratch,s=tmp_sub)
add_case('m0183',p=tmp_scratch,s=tmp_sub)
add_case('m0185',p=tmp_scratch,s=tmp_sub)
add_case('m0186',p=tmp_scratch,s=tmp_sub)
add_case('m0190',p=tmp_scratch,s=tmp_sub)
add_case('m0193',p=tmp_scratch,s=tmp_sub)
add_case('m0197',p=tmp_scratch,s=tmp_sub)
add_case('m0199',p=tmp_scratch,s=tmp_sub)
add_case('m0200',p=tmp_scratch,s=tmp_sub)
add_case('m0201',p=tmp_scratch,s=tmp_sub)
add_case('m0203',p=tmp_scratch,s=tmp_sub)
add_case('m0205',p=tmp_scratch,s=tmp_sub)
add_case('m0207',p=tmp_scratch,s=tmp_sub)
add_case('m0213',p=tmp_scratch,s=tmp_sub)
add_case('m0216',p=tmp_scratch,s=tmp_sub)
add_case('m0217',p=tmp_scratch,s=tmp_sub)
add_case('m0218',p=tmp_scratch,s=tmp_sub)
add_case('m0224',p=tmp_scratch,s=tmp_sub)
add_case('m0225',p=tmp_scratch,s=tmp_sub)
add_case('m0226',p=tmp_scratch,s=tmp_sub)
add_case('m0227',p=tmp_scratch,s=tmp_sub)
add_case('m0228',p=tmp_scratch,s=tmp_sub)
add_case('m0229',p=tmp_scratch,s=tmp_sub)
add_case('m0230',p=tmp_scratch,s=tmp_sub)
add_case('m0231',p=tmp_scratch,s=tmp_sub)
add_case('m0232',p=tmp_scratch,s=tmp_sub)
add_case('m0237',p=tmp_scratch,s=tmp_sub)
add_case('m0238',p=tmp_scratch,s=tmp_sub)
add_case('m0240',p=tmp_scratch,s=tmp_sub)
add_case('m0243',p=tmp_scratch,s=tmp_sub)
add_case('m0244',p=tmp_scratch,s=tmp_sub)
add_case('m0246',p=tmp_scratch,s=tmp_sub)
add_case('m0248',p=tmp_scratch,s=tmp_sub)
add_case('m0251',p=tmp_scratch,s=tmp_sub)
add_case('m0252',p=tmp_scratch,s=tmp_sub)
add_case('m0254',p=tmp_scratch,s=tmp_sub)
add_case('m0255',p=tmp_scratch,s=tmp_sub)
add_case('m0257',p=tmp_scratch,s=tmp_sub)
add_case('m0259',p=tmp_scratch,s=tmp_sub)
add_case('m0260',p=tmp_scratch,s=tmp_sub)
add_case('m0280',p=tmp_scratch,s=tmp_sub)
add_case('m0281',p=tmp_scratch,s=tmp_sub)
add_case('m0283',p=tmp_scratch,s=tmp_sub)
add_case('m0284',p=tmp_scratch,s=tmp_sub)
add_case('m0285',p=tmp_scratch,s=tmp_sub)
add_case('m0286',p=tmp_scratch,s=tmp_sub)


# add_case('m0123',p=tmp_scratch,s=tmp_sub)
# add_case('m0138',p=tmp_scratch,s=tmp_sub)


# clr = hapy.fill_color_list(clr,cmap=None,randomize=True)
# clr = hapy.fill_color_list(clr,cmap='gist_ncar',randomize=True)
# clr = hapy.fill_color_list(clr,cmap='plasma',randomize=True)
clr = hapy.fill_color_list(clr,cmap='rainbow',randomize=True)
# clr = hapy.fill_color_list(clr,cmap='brg',randomize=True)

# for c in clr:
#    print(c)
# exit()

#-------------------------------------------------------------------------------
# grid_root = '/lustre/orion/cli115/proj-shared/hannah6/HICCUP/data'
# grid_root = os.getenv('HOME')+'/E3SM/data_grid'
grid_root = '/pscratch/sd/w/whannah/e3sm_scratch/files_grid'

scrip_file_sim = f'{grid_root}/ne30pg2_scrip.nc'
tmp_file_type = 'output.scream.AutoCal.hourly_inst_ne30pg2.INSTANT.nhours_x1.'
# tmp_file_type = 'output.scream.AutoCal.daily_avg_ne30pg2.AVERAGE.nhours_x24'

# scrip_file_sim = f'{grid_root}/scrip_ne120pg2.nc'
# tmp_file_type = 'output.scream.AutoCal.3hourly_avg_ne120pg2.AVERAGE.nhours_x3.'

# scrip_file_sim = f'{grid_root}/ne256pg2_scrip.nc'
# tmp_file_type = 'output.scream.AutoCal.hourly_inst_ne256pg2.INSTANT.nhours_x1'

# scrip_file_sim = f'{grid_root}/ne1024pg2_scrip.nc'
# tmp_file_type = 'output.scream.AutoCal.hourly_inst_ne1024pg2.INSTANT.nhours_x1.'

# scrip_file_sim = f'{grid_root}/ne256pg2_scrip.nc'
# tmp_file_type = 'output.scream.2D.1hr.ne256pg2.INSTANT.nhours_x1'

#-------------------------------------------------------------------------------

# tmp_file_type  = 'output.scream.2D.1hr.inst.INSTANT.nhours_x1'
# tmp_file_type2 = 'output.scream.2D.1hr.ne256pg2.INSTANT.nhours_x1'

# add_var('precip_total_surf_mass_flux',tmp_file_type,name=None,unit_fac=86400*1e3)
# add_var('zm_prec',tmp_file_type,unit_fac=86400*1e3)

# add_var('precip_total_surf_mass_flux', tmp_file_type, name='precip')
# add_var('LW_flux_up_at_model_top',     tmp_file_type, name='LW Flux up @ TOA')
# add_var('LW_flux_dn_at_model_top',     tmp_file_type, name='LW Flux dn @ TOA')
# add_var('SW_flux_up_at_model_top',     tmp_file_type, name='SW Flux up @ TOA')
# add_var('SW_flux_dn_at_model_top',     tmp_file_type, name='SW Flux dn @ TOA')

# add_var('LW_flux_net_at_model_top',     tmp_file_type, name='Net LW Flux @ TOA')
# add_var('SW_flux_net_at_model_top',     tmp_file_type, name='Net SW Flux @ TOA')

# add_var('ps',                          tmp_file_type, name='sfc pressure',unit_fac=1e-2, unit='mb')
# add_var('VapWaterPath',                tmp_file_type, name='VapWP')
add_var('LiqWaterPath',                tmp_file_type, name='TLWP (kg/m2)')
# add_var('IceWaterPath',                tmp_file_type, name='TIWP [kg/m2]', unit_fac=1e3)

# add_var('T_mid',                       tmp_file_type, name='T500',lev=-21)
# add_var('T_mid',                       tmp_file_type, name='T850',lev=-30)
# add_var('T_2m',                        tmp_file_type,obs_var='t2m', obs_file_type='daily.sfc', name='2m Temperature')

# add_var('U',                           tmp_file_type,name='U100',lev=-10)
# add_var('U',                           tmp_file_type,name='U500',lev=-21)
# add_var('U',                           tmp_file_type,name='U850',lev=-30)


# add_var('SW_flux_up@tom',  'output.scream.TOMVars.INSTANT',name='TOA SW Up')
# add_var('LW_flux_up@tom',  'output.scream.TOMVars.INSTANT',name='TOA LW Up')
# add_var('VapWaterPath',    'output.scream.VertIntegrals.INSTANT',name='water vapor path')
# add_var('IceWaterPath',    'output.scream.VertIntegrals.INSTANT',name='ice water path')
# add_var('T_2m',            'output.scream.SurfVars.INSTANT',name='2m Temperature')
# add_var('ps',              'output.scream.SurfVars.INSTANT',name='Psfc')
# add_var('wind_speed_10m',  'output.scream.SurfVars.INSTANT',name='') # inconsistent trends?
# add_var('qv_2m',           'output.scream.SurfVars.INSTANT',name='') # error with data?



# add_var('T_2m',                        'scream.h1.INSTANT.nhours_x1',obs_var='t2m', obs_file_type='daily.sfc', name='2m Temperature')
# add_var('VapWaterPath',                'scream.h1.INSTANT.nhours_x1',obs_var='tcw', obs_file_type='daily.sfc', name='Column Water Vapor')
# add_var('LiqWaterPath',                'scream.h1.INSTANT.nhours_x1',obs_var='tclw',obs_file_type='daily.sfc', name='Liquid Water Path')
# add_var('IceWaterPath',                'scream.h1.INSTANT.nhours_x1',obs_var='tciw',obs_file_type='daily.sfc', name='Ice Water Path')
# add_var('ps',                          'scream.h1.INSTANT.nhours_x1',obs_var='sp',  obs_file_type='daily.sfc', name='Sfc Pressure')
# add_var('precip_total_surf_mass_flux', 'scream.h1.INSTANT.nhours_x1',obs_var='tp',  obs_file_type='daily.sfc', name='Total Precipitation')
# add_var('horiz_winds_at_model_bot_u',      'scream.h1.INSTANT.nhours_x1',obs_var='u10',  obs_file_type='daily.sfc', name='U "near sfc"')

# add_var('U_at_500hPa',      'scream.h1.INSTANT.nhours_x1',obs_var='u',  obs_file_type='daily.atm', name='U500')

# add_var('precip','output.scream.SurfVars.INSTANT','precipitation','3B-HHR.MS.MRG.3IMERG',name='Precipitation')
# add_var('precip','output.scream.SurfVars.INSTANT','precipAvg','3B-DAY.MS.MRG.3IMERG',name='Precipitation')

# add_var('precip_liq_surf_mass','output.scream.SurfVars.INSTANT','precipAvg','3B-DAY.MS.MRG.3IMERG.2016')
# add_var('precip_ice_surf_mass','output.scream.SurfVars.INSTANT')
# add_var('surf_evap','output.scream.SurfVars.INSTANT')
# add_var('horiz_winds@bot','output.scream.SurfVars.INSTANT')

#-------------------------------------------------------------------------------
# set lat bounds for regional subsets

# lat1,lat2 = -30,30 # tropical belt

# xlat,xlon,dy,dx =  40,360-100,10,10;
# xlat,xlon,dy,dx =  40,360-90,5,5;

# xlat,xlon,dy,dx =  0,180,10,10;
# xlat,xlon,dy,dx = 40,180,5,5;
# xlat,xlon,dy,dx = 80,180,5,5;
# xlat,xlon,dy,dx = 15,360-40,10,10;
# xlat,xlon,dy,dx = 60,60,1,1;
# xlat,xlon,dy,dx = 60,60,5,5;
# xlat,xlon,dy,dx = 60,60,10,10;

# xlat,xlon,dy,dx = 0,75,1,1;
# xlat,xlon,dy,dx = 0,75,20,20;

if 'xlat' in locals(): lat1,lat2,lon1,lon2 = xlat-dy/2,xlat+dy/2,xlon-dx/2,xlon+dx/2

#-------------------------------------------------------------------------------
tmp_data_path = '/global/homes/w/whannah/publication_scripts/2026_scream_autocalibration/data'

nday_data = 5
nday_plot = 5

recalculate = True

plot_as_anomaly = False

# -------------------------------------------------------------------
# show actual dates (mm-dd-yyyy) on the x-axis instead of day-of-year
use_date_axis = True

# var_x_case = False
# num_plot_col = len(var)
num_plot_col = 1

subtitle_font_height = 0.015
# subtitle_font_height = 0.005

print_stats = False
#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_var = len(case),len(var)

#---------------------------------------------------------------------------------------------------
# determine panel layout and create figure and axes
if 'num_plot_col' in locals():
   panel_cols = num_plot_col
   panel_rows = int(np.ceil(num_var/float(num_plot_col)))
else:
   panel_rows,panel_cols = 1,num_var

fig,axs = plt.subplots( panel_rows, panel_cols,
                        figsize=(8*panel_cols,4*panel_rows),
                        squeeze=False )
axs = axs.flatten()

#---------------------------------------------------------------------------------------------------
# plot style parameters
line_thickness       = 1.5   # main time series line width
ref_line_thickness   = 1.0   # reference/anomaly line width
label_font_size      = 12
title_font_size      = 14

#---------------------------------------------------------------------------------------------------
# map NGL dash patterns to matplotlib line styles
def dash_to_linestyle(d):
   return {0:'-',1:'--',2:':',3:'-.'}.get(d,'-')

# print();print(clr)

# # use random colors for all but first case
# ngl.define_colormap(wks,'ncl_default')
# clr[1:] = np.linspace(2,len( ngl.retrieve_colormap(wks) )-1,num_case-1,dtype=int)

# print();print(clr)
# exit()

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

sim_grid_ds = xr.open_dataset( scrip_file_sim ).rename({'grid_size':'ncol'})
# sim_grid_ds['grid_center_lon'] = xr.where(sim_grid_ds['grid_center_lon']<0,sim_grid_ds['grid_center_lon']+360,sim_grid_ds['grid_center_lon'])
# sim_grid_ds['grid_corner_lon'] = xr.where(sim_grid_ds['grid_corner_lon']<0,sim_grid_ds['grid_corner_lon']+360,sim_grid_ds['grid_corner_lon'])
sim_area = sim_grid_ds['grid_area']

# obs_grid_ds = xr.open_dataset( scrip_file_obs ).rename({'grid_size':'ncol'})
# obs_grid_ds['grid_center_lon'] = xr.where(obs_grid_ds['grid_center_lon']<0,obs_grid_ds['grid_center_lon']+360,obs_grid_ds['grid_center_lon'])
# obs_grid_ds['grid_corner_lon'] = xr.where(obs_grid_ds['grid_corner_lon']<0,obs_grid_ds['grid_corner_lon']+360,obs_grid_ds['grid_corner_lon'])
# obs_area = obs_grid_ds['grid_area']

# era_grid_ds = xr.open_dataset( scrip_file_era ).rename({'grid_size':'ncol'})
# era_grid_ds['grid_center_lon'] = xr.where(era_grid_ds['grid_center_lon']<0,era_grid_ds['grid_center_lon']+360,era_grid_ds['grid_center_lon'])
# era_grid_ds['grid_corner_lon'] = xr.where(era_grid_ds['grid_corner_lon']<0,era_grid_ds['grid_corner_lon']+360,era_grid_ds['grid_corner_lon'])

#---------------------------------------------------------------------------------------------------
if 'lat1' in locals():

   scrip_lat_name,scrip_lon_name = 'grid_center_lat','grid_center_lon'

   print('Creating obs mask...')
   obs_mask = xr.DataArray( np.ones(len(obs_grid_ds[scrip_lat_name]),dtype=bool), coords=obs_grid_ds[scrip_lat_name].coords )
   if 'lat1' in locals(): obs_mask = obs_mask & (obs_grid_ds[scrip_lat_name]>=lat1) & (obs_grid_ds[scrip_lat_name]<=lat2)
   if 'lon1' in locals(): obs_mask = obs_mask & (obs_grid_ds[scrip_lon_name]>=lon1) & (obs_grid_ds[scrip_lon_name]<=lon2)
   obs_area = obs_area.where(obs_mask,drop=True)

   print('Creating sim mask...')
   sim_mask = xr.DataArray( np.ones(len(sim_grid_ds[scrip_lat_name]),dtype=bool), coords=sim_grid_ds[scrip_lat_name].coords )
   if 'lat1' in locals(): sim_mask = sim_mask & (sim_grid_ds[scrip_lat_name]>=lat1) & (sim_grid_ds[scrip_lat_name]<=lat2)
   if 'lon1' in locals(): sim_mask = sim_mask & (sim_grid_ds[scrip_lon_name]>=lon1) & (sim_grid_ds[scrip_lon_name]<=lon2)
   sim_area = sim_area.where(sim_mask,drop=True)
   
   print('done.')

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_file_list(case,case_dir,case_sub,file_type):
   file_path = f'{case_dir}/{case}/{case_sub}/*{file_type}*'
   file_list = sorted(glob.glob(file_path))
   if 'first_file' in globals(): file_list = file_list[first_file:]
   if 'num_files'  in globals(): file_list = file_list[:num_files]

   if case!='seed': file_list = file_list[:2]

   # print(); print(file_list); print(); exit()
   if file_list==[]:
      print('ERROR: Empty file list:')
      print(); print(file_path)
      print(); print(file_list)
      exit()
   return file_list
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
def get_ctr_str(glb_avg=None):
   ctr_str = ''
   if 'lat1' in globals():
      lat1_str = f'{lat1}N' if lat1>=0 else f'{(lat1*-1)}S'
      lat2_str = f'{lat2}N' if lat2>=0 else f'{(lat2*-1)}S'
      ctr_str += f' {lat1_str}:{lat2_str} '
   if 'lon1' in globals():
      lon1_str = f'{lon1}E' #if lon1>=0 and lon1<=360 else f'{(lon1*-1)}S'
      lon2_str = f'{lon2}E' #if lon2>=0 and lon2<=360 else f'{(lon2*-1)}S'
      ctr_str += f' {lon1_str}:{lon2_str} '
   # if glb_avg is not None:
   #    # add logic here t display global average value?
   return ctr_str
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   time_list = []
   data_list = []
   date_list = []
   for c in range(num_case):
      if obs_flag[c]:
         tfile_type = obs_file_type_list[v]
         tvar = obs_var_list[v]
      else:
         tfile_type = file_type_list[v]
         tvar = var[v]
         if var[v]=='horiz_winds_at_model_bot_u': tvar = 'horiz_winds_at_model_bot'
         if var[v]=='horiz_winds_at_model_bot_v': tvar = 'horiz_winds_at_model_bot'
      if c==0: print(' '*2+'var: '+hapy.tclr.GREEN+tvar+hapy.tclr.END)
      print('\n'+' '*4+'case: '+hapy.tclr.CYAN+case[c]+hapy.tclr.END)
      #-------------------------------------------------------------------------
      tmp_file = f'{tmp_data_path}/global-mean-timeseries.v1.tmp.nday_{nday_data}.{tvar}.{case[c]}.nc'

      if recalculate :
         print(' '*6+'recalculating...')
         #----------------------------------------------------------------------
         # idenfity the files to load
         tcase = case[c]
         # if 'IMERG' in tcase: tcase = 'IMERG'
         # if 'ERA5'  in tcase: tcase = 'ERA5'
         # file_path = f'{case_dir[c]}/{tcase}/{case_sub[c]}/{tcase}.{tfile_type}*'
         tmp_file_type = file_type_list[v]
         if htype_alt_list[c]: tmp_file_type = tmp_file_type2
         file_list = get_file_list(case[c],case_dir[c],case_sub[c],tmp_file_type)

         # file_path = f'{case_dir[c]}/{tcase}/{case_sub[c]}/{tfile_type}*'
         # file_list = sorted(glob.glob(file_path))

         # trim down file list
         # nf = nday_data
         # if obs_flag[c]: nf = nday_data*48 # IMERG frequency is 30min
         # file_list = file_list[:nf] # use initial files
         if len(file_list)>2:
            file_list = file_list[:len(file_list)-1] # skip last (incomplete) file
         #----------------------------------------------------------------------
         # if file_list==[]:
         #    print('ERROR: Empty file list:')
         #    # print(); print(file_path)
         #    print(); print(file_list)
         #    exit()
         #----------------------------------------------------------------------
         print(' '*6+f'Loading data ({tvar})...')
         for f in range(len(file_list)):
            print(f' '*8+f'f: {f:03d}  file: {hapy.tclr.YELLOW}{file_list[f]}{hapy.tclr.END}')
            #-------------------------------------------------------------------
            group_name = None # group_name = 'Grid' if case[c]=='IMERG' else None
            ds = xr.open_dataset( file_list[f], group=group_name)
            #-------------------------------------------------------------------
            # Load the data
            if tvar=='precip':
               data = ds['precip_liq_surf_mass'] + ds['precip_ice_surf_mass']
            # elif tvar=='LW_flux_net_at_model_top':
            #    data = ds['LW_flux_up_at_model_top'] + ds['']
            # elif tvar=='SW_flux_net_at_model_top':
            #    data = ds['SW_flux_up_at_model_top'] + ds['']
            else:
               data = ds[tvar]
            if obs_flag[c]: 
               if 'latitude'  in data.dims: data = data.rename({'latitude':'lat'})
               if 'longitude' in data.dims: data = data.rename({'longitude':'lon'})
               data = data.stack(ncol=('lat','lon'))
            #-------------------------------------------------------------------
            if obs_flag[c]:
               if var[v]=='U_at_850hPa': data = data.sel(level=850)
               if var[v]=='V_at_850hPa': data = data.sel(level=850)
               if var[v]=='U_at_500hPa': data = data.sel(level=500)
               if var[v]=='V_at_500hPa': data = data.sel(level=500)
            else:
               if var[v]=='horiz_winds_at_model_bot_u': data = data.isel(dim2=0)
               if var[v]=='horiz_winds_at_model_bot_v': data = data.isel(dim2=1)
            #-------------------------------------------------------------------
            if 'lev' in data.dims and lev_list[v] is not None: 
               if lev_list[v]<0: data = data.isel({'lev':np.absolute(lev_list[v])})
            #-------------------------------------------------------------------
            # convert units
            # if var[v]=='precip_total_surf_mass_flux':
            #    data = data*86400*1e3 # m/s to mm/day
            if unit_fac_list[v] is not None : data = data * unit_fac_list[v]
            #-------------------------------------------------------------------
            # print('-'*80)
            # print()
            # print(data)
            # print()
            # print('-'*80)
            # exit('stopping to sanity check the data')
            #-------------------------------------------------------------------
            # build time coodinate
            time_tmp = data['time.dayofyear'].values   \
                      +data['time.hour'].values/24     \
                      +data['time.minute'].values/24/60
            # time_tmp = data['time.year'].values*365    \
            #           +data['time.dayofyear'].values   \
            #           +data['time.hour'].values/24     \
            #           +data['time.minute'].values/24/60
            #-------------------------------------------------------------------
            # build datetime coordinate (for optional date x-axis)
            date_tmp = np.array([ datetime.datetime(int(yr),int(mo),int(dd),int(hr),int(mn))
                                  for yr,mo,dd,hr,mn in zip( data['time.year'].values,
                                                             data['time.month'].values,
                                                             data['time.day'].values,
                                                             data['time.hour'].values,
                                                             data['time.minute'].values ) ])
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
            # apply mask for spatial subset
            if 'lat1' in locals():
               print(' '*8+f'Applying spatial mask...')
               if obs_flag[c]:
                  data = data.where(obs_mask,drop=True)
               else:
                  data = data.where(sim_mask,drop=True)
            #-------------------------------------------------------------------
            # print('-'*80)
            # print()
            # print(data)
            # print()
            # print('-'*80)
            # exit('stopping to sanity check the data')
            #-------------------------------------------------------------------
            # print(' '*8+f'Calculating spatial mean...')
            # gbl_mean = ( (data*area).sum(dim=xy_dims) / area.sum(dim=xy_dims) )
            gbl_mean_tmp = (data*area.values).sum(dim=xy_dims) / np.sum(area.values)
            # np.ma.masked_invalid(var_out.values)
            gbl_mean_tmp['time'] = time_tmp
            # gbl_mean_tmp = xr.DataArray( np.zeros( (ntime,nbins) ), \
            #                         coords=[('time',time_tmp),('lon',lon_bins)], \
            #                         dims=['time','lon'] )
            #-------------------------------------------------------------------
            # if var[v]=='LiqWaterPath':
            #    if gbl_mean_tmp[0]>90:
            #       print(f'DEBUG - case: {case[c]}   gbl_mean_tmp[0]: {gbl_mean_tmp[0]}')
            #-------------------------------------------------------------------
            # if not obs_flag[c]:
            #    # first hourly data point being zero?
            #    # just set it to the next point and get on with it!
            #    gbl_mean_tmp[0] = gbl_mean_tmp[1]
            #-------------------------------------------------------------------------
            # dy = data['time.dayofyear'].values
            # hr = data['time.hour'].values
            # mn = data['time.minute'].values
            # print()
            # for t in range(len(time_tmp)): 
            #    # print(f'  {time_tmp[t]:8.4f}  {gbl_mean_tmp[t].values:8.2f}')
            #    print(f'  {dy[t]:03d}  {hr[t]:02d}:{mn[t]:02d} {gbl_mean_tmp[t].values:20.6f}')
            # print()
            # for t in range(len(time_tmp)): 
            #    tmp = (data[t,:]*area.values).sum(dim=xy_dims) / np.sum(area.values)
            #    print(f'  {dy[t]:03d}  {hr[t]:02d}:{mn[t]:02d} {tmp.values:20.6f}')
            # print()
            # # for t in range(len(time_tmp)): 
            # #    print(f'  {dy[t]:03d}  {hr[t]:02d}:{mn[t]:02d} {data[t,:].sum(dim=xy_dims).values:20.6f}')
            # exit()
            #-------------------------------------------------------------------
            if f==0:
               gbl_mean = gbl_mean_tmp
               time     = time_tmp
               date     = date_tmp
            else:
               gbl_mean = xr.concat([gbl_mean,gbl_mean_tmp], dim='time')
               time     = np.append(time,time_tmp)
               date     = np.append(date,date_tmp)
         # exit()
         #-------------------------------------------------------------------
         # convert units
         if var[v]=='precip':
            if obs_flag[c]:
               gbl_mean = gbl_mean*24. # mm/hr to mm/day
            else:
               gbl_mean = (gbl_mean/100)*86400 # kg/m2 to mm/day using dtime=100 (ne1024)
         #----------------------------------------------------------------------
         # print some summary stats after averaging globally
         if print_stats: hapy.print_stat(gbl_mean,name=tvar,compact=True,indent=' '*6)
         #----------------------------------------------------------------------
         # Write to file
         print(' '*6+f'Writing data to file: {tmp_file}')
         gbl_mean.name = tvar
         tmp_ds = xr.Dataset()
         # tmp_ds['time'] = time
         tmp_ds[tvar]   = gbl_mean
         tmp_ds['date'] = date
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         print(' '*6+f'Reading pre-calculated data from file: {hapy.tclr.MAGENTA}{tmp_file}{hapy.tclr.END}')
         tmp_ds = xr.open_dataset( tmp_file )
         # time_tmp = tmp_ds['time'].values
         gbl_mean = tmp_ds[tvar]
         date     = tmp_ds['date'].values
      #-------------------------------------------------------------------------
      time_tmp = gbl_mean['time'].values
      # time_tmp = time_tmp/(24*60)
      # time_tmp = time_tmp - time_tmp[0]# + 20
      # print(time_tmp)
      time_list.append( time_tmp )
      data_list.append( gbl_mean.values )
      date_list.append( date )
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

   # print()
   # print(f'  data_min: {data_min}')
   # print(f'  data_max: {data_max}')
   # dy = (data_max-data_min); print(); print(f'  dy: {dy}')

   # use fixed axis width for some variables
   if var[v] in ['T_mid','T_2m']:
      dy_fixed = 20
      if ( (data_max-data_min)!=dy_fixed ):
         delta = dy_fixed - (data_max-data_min)
         data_max = data_max + delta/2
         data_min = data_min - delta/2

   # print()
   # print(f'  data_min: {data_min}')
   # print(f'  data_max: {data_max}')
   # dy = (data_max-data_min); print(); print(f'  dy: {dy}')

   ax = axs[v]
   #----------------------------------------------------------------------------
   # select x-axis values: actual dates or day-of-year
   if use_date_axis:
      xvals_list = date_list
      time_min = np.min([np.min(d) for d in date_list])
      time_max = np.max([np.max(d) for d in date_list])
   else:
      xvals_list = time_list
   ax.set_ylim(data_min,data_max)
   ax.set_xlim(time_min,time_max)
   # if 'nday_plot' in locals():
   #    ax.set_xlim(time_min,nday_plot)
   # else:
   #    ax.set_xlim(time_min,time_max)
   if use_date_axis:
      ax.set_xlabel('Date',fontsize=label_font_size)
   else:
      ax.set_xlabel('Time [days]',fontsize=label_font_size)
   ax.set_ylabel(f'{var_str[v]}',fontsize=label_font_size)
   # ax.set_ylabel(f'[{var_unit[v]}]',fontsize=label_font_size)
   #----------------------------------------------------------------------------
   # Calculate anomalies from a reference value
   if plot_as_anomaly:
      for c in range(num_case):
         # calculate ref value as mean of last 2 days
         sample_per_day = int(24*60/15) # 15 minute sampling fro SCREAM
         nsample = 2*sample_per_day
         ref_value = np.mean(data_list[c][-1*nsample:])
         # subtract reference value to collapse curves
         data_list[c] = data_list[c] - ref_value
   #----------------------------------------------------------------------------
   # data_min = np.min([np.min(d) for d in data_list])
   # data_max = np.max([np.max(d) for d in data_list])
   # res.trYMinF = data_min
   # res.trYMaxF = data_max
   for c in range(num_case):
      #-------------------------------------------------------------------------
      tlabel = case_name[c] if case_name[c]!='' else case[c]
      ax.plot( xvals_list[c], data_list[c],
               color=clr[c], linestyle=dash_to_linestyle(dsh[c]),
               linewidth=line_thickness, label=tlabel )
      #-------------------------------------------------------------------------
      # add line to indicate reference value
      if plot_as_anomaly:
         ref_value = np.mean(data_list[c][-1*nsample:])
         ref_array = np.ones(time_list[c].shape) * ref_value
         ax.plot( xvals_list[c], ref_array,
                  color=clr[c], linestyle='--', linewidth=ref_line_thickness )
   #----------------------------------------------------------------------------
   # subtitles: region string (left) and variable name (right)
   # ax.set_title( get_ctr_str(), loc='left',  fontsize=label_font_size )
   # ax.set_title( var_str[v],    loc='right', fontsize=title_font_size )
   ax.set_title( "Global Average Spin Up and PPE Members", fontsize=title_font_size )
   ax.tick_params(labelsize=label_font_size)
   ax.grid(True,alpha=0.3)
   #----------------------------------------------------------------------------
   # format x-axis as dates (mm-dd-yyyy)
   if use_date_axis:
      ax.xaxis.set_major_formatter( mdates.DateFormatter('%m-%d-%Y') )
      for lbl in ax.get_xticklabels():
         lbl.set_rotation(30)
         lbl.set_ha('right')
   if any(n!='' for n in case_name):
      ax.legend(fontsize=label_font_size*0.8)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# hide any unused panels
for ax in axs[num_var:]: ax.set_visible(False)

fig.tight_layout()

fig_path = f'{fig_file}.{fig_type}'
os.makedirs(os.path.dirname(fig_path),exist_ok=True)
fig.savefig( fig_path, dpi=150, bbox_inches='tight' )
plt.close(fig)
print(f'\n{fig_path}\n')
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------