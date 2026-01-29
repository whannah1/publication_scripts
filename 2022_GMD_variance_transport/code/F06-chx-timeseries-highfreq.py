# plot the fractional occurence of all possible patterns
import os, ngl, sys, numba, copy, xarray as xr, numpy as np, re, string
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import pg_checkerboard_utilities as pg
import cmocean

case,name,clr,dsh,pat = [],[],[],[],[]
var,var_name,lev_list = [],[],[]
subset_min_length_list = []
def add_case(case_in,sml=4,n='',c='black',d=0,p=0):
   global name,case
   case.append(case_in); name.append(n); 
   clr.append(c); dsh.append(d); pat.append(p)
   subset_min_length_list.append(sml)

def add_var(v,n=None,lev=-1): var.append(v); var_name.append(n); lev_list.append(lev)
#-------------------------------------------------------------------------------
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.VT_1.00',sml=4,n='E3SM-MMF+FVT1',c='blue1',d=0)
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.VT_2.00',sml=4,n='E3SM-MMF+FVT2',c='blue2',d=0)
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.VT_4.00',sml=4,n='E3SM-MMF+FVT4',c='blue3',d=0)
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.VT_8.00',sml=4,n='E3SM-MMF+FVT8',c='blue4',d=0)
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.VT_0.00',sml=4,n='E3SM-MMF+BVT' ,c='green',d=0)

add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.00',n='E3SM-MMF+BVT' ,c='green',d=0)
add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_1.00',n='E3SM-MMF+FVT' ,c='blue',d=0)
add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.00',           n='E3SM-MMF'     ,c='red' ,d=0)

# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.TEST_QU.00',n='MMF+BVT only Q+U', c='magenta')
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.TEST_TQ.00',n='MMF+BVT only T+Q', c='purple')
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.TEST_TU.00',n='MMF+BVT only T+U', c='pink')
# # add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.TEST_T.00', n='MMF+BVT only T',   c='magenta')
# # add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.TEST_Q.00', n='MMF+BVT only Q',   c='purple')
# # add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.TEST_U.00', n='MMF+BVT only U',   c='pink')
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.00',        n='MMF+BVT',          c='green')
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.00',                   n='MMF',              c='red')
#-------------------------------------------------------------------------------

var = []
add_var('TGCLDLWP',n='Liq Water Path')
add_var('PRECT',n='Precipitation')


fig_type = 'png'
fig_file = 'figs/F06-chx-timeseries-highfreq'
tmp_file_head = 'data/timeseries-highfreq.v1'

# scratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl/chx-detection-data'
scratch = '/gpfs/alpine/scratch/hannah6/cli115/analysis_files/vt_validation/chx-detection-data'
chx_data_file_head=f'{scratch}/occurence-partial-chx-over-time'

recalculate = False

ocean_only = True

add_obs = True

var_x_case           = True    # controls plot panel arrangement
print_stats          = True    #
verbose              = True    #

lat1,lat2 = -60,60
# lat1,lat2,lon1,lon2 = 0,30,140,180+40

# ilat,ilon,idx = -20,360-115,20
# lat1,lat2 = ilat-idx,ilat+idx
# lon1,lon2 = ilon-idx, ilon+idx

# subset_min_length = 4       # threshold for partial checkerboard when combine_mode=1

#---------------------------------------------------------------------------------------------------
# generate sets of possible neighbor states
#---------------------------------------------------------------------------------------------------

### only pure checkboard
# tmp_file='data/occurrence-chx-only'; rotate_sets=False; sets=pg.chx_only_sets

### find index of pure checkerboard
# chx_idx = []
# for s,tset in enumerate(sets):
#    for chx in chx_sets:
#       if all(tset==chx): chx_idx.append(s)

# (num_set,set_len) = sets.shape
# set_coord,nn_coord = np.arange(num_set),np.arange(set_len)
# sets.assign_coords(coords={'set':set_coord,'neighbors':nn_coord})

# set_labels = pg.get_set_labels(sets)

set_labels = ['no checkerboard','partial checkerboard']

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_var = len(case),len(var)

wkres = ngl.Resources()
# npix=2**13; wkres.wkWidth,wkres.wkHeight=npix,npix # use this for plotting all patterns w/ rotation
wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*(num_var)

res = hs.res_xy()
# res.tmXBOn = False
# res.tmYLOn = False

res.tiXAxisFontHeightF = 0.02
res.tiYAxisFontHeightF = 0.02

res.tiXAxisString = 'Time [days]'
res.tiYAxisString = 'Fractional Occurrence of Partial Checkerboard'

# res.xyXStyle = 'Log'
# res.xyYStyle = 'Log'

# res.xyMonoLineThickness = True
# res.xyLineThicknessF = 20.


scrip_file_path = 'scrip_files/ne30pg2_scrip.nc'
scripfile = xr.open_dataset(scrip_file_path)
lat = scripfile['grid_center_lat'].rename({'grid_size':'ncol'})
lon = scripfile['grid_center_lon'].rename({'grid_size':'ncol'})

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_alt_name(case,var):
   tname = case
   if 'OBS' in case and var=='TGCLDLWP': tname = 'MAC'
   if 'OBS' in case and var=='PRECT'   : tname = 'GPM'
   alt_name = case
   if 'OBS' in case: alt_name = alt_name.replace('OBS',tname)
   return alt_name

def get_data_file_name(case,var,sml):
   case_tmp = get_alt_name(case,var)
   data_file = f'{chx_data_file_head}.high-freq.{case_tmp}.{var}.sml_{sml}.nc'
   return data_file

def print_file_status(tfile):
   msg = tfile
   msg = msg.replace('/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl/chx-detection-data','<...>')
   msg = f'  {msg:140}  '
   if os.path.exists(tfile): 
      print(msg+hc.tcolor.GREEN+'OK'+hc.tcolor.ENDC)
   else:
      print(msg+hc.tcolor.RED+'MISSING!'+hc.tcolor.ENDC)

def get_tmp_file_name(case,var,sml):
   case_tmp = get_alt_name(case,var)
   tmp_tmp_file = f'{tmp_file_head}.{case_tmp}.{var}'
   if 'lat1' in globals(): tmp_tmp_file+= f'.lat1_{lat1}.lat2_{lat2}'
   if 'lon1' in globals(): tmp_tmp_file+= f'.lon1_{lat1}.lon2_{lat2}'
   tmp_tmp_file+= f'.sml_{sml}.nc'
   return tmp_tmp_file

def get_mask(ds):
   mask = xr.DataArray( np.ones(len(ds['ncol']),dtype=bool), coords=[ds['ncol']] )
   if 'lat1' in locals(): mask = mask & (lat>=lat1) & (lat<=lat2)
   if 'lon1' in locals(): mask = mask & (lon>=lon1) & (lon<=lon2)
   return mask
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# check if any files are missing
print()
for i in range(2):
   print()
   for v in range(num_var):
      for c in range(num_case):
         if i==0: print_file_status( get_data_file_name(case[c],var[v],subset_min_length_list[c]) )
         if i==1: print_file_status( get_tmp_file_name(case[c],var[v],subset_min_length_list[c]) )
print()

# load data for ocean mask
if ocean_only:    
   landfrac_ds = xr.open_dataset('data/E3SM_landfrac_ne30pg2.nc').isel(time=0)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
list_cnt = 0
# plot = [None]*(num_var*num_case*num_set)
plot = [None]*(num_var)
# plot = [None]

for v in range(num_var):
   print('\n    var: '+hc.tcolor.GREEN+f'{var[v]:10}'+hc.tcolor.ENDC)

   #----------------------------------------------------------------------------
   # load obs data
   #----------------------------------------------------------------------------
   if add_obs:
      if var[v]=='TGCLDLWP': alt_case = f'MAC-PG'
      if var[v]=='PRECT'   : alt_case = f'GPM-PG'

      chx_data_file = f'data/chx-occurrence.daily.{alt_case}.{var[v]}.nc'
      # chx_data_file = get_data_file_name(case[c],var[v],subset_min_length_list[c])
      cnt_ds = xr.open_dataset( chx_data_file )
      mask = get_mask(cnt_ds)
      cnt_ds = cnt_ds.where(mask,drop=True)
      if ocean_only:
         landfrac = landfrac_ds['LANDFRAC']
         landfrac = landfrac.where(mask,drop=True)    
         ocn_mask = xr.DataArray( np.ones(len(cnt_ds['ncol']),dtype=bool), coords=[cnt_ds['ncol']] )
         ocn_mask = ocn_mask & (landfrac.values<0.5)
         cnt_ds = cnt_ds.where(ocn_mask,drop=True)
      obs_pchx = ( cnt_ds['cnt'].isel(set=1) / cnt_ds['num_valid'] ).mean().values
      print(f'    obs climatology value: {obs_pchx}')
   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------

   data_list,time_list = [],[]

   for c in range(num_case):

      chx_data_file = get_data_file_name(case[c],var[v],subset_min_length_list[c])
      tmp_file      = get_tmp_file_name(case[c],var[v],subset_min_length_list[c])
      
      # print('\n  case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC+'    '+chx_data_file)
      print('\n      case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)

      if recalculate:
         #----------------------------------------------------------------------
         # Load the teporally resolved chx identification data
         #----------------------------------------------------------------------
         cnt_ds = xr.open_dataset( chx_data_file )

         # Apply regional mask
         mask = get_mask(cnt_ds)
         cnt_ds = cnt_ds['cnt'].where(mask,drop=True).to_dataset()

         # just use the partial chx flag
         tmp_data = cnt_ds['cnt'][1,:,:]

         #----------------------------------------------------------------------
         # Ocean mask
         #----------------------------------------------------------------------
         if ocean_only:
            landfrac = landfrac_ds['LANDFRAC']
            landfrac = landfrac.where(mask,drop=True)    
            ocn_mask = xr.DataArray( np.ones(len(cnt_ds['ncol']),dtype=bool), coords=[cnt_ds['ncol']] )
            ocn_mask = ocn_mask & (landfrac.values<0.5)
            # cnt_ds['cnt'] = cnt_ds['cnt'].where(ocn_mask,drop=True)
            tmp_data = tmp_data.where(ocn_mask,drop=True)
         #----------------------------------------------------------------------
         # misc data prep
         #----------------------------------------------------------------------
         list_cnt += 1

         ### average over all columns
         tmp_data.load()
         tmp_data = tmp_data.mean(dim='ncol')

         ### Make time start at zero
         tmp_time = ( cnt_ds['time'] - cnt_ds['time'][0] ).astype('float') / 86400e9

         #----------------------------------------------------------------------
         # write to file
         #----------------------------------------------------------------------
         tmp_ds = xr.Dataset()
         tmp_ds['time'] = ( ('time'), tmp_time )
         tmp_ds['data'] = ( ('time'), tmp_data )
         print(f'      writing to file: {tmp_file}')
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         tmp_ds = xr.open_dataset( tmp_file )
         tmp_data = tmp_ds['data']
         tmp_time = tmp_ds['time']

      ### Convert to hourly
      # print(); print(tmp_data); print()
      # tmp_data = tmp_data.resample(time='H').mean(dim='time')
      # tmp_time = tmp_time.resample(time='H').mean(dim='time')

      if print_stats: hc.print_stat(tmp_data,name='avg pchx occurrence',stat='nax',indent=' '*6,compact=True)

      data_list.append(tmp_data.values)
      time_list.append(tmp_time.values)
         
   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   ip = v

   tres = copy.deepcopy(res)
   tres.xyLineColors = clr
   tres.xyDashPatterns = dsh

   plot[ip] = ngl.xy(wks, np.stack(time_list), np.stack(data_list), tres) 

   if add_obs:
      ox = np.array([np.min(time_list),np.max(time_list)])
      oy = np.array([obs_pchx,obs_pchx])
      tres.xyLineColors = 'black'
      tres.xyDashPatterns = 1
      ngl.overlay( plot[ip], ngl.xy(wks, ox, oy, tres)  )
   

   var_name = var[v]
   if var_name=='TGCLDLWP': var_name = 'Liq Water Path'
   if var_name=='PRECT':    var_name = 'Precipitation'
   hs.set_subtitles(wks, plot[ip], var_name, '', '', font_height=0.020)

#---------------------------------------------------------------------------------------------------
# Add legend
#---------------------------------------------------------------------------------------------------
if add_obs:
   name.append('OBS')
   clr.append('black')

lgres = ngl.Resources()
lgres.vpWidthF, lgres.vpHeightF  = 0.09, 0.12
if num_case>4: lgres.vpHeightF  = 0.18
lgres.lgLabelFontHeightF = 0.012
lgres.lgLineThicknessF   = 20
lgres.lgMonoLineColor,lgres.lgLineColors  = False, clr
lgres.lgMonoDashIndex,lgres.lgDashIndexes = True, 0
lgres.lgLabelJust    = 'CenterLeft'

lname = [f' {n}' for n in name]
for (n,l) in enumerate(lname):
   if 'OBS' in lname[n]: 
      lname[n] = lname[n].replace('OBS','MAC')+' climatology'
xpos = 0.25
ypos = 0.43
pid = ngl.legend_ndc(wks, len(name), lname, xpos, ypos, lgres)

lname = [f' {n}' for n in name]
for (n,l) in enumerate(lname):
   if 'OBS' in lname[n]: 
      lname[n] = lname[n].replace('OBS','IMERG')+' climatology'
xpos = 0.75
pid = ngl.legend_ndc(wks, len(name), lname, xpos, ypos, lgres)
#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
layout = [1,num_var]
# layout = [1,len(plot)]

pnl_res = hs.setres_panel()
# pnl_res.nglPanelYWhiteSpacePercent = 10
# pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
# pnl_res.nglPanelFigureStringsJust        = "TopLeft"
# pnl_res.nglPanelFigureStringsFontHeightF = 0.015   

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------