import os, glob, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string, sys, cmocean
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
#-------------------------------------------------------------------------------
case_dir,case_sub = [],[]
case,name = [],[]
htype_list = []
comp_list = []
clr,dsh,mrk = [],[],[]
def add_case(case_in,n=None,comp=None,htype=None,p=None,s=None,d=0,c='black',m=0):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); name.append(n); 
   comp_list.append(comp); htype_list.append(htype)
   case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
#-------------------------------------------------------------------------------
var_list = []
eam_var_list = []
eamxx_var_list = []
lev_list = []
unit_list = []
def add_var(var_name,eam_var=None,eamxx_var=None,var_str=None,lev=None,unit=None): 
   var_list.append(var_name)
   eam_var_list.append(eam_var)
   eamxx_var_list.append(eamxx_var)
   lev_list.append(lev)
   unit_list.append(unit)
#-------------------------------------------------------------------------------
tmp_path_ne1024 = '/global/cfs/cdirs/e3smdata/simulations/scream-decadal'
# tmp_sub_ne1024  = 'data_remap_90x180/output.scream.decadal.3hourlyINST_ne30pg2'
# tmp_sub_ne1024  = 'data_remap_90x180/output.scream.decadal.6hourlyAVG_ne30pg2'
tmp_path_hst_v3 = '/global/cfs/cdirs/m3312/whannah/e3smv3_amip'
# tmp_path_qbo_bm = '/global/cfs/cdirs/m4310/data/sims'
tmp_path_qbo_bm = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu'

# add_case('Obs')
add_case('v3.LR.amip_0101.QBObenchmark.20241008', n='EAMv3 AMIP',      comp='eam',  htype='eam.h3', p=tmp_path_qbo_bm, s='data_remap_90x180')
# add_case('v3.LR.amip_0201',                       n='v3.LR.amip_0201', comp='eam',  htype='eam.h1', p=tmp_path_hst_v3, s='data_remap_90x180')
# add_case('decadal-production-run6', n='SCREAM ne1024',   comp='eamxx',c='blue', htype='output.scream.decadal.3hourlyINST_ne30pg2',p=tmp_path_ne1024,s='data_remap_90x180/output.scream.decadal.3hourlyINST_ne30pg2')
# add_case('decadal-production-run6', n='SCREAM ne1024',   comp='eamxx',htype='output.scream.decadal.6hourlyAVG_ne30pg2',p=tmp_path_ne1024,s='data_remap_90x180/output.scream.decadal.6hourlyAVG_ne30pg2')


'''
ERA5 hourly data:
/global/cfs/projectdirs/m3522/cmip6/ERA5/e5.oper.an.pl/201111/e5.oper.an.pl.128_130_t.ll025sc.2011111900_2011111923.nc
/global/cfs/cdirs/m3312/whannah/ERA5/daily
'''
#-------------------------------------------------------------------------------

# add_var('Precipitation',    eam_var='PRECT',eamxx_var='precip_liq_surf_mass_flux',unit='mm/day')
# add_var('OLR',              eam_var='FLUT', eamxx_var='LW_flux_up_at_model_top',  unit='W/m2')
add_var('850 mb Zonal wind',eam_var='U850', eamxx_var='U_at_850hPa',              unit='m/s')
# add_var('200 mb Zonal wind',eam_var='U200', eamxx_var='U_at_200hPa',              unit='m/s')

num_plot_col = len(var_list)

#-------------------------------------------------------------------------------

fig_type,fig_file = 'png',f'figs/FXX-diurnal-composite-hovmoller-v1'

tmp_file_head = 'data/diurnal-composite-hovmoller-v1'


# yr1,yr2 = 1975,2020
yr1,yr2 = 1995,2004
# yr1,yr2 = 1995,1995

recalculate = True

lat1,lat2 = -15,15

var_x_case  = True

# num_plot_col = 2

use_common_label_bar = True

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
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_var = len(case),len(var_list)

subtitle_font_height = 0.01

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

plot = [None]*(num_var)

res = hs.res_contour_fill()
res.vpHeightF = 0.3
res.lbLabelFontHeightF     = 0.022
res.tiXAxisFontHeightF     = 0.025
res.tiYAxisFontHeightF     = 0.025
res.tmYLLabelFontHeightF   = 0.02
res.tmXBLabelFontHeightF   = 0.02
res.tiXAxisString = ''
res.tiYAxisString = ''

if use_common_label_bar:  res.lbLabelBarOn = False

# res.trYMinF,res.trYMaxF = 0.0,0.5
res.trXMinF,res.trXMaxF = 0,100


lres = hs.res_xy()
lres.xyLineThicknessF = 1
lres.xyLineColor      = 'black'
lres.xyDashPattern    = 0

#---------------------------------------------------------------------------------------------------  
#---------------------------------------------------------------------------------------------------
def apply_mask(data):
   global lat1, lat2, lon1, lon2
   mask = xr.DataArray( np.ones([len(data['lat']),len(data['lon'])],dtype=bool), dims=('lat','lon') )
   if 'lat1' in globals(): mask = mask & (data['lat']>=lat1) & (data['lat']<=lat2)
   if 'lon1' in globals(): mask = mask & (data['lon']>=lon1) & (data['lon']<=lon2)
   data = data.where( mask, drop=True)
   return data
#---------------------------------------------------------------------------------------------------
# def get_obs_name(case,var):
#    if case[c]=='Obs' and var[v] in ['OLR','FLNT','FLUT']: obs_var = 'olr'  ; obs_name = 'NOAA'
#    if case[c]=='Obs' and var[v] in ['PRECT']            : obs_var = 'PRECT'; obs_name = 'IMERG'
#    if case[c]=='Obs' and var[v] in ['U850']             : obs_var = 'U850' ; obs_name = 'ERA5'
#    return obs_name, obs_var
#---------------------------------------------------------------------------------------------------

for v in range(num_var):
   print(' '*4+f'var: '+hc.tcolor.GREEN+f'{var_list[v]}'+hc.tcolor.ENDC)

   data_list = []
   time_list = []
   lon_list = []

   for c in range(num_case):
      tcase = case[c]
      # if case[c]=='Obs': tcase,tvar = get_obs_name(case,var)
      print(' '*2+f'case: '+hc.tcolor.CYAN+f'{tcase}'+hc.tcolor.ENDC)
      #-------------------------------------------------------------------------
      tvar = None
      if comp_list[c]=='eam'  : tvar = eam_var_list[v]
      if comp_list[c]=='eamxx': tvar = eamxx_var_list[v]
      if tvar is None: raise ValueError(f'ERROR - variable name not found?!?!?  comp: {comp_list[c]} ')
      #-------------------------------------------------------------------------
      tmp_file = f'{tmp_file_head}.{tcase}.{var_list[v]}'
      tmp_file+= f'.yr_{yr1}_{yr2}'
      tmp_file+= f'.lat_{lat1}_{lat2}'
      # tmp_file+= f'.lat_{lon1}_{lon2}'
      tmp_file+= f'.nc'

      ##########################################################################
      ##########################################################################
      if recalculate :
         #----------------------------------------------------------------------
         # identify files
         # if tcase=='NOAA':
         #    file_list = ['/global/cfs/cdirs/m3312/whannah/obs_data/OLR/olr.day.mean.nc']
         # elif tcase=='IMERG':
         #    file_list = ['/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/IMERG_Daily/PRECT_200101_202012.nc']
         if tcase=='ERA5':
            file_list = ['/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/ERA5_Daily/U850_198001_202212.nc']
         else:
            file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/*{htype_list[c]}*'
            file_list = sorted(glob.glob(file_path))
         if file_list==[]:
            print()
            print(f'file_path: {file_path}')
            print(f'file_list: {file_list}')
            exit('ERROR: no files found')
         #----------------------------------------------------------------------
         # print()
         # for f in file_list: print(f)
         # print()
         #----------------------------------------------------------------------
         # read the data
         ds = xr.open_mfdataset(file_list)
         data = ds[tvar]
         #----------------------------------------------------------------------
         if comp_list[c]=='eamxx':
            if eam_var_list[v]=='PRECT': data = data + ds['precip_ice_surf_mass_flux']
         #----------------------------------------------------------------------
         # subset the data based on year
         data = data.where( data['time.year']>=yr1, drop=True)
         data = data.where( data['time.year']<=yr2, drop=True)
         #----------------------------------------------------------------------
         # change unit of time coordinate
         data['time'] = data['time.year']*365 + data['time.day'] + data['time.hour']/24
         data['time'] = data['time'] - data['time'][0]
         #----------------------------------------------------------------------
         # reduce to equatorial region
         data = apply_mask(data).mean(dim='lat')
         #----------------------------------------------------------------------
         # adjust units
         if eam_var_list[v]=='PRECT' and comp_list[c] in ['eam','eamxx']: data = data*86400.*1e3
         #----------------------------------------------------------------------
         # # Convert to daily mean
         # data = data.resample(time='D').mean(dim='time')
         #----------------------------------------------------------------------
         hc.print_stat(data,name=var_list[v],compact=True,indent=' '*6)
         #----------------------------------------------------------------------
         # diurnal composite
         comp = data[:8,:]*0.
         for n in range(4):
            comp[n,:] = data[n:-4+n:4,:].mean(dim='time') - data.mean(dim='time')
            comp[n+4,:] = comp[n,:]
         data = comp
         #----------------------------------------------------------------------
         hc.print_stat(data,name=var_list[v],compact=True,indent=' '*6)
         # exit()
         #----------------------------------------------------------------------
         # print(); print(data)
         # print()
         # exit()
         #----------------------------------------------------------------------
         # Write to file 
      #    print(' '*6+f'writing to file - {tmp_file}')
      #    ds = xr.Dataset()
      #    ds[eam_var_list[v]] = data
      #    ds.to_netcdf(path=tmp_file,mode='w')
      # else:
      #    print(f'      loading pre-calculated spectra... {tmp_file}')
      #    ds = xr.open_mfdataset( tmp_file )
      #    data = ds[eam_var_list[v]]


      data_list.append(data.values)
      time_list.append(data['time'].values)
      lon_list.append(data['lon'].values)

   #----------------------------------------------------------------------------
   # calculate limits for common color bar
   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])
   #----------------------------------------------------------------------------
   # Create plot
   for c in range(num_case):

      tres = copy.deepcopy(res)

      ip = v*num_case+c if var_x_case else c*num_var+v

      tres.sfXArray = lon_list[c]
      tres.sfYArray = time_list[c]

      # print()
      # print(data_list[c].shape)
      # for t in range(len(time_list[c])):
      #    for x in range(len(lon_list[c])):
      #       print(type(data_list[c][t,x]))
      # print()
      # exit()

      plot[ip] = ngl.contour(wks, data_list[c], tres)
      
      hs.set_subtitles(wks, plot[ip], 
                       left_string=name[c], 
                       center_string='', 
                       right_string=var_list[v], 
                       font_height=subtitle_font_height)

# #---------------------------------------------------------------------------------------------------
# # Finalize plot
# #---------------------------------------------------------------------------------------------------

layout = [num_var,num_case] if var_x_case else [num_case,num_var]

# layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

pnl_res = hs.setres_panel()
# pnl_res.nglPanelYWhiteSpacePercent = 5
# pnl_res.nglPanelXWhiteSpacePercent = 5
# pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
# pnl_res.nglPanelFigureStringsJust        = "TopLeft"
# pnl_res.nglPanelFigureStringsFontHeightF = 0.015

ngl.panel(wks,plot,layout,pnl_res)
hc.trim_png(fig_file)

ngl.end()
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
