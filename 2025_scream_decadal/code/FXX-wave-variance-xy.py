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
   # if comp  is None: raise ValueError(f'ERROR - add_case: comp argument cannot be None')
   # if htype is None: raise ValueError(f'ERROR - add_case: htype argument cannot be None')
   # if p     is None: raise ValueError(f'ERROR - add_case: p argument cannot be None')
   # if s     is None: raise ValueError(f'ERROR - add_case: s argument cannot be None')
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
tmp_path_ne256  = '/global/cfs/cdirs/e3smdata/simulations/scream-decadal-ne256'
# tmp_sub_ne1024  = 'data_remap_90x180/output.scream.decadal.3hourlyINST_ne30pg2'
# tmp_sub_ne1024  = 'data_remap_90x180/output.scream.decadal.6hourlyAVG_ne30pg2'
tmp_path_hst_v3 = '/global/cfs/cdirs/m3312/whannah/e3smv3_amip'
# tmp_path_qbo_bm = '/global/cfs/cdirs/m4310/data/sims'
tmp_path_qbo_bm = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu/'

add_case('Obs')
add_case('v3.LR.amip_0101.QBObenchmark.20241008',                                        n='EAMv3 AMIP',          comp='eam',  d=0,c='red',  htype='eam.h2',                                   p=tmp_path_qbo_bm,s='data_remap_73x144')
# add_case('decadal-production-run6',                                                      n='SCREAM 3km',          comp='eamxx',d=0,c='blue', htype='output.scream.decadal.6hourlyINST_ne30pg2',p=tmp_path_ne1024,s='data_remap_73x144/output.scream.decadal.6hourlyINST_ne30pg2')
add_case('decadal-production-run6',                                                      n='SCREAM 3km',          comp='eamxx',d=0,c='blue', htype='output.scream.decadal.6hourlyAVG_ne30pg2', p=tmp_path_ne1024,s='data_remap_73x144/output.scream.decadal.6hourlyAVG_ne30pg2')
add_case('ne256pg2_ne256pg2.F20TR-SCREAMv1.May-12.with.rain.frac.n0128',                 n='SCREAM 13km control', comp='eamxx',d=0,c='green',htype='6ha_ne30pg2.AVERAGE.nhours_x6',            p=tmp_path_ne256, s='data_remap_73x144/6ha_ne30pg2.AVERAGE.nhours_x6')
add_case('ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1', n='SCREAM 13km tuned',   comp='eamxx',d=1,c='green',htype='6ha_ne30pg2.AVERAGE.nhours_x6',            p=tmp_path_ne256, s='data_remap_73x144/6ha_ne30pg2.AVERAGE.nhours_x6')
### add_case('v3.LR.amip_0201',                       n='v3.LR.amip_0201', comp='eam',  htype='eam.h2', p=tmp_path_hst_v3, s='data_remap_90x180')

'''
ERA5 hourly data:
/global/cfs/projectdirs/m3522/cmip6/ERA5/e5.oper.an.pl/201111/e5.oper.an.pl.128_130_t.ll025sc.2011111900_2011111923.nc
/global/cfs/cdirs/m3312/whannah/ERA5/daily
'''
#-------------------------------------------------------------------------------

add_var('OLR', eam_var='FLUT', eamxx_var='LW_flux_up_at_model_top', unit='W/m2')
add_var('U850',eam_var='U850', eamxx_var='U_at_850hPa',             unit='m/s')
# add_var('U',eam_var='U', eamxx_var='U',     unit='m/s')
# add_var('T',eam_var='T', eamxx_var='T_mid', unit='K')

num_plot_col = 1

#-------------------------------------------------------------------------------

wave_type_list = []
wave_type_list.append('MJO')
wave_type_list.append('Kelvin')
wave_type_list.append('Rossby')
# wave_type_list.append('MRG')

fig_type,fig_file = 'png',f'figs/FXX-wave-variance-xy'
# fig_file_diff = f'{fig_file}.diff'

# tmp_file_head = 'data/kelvin-variance-map'

# yr1,yr2 = 1975,2020
yr1,yr2 = 1995,2004
# yr1,yr2 = 1995,1995
# yr1,yr2 = 1995,1999

# recalculate = True

var_x_case  = False

# num_plot_col = 2

use_common_label_bar = True

#---------------------------------------------------------------------------------------------------
def run_cmd(cmd,verbose=True,indent='  '):
   cmd_str = indent + hc.tclr.GREEN + cmd + hc.tclr.END
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
num_wave = len(wave_type_list)

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

plot = [None]*(num_var*num_wave)

res = hs.res_xy()
res.vpHeightF = 0.2
res.tmYLLabelFontHeightF         = 0.008
res.tmXBLabelFontHeightF         = 0.008
res.tiXAxisFontHeightF           = 0.01
res.tiYAxisFontHeightF           = 0.01

res.tmXTOn = False
res.tmYROn = False

lres = hs.res_xy()
lres.xyLineThicknessF = 1
lres.xyLineColor      = 'black'
lres.xyDashPattern    = 0

#---------------------------------------------------------------------------------------------------  
max_case_len = 0
for c in range(num_case):
   max_case_len = max(max_case_len,len(case[c]))
#---------------------------------------------------------------------------------------------------  
#---------------------------------------------------------------------------------------------------
def get_obs_name(case,var):
   obs_name, obs_var = None, None
   if case=='Obs' and var in ['OLR','FLNT','FLUT']: obs_var = 'olr'  ; obs_name = 'NOAA'
   if case=='Obs' and var in ['PRECT']            : obs_var = 'PRECT'; obs_name = 'IMERG'
   if case=='Obs' and var in ['U850']             : obs_var = 'U850' ; obs_name = 'ERA5'
   if obs_name is None: raise ValueError(f'get_obs_name(): obs_name cannot be None!  case: {case}  var: {var}')
   return obs_name, obs_var
#---------------------------------------------------------------------------------------------------

for w in range(num_wave):
   wave_type = wave_type_list[w]
   for v in range(num_var):
      print(' '*2+f'var: '+hc.tclr.GREEN+f'{var_list[v]}'+hc.tclr.END)

      data_list = []
      lat_list = []
      lon_list = []

      for c in range(num_case):
         #-------------------------------------------------------------------------
         tcase,tvar = case[c],None
         dst_root = f'{case_dir[c]}/{case[c]}/data_wave'
         if case[c]=='Obs':
            tcase,tvar = get_obs_name(case[c],eam_var_list[v])
            dst_root = None
            if tcase=='NOAA':  dst_root  = '/global/cfs/cdirs/m3312/whannah/obs_data/OLR'
            if tcase=='IMERG': dst_root  = '/global/cfs/cdirs/m3312/whannah/obs_data/IMERG'
            if tcase=='ERA5':  dst_root  = '/global/cfs/cdirs/m3312/whannah/obs_data/ERA5'
            if dst_root is None: raise ValueError('dst_root cannot be None!')
         # if case[c]=='Obs':
         #    if var_list[v]=='OLR': tcase, tvar, dst_root = 'NOAA', 'olr', '/global/cfs/cdirs/m3312/whannah/obs_data/OLR'
         #    name[c] = tcase
         #-------------------------------------------------------------------------
         if tvar is None:
            if comp_list[c]=='eam'  : tvar = eam_var_list[v]
            if comp_list[c]=='eamxx': tvar = eamxx_var_list[v]
            if tvar is None: raise ValueError(f'ERROR - variable name not found?!?!?  comp: {comp_list[c]} ')
         #-------------------------------------------------------------------------
         wave_file = f'{dst_root}/{tcase}.{tvar}.{wave_type}.nc'
         #-------------------------------------------------------------------------
         print(' '*4+f'case: '+hc.tclr.CYAN+f'{case[c]}'+hc.tclr.END)
         print(' '*4+f'wave_file: {hc.tclr.YELLOW}{wave_file}{hc.tclr.END}')
         #-------------------------------------------------------------------------
         # tcase,tvar = case[c],None
         # if case[c]=='Obs': tcase,tvar = get_obs_name(case[c],var_list[v])
         #-------------------------------------------------------------------------
         # print(' '*4+f'case: {hc.tclr.CYAN}{tcase:{max_case_len}}{hc.tclr.END}   wave_file: {hc.tclr.YELLOW}{wave_file}{hc.tclr.END}')
         # #-------------------------------------------------------------------------
         # if tvar is None:
         #    if comp_list[c]=='eam'  : tvar = eam_var_list[v]
         #    if comp_list[c]=='eamxx': tvar = eamxx_var_list[v]
         #    if tvar is None: raise ValueError(f'ERROR - variable name not found?!?!?  comp: {comp_list[c]} ')
         #-------------------------------------------------------------------------
         # if comp_list[c]=='eam'  : tvar = 'FLUT'
         # if comp_list[c]=='eamxx': tvar = 'LW_flux_up_at_model_top'
         # wave_file = f'{case_dir[c]}/{case[c]}/data_wave/{case[c]}.{tvar}.{wave_type}.nc'
         #-------------------------------------------------------------------------
         # load wave filtered data
         wave_ds = xr.open_dataset(wave_file)
         data = wave_ds[tvar]
         data = data.where( data['time.year']>=yr1, drop=True)
         data = data.where( data['time.year']<=yr2, drop=True)
         #-------------------------------------------------------------------------
         # hc.print_stat(data,name=tvar,compact=True,indent=' '*6)
         #-------------------------------------------------------------------------
         data_avg = data.mean(dim='lat').var(dim='time')
         #-------------------------------------------------------------------------
         data_list.append(np.ma.masked_invalid(data_avg.values))
         # lat_list.append(wave_ds['lat'].values)
         lon_list.append(np.ma.masked_invalid(wave_ds['lon'].values))
      #----------------------------------------------------------------------------
      # calculate limits for common color bar
      data_min = np.min([np.nanmin(d) for d in data_list])
      data_max = np.max([np.nanmax(d) for d in data_list])
      if data_min==data_max: raise ValueError(hc.tclr.RED+'WARNING: Difference is zero!'+hc.tclr.END)
      #----------------------------------------------------------------------------
      tres = copy.deepcopy(res)
      tres.trXMinF = 0
      tres.trXMaxF = 360
      tres.trYMinF = data_min
      tres.trYMaxF = data_max
      #----------------------------------------------------------------------------
      # # tres.cnLevels = np.arange(-30,30+3,3)/1e3
      # tres.cnLevelSelectionMode = "ExplicitLevels"
      # if not hasattr(tres,'cnLevels') : 
      #    cmin,cmax,cint,clev = ngl.nice_cntr_levels(data_min, data_max, cint=None, max_steps=21, \
      #                                               returnLevels=True, aboutZero=False )
      #    tres.cnLevels = np.linspace(cmin,cmax,num=21)
      #----------------------------------------------------------------------------
      # Create plot

      ip = v*num_wave+w if var_x_case else w*num_var+v

      tres.tiXAxisString = f'Longitude'
      tres.tiYAxisString = f'{var_list[v]} Variance [{unit_list[v]}]'

      for c in range(num_case):
         tres.xyLineColor = clr[c]
         tres.xyDashPattern = dsh[c]
         tplot = ngl.xy(wks,lon_list[c],data_list[c],tres)
         if c==0:
            plot[ip] = tplot
         else:
            ngl.overlay(plot[ip],tplot)

      ctr_str = 'MJO' if wave_type=='MJO' else f'{wave_type} Waves'
      hs.set_subtitles(wks, plot[ip], 
                       left_string='', 
                       center_string=ctr_str, 
                       right_string='',
                       font_height=0.01)

#---------------------------------------------------------------------------------------------------
# Add legend
for v in range(num_var):
   lgres = ngl.Resources()
   # lgres.vpWidthF           = 0.1
   # lgres.vpHeightF          = 0.15
   lgres.vpWidthF           = 0.05
   lgres.vpHeightF          = 0.08
   lgres.lgLabelFontHeightF = 0.006
   lgres.lgLineLabelsOn     = False
   lgres.lgLineThicknessF   = 10
   lgres.lgLabelJust        = 'CenterLeft'
   lgres.lgLineColors       = clr
   lgres.lgDashIndexes      = dsh
   # lgres.lgMonoDashIndex    = True

   lgd_lbl = name
   for c in range(num_case): lgd_lbl[c] = f'  {lgd_lbl[c]}'

   if eam_var_list[v]=='FLUT': lgd_lbl[0] = '  NOAA'
   if eam_var_list[v]=='U850': lgd_lbl[0] = '  ERA5'

   if eam_var_list[v]=='FLUT': pid = ngl.legend_ndc(wks, num_case, lgd_lbl, 0.35, 0.77, lgres)
   if eam_var_list[v]=='U850': pid = ngl.legend_ndc(wks, num_case, lgd_lbl, 0.85, 0.77, lgres)
#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

layout = [num_var,num_wave] if var_x_case else [num_wave,num_var]

pnl_res = hs.setres_panel()
# pnl_res.nglPanelYWhiteSpacePercent = 5
# pnl_res.nglPanelXWhiteSpacePercent = 5
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.01

ngl.panel(wks,plot,layout,pnl_res)
hc.trim_png(fig_file)

ngl.end()
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
