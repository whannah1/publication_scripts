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
# tmp_sub_ne1024  = 'data_remap_90x180/output.scream.decadal.3hourlyINST_ne30pg2'
# tmp_sub_ne1024  = 'data_remap_90x180/output.scream.decadal.6hourlyAVG_ne30pg2'
tmp_path_hst_v3 = '/global/cfs/cdirs/m3312/whannah/e3smv3_amip'
tmp_path_qbo_bm = '/global/cfs/cdirs/m4310/data/sims'

add_case('Obs')
add_case('v3.LR.amip_0101.QBObenchmark.20241008', n='EAMv3 AMIP',      comp='eam',  htype='eam.h2', p=tmp_path_qbo_bm, s='data_remap_90x180')
# add_case('v3.LR.amip_0201',                       n='v3.LR.amip_0201', comp='eam',  htype='eam.h2', p=tmp_path_hst_v3, s='data_remap_90x180')
add_case('decadal-production-run6',               n='SCREAM ne1024',   comp='eamxx',htype='output.scream.decadal.6hourlyINST_ne30pg2',p=tmp_path_ne1024,s='data_remap_90x180/output.scream.decadal.6hourlyINST_ne30pg2')


'''
ERA5 hourly data:
/global/cfs/projectdirs/m3522/cmip6/ERA5/e5.oper.an.pl/201111/e5.oper.an.pl.128_130_t.ll025sc.2011111900_2011111923.nc
/global/cfs/cdirs/m3312/whannah/ERA5/daily
'''
#-------------------------------------------------------------------------------

add_var('OLR',eam_var='FLUT', eamxx_var='LW_flux_up_at_model_top',  unit='W/m2')
# add_var('U',eam_var='U', eamxx_var='U',     unit='m/s')
# add_var('T',eam_var='T', eamxx_var='T_mid', unit='K')

num_plot_col = 1

#-------------------------------------------------------------------------------

# wave_type = 'Kelvin'
# wave_type = 'Rossby'
wave_type = 'MRG'

fig_type,fig_file = 'png',f'figs/FXX-wave-variance-map.{wave_type}'
# fig_file_diff = f'{fig_file}.diff'

# tmp_file_head = 'data/kelvin-variance-map'

lat1,lat2 = -15,15

# yr1,yr2 = 1975,2020
# yr1,yr2 = 1995,2004
yr1,yr2 = 1995,1995

recalculate = True

var_x_case  = False

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

plot = [None]*(num_var*num_case)

res = hs.res_contour_fill_map()
res.tmYLLabelFontHeightF         = 0.008
res.tmXBLabelFontHeightF         = 0.008
res.tmXBOn                       = False
res.tmYLOn                       = False
res.lbLabelFontHeightF           = 0.01
res.lbLabelBarOn                 = True
# res.mpGeophysicalLineColor       = 'white'
# res.mpProjection = 'Robinson'

res.mpMinLatF = lat1; res.mpMaxLatF = lat2
# res.mpMinLonF = lon1; res.mpMaxLonF = lon2

lres = hs.res_xy()
lres.xyLineThicknessF = 1
lres.xyLineColor      = 'black'
lres.xyDashPattern    = 0

#---------------------------------------------------------------------------------------------------  
#---------------------------------------------------------------------------------------------------
def get_obs_name(case,var):
   if case=='Obs' and var=='OLR': obs_var = 'olr'  ; obs_name = 'NOAA'
   return obs_name, obs_var
#---------------------------------------------------------------------------------------------------

for v in range(num_var):
   print(' '*4+f'var: '+hc.tcolor.GREEN+f'{var_list[v]}'+hc.tcolor.ENDC)

   data_list = []
   lat_list = []
   lon_list = []

   for c in range(num_case):
      #-------------------------------------------------------------------------
      tcase,tvar = case[c],None
      dst_root = f'{case_dir[c]}/{case[c]}/data_wave'
      # if case[c]=='Obs': tcase,tvar = get_obs_name(case[c],var_list[v])
      if case[c]=='Obs':
         if var_list[v]=='OLR': tcase, tvar, dst_root = 'NOAA', 'olr', '/global/cfs/cdirs/m3312/whannah/obs_data/OLR'
         name[c] = tcase
      #-------------------------------------------------------------------------
      if tvar is None:
         if comp_list[c]=='eam'  : tvar = eam_var_list[v]
         if comp_list[c]=='eamxx': tvar = eamxx_var_list[v]
         if tvar is None: raise ValueError(f'ERROR - variable name not found?!?!?  comp: {comp_list[c]} ')
      #-------------------------------------------------------------------------
      wave_file = f'{dst_root}/{tcase}.{tvar}.{wave_type}.nc'
      print(f'  wave_file: {wave_file}')
      #-------------------------------------------------------------------------
      # tcase,tvar = case[c],None
      # if case[c]=='Obs': tcase,tvar = get_obs_name(case[c],var_list[v])
      # print(' '*2+f'case: '+hc.tcolor.CYAN+f'{tcase}'+hc.tcolor.ENDC)
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
      # create kelvin wave indx
      wave_ds = xr.open_dataset(wave_file)
      data = wave_ds[tvar]
      data = data.where( data['time.year']>=yr1, drop=True)
      data = data.where( data['time.year']<=yr2, drop=True)
      #-------------------------------------------------------------------------
      # print()
      # print(data.var(dim='time'))
      # print()
      # hc.print_stat(data,name='KW Filtered OLR',compact=True,indent=' '*6)
      # exit()
      #-------------------------------------------------------------------------
      data_list.append(data.var(dim='time').values)
      lat_list.append(wave_ds['lat'].values)
      lon_list.append(wave_ds['lon'].values)
   #----------------------------------------------------------------------------
   # calculate limits for common color bar
   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])
   if data_min==data_max: raise ValueError(hc.tcolor.RED+'WARNING: Difference is zero!'+hc.tcolor.ENDC)
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)

   # tres.cnLevels = np.arange(-30,30+3,3)/1e3
   tres.cnLevelSelectionMode = "ExplicitLevels"
   if not hasattr(tres,'cnLevels') : 
      cmin,cmax,cint,clev = ngl.nice_cntr_levels(data_min, data_max, cint=None, max_steps=21, \
                                                 returnLevels=True, aboutZero=False )
      tres.cnLevels = np.linspace(cmin,cmax,num=21)
   #----------------------------------------------------------------------------
   # Create plot
   for c in range(num_case):

      tres.cnFillMode = "AreaFill"
      # tres.cnFillMode = "RasterFill"
      tres.sfXArray = lon_list[c]
      tres.sfYArray = lat_list[c]

      ip = v*num_case+c if var_x_case else c*num_var+v

      plot[ip] = ngl.contour_map(wks,np.ma.masked_invalid(data_list[c]),tres)
   
      hs.set_subtitles(wks, plot[ip], 
                       left_string=name[c], 
                       center_string='', 
                       right_string=f'{var_list[v]} Kelvin Wave Variance',
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
