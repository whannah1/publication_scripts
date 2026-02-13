import os, glob, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string, sys, cmocean
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
#-------------------------------------------------------------------------------
case,opts_list = [],[]
def add_case(case_in,**kwargs):
   case.append(case_in)
   case_opts = {}
   for k, val in kwargs.items(): case_opts[k] = val
   opts_list.append(case_opts)
#-------------------------------------------------------------------------------
tmp_path_ne1024 = '/global/cfs/cdirs/e3smdata/simulations/scream-decadal'
tmp_path_ne256  = '/global/cfs/cdirs/e3smdata/simulations/scream-decadal-ne256'
tmp_path_hst_v3 = '/global/cfs/cdirs/m3312/whannah/e3smv3_amip'
tmp_path_qbo_bm = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu/'
eamxx_var_pcp, eamxx_var_olr = 'precip_liq_surf_mass_flux','LW_flux_up_at_model_top'

add_case('Obs',name='NOAA/IMERG',var_pcp='PRECT',var_olr='olr',dsh=0,clr='black')
add_case('v3.LR.amip_0101.QBObenchmark.20241008',                                        name='EAMv3 AMIP',          dsh=0,clr='red',  var_pcp='PRECT',var_olr='FLUT',htype='eam.h1',                                  p=tmp_path_qbo_bm,s='data_remap_73x144')
add_case('decadal-production-run6',                                                      name='SCREAM 3-km',          dsh=0,clr='blue', var_pcp=eamxx_var_pcp,var_olr=eamxx_var_olr,htype='output.scream.decadal.6hourlyAVG_ne30pg2',p=tmp_path_ne1024,s='data_remap_73x144/output.scream.decadal.6hourlyAVG_ne30pg2')
add_case('ne256pg2_ne256pg2.F20TR-SCREAMv1.May-12.with.rain.frac.n0128',                 name='SCREAM 13-km control', dsh=0,clr='green',var_pcp=eamxx_var_pcp,var_olr=eamxx_var_olr,htype='6ha_ne30pg2.AVERAGE.nhours_x6',           p=tmp_path_ne256, s='data_remap_73x144/6ha_ne30pg2.AVERAGE.nhours_x6')
add_case('ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1', name='SCREAM 13-km tuned',   dsh=1,clr='green',var_pcp=eamxx_var_pcp,var_olr=eamxx_var_olr,htype='6ha_ne30pg2.AVERAGE.nhours_x6',           p=tmp_path_ne256, s='data_remap_73x144/6ha_ne30pg2.AVERAGE.nhours_x6')

#-------------------------------------------------------------------------------

num_plot_col = 1

#-------------------------------------------------------------------------------

fig_type,fig_file = 'png',f'figs/FXX-GEF-xy'
tmp_file_head = 'data/GEF-xy'

# yr1,yr2 = 2001,2002
yr1,yr2 = 2001,2005
# yr1,yr2 = 2001,2020

# lat1,lat2 = -15,15; lon1,lon2 = 60,180
lat1,lat2 = -30,30


recalculate = False

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
num_case = len(case)

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

plot = [None]*3

res = hs.res_xy()
res.vpHeightF = 0.2
res.tmYLLabelFontHeightF         = 0.009
res.tmXBLabelFontHeightF         = 0.009
res.tiXAxisFontHeightF           = 0.012
res.tiYAxisFontHeightF           = 0.012

res.xyLineThicknessF = 10

# res.tmXTOn = False; tmXUseBottom = False
# res.tmYROn = False; tmYUseLeft = False

res.xyXStyle = 'Log'

# res.tiXAxisString = 'Precipitation [mm/day]'

lres = hs.res_xy()
lres.xyLineThicknessF = 1
lres.xyLineColor      = 'black'
lres.xyDashPattern    = 0

#---------------------------------------------------------------------------------------------------
# define precip bins - see Kim et al. (2015) - https://doi.org/10.1175/JCLI-D-14-00767.1
num_pcp_bins = 51
pcp_bin_bot = np.zeros(num_pcp_bins)
pcp_bin_top = np.zeros(num_pcp_bins)
pcp_bin_bot[0] = 0
pcp_bin_top[0] = 0.09797

pcp_bin_bot[0] = np.power( 10., np.log10(pcp_bin_top[0]) - 0.065 ) # use this to avoid "kink" in lowest bin

for i in range(1,num_pcp_bins-1):
   pcp_bin_bot[i] = pcp_bin_top[i-1]
   pcp_bin_top[i] = np.power( 10., 0.065 + np.log10(pcp_bin_bot[i]) )
pcp_bin_top[num_pcp_bins-2] = 150
pcp_bin_bot[num_pcp_bins-1] = 150
pcp_bin_top[num_pcp_bins-1] = 1000

bin_ctr = np.zeros(num_pcp_bins)
for i in range(num_pcp_bins):
   bin_ctr[i] = ( pcp_bin_top[i] + pcp_bin_bot[i] )/2.


# for i in range(num_pcp_bins): print(f'{i:3}  {pcp_bin_bot[i]:12.4f}  /  {bin_ctr[i]:12.4f}  /  {pcp_bin_top[i]:12.4f}')
# for i in range(5): print(f'{i:3}  {pcp_bin_bot[i]:12.4f}  /  {bin_ctr[i]:12.4f}  /  {pcp_bin_top[i]:12.4f}')
# exit()

#---------------------------------------------------------------------------------------------------
bin_list = []
gef_list = []
olr_list = []
pct_list = []

for c in range(num_case):
   #----------------------------------------------------------------------------
   tcase,tvar = case[c],None
   case_opts = opts_list[c]
   #----------------------------------------------------------------------------
   # if case[c]=='Obs':
   #    tcase,tvar = get_obs_name(case[c],eam_var_list[v])
   #    dst_root = None
   #    if tcase=='NOAA':  dst_root  = '/global/cfs/cdirs/m3312/whannah/obs_data/OLR'
   #    if tcase=='IMERG': dst_root  = '/global/cfs/cdirs/m3312/whannah/obs_data/IMERG'
   #    if tcase=='ERA5':  dst_root  = '/global/cfs/cdirs/m3312/whannah/obs_data/ERA5'
   #    if dst_root is None: raise ValueError('dst_root cannot be None!')
   #----------------------------------------------------------------------------
   # if tvar is None:
   #    if comp_list[c]=='eam'  : tvar = eam_var_list[v]
   #    if comp_list[c]=='eamxx': tvar = eamxx_var_list[v]
   #    if tvar is None: raise ValueError(f'ERROR - variable name not found?!?!?  comp: {comp_list[c]} ')
   #----------------------------------------------------------------------------
   print(' '*2+f'case: '+hc.tclr.CYAN+f'{case[c]}'+hc.tclr.END)
   #----------------------------------------------------------------------------
   # if tvar is None:
   #    if comp_list[c]=='eam'  : tvar = eam_var_list[v]
   #    if comp_list[c]=='eamxx': tvar = eamxx_var_list[v]
   #    if tvar is None: raise ValueError(f'ERROR - variable name not found?!?!?  comp: {comp_list[c]} ')
   #----------------------------------------------------------------------------
   tmp_file = f'{tmp_file_head}.{tcase}'
   tmp_file+= f'.yr_{yr1}_{yr2}'
   if 'lat1' in locals(): tmp_file+= f'.lat_{lat1}_{lat2}'
   if 'lon1' in locals(): tmp_file+= f'.lat_{lon1}_{lon2}'
   tmp_file+= f'.nc'
   #----------------------------------------------------------------------------
   if recalculate:
      #-------------------------------------------------------------------------
      if case[c]=='Obs':
         file_list_olr = ['/global/cfs/cdirs/m3312/whannah/obs_data/OLR/olr.day.mean.nc']
         file_list_pcp = ['/global/cfs/cdirs/m3312/whannah/obs_data/IMERG/IMERG_Daily_PRECT_200101_202012.remap_73x144.nc']
      else:
         case_dir = case_opts['p']
         case_sub = case_opts['s']
         htype    = case_opts['htype']
         file_path_pcp = f'{case_dir}/{case[c]}/{case_sub}/*{htype}*'; file_list_pcp = sorted(glob.glob(file_path_pcp))
         file_path_olr = f'{case_dir}/{case[c]}/{case_sub}/*{htype}*'; file_list_olr = sorted(glob.glob(file_path_olr))
      #-------------------------------------------------------------------------
      if file_list_pcp==[]: print();print(f'file_path_pcp: {file_path_pcp}');exit('ERROR: file_path_pcp is empty!')
      if file_list_olr==[]: print();print(f'file_path_olr: {file_path_olr}');exit('ERROR: file_path_olr is empty!')
      #-------------------------------------------------------------------------
      ds_pcp = xr.open_mfdataset(file_list_pcp); data_pcp = ds_pcp[case_opts['var_pcp']]
      ds_olr = xr.open_mfdataset(file_list_olr); data_olr = ds_olr[case_opts['var_olr']]
      #-------------------------------------------------------------------------
      # print(); print(data_pcp['time.year'])
      # print(); print(data_olr['time.year'])
      # exit()
      #-------------------------------------------------------------------------
      data_pcp = data_pcp.where( data_pcp['time.year']>=yr1, drop=True).where( data_pcp['time.year']<=yr2, drop=True)
      data_olr = data_olr.where( data_olr['time.year']>=yr1, drop=True).where( data_olr['time.year']<=yr2, drop=True)
      #-------------------------------------------------------------------------
      if len(data_pcp.time)==0: raise ValueError('ERROR - something went wrong - data_pcp time dimension is zero!')
      if len(data_olr.time)==0: raise ValueError('ERROR - something went wrong - data_olr time dimension is zero!')
      #-------------------------------------------------------------------------
      # Convert to daily mean
      if case[c]!='Obs':
         data_pcp = data_pcp.resample(time='D').mean(dim='time')
         data_olr = data_olr.resample(time='D').mean(dim='time')
      #-------------------------------------------------------------------------
      # flip latitude dimension and make time coords consistent
      if case[c]=='Obs':
         # print(); print(data_pcp['time'].coords)
         # print(); print(data_olr['time'].coords)
         data_olr = data_olr.reindex(lat=data_olr.lat[::-1])
         data_pcp['time'] = data_olr['time']
         # print(); print(data_pcp['time'].coords)
         # print(); print(data_olr['time'].coords)
      # exit()
      #-------------------------------------------------------------------------
      def remove_leap(da):
         # Create a boolean mask for all dates that are NOT February 29th
         is_not_leap_day = ~((da.time.dt.month == 2) & (da.time.dt.day == 29))
         da = da.sel(time=is_not_leap_day)
         return da
      #-------------------------------------------------------------------------
      # remove leap days for obs data
      if case[c]=='Obs':
         data_pcp = remove_leap(data_pcp)
         data_olr = remove_leap(data_olr)
      #-------------------------------------------------------------------------
      data_pcp = data_pcp*86400.*1e3 # units = mm/day
      #-------------------------------------------------------------------------
      # reduce to equatorial region
      mask = xr.DataArray( np.ones([len(data_pcp['lat']),len(data_pcp['lon'])],dtype=bool), dims=('lat','lon') )
      if 'lat1' in locals(): mask = mask & (data_pcp['lat']>=lat1) & (data_pcp['lat']<=lat2)
      if 'lon1' in locals(): mask = mask & (data_pcp['lon']>=lon1) & (data_pcp['lon']<=lon2)
      data_pcp = data_pcp.where( mask,drop=True)
      data_olr = data_olr.where( mask,drop=True)
      #-------------------------------------------------------------------------
      data_pcp.load()
      data_olr.load()
      #-------------------------------------------------------------------------
      # hc.print_stat(data_pcp,name=case_opts['var_pcp'],compact=True,indent=' '*6)
      # hc.print_stat(data_olr,name=case_opts['var_olr'],compact=True,indent=' '*6)
      #-------------------------------------------------------------------------
      # remove daily climatology
      print(' '*4+'Removing daily climatology...')
      data_pcp_b = data_pcp.groupby('time.dayofyear').mean(dim='time',keep_attrs=True,skipna=True)
      data_olr_b = data_olr.groupby('time.dayofyear').mean(dim='time',keep_attrs=True,skipna=True)
      data_pcp_p = data_pcp.copy(deep=True)
      data_olr_p = data_olr.copy(deep=True)
      for d in range(365):
         data_pcp_p[d::365,:,:] = data_pcp[d::365,:,:] - data_pcp_b[d,:,:]
         data_olr_p[d::365,:,:] = data_olr[d::365,:,:] - data_olr_b[d,:,:]
      #-------------------------------------------------------------------------
      # hc.print_stat(data_pcp_p,name=case_opts['var_pcp']+'_p',compact=True,indent=' '*6)
      # hc.print_stat(data_olr_p,name=case_opts['var_olr']+'_p',compact=True,indent=' '*6)
      #-------------------------------------------------------------------------
      # bin the data
      print(' '*4+'Binning the data...')
      bin_coord = xr.DataArray( bin_ctr )
      dims,coord = 'bins',[('bins', bin_ctr)]
      bin_cnt = xr.DataArray( np.full(num_pcp_bins,np.nan), coords=coord, dims=dims )
      bin_olr = xr.DataArray( np.full(num_pcp_bins,np.nan), coords=coord, dims=dims )
      bin_pcp = xr.DataArray( np.full(num_pcp_bins,np.nan), coords=coord, dims=dims )
      # bin_gef = xr.DataArray( np.zeros(shape,dtype=data.dtype), coords=coord, dims=dims )
      for b in range(num_pcp_bins):
         condition = ( data_pcp_p.values>=pcp_bin_bot[b] ) \
                    &( data_pcp_p.values <pcp_bin_top[b] )
         condition = xr.DataArray( condition, coords=data_pcp.coords )
         # print()
         # print(); print(condition.coords)
         # print(); print(data_olr_p.coords)
         # print()
         # exit()
         bin_cnt[b] = condition.sum()
         if bin_cnt[b]>0:
            bin_olr[b] = data_olr_p.where(condition,drop=True).mean(skipna=True)#.sum() / (bin_spc_pct/1e2)
            bin_pcp[b] = data_pcp_p.where(condition,drop=True).mean(skipna=True)#.sum() / (bin_spc_pct/1e2)
            # print()
            # print();print(f'bin_cnt[b]: {bin_cnt[b].values}')
            # print();print(f'bin_pcp[b]: {bin_pcp[b].values}')
            # print();print(f'bin_olr[b]: {bin_olr[b].values}')
            # olr_sum = data_olr_p.where(condition,drop=True).sum().values
            # nan_cnt = np.sum( np.isnan(   data_olr_p.where(condition,drop=True).values) )
            # fin_cnt = np.sum( np.isfinite(data_olr_p.where(condition,drop=True).values) )
            # print();print(f'olr_sum: {olr_sum}')
            # print();print(f'nan_cnt: {nan_cnt}')
            # print();print(f'fin_cnt: {fin_cnt}')
            # print()
            # print(); print( data_olr_p.where(condition,drop=True) )
            # print(); print( data_olr_p.where(condition,drop=True).values )
            # print()
            # print(); print(condition)
            # print(); print(np.sum(condition.values))
            # exit()
      
      bin_pct = bin_cnt.values/bin_cnt.sum().values*1e2
   #----------------------------------------------------------------------------
   # Write to file 
      print(' '*4+f'writing to file - {tmp_file}')
      ds = xr.Dataset()
      ds['bin_ctr'] = bin_ctr
      ds['bin_pcp'] = bin_pcp
      ds['bin_olr'] = bin_olr
      ds['bin_pct'] = bin_pct
      ds.to_netcdf(path=tmp_file,mode='w')
   else:
      print(' '*4+f'loading pre-calculated spectra... {tmp_file}')
      ds = xr.open_mfdataset( tmp_file )
      bin_ctr = ds['bin_ctr']
      bin_pcp = ds['bin_pcp']
      bin_olr = ds['bin_olr']
      bin_pct = ds['bin_pct']
   #----------------------------------------------------------------------------
   rho = 1e3      # kg/m3
   Lv  = 2.5104e6 # J/kg
   bin_pcp = bin_pcp / (86400.*1e3) * rho * Lv # need to convert precip from mm/day to W/m2
   # bin_pcp = bin_pcp * rho * Lv # need to convert precip from m/s to W/m2
   bin_gef = -1*bin_olr / bin_pcp
   #----------------------------------------------------------------------------
   # hc.print_stat(bin_ctr,name='bin_ctr',compact=True,indent=' '*6)
   # hc.print_stat(bin_gef,name='bin_gef',compact=True,indent=' '*6)
   # hc.print_stat(bin_pcp,name='bin_pcp',compact=True,indent=' '*6)
   # hc.print_stat(bin_olr,name='bin_olr',compact=True,indent=' '*6)
   # hc.print_stat(bin_pct,name='bin_pct',compact=True,indent=' '*6)
   #----------------------------------------------------------------------------
   bin_list.append( np.ma.masked_invalid(bin_ctr) )
   gef_list.append( np.ma.masked_invalid(bin_gef) )
   olr_list.append( np.ma.masked_invalid(bin_olr) )
   pct_list.append( np.ma.masked_invalid(bin_pct) )

#-------------------------------------------------------------------------------
# calculate limits for common color bar
# data_min = np.min([np.nanmin(d) for d in data_list])
# data_max = np.max([np.nanmax(d) for d in data_list])
# if data_min==data_max: raise ValueError(hc.tclr.RED+'WARNING: Difference is zero!'+hc.tclr.END)
#-------------------------------------------------------------------------------
tres = copy.deepcopy(res)
tres.trXMinF = 0.1 # np.min(bin_ctr)
# tres.trXMaxF = 100 # 
tres.trXMaxF = 500

clr,dsh = [None]*num_case,[None]*num_case
for c in range(num_case):
   clr[c] = opts_list[c]['clr']
   dsh[c] = opts_list[c]['dsh']

tres.xyLineColors   = clr
tres.xyDashPatterns = dsh


gres = copy.deepcopy(tres); gres.trYMinF,gres.trYMaxF = -1,4
ores = copy.deepcopy(tres); ores.trYMinF,ores.trYMaxF = -150,0
pres = copy.deepcopy(tres); pres.trYMinF,pres.trYMaxF = 0,5


gres = copy.deepcopy(tres); gres.trYMinF,gres.trYMaxF = 0.01,4; gres.xyYStyle='Log'

#-------------------------------------------------------------------------------
# Create plot

# ip = v*num_wave+w if var_x_case else w*num_var+v

# tres.tiXAxisString = f'Longitude'
# tres.tiYAxisString = f'{var_list[v]} Variance [{unit_list[v]}]'


gres.tiYAxisString = 'GEF'
ores.tiYAxisString = 'OLR [W/m2]'
pres.tiYAxisString = 'Frequency [%]'
pres.tiXAxisString = 'Precipitation [mm/day]'

# adjust panel widths so they are consistent
ores.vpWidthF = 0.6*1.15

plot[0] = ngl.xy(wks, np.stack(bin_list) , np.stack(gef_list) ,gres)
plot[1] = ngl.xy(wks, np.stack(bin_list) , np.stack(olr_list) ,ores)
plot[2] = ngl.xy(wks, np.stack(bin_list) , np.stack(pct_list) ,pres)

# for c in range(num_case):
#    tres.xyLineColor = clr[c]
#    tres.xyDashPattern = dsh[c]
#    tplot = ngl.xy(wks,lon_list[c],data_list[c],tres)
#    if c==0:
#       plot[ip] = tplot
#    else:
#       ngl.overlay(plot[ip],tplot)

# ctr_str = 'MJO' if wave_type=='MJO' else f'{wave_type} Waves'

# hs.set_subtitles(wks, plot[0], left_string='', center_string='', right_string='',font_height=0.01)
# hs.set_subtitles(wks, plot[1], left_string='', center_string='', right_string='',font_height=0.01)
# hs.set_subtitles(wks, plot[2], left_string='', center_string='', right_string='',font_height=0.01)

#---------------------------------------------------------------------------------------------------
# Add legend

lgres = ngl.Resources()
lgres.vpWidthF,lgres.vpHeightF           = 0.08,0.12
lgres.lgLabelFontHeightF = 0.01
# lgres.vpWidthF,lgres.vpHeightF           = 0.18,0.15
# lgres.lgLabelFontHeightF = 0.014
lgres.lgLineLabelsOn     = False
lgres.lgLineThicknessF   = 12
lgres.lgLabelJust        = 'CenterLeft'
lgres.lgLineColors       = clr
lgres.lgDashIndexes      = dsh
# lgres.lgMonoDashIndex    = True

lgd_lbl = [None]*num_case
for c in range(num_case): lgd_lbl[c] = f'  {opts_list[c]["name"]}'

# pid = ngl.legend_ndc(wks, num_case, lgd_lbl, 0.6, 0.97, lgres)
pid = ngl.legend_ndc(wks, num_case, lgd_lbl, 0.62, 0.89, lgres)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

# layout = [num_var,num_wave] if var_x_case else [num_wave,num_var]
num_plot_col = 1
layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

pnl_res = hs.setres_panel()
# pnl_res.nglPanelYWhiteSpacePercent = 2
# pnl_res.nglPanelXWhiteSpacePercent = 5
pnl_res.nglPanelTop = 0.9
pnl_res.nglPanelBottom = 0.1
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "BottomLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.012

ngl.panel(wks,plot,layout,pnl_res)
hc.trim_png(fig_file)

ngl.end()
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
