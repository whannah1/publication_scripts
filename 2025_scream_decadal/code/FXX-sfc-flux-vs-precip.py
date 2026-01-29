import os, glob, subprocess as sp, numpy as np, xarray as xr, copy, string, sys, cmocean
import hapy_common as hc#, hapy_E3SM   as he#, hapy_setres as hs
import matplotlib.pyplot as plt
#-------------------------------------------------------------------------------
case,opts_list = [],[]
def add_case(case_in,**kwargs):
   case.append(case_in)
   case_opts = {}
   for k, val in kwargs.items(): case_opts[k] = val
   opts_list.append(case_opts)
#-------------------------------------------------------------------------------
hx_data_root = '/global/cfs/cdirs/e3smdata/simulations/scream-decadal/DYNAMO_MJO_hindcasts'
# hx_data_sub  = 'data_remap_73x144'
hx_data_sub  = 'run'

htype = 'output.scream.2D.1hr.ne30pg2.AVERAGE.nhours_x1'
# htype = 'output.scream.2D.6hr.INSTANT.nhours_x6' # no flux data in this stream
htype_3D = 'output.scream.3D.6hr.ne30pg2.AVERAGE.nhours_x6'

eamxx_var_pcp = 'precip_total_surf_mass_flux'
eamxx_var_shf = 'surf_sens_flux'
eamxx_var_lhf = 'surface_upward_latent_heat_flux'



add_case('SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-15.rfrac_fix_0', name='SCREAMv1 3-km control',           dsh=0,clr='blue', var_pcp=eamxx_var_pcp,var_shf=eamxx_var_shf,var_lhf=eamxx_var_lhf,htype=htype,p=hx_data_root,s=hx_data_sub)
add_case('SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-15.rfrac_fix_1', name='SCREAMv1 3-km w/ cold pool fix',  dsh=1,clr='blue', var_pcp=eamxx_var_pcp,var_shf=eamxx_var_shf,var_lhf=eamxx_var_lhf,htype=htype,p=hx_data_root,s=hx_data_sub)

scrip_file = 'files_grid/73x144_scrip.nc'

#-------------------------------------------------------------------------------

num_plot_col = 1

#-------------------------------------------------------------------------------

fig_type,fig_file = 'png',f'figs/FXX-sfc-flux-vs-precip'
tmp_file_head = 'data/sfc-flux-vs-precip'

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

# wkres = ngl.Resources()
# npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
# wks = ngl.open_wks(fig_type,fig_file,wkres)

# plot = [None]*3

# res = hs.res_xy()
# res.vpwidthF = 0.5
# res.vpHeightF = 0.2
# res.tmYLLabelFontHeightF         = 0.009
# res.tmXBLabelFontHeightF         = 0.009
# res.tiXAxisFontHeightF           = 0.012
# res.tiYAxisFontHeightF           = 0.012

# res.xyLineThicknessF = 10

# res.trXMinF = 0.1
# res.trXMaxF = 500

# # res.tmXTOn = False; tmXUseBottom = False
# # res.tmYROn = False; tmYUseLeft = False

# res.xyXStyle = 'Log'

# # res.tiXAxisString = 'Precipitation [mm/day]'

# lres = hs.res_xy()
# lres.xyLineThicknessF = 1
# lres.xyLineColor      = 'black'
# lres.xyDashPattern    = 0

#---------------------------------------------------------------------------------------------------
# define precip bins - see Kim et al. (2015) - https://doi.org/10.1175/JCLI-D-14-00767.1
num_pcp_bins = 57#51
pcp_bin_bot = np.zeros(num_pcp_bins)
pcp_bin_top = np.zeros(num_pcp_bins)
pcp_bin_bot[0] = 0
pcp_bin_top[0] = 0.09797

pcp_bin_bot[0] = np.power( 10., np.log10(pcp_bin_top[0]) - 0.065 ) # use this to avoid "kink" in lowest bin

# for i in range(1,num_pcp_bins-1):
for i in range(1,num_pcp_bins):
   pcp_bin_bot[i] = pcp_bin_top[i-1]
   pcp_bin_top[i] = np.power( 10., 0.065 + np.log10(pcp_bin_bot[i]) )
# pcp_bin_top[num_pcp_bins-2] = 150
# pcp_bin_bot[num_pcp_bins-1] = 150
# pcp_bin_top[num_pcp_bins-1] = 1000

bin_ctr = np.zeros(num_pcp_bins)
for i in range(num_pcp_bins):
   bin_ctr[i] = ( pcp_bin_top[i] + pcp_bin_bot[i] )/2.


# for i in range(num_pcp_bins): print(f'{i:3}  {pcp_bin_bot[i]:12.4f}  /  {bin_ctr[i]:12.4f}  /  {pcp_bin_top[i]:12.4f}')
# exit()

#---------------------------------------------------------------------------------------------------
bin_list = []
shf_list = []
lhf_list = []
pct_list = []
wnd_list = []
tmp_list = []

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
   # tmp_file+= f'.yr_{yr1}_{yr2}'
   if 'lat1' in locals(): tmp_file+= f'.lat_{lat1}_{lat2}'
   if 'lon1' in locals(): tmp_file+= f'.lat_{lon1}_{lon2}'
   tmp_file+= f'.nc'
   #----------------------------------------------------------------------------
   if recalculate:
      #-------------------------------------------------------------------------
      if case[c]=='Obs':
         file_list_pcp = ['/global/cfs/cdirs/m3312/whannah/obs_data/IMERG/IMERG_Daily_PRECT_200101_202012.remap_73x144.nc']
         # file_list_flx = ['/global/cfs/cdirs/m3312/whannah/obs_data/OLR/olr.day.mean.nc']
      else:
         case_dir = case_opts['p']
         case_sub = case_opts['s']
         htype    = case_opts['htype']
         file_path_pcp = f'{case_dir}/{case[c]}/{case_sub}/*{htype}*'; file_list_pcp = sorted(glob.glob(file_path_pcp))
         file_path_flx = f'{case_dir}/{case[c]}/{case_sub}/*{htype}*'; file_list_flx = sorted(glob.glob(file_path_flx))
         file_path_wnd = f'{case_dir}/{case[c]}/{case_sub}/*{htype_3D}*'; file_list_wnd = sorted(glob.glob(file_path_wnd))
      #-------------------------------------------------------------------------
      if file_list_pcp==[]: print();print(f'file_path_pcp: {file_path_pcp}');exit('ERROR: file_path_pcp is empty!')
      if file_list_flx==[]: print();print(f'file_path_flx: {file_path_flx}');exit('ERROR: file_path_flx is empty!')
      if file_list_wnd==[]: print();print(f'file_path_flx: {file_path_wnd}');exit('ERROR: file_path_flx is empty!')
      #-------------------------------------------------------------------------
      ds_pcp = xr.open_mfdataset(file_list_pcp,data_vars='all')
      ds_flx = xr.open_mfdataset(file_list_flx,data_vars='all')
      data_pcp = ds_pcp[case_opts['var_pcp']]
      data_shf = ds_flx[case_opts['var_shf']]
      data_lhf = ds_flx[case_opts['var_lhf']]
      #-------------------------------------------------------------------------
      ds_wnd = xr.open_mfdataset(file_list_wnd,data_vars='all')
      data_wnd = ds_wnd['horiz_winds'].isel(lev=-1)
      data_wnd = np.sqrt( np.square(data_wnd.isel(dim2=0)) + np.square(data_wnd.isel(dim2=1)) )
      data_tmp = ds_wnd['T_mid'].isel(lev=-1)
      #-------------------------------------------------------------------------
      if len(data_pcp.time)==0: raise ValueError('ERROR - something went wrong - data_pcp time dimension is zero!')
      if len(data_shf.time)==0: raise ValueError('ERROR - something went wrong - data_shf time dimension is zero!')
      if len(data_lhf.time)==0: raise ValueError('ERROR - something went wrong - data_lhf time dimension is zero!')
      #-------------------------------------------------------------------------
      # Convert to 6-hr mean for comparison with 3D data - throw out first time for consistency
      if case[c]!='Obs':
         data_pcp = data_pcp.resample(time='6h').mean(dim='time')[1:,:]
         data_shf = data_shf.resample(time='6h').mean(dim='time')[1:,:]
         data_lhf = data_lhf.resample(time='6h').mean(dim='time')[1:,:]
      # print()
      # print(data_pcp.shape)
      # print(data_shf.shape)
      # print(data_lhf.shape)
      # print(data_wnd.shape)
      # print()
      # for t in data_lhf.time[0:4]: print(t.values)
      # print()
      # for t in data_wnd.time[0:4]: print(t.values)
      # print()
      # exit()
      #-------------------------------------------------------------------------
      # # Convert to daily mean
      # if case[c]!='Obs':
      #    data_pcp = data_pcp.resample(time='D').mean(dim='time')
      #    data_flx = data_flx.resample(time='D').mean(dim='time')
      #-------------------------------------------------------------------------
      # # flip latitude dimension and make time coords consistent
      # if case[c]=='Obs':
      #    data_flx = data_flx.reindex(lat=data_flx.lat[::-1])
      #    data_pcp['time'] = data_flx['time']
      #-------------------------------------------------------------------------
      # def remove_leap(da):
      #    # Create a boolean mask for all dates that are NOT February 29th
      #    is_not_leap_day = ~((da.time.dt.month == 2) & (da.time.dt.day == 29))
      #    da = da.sel(time=is_not_leap_day)
      #    return da
      #-------------------------------------------------------------------------
      # # remove leap days for obs data
      # if case[c]=='Obs':
      #    data_pcp = remove_leap(data_pcp)
      #    data_flx = remove_leap(data_flx)
      #-------------------------------------------------------------------------
      data_pcp = data_pcp*86400.*1e3 # units = mm/day
      #-------------------------------------------------------------------------
      # reduce to equatorial region
      # mask = xr.DataArray( np.ones([len(ds_pcp['lat']),len(ds_pcp['lon'])],dtype=bool), dims=('lat','lon') )
      mask = xr.DataArray( np.ones(data_pcp.shape,dtype=bool), dims=data_pcp.dims )
      if 'lat1' in locals(): mask = mask & (ds_pcp['lat'].values>=lat1) & (ds_pcp['lat'].values<=lat2)
      if 'lon1' in locals(): mask = mask & (ds_pcp['lon'].values>=lon1) & (ds_pcp['lon'].values<=lon2)
      data_pcp = data_pcp.where( mask,drop=True)
      data_shf = data_shf.where( mask,drop=True)
      data_lhf = data_lhf.where( mask,drop=True)
      data_wnd = data_wnd.where( mask,drop=True)
      data_tmp = data_tmp.where( mask,drop=True)
      #-------------------------------------------------------------------------
      data_pcp.load()
      data_shf.load()
      data_lhf.load()
      data_wnd.load()
      data_tmp.load()
      #-------------------------------------------------------------------------
      # hc.print_stat(data_pcp,name=case_opts['var_pcp'],compact=True,indent=' '*6)
      # hc.print_stat(data_flx,name=case_opts['var_flx'],compact=True,indent=' '*6)
      #-------------------------------------------------------------------------
      # # remove daily climatology
      # print(' '*4+'Removing daily climatology...')
      # data_pcp_b = data_pcp.groupby('time.dayofyear').mean(dim='time',keep_attrs=True,skipna=True)
      # data_flx_b = data_flx.groupby('time.dayofyear').mean(dim='time',keep_attrs=True,skipna=True)
      # data_pcp_p = data_pcp.copy(deep=True)
      # data_flx_p = data_flx.copy(deep=True)
      # for d in range(365):
      #    data_pcp_p[d::365,:,:] = data_pcp[d::365,:,:] - data_pcp_b[d,:,:]
      #    data_flx_p[d::365,:,:] = data_flx[d::365,:,:] - data_flx_b[d,:,:]
      #-------------------------------------------------------------------------
      # hc.print_stat(data_pcp_p,name=case_opts['var_pcp']+'_p',compact=True,indent=' '*6)
      # hc.print_stat(data_flx_p,name=case_opts['var_flx']+'_p',compact=True,indent=' '*6)
      #-------------------------------------------------------------------------
      # bin the data
      print(' '*4+'Binning the data...')
      bin_coord = xr.DataArray( bin_ctr )
      dims,coord = 'bins',[('bins', bin_ctr)]
      bin_cnt = xr.DataArray( np.full(num_pcp_bins,np.nan), coords=coord, dims=dims )
      bin_shf = xr.DataArray( np.full(num_pcp_bins,np.nan), coords=coord, dims=dims )
      bin_lhf = xr.DataArray( np.full(num_pcp_bins,np.nan), coords=coord, dims=dims )
      bin_pcp = xr.DataArray( np.full(num_pcp_bins,np.nan), coords=coord, dims=dims )
      bin_wnd = xr.DataArray( np.full(num_pcp_bins,np.nan), coords=coord, dims=dims )
      bin_tmp = xr.DataArray( np.full(num_pcp_bins,np.nan), coords=coord, dims=dims )
      for b in range(num_pcp_bins):
         condition = ( data_pcp.values>=pcp_bin_bot[b] ) \
                    &( data_pcp.values <pcp_bin_top[b] )
         condition = xr.DataArray( condition, coords=data_pcp.coords )
         bin_cnt[b] = condition.sum()
         if bin_cnt[b]>0:
            bin_shf[b] = data_shf.where(condition,drop=True).mean(skipna=True)#.sum() / (bin_spc_pct/1e2)
            bin_lhf[b] = data_lhf.where(condition,drop=True).mean(skipna=True)#.sum() / (bin_spc_pct/1e2)
            bin_pcp[b] = data_pcp.where(condition,drop=True).mean(skipna=True)#.sum() / (bin_spc_pct/1e2)
            bin_wnd[b] = data_wnd.where(condition,drop=True).mean(skipna=True)#.sum() / (bin_spc_pct/1e2)
            bin_tmp[b] = data_tmp.where(condition,drop=True).mean(skipna=True)#.sum() / (bin_spc_pct/1e2)
      
      bin_pct = bin_cnt.values/bin_cnt.sum().values*1e2
   #----------------------------------------------------------------------------
   # Write to file 
      print(' '*4+f'writing to file - {tmp_file}')
      ds = xr.Dataset()
      ds['bin_ctr'] = bin_ctr
      ds['bin_pcp'] = bin_pcp
      ds['bin_shf'] = bin_shf
      ds['bin_lhf'] = bin_lhf
      ds['bin_pct'] = bin_pct
      ds['bin_wnd'] = bin_wnd
      ds['bin_tmp'] = bin_tmp
      ds.to_netcdf(path=tmp_file,mode='w')
      print(' '*4+f'done')
   else:
      print(' '*4+f'loading pre-calculated spectra... {tmp_file}')
      ds = xr.open_mfdataset( tmp_file )
      bin_ctr = ds['bin_ctr']
      bin_pcp = ds['bin_pcp']
      bin_shf = ds['bin_shf']
      bin_lhf = ds['bin_lhf']
      bin_pct = ds['bin_pct']
      bin_wnd = ds['bin_wnd']
      bin_tmp = ds['bin_tmp']
   #----------------------------------------------------------------------------
   # hc.print_stat(bin_ctr,name='bin_ctr',compact=True,indent=' '*6)
   hc.print_stat(bin_pcp,name='bin_pcp',compact=True,indent=' '*6)
   hc.print_stat(bin_shf,name='bin_shf',compact=True,indent=' '*6)
   hc.print_stat(bin_lhf,name='bin_lhf',compact=True,indent=' '*6)
   hc.print_stat(bin_wnd,name='bin_wnd',compact=True,indent=' '*6)
   hc.print_stat(bin_tmp,name='bin_tmp',compact=True,indent=' '*6)
   # hc.print_stat(bin_pct,name='bin_pct',compact=True,indent=' '*6)
   #----------------------------------------------------------------------------
   bin_list.append( np.ma.masked_invalid(bin_ctr) )
   shf_list.append( np.ma.masked_invalid(bin_shf) )
   lhf_list.append( np.ma.masked_invalid(bin_lhf) )
   pct_list.append( np.ma.masked_invalid(bin_pct) )
   wnd_list.append( np.ma.masked_invalid(bin_wnd) )
   tmp_list.append( np.ma.masked_invalid(bin_tmp) )

#-------------------------------------------------------------------------------
# calculate limits for common color bar
# data_min = np.min([np.nanmin(d) for d in data_list])
# data_max = np.max([np.nanmax(d) for d in data_list])
# if data_min==data_max: raise ValueError(hc.tclr.RED+'WARNING: Difference is zero!'+hc.tclr.END)
#-------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
# use matplotlib

fig, (ax1,ax3,ax_pct) = plt.subplots(3,1,figsize=(10,12))

kwargs = []
kwargs.append({'linestyle':'solid',})
kwargs.append({'linestyle':'dashed',})

# plt.rc('axes',  titlesize=14)
# plt.rc('axes',  labelsize=14)
# plt.rc('xtick', labelsize=14)
# plt.rc('ytick', labelsize=14)
title_fontsize = 14
label_fontsize = 12


color = 'red'#'tab:red'
ax1.set_xlabel('Precipitation [mm/day]', fontsize=title_fontsize)
ax1.set_ylabel('Sensible Heat Flux [W/m2]', color=color, fontsize=title_fontsize)
for c in range(num_case):
   ax1.plot(bin_list[c], np.ma.masked_invalid(shf_list[c]), color=color, **kwargs[c])
ax1.tick_params(axis='x', labelsize=label_fontsize)
ax1.tick_params(axis='y', labelsize=label_fontsize, labelcolor=color)
ax1.set_xlim(xmin=0.1)
ax1.set_xscale('log')


ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

color = 'blue'#'tab:blue'
ax2.set_ylabel('Latent Heat Flux [W/m2]', color=color, fontsize=title_fontsize)  # we already handled the x-label with ax1
for c in range(num_case):
   ax2.plot(bin_list[c], np.ma.masked_invalid(lhf_list[c]), color=color, **kwargs[c])
ax2.tick_params(axis='y', labelsize=label_fontsize, labelcolor=color)
ax2.set_ylim(ymin=-50)



color = 'green'
ax3.set_xlabel('Precipitation [mm/day]', fontsize=title_fontsize)
ax3.set_ylabel('Lowest Level Wind Speed [m/s]', color=color, fontsize=title_fontsize)
for c in range(num_case):
   ax3.plot(bin_list[c], np.ma.masked_invalid(wnd_list[c]), color=color, **kwargs[c])
ax3.tick_params(axis='x', labelsize=label_fontsize)
ax3.tick_params(axis='y', labelsize=label_fontsize, labelcolor=color)
ax3.set_xlim(xmin=0.1)
# ax3.set_ylim(ymin=0)
ax3.set_xscale('log')

ax4 = ax3.twinx()  # instantiate a second Axes that shares the same x-axis

color = 'magenta'
ax4.set_xlabel('Precipitation [mm/day]', fontsize=title_fontsize)
ax4.set_ylabel('Lowest Level Temperature [K]', color=color, fontsize=title_fontsize)
for c in range(num_case):
   ax4.plot(bin_list[c], np.ma.masked_invalid(tmp_list[c]), color=color, **kwargs[c])
ax4.tick_params(axis='x', labelsize=label_fontsize)
ax4.tick_params(axis='y', labelsize=label_fontsize, labelcolor=color)
ax4.set_xlim(xmin=0.1)
# ax4.set_ylim(ymin=0)
ax4.set_xscale('log')



color = 'black'
ax_pct.set_xlabel('Precipitation [mm/day]', fontsize=title_fontsize)
ax_pct.set_ylabel('Frequency [%]', fontsize=title_fontsize)
for c in range(num_case):
   ax_pct.plot(bin_list[c], np.ma.masked_invalid(pct_list[c]), color=color, **kwargs[c])
ax_pct.tick_params(axis='x', labelsize=label_fontsize)
ax_pct.tick_params(axis='y', labelsize=label_fontsize)
ax_pct.set_xlim(xmin=0.1)
ax_pct.set_ylim(ymin=0)
ax_pct.set_xscale('log')

# Create proxy handles for legend
from matplotlib.lines import Line2D
p1 = Line2D([0], [0], color='black', lw=2, linestyle='solid')
p2 = Line2D([0], [0], color='black', lw=2, linestyle='dashed')
ax1.legend( [p1, p2], [ opts_list[0]['name'], opts_list[1]['name'] ], \
            loc='upper left', borderpad=1.0, fontsize=title_fontsize )
            # bbox_to_anchor=(0.05, 0.95)

ax1   .annotate('a',xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize')
ax3   .annotate('b',xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize')
ax_pct.annotate('c',xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize')
        # fontsize='medium', verticalalignment='top', fontfamily='serif',
        # bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))

fig.tight_layout()  # otherwise the right y-label is slightly clipped
fig.savefig(f'{fig_file}.png', dpi=200)

run_cmd(f'magick convert -trim +repage {fig_file}.png {fig_file}.png')

#---------------------------------------------------------------------------------------------------

