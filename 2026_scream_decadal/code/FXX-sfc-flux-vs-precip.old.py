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
hx_data_root = '/global/cfs/cdirs/e3smdata/simulations/scream-decadal/DYNAMO_MJO_hindcasts'
# hx_data_sub  = 'data_remap_73x144'
hx_data_sub  = 'run'

htype = 'output.scream.2D.1hr.ne30pg2.AVERAGE.nhours_x1'

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

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

plot = [None]*3

res = hs.res_xy()
res.vpwidthF = 0.5
res.vpHeightF = 0.2
res.tmYLLabelFontHeightF         = 0.009
res.tmXBLabelFontHeightF         = 0.009
res.tiXAxisFontHeightF           = 0.012
res.tiYAxisFontHeightF           = 0.012

res.xyLineThicknessF = 10

res.trXMinF = 0.1
res.trXMaxF = 500

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
num_pcp_bins = 65#51
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
      #-------------------------------------------------------------------------
      if file_list_pcp==[]: print();print(f'file_path_pcp: {file_path_pcp}');exit('ERROR: file_path_pcp is empty!')
      if file_list_flx==[]: print();print(f'file_path_flx: {file_path_flx}');exit('ERROR: file_path_flx is empty!')
      #-------------------------------------------------------------------------
      ds_pcp = xr.open_mfdataset(file_list_pcp)
      ds_flx = xr.open_mfdataset(file_list_flx)
      data_pcp = ds_pcp[case_opts['var_pcp']]
      data_shf = ds_flx[case_opts['var_shf']]
      data_lhf = ds_flx[case_opts['var_lhf']]
      #-------------------------------------------------------------------------
      if len(data_pcp.time)==0: raise ValueError('ERROR - something went wrong - data_pcp time dimension is zero!')
      if len(data_shf.time)==0: raise ValueError('ERROR - something went wrong - data_shf time dimension is zero!')
      if len(data_lhf.time)==0: raise ValueError('ERROR - something went wrong - data_lhf time dimension is zero!')
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
      #-------------------------------------------------------------------------
      data_pcp.load()
      data_shf.load()
      data_lhf.load()
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
      for b in range(num_pcp_bins):
         condition = ( data_pcp.values>=pcp_bin_bot[b] ) \
                    &( data_pcp.values <pcp_bin_top[b] )
         condition = xr.DataArray( condition, coords=data_pcp.coords )
         bin_cnt[b] = condition.sum()
         if bin_cnt[b]>0:
            bin_shf[b] = data_shf.where(condition,drop=True).mean(skipna=True)#.sum() / (bin_spc_pct/1e2)
            bin_lhf[b] = data_lhf.where(condition,drop=True).mean(skipna=True)#.sum() / (bin_spc_pct/1e2)
            bin_pcp[b] = data_pcp.where(condition,drop=True).mean(skipna=True)#.sum() / (bin_spc_pct/1e2)
      
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
   #----------------------------------------------------------------------------
   # hc.print_stat(bin_ctr,name='bin_ctr',compact=True,indent=' '*6)
   hc.print_stat(bin_pcp,name='bin_pcp',compact=True,indent=' '*6)
   hc.print_stat(bin_shf,name='bin_shf',compact=True,indent=' '*6)
   hc.print_stat(bin_lhf,name='bin_lhf',compact=True,indent=' '*6)
   # hc.print_stat(bin_pct,name='bin_pct',compact=True,indent=' '*6)
   #----------------------------------------------------------------------------
   bin_list.append( np.ma.masked_invalid(bin_ctr) )
   shf_list.append( np.ma.masked_invalid(bin_shf) )
   lhf_list.append( np.ma.masked_invalid(bin_lhf) )
   pct_list.append( np.ma.masked_invalid(bin_pct) )

#-------------------------------------------------------------------------------
# calculate limits for common color bar
# data_min = np.min([np.nanmin(d) for d in data_list])
# data_max = np.max([np.nanmax(d) for d in data_list])
# if data_min==data_max: raise ValueError(hc.tclr.RED+'WARNING: Difference is zero!'+hc.tclr.END)
#-------------------------------------------------------------------------------
tres = copy.deepcopy(res)

clr,dsh = [None]*num_case,[None]*num_case
for c in range(num_case):
   clr[c] = opts_list[c]['clr']
   dsh[c] = opts_list[c]['dsh']

# tres.xyLineColors   = clr
tres.xyDashPatterns = dsh

# fres = copy.deepcopy(tres); #ores.trYMinF,ores.trYMaxF = -150,0
pres = copy.deepcopy(tres); #pres.trYMinF,pres.trYMaxF = 0,5

#-------------------------------------------------------------------------------
# Create plot

# ip = v*num_wave+w if var_x_case else w*num_var+v

# tres.tiXAxisString = f'Longitude'
# tres.tiYAxisString = f'{var_list[v]} Variance [{unit_list[v]}]'


pres.tiYAxisString = 'Frequency [%]'
pres.tiXAxisString = 'Precipitation [mm/day]'

# adjust panel widths so they are consistent
# ores.vpWidthF = 0.6*1.15

sres = copy.deepcopy(tres);
sres.xyLineColor = 'red';
# sres.tiYAxisString = '[W/m2]'
# sres.tmXTOn = False
# sres.tmXUseBottom = False
sres.tmYROn = True
sres.tmYLOn = False
# sres.tmYUseLeft = False
sres.tmYRLabelFontColor = 'red'
sres.tiYAxisFontColor = 'red'
sres.tmYRLabelsOn = True
sres.tmYRMajorOutwardLengthF      = 0.
sres.tmYRMinorOutwardLengthF      = 0.
sres.tiYAxisString = 'Sensible Heat Flux [W/m2]'
sres.tiYAxisSide = "Right"

lres = copy.deepcopy(tres);
lres.xyLineColor = 'blue';
lres.tiYAxisString = 'Latent Heat Flux [W/m2]'
lres.tmYLLabelFontColor = 'blue'
lres.tiYAxisFontColor = 'blue'
lres.tmYROn = False
lres.tmYUseLeft = False

# shf_max = np.max([np.nanmax(d) for d in shf_list])
# lhf_max = np.max([np.nanmax(d) for d in lhf_list])
# for i in range(len(shf_list)): shf_list[i] = shf_list[i]*lhf_max/shf_max

# plot_shf = ngl.xy(wks, np.stack(bin_list) , np.ma.masked_invalid(np.stack(shf_list)) ,sres)
# plot_lhf = ngl.xy(wks, np.stack(bin_list) , np.ma.masked_invalid(np.stack(lhf_list)) ,lres)

# plot[0] = plot_lhf
# ngl.overlay(plot[0], plot_shf)
# plot[1] = ngl.xy(wks, np.stack(bin_list) , np.ma.masked_invalid(np.stack(shf_list)) ,sres)
# plot[2] = ngl.xy(wks, np.stack(bin_list) , np.ma.masked_invalid(np.stack(pct_list)) ,pres)

# #---------------------------------------------------------------------------------------------------
# # Add legend

# lgres = ngl.Resources()
# lgres.vpWidthF,lgres.vpHeightF           = 0.08,0.12
# lgres.lgLabelFontHeightF = 0.01
# # lgres.vpWidthF,lgres.vpHeightF           = 0.18,0.15
# # lgres.lgLabelFontHeightF = 0.014
# lgres.lgLineLabelsOn     = False
# lgres.lgLineThicknessF   = 12
# lgres.lgLabelJust        = 'CenterLeft'
# lgres.lgLineColors       = clr
# lgres.lgDashIndexes      = dsh
# # lgres.lgMonoDashIndex    = True

# lgd_lbl = [None]*num_case
# for c in range(num_case): lgd_lbl[c] = f'  {opts_list[c]["name"]}'

# # pid = ngl.legend_ndc(wks, num_case, lgd_lbl, 0.6, 0.97, lgres)
# pid = ngl.legend_ndc(wks, num_case, lgd_lbl, 0.62, 0.89, lgres)

# #---------------------------------------------------------------------------------------------------
# # Finalize plot
# #---------------------------------------------------------------------------------------------------

# # layout = [num_var,num_wave] if var_x_case else [num_wave,num_var]
# num_plot_col = 1
# layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

# pnl_res = hs.setres_panel()
# # pnl_res.nglPanelYWhiteSpacePercent = 2
# # pnl_res.nglPanelXWhiteSpacePercent = 5
# pnl_res.nglPanelTop = 0.9
# pnl_res.nglPanelBottom = 0.1
# pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
# pnl_res.nglPanelFigureStringsJust        = "TopLeft"
# pnl_res.nglPanelFigureStringsFontHeightF = 0.012

# ngl.panel(wks,plot,layout,pnl_res)
# hc.trim_png(fig_file)

# ngl.end()

#-------------------------------------------------------------------------------
# Create plots separately

fig_file_shf = f'{fig_file}.shf'
fig_file_lhf = f'{fig_file}.lhf'

# wkres.wkBackgroundColor = 'transparent'
wks_shf = ngl.open_wks(fig_type,fig_file_shf,wkres)
# wkres.wkBackgroundColor = 'white'
wks_lhf = ngl.open_wks(fig_type,fig_file_lhf,wkres)


plot_shf = ngl.xy(wks_shf, np.stack(bin_list) , np.ma.masked_invalid(np.stack(shf_list)) ,sres)
plot_lhf = ngl.xy(wks_lhf, np.stack(bin_list) , np.ma.masked_invalid(np.stack(lhf_list)) ,lres)

ngl.draw(wks_shf); ngl.frame(wks_shf);
ngl.draw(wks_lhf); ngl.frame(wks_lhf);

# hc.trim_png(fig_file_shf)
# hc.trim_png(fig_file_lhf)


run_cmd(f'magick convert {fig_file_shf}.png -transparent "white" {fig_file_shf}.png')
run_cmd(f'magick convert {fig_file_lhf}.png  {fig_file_shf}.png -gravity Center -geometry +50+0 -composite {fig_file}.png')
# run_cmd(f'magick convert -trim +repage {fig_file}.png {fig_file}.png')

ngl.end()

'''
-geometry 2048x2048+30+5
-resize 64x64
'''


#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
