# v2 is a cleaned up version of v1 - functionality is mostly the same
import os, glob, ngl, subprocess as sp, numpy as np, xarray as xr, dask, copy, string, cmocean, pandas as pd
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
host = hc.get_host()
#-------------------------------------------------------------------------------
''' mapping commands

ncremap -G ttl='BEST grid'#latlon==180,360#lat_typ=uni#lon_typ=180_wst -g grid_files/BEST_180x360_scrip.20241001.nc
SRC_GRID_FILE=grid_files/ne30pg2_scrip.nc
DST_GRID_FILE=grid_files/BEST_180x360_scrip.20241001.nc
MAP_FILE=map_files/map_ne30pg2_to_180x360_BEST_traave_20241001.nc
ncremap -a traave --src_grd=${SRC_GRID_FILE} --dst_grd=${DST_GRID_FILE} --map_file=${MAP_FILE}

SRC_GRID_FILE=grid_files/BEST_180x360_scrip.20241001.nc
DST_GRID_FILE=grid_files/ne30pg2_scrip.nc
MAP_FILE=map_files/map_180x360_BEST_to_ne30pg2_traave_20241001.nc
ncremap -a traave --src_grd=${SRC_GRID_FILE} --dst_grd=${DST_GRID_FILE} --map_file=${MAP_FILE}

SRC_FILE=/global/cfs/cdirs/m3312/whannah/obs_data/BEST/Land_and_Ocean_LatLong1.nc
DST_FILE=/global/cfs/cdirs/m3312/whannah/obs_data/BEST/Land_and_Ocean_LatLong1.remap_ne30pg2.nc
MAP_FILE=map_files/map_180x360_BEST_to_ne30pg2_traave_20241001.nc
ncremap -m ${MAP_FILE} -i ${SRC_FILE} -o ${DST_FILE}

ncremap -G ttl='HadCRU grid'#latlon==36,72#lat_typ=uni#lon_typ=grn_ctr -g grid_files/HadCRU_36x72_scrip.20241001.nc
SRC_GRID_FILE=grid_files/ne30pg2_scrip.nc
DST_GRID_FILE=grid_files/HadCRU_36x72_scrip.20241001.nc
MAP_FILE=map_files/map_ne30pg2_to_36x72_traave_20241001.nc
ncremap -a traave --src_grd=${SRC_GRID_FILE} --dst_grd=${DST_GRID_FILE} --map_file=${MAP_FILE}

DATA_ROOT=/global/cfs/cdirs/m3312/whannah/obs_data/HadCRU
SRC_DATA=${DATA_ROOT}/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc
DST_DATA=${DATA_ROOT}/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.remap_ne30pg2.nc
MAP_FILE=map_files/map_36x72_to_ne30pg2_traave_20241001.nc
ncremap -m ${MAP_FILE} -i ${SRC_DATA} -o ${DST_DATA}
'''
#-------------------------------------------------------------------------------
def run_cmd(cmd):
   msg = hc.tcolor.GREEN + cmd + hc.tcolor.ENDC ; print(f'\n{msg}')
   os.system(cmd); return
#-------------------------------------------------------------------------------
name,case,case_dir,case_sub = [],[],[],[]
obs_flag = []
def add_case(case_in,n=None,p=None,s=None,g=None,c=None,obs=False):
   global name,case,case_dir,case_sub
   tmp_name = case_in if n is None else n
   case.append(case_in); name.append(tmp_name); 
   case_dir.append(p); case_sub.append(s);
   obs_flag.append(obs)
#-------------------------------------------------------------------------------
var,lev_list,var_str = [],[],[]
def add_var(var_in,lev=-1,name=None): 
   var.append(var_in); lev_list.append(lev)
   var_str.append(var_in if name is None else name)
#-------------------------------------------------------------------------------
scrip_file_path  = os.getenv('HOME')+'/E3SM/data_grid/ne30pg2_scrip.nc'
tmp_path_hst_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_v2  = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
tmp_path_co2_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_sub = 'archive/atm/hist'


# add_case('ERA5',obs=True)
add_case('BEST',obs=True)

add_case('v2.LR.historical_0101',                                  n='E3SMv2',  c='cyan', p=tmp_path_hst_v2, s='archive/atm/hist')
add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',c='blue', p=tmp_path_hst_mmf,s='archive/atm/hist')

#-------------------------------------------------------------------------------

fig_file,fig_type = 'figs/F02-map-sfc-temperature','png'
tmp_file_head = 'map-TS'

# add_var('TS',name='Tsfc')
add_var('TS',name='Sfc Air Temp')

#-------------------------------------------------------------------------------
htype,yr1,yr2 = 'ha',1979,2014
date_str1 = f'{yr1}-{yr2} mean'
date_str2 = f'{yr2}-{yr1} diff'

plot_diff = True

recalculate = False

print_stats          = True
var_x_case           = False
num_plot_col         = 1
use_common_label_bar = False

#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)
diff_base = 0
if 'first_file'      not in locals(): first_file = 0
if 'num_files'       not in locals(): num_files = 0
if 'scrip_file_path' not in locals(): scrip_file_path = None
#---------------------------------------------------------------------------------------------------
wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*(num_var*2*num_case)

res = hs.res_contour_fill_map()
res.lbLabelFontHeightF           = 0.012
res.tmXBOn                       = False
res.tmYLOn                       = False
res.mpCenterLonF                 = 0
res.pmTickMarkDisplayMode        = 'Never'
res.mpProjection                 = 'Robinson'

subtitle_font_height = 0.01

tval_crit = 2.042 # 95% 2-tail t-statistic critical value w/ 30 dof
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   data_list = []
   trend_list = []
   trend_tval_list = []
   stddev_list = []
   glb_avg_list = []
   lat_list,lon_list = [],[]
   if 'lev_list' in locals(): lev = lev_list[v]
   #----------------------------------------------------------------------------
   scrip_ds = xr.open_dataset('grid_files/ne30pg2_scrip.nc')
   area = scrip_ds['grid_area'].rename({'grid_size':'ncol'})
   #----------------------------------------------------------------------------
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)
      #-------------------------------------------------------------------------
      tmp_file = f'data/{tmp_file_head}.{case[c]}.{var[v]}.nc'
      if recalculate:
         #----------------------------------------------------------------------
         # if case[c]=='ERA5':
         #    obs_root = '/global/cfs/cdirs/m3312/whannah/2023-CPL/ERA5'
         #    obs_file = f'{obs_root}/ERA5_ANN_198501_201412_climo.ne30pg2.nc'
         #    ds = xr.open_dataset(obs_file)
         #    data = ds['ts'] 
         #    print()
         #    print(data)
         #    print()
         #    exit()
         if case[c]=='BEST':
            obs_file = '/global/cfs/cdirs/m3312/whannah/obs_data/BEST/Land_and_Ocean_LatLong1.remap_ne30pg2.nc'
            ds = xr.open_dataset(obs_file)#.rename({'latitude':'lat','longitude':'lon'})
            data = ds['temperature']
            time_list = [None]*len(data['time'])
            # convert anomalies to absolute temperature and build time coord
            for t,tt in enumerate(data['time'].values):
               yr_val = int(np.floor(tt))
               mn_ind = int(np.floor((tt-yr_val)*12))
               time_list[t] = np.datetime64(f'{yr_val}-{(mn_ind+1):02d}')
               data[t,:] = data[t,:] + ds['climatology'].isel(month_number=mn_ind)
            time = xr.DataArray(time_list,dims=('time'))
            time = time.assign_coords(time=time)
            month_length = time.dt.days_in_month
            wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
            data = data.assign_coords(time=time)
            data = data.resample(time='Y').mean(dim='time')
            time = time.resample(time='Y').mean(dim='time')
            # subset in time
            for t,y in enumerate(data['time.year'].values):
               if y==yr1: t_beg=t
               if y==yr2: t_end=t
            data = data.isel(time=slice(t_beg,t_end+1))
         else:
            file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/*.eam.{htype}.*'
            file_list_all = sorted(glob.glob(file_path))
            file_list1 = [] # data for anomalies
            file_list2 = [] # data for baseline
            for f in range(len(file_list_all)):
               yr = int(file_list_all[f][-7:-7+4])
               if yr>= yr1 and yr<= yr2: file_list1.append(file_list_all[f])
            if file_list1==[]: exit(f'\nERROR: no files found for file_path:\n{file_path}\n')   
            ds_main = xr.open_mfdataset( file_list1 )
            # ds_base = xr.open_mfdataset( file_list2 )
            # data = ds_main[var[v]] - ds_base[var[v]].mean(dim='time')
            data = ds_main[var[v]] - 273.15
         #----------------------------------------------------------------------
         # simple and fast method for regression coeff and intercept
         px = data['time.year']
         py = data
         a = xr.cov( px, py, dim='time' ) / px.var(dim='time')
         b = py.mean(dim='time') - a*px.mean(dim='time')
         #----------------------------------------------------------------------
         # print()
         # print(f'a: {     a.min().values}  /  {     a.max().values}')
         # print(f'b: {     b.min().values}  /  {     b.max().values}')
         # print(f'y: {    py.min().values}  /  {    py.max().values}  /  {    py.mean().values}')
         # print(f'x: {(a*px).min().values}  /  {(a*px).max().values}  /  {(a*px).mean().values}')
         # print()
         # exit()
         #----------------------------------------------------------------------
         # calculate std error and t-statistic for trend
         dof = yr2 - yr1 - 2
         y_r = (a*px+b) - py
         x_r = px - px.mean()
         yr_variance = (np.square(y_r)).sum(dim='time')
         xr_variance = (np.square(x_r)).sum(dim='time')
         stderr_sq = (1/dof) * yr_variance / xr_variance
         stderr = np.sqrt( stderr_sq )
         trend_tval = a / stderr
         #----------------------------------------------------------------------
         print()
         print(' '*6+f'py: {np.min(py.values):8.4f}  /  {np.max(py.values):8.4f}')
         print(' '*6+f'px: {np.min(px.values):8.4f}  /  {np.max(px.values):8.4f}')
         print(' '*6+f'yr: {np.min(y_r.values):8.4f}  /  {np.max(y_r.values):8.4f}')
         print(' '*6+f'xr: {np.min(x_r.values):8.4f}  /  {np.max(x_r.values):8.4f}')
         print(' '*6+f'trnd min/max: {np.min(a.values)         :8.4f}  /  {np.max(a.values)         :8.4f}')
         print(' '*6+f'serr min/max: {np.min(stderr.values)    :8.4f}  /  {np.max(stderr.values)    :8.4f}')
         print(' '*6+f'tval min/max: {np.min(trend_tval.values):8.4f}  /  {np.max(trend_tval.values):8.4f}')
         print()
         # exit()
         #----------------------------------------------------------------------
         # average over time dimension
         # hc.print_time_length(data.time,indent=' '*6)
         stddev = data.std(dim='time')
         data = data.mean(dim='time')
         #----------------------------------------------------------------------
         # Write to file 
         if os.path.isfile(tmp_file) : os.remove(tmp_file)
         tmp_ds = xr.Dataset()
         tmp_ds[var[v]]   = data
         tmp_ds['trend']  = a
         tmp_ds['stddev'] = stddev
         tmp_ds['dof']    = dof
         tmp_ds['stderr'] = stderr
         tmp_ds['trend_tval'] = trend_tval
         
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      #-------------------------------------------------------------------------
      tmp_ds = xr.open_dataset( tmp_file )
      data       = tmp_ds[var[v]]
      stddev     = tmp_ds['stddev']
      trend      = tmp_ds['trend']
      trend_tval = tmp_ds['trend_tval']
      #-------------------------------------------------------------------------
      data = data + 273.15
      #-------------------------------------------------------------------------
      if print_stats: 
         hc.print_stat(data,name=var[v],stat='naxsh',indent=' '*6,compact=True)
         hc.print_stat(trend,name='trend',stat='naxsh',indent=' '*6,compact=True)
      #-------------------------------------------------------------------------
      # Calculate area weighted global mean
      gbl_mean = ( (data*area).sum() / area.sum() ).values 
      print(hc.tcolor.CYAN+' '*6+f'Area Weighted Global Mean : {gbl_mean:6.4}'+hc.tcolor.ENDC)
      glb_avg_list.append(gbl_mean)
      #-------------------------------------------------------------------------
      # append to data lists
      data_list.append( data.values )
      trend_list.append( trend.values )
      trend_tval_list.append( trend_tval.values )
      stddev_list.append( stddev.values )
      #-------------------------------------------------------------------------
      # save baseline for diff map
      if plot_diff :
         if c==diff_base:
            data_baseline = data.copy().values
   #----------------------------------------------------------------------------
   # calculate common limits for consistent contour levels
   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])

   if plot_diff:
      rmse_list = []
      tmp_data_list = []
      for c in range(num_case): 
         rmse_tmp = np.sqrt( np.mean( np.square( data_list[c] - data_baseline )))
         rmse_list.append(rmse_tmp)
         tmp_data_list.append( data_list[c] - data_baseline )
      diff_data_min = np.min([np.nanmin(d) for d in tmp_data_list])
      diff_data_max = np.max([np.nanmax(d) for d in tmp_data_list])
   #----------------------------------------------------------------------------
   # Plot averaged data
   for c in range(num_case):
      #-------------------------------------------------------------------------
      # Set color palette and levels
      tres = copy.deepcopy(res)
      main_colormap = np.array( cmocean.cm.thermal(np.linspace(0,1,256)) )
      # main_colormap = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
      diff_colormap = np.array( cmocean.cm.diff(np.linspace(0,1,256)) )
      # diff_colormap = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )
      # diff_colormap = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
      # trend_colormap = np.array( cmocean.cm.curl(np.linspace(0,1,256)) )
      # trend_colormap = np.array( cmocean.cm.diff(np.linspace(0,1,256)) )
      trend_colormap = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
      # trend_colormap = 'BlueWhiteOrangeRed'

      main_levels  = np.arange(-60,34+6,6) + 273#.15
      # diff_levels  = np.arange(-11,11+2,2)
      diff_levels  = np.arange(-11,11+1,1)
      # diff_levels  = np.arange(-10.5,10.5+1.5,1.5)
      # trend_levels = np.arange(-0.13,0.13+0.02,0.02)[:-1]
      trend_levels = np.arange(-0.14,0.14+0.02,0.02)[:-1]
      #-------------------------------------------------------------------------
      tres.sfXArray = scrip_ds['grid_center_lon'].values
      tres.sfYArray = scrip_ds['grid_center_lat'].values

      ip1 = v*1*num_case+c if var_x_case else c*(num_var*2)+v+0
      ip2 = v*2*num_case+c if var_x_case else c*(num_var*2)+v+1

      if not plot_diff  or (plot_diff and c==diff_base) : 

         tres.cnFillPalette = main_colormap
         tres.cnLevelSelectionMode = 'ExplicitLevels'
         tres.cnLevels = main_levels
         plot[ip1] = ngl.contour_map(wks,np.ma.masked_invalid(data_list[c]),tres)
            
         tres.cnFillPalette = trend_colormap
         tres.cnLevels      = trend_levels
         plot[ip2] = ngl.contour_map(wks,np.ma.masked_invalid(trend_list[c]),tres)
         
         hs.set_subtitles( wks, plot[ip1],
                           left_string=name[c], 
                           center_string=f'Mean: {glb_avg_list[c]:5.2f} K',
                           right_string=f'{var_str[v]} Mean',
                           font_height=subtitle_font_height,
                           right_sub_string=f'{yr1}-{yr2}')

         gbl_mean_trend = ( (trend*area).sum() / area.sum() ).values 
         hs.set_subtitles( wks, plot[ip2], 
                           left_string=name[c],
                           center_string=f'Mean: {gbl_mean_trend:5.2f} K/yr', 
                           right_string=f'{var_str[v]} Trend',
                           font_height=subtitle_font_height,
                           right_sub_string=f'{yr1}-{yr2}')

      #-------------------------------------------------------------------------
      # create difference plot
      if plot_diff and c!=diff_base :
         
         data_list[c] = data_list[c] - data_baseline
         
         tres.cnLevelSelectionMode = "ExplicitLevels"
         if hasattr(tres,'cnLevels') : del tres.cnLevels

         tres.cnFillPalette = diff_colormap
         tres.cnLevelSelectionMode = "ExplicitLevels"
         tres.cnLevels = diff_levels
         plot[ip1] = ngl.contour_map(wks,np.ma.masked_invalid(data_list[c]),tres)

         tres.cnFillPalette = trend_colormap
         tres.cnLevels      = trend_levels
         plot[ip2] = ngl.contour_map(wks,np.ma.masked_invalid(trend_list[c]),tres)
         
         glb_mean_diff = glb_avg_list[c] - glb_avg_list[diff_base]
         hs.set_subtitles( wks, plot[ip1],
                           left_string=name[c], 
                           center_string=f'Mean: {glb_mean_diff:5.2f} K',
                           right_string=f'{var_str[v]} Bias',
                           font_height=subtitle_font_height,
                           right_sub_string=f'{yr1}-{yr2}')

         gbl_mean_trend = ( (trend*area).sum() / area.sum() ).values 
         hs.set_subtitles( wks, plot[ip2], 
                           left_string=name[c],
                           center_string=f'Mean: {gbl_mean_trend:5.2f} K/yr', 
                           right_string=f'{var_str[v]} Trend',
                           font_height=subtitle_font_height,
                           right_sub_string=f'{yr1}-{yr2}')

      #-------------------------------------------------------------------------
      sres = hs.res_stippling()
      sres.cnFillScaleF   = 0.5
      sres.cnFillDotSizeF = 0.0015
      sres.sfXArray       = scrip_ds['grid_center_lon'].values
      sres.sfYArray       = scrip_ds['grid_center_lat'].values

      # NOTE - when I checked the bias against the confidence interval it was
      # pretty much all significant, which seems suspicious...

      # # overlay stippling to indicate bias confidence interval
      # if plot_diff and c!=diff_base :
      #    var_test = np.square(stddev_list[c])
      #    var_base = np.square(stddev_list[0])
      #    bias_conf = tval_crit * np.sqrt( np.absolute(var_test-var_base)/(yr2-yr1) )
      #    sig_mask = np.where( np.absolute(data_list[c]) > bias_conf, 1, 0 )
      #    ngl.overlay( plot[ip1], ngl.contour(wks,sig_mask,sres) )

      # overlay stippling to indicate where trend is significant
      sig_mask = np.where( np.absolute(trend_tval_list[c]) > tval_crit, 1, 0 )
      ngl.overlay( plot[ip2], ngl.contour(wks,sig_mask,sres) )

#---------------------------------------------------------------------------------------------------
# Finalize plot
pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent       = 5
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.012

layout = [num_var*2,num_case] if var_x_case else [num_case,num_var*2]

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
