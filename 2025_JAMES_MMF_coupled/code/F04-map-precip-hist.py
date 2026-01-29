# v2 is a cleaned up version of v1 - functionality is mostly the same
import os, glob, ngl, subprocess as sp, numpy as np, xarray as xr, dask, copy, string, cmocean, pandas as pd
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
host = hc.get_host()
#-------------------------------------------------------------------------------
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


add_case('IMERG',obs=True)

add_case('v2.LR.historical_0101',                                  n='E3SMv2',  c='cyan', p=tmp_path_hst_v2, s='archive/atm/hist')
add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',c='blue', p=tmp_path_hst_mmf,s='archive/atm/hist')

#-------------------------------------------------------------------------------

fig_file,fig_type = 'figs/F04-map-precip','png'
tmp_file_head = 'map-PRECT'

add_var('PRECT',name='Precip')

unit_str = 'mm/day'

#-------------------------------------------------------------------------------
htype,yr1,yr2 = 'ha',2001,2014
date_str1 = f'{yr1}-{yr2} mean'
date_str2 = f'{yr2}-{yr1} diff'

plot_diff = True

recalculate = False

print_stats          = True
var_x_case           = True
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
plot = [None]*(num_var*1*num_case)

res = hs.res_contour_fill_map()
res.lbLabelFontHeightF           = 0.015
res.tmXBOn                       = False
res.tmYLOn                       = False
res.mpCenterLonF                 = 180
res.pmTickMarkDisplayMode        = 'Never'
res.mpProjection                 = 'Robinson'

res.cnLevelSelectionMode = 'ExplicitLevels'
# res.cnLevelSelectionMode = 'AutomaticLevels'

subtitle_font_height = 0.008
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   data_list = []
   delta_list = []
   trend_list = []
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
         if case[c]=='IMERG':
            file_list = sorted(glob.glob('/pscratch/sd/w/whannah/Obs/IMERG/monthly_remap_ne30pg2/*'))
            ds = xr.open_mfdataset(file_list)#.rename({'latitude':'lat','longitude':'lon'})
            
            time_index_list = []
            # convert anomalies to absolute temperature and build time coord
            for t,yr in enumerate(ds['time.year'].values):
               if yr>=yr1 and yr<=yr2: time_index_list.append(t)
            data = ds['precipitation'].isel(time=time_index_list)
            time = data['time']
            data = data.resample(time='Y').mean(dim='time')
            time = time.resample(time='Y').mean(dim='time')
            data = data * 24.
         else:
            file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/*.eam.{htype}.*'
            file_list_all = sorted(glob.glob(file_path))
            file_list = [] # data for anomalies
            for f in range(len(file_list_all)):
               yr = int(file_list_all[f][-7:-7+4])
               if yr>= yr1 and yr<= yr2: file_list.append(file_list_all[f])
            if file_list==[]: exit(f'\nERROR: no files found for file_path:\n{file_path}\n')   
            ds = xr.open_mfdataset( file_list )
            data = ds['PRECC'] + ds['PRECL']
            data = data * 86400.*1000.
         #-------------------------------------------------------------------------
         # simple and fast method for regression coeff and intercept
         px = data['time.year']
         py = data
         trend = xr.cov( px, py, dim='time' ) / px.var()
         #----------------------------------------------------------------------
         # average over time dimension
         hc.print_time_length(data.time,indent=' '*6)
         delta = data.isel(time=-1) - data.isel(time=0)
         data = data.mean(dim='time')
         #----------------------------------------------------------------------
         # Write to file 
         if os.path.isfile(tmp_file) : os.remove(tmp_file)
         tmp_ds = xr.Dataset()
         tmp_ds[var[v]] = data
         tmp_ds['delta'] = delta
         tmp_ds['trend'] = trend
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      #-------------------------------------------------------------------------
      tmp_ds = xr.open_dataset( tmp_file )
      data = tmp_ds[var[v]]
      delta = tmp_ds['delta']
      trend = tmp_ds['trend']
      #-------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------
      if print_stats: 
         hc.print_stat(data,name=var[v],stat='naxsh',indent=' '*6,compact=True)
         hc.print_stat(delta,name='delta',stat='naxsh',indent=' '*6,compact=True)
         hc.print_stat(trend,name='trend',stat='naxsh',indent=' '*6,compact=True)
      #-------------------------------------------------------------------------
      # Calculate area weighted global mean
      gbl_mean = ( (data*area).sum() / area.sum() ).values 
      print(hc.tcolor.CYAN+' '*6+f'Area Weighted Global Mean : {gbl_mean:6.4}'+hc.tcolor.ENDC)
      glb_avg_list.append(gbl_mean)
      #-------------------------------------------------------------------------
      # append to data lists
      data_list.append( data.values )
      delta_list.append( delta.values )
      trend_list.append( trend.values )
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
      
      main_colormap = 'MPL_viridis'
      # main_colormap = np.array( cmocean.cm.thermal(np.linspace(0,1,256)) )
      # main_colormap = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
      
      # diff_colormap = np.array( cmocean.cm.diff(np.linspace(0,1,256)) )
      # diff_colormap = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )
      diff_colormap = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
      

      main_levels  = np.arange(2,16+2,2)
      # diff_levels  = np.array([-5,-4,-3,-2,-1,1,2,3,4,5])
      diff_levels  = np.arange(-5,5+1,1)
      # diff_levels  = np.arange(-5.5,5.5+1,1)
      # diff_levels  = np.arange(-4.4,4.4+0.8,0.8)[:-1]
      
      # dc=2.0; cmax=6*dc-dc/2.; diff_levels  = np.arange(-cmax,cmax+dc,dc)
      # dc=0.8; cmax=7*dc-dc/2.; diff_levels  = np.arange(-cmax,cmax+dc,dc)
      # if len(diff_levels)%2 != 0 : diff_levels = diff_levels[:-1]
      
      # delta_levels = np.arange(-0.5,0.5+0.05,0.05)#[:-1]
      # print(diff_levels)
      # exit()
      #-------------------------------------------------------------------------
      tres.sfXArray = scrip_ds['grid_center_lon'].values
      tres.sfYArray = scrip_ds['grid_center_lat'].values

      # ip1 = v*1*num_case+c if var_x_case else c*(num_var*2)+v+0
      # ip2 = v*2*num_case+c if var_x_case else c*(num_var*2)+v+1
      ip1 = v*1*num_case+c if var_x_case else c*(num_var*1)+v+0

      if not plot_diff  or (plot_diff and c==diff_base) : 

         tres.cnFillPalette = main_colormap
         tres.cnLevels = main_levels
         plot[ip1] = ngl.contour_map(wks,np.ma.masked_invalid(data_list[c]),tres)
            
         # tres.cnFillPalette = delta_colormap
         # # tres.cnLevels = delta_levels
         # # plot[ip2] = ngl.contour_map(wks,np.ma.masked_invalid(delta_list[c]),tres)
         # plot[ip2] = ngl.contour_map(wks,np.ma.masked_invalid(trend_list[c]),tres)
         
         hs.set_subtitles( wks, plot[ip1],
                           left_string=name[c], 
                           center_string=f'Mean: {glb_avg_list[c]:5.2f} {unit_str}',
                           right_string=f'{var_str[v]} Mean',
                           font_height=subtitle_font_height,
                           right_sub_string=f'{yr1}-{yr2}')

         # gbl_mean_trend = ( (trend*area).sum() / area.sum() ).values 
         # hs.set_subtitles( wks, plot[ip2], 
         #                   left_string=name[c],
         #                   center_string=f'Mean: {gbl_mean_trend:5.2f} {unit_str}', 
         #                   right_string=f'{var_str[v]} Trend',
         #                   font_height=subtitle_font_height,
         #                   right_sub_string=f'{yr1}-{yr2}')

      #-------------------------------------------------------------------------
      # create difference plot
      if plot_diff and c!=diff_base :

         data_list[c] = data_list[c] - data_baseline
         
         if hasattr(tres,'cnLevels') : del tres.cnLevels

         tres.cnFillPalette = diff_colormap
         tres.cnLevels = diff_levels
         plot[ip1] = ngl.contour_map(wks,np.ma.masked_invalid(data_list[c]),tres)

         # tres.cnFillPalette = delta_colormap
         # tres.cnLevels = delta_levels
         # # plot[ip2] = ngl.contour_map(wks,np.ma.masked_invalid(delta_list[c]),tres) 
         # plot[ip2] = ngl.contour_map(wks,np.ma.masked_invalid(trend_list[c]),tres)
         
         glb_mean_diff = glb_avg_list[c] - glb_avg_list[diff_base]
         hs.set_subtitles( wks, plot[ip1],
                           left_string=name[c], 
                           center_string=f'Mean: {glb_mean_diff:5.3f} {unit_str}',
                           right_string=f'{var_str[v]} Bias',
                           font_height=subtitle_font_height,
                           right_sub_string=f'{yr1}-{yr2}')

         # gbl_mean_trend = ( (trend*area).sum() / area.sum() ).values 
         # hs.set_subtitles( wks, plot[ip2], 
         #                   left_string=name[c],
         #                   center_string=f'Mean: {gbl_mean_trend:5.2f} {unit_str}/yr', 
         #                   right_string=f'{var_str[v]} Trend',
         #                   font_height=subtitle_font_height,
         #                   right_sub_string=f'{yr1}-{yr2}')

#---------------------------------------------------------------------------------------------------
# Finalize plot
pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent       = 5
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.008

layout = [num_var*1,num_case] if var_x_case else [num_case,num_var*1]

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
