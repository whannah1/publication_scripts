# v2 is a cleaned up version of v1 - functionality is mostly the same
import os, glob, ngl, subprocess as sp, numpy as np, xarray as xr, dask, copy, string, cmocean
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
host = hc.get_host()
#-------------------------------------------------------------------------------
name,case,case_dir,case_sub = [],[],[],[]
clr,dsh = [],[]
def add_case(case_in,n=None,p=None,s=None,g=None,c='black',d=0):
   global name,case,case_dir,case_sub
   tmp_name = case_in if n is None else n
   case.append(case_in); name.append(tmp_name); 
   case_dir.append(p); case_sub.append(s); 
   clr.append(c); dsh.append(d)
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
add_case('BEST',n='BEST',c='black')
add_case('v2.LR.historical_0101',                                  n='E3SMv2',  c='red',    p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0101',                                  n='E3SMv2',  c='cyan',    p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0151',                                  n='E3SMv2',  c='orange', p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0201',                                  n='E3SMv2',  c='green',  p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0251',                                  n='E3SMv2',  c='purple',   p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0301',                                  n='E3SMv2',  c='pink', p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical',                                       n='E3SMv2 ens',  c='red',  p=None, s=None)
add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',c='blue', p=tmp_path_hst_mmf,s='archive/atm/hist')
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2',n='1xCO2',c='blue' ,p=tmp_path_co2_mmf,s=tmp_sub)
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.2xCO2',n='2xCO2',c='green',p=tmp_path_co2_mmf,s=tmp_sub)
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.4xCO2',n='4xCO2',c='red'  ,p=tmp_path_co2_mmf,s=tmp_sub)
#-------------------------------------------------------------------------------

fig_file,fig_type = 'figs/FXX-zonal-mean-Tsfc','png'

# add_var('PRECT',name='Precipitation')
# add_var('TS',name='Tsfc')
add_var('TS',name='Tsfc Trend (1979-2014)')
# add_var('TMQ')
# add_var('LHFLX')
# add_var('SHFLX')
# add_var('TGCLDLWP')
# add_var('TGCLDIWP')
# add_var('NET_TOA_RAD')
# add_var('FSNT'); add_var('FLNT')
# add_var('FSNS'); add_var('FLNS')
# add_var('LWCF'); add_var('SWCF')

plot_diff = False

bin_dlat = 2

#-------------------------------------------------------------------------------
# lat1,lat2 = -60,60

# htype,first_file,num_files = 'ha',100,20
htype,yr1,yr2 = 'ha',1979,2014

use_remap,remap_grid = False,'90x180'

print_stats          = True
var_x_case           = False
num_plot_col         = 1 # len(case)
use_common_label_bar = False

#---------------------------------------------------------------------------------------------------
# Set up plot resources
if case==[]: raise ValueError('ERROR - case list is empty!')
num_var,num_case = len(var),len(case)

subtitle_font_height = 0.015


if 'scrip_file_path' not in locals(): scrip_file_path = None

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix

wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*(num_var)
   
res = hs.res_xy()
res.vpHeightF = 0.3
res.xyLineThicknessF = 10
res.tiYAxisString = '[deg C]'
res.tiXAxisString = 'Latitude'
# res.tiXAxisString = 'sin( Latitude )'
res.tmYLAutoPrecision = False
res.tmYLPrecision = 2

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   data_list = []
   bin_list = []
   glb_avg_list = []
   lat_list,lon_list = [],[]
   if 'lev_list' in locals(): lev = lev_list[v]
   #----------------------------------------------------------------------------
   scrip_ds = xr.open_dataset('grid_files/ne30pg2_scrip.nc')
   area = scrip_ds['grid_area'].rename({'grid_size':'ncol'})
   lat  = scrip_ds['grid_center_lat'].rename({'grid_size':'ncol'})
   #----------------------------------------------------------------------------
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)
      #-------------------------------------------------------------------------
      # data_dir_tmp,data_sub_tmp = None, None
      # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
      # if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      # if case_sub[c] is not None: data_sub_tmp = case_sub[c]
      # case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp )
      # case_obj.set_coord_names(var[v])
      # #-------------------------------------------------------------------------
      # # read the data
      # with dask.config.set(**{'array.slicing.split_large_chunks': True}):
      #    tmp_first_file = first_file
      #    if 'v2.LR.historical' in case[c]: tmp_first_file = 50 + first_file
      #    lat  = case_obj.load_data('lat', htype=htype,num_files=1)
      #    area = case_obj.load_data('area',htype=htype,num_files=1).astype(np.double)
      #    data = case_obj.load_data(var[v],htype=htype,ps_htype=htype,lev=lev,
      #                                    first_file=tmp_first_file,num_files=num_files,
      #                                    use_remap=use_remap,remap_str=f'remap_{remap_grid}')
      # #-------------------------------------------------------------------------
      # # Special handling of various specific circumstances
      # if 'lev' in data.dims : data = data.isel(lev=0)
      # if var[v]=='TS': data = data-273

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
         data = data.resample(time='YE').mean(dim='time')
         time = time.resample(time='YE').mean(dim='time')
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
      # # calculate std error and t-statistic for trend
      # dof = yr2 - yr1 - 2
      # y_r = (a*px+b) - py
      # x_r = px - px.mean()
      # yr_variance = (np.square(y_r)).sum(dim='time')
      # xr_variance = (np.square(x_r)).sum(dim='time')
      # stderr_sq = (1/dof) * yr_variance / xr_variance
      # stderr = np.sqrt( stderr_sq )
      # trend_tval = a / stderr

      #-------------------------------------------------------------------------
      # print stats before time averaging
      if print_stats: hc.print_stat(data,name=var[v],stat='naxsh',indent='    ',compact=True)
      #-------------------------------------------------------------------------
      # average over time dimension
      if 'time' in data.dims : 
         hc.print_time_length(data.time,indent=' '*6)
         data = data.mean(dim='time')
      #-------------------------------------------------------------------------
      # Calculate area weighted global mean
      if 'area' in locals() :
         gbl_mean = ( (data*area).sum() / area.sum() ).values 
         print(hc.tcolor.CYAN+f'      Area Weighted Global Mean : {gbl_mean:6.4}'+hc.tcolor.ENDC)
         # glb_avg_list.append(gbl_mean)
      #-------------------------------------------------------------------------
      # Calculate time and zonal mean
      # bin_ds = hc.bin_YbyX( data, lat, 
      bin_ds = hc.bin_YbyX( a, lat, 
                           bin_min=-90, bin_max=90, 
                           bin_spc=bin_dlat, wgt=area )

      data_list.append( bin_ds['bin_val'].values )
      bin_list.append(bin_ds['bins'].values)
      
      # std_list.append( bin_ds['bin_std'].values )
      # cnt_list.append( bin_ds['bin_cnt'].values )

      # lat_bins = bin_ds['bins'].values
      # sin_lat_bins = np.sin(lat_bins*np.pi/180.)
      # bin_list.append(sin_lat_bins)
      #-------------------------------------------------------------------------
      # # append to data lists
      # if case[c]=='TRMM' and 'lon1' not in locals(): data = ngl.add_cyclic(data.values)
      
      # data_list.append( data.values )
   #----------------------------------------------------------------------------
   # calculate common limits for consistent contour levels
   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])
   #----------------------------------------------------------------------------
   # Plot averaged data
   
   res.trXMinF = np.min([np.nanmin(d) for d in  bin_list])
   res.trXMaxF = np.max([np.nanmax(d) for d in  bin_list])
   res.trYMinF = np.min([np.nanmin(d) for d in data_list])
   res.trYMaxF = np.max([np.nanmax(d) for d in data_list])

   lat_tick = np.array([-90,-60,-30,0,30,60,90])
   res.tmXBMode = "Explicit"
   # res.tmXBValues = np.sin( lat_tick*3.14159/180. )
   # res.trXMinF,res.trXMaxF = -1.,1. 
   res.tmXBLabels = lat_tick

   # plot[v] = ngl.xy(wks, np.stack(bin_list), np.ma.masked_invalid(  np.stack(data_list) ), res)

   for c in range(num_case):
      res.xyLineColor   = clr[c]
      res.xyDashPattern = dsh[c]
      tplot = ngl.xy(wks, bin_list[c], np.ma.masked_invalid( data_list[c] ), res)
      if c==0:
         plot[v] = tplot
      else:
         ngl.overlay( plot[v], tplot )

   hs.set_subtitles(wks, plot[v], '', '', 'Tsfc', font_height=subtitle_font_height)
      #-------------------------------------------------------------------------
      # Create plot
      # if use_remap \
      # or case_obj.obs \
      # or ('CESM' in case[c] and 'ne30' not in case[c]) :
      #    if case[c]=='ERAi': lat,lon = erai_lat,erai_lon
      #    hs.set_cell_fill(tres,case_obj=case_obj,lat=lat,lon=lon)
      # else:
      #    if case[c]=='v2.LR.amip_0101' : case_obj.grid = 'ne30pg2'
      #    hs.set_cell_fill(tres,case_obj=case_obj,htype=htype,scrip_file_path=scrip_file_path)

      # tres.lbLabelBarOn = False if use_common_label_bar else True      
         
      # # if plot_diff and c==diff_base : base_name = name[c]

      # num_case_alt = num_case*2-1 if (plot_diff and add_diff) else num_case
      # ip = v*num_case_alt+c if var_x_case else c*num_var+v

      # if not plot_diff  or (plot_diff and add_diff) or (plot_diff and c==diff_base) : 

      #    plot[ip] = ngl.contour_map(wks,np.ma.masked_invalid(data_list[c]),tres) 
         
      #    # set plot subtitles
      #    ctr_str = f'{glb_avg_list[c]:6.4}' if glb_avg_list!=[] else ''
      #    hs.set_subtitles(wks, plot[ip], name[c], ctr_str, var_str[v], font_height=subtitle_font_height)

      
#---------------------------------------------------------------------------------------------------
# Finalize plot

num_case_alt = num_case*2-1 if (plot_diff and add_diff) else num_case
layout = [num_var,num_case_alt] if var_x_case else [num_case_alt,num_var]


if not (plot_diff and add_diff):
   if num_case==1 or num_var==1:
      layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
   
pnl_res = hs.setres_panel()
if use_common_label_bar: pnl_res.nglPanelLabelBar = True
### add panel labels
# pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
# pnl_res.nglPanelFigureStringsJust        = "TopLeft"
# pnl_res.nglPanelFigureStringsFontHeightF = 0.01
# if layout==[3,2] : pnl_res.nglPanelFigureStringsFontHeightF = 0.015

pnl_res.nglPanelYWhiteSpacePercent = 5

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
