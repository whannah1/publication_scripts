import os, copy, string, ngl, xarray as xr, numpy as np, warnings, glob
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
#---------------------------------------------------------------------------------------------------
case_name,case,case_dir,case_sub,case_grid,clr,dsh,mrk = [],[],[],[],[],[],[],[]
def add_case(case_in,n=None,p=None,s=None,g=None,d=0,c='black',m=0):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   if n is None:
      tmp_name = ''
   else:
      tmp_name = n
   case.append(case_in); case_name.append(tmp_name)
   case_dir.append(p); case_sub.append(s); case_grid.append(g)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
#---------------------------------------------------------------------------------------------------
var,lev_list,mask_flag,var_str = [],[],[],[]
comp_list = []
def add_var(var_name,lev=-1,mask='glb',vstr=None,comp='atm'): 
   if vstr is None: vstr = var_name
   var.append(var_name)
   lev_list.append(lev)
   mask_flag.append(mask)
   var_str.append(vstr)
   comp_list.append(comp)
#---------------------------------------------------------------------------------------------------
tmp_path_co2_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_v2  = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
tmp_sub = 'archive/atm/hist'

add_case('v2.LR.historical_0101',                                  n='E3SMv2',  c='cyan',     p=tmp_path_hst_v2, s=tmp_sub)
# add_case('v2.LR.historical_0151',                                  n='E3SMv2',  c='orange',   p=tmp_path_hst_v2, s=tmp_sub)
# add_case('v2.LR.historical_0201',                                  n='E3SMv2',  c='green',    p=tmp_path_hst_v2, s=tmp_sub)
# add_case('v2.LR.historical_0251',                                  n='E3SMv2',  c='purple',   p=tmp_path_hst_v2, s=tmp_sub)
# add_case('v2.LR.historical_0301',                                  n='E3SMv2',  c='pink',     p=tmp_path_hst_v2, s=tmp_sub)
add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',c='blue',p=tmp_path_hst_mmf,s=tmp_sub)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# anomaly_yr1,anomaly_yr2 = 1961,1990 # match HadCRU
anomaly_yr1,anomaly_yr2 = 1950,1980

fig_file,fig_type = 'figs/FXX-timeseries-mpas','png'
add_var('meridionalHeatTransportLat',comp='ocn')

# fig_file,fig_type = 'figs/F01-timeseries-global-hist-chk','png'
# add_var('oceanHeatContentSfcToBot', comp='ocn')
# add_var('oceanHeatContentSfcTo700m', comp='ocn')
# add_var('oceanHeatContent700mTo2000m', comp='ocn')
# add_var('oceanHeatContent2000mToBot', comp='ocn')
# num_plot_col = 2 

#---------------------------------------------------------------------------------------------------
'''
   double oceanHeatContentSfcToBot(time, nCells) ;
   double oceanHeatContentSfcTo700m(time, nCells) ;
   double oceanHeatContent700mTo2000m(time, nCells) ;
   double oceanHeatContent2000mToBot(time, nCells) ;
   double binBoundaryMerHeatTrans(time, nMerHeatTransBinsP1) ;
   double meridionalHeatTransportLatZ(time, nMerHeatTransBinsP1, nVertLevels) ;
   double meridionalHeatTransportLat(time, nMerHeatTransBinsP1) ;
   double refZMid(time, nVertLevels) ;
   double refBottomDepth(time, nVertLevels) ;
'''
#---------------------------------------------------------------------------------------------------

# htype,years,months,first_file,num_files = 'ha',[],[],0,65 # pre-calculated annual means
htype,yr1,yr2 = 'ha',1950,2014

# fig_file,fig_type = 'figs/F01-timeseries-global-hist','png'
tmp_file_head = 'data/timeseries-global-hist'

#---------------------------------------------------------------------------------------------------
write_file    = False
print_stats   = True

recalculate = True

add_trend = False

if 'num_plot_col' not in locals(): num_plot_col  = 3

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)

plot = [None]*(num_var)

wkres = ngl.Resources()
npix=1024*4; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

res = hs.res_xy()
res.vpHeightF = 0.4
# res.vpHeightF = 0.2
res.tmYLLabelFontHeightF   = 0.022 # 0.015
res.tmXBLabelFontHeightF   = 0.022 # 0.015
res.tiXAxisFontHeightF     = 0.022 # 0.015
res.tiYAxisFontHeightF     = 0.022 # 0.015
res.xyLineThicknessF       = 20
# res.tiYAxisString          = 'Temperature Anomaly [K]'
res.tiXAxisString          = 'Time [years]'

# res.tmYLPrecision = 4

lres = hs.res_xy()
lres.xyDashPattern    = 1
lres.xyLineThicknessF = 1
lres.xyLineColor      = "black"

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
land_frac_ds = xr.open_dataset('/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4pg2_16xdel2.c20200108.nc')
land_frac = land_frac_ds['LANDFRAC']
land_frac = land_frac / land_frac.max()
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)

   if 'lev_list' in locals(): lev = lev_list[v]

   time_list,data_list = [],[]
   for c in range(num_case):

      if comp_list[v]=='ocn' and case[c]=='BEST': 
         data_list.append( None )
         time_list.append( None )
         continue

      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}.{mask_flag[v]}.nc'

      print('    case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)
      print('    time series file: '+tmp_file)

      if recalculate:
         #----------------------------------------------------------------------
         scrip_file_path  = os.getenv('HOME')+'/E3SM/data_grid/ne30pg2_scrip.nc'
         scrip_ds = xr.open_mfdataset(scrip_file_path).rename({'grid_size':'ncol'})
         area = scrip_ds['grid_area']
         #----------------------------------------------------------------------
         if comp_list[v]=='ocn':
            if case[c]=='BEST': 
               data = None
            else:
               ocn_scrip_file_path = '/global/cfs/cdirs/e3sm/inputdata/ocn/mpas-o/EC30to60E2r2/EC30to60E2r2.scrip.201005.nc'
               ocn_scrip_ds = xr.open_mfdataset(ocn_scrip_file_path).rename({'grid_size':'ncol'})
               ocn_area = ocn_scrip_ds['grid_area']
               tmp_sub = case_sub[c].replace("atm","ocn")
               data_root = f'{case_dir[c]}/{case[c]}/{tmp_sub}'
               file_path = f'{data_root}/*.mpaso.hist.annual.*'
               file_list = sorted(glob.glob(file_path))
               ds = xr.open_mfdataset( file_list ).rename({'nCells':'ncol'})
               if var[v]=='meridionalHeatTransportLat':
                  ds['binBoundaryMerHeatTrans'] = ds['binBoundaryMerHeatTrans']*180./np.pi
                  latbins = ds['binBoundaryMerHeatTrans'].isel(time=0).values
                  bval = None
                  ctr_lat = None
                  for i,lat in enumerate(latbins):
                     if lat>26.0 and lat<26.9:
                        bval = i
                        break
                  ctr_lat = ( ( latbins[bval] + latbins[bval+1] ) / 2. )
                  # print(f'bval: {bval}')
                  # print(f'bot_lat : {latbins[bval]}')
                  # print(f'ctr_lat : {ctr_lat}')
                  data = ds['meridionalHeatTransportLat'].isel(nMerHeatTransBinsP1=bval)
                  # print()
                  # print(data)
                  # print()
                  # exit()
               else:
                  data = ( (ds[var[v]]*ocn_area).sum(dim='ncol') / ocn_area.sum(dim='ncol') )

         if comp_list[v]=='atm':
            if case[c]=='BEST':
               obs_file = '/global/cfs/cdirs/m3312/whannah/obs_data/BEST/Land_and_Ocean_LatLong1.remap_ne30pg2.nc'
               ds = xr.open_dataset(obs_file)
               data = ds['temperature']
               file_time_list = [None]*len(data['time'])
               # convert anomalies to absolute temperature and create time coord
               for t,tt in enumerate(data['time'].values):
                  yr_val = int(np.floor(tt))
                  mn_ind = int(np.floor((tt-yr_val)*12))
                  file_time_list[t] = np.datetime64(f'{yr_val}-{(mn_ind+1):02d}')
                  # climatology = ds['climatology'].isel(month_number=mn_ind)
                  # climatology = climatology.where(np.isfinite(data[t,:]),np.nan)
                  # data[t,:] = data[t,:] + climatology
               with warnings.catch_warnings():
                  warnings.simplefilter("ignore", category=UserWarning)
                  time = xr.DataArray(file_time_list,dims=('time'))
               data = data.assign_coords(time=time)
               #----------------------------------------------------------------
               # convert Celsius to Kelvin
               data = data + 273.15
               #----------------------------------------------------------------
               # mask area to deal with missing data
               area = area.expand_dims(dim={'time':len(data['time'])}, axis=0)
               area = area.where(np.isfinite(data),np.nan)
               #----------------------------------------------------------------
               # land/ocn mask and area weighted spatial averaging
               if mask_flag[v] in ['lnd','ocn']:
                  mask_data = xr.DataArray( np.ones(land_frac.shape,dtype=bool), dims=land_frac.dims )
                  if mask_flag[v]=='lnd': mask_data = mask_data & (land_frac.values>0.5)
                  if mask_flag[v]=='ocn': mask_data = mask_data & (land_frac.values<0.5)
                  data = data.where( mask_data, drop=True)
                  area = area.where( mask_data, drop=True)
               #----------------------------------------------------------------
               # area weighted global mean
               data = ( (data*area).sum(dim='ncol',skipna=True) / area.sum(dim='ncol',skipna=True) )
               #----------------------------------------------------------------
               # convert monthly to yearly mean
               month_length = data['time'].dt.days_in_month
               mn_wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
               data = (data*mn_wgts).resample(time='A').sum('time') / (mn_wgts).resample(time='A').sum(dim='time')
               #----------------------------------------------------------------
            else:
               file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/*.eam.{htype}.*'
               file_list = sorted(glob.glob(file_path))
               ds = xr.open_mfdataset( file_list )
               data = ds['TS']
               #----------------------------------------------------------------
               # land/ocn mask and area weighted spatial averaging
               if mask_flag[v] in ['lnd','ocn']:
                  mask_data = xr.DataArray( np.ones(land_frac.shape,dtype=bool), dims=land_frac.dims )
                  if mask_flag[v]=='lnd': mask_data = mask_data & (land_frac.values>0.5)
                  if mask_flag[v]=='ocn': mask_data = mask_data & (land_frac.values<0.5)
                  data = data.where( mask_data, drop=True)
                  area = area.where( mask_data, drop=True)
               #----------------------------------------------------------------
               # area weighted global mean
               data = ( (data*area).sum(dim='ncol',skipna=True) / area.sum(dim='ncol',skipna=True) )
         #----------------------------------------------------------------------
         # if data is not None
         #-------------------------------------------------------------------
         # subset in time
         data = data.where( data['time.year']>=yr1, drop=True)
         data = data.where( data['time.year']<=yr2, drop=True)
         time = data['time']
         #-------------------------------------------------------------------
         # write to file
         tmp_ds = xr.Dataset( coords=data.coords )
         tmp_ds[var[v]] = data
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      #-------------------------------------------------------------------------
      else:
         tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
         data = tmp_ds[var[v]]
      #-------------------------------------------------------------------------
      # redefine values as anomalies

      # if comp_list[v]=='atm':
      #    mean_for_anomalies = data.copy(deep=True)
      #    mean_for_anomalies = data.where( data['time.year']>=anomaly_yr1, drop=True)
      #    mean_for_anomalies = data.where( data['time.year']<=anomaly_yr2, drop=True)
      #    mean_for_anomalies = mean_for_anomalies.mean()
      #    data = data - mean_for_anomalies
         
      #-------------------------------------------------------------------------
      if var[v]=='oceanHeatContentSfcTo700m':
         # fac = 1e-22
         fac = 1e-20
         data = data*fac
      #-------------------------------------------------------------------------
      time_mean = data.mean(dim='time').values
      # print('      Area Weighted Time Mean : '+hc.tcolor.GREEN+f'{time_mean:10.6f}'+hc.tcolor.ENDC)
      print('      Area Weighted Time Mean : '+hc.tcolor.GREEN+f'{time_mean}'+hc.tcolor.ENDC)

      if print_stats: hc.print_stat(data,name='',stat='naxs',indent=' '*6,compact=True)

      data_list.append( data.values )
      # time_list.append( data['time'].values )
      time_list.append( data['time.year'].values )

   #----------------------------------------------------------------------------
   # set plot bounds here before loading obs
   
   tmp_time_list = []
   tmp_data_list = []
   for c in range(num_case):
      if time_list[c] is not None: tmp_time_list.append(time_list[c])
      if data_list[c] is not None: tmp_data_list.append(data_list[c])
   
   res.trXMinF = np.min([np.nanmin(d) for d in tmp_time_list])
   res.trXMaxF = np.max([np.nanmax(d) for d in tmp_time_list])
   data_min    = np.min([np.nanmin(d) for d in tmp_data_list])
   data_max    = np.max([np.nanmax(d) for d in tmp_data_list])
   res.trYMinF = data_min - (data_max-data_min)*0.04
   res.trYMaxF = data_max + (data_max-data_min)*0.04

   res.trXMinF = yr1

   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)

   tres.xyCurveDrawOrder = 'PostDraw'

   for c in range(num_case):

      if comp_list[v]=='ocn' and case[c]=='BEST': continue

      ip = v

      tres.tiYAxisString = 'Temperature Anomaly [K]'
      if var[v]=='oceanHeatContentSfcTo700m':
         tres.tiYAxisString = 'Ocean Heat Content [J x 10~S~20~N~]'
      if var[v]=='meridionalHeatTransportLat':
         tres.tiYAxisString = 'Meridional Heat Transport []'
      
      tres.xyLineColor   = clr[c]
      tres.xyDashPattern = dsh[c]

      tplot = ngl.xy(wks, time_list[c], data_list[c], tres)
      
      if plot[ip] is None: 
         plot[ip] = tplot
      else:
         ngl.overlay(plot[ip],tplot)

   #----------------------------------------------------------------------------
   # Set strings at top of plot
   #----------------------------------------------------------------------------
   hs.set_subtitles(wks, plot[ip], '', var_str[v], '', font_height=0.010)#0.008)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
hc.printline()

if 'num_plot_col' in locals():
   layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
else:
   layout = [num_var,num_case]

pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent       = 5
pnl_res.nglPanelXWhiteSpacePercent       = 5
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.012
ngl.panel(wks,plot[0:len(plot)],layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
