import os, copy, ngl, xarray as xr, numpy as np, warnings
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
def add_var(var_name,lev=-1,mask=None,vstr=None): 
   var.append(var_name); lev_list.append(lev),mask_flag.append(mask)
   if vstr is None: vstr = var_name
   var_str.append(vstr)
#---------------------------------------------------------------------------------------------------
tmp_path_co2_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_v2  = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
tmp_sub = 'archive/atm/hist'
add_case('v2.LR.historical',                                       n='E3SMv2 Ens Mean',  c='red', p=tmp_path_hst_v2, s=tmp_sub)
add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',c='blue',p=tmp_path_hst_mmf,s=tmp_sub)

ens_spread_clr = 'lightpink'
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
anomaly_yr1,anomaly_yr2 = 1961,1990 # match HadCRU

add_var('TS',vstr=f'Global Surface Temperature Anomaly from {anomaly_yr1}-{anomaly_yr2} Mean')

htype,years,months,first_file,num_files = 'ha',[],[],0,65 # pre-calculated annual means

# lat1,lat2 = -30,30

fig_file,fig_type = 'figs/F01-timeseries-global-hist','png'
tmp_file_head = 'data/timeseries-global-hist'

#---------------------------------------------------------------------------------------------------
write_file    = False
print_stats   = True
overlay_cases = True

recalculate = False

add_obs_TS = True

add_trend = False

num_plot_col  = 2

# year_start = 1950 + first_file/12
# year_start = 00 + first_file/12

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
num_var = len(var)
num_case = len(case)

obs_pos = num_case

# # overlay_cases with single case causes segfault
# if num_case==1 : overlay_cases = False

if 'clr' not in vars() or clr==[]: clr = ['black']*num_case
if 'dsh' not in vars() or dsh==[]: dsh = np.zeros(num_case)

if 'lev' not in vars(): lev = np.array([0])

if add_obs_TS:
   clr.insert(obs_pos,'black')
   dsh.insert(obs_pos,0)

#-------------------------------------------------------------------------------
# plot legend in separate file
#-------------------------------------------------------------------------------
# if num_case>1:
#    legend_file = fig_file+'.legend'
#    wkres = ngl.Resources() #; npix = 1024 ; wkres.wkWidth,wkres.wkHeight=npix,npix
#    lgd_wks = ngl.open_wks('png',legend_file,wkres)
#    lgres = ngl.Resources()
#    lgres.vpWidthF           = 0.06
#    lgres.vpHeightF          = 0.06#*num_case
#    lgres.lgLabelFontHeightF = 0.008
#    lgres.lgLabelFont        = "courier"
#    lgres.lgMonoDashIndex    = False
#    lgres.lgLineLabelsOn     = False
#    lgres.lgLineThicknessF   = 8#16
#    lgres.lgLabelJust        = 'CenterLeft'
#    lgres.lgLineColors       = clr
#    lgres.lgDashIndexes      = dsh

#    indent = ' '*4
#    labels = case_name
#    for i in range(len(labels)): labels[i] = indent+labels[i] 

#    if add_obs_TS: labels.insert(0,indent+'HadCRU')

#    pid = ngl.legend_ndc(lgd_wks, len(labels), labels, 0.5, 0.65, lgres)

#    ngl.frame(lgd_wks)
#    hc.trim_png(legend_file)

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
if overlay_cases:
   plot = [None]*(num_var)
else:
   plot = [None]*(num_case*num_var)

wkres = ngl.Resources()
npix=1024*4; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

if 'legend_file' in locals(): lgd_wks = ngl.open_wks('png',legend_file,wkres)


res = hs.res_xy()
res.vpHeightF = 0.4
# res.vpHeightF = 0.2
res.tmYLLabelFontHeightF   = 0.015
res.tmXBLabelFontHeightF   = 0.015
res.tiXAxisFontHeightF     = 0.015
res.tiYAxisFontHeightF     = 0.015
res.xyLineThicknessF       = 20
res.tiYAxisString          = 'Temperature Anomaly [K]'
res.tiXAxisString          = 'Time [years]'

lres = hs.res_xy()
lres.xyDashPattern    = 1
lres.xyLineThicknessF = 1
lres.xyLineColor      = "black"

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)

   if 'lev_list' in locals(): lev = lev_list[v]

   time_list,data_list = [],[]
   for c in range(num_case):

      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}.nc'

      print('    case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)
      print('    time series file: '+tmp_file)

      if recalculate:

         #----------------------------------------------------------------------
         # set up the case object
         data_dir_tmp,data_sub_tmp = None, None
         # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
         if case_dir[c] is not None: data_dir_tmp = case_dir[c]
         if case_sub[c] is not None: data_sub_tmp = case_sub[c]

         case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp  )

         tvar = var[v]

         if 'lat1' in vars() : case_obj.lat1,case_obj.lat2 = lat1,lat2
         if 'lon1' in vars() : case_obj.lon1,case_obj.lon2 = lon1,lon2

         #----------------------------------------------------------------------
         # read the data
         tmp_first_file = first_file
         if 'v2.LR.historical' in case[c]: tmp_first_file = 50*12 + first_file

         lat = case_obj.load_data('lat',  htype=htype)
         lon = case_obj.load_data('lon',  htype=htype)
         area = case_obj.load_data('area',htype=htype,num_files=1).astype(np.double)
         data = case_obj.load_data(tvar,  htype=htype, \
                                   first_file=tmp_first_file,\
                                   num_files=num_files,lev=lev)

         time_bnds = case_obj.load_data('time_bnds',htype=htype,first_file=tmp_first_file,num_files=num_files)
         time = time_bnds.isel(nbnd=0)

         # print(); print(time)

         #----------------------------------------------------------------------
         # deal with land mask and area weighted spatial averaging

         land_frac = case_obj.load_data('LANDFRAC',htype='h0',first_file=0,num_files=1).astype(np.double)
         if 'time' in land_frac.dims: land_frac = land_frac.isel(time=0)
         land_frac = land_frac / land_frac.max()

         mask = xr.DataArray( np.ones(land_frac.shape,dtype=bool), dims=land_frac.dims )
         if mask_flag[v]=='lnd': mask = mask & (land_frac.values>0.5)
         if mask_flag[v]=='ocn': mask = mask & (land_frac.values<0.5)

         data = data.where( mask, drop=True)
         data = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )

         # if print_stats: hc.print_stat(data,name='',stat='naxs',indent=' '*6,compact=True)
    
         #----------------------------------------------------------------------
         # Get rid of lev dimension
         if np.all( lev < 0 ) and 'lev' in data.coords : print(f'    lev value: {data.lev.values}')
         if 'lev' in data.dims : data = data.isel(lev=0)
         #----------------------------------------------------------------------
         tmp_ds = xr.Dataset( coords=data.coords )
         tmp_ds[var[v]] = data
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      #-------------------------------------------------------------------------
      else:
         if case[c]=='v2.LR.historical':
            v2_amip_ens_list = []
            v2_amip_ens_list.append('v2.LR.historical_0101')
            v2_amip_ens_list.append('v2.LR.historical_0151')
            v2_amip_ens_list.append('v2.LR.historical_0201')
            v2_amip_ens_list.append('v2.LR.historical_0251')
            v2_amip_ens_list.append('v2.LR.historical_0301')
            cnt = 0
            for e,ens_member in enumerate(v2_amip_ens_list):
               tmp_file = f'{tmp_file_head}.{ens_member}.{var[v]}.nc'
               tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
               ens_member_data = tmp_ds[var[v]]
               if cnt==0: 
                  data = xr.zeros_like(ens_member_data)
                  ens_min = ens_member_data.copy()
                  ens_max = ens_member_data.copy()
               else:
                  for t in range(len(ens_member_data)):
                     ens_min[t] = np.min([ens_min[t].values,ens_member_data[t].values])
                     ens_max[t] = np.max([ens_max[t].values,ens_member_data[t].values])
               data = ( data*cnt + ens_member_data ) / (cnt+1)
               cnt += 1
         else:
            tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
            data = tmp_ds[var[v]]

      #-------------------------------------------------------------------------
      # Make time start at zero
      # year_start = 1950 + first_file/12
      year_start = data['time.year'].values[0]

      data['time'] = ( data['time'] - data['time'][0] ).astype('float') / 86400e9 / 365 + year_start
      
      #-------------------------------------------------------------------------
      # redefine values as anomalies
      if var[v]=='TS' and add_obs_TS:
         # use these years to redefine simulation data as anomalies relative to a climatology
         clim_num_years = anomaly_yr2 - anomaly_yr1 + 1
         yr = data['time'].values
         for t,y in enumerate(yr):
            if y>=year_start: t_start=t;break

         mean_for_anomalies = data.isel(time=slice(t_start,t_start+clim_num_years)).mean()
         data = data - mean_for_anomalies
         
         if case[c]=='v2.LR.historical':
            ens_min = ens_min - mean_for_anomalies
            ens_max = ens_max - mean_for_anomalies

      #-------------------------------------------------------------------------
      time_mean = data.mean(dim='time').values
      # print('      Area Weighted Time Mean : '+hc.tcolor.GREEN+f'{time_mean:10.6f}'+hc.tcolor.ENDC)
      print('      Area Weighted Time Mean : '+hc.tcolor.GREEN+f'{time_mean}'+hc.tcolor.ENDC)

      if print_stats: hc.print_stat(data,name='',stat='naxs',indent=' '*6,compact=True)

      data_list.append( data.values )
      time_list.append( data['time'].values )

      #-------------------------------------------------------------------------
      # write to file
      #-------------------------------------------------------------------------
      # if write_file : 
      #    tfile = f'/global/homes/w/whannah/E3SM/scratch/{case[0]}/run/{case[0]}.time_series.{var[v]}.nc'
      #    print('writing to file: '+tfile)
      #    avg_X.name = var[v]
      #    avg_X.to_netcdf(path=tfile,mode='w')
      #    exit()

   #----------------------------------------------------------------------------
   # set plot bounds here before loading obs
   res.trXMinF = np.min([np.nanmin(d) for d in time_list])
   res.trXMaxF = np.max([np.nanmax(d) for d in time_list])
   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])
   res.trYMinF = data_min - (data_max-data_min)*0.04
   res.trYMaxF = data_max + (data_max-data_min)*0.04

   #----------------------------------------------------------------------------
   # Add Obs data
   #----------------------------------------------------------------------------
   # if var[v]=='TS' and add_obs_TS:
   #    print(' '*4+'Loading HadCRU data...')
   #    file_obs = '/global/cfs/cdirs/m3312/whannah/obs_data/HadCRU/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc'
   #    # file_obs = '/gpfs/alpine/scratch/hannah6/cli115/Obs/HadCRU/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc'
   #    ds_obs = xr.open_dataset(file_obs)
   #    data_obs = ds_obs['tas_mean']
   #    # hc.print_stat(data_obs,name='HadCRU',stat='naxs',indent=' '*6,compact=True)
   #    xy_dims = ('longitude','latitude')
   #    xlon, ylat = np.meshgrid(ds_obs['latitude'],ds_obs['longitude'])
   #    R = hc.earth_radius(ylat)
   #    dlat = np.deg2rad(np.gradient(ylat, axis=0))
   #    dlon = np.deg2rad(np.gradient(xlon, axis=1))
   #    dy,dx = dlat * R , dlon * R * np.cos(np.deg2rad(ylat))
   #    area_obs = np.absolute(dy*dx) / np.square(R) # calculate area and convert to steridians
   #    area_obs = xr.DataArray(area_obs,dims=xy_dims).transpose()
   #    gbl_mean_obs = ( (data_obs*area_obs).sum(dim=xy_dims) / area_obs.sum(dim=xy_dims) )

   #    # convert to annual mean
   #    time_obs = ds_obs['time_bnds'].isel(bnds=0)
   #    month_length = time_obs.dt.days_in_month
   #    # print(); print(time_obs)
   #    # print(); print(month_length)
   #    month_length['time'] = time_obs
   #    # print(); print(month_length)
   #    gbl_mean_obs['time'] = time_obs
   #    # print(); print(month_length.groupby("time.year"))
   #    # exit()
   #    wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
   #    gbl_mean_obs = (gbl_mean_obs*wgts).resample(time='A').sum('time') / (wgts).resample(time='A').sum(dim='time')

   #    hc.print_stat(gbl_mean_obs,name='HadCRU',stat='naxs',indent=' '*6,compact=True)

   #    time_obs = ds_obs['time.year'].resample(time='Y').mean(dim='time')

   #    # limit extent of obs data to match size of model data
   #    sim_num_t = np.max([len(d) for d in time_list])

   #    for t,y in enumerate(time_obs):
   #       # if y>=0: t_start=t;break
   #       if y>=year_start: t_start=t;break

   #    # print(); print(f't_start: {t_start}'); print()

   #    gbl_mean_obs = gbl_mean_obs.isel(time=slice(t_start,t_start+sim_num_t))
   #    time_obs     = time_obs    .isel(time=slice(t_start,t_start+sim_num_t))

   #    data_list.insert(obs_pos, gbl_mean_obs.values )
   #    time_list.insert(obs_pos, time_obs.values )

   if var[v]=='TS' and add_obs_TS:
      print(' '*4+'Loading BEST data...')
      obs_file = '/global/cfs/cdirs/m3312/whannah/obs_data/BEST/Land_and_Ocean_LatLong1.remap_ne30pg2.nc'
      ds = xr.open_dataset(obs_file)#.rename({'latitude':'lat','longitude':'lon'})
      data = ds['temperature']
      file_time_list = [None]*len(data['time'])
      # convert anomalies to absolute temperature and build time coord
      for t,tt in enumerate(data['time'].values):
         yr_val = int(np.floor(tt))
         mn_ind = int(np.floor((tt-yr_val)*12))
         file_time_list[t] = np.datetime64(f'{yr_val}-{(mn_ind+1):02d}')
         data[t,:] = data[t,:] #+ ds['climatology'].isel(month_number=mn_ind)
      with warnings.catch_warnings():
         warnings.simplefilter("ignore", category=UserWarning)
         time = xr.DataArray(file_time_list,dims=('time'))
      time = time.assign_coords(time=time)
      month_length = time.dt.days_in_month
      wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
      data = data.assign_coords(time=time)
      # convert to annual mean
      data = data.resample(time='Y').mean(dim='time')
      time = time.resample(time='Y').mean(dim='time')
      # calculate global mean
      scrip_ds = xr.open_dataset('grid_files/ne30pg2_scrip.nc')
      area = scrip_ds['grid_area'].rename({'grid_size':'ncol'})

      area = area.expand_dims(dim={'time':len(data['time'])}, axis=0)
      area = area.where(np.isfinite(data),np.nan)
      
      gbl_mean_obs = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )

      # gbl_mean_obs = np.nansum(data*area,axis=1) / np.nansum(area,axis=1)
      # gbl_mean_obs = xr.DataArray(gbl_mean_obs,coords={'time':data['time']})

      for t,y in enumerate(gbl_mean_obs['time.year'].values):
         if y>=year_start: t_start=t;break
      clim_num_years = anomaly_yr2 - anomaly_yr1 + 1
      mean_for_anomalies = gbl_mean_obs.isel(time=slice(t_start,t_start+clim_num_years)).mean()
      gbl_mean_obs = gbl_mean_obs - mean_for_anomalies

      hc.print_stat(gbl_mean_obs,name='BEST',stat='naxs',indent=' '*6,compact=True)

      # time_obs = ds_obs['time.year'].resample(time='Y').mean(dim='time')
      # # limit extent of obs data to match size of model data
      # sim_num_t = np.max([len(d) for d in time_list])
      # for t,y in enumerate(time_obs):
      #    # if y>=0: t_start=t;break
      #    if y>=year_start: t_start=t;break
      # # print(); print(f't_start: {t_start}'); print()
      # gbl_mean_obs = gbl_mean_obs.isel(time=slice(t_start,t_start+sim_num_t))
      # time_obs     = time_obs    .isel(time=slice(t_start,t_start+sim_num_t))

      time = np.array(gbl_mean_obs['time.year'].values,dtype=float)

      data_list.insert(obs_pos, gbl_mean_obs.values )
      time_list.insert(obs_pos, time )

   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)

   tres.xyCurveDrawOrder = 'PostDraw'

   ### Make sure plot bounds are consistent
   # tres.trYMinF = 14.4

   # data_min = np.min([np.nanmin(d) for d in data_list])
   # data_max = np.max([np.nanmax(d) for d in data_list])
   # tres.trYMinF = data_min - (data_max-data_min)*0.04
   # tres.trYMaxF = data_max + (data_max-data_min)*0.04
   # tres.trXMinF = np.min([np.nanmin(d) for d in time_list])
   # tres.trXMaxF = np.max([np.nanmax(d) for d in time_list])

   if var[v]=='NET_TOA_RAD':
      tres.trYMinF = -20
      tres.trYMaxF =  20

   tmp_num_case = num_case
   if var[v]=='TS' and add_obs_TS: tmp_num_case = num_case + 1

   for c in range(tmp_num_case):
      ip = c*num_var + v
      if overlay_cases: ip = v
      
      tres.xyLineColor   = clr[c]
      tres.xyDashPattern = dsh[c]

      # print()
      # print()
      # print()
      # print()
      # print(); print(time_list)
      # # print(); print(data_list[c])
      # print()
      # exit()

      tplot = ngl.xy(wks, time_list[c], data_list[c], tres)


      
      if overlay_cases: 
         if c==0: 
            plot[ip] = tplot
         else:
            ngl.overlay(plot[ip],tplot)
      else:
         plot[ip] = tplot

      #------------------------------------------------
      # add ensemble spread
      #------------------------------------------------
      ens_case = 'v2.LR.historical'
      if (add_obs_TS and c>0 and case[c-1]==ens_case) or (not add_obs_TS and case[c]==ens_case):
         # eres = copy.deepcopy(tres)
         # ens_spread_data      = np.zeros([2,len(ens_min)])
         # ens_spread_data[0,:] = ens_min.values
         # ens_spread_data[1,:] = ens_max.values
         # eres.xyLineColor = [0,0,0,0]
         # eres.nglXYAboveFillColors = 'orange'
         # eres.nglXYBelowFillColors = 'orange'
         # ngl.overlay(plot[ip], ngl.xy(wks, time_list[c], ens_spread_data, eres) )

         # tmp = np.array([0,1,2,3,4])
         # n = len(tmp)
         # print(); print(tmp[:: 1])
         # print(); print(tmp[::-1])

         # n = len(time_list[c])
         # print(); print(time_list[c][:: 1])
         # print(); print(time_list[c][::-1])
         # exit()

         n = len(time_list[c])
         ens_spread_data = np.zeros( 2*n+1 )
         ens_spread_data[n*0:n*1] = ens_min[:: 1].values
         ens_spread_data[n*1:n*2] = ens_max[::-1].values
         ens_spread_data[n*2]       = ens_min[0].values
         ens_spread_time = np.zeros( 2*n+1 )
         ens_spread_time[n*0:n*1] = time_list[c][:: 1]
         ens_spread_time[n*1:n*2] = time_list[c][::-1]
         ens_spread_time[n*2]       = time_list[c][0]
         eres = ngl.Resources()
         eres.gsFillColor = ens_spread_clr
         # eres.gsFillOpacityF = 0.1
         # eres.tfPolyDrawOrder = 'Draw'
         dum = ngl.add_polygon(wks, plot[ip], ens_spread_time, ens_spread_data, eres)

      #------------------------------------------------
      # add linear trend
      #------------------------------------------------
      if add_trend:
         px = time_list[c]
         py = data_list[c]
         # simple and fast method for regression coeff and intercept
         a = np.cov( px.flatten(), py.flatten() )[1,0] / np.var( px )
         b = np.mean(py) - a*np.mean(px)

         # print regression info
         # if c==0: print()
         # print(' '*4+f'linear regression a: {a}    b: {b}')
         # if c==(num_case-1): print()

         px_range = np.abs( np.max(px) - np.min(px) )
         lx = np.array([-1e2*px_range,1e2*px_range])

         lres.xyLineColor = clr[c]
         ngl.overlay( plot[ip], ngl.xy(wks, lx, lx*a+b , lres) )

   #----------------------------------------------------------------------------
   # Set strings at top of plot
   #----------------------------------------------------------------------------
   # lft_str,ctr_str = '',''
   # lat_chk,lon_chk = 'lat1' in locals(), 'lon1' in locals()
   # if not lat_chk and not lon_chk : 
   #    if mask_flag[v] is None: ctr_str = 'Global'
   #    if mask_flag[v]=='lnd' : ctr_str = 'Land Only'
   #    if mask_flag[v]=='ocn' : ctr_str = 'Ocean Only'
   #    if mask_flag[v]=='lnd/ocn ratio':ctr_str = 'Land / Ocean Ratio'
   # else:
   #    if lat_chk:      ctr_str += f' {lat1}:{lat2}N '
   #    if lon_chk:      ctr_str += f' {lon1}:{lon2}E '
   
   if overlay_cases:
      hs.set_subtitles(wks, plot[ip], '', var_str[v], '', font_height=0.01)
   else:
      hs.set_subtitles(wks, plot[ip], case_name[c], '', var_str[v], font_height=0.01)

   #----------------------------------------------------------------------------
   # Ad legend
   #----------------------------------------------------------------------------
   lgres = ngl.Resources()
   lgres.vpWidthF, lgres.vpHeightF  = 0.08, 0.08  
   lgres.lgLabelFontHeightF = 0.01
   lgres.lgLineThicknessF   = 4
   lgres.lgMonoLineColor    = False
   # lgres.lgMonoDashIndex    = True
   lgres.lgLineColors       = clr
   lgres.lgDashIndexes      = dsh
   lgres.lgLabelJust    = 'CenterLeft'
   # pid = ngl.legend_ndc(wks, len(name), name, 0.5, 0.4, lgres)  # 3x2
   # pid = ngl.legend_ndc(wks, len(name), name, 0.5, 0.1, lgres)  # 3x2
   # pid = ngl.legend_ndc(wks, len(name), name, 0.3, 0.5, lgres)  # 1x2

#-------------------------------------------------------------------------------
# Add legend
#-------------------------------------------------------------------------------

lgres = ngl.Resources()
lgres.vpWidthF           = 0.06
lgres.vpHeightF          = 0.08#*num_case
lgres.lgLabelFontHeightF = 0.008
lgres.lgLabelFont        = "courier"
lgres.lgMonoDashIndex    = False
lgres.lgLineLabelsOn     = False
lgres.lgLineThicknessF   = 30
lgres.lgLabelJust        = 'CenterLeft'
lgres.lgLineColors       = clr
lgres.lgDashIndexes      = dsh

indent = ' '*4
labels = case_name
for i in range(len(labels)): labels[i] = indent+labels[i] 

# if add_obs_TS: labels.insert(obs_pos,indent+'HadCRU')
if add_obs_TS: labels.insert(obs_pos,indent+'BEST')
if 'v2.LR.historical' in case:
   labels.insert(0,indent+'E3SMv2 Ens Min/Max')
   lgres.lgLineColors.insert(0,ens_spread_clr)
   lgres.lgDashIndexes.insert(0,0)
   # xyLineOpacities

pid = ngl.legend_ndc(wks, len(labels), labels, 0.38, 0.63, lgres)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
hc.printline()

if 'num_plot_col' in locals():
   if overlay_cases :
      layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
   else:
      if num_case==1 or num_var==1:
         layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
      else:
         layout = [num_var,num_case]
         # layout = [num_case,num_var]
else:
   layout = [num_var,num_case]

ngl.panel(wks,plot[0:len(plot)],layout,hs.setres_panel())
ngl.end()

hc.trim_png(fig_file)

if 'legend_file' in locals(): hc.trim_png(legend_file)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
