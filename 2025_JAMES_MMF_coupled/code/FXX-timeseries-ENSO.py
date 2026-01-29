import os, copy, ngl, xarray as xr, numpy as np
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
var,lev_list,mask_flag = [],[],[]
def add_var(var_name,lev=-1,mask=None): 
   var.append(var_name); lev_list.append(lev),mask_flag.append(mask)
#---------------------------------------------------------------------------------------------------
tmp_path_co2_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_v2  = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
tmp_sub = 'archive/atm/hist'
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2',n='MMF PI 1xCO2',c='blue' ,p=tmp_path_co2_mmf,s=tmp_sub)
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.2xCO2',n='MMF PI 2xCO2',c='green',p=tmp_path_co2_mmf,s=tmp_sub)
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.4xCO2',n='MMF PI 4xCO2',c='red'  ,p=tmp_path_co2_mmf,s=tmp_sub)
# add_case('v2.LR.historical_0101',                                  n='E3SMv2',  p=tmp_path_hst_v2, s=tmp_sub)
# add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',p=tmp_path_hst_mmf,s=tmp_sub)
#---------------------------------------------------------------------------------------------------
# scrip_file_path = os.getenv('HOME')+'/E3SM/data_grid/ne30pg2_scrip.nc'

add_var('TS')

# htype,first_file,num_files = 'h0',    0,35*12
htype,first_file,num_files = 'h0',50*12,70*12
# htype,first_file,num_files = 'h0',50*12,5*12


# region,lat1,lat2,lon1,lon2 = 'Nino3'  ,-5,5,360-150,360-90
# region,lat1,lat2,lon1,lon2 = 'Nino4'  ,-5,5,160,360-150
region,lat1,lat2,lon1,lon2 = 'Nino3.4',-5,5,190,240


fig_type,fig_file = 'png','figs/FXX-timeseries-ENSO'
tmp_file_head = 'data/timeseries-ENSO'

write_file    = False
print_stats   = True
overlay_cases = False

convert_to_annual_mean = False
remove_annual_cycle    = True

recalculate = True
add_obs_TS  = False
add_trend   = False

num_plot_col  = 1

# year_start = 1950+first_file/12
year_start = 50

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
num_var = len(var)
num_case = len(case)

# overlay_cases with single case causes segfault
if num_case==1 and not add_obs_TS : overlay_cases = False

if 'clr' not in vars() or clr==[]: clr = ['black']*num_case
if 'dsh' not in vars() or dsh==[]: dsh = np.zeros(num_case)

if 'lev' not in vars(): lev = np.array([0])

#-------------------------------------------------------------------------------
# plot legend in separate file
#-------------------------------------------------------------------------------
if num_case>1:
   legend_file = fig_file+'.legend'
   wkres = ngl.Resources() #; npix = 1024 ; wkres.wkWidth,wkres.wkHeight=npix,npix
   lgd_wks = ngl.open_wks('png',legend_file,wkres)
   lgres = ngl.Resources()
   lgres.vpWidthF           = 0.05
   lgres.vpHeightF          = 0.03*num_case
   lgres.lgLabelFontHeightF = 0.008
   lgres.lgLabelFont        = "courier"
   lgres.lgMonoDashIndex    = False
   lgres.lgLineLabelsOn     = False
   lgres.lgLineThicknessF   = 2
   lgres.lgLabelJust        = 'CenterLeft'
   lgres.lgLineColors       = clr
   lgres.lgDashIndexes      = dsh

   labels = case_name
   for i in range(len(labels)): labels[i] = ' '*4+labels[i] 

   pid = ngl.legend_ndc(lgd_wks, len(labels), labels, 0.5, 0.65, lgres)

   ngl.frame(lgd_wks)
   hc.trim_png(legend_file)
   # exit()

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
if overlay_cases:
   plot = [None]*(num_var)
else:
   if add_obs_TS:
      plot = [None]*((num_case+1)*num_var)
   else:
      plot = [None]*(num_case*num_var)

wkres = ngl.Resources()
npix=1024; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
if 'legend_file' in locals(): lgd_wks = ngl.open_wks('png',legend_file,wkres)

res = hs.res_xy()
# res.vpHeightF = 0.5
res.vpHeightF = 0.2
res.tmYLLabelFontHeightF         = 0.015
res.tmXBLabelFontHeightF         = 0.015
res.tiXAxisFontHeightF           = 0.015
res.tiYAxisFontHeightF           = 0.015
res.xyLineThicknessF = 6

lres = hs.res_xy()
lres.xyDashPattern    = 1
lres.xyLineThicknessF = 1
lres.xyLineColor      = "black"

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def CESM_FV(case):
   CESM_FV = False
   if 'CESM' in case and any([g in case for g in ['f09','f19']]): CESM_FV = True
   return CESM_FV

def get_comp(case):
   comp = 'eam'
   if case=='ERA5': comp = None
   if case=='MAC': comp = None
   if 'CESM' in case: comp = 'cam'
   if case=='E3SM.RGMA.ne30pg2_r05_oECv3.FC5AV1C-L.00': comp = 'cam'
   if 'E3SM.PI-CPL.v1.' in case: comp = 'cam'
   return comp

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)

   if 'lev_list' in locals(): lev = lev_list[v]

   time_list,data_list = [],[]
   for c in range(num_case):

      print('    case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)

      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}.nc'

      if recalculate:

         data_dir_tmp,data_sub_tmp = None, None
         # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
         if case_dir[c] is not None: data_dir_tmp = case_dir[c]
         if case_sub[c] is not None: data_sub_tmp = case_sub[c]

         case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp  )

         tvar = var[v]

         if 'lat1' in vars() : case_obj.lat1,case_obj.lat2 = lat1,lat2
         if 'lon1' in vars() : case_obj.lon1,case_obj.lon2 = lon1,lon2

         #-------------------------------------------------------------------------
         # read the data
         #-------------------------------------------------------------------------
         tnum_files = num_files
         tfirst_file = first_file
         # if case=='E3SM.INCITE2022-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR':
         #    tnum_files = 12*40
         if case[c]=='v2.LR.historical_0101': tfirst_file = 50*12


         lat  = case_obj.load_data('lat',  component=get_comp(case[c]),htype=htype)
         lon  = case_obj.load_data('lon',  component=get_comp(case[c]),htype=htype)
         area = case_obj.load_data('area',component=get_comp(case[c]),htype=htype,num_files=1).astype(np.double)
         data = case_obj.load_data(tvar,  component=get_comp(case[c]),htype=htype,
                                    first_file=tfirst_file,num_files=tnum_files,lev=lev)

         # if htype in ['h0'] and convert_to_annual_mean:
         #    num_time = len(data.time)
         #    #print(num_time); print(np.floor(num_time/12)*12); exit()
         #    data = data.isel(time=slice(0,int(np.floor(num_time/12)*12)))

         # land or ocean mask
         if mask_flag[v] is not None :
            land_frac = case_obj.load_data('LANDFRAC',component=get_comp(case[c]),htype='h0',
                                             first_file=first_file,num_files=tnum_files).astype(np.double)
            land_frac = land_frac.isel(time=0)
         if mask_flag[v]=='lnd': data = data*land_frac
         if mask_flag[v]=='ocn': data = data*(1-land_frac)

         # if print_stats: hc.print_stat(data,name='',stat='naxs',indent=' '*6,compact=True)

         if htype in ['h0'] and convert_to_annual_mean: 
            # truncate months past the last full year
            extra_months = len(data.time)%12
            print(f'  extra_months: {extra_months}')
            if extra_months !=0: data = data.isel(time=slice(0,-1-extra_months))
            # data = data.resample(time='Y').mean(dim='time')
            month_length = data.time.dt.days_in_month
            wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
            data = (data*wgts).resample(time='A').sum('time') / (wgts).resample(time='A').sum(dim='time')

         if np.all( lev < 0 ) and 'lev' in data.coords : print(f'    lev value: {data.lev.values}')

         #-------------------------------------------------------------------------
         #reset time index to start at zero and convert to days
         dtime = ( data['time'][-1] - data['time'][0] ).values.astype('timedelta64[D]')
         print('      Time length: '+str(dtime)+'  ('+str(dtime.astype('timedelta64[M]'))+')')

         #-------------------------------------------------------------------------
         # time/space means
         avg_X = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )

         if remove_annual_cycle:
            for n in range(12):
               month_index = slice(n,len(avg_X.time),12)
               avg_X[month_index] = avg_X[month_index] - np.mean(avg_X[month_index])
         
         # time_mean = avg_X.mean(dim='time').values
         # print('      Area Weighted Time Mean : '+hc.tcolor.GREEN+f'{time_mean:10.6f}'+hc.tcolor.ENDC)
         
         # Make time start at zero
         avg_X['time'] = ( avg_X['time'] - avg_X['time'][0] ).astype('float') / 86400e9
         
         avg_X = avg_X - avg_X.mean()

         # detrend in time
         fit = xr.polyval(avg_X['time'], avg_X.polyfit(dim='time', deg=1).polyfit_coefficients)
         avg_X = avg_X - fit

         avg_X.load()

         time = avg_X['time']

         #-------------------------------------------------------------------------
         # Write to file 
         #-------------------------------------------------------------------------
         if os.path.isfile(tmp_file) : os.remove(tmp_file)
         tmp_ds = xr.Dataset()
         tmp_ds[var[v]] = avg_X
         tmp_ds['time'] = time
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         tmp_ds = xr.open_dataset( tmp_file )
         avg_X = tmp_ds[var[v]]
         time  = tmp_ds['time']

      # print(); print(time)
      # print(); print(avg_X.values)
      # exit()

      
      #-------------------------------------------------------------------------
      if print_stats: hc.print_stat(avg_X,name='',stat='naxs',indent=' '*6,compact=True)
         
      data_list.append( avg_X.values )
      time_list.append( time.values/365 )

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
   # Add Obs data
   #----------------------------------------------------------------------------
   if var[v]=='TS' and add_obs_TS:
      obs_name = 'HadSST' # HadCRU / HadSST / NOAA

      if obs_name=='HadCRU':
         lat_name,lon_name = 'latitude','longitude'
         print(' '*4+'Loading HadCRU data...')
         file_obs = '/global/cfs/cdirs/m3312/whannah/obs_data/HadCRU/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc'
         ds_obs = xr.open_dataset(file_obs)
         data_obs = ds_obs['tas_mean']

      if obs_name=='HadSST':
         lat_name,lon_name = 'latitude','longitude'
         print(' '*4+'Loading HadSST data...')
         file_obs = os.getenv('HOME')+'/Data/Obs/HadSST/HadISST_sst.nc'
         ds_obs = xr.open_dataset(file_obs)
         data_obs = ds_obs['sst']

      if obs_name=='NOAA':
         lat_name,lon_name = 'lat','lon'
         print(' '*4+'Loading NOAA OI SST data...')
         file_obs = os.getenv('HOME')+'/Data/Obs/NOAA/monthly/sst.mnmean.nc'
         ds_obs = xr.open_dataset(file_obs)
         data_obs = ds_obs['sst'] * ds_obs['scale_factor'] + ds_obs['add_offset']

      xy_dims = (lon_name,lat_name)
      xlon, ylat = np.meshgrid(ds_obs[lat_name],ds_obs[lon_name])
      R = hc.earth_radius(ylat)
      drlat = np.deg2rad(np.gradient(ylat, axis=0))
      drlon = np.deg2rad(np.gradient(xlon, axis=1))

      dy,dx = drlat * R , drlon * R * np.cos(np.deg2rad(ylat))
      area_obs = np.absolute(dy*dx) / np.square(R) # calculate area and convert to steridians
      area_obs = xr.DataArray(area_obs,dims=xy_dims).transpose()

      xlon = np.transpose(xlon)
      ylat = np.transpose(ylat)

      tlon1,tlon2 = lon1,lon2
      if tlon1>180: tlon1 = tlon1 - 360
      if tlon2>180: tlon2 = tlon2 - 360
      lon1,lon2 = tlon1,tlon2

      tmp_data = np.ones([len(ds_obs[lat_name]),len(ds_obs[lon_name])],dtype=bool)
      tmp_coords = {lat_name:ds_obs[lat_name],lon_name:ds_obs[lon_name]}
      mask = xr.DataArray( tmp_data, coords=tmp_coords, dims=(lat_name,lon_name) )
      mask = mask & (ds_obs[lat_name]>=lat1) & (ds_obs[lat_name]<=lat2)
      mask = mask & (ds_obs[lon_name]>=lon1) & (ds_obs[lon_name]<=lon2)
      # mask = mask & (ylat>=lat1) & (ylat<=lat2) & (xlon>=lon1) & (xlon<=lon2)


      # print(); print(ds_obs[lon_name])
      # print(); print(lon1)
      # print(); print(lon2)
      # exit()

      # data_obs.load()
      # area_obs.load()

      data_obs = data_obs.where( mask, drop=True)
      area_obs = area_obs.where( mask, drop=True)

      # print(); print(area_obs)

      # print();hc.print_stat(area_obs,name='area_obs',stat='naxs',indent=' '*6,compact=True)
      # print();hc.print_stat(data_obs,name='data_obs',stat='naxs',indent=' '*6,compact=True)
      # exit()

      gbl_mean_obs = ( (data_obs*area_obs).sum(dim=xy_dims) / area_obs.sum(dim=xy_dims) )

      if remove_annual_cycle:
         for n in range(12):
            month_index = slice(n,len(gbl_mean_obs.time),12)
            gbl_mean_obs[month_index] = gbl_mean_obs[month_index] - np.mean(gbl_mean_obs[month_index])

      if convert_to_annual_mean: 
         gbl_mean_obs = gbl_mean_obs.resample(time='Y').mean(dim='time')
         # month_length = data.time.dt.days_in_month
         # wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
         # data = (data*wgts).resample(time='A').sum('time') / (wgts).resample(time='A').sum(dim='time')

      hc.print_stat(gbl_mean_obs,name=obs_name,stat='naxs',indent=' '*6,compact=True)
      # exit()

      if convert_to_annual_mean:
         # time_obs = ds_obs['time.year'].isel[slice(0,0,365)]
         time_obs = ds_obs['time.year'].resample(time='Y').mean(dim='time')
      else:
         time_obs = ds_obs['time.year'] + ds_obs['time.dayofyear']/365

      time_obs = time_obs - year_start

      # limit extent of obs data to match size of model data
      sim_num_t = np.max([len(d) for d in time_list])

      t_start = None
      for t,y in enumerate(time_obs):
         if y>=0: t_start=t;break

      if t_start is None: raise ValueError('ERROR: no appropriate value found for t_start')

      gbl_mean_obs = gbl_mean_obs.isel(time=slice(t_start,t_start+sim_num_t))
      time_obs     = time_obs    .isel(time=slice(t_start,t_start+sim_num_t))
      case_name.insert(0, obs_name )
      data_list.insert(0, gbl_mean_obs.values )
      time_list.insert(0, time_obs.values )

      clr.insert(0,'black')
      dsh.insert(0,0)

      # print();print(ds_obs)
      # print();print(ds_obs.time)
      # print();print(time_obs)
      # print();print(gbl_mean_obs)

      # for i,y in enumerate(time_obs.values):
      #    if i%12==0:
      #       print(f'  {i:4d}  {y}  {gbl_mean_obs.values[i]}')
      #       # if i==12: exit()

      # tres.xyLineColor   = 'black'
      # tres.xyDashPattern = 0
      # ngl.overlay(plot[ip], ngl.xy(wks, time_obs.values, gbl_mean_obs.values , tres) )

      # replace plot with obs-obly plot for debugging
      # tres = copy.deepcopy(res)
      # tres.trXMinF = np.min(time_obs.values)
      # tres.trXMaxF = np.max(time_obs.values)
      # plot[ip] = ngl.xy(wks, time_obs.values, gbl_mean_obs.values , tres)

   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)
   # tres.tiXAxisString = 'Time [days]'
   tres.tiXAxisString = 'Time [years]'
   if convert_to_annual_mean: tres.tiXAxisString = 'Time [years]'

   # reset start year for all data
   for t in range(len(time_list)):
      time_list[t] = time_list[t] + year_start


   ### Make sure plot bounds are consistent
   tres.trYMinF = np.min([np.nanmin(d) for d in data_list])
   tres.trYMaxF = np.max([np.nanmax(d) for d in data_list])
   tres.trXMinF = np.min([np.nanmin(d) for d in time_list])
   tres.trXMaxF = np.max([np.nanmax(d) for d in time_list])


   tmp_num_case = num_case
   if var[v]=='TS' and add_obs_TS: tmp_num_case = num_case + 1

   for c in range(tmp_num_case):
      ip = c*num_var + v
      if overlay_cases: ip = v
      
      tres.xyLineColor   = clr[c]
      tres.xyDashPattern = dsh[c]
      tplot = ngl.xy(wks, time_list[c], data_list[c], tres)
      
      if overlay_cases: 
         if c==0: 
            plot[ip] = tplot
         else:
            ngl.overlay(plot[ip],tplot)
      else:
         plot[ip] = tplot

      cres = copy.deepcopy(tres)
      fill_data      = np.zeros([2,len(data_list[c])])
      fill_data[0,:] = data_list[c]
      cres.xyLineColor = [0,0,0,0]

      fill_data[1,:] = 0.5
      cres.nglXYAboveFillColors = 'red'
      cres.nglXYBelowFillColors = -1
      ngl.overlay(plot[ip], ngl.xy(wks, time_list[c], fill_data, cres) )
      
      fill_data[1,:] = -0.5
      cres.nglXYAboveFillColors = -1
      cres.nglXYBelowFillColors = 'blue'
      ngl.overlay(plot[ip], ngl.xy(wks, time_list[c], fill_data, cres) )

      xx = np.array([-1e8,1e8])
      yy = np.array([0,0])
      ngl.overlay( plot[ip], ngl.xy(wks, xx, yy, lres) )

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
   
   #------------------------------------------------
   #------------------------------------------------


   # if overlay_cases:
   #    ip = v

   #    ### use this for overlaying variables on same plot
   #    for c in range(num_case):
   #       tres.xyLineColor   = clr[c]
   #       tres.xyDashPattern = dsh[c]
   #       tplot = ngl.xy(wks, time_list[c], data_list[c], tres)
   #       if c==0: 
   #          plot[ip] = tplot
   #       else:
   #          ngl.overlay(plot[ip],tplot)

   #       #------------------------------------------------
   #       # add linear trend
   #       #------------------------------------------------
   #       if add_trend:
   #          px = time_list[c]
   #          py = data_list[c]
   #          # simple and fast method for regression coeff and intercept
   #          a = np.cov( px.flatten(), py.flatten() )[1,0] / np.var( px )
   #          b = np.mean(py) - a*np.mean(px)

   #          print(f'\n    linear regression a: {a}    b: {b}\n')

   #          px_range = np.abs( np.max(px) - np.min(px) )
   #          lx = np.array([-1e2*px_range,1e2*px_range])

   #          lres.xyLineColor = clr[c]
   #          ngl.overlay( plot[ip], ngl.xy(wks, lx, lx*a+b , lres) )
   #       #------------------------------------------------
   #       #------------------------------------------------
   # else:
   #    for c in range(num_case):
   #       ip = c*num_var + v
   #       # ip = v*num_case + c
   #       tres.xyLineColor   = clr[c]
   #       tres.xyDashPattern = dsh[c]
   #       plot[ip] = ngl.xy(wks, time_list[c], data_list[c], tres)


   #----------------------------------------------------------------------------
   # Set strings at top of plot
   #----------------------------------------------------------------------------
   var_str = var[v]
   # if var[v]=="PRECT" : var_str = "Precipitation [mm/day]"
   # if var[v]=="TMQ"   : var_str = "Column Water Vapor [mm]"

   lft_str = ''
   ctr_str = ''
   # if var[v] in ['PRECT','PRECC','PRECL'] : ctr_str = 'Mean: '+'%.2f'%avg_X+' [mm/day]'


   lat_chk,lon_chk = 'lat1' in locals(), 'lon1' in locals()


   if not lat_chk and not lon_chk : 
      ctr_str = 'Global'
   else:
      ctr_str += f' {region} '
      # if lat_chk:      ctr_str += f' {lat1}:{lat2}N '
      # if lon_chk:      ctr_str += f' {lon1}:{lon2}E '

   
   if mask_flag[v]=='lnd': ctr_str += ' (Land Only)'
   if mask_flag[v]=='ocn': ctr_str += ' (Ocean Only)'

   if overlay_cases:
      hs.set_subtitles(wks, plot[ip], var_str, ctr_str, '', font_height=0.01)
   else:
      for c in range(tmp_num_case):
         ip = c*num_var + v
         hs.set_subtitles(wks, plot[ip], case_name[c], ctr_str, var_str, font_height=0.02)

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
