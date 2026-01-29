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

var,lev_list,mask_flag = [],[],[]
def add_var(var_name,lev=-1,mask=None): 
   var.append(var_name); lev_list.append(lev),mask_flag.append(mask)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

### coupled runs for Jim Benedict
# add_case('E3SM.PI-CPL.v1.ne30.01',  n='E3SMv1')
# add_case('E3SM.PI-CPL.v2.ne30.01',  n='E3SMv2')

### AMIP test on perlmutter
# add_case('E3SM.AMIP.GNUGPU.ne30pg2_r05_oECv3.F20TR-MMFXX-CMIP6.MOMFB.BVT.00',  n='E3SM-MMF AMIP', p='/pscratch/sd/w/whannah/e3sm_scratch/perlmutter',s='run')


### INCITE 2022 coupled runs
# scrip_file_path = os.getenv('HOME')+'/E3SM/data_grid/ne30pg2_scrip.nc'
# # add_case('CERES-EBAF',n='CERES-EBAF',c='gray',d=1,ref=True); obs_remap_sub='clim_ne30pg2'
# # add_case('ERA5',n='ERA5',c='gray',d=1,ref=True); obs_remap_sub='monthly_ne30pg2'; obs_data_path = '/gpfs/alpine/scratch/hannah6/cli115/Obs/ERA5'
# add_case('E3SM.INCITE2022-CPL.ne30pg2_EC30to60E2r2.WCYCL1950',     n='E3SM (1950)',c='orange')
# add_case('E3SM.INCITE2022-CPL.ne30pg2_EC30to60E2r2.WCYCL1950-MMF1',n='MMF  (1950)',c='cyan')
# add_case('E3SM.INCITE2022-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR',     n='E3SM (20TR)',c='red' )
# add_case('E3SM.INCITE2022-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='MMF  (20TR)',c='blue')

### 4xCO2 tests on Summit
# add_case('E3SM.2023-CO2-TEST-00.GNUGPU.ne30pg2_oECv3.WCYCL1850-MMF1.1xCO2',n='MMF 1xCO2',c='blue')
# add_case('E3SM.2023-CO2-TEST-00.GNUGPU.ne30pg2_oECv3.WCYCL1850-MMF1.4xCO2',n='MMF 4xCO2',c='red')
tmp_path,tmp_sub = '/gpfs/alpine/cli115/proj-shared/hannah6/e3sm_scratch','archive/atm/hist'
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2',n='E3SM-MMF PI 1xCO2',c='blue' ,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.2xCO2',n='E3SM-MMF PI 2xCO2',c='green',p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.4xCO2',n='E3SM-MMF PI 4xCO2',c='red'  ,p=tmp_path,s=tmp_sub)


### validation for v3+L80
# add_case('20230629.v3alpha02.amip.chrysalis.L80',             n='AMIP L80', c='red',p='/lcrc/group/e3sm/ac.whannah/E3SMv3_dev',s='archive/atm/hist')
# add_case('20231002.v3alpha04_bigrid_L80_QBO1.F2010.chrysalis',n='F2010 L80',c='blue',p='/lcrc/group/acme/ac.benedict/E3SMv3_dev',s='run')

### SCREAM LR Cess tests
# add_case('E3SM.SCREAM-CESS-LR.ne30pg2_oECv3.F2010-SCREAM-LR.SSTP_0K.NN_64',n='SCREAMv0 +0K',d=0,c='blue')
# add_case('E3SM.SCREAM-CESS-LR.ne30pg2_oECv3.F2010-SCREAM-LR.SSTP_4K.NN_64',n='SCREAMv0 +4K',d=0,c='red')
# add_case('E3SM.SCREAM-CESS-LR.ne30pg2_oECv3.F2010.SSTP_0K.NN_64',          n='E3SMv2 +0K',  d=1,c='blue')
# add_case('E3SM.SCREAM-CESS-LR.ne30pg2_oECv3.F2010.SSTP_4K.NN_64',          n='E3SMv2 +4K',  d=1,c='red')

### 2022 coupled historical runs
# add_case('E3SM.INCITE2022-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR',     n='E3SM (20TR)',c='red' ,p='/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu',s='run')
# add_case('E3SM.INCITE2022-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='MMF  (20TR)',c='blue')
# add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='MMF  (20TR)',c='blue',s='archive/atm/hist')

### PAM dev tests
# add_case('E3SM.PAM-DEV-30.GNUGPU.F2010-MMF-PAM-C.ne30pg2_oECv3.CDT_05.F-QTOT',n='MMF+PAM (crm_dt=5)')
# add_case('E3SM.PAM-DEV-2023-35.GNUGPU.ne4pg2_ne4pg2.FRCE-MMF1',             n='RCE-MMF1',c='red')
# add_case('E3SM.PAM-DEV-2023-35.GNUGPU.ne4pg2_ne4pg2.FRCE-MMF2.CDT_10.DPP_4',n='RCE-MMF2',c='blue')

### L80 ensemble
# tmp_path,tmp_sub = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu','archive/atm/hist'
# add_case('E3SM.2023-SCIDAC.ne30pg2_EC30to60E2r2.AMIP.EF_0.35.CF_10.HD_1.00',c='red',n='',p=tmp_path,s=tmp_sub) # top-tier
# add_case('E3SM.2023-SCIDAC.ne30pg2_EC30to60E2r2.AMIP.EF_0.35.CF_10.HD_0.50',c='orange',n='',p=tmp_path,s=tmp_sub) # top-tier
# add_case('E3SM.2023-SCIDAC.ne30pg2_EC30to60E2r2.AMIP.EF_0.20.CF_15.HD_0.50',c='green',n='',p=tmp_path,s=tmp_sub) # top-tier
# add_case('E3SM.2023-SCIDAC.ne30pg2_EC30to60E2r2.AMIP.EF_0.05.CF_25.HD_0.50',c='blue',n='',p=tmp_path,s=tmp_sub) # top-tier
# add_case('E3SM.2023-SCIDAC.ne30pg2_EC30to60E2r2.AMIP.EF_0.90.CF_20.HD_0.25',c='cyan',n='',p=tmp_path,s=tmp_sub) # top-tier
# add_case('E3SM.2023-SCIDAC.ne30pg2_EC30to60E2r2.AMIP.EF_0.40.CF_10.HD_1.00',c='purpl',n='',p=tmp_path,s=tmp_sub) # top-tier
# add_case('E3SM.2023-SCIDAC.ne30pg2_EC30to60E2r2.AMIP.EF_0.60.CF_07.HD_1.35',c='pink',n='',p=tmp_path,s=tmp_sub) # top-tier

# add_var('TS')
# add_var('TS',mask='lnd')
# add_var('TS',mask='ocn')

# add_var('PS')
# add_var('TMQ')
# add_var('TCO')
# add_var('PRECT')
# add_var('SHFLX')
# add_var('LHFLX')
# add_var('P-E')

# add_var('TGCLDLWP')
# add_var('TGCLDIWP')

# add_var('TBOT'); 
# add_var('QBOT'); 
# add_var('SHFLX')
# add_var('LHFLX')

# add_var('PRECT',mask='lnd/ocn ratio')
# add_var('TGCLDLWP',mask='lnd/ocn ratio')

# add_var('TMQ'  ,mask='lnd')
# add_var('TMQ'  ,mask='ocn')
# add_var('PRECT',mask='lnd')
# add_var('PRECT',mask='ocn')
# add_var('LHFLX',mask='lnd')
# add_var('LHFLX',mask='ocn')
# add_var('TGCLDLWP',mask='lnd')
# add_var('TGCLDLWP',mask='ocn')
# add_var('TGCLDIWP',mask='lnd')
# add_var('TGCLDIWP',mask='ocn')

# add_var('PRECT',mask='lnd')
# add_var('PRECT',mask='ocn')

# tlev = -59
# add_var('T',lev=tlev,mask='lnd')
# add_var('Q',lev=tlev,mask='ocn')








# add_var('UBOT'); 
# add_var('TAUX')
# add_var('VBOT'); 
# add_var('TAUY')
# add_var('TMQ')

# add_var('SOLIN')
# add_var('FLNT')
# add_var('FSNT')
# add_var('NET_TOA_RAD')

# add_var('CLDLOW')
# add_var('CLDHGH')
# add_var('SWCF')
# add_var('LWCF')

# add_var('U',lev=850)
# add_var('U',lev=200)

# tlev = -56
# add_var('T',lev=tlev)
# add_var('Q',lev=tlev)
# tlev = -57
# add_var('T',lev=tlev)
# add_var('Q',lev=tlev)
# tlev = -58
# add_var('T',lev=tlev)
# add_var('Q',lev=tlev)
# tlev = -59
# add_var('T',lev=tlev)
# add_var('Q',lev=tlev)
# add_var('U',lev=tlev)
# add_var('V',lev=tlev)
# add_var('CLOUD',lev=tlev)
# add_var('CLDLIQ',lev=tlev)
# add_var('MMF_TLS',lev=tlev); add_var('MMF_QTLS',lev=tlev)
# add_var('MMF_DT' ,lev=tlev); add_var('MMF_DQ'  ,lev=tlev)
# add_var('MMF_DU' ,lev=tlev); add_var('MMF_DV'  ,lev=tlev)
# add_var('U' ,lev=tlev); add_var('MMF_DU' ,lev=tlev); add_var('TAUX');

# add_var('T',lev=-71); add_var('SHFLX')
# add_var('Q',lev=-71); add_var('LHFLX')
# add_var('U',lev=-71); add_var('TAUX')
# add_var('V',lev=-71); add_var('TAUY')


htype,years,months,first_file,num_files = 'ha',[],[],0,0 # pre-calculated annual means
# htype,years,months,first_file,num_files = 'h0',[],[],0,0
# htype,years,months,first_file,num_files = 'h1',[],[],0,0
# htype,years,months,first_file,num_files = 'h0',[],[],0*12,10*12
# htype,years,months,first_file,num_files = 'h1',[],[],0,int(365/2)

# lat1,lat2 = -30,30


fig_type = "png"
fig_file = os.getenv('HOME')+'/Research/E3SM/figs_clim/clim.timeseries.global_mean.v1'

write_file    = False
print_stats   = True
overlay_cases = True

convert_to_daily_mean  = False
convert_to_annual_mean = True

add_obs_TS = False

add_trend = True

num_plot_col  = 2

# year_start = 1950 + first_file/12
# year_start = 1985 + first_file/12
year_start = 00 + first_file/12

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
num_var = len(var)
num_case = len(case)

# # overlay_cases with single case causes segfault
# if num_case==1 : overlay_cases = False

if 'clr' not in vars() or clr==[]: clr = ['black']*num_case
if 'dsh' not in vars() or dsh==[]: dsh = np.zeros(num_case)

if 'lev' not in vars(): lev = np.array([0])

#-------------------------------------------------------------------------------
# plot legend in separate file
#-------------------------------------------------------------------------------
if num_case>1:
   legend_file = fig_file+'.legend'
   wkres = ngl.Resources() ; npix = 2048 ; wkres.wkWidth,wkres.wkHeight=npix,npix
   lgd_wks = ngl.open_wks('png',legend_file,wkres)
   lgres = ngl.Resources()
   lgres.vpWidthF           = 0.06
   lgres.vpHeightF          = 0.06#*num_case
   lgres.lgLabelFontHeightF = 0.008
   lgres.lgLabelFont        = "courier"
   lgres.lgMonoDashIndex    = False
   lgres.lgLineLabelsOn     = False
   lgres.lgLineThicknessF   = 4#16
   lgres.lgLabelJust        = 'CenterLeft'
   lgres.lgLineColors       = clr
   lgres.lgDashIndexes      = dsh

   labels = case_name
   for i in range(len(labels)): labels[i] = ' '*4+labels[i] 

   pid = ngl.legend_ndc(lgd_wks, len(labels), labels, 0.5, 0.65, lgres)

   ngl.frame(lgd_wks)
   hc.trim_png(legend_file)
   exit()

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
if overlay_cases:
   plot = [None]*(num_var)
else:
   plot = [None]*(num_case*num_var)

wkres = ngl.Resources()
npix=1024; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

if 'legend_file' in locals(): lgd_wks = ngl.open_wks('png',legend_file,wkres)


res = hs.res_xy()
res.vpHeightF = 0.5
# res.vpHeightF = 0.2
res.tmYLLabelFontHeightF         = 0.015
res.tmXBLabelFontHeightF         = 0.015
res.tiXAxisFontHeightF           = 0.015
res.tiYAxisFontHeightF           = 0.015
res.xyLineThicknessF = 5

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

      print('    case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)

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
      lat = case_obj.load_data('lat',  htype=htype)
      lon = case_obj.load_data('lon',  htype=htype)
      area = case_obj.load_data('area',htype=htype,num_files=1).astype(np.double)
      data = case_obj.load_data(tvar,  htype=htype, \
                                first_file=first_file,\
                                num_files=num_files,lev=lev)
      
      # time = data.time - data.time[0]

      time_bnds = case_obj.load_data('time_bnds',htype=htype,first_file=first_file,num_files=num_files)
      time = time_bnds.isel(nbnd=0)

      # if htype in ['h0'] and convert_to_annual_mean:
      #    num_time = len(data.time)
      #    #print(num_time); print(np.floor(num_time/12)*12); exit()
      #    data = data.isel(time=slice(0,int(np.floor(num_time/12)*12)))

      # deal with land data dimension names
      if 'levgrnd' in data.dims: 
         data = data.rename({'lndgrid':'ncol'})
         area = area.rename({'lndgrid':'ncol'})

      # land or ocean mask
      if mask_flag[v] is not None :
         land_frac = case_obj.load_data('LANDFRAC',component=get_comp(case[c]),htype='rh0',first_file=first_file,num_files=num_files).astype(np.double)
         if 'time' in land_frac.dims: land_frac = land_frac.isel(time=0)
         land_frac = land_frac / land_frac.max()
         # if mask_flag[v]=='lnd': data = data*land_frac
         # if mask_flag[v]=='ocn': data = data*(1-land_frac)

         mask = xr.DataArray( np.ones(land_frac.shape,dtype=bool), dims=land_frac.dims )
         if mask_flag[v]=='lnd': mask = mask & (land_frac.values>0.5)
         if mask_flag[v]=='ocn': mask = mask & (land_frac.values<0.5)
         data = data.where( mask, drop=True)


      data = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )

      # if print_stats: hc.print_stat(data,name='',stat='naxs',indent=' '*6,compact=True)

      # Convert to daily mean
      if htype in ['h1','h2'] and convert_to_daily_mean: 
         data = data.resample(time='D').mean(dim='time')

      if htype in ['h0'] and convert_to_annual_mean: 
         # truncate months past the last full year
         extra_months = len(time)%12
         print(f'  extra_months: {extra_months}')
         if extra_months !=0: 
            data = data.isel(time=slice(0,-1-extra_months))
            time = time.isel(time=slice(0,-1-extra_months))
         month_length = time.dt.days_in_month
         month_length['time'] = time
         data['time'] = time
         mn_wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
         data = (data*mn_wgts).resample(time='A').sum('time') / (mn_wgts).resample(time='A').sum(dim='time')


      if np.all( lev < 0 ) and 'lev' in data.coords : print(f'    lev value: {data.lev.values}')


      #-------------------------------------------------------------------------
      if 'levgrnd' in data.dims : data = data.sum(dim='levgrnd')
      # if 'levgrnd' in data.dims : data = data.isel(levgrnd=9)     # deepest layer
      
      # Get rid of lev dimension
      if 'lev' in data.dims : data = data.isel(lev=0)

      #-------------------------------------------------------------------------
      # #reset time index to start at zero and convert to days
      # dtime = ( data['time'][-1] - data['time'][0] ).values.astype('timedelta64[D]')
      # print('      Time length: '+str(dtime)+'  ('+str(dtime.astype('timedelta64[M]'))+')')

      #-------------------------------------------------------------------------
      # Make time start at zero
      data['time'] = ( data['time'] - data['time'][0] ).astype('float') / 86400e9 / 365 + year_start
      
      # # convert to anomaly
      # avg_X = avg_X - avg_X.mean()

      #-------------------------------------------------------------------------
      # redefine values as anomalies
      if var[v]=='TS' and add_obs_TS:
         # use these years to redefine simulation data as anomalies relative to a climatology
         anomaly_yr1,anomaly_yr2 = 1961,1990 # match HadCRU
         clim_num_years = anomaly_yr2 - anomaly_yr1 + 1
         yr = data['time'].values
         for t,y in enumerate(yr):
            if y>=year_start: t_start=t;break

         data = data - data.isel(time=slice(t_start,t_start+clim_num_years)).mean()
         


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
   # Add Obs data
   #----------------------------------------------------------------------------
   if var[v]=='TS' and add_obs_TS:
      print(' '*4+'Loading HadCRU data...')
      # file_obs = '/global/cfs/cdirs/m3312/whannah/obs_data/HadCRU/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc'
      file_obs = '/gpfs/alpine/scratch/hannah6/cli115/Obs/HadCRU/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc'
      ds_obs = xr.open_dataset(file_obs)
      data_obs = ds_obs['tas_mean']
      # hc.print_stat(data_obs,name='HadCRU',stat='naxs',indent=' '*6,compact=True)
      xy_dims = ('longitude','latitude')
      xlon, ylat = np.meshgrid(ds_obs['latitude'],ds_obs['longitude'])
      R = hc.earth_radius(ylat)
      dlat = np.deg2rad(np.gradient(ylat, axis=0))
      dlon = np.deg2rad(np.gradient(xlon, axis=1))
      dy,dx = dlat * R , dlon * R * np.cos(np.deg2rad(ylat))
      area_obs = np.absolute(dy*dx) / np.square(R) # calculate area and convert to steridians
      area_obs = xr.DataArray(area_obs,dims=xy_dims).transpose()
      gbl_mean_obs = ( (data_obs*area_obs).sum(dim=xy_dims) / area_obs.sum(dim=xy_dims) )

      if convert_to_annual_mean: 
         time_obs = ds_obs['time_bnds'].isel(bnds=0)
         month_length = time_obs.dt.days_in_month
         month_length['time'] = time_obs
         gbl_mean_obs['time'] = time_obs
         wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
         gbl_mean_obs = (gbl_mean_obs*wgts).resample(time='A').sum('time') / (wgts).resample(time='A').sum(dim='time')

      hc.print_stat(gbl_mean_obs,name='HadCRU',stat='naxs',indent=' '*6,compact=True)

      if convert_to_annual_mean:
         # time_obs = ds_obs['time.year'].isel[slice(0,0,365)]
         time_obs = ds_obs['time.year'].resample(time='Y').mean(dim='time')
         # time_obs = time_obs.resample(time='Y').mean(dim='time')
         # time_obs_yr = ds_obs['time.year']
      else:
         time_obs = ds_obs['time.year'] + ds_obs['time.dayofyear']/365
      
      # time_obs = time_obs - year_start

      # limit extent of obs data to match size of model data
      sim_num_t = np.max([len(d) for d in time_list])


      # print()
      # print(time_obs.values)
      # print()
      # print(time_obs_yr.values)
      # print()
      # print(ds_obs['time.year'].values)
      # print()
      # exit()

      for t,y in enumerate(time_obs):
         # if y>=0: t_start=t;break
         if y>=year_start: t_start=t;break

      # print(); print(f't_start: {t_start}'); print()

      gbl_mean_obs = gbl_mean_obs.isel(time=slice(t_start,t_start+sim_num_t))
      time_obs     = time_obs    .isel(time=slice(t_start,t_start+sim_num_t))

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
   tres.tiXAxisString = 'Time [years]'
   # tres.tiXAxisString = 'Time [years]'
   if convert_to_annual_mean: tres.tiXAxisString = 'Time [years]'

   # reset start year for all data
   if 'year_start' in locals():
      for t in range(len(time_list)):
         time_list[t] = time_list[t] + year_start

   ### Make sure plot bounds are consistent
   # tres.trYMinF = 14.4
   tres.trYMinF = np.min([np.nanmin(d) for d in data_list])
   tres.trYMaxF = np.max([np.nanmax(d) for d in data_list])
   tres.trXMinF = np.min([np.nanmin(d) for d in time_list])
   tres.trXMaxF = np.max([np.nanmax(d) for d in time_list])

   print()
   print(f'tres.trYMinF = {tres.trYMinF}')
   print(f'tres.trYMaxF = {tres.trYMaxF}')
   print(f'tres.trXMinF = {tres.trXMinF}')
   print(f'tres.trXMaxF = {tres.trXMaxF}')
   print()

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
      tplot = ngl.xy(wks, time_list[c], data_list[c], tres)
      
      if overlay_cases: 
         if c==0: 
            plot[ip] = tplot
         else:
            ngl.overlay(plot[ip],tplot)
      else:
         plot[ip] = tplot

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
      if mask_flag[v] is None: ctr_str = 'Global'
      if mask_flag[v]=='lnd' : ctr_str = 'Land Only'
      if mask_flag[v]=='ocn' : ctr_str = 'Ocean Only'
      if mask_flag[v]=='lnd/ocn ratio':ctr_str = 'Land / Ocean Ratio'
   else:
      if lat_chk:      ctr_str += f' {lat1}:{lat2}N '
      if lon_chk:      ctr_str += f' {lon1}:{lon2}E '
   
   

   if overlay_cases:
      hs.set_subtitles(wks, plot[ip], var_str, ctr_str, '', font_height=0.012)
   else:
      hs.set_subtitles(wks, plot[ip], case_name[c], ctr_str, var_str, font_height=0.015)

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
