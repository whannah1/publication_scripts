import os, copy, ngl, xarray as xr, numpy as np, cmocean
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
tmp_path_hst_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_v2  = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
scrip_file_path  = os.getenv('HOME')+'/E3SM/data_grid/ne30pg2_scrip.nc'
# add_case('HadSST',                                                 n='HadSST',  c='black')
add_case('v2.LR.historical_0101',                                  n='E3SMv2',  c='cyan',    p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0151',                                  n='E3SMv2',  c='orange', p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0201',                                  n='E3SMv2',  c='green',  p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0251',                                  n='E3SMv2',  c='purple',   p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0301',                                  n='E3SMv2',  c='pink', p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical',                                       n='E3SMv2',  c='red',  p=None, s=None)
add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',c='blue', p=tmp_path_hst_mmf,s='archive/atm/hist')
#---------------------------------------------------------------------------------------------------

add_var('TS')

# htype,num_files = 'h0',12*30 ; first_file,first_file_v2 = 12*30,12*80
htype,num_files = 'h0',12*65 ; first_file,first_file_v2 = 12*00,12*50

# htype,yr1,yr2 = 'ha',1950,2014

fig_file,fig_type = 'figs_ENSO/FXX-ENSO-sst-variance.v1','png'

write_file    = False
print_stats   = True
overlay_cases = False

remove_annual_cycle    = True

recalculate = False

add_obs_TS = True

year_start = 1950+first_file/12

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
num_var = len(var)
num_case = len(case)

# overlay_cases with single case causes segfault
if num_case==1 and not add_obs_TS : overlay_cases = False

if 'clr' not in vars() or clr==[]: clr = ['black']*num_case
if 'dsh' not in vars() or dsh==[]: dsh = np.zeros(num_case)

if 'lev' not in vars(): lev = np.array([0])


#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
if add_obs_TS:
   plot = [None]*((num_case+1)*num_var)
else:
   plot = [None]*(num_case*num_var)

wkres = ngl.Resources()
npix=1024*2; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
if 'legend_file' in locals(): lgd_wks = ngl.open_wks('png',legend_file,wkres)

res = hs.res_contour_fill_map()
res.tmYLLabelFontHeightF         = 0.008
res.tmXBLabelFontHeightF         = 0.008
res.lbLabelFontHeightF           = 0.015
res.tmXBOn                       = False
res.tmYLOn                       = False

res.mpMinLatF = -20
res.mpMaxLatF = 20
res.mpMinLonF = 130
res.mpMaxLonF = 280




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

      tmp_file = os.getenv('HOME')+f'/Research/E3SM/data_temp/ENSO.sst-variance.v1.{case[c]}.{var[v]}.nc'
      # '.lat_{lat1}_{lat2}.lon_{lon1}_{lon2}.nc'

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
         first_file_tmp = first_file_v2 if 'v2.LR.historical' in case[c] else first_file

         lat  = case_obj.load_data('lat',  component=get_comp(case[c]),htype=htype)
         lon  = case_obj.load_data('lon',  component=get_comp(case[c]),htype=htype)
         area = case_obj.load_data('area',component=get_comp(case[c]),htype=htype,num_files=1).astype(np.double)
         data = case_obj.load_data(tvar,  component=get_comp(case[c]),htype=htype,
                                    first_file=first_file_tmp,
                                    num_files=num_files,
                                    lev=lev)

         #-------------------------------------------------------------------------
         #reset time index to start at zero and convert to days
         dtime = ( data['time'][-1] - data['time'][0] ).values.astype('timedelta64[D]')
         print('      Time length: '+str(dtime)+'  ('+str(dtime.astype('timedelta64[M]'))+')')

         #-------------------------------------------------------------------------
         # time/space means
         # avg_X = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )

         # print(); print(data)

         if remove_annual_cycle:
            for n in range(12):
               month_index = slice(n,len(data.time),12)
               data[month_index,:] = data.isel(time=month_index) - data.isel(time=month_index).mean(dim='time')
         
         # convert to anomalies
         data = data - data.mean()

         # detrend in time
         fit = xr.polyval(data['time'], data.polyfit(dim='time', deg=1).polyfit_coefficients)
         data = data - fit

         data = data.var(dim='time')

         #-------------------------------------------------------------------------
         # Write to file 
         #-------------------------------------------------------------------------
         if os.path.isfile(tmp_file) : os.remove(tmp_file)
         tmp_ds = xr.Dataset( coords=data.coords )
         tmp_ds[var[v]] = data
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         tmp_ds = xr.open_dataset( tmp_file )
         data = tmp_ds[var[v]]

      
      
      #-------------------------------------------------------------------------
      if print_stats: hc.print_stat(data,name='',stat='naxs',indent=' '*6,compact=True)
         
      data_list.append( data.values )


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
         obs_scrip_file = os.getenv('HOME')+'/E3SM/data_grid/180x360_HadSST_scrip.nc'
         # ncremap -G ttl='Equi-Angular grid 180x360'#latlon=180,360#lat_typ=uni#lat_drc=n2s#lon_typ=180_wst  -g ~/E3SM/data_grid/180x360_HadSST_scrip.nc

      if obs_name=='NOAA':
         lat_name,lon_name = 'lat','lon'
         print(' '*4+'Loading NOAA OI SST data...')
         file_obs = os.getenv('HOME')+'/Data/Obs/NOAA/monthly/sst.mnmean.nc'
         ds_obs = xr.open_dataset(file_obs)
         data_obs = ds_obs['sst'] * ds_obs['scale_factor'] + ds_obs['add_offset']

      # print(ds_obs)
      # exit()

      # xy_dims = (lon_name,lat_name)
      # xlon, ylat = np.meshgrid(ds_obs[lat_name],ds_obs[lon_name])
      # R = hc.earth_radius(ylat)
      # drlat = np.deg2rad(np.gradient(ylat, axis=0))
      # drlon = np.deg2rad(np.gradient(xlon, axis=1))

      # dy,dx = drlat * R , drlon * R * np.cos(np.deg2rad(ylat))
      # area_obs = np.absolute(dy*dx) / np.square(R) # calculate area and convert to steridians
      # area_obs = xr.DataArray(area_obs,dims=xy_dims).transpose()

      # xlon = np.transpose(xlon)
      # ylat = np.transpose(ylat)

      # tlon1,tlon2 = lon1,lon2
      # if tlon1>180: tlon1 = tlon1 - 360
      # if tlon2>180: tlon2 = tlon2 - 360
      # lon1,lon2 = tlon1,tlon2

      # tmp_data = np.ones([len(ds_obs[lat_name]),len(ds_obs[lon_name])],dtype=bool)
      # tmp_coords = {lat_name:ds_obs[lat_name],lon_name:ds_obs[lon_name]}
      # mask = xr.DataArray( tmp_data, coords=tmp_coords, dims=(lat_name,lon_name) )
      # mask = mask & (ds_obs[lat_name]>=lat1) & (ds_obs[lat_name]<=lat2)
      # mask = mask & (ds_obs[lon_name]>=lon1) & (ds_obs[lon_name]<=lon2)
      # # mask = mask & (ylat>=lat1) & (ylat<=lat2) & (xlon>=lon1) & (xlon<=lon2)


      # print(); print(ds_obs[lon_name])
      # print(); print(lon1)
      # print(); print(lon2)
      # exit()

      # data_obs.load()
      # area_obs.load()

      # data_obs = data_obs.where( mask, drop=True)
      # area_obs = area_obs.where( mask, drop=True)

      # print(); print(area_obs)

      # print();hc.print_stat(area_obs,name='area_obs',stat='naxs',indent=' '*6,compact=True)
      # print();hc.print_stat(data_obs,name='data_obs',stat='naxs',indent=' '*6,compact=True)
      # exit()

      # gbl_mean_obs = ( (data_obs*area_obs).sum(dim=xy_dims) / area_obs.sum(dim=xy_dims) )

      # if remove_annual_cycle:
      #    for n in range(12):
      #       month_index = slice(n,len(gbl_mean_obs.time),12)
      #       gbl_mean_obs[month_index] = gbl_mean_obs[month_index] - np.mean(gbl_mean_obs[month_index])

      # hc.print_stat(gbl_mean_obs,name=obs_name,stat='naxs',indent=' '*6,compact=True)
      # exit()

      time_obs = ds_obs['time.year'] + ds_obs['time.dayofyear']/365
      time_obs = time_obs - year_start

      # limit extent of obs data to match size of model data
      sim_num_t = num_files # np.max([len(d) for d in time_list])

      t_start = None
      for t,y in enumerate(time_obs):
         if y>=0: t_start=t;break

      if t_start is None: raise ValueError('ERROR: no appropriate value found for t_start')

      # gbl_mean_obs = gbl_mean_obs.isel(time=slice(t_start,t_start+sim_num_t))
      # time_obs     = time_obs    .isel(time=slice(t_start,t_start+sim_num_t))

      # data_list.insert(0, gbl_mean_obs.values )
      # # time_list.insert(0, time_obs.values )

      data_obs = data_obs.isel(time=slice(t_start,t_start+sim_num_t))

      if remove_annual_cycle:
            for n in range(12):
               month_index = slice(n,len(data_obs.time),12)
               data_obs[month_index,:] = data_obs.isel(time=month_index) - data_obs.isel(time=month_index).mean(dim='time')

      # convert to anomalies
      data_obs = data_obs - data_obs.mean()

      # detrend in time
      fit = xr.polyval(data_obs['time'], data_obs.polyfit(dim='time', deg=1).polyfit_coefficients)
      data_obs = data_obs - fit
         
      data_obs = data_obs.var(dim='time')

      data_obs = data_obs.stack(ncol=('latitude','longitude'))
      # data_obs = data_obs.stack(ncol=('longitude','latitude'))

      data_list.insert(0, data_obs.values )
      case_name.insert(0,obs_name)
      clr.insert(0,'black')
      dsh.insert(0,0)

   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   # data_min = np.min([np.nanmin(d) for d in data_list])
   # data_max = np.max([np.nanmax(d) for d in data_list])

   tres = copy.deepcopy(res)
   # tres.tiXAxisString = 'Time [days]'

   # tres.cnFillPalette = "MPL_viridis"
   tres.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
   # tres.cnFillPalette = np.array( cmocean.cm.amp(np.linspace(0,1,256)) )
   tres.cnLevelSelectionMode = "ExplicitLevels"
   tres.cnLevels = np.linspace(0,3,21)

   tmp_num_case = num_case
   if var[v]=='TS' and add_obs_TS: tmp_num_case = num_case + 1

   pdum = [None]*tmp_num_case*3

   for c in range(tmp_num_case):
      
      # hs.set_cell_fill(tres,case_obj)

      if (add_obs_TS and c>0) or (not add_obs_TS):
         scrip_file_name = os.getenv('HOME')+'/Research/E3SM/data_grid/ne30pg2_scrip.nc'
      else:
         scrip_file_name = obs_scrip_file
         # scrip_file_name = os.getenv('HOME')+'/E3SM/data_grid/cmip6_180x360_scrip.20181001.nc'

      scrip_ds = xr.open_dataset(scrip_file_name)
      tres.cnFillMode       = "CellFill"
      tres.sfXArray         = scrip_ds.variables['grid_center_lon'].values
      tres.sfYArray         = scrip_ds.variables['grid_center_lat'].values
      tres.sfXCellBounds    = scrip_ds.variables['grid_corner_lon'].values
      tres.sfYCellBounds    = scrip_ds.variables['grid_corner_lat'].values

      ip = v*num_case+c
      # ip = c*num_var+v

      plot[ip] = ngl.contour_map(wks,np.ma.masked_invalid(data_list[c]),tres)


      ### draw box around ENSO regions
      pgres = ngl.Resources()
      pgres.nglDraw,pgres.nglFrame = False,False
      

      pgres.gsLineColor = 'blue'
      pgres.gsLineDashPattern = 0
      pgres.gsLineThicknessF = 4
      region,lat1,lat2,lon1,lon2 = 'Nino3'  ,-5,5,360-150,360-90
      bx = np.array([lon1,lon2,lon2,lon1,lon1])
      by = np.array([lat1,lat1,lat2,lat2,lat1])
      pdum[ip+1] = ngl.add_polyline(wks,plot[ip],bx,by,pgres)

      pgres.gsLineColor = 'red'
      pgres.gsLineDashPattern = 0
      pgres.gsLineThicknessF = 4
      region,lat1,lat2,lon1,lon2 = 'Nino4'  ,-5,5,160,360-150
      bx = np.array([lon1,lon2,lon2,lon1,lon1])
      by = np.array([lat1,lat1,lat2,lat2,lat1])
      pdum[ip+1] = ngl.add_polyline(wks,plot[ip],bx,by,pgres)

      pgres.gsLineColor = 'green'
      pgres.gsLineDashPattern = 2
      pgres.gsLineThicknessF = 6
      region,lat1,lat2,lon1,lon2 = 'Nino3.4',-5,5,190,240
      bx = np.array([lon1,lon2,lon2,lon1,lon1])
      by = np.array([lat1,lat1,lat2,lat2,lat1])
      pdum[ip+0] = ngl.add_polyline(wks,plot[ip],bx,by,pgres)

      #-------------------------------------------------------------------------
      # Set strings at top of plot
      #-------------------------------------------------------------------------
      # var_str = f'{var[v]} variance'
      var_str = f'SST variance'
      # if var[v]=="PRECT" : var_str = "Precipitation [mm/day]"
      # if var[v]=="TMQ"   : var_str = "Column Water Vapor [mm]"

      lft_str = ''
      ctr_str = ''
      # lat_chk,lon_chk = 'lat1' in locals(), 'lon1' in locals()

      # if not lat_chk and not lon_chk : 
      #    ctr_str = 'Global'
      # else:
      #    ctr_str += f' {region} '
      #    if lat_chk:      ctr_str += f' {lat1}:{lat2}N '
      #    if lon_chk:      ctr_str += f' {lon1}:{lon2}E '
      
      if mask_flag[v]=='lnd': ctr_str += ' (Land Only)'
      if mask_flag[v]=='ocn': ctr_str += ' (Ocean Only)'

      hs.set_subtitles(wks, plot[ip], case_name[c], ctr_str, var_str, font_height=0.02)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
hc.printline()

if add_obs_TS:
   # layout = [num_var,num_case+1]
   layout = [num_case+1,num_var]
else:
   # layout = [num_var,num_case]
   layout = [num_case,num_var]


ngl.panel(wks,plot[0:len(plot)],layout,hs.setres_panel())
ngl.end()

hc.trim_png(fig_file)
if 'legend_file' in locals(): hc.trim_png(legend_file)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
