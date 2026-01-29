# v2 is a cleaned up version of v1 - functionality is mostly the same
import os, ngl, subprocess as sp, numpy as np, xarray as xr, dask, copy, string, cmocean, numba, glob
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
host = hc.get_host()
#-------------------------------------------------------------------------------
name,case,case_dir,case_sub = [],[],[],[]
def add_case(case_in,n=None,p=None,s=None,g=None,c=None):
   global name,case,case_dir,case_sub
   tmp_name = case_in if n is None else n
   case.append(case_in); name.append(tmp_name); 
   case_dir.append(p); case_sub.append(s); 
#-------------------------------------------------------------------------------
var,ovar,var_str = [],[],[]
def add_var(var_in,ovar_in,name=None): 
   var.append(var_in)
   ovar.append(ovar_in)
   var_str.append(var_in if name is None else name)
#-------------------------------------------------------------------------------
scrip_file_path  = os.getenv('HOME')+'/E3SM/data_grid/ne30pg2_scrip.nc'
tmp_path_hst_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_v2  = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
tmp_path_co2_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
# obs_file_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology_1985-2014/ERA5/ERA5_ANN_198501_201412_climo.nc'
# obs_file_path = '/global/cfs/cdirs/m3312/whannah/2023-CPL/ERA5/ERA5_ANN_198501_201412_climo.ne30pg2.nc'
tmp_sub = 'archive/atm/hist'
add_case('ERA5',                                                   n='ERA5')
add_case('v2.LR.historical_0101',                                  n='E3SMv2',  c='cyan',    p=tmp_path_hst_v2, s='archive/atm/hist')
add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',c='blue', p=tmp_path_hst_mmf,s='archive/atm/hist')
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2',n='1xCO2',c='blue' ,p=tmp_path_co2_mmf,s=tmp_sub)
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.2xCO2',n='2xCO2',c='green',p=tmp_path_co2_mmf,s=tmp_sub)
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.4xCO2',n='4xCO2',c='red'  ,p=tmp_path_co2_mmf,s=tmp_sub)
#-------------------------------------------------------------------------------

fig_file,fig_type = 'figs/F03-zonal-mean-UTQ','png'
tmp_file_head = 'data/zonal-mean'

add_var('U',None,name='Zonal Wind')
# add_var('V',None,name='Meridional Wind')
add_var('T','ta',name='Temperature')
add_var('Q','hus',name='Specific Humidity')

#-------------------------------------------------------------------------------
# plev = np.array([   1.,    2.,    3.,    5.,    7.,   10.,   20.,   30.,   50.,   70.,  100.,  125.,  150.,])
# plev = np.array([10,30,50,75,100,125,150,200,250,300,350,400,450,500,550,600,650,700,750,800,825,850,875,900,925,950,975,1000])
# plev = np.array([1,2,4,6,10,30,50,100,150,200,300,400,500,600,700,800,850,925,950])
# plev = np.array([5,10,30,50,100,150,200,300,400,500,600,700,800,850,925,975,1000])
# plev = np.array([10,20,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900])
# plev = np.array([1,2,3,5,7,10,20,30,50,70,100,125,150,175,200,225,250,300,350,400,450,500,550,600,650,700,750,775,800,825,850,875,900,925,950,975,1000.])
plev = np.array([1,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950])
#-------------------------------------------------------------------------------

# htype,first_file,num_files = 'ha',35,30
htype,yr1,yr2 = 'ha',1980,2014
# htype,yr1,yr2 = 'ha',1980,1981

recalculate = False

plot_diff,add_diff   = True,False
print_stats          = True
var_x_case           = False
num_plot_col         = len(var)
use_common_label_bar = False

lat1, lat2, dlat = -88., 88., 2

#---------------------------------------------------------------------------------------------------
# Set up plot resources
if case==[]: raise ValueError('ERROR - case list is empty!')
num_var,num_case = len(var),len(case)

subtitle_font_height = 0.01


if 'scrip_file_path' not in locals(): scrip_file_path = None

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix

wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*(num_var*num_case)
   
res = hs.res_contour_fill()
res.vpHeightF = 0.3
res.trYReverse = True
res.tiXAxisString = 'Latitude'
res.tiYAxisString = 'Pressure [hPa]'

res.tmYLLabelFontHeightF   = 0.015
res.tmXBLabelFontHeightF   = 0.015
res.lbLabelFontHeightF     = 0.02
res.lbRightMarginF         = -1
res.lbLeftMarginF          = -1
res.lbBottomMarginF        = 0.5

cres = hs.res_contour()
cres.cnLineThicknessF       = 1
# cres.cnLineColor            = 'gray'

# res.trYMaxF = 900

# lev_tick = np.array([850,500,200,100,50,10])
lev_tick = np.array([900,700,500,300,100])
res.tmYLMode = "Explicit"
res.tmYLValues = lev_tick
res.tmYLLabels = lev_tick

#---------------------------------------------------------------------------------------------------
def calculate_obs_area(lon,lat,lon_bnds,lat_bnds):
   re = 6.37122e06  # radius of earth
   nlat,nlon = len(lat),len(lon)
   area = np.empty((nlat,nlon),np.float64)
   for j in range(nlat):
      for i in range(nlon):
         dlon = np.absolute( lon_bnds[j,1] - lon_bnds[j,0] )
         dlat = np.absolute( lat_bnds[j,1] - lat_bnds[j,0] )
         dx = re*dlon*np.pi/180.
         dy = re*dlat*np.pi/180.
         area[j,i] = dx*dy
   return area
#---------------------------------------------------------------------------------------------------
@numba.njit()
def apply_pressure_mask(data,lev,psfc3D,ntime,nlev,ncol):
   data_masked = data
   for t in range(ntime):
      for i in range(ncol):
         for k in range(nlev):
            if lev[k] > psfc3D[t,k,i]:
               data_masked[t,k,i] = np.nan
   return data_masked
#---------------------------------------------------------------------------------------------------
@numba.njit()
def bin_numba(data, lat, lat_bins, nbin, nlev, ncol, bin_val, bin_cnt):
   for k in range( nlev ):
      for b in range( nbin ):
         bin_bot = lat_bins[b] - dlat/2. 
         bin_top = lat_bins[b] - dlat/2. + dlat
         for n in range( ncol ):
            if lat[n]>=bin_bot and lat[n]<bin_top :
               if np.isfinite(data[k,n]):
                  bin_cnt[b,k] = bin_cnt[b,k]+1
                  bin_val[b,k] = ( bin_val[b,k]*bin_cnt[b,k] + data[k,n] ) / ( bin_cnt[b,k] + 1 )
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   lev_list = []
   lat_list = []
   data_list = []
   psfc_list = []
   glb_avg_list = []
   for c in range(num_case):
      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}.nc'
      print(' '*4+f'case: {hc.tclr.CYAN}{case[c]}{hc.tclr.END}  =>  {tmp_file}')
      if recalculate:
         scrip_ds = xr.open_mfdataset(scrip_file_path).rename({'grid_size':'ncol'})
         area = scrip_ds['grid_area']
         lat  = scrip_ds['grid_center_lat']
         #----------------------------------------------------------------------
         # if case[c]=='ERA5':
         #    ds = xr.open_dataset(obs_file_path)
         #    area = calculate_obs_area(ds['lon'].values,ds['lat'].values,ds['lon_bnds'].values,ds['lat_bnds'].values)
         #    area = xr.DataArray( area, coords=[ds['lat'],ds['lon']] ).stack(ncol=('lat','lon'))
         #    data = ds['ta'].sel(plev=plev*1e2).isel(time=0).stack(ncol=('lat','lon')).rename({'plev':'lev'})
         #    lat2D = np.repeat( ds['lat'].values[...,None],len(ds['lon']),axis=1)
         #    lat   = xr.DataArray(lat2D,dims=('lat','lon')).stack(ncol=('lat','lon'))
         #    print(); print(data); print(); print(area); print(); print(lat); exit()
         # if case[c]=='ERA5':
         #    scrip_ds = xr.open_dataset(scrip_file_path)
         #    area = scrip_ds['grid_area'].rename({'grid_size':'ncol'})
         #    ds = xr.open_dataset(obs_file_path)
         #    lat = ds['lat']
         #    data = ds[ovar[v]].isel(time=0).sel(plev=plev*1e2).rename({'plev':'lev'})
         #    data['lev'] = data['lev']/1e2
         #    psfc = ds['ps'].isel(time=0)/1e2
         if case[c]=='ERA5':
            era5_data_root = '/global/cfs/cdirs/m3312/whannah/2023-CPL/ERA5/annual'
            #-------------------------------------------------------------------
            if 'file_list' in locals(): del file_list
            if var[v]=='T' : file_list = sorted(glob.glob(f'{era5_data_root}/e5.oper.an.pl.128_130_t.*'))
            if var[v]=='U' : file_list = sorted(glob.glob(f'{era5_data_root}/e5.oper.an.pl.128_131_u.*'))
            if var[v]=='V' : file_list = sorted(glob.glob(f'{era5_data_root}/e5.oper.an.pl.128_132_v.*'))
            if var[v]=='Q' : file_list = sorted(glob.glob(f'{era5_data_root}/e5.oper.an.pl.128_133_q.*'))
            ds = xr.open_mfdataset(file_list)
            data = ds[var[v]].rename({'level':'lev'})
            data = data.sel(lev=plev) # select levels specified above (no interpolation)
            #-------------------------------------------------------------------
            file_list = sorted(glob.glob(f'{era5_data_root}/e5.oper.an.sfc.128_134_sp.*'))
            ds = xr.open_mfdataset(file_list)
            psfc = ds['SP']/1e2
            #-------------------------------------------------------------------
            data = data.where( data['time.year']>=yr1, drop=True)
            data = data.where( data['time.year']<=yr2, drop=True)
            psfc = psfc.where( psfc['time.year']>=yr1, drop=True)
            psfc = psfc.where( psfc['time.year']<=yr2, drop=True)
            #-------------------------------------------------------------------
            # print(); print(data)
            # print(); print(psfc)
            # hc.print_stat(data,name=var[v],stat='naxsh',indent='    ',compact=True)
            # hc.print_stat(psfc,name='PS',  stat='naxsh',indent='    ',compact=True)
            #-------------------------------------------------------------------
            psfc3D = psfc.expand_dims(dim={'lev':data['lev']},axis=1)
            data_masked = apply_pressure_mask(data.values,plev,psfc3D.values,
                                              len(data['time'].values),
                                              len(data[ 'lev'].values),
                                              len(data['ncol'].values) )
            data = xr.DataArray( np.ma.masked_invalid(data_masked) , coords=data.coords )
            #-------------------------------------------------------------------
            data = data.mean(dim='time',skipna=True)
            # psfc = psfc.mean(dim='time',skipna=True)
         else:
            # data_dir_tmp,data_sub_tmp = None, None
            # if case_dir[c] is not None: data_dir_tmp = case_dir[c]
            # if case_sub[c] is not None: data_sub_tmp = case_sub[c]
            # case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp )
            # case_obj.set_coord_names(var[v])
            # #-------------------------------------------------------------------
            # # read the data
            # with dask.config.set(**{'array.slicing.split_large_chunks': True}):
            #    lat  = case_obj.load_data('lat',  htype=htype)
            #    area = case_obj.load_data('area', htype=htype).astype(np.double)
            #    data = case_obj.load_data(var[v],htype=htype,ps_htype=htype,lev=plev,
            #                                    first_file=first_file,num_files=num_files,
            #                                    extrap_flag=False)
            #    hc.print_time_length(data.time,indent=' '*6)
            #    data = data.mean(dim='time',skipna=True)
            #    psfc = case_obj.load_data('PS',htype=htype,ps_htype=htype,first_file=first_file,num_files=num_files)
            #    psfc = psfc.mean(dim='time',skipna=True)/1e2
            #    #----------------------------------------------------------------
            #    lev2D  =               np.repeat( plev       [...,None],len(psfc),axis=1)
            #    psfc2D = np.transpose( np.repeat( psfc.values[...,None],len(plev),axis=1) )
            #    data_masked = apply_pressure_mask(data.values,lev2D,psfc2D)
            #    data = xr.DataArray( np.ma.masked_invalid(data_masked) , coords=data.coords )
            #    # data = np.ma.masked_invalid(data_masked)
            #-------------------------------------------------------------------
            file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/*.eam.{htype}.*'
            file_list = sorted(glob.glob(file_path))
            ds = xr.open_mfdataset( file_list )
            #-------------------------------------------------------------------
            ds = ds.where( ds['time.year']>=yr1, drop=True)
            ds = ds.where( ds['time.year']<=yr2, drop=True)
            #-------------------------------------------------------------------
            data = ds[var[v]]
            psfc = ds['PS']
            data = he.interpolate_to_pressure(ds,data_mlev=data,lev=plev,ds_ps=ds,ps_var='PS'
                                          ,interp_type=2,extrap_flag=True)
            #-------------------------------------------------------------------
            psfc3D = psfc.expand_dims(dim={'lev':data['lev']},axis=1)
            data_masked = apply_pressure_mask(data.values,plev,psfc3D.values,
                                              len(data['time'].values),
                                              len(data[ 'lev'].values),
                                              len(data['ncol'].values) )
            data = xr.DataArray( np.ma.masked_invalid(data_masked) , coords=data.coords )
            #-------------------------------------------------------------------
            data = data.mean(dim='time',skipna=True)
            # psfc = psfc.mean(dim='time',skipna=True)/1e2
            
         #----------------------------------------------------------------------
         # Special handling of various specific circumstances
         if var[v]=='Q': data = data*1e3
         #----------------------------------------------------------------------
         # print stats after time averaging
         if print_stats: 
            hc.print_stat(data,name=var[v],stat='naxsh',indent='    ',compact=True)
            # hc.print_stat(psfc,name='PS',  stat='naxsh',indent='    ',compact=True)
         #----------------------------------------------------------------------
         # # Calculate time and zonal mean
         # bin_ds = hc.bin_YbyX( data, lat, bin_min=lat1, bin_max=lat2, bin_spc=dlat, wgt=area, keep_lev=True )
         #----------------------------------------------------------------------
         lat_bins = np.arange(lat1,lat2+dlat/2,dlat)
         nbin = len(lat_bins)
         nlev = len(data['lev'].values)
         ncol = len(data['ncol'].values)
         bin_val = np.zeros([nbin,nlev])
         bin_cnt = np.zeros([nbin,nlev])
         bin_numba(data=data.values, lat=lat.values, 
                   lat_bins=lat_bins, nbin=nbin,nlev=nlev, ncol=ncol, 
                   bin_val=bin_val, bin_cnt=bin_cnt)
         bin_ds = xr.Dataset()
         coords = {'bins':lat_bins,'lev':plev}
         bin_ds['bin_val'] = xr.DataArray(bin_val,coords=coords)
         bin_ds['bin_cnt'] = xr.DataArray(bin_cnt,coords=coords)
         #----------------------------------------------------------------------
         # bin_psfc_ds = hc.bin_YbyX( psfc, lat, bin_min=lat1, bin_max=lat2, bin_spc=dlat, wgt=area )
         # bin_ds['psfc'] = bin_psfc_ds['bin_val']
         #----------------------------------------------------------------------
         bin_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         bin_ds = xr.open_dataset( tmp_file, use_cftime=True  )
      #-------------------------------------------------------------------------
      # if case[c]=='ERA5':
      #    cnt_mask = np.where( bin_ds['bin_cnt'].transpose().values<50, True, False )
      # bin_ds['bin_val'] = np.ma.masked_where( cnt_mask, data_list[c] )
      #-------------------------------------------------------------------------
      lev_list.append( bin_ds['lev'].values )
      lat_list.append( bin_ds['bins'].values )
      data_list.append( bin_ds['bin_val'].transpose().values )
      # psfc_list.append( bin_ds['psfc'].values )
   #----------------------------------------------------------------------------
   # Plot averaged data

   tres = copy.deepcopy(res)
   
   # res.trXMinF = np.min( lat_bins )
   # res.trXMaxF = np.max( lat_bins )
   # res.trYMinF = np.min([np.nanmin(d) for d in data_list])
   # res.trYMaxF = np.max([np.nanmax(d) for d in data_list])

   # lat_tick = np.array([-90,-60,-30,0,30,60,90])
   # res.tmXBMode = "Explicit"
   # res.tmXBValues = np.sin( lat_tick*3.14159/180. )
   # res.tmXBLabels = lat_tick

   # plot[v] = ngl.xy(wks, np.stack(bin_list), np.ma.masked_invalid(  np.stack(data_list) ), res)

   tres.cnLevelSelectionMode = "ExplicitLevels"
   cres.cnLevelSelectionMode = "ExplicitLevels"

   if var[v]=='U'     : clev_main = np.arange(-32,32+4,4)
   if var[v]=='V'     : clev_main = np.arange(-2,2+0.4,0.4)
   if var[v]=='OMEGA' : clev_main = np.arange(-0.075,0.075+0.015,0.015)
   if var[v]=='T'     : clev_main = np.arange(200,300+10,10)
   if var[v]=='Q'     : clev_main = np.arange(1,13+2,2)


   if var[v]=='U' : tres.cnFillPalette = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )
   if var[v]=='T' : tres.cnFillPalette = np.array( cmocean.cm.thermal(np.linspace(0,1,256)) )
   if var[v]=='Q' : tres.cnFillPalette = np.array( cmocean.cm.dense(np.linspace(0,1,256)) )
   

   # tres.lbTitlePosition = 'Bottom'
   # tres.lbTitleFontHeightF = 0.02
   # if var[v]=='U'     : tres.lbTitleString = '[m/s]'
   # if var[v]=='OMEGA' : tres.lbTitleString = '[Pa/s]'
   # if var[v]=='Q'     : tres.lbTitleString = '[g/kg]'
   # if var[v]=='T'     : tres.lbTitleString = '[K]'


   # # calculate common limits for consistent contour levels
   # data_min = np.min([np.nanmin(d) for d in data_list])
   # data_max = np.max([np.nanmax(d) for d in data_list])
   # cmin,cmax,cint = ngl.nice_cntr_levels(data_min, data_max, aboutZero=False )
   # clev = np.linspace(cmin,cmax,num=21)

   if plot_diff and num_case>1:
      if 'dlev' in locals(): del dlev
      # if var[v]=='T': dlev = np.linspace(-7.0,7.0,num=21)
      # if var[v]=='Q': dlev = np.linspace(-1.5,1.5,num=21)
      if var[v]=='U'     : dlev = np.arange(-9,9+2,2)
      if var[v]=='T'     : dlev = np.arange(-5.5,5.5+1,1)
      if var[v]=='Q'     : dlev = np.arange(-18,18+4,4)/1e1
      # if 'dlev' not in locals():
      #    tmp_data_list = [0]*len(data_list)
      #    for c in range(1,num_case): tmp_data_list[c] = data_list[c] - data_list[0]
      #    diff_min = np.min([np.min(d) for d in tmp_data_list])
      #    diff_max = np.max([np.max(d) for d in tmp_data_list])
      #    cmin,cmax,clev = ngl.nice_cntr_levels(diff_min, diff_max, aboutZero=True )
      #    dlev = np.linspace(cmin,cmax,num=21)
      #    # print(f'data_min: {data_min}')
      #    # print(f'data_max: {data_max}')
      #    print(f'diff_min: {diff_min}')
      #    print(f'diff_max: {diff_max}')

   

   # # apply mask based on time-mean surface pressure
   # for c in range(num_case):
   #    lev2D  =               np.repeat(  lev_list[c][...,None],len(psfc_list[c]),axis=1)
   #    psfc2D = np.transpose( np.repeat( psfc_list[c][...,None],len( lev_list[c]),axis=1) )
   #    psfc_mask = np.where( psfc2D<lev2D, True, False )
   #    data_list[c] = np.ma.masked_where( psfc_mask, data_list[c] )

   for c in range(num_case):
      tres.sfXArray,tres.sfYArray = lat_list[c],lev_list[c]
      

      # lev_tick = np.array([900,700,500,300,200,100,50])
      # tres.tmYLMode,tres.tmYLValues = 'Explicit',lev_tick

      # tres.cnLevels = clev

      tmp_data = data_list[c]
      if plot_diff:
         tres.cnLevels = clev_main
         if c==0:
            baseline = tmp_data
         if c>=1: 
            tmp_data = tmp_data - baseline
            tres.cnLevels = dlev
            tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
      
      ip = v*num_case+c if var_x_case else c*num_var+v

      plot[ip] = ngl.contour(wks, tmp_data, tres)

      # overlay contour lines from first case
      if plot_diff and c>=1:
         cres.sfXArray,cres.sfYArray = lat_list[c],lev_list[c]

         cres.cnLineThicknessF  = 1
         cres.cnLineDashPattern = 1
         cres.cnLevels          = clev_main
         ngl.overlay( plot[ip], ngl.contour(wks, data_list[0], cres) )

         cres.cnLineThicknessF  = 1
         cres.cnLineDashPattern = 0
         cres.cnLevels          = clev_main[ np.where(clev_main>=0) ]
         ngl.overlay( plot[ip], ngl.contour(wks, data_list[0], cres) )

         cres.cnLineThicknessF  = 3
         cres.cnLineDashPattern = 0
         cres.cnLevels          = [0]
         ngl.overlay( plot[ip], ngl.contour(wks, data_list[0], cres) )

      if var[v]=='U'     : unit_str = '[m/s]'
      if var[v]=='OMEGA' : unit_str = '[Pa/s]'
      if var[v]=='Q'     : unit_str = '[g/kg]'
      if var[v]=='T'     : unit_str = '[K]'

      rstr = f'{var_str[v]} {unit_str}' if not (plot_diff and c>=1) else f'{var_str[v]} Bias {unit_str}'
      # cstr = f'{yr1}-{yr2}'
      hs.set_subtitles(wks, plot[ip], name[c], '', rstr, font_height=subtitle_font_height)
      
#---------------------------------------------------------------------------------------------------
# Finalize plot

num_case_alt = num_case*2-1 if (plot_diff and add_diff) else num_case
layout = [num_var,num_case_alt] if var_x_case else [num_case_alt,num_var]


if not (plot_diff and add_diff):
   if num_case==1 or num_var==1:
      layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
   
pnl_res = hs.setres_panel()
if use_common_label_bar: pnl_res.nglPanelLabelBar = True
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.01

# pnl_res.nglPanelYWhiteSpacePercent = 5

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
