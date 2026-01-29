#---------------------------------------------------------------------------------------------------
# Plot the zonal mean of the specified variables
#---------------------------------------------------------------------------------------------------
import os, ngl, xarray as xr, numpy as np, glob
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
#---------------------------------------------------------------------------------------------------
use_obs_remap = True
# obs_data_path = '/gpfs/alpine/scratch/hannah6/cli115/Obs/CERES-EBAF'; obs_remap_sub='clim_ne30pg2'
obs_data_path = os.getenv('HOME')+'/Data/Obs/CERES'; obs_remap_sub='clim_ne30pg2'
scrip_file_path = os.getenv('HOME')+'/E3SM/data_grid/ne30pg2_scrip.nc'
scrip_ds = xr.open_dataset(scrip_file_path)
#---------------------------------------------------------------------------------------------------
case,name,case_dir,case_sub = [],[],[],[]
clr,dsh,mrk = [],[],[]
def add_case(case_in,n=None,c='black',d=0,m=1,p=None,s=None):
   global name,case,clr,dsh,mrk
   if n is None: n = '' 
   case.append(case_in); name.append(n)
   case_dir.append(p); case_sub.append(s)
   clr.append(c); dsh.append(d); mrk.append(m)
#---------------------------------------------------------------------------------------------------
pvar,lev_list = [],[]
def add_var(var_name,lev=-1): pvar.append(var_name); lev_list.append(lev)
#---------------------------------------------------------------------------------------------------
# add_case('CERES-EBAF', n='CERES-EBAF' ,c='black', d=1)
# add_case('E3SM.RAD-SENS.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_128.00', n='nx_rad=128',c='red')
# add_case('E3SM.RAD-SENS.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_64.00',  n='nx_rad=64' ,c='orange')
# add_case('E3SM.RAD-SENS.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_32.00',  n='nx_rad=32' ,c='gold')
# add_case('E3SM.RAD-SENS.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_16.00',  n='nx_rad=16' ,c='palegreen')
# add_case('E3SM.RAD-SENS.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_8.00',   n='nx_rad=8'  ,c='green')
# add_case('E3SM.RAD-SENS.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_4.00',   n='nx_rad=4'  ,c='cyan')
# add_case('E3SM.RAD-SENS.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_2.00',   n='nx_rad=2'  ,c='blue')
# add_case('E3SM.RAD-SENS.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_1.00',   n='nx_rad=1'  ,c='purple')

add_case('CERES-EBAF', n='CERES-EBAF' ,c='black', d=2)
tmp_path,tmp_sub = '/pscratch/sd/w/whannah/e3sm_scratch/pm-gpu','archive/atm/hist'
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_1.MCICA_OFF',   n='nx/rx= 64/ 1 mcica OFF',c='red',      d=1,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_2.MCICA_OFF',   n='nx/rx= 64/ 2 mcica OFF',c='orange',   d=1,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_4.MCICA_OFF',   n='nx/rx= 64/ 4 mcica OFF',c='gold',     d=1,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_8.MCICA_OFF',   n='nx/rx= 64/ 8 mcica OFF',c='palegreen',d=1,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_16.MCICA_OFF',  n='nx/rx= 64/16 mcica OFF',c='green',    d=1,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_32.MCICA_OFF',  n='nx/rx= 64/32 mcica OFF',c='blue',     d=1,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_64.MCICA_OFF',  n='nx/rx= 64/64 mcica OFF',c='purple',   d=1,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_1.MCICA_ON',    n='nx/rx= 64/ 1 mcica ON ',c='red',      d=0,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_2.MCICA_ON',    n='nx/rx= 64/ 2 mcica ON ',c='orange',   d=0,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_4.MCICA_ON',    n='nx/rx= 64/ 4 mcica ON ',c='gold',     d=0,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_8.MCICA_ON',    n='nx/rx= 64/ 8 mcica ON ',c='palegreen',d=0,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_16.MCICA_ON',   n='nx/rx= 64/16 mcica ON ',c='green',    d=0,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_32.MCICA_ON',   n='nx/rx= 64/32 mcica ON ',c='blue',     d=0,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_64.MCICA_ON',   n='nx/rx= 64/64 mcica ON ',c='purple',   d=0,p=tmp_path,s=tmp_sub)

# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_1.MCICA_OFF',  n='nx/rx= 128/  1 mcica OFF',c='red',      d=1,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_2.MCICA_OFF',  n='nx/rx= 128/  2 mcica OFF',c='orange',   d=1,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_4.MCICA_OFF',  n='nx/rx= 128/  4 mcica OFF',c='gold',     d=1,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_8.MCICA_OFF',  n='nx/rx= 128/  8 mcica OFF',c='palegreen',d=1,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_16.MCICA_OFF', n='nx/rx= 128/ 16 mcica OFF',c='green',    d=1,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_32.MCICA_OFF', n='nx/rx= 128/ 32 mcica OFF',c='blue',     d=1,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_64.MCICA_OFF', n='nx/rx= 128/ 64 mcica OFF',c='purple',   d=1,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_128.MCICA_OFF',n='nx/rx= 128/128 mcica OFF',c='purple',d=1,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_1.MCICA_ON',   n='nx/rx= 128/  1 mcica ON ',c='red',      d=0,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_2.MCICA_ON',   n='nx/rx= 128/  2 mcica ON ',c='orange',   d=0,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_4.MCICA_ON',   n='nx/rx= 128/  4 mcica ON ',c='gold',     d=0,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_8.MCICA_ON',   n='nx/rx= 128/  8 mcica ON ',c='palegreen',d=0,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_16.MCICA_ON',  n='nx/rx= 128/ 16 mcica ON ',c='green',    d=0,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_32.MCICA_ON',  n='nx/rx= 128/ 32 mcica ON ',c='blue',     d=0,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_64.MCICA_ON',  n='nx/rx= 128/ 64 mcica ON ',c='purple',   d=0,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_128.MCICA_ON', n='nx/rx= 128/128 mcica ON ',c='purple',d=0,p=tmp_path,s=tmp_sub)

### extended runs
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_1.MCICA_OFF',   n='nx/rx= 64/ 1 mcica OFF',c='red',      d=1,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_8.MCICA_OFF',   n='nx/rx= 64/ 8 mcica OFF',c='palegreen',d=1,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_64.MCICA_OFF',  n='nx/rx= 64/64 mcica OFF',c='purple',   d=1,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_1.MCICA_ON',    n='nx/rx= 64/ 1 mcica ON ',c='red',      d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_8.MCICA_ON',    n='nx/rx= 64/ 8 mcica ON ',c='palegreen',d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_64.MCICA_ON',   n='nx/rx= 64/64 mcica ON ',c='purple',   d=0,p=tmp_path,s=tmp_sub)

#---------------------------------------------------------------------------------------------------
# add_var('PRECT')
# # add_var('TMQ')
# add_var('LHFLX')
# # add_var('SHFLX')
# add_var('P-E')
# add_var('TS')
# add_var('TGCLDLWP')
# add_var('TGCLDIWP')

# add_var('NET_TOA_RAD')

add_var('FSNTOA');add_var('FLUT')
# add_var('FSNS');add_var('FLNS')
add_var('LWCF');add_var('SWCF')
# add_var('CLDLOW')
# add_var('CLDHGH')
# add_var('CLDTOT')

#---------------------------------------------------------------------------------------------------

num_plot_col = 2

htype,years,months,first_file,num_files = 'h0',[],[],0,5*12

fig_file = 'figs/FXX-zonal-mean'

bin_dlat = 2

plot_diff = False

print_rmse  = False
print_stats = True

chk_significance = False # use with plot_diff

add_legend = False

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_pvar,num_case = len(pvar),len(case)

wkres = ngl.Resources()
npix=4096; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks('png',fig_file,wkres)
plot = [None]*num_pvar
res = hs.res_xy()
res.vpHeightF = 0.3
res.xyLineThicknessF = 4

if 'clr' not in locals(): 
   if num_case>1 : clr = np.linspace(2,len( ngl.retrieve_colormap(wks) )-1,num_case,dtype=int)
   else : clr = ['black']

# if num_case>1 and 'dsh' not in vars(): dsh = np.arange(0,num_case,1)
if 'dsh' not in locals(): 
   if num_case>1 : dsh = np.zeros(num_case)
   else : dsh = [0]
# res.xyLineColors   = clr
# res.xyDashPatterns = dsh

# res.tiXAxisString = 'Latitude'
res.tiXAxisString = 'sin( Latitude )'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
msg_list = []
for v in range(num_pvar):
   print('  var: '+hc.tcolor.MAGENTA+pvar[v]+hc.tcolor.ENDC)
   data_list,std_list,cnt_list = [],[],[]
   bin_list = []
   glb_avg_list = []
   if 'lev_list' in locals(): lev = lev_list[v]
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)
      #-------------------------------------------------------------------------
      if case[c]!='CERES-EBAF':
         case_obj = he.Case( name=case[c], data_dir=case_dir[c], data_sub=case_sub[c], populate_files=True )
      #-------------------------------------------------------------------------
      if case[c]=='CERES-EBAF':
         if use_obs_remap:
            file_list = sorted(glob.glob(f'{obs_data_path}/{obs_remap_sub}/*'))
         else:
            file_list = sorted(glob.glob(f'{obs_data_path}/clim_180x360/*'))
         file_list = file_list[:num_files]
         # ds = xr.open_mfdataset(file_list,combine='by_coords',concat_dim='time').isel(nv=0)
         ds = xr.open_mfdataset(file_list,combine='by_coords').isel(nv=0)

         tvar = pvar[v]
         if tvar=='NET_TOA_RAD': tvar = 'RESTOA'
         if pvar[v]=='NET_CF': tvar = 'SWCF'
         
         data = ds[tvar]
         
         if pvar[v]=='NET_CF': data = data + ds['LWCF']

         if use_obs_remap:
            lat,lon = ds['lat'],ds['lon']
            # scrip_ds = xr.open_dataset(os.getenv('HOME')+f'/E3SM/data_grid/{CERES_remap_scrip}.nc')
            area = scrip_ds['grid_area'].rename({'grid_size':'ncol'})
         else:
            data = data.stack(ncol=('lat','lon'))
            lon2D = np.transpose( np.repeat( ds['lon'].values[...,None],len(ds['lat']),axis=1) )
            lat2D =               np.repeat( ds['lat'].values[...,None],len(ds['lon']),axis=1)
            lon = xr.DataArray(lon2D,dims=('lat','lon')).stack(ncol=('lat','lon'))
            lat = xr.DataArray(lat2D,dims=('lat','lon')).stack(ncol=('lat','lon'))
            fv_scrip_ds = xr.open_dataset(os.getenv('HOME')+'/E3SM/data_grid/cmip6_180x360_scrip.20181001.nc')
            area = fv_scrip_ds['grid_area'].rename({'grid_size':'ncol'})
            # convert ncol to DataArray - multiindex causes problems
            data['ncol'] = np.arange(data.shape[1])
            area['ncol'] = np.arange(data.shape[1])
            lat['ncol']  = np.arange(data.shape[1])
            lon['ncol']  = np.arange(data.shape[1])

         # scrip_ds = xr.open_dataset(os.getenv('HOME')+'/E3SM/data_grid/cmip6_180x360_scrip.20181001.nc')
         area = scrip_ds['grid_area'].rename({'grid_size':'ncol'})

         # avoid dask array
         data.load()

         # convert ncol to DataArray - multiindex causes problems
         data['ncol'] = np.arange(data.shape[1])
         area['ncol'] = np.arange(data.shape[1])
         lat['ncol']  = np.arange(data.shape[1])
      #-------------------------------------------------------------------------
      else:
         if 'lon1' in locals() : case_obj.lon1 = lon1
         if 'lon2' in locals() : case_obj.lon2 = lon2
         if 'lev'  in locals() : case_obj.lev  = lev

         lat  = case_obj.load_data('lat',  htype=htype,num_files=1)
         area = case_obj.load_data('area', htype=htype,num_files=1).astype(np.double)
         data = case_obj.load_data(pvar[v],htype=htype,years=years,months=months,
                                   first_file=first_file,num_files=num_files)
      #-------------------------------------------------------------------------
      hc.print_time_length(data.time,indent=(' '*6),print_span=True, print_length=False)

      
      if 'area' in locals() :
         gbl_mean = ( (data.mean(dim='time')*area).sum() / area.sum() ).values 
         glb_avg_list.append(gbl_mean)
      else:
         lat1D,lon1D = ds['lat'],ds['lon']
         lat2D = np.repeat( lat1D.values[...,None],len(lon1D),axis=1)
         weights = np.cos( np.deg2rad( lat2D ) )
         gbl_mean = np.average( np.multiply( data.mean(dim='time').values, weights ) )
         glb_avg_list.append(gbl_mean)

      if print_stats:
         msg = hc.print_stat(data,name=pvar[v],stat='naxsh',indent=(' '*6),compact=True)
         msg_list.append('  case: '+case[c]+'\n'+msg)
         if 'area' in locals() :
            print(f'      Area Weighted Global Mean : '+hc.tcolor.CYAN+f'{gbl_mean:6.4}'+hc.tcolor.ENDC)
            

      if print_rmse:
         if c==0:baseline = data
         if c>0:
            rmse = np.sqrt( np.mean( np.square( data.to_masked_array() - baseline.to_masked_array() )))
            print(f'      Root Mean Square Error    : {rmse:6.4}')
            # exit()
      #-------------------------------------------------------------------------
      # Calculate time and zonal mean
      #-------------------------------------------------------------------------
      bin_ds = hc.bin_YbyX( data.mean(dim='time'), lat, 
                           bin_min=-90, bin_max=90, 
                           bin_spc=bin_dlat, wgt=area )

      data_list.append( bin_ds['bin_val'].values )
      lat_bins = bin_ds['bins'].values

      sin_lat_bins = np.sin(lat_bins*np.pi/180.)

      bin_list.append(sin_lat_bins)

      if 'area' in locals(): del area

   #----------------------------------------------------------------------------
   # Take difference from first case
   #----------------------------------------------------------------------------
   if plot_diff :
      data_tmp = data_list
      data_baseline = data_list[0]
      for c in range(num_case): data_list[c] = data_list[c] - data_baseline
   #----------------------------------------------------------------------------
   # Check significance using t-test
   # https://stattrek.com/hypothesis-test/difference-in-means.aspx
   #----------------------------------------------------------------------------
   if plot_diff and chk_significance :
      for c in range(1,num_case):

         N0,N1 = cnt_list[0],cnt_list[c]
         # using number of months for N might make more sense?
         if num_files>0: N0,N1 = num_files,num_files  
         if len(years)>0: N0,N1 = len(years)*12,len(years)*12
         S0,S1 = std_list[0],std_list[c]
         
         # Standard error
         SE = np.sqrt( S0**2/N0 + S1**2/N1 )

         hc.print_stat(SE,name='SE',indent='    ')

         # Degrees of freedom
         DF = (S0**2/N0 + S1**2/N1)**2       \
             /( ( (S0**2/N0)**2 / (N0-1) )   \
               +( (S1**2/N1)**2 / (N1-1) ) )

         # t-statistic - aX is the difference now
         t_stat = data_list[c] / SE

         hc.print_stat(t_stat,name='t statistic',indent='    ')

         # Critical t-statistic
         t_crit = 2.24   # 2-tail test w/ inf dof & P=0.05

         for i in range(len(lat_bins)):
            msg = f'  lat: {lat_bins[i]}   t_stat: {t_stat[i]}   '
            if np.absolute(t_stat[i])>t_crit: msg = msg+tcolor.RED+'SIGNIFICANT'+tcolor.ENDC
            print(msg)

   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   unit_str = ''
   if pvar[v] in ['PRECT','PRECC','PRECL']   : unit_str = '[mm/day]'
   if pvar[v] in ['LHFLX','SHFLX']           : unit_str = '[W/m2]'
   res.tiYAxisString = unit_str

   res.trXMinF = -1. #np.min( sin_lat_bins )
   res.trXMaxF =  1. #np.max( sin_lat_bins )

   res.trYMinF = np.min([np.nanmin(d) for d in data_list])
   res.trYMaxF = np.max([np.nanmax(d) for d in data_list])

   lat_tick = np.array([-90,-60,-30,0,30,60,90])
   res.tmXBMode = "Explicit"
   res.tmXBValues = np.sin( lat_tick*3.14159/180. )
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

   var_str = pvar[v]
   if pvar[v]=="PRECT" : var_str = "Precipitation"
   hs.set_subtitles(wks, plot[v], "", "", var_str, font_height=0.005)


#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
# layout = [num_pvar,1]
layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
# layout = [1,num_pvar]
# if num_pvar==4 : layout = [2,2]
# if num_pvar==6 : layout = [3,2]
pres = hs.setres_panel()
pres.nglFrame = False
pres.nglPanelRight = 0.5
ngl.panel(wks,plot,layout,pres)
#-------------------------------------------------------------------------------
# Add legend
#-------------------------------------------------------------------------------
if num_case>1 and add_legend:
   lgres = ngl.Resources()
   lgres.vpWidthF           = 0.05
   lgres.vpHeightF          = 0.01*num_case
   lgres.lgLabelFontHeightF = 0.005
   lgres.lgLabelFont        = "courier"
   lgres.lgMonoDashIndex    = False
   lgres.lgLineLabelsOn     = False
   lgres.lgLineThicknessF   = 16
   lgres.lgLabelJust        = 'CenterLeft'
   lgres.lgLineColors       = clr
   lgres.lgDashIndexes      = dsh

   labels = name
   max_label_len = max([len(n) for n in name])+2
   for n,nn in enumerate(name):
      labels[n] = f'  {nn:{max_label_len}}'
      if num_pvar==1: labels[n] += f'  {glb_avg_list[n]:6.1f}'

   # ndc_T,ndc_B,ndc_L,ndc_R = ngl.get_bounding_box(plot[0])
   # ndcx = ndc_R + 0.02
   # ndcy = np.average(np.array([ndc_T,ndc_B]))

   ndcx,ndcy = 0.5,0.55

   pid = ngl.legend_ndc(wks, len(labels), labels, ndcx, ndcy, lgres)

   # legend_id_list = [None]*num_pvar
   # for v in range(num_pvar):
   #    # ndcx, ndcy = ngl.datatondc(plot[v], 1., 0. )
   #    ndc_T,ndc_B,ndc_L,ndc_R = ngl.get_bounding_box(plot[v])
   #    ndcx = ndc_R + 0.02
   #    ndcy = np.average(np.array([ndc_T,ndc_B]))
   #    print(f'v: {v}   ndcx / ndcy: {ndcx} / {ndcy}')
   #    legend_id_list[v] = ngl.legend_ndc(wks, len(labels), labels, ndcx, ndcy, lgres)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
ngl.frame(wks)
ngl.end()

# print()
# for msg in msg_list: print(msg)

hc.trim_png(fig_file)
# print(f'\n{fig_file}.png\n')
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------