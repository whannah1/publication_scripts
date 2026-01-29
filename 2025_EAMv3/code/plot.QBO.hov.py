import os, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string, dask
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
import pandas as pd
import cmocean
#-------------------------------------------------------------------------------
case,case_name,case_dir,case_sub = [],[],[],[]
def add_case(case_in,n=None,p=None,s=None,):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   tmp_name = '' if n is None else n
   case.append(case_in); case_name.append(tmp_name)
   case_dir.append(p); case_sub.append(s)
#-------------------------------------------------------------------------------
# Create tarball of hov data:
# tar -czvf QBO.hov.2023_L80_ensemble.tar.gz data_temp/QBO.hov.* 
# tar -czvf QBO.hov.2023_AMIP.tar.gz data_temp/QBO.hov.v1.E3SM.2023-SCIDAC-v2-AMIP* 
#-------------------------------------------------------------------------------

tmp_path,tmp_sub = '/global/cfs/cdirs/m3312/whannah/e3smv3_amip','archive/atm/hist'
add_case('v3.LR.amip_0101', n='EAMv3', p=tmp_path,s=tmp_sub)

tmp_path,tmp_sub = '/global/cfs/cdirs/e3smdata/simulations/','archive/atm/hist'
add_case('v2.LR.amip_0101', n='EAMv2', p=tmp_path,s=tmp_sub) 

#-------------------------------------------------------------------------------
lev = np.array([ 0.1, 0.2, 0.5, 1.,  2.,  3.,  5.,  7., 
                 10., 20., 30., 50.,70.,100.,125.,150.,])
#-------------------------------------------------------------------------------
var = ['U']
#-------------------------------------------------------------------------------

fig_file,fig_type = 'figs/QBO-hov','png'

lat1,lat2 = -5,5
# lat1,lat2 = -10,10

htype,first_file,num_files = 'h0', (1979-1870)*12, 12*35#65


print_stats = True
print_time  = False

var_x_case = False

use_common_label_bar = True

num_plot_col = len(var)#1

add_obs = True
obs_case = 'ERA5'

recalculate_sim = False
recalculate_obs = False

# write_to_file = True

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)

if 'lev' not in vars(): lev = np.array([0])

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
num_plot = (num_case+1)*num_var if add_obs else num_case*num_var
plot = [None]*(num_plot)
res = hs.res_contour_fill()
res.vpHeightF = 0.1
# res.vpHeightF = 0.15
res.tmYLLabelFontHeightF         = 0.008
res.tmXBLabelFontHeightF         = 0.008
res.tiXAxisFontHeightF           = 0.008
res.tiYAxisFontHeightF           = 0.008
if use_common_label_bar:
   res.lbLabelBarOn = False

# # disable these by default - turn back on for bottom panel(s)
# res.tmXBOn = False
# res.tiXAxisOn = False

res.tiXAxisString = 'Time'
# if htype=='h0': res.tiXAxisString = 'Time [months]'
if htype=='h0': res.tiXAxisString = 'Time [years]'

res.tiYAxisString = 'Pressure [hPa]'
res.trYReverse = True
# res.trYLog = True # doesn't work due to irregular spacing :(
res.trYMinF = 5
res.trYMaxF = 100

tm_vals = [5,10,50,100,200]
# tm_vals = [2,4,8,10,20,40,80,100]
res.tmYLMode = 'Explicit'
res.tmYLValues = tm_vals
res.tmYLLabels = tm_vals
#---------------------------------------------------------------------------------------------------
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
def get_tmp_file_name(case_in,var_in):
   tmp_file = f'data/QBO.hov.v1.{case_in}.{var_in}.lat_{lat1}_{lat2}.nc'
   return tmp_file
#---------------------------------------------------------------------------------------------------
def write_file(case_in,var_in,data_avg,Z_avg=None):
   tmp_file = get_tmp_file_name(case_in,var_in)
   print('    writing to file: '+tmp_file)
   ds_out = xr.Dataset( coords=data_avg.coords )
   ds_out[var[v]] = data_avg
   ds_out.to_netcdf(path=tmp_file,mode='w')
#---------------------------------------------------------------------------------------------------
def load_file(case_in,var_in):
   tmp_file = get_tmp_file_name(case_in,var_in)
   ds = xr.open_dataset( tmp_file )
   data_avg = ds[var_in]
   return data_avg
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+var[v])
   tvar = var[v]
   data_list,time_list,lev_list = [],[],[]
   
   #----------------------------------------------------------------------------
   # read the model data
   #----------------------------------------------------------------------------
   
   for c in range(num_case):
      print(f'    case: {hc.tclr.CYAN}{case[c]}{hc.tclr.END}')

      if recalculate_sim:

         data_dir_tmp,data_sub_tmp = None, None
         # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
         if case_dir[c] is not None: data_dir_tmp = case_dir[c]
         if case_sub[c] is not None: data_sub_tmp = case_sub[c]

         case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp, time_freq=None )

         if 'lat1' in vars() : case_obj.lat1,case_obj.lat2 = lat1,lat2
         if 'lon1' in vars() : case_obj.lon1,case_obj.lon2 = lon1,lon2

         # avoid creating large chunks
         with dask.config.set(**{'array.slicing.split_large_chunks': True}):  

            data = case_obj.load_data(tvar,   htype=htype, first_file=first_file, num_files=num_files, lev=lev)
            area = case_obj.load_data('area', htype=htype, first_file=first_file, num_files=num_files).astype(np.double)

            hc.print_time_length(data.time,indent=' '*4)

            avg_dims= ('ncol')
            data_avg = ( (data*area).sum(dim=avg_dims) / area.sum(dim=avg_dims) )

            ### adjust time to represent the middle of the month instead of the end
            time = data_avg.time.values
            time_orig = copy.deepcopy(time)
            for i,t in enumerate(time):
               if i==0:
                  dt = pd.Timedelta('15 days')
               else:
                  dt = ( time_orig[i] - time_orig[i-1] ) / 2
               time[i] = time_orig[i] - dt
            data_avg['time'] = time
            
            time = data_avg['time.year'].values + data_avg['time.month'].values/12.

            if print_time: print(); print(time); print()

            write_file(case[c],var[v],data_avg)
      else:
         data_avg = load_file(case[c],var[v])

      data_list.append( data_avg.transpose().values )
      lev_list.append( data_avg['lev'].values )
      time = data_avg['time.year'].values + data_avg['time.month'].values/12.
      time_list.append( time )
   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   if var[v] in ['U','V','T','Q'] and add_obs:
      if os.path.exists( '/lcrc/group/e3sm' ):
         obs_scratch_path = '/lcrc/group/e3sm/diagnostics/observations/Atm/time-series'
      if os.path.exists( '/global/cfs/cdirs/e3sm' ):
         obs_scratch_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series'
      if obs_case=='ERA5': 
         if var[v]=='U': obs_var,obs_data_file =  'ua',f'{obs_scratch_path}/ERA5/ua_197901_201912.nc'
         if var[v]=='V': obs_var,obs_data_file =  'va',f'{obs_scratch_path}/ERA5/va_197901_201912.nc'
         if var[v]=='T': obs_var,obs_data_file =  'ta',f'{obs_scratch_path}/ERA5/ta_197901_201912.nc'
         if var[v]=='Q': obs_var,obs_data_file = 'hus',f'{obs_scratch_path}/ERA5/hus_197901_201912.nc'

      print(f'\n    case: {obs_case:10}  obs_data_file: {obs_data_file}')

      if recalculate_obs:

         ds = xr.open_dataset(obs_data_file)
         area = calculate_obs_area(ds['lon'].values,ds['lat'].values,ds['lon_bnds'].values,ds['lat_bnds'].values)
         area = xr.DataArray( area, coords=[ds['lat'],ds['lon']] )
         data = ds[obs_var]
         data = data.sel(lat=slice(lat1,lat2))
         area = area.sel(lat=slice(lat1,lat2))
         data_avg = (data*area).sum(dim=('lon','lat')) / area.sum(dim=('lon','lat'))

         time = data_avg['time.year'].values + data_avg['time.month'].values/12.

         ### truncate obs to match model data
         # time = data_avg['time.year'].values + data_avg['time.month'].values/12.
         sim_t1 = time_list[0][ 0]
         sim_t2 = time_list[0][-1]
         t_beg,t_end = None,None
         if sim_t1>1850:
            for i,t in enumerate(time):
               if t_beg is None and t==sim_t1: t_beg = i
               if t_end is None and t==sim_t2: t_end = i
         else:
            # looks like we have an F-compset - so just truncate length to match
            t_beg,t_end = 0,int(sim_t2*12-sim_t1*12+1)

         print(f't_beg / t_end: {t_beg} / {t_end}')

         if t_beg is not None and t_end is not None:
            data_avg = data_avg.isel(time=slice(t_beg,t_end))
            # time = data_avg['time.year'].values + data_avg['time.month'].values/12.
         else:
            exit(f'Something went wrong? t_beg: {t_beg}  t_end: {t_end}')

         # if print_time: print(); print(time); print()

         k_list = []
         for k,p in enumerate(data_avg['plev'].values/1e2):
            if p in lev: k_list.append(k)
         data_avg = data_avg.isel(plev=k_list)

         write_file(obs_case,var[v],data_avg)
      else:
         data_avg = load_file(obs_case,var[v])

      time = data_avg['time.year'].values + data_avg['time.month'].values/12.

      # if print_time: print(); print(time); print()

      if var[v]=='Q': data_avg = data_avg*1e3

      lev_list.insert(0, data_avg['plev'].values/1e2 )
      time_list.insert(0, time )
      data_list.insert(0, data_avg.transpose().values )
      if v==0: case_name.insert(0, obs_case )
   #----------------------------------------------------------------------------
   # print stats after time averaging
   print()
   num_case_tmp = (num_case+1) if add_obs else num_case
   for c in range(num_case_tmp):
      hc.print_stat(data_list[c],name=f'{case_name[c]} {var[v]}',stat='naxsh',indent='    ',compact=True)
   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)

   # tres.cnFillPalette = np.array( cmocean.cm.diff(np.linspace(0,1,256)) )
   # tres.cnFillPalette = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )
   tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )

   if var[v]=='Q': tres.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )

   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])

   tres.cnLevelSelectionMode = 'ExplicitLevels'

   if var[v]=='U': tres.cnLevels = np.linspace(-50,50,num=21)
   if var[v]=='T': tres.cnLevels = np.arange(190,260+5,5)
   if var[v]=='Q': tres.cnLevels = np.arange(2,40+2,2)*1e-4

   num_case_tmp = (num_case+1) if add_obs else num_case
   for c in range(num_case_tmp):
      
      time = time_list[c]
      # time = ( time - time[0] ).astype('float') / 86400e9
      # time = ( time - time[0] ).astype('float') / (60*60*24*365)

      tres.sfYArray = lev_list[c]#.astype(int)
      tres.sfXArray = time
      # tres.sfXArray = np.linspace( 1./12., float(num_files)/12., num=len(time) )

      if use_common_label_bar and c==(num_case_tmp-1): 
         tres.lbLabelBarOn = True
         # tres.lbTopMarginF       =  0.2
         tres.lbBottomMarginF    = 0.5
         tres.lbLeftMarginF      = -2
         tres.lbRightMarginF     = -2

      ip = c*num_var + v

      plot[ip] = ngl.contour(wks, data_list[c], tres)

      #-------------------------------------------------------------------------
      # Set strings at top of plot
      #-------------------------------------------------------------------------
      var_str = var[v]
      if var[v]=='PRECT' : var_str = 'Precipitation [mm/day]'
      if var[v]=='U'     : var_str = 'Zonal Wind [m/s]'

      ctr_str = ''
      # if var[v] in ['PRECT','PRECC','PRECL'] : ctr_str = 'Mean: '+'%.2f'%avg_X+' [mm/day]'
      # hs.set_subtitles(wks, plot[len(plot)-1], case_name[c], ctr_str, var_str, font_height=0.01)
      # hs.set_subtitles(wks, plot[len(plot)-1], case_name[c], '', '', font_height=0.008)
      hs.set_subtitles(wks, plot[ip], case_name[c], f'{lat1}:{lat2}N', var_str, font_height=0.015)

      #-------------------------------------------------------------------------

      # lres = hs.res_xy()
      # lres.xyLineThicknessF = 1
      # # lres.xyLineColor = 'gray'
      # lres.xyDashPattern = 2
      # xx = np.array([-1e8,1e8])
      # for k in [50]:
      #    yy = np.array([1,1]) * k
      #    ngl.overlay(plot[ip], ngl.xy(wks, xx, yy, lres) )

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

layout = [num_case+1,num_var] if add_obs else [num_case,num_var]

pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent = 5

# if use_common_label_bar: 
#    pnl_res.nglPanelLabelBar   = True
#    # pnl_res.lbTopMarginF       =  0.2
#    # pnl_res.lbBottomMarginF    = -0.2
#    pnl_res.lbLeftMarginF      = 0.5+0.5
#    pnl_res.lbRightMarginF     = 0.5
#    # pnl_res.lbTitleString      = "m/s"
#    # pnl_res.lbTitlePosition    = "bottom"
#    # pnl_res.lbLabelFontHeightF = 0.001
#    # pnl_res.lbTitleFontHeightF = 0.01

### add panel labels
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.01

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
