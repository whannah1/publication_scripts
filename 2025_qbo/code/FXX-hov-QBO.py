import os, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string, dask, pandas as pd
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
import cmocean
cscratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl'
gscratch = '/global/cfs/cdirs/m3312/whannah/2022-QBO-TEST'
pscratch = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu'
#-------------------------------------------------------------------------------
case_dir,case_sub = [],[]
case,case_name = [],[]
clr,dsh,mrk = [],[],[]
def add_case(case_in,n=None,p=None,s=None,d=0,c='black',m=0):
   global case_name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); case_name.append(n)
   case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
var,lev_list,vclr,vdsh = [],[],[],[]
def add_var(var_name,lev=-1,c='black',d=0): 
   var.append(var_name); lev_list.append(lev)
   vclr.append(c); vdsh.append(d)
##------------------------------------------------------------------------------
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01',       n='E3SM control',     c='red'  ,p=gscratch,s='run')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01', n='E3SM L72 smoothed',c='green',p=gscratch,s='run')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01', n='E3SM L80 refined', c='blue' ,p=scratch,s='run')

# add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72.GW-MOD-0', n='E3SM L72',      p=pscratch,s='run')
# add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72.GW-MOD-1', n='E3SM L72 GW-mod',   p=pscratch,s='run')
# add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-nsu40',    n='E3SM L72-nsu40',p=pscratch,s='run')
# add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-rscl',     n='E3SM L72-rscl', p=pscratch,s='run')
# add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-rlim',     n='E3SM L72-rlim', p=pscratch,s='run')

# add_case('E3SM.QBO-TEST-03.F2010.ne30pg2.L72',       n='E3SM L72',      p=pscratch,s='run')
# add_case('E3SM.QBO-TEST-03.F2010.ne30pg2.L72-nsu40', n='E3SM L72-smth', p=pscratch,s='run')
# add_case('E3SM.QBO-TEST-03.F2010.ne30pg2.L72-rlim',  n='E3SM L72-rlim', p=pscratch,s='run')
# add_case('E3SM.QBO-TEST-03.F2010.ne30pg2.L72-rscl',  n='E3SM L72-rscl', p=pscratch,s='run')

scratch = '/global/cfs/cdirs/m4310/whannah/E3SM'
add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L72' ,n='E3SM L72', p=scratch,s='archive/atm/hist')
add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L80' ,n='E3SM L80', p=scratch,s='archive/atm/hist')
add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L128',n='E3SM L128',p=scratch,s='run')

scratch = '/global/cfs/cdirs/m4310/wandiyu/E3SM/'
add_case('E3SM.2023-SCIDAC.ne30pg2_EC30to60E2r2.F20TR-MMF1.L64',n='MMF L64', p=scratch,s='archive/atm/hist')
add_case('E3SM.2023-SCIDAC.ne30pg2_EC30to60E2r2.F20TR-MMF1.L72',n='MMF L72', p=scratch,s='archive/atm/hist')


# remap_str,search_str = 'remap_90x180','h0.tem.'; first_file,num_files = 12*0,12*20

#-------------------------------------------------------------------------------
lev = np.array([   1.,    2.,    3.,    5.,    7.,   10.,   20.,   30.,   50.,   70.,  100.,  125.,  150.,])

var = ['U']

# add_var('dudt'     ,c='black',d=1)
# add_var('RES'      ,c='red')
# add_var('TEND'     ,c='blue')

# add_var('dudt'     )
# add_var('utendepfd')
# add_var('utendvtem')
# add_var('utendwtem')
# add_var('BUTGWSPEC')
# add_var('UTGWSPEC' )

#-------------------------------------------------------------------------------

fig_type = 'png'
fig_file = 'figs/FXX-hov-QBO'

lat1,lat2 = -5,5
# lat1,lat2 = -10,10

htype,num_files = 'h0',12*20

use_height_coord = False

print_stats = True
print_time  = False

var_x_case = False

use_common_label_bar = True

num_plot_col = 1

add_obs = False
obs_case = 'ERA5' # ERAi / ERA5

write_to_file = True

#-------------------------------------------------------------------------------
# Set up plot resources
#-------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)

if 'lev' not in vars(): lev = np.array([0])

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = []
res = hs.res_contour_fill()
res.vpHeightF = 0.1
res.tmYLLabelFontHeightF         = 0.01
res.tmXBLabelFontHeightF         = 0.01
res.tiXAxisFontHeightF           = 0.01
res.tiYAxisFontHeightF           = 0.01
res.lbLabelBarOn = False

# disable these by default - turn back on for bottom panel(s)
res.tmXBOn = False
res.tiXAxisOn = False

res.tiXAxisString = 'Time'
# if htype=='h0': res.tiXAxisString = 'Time [months]'
if htype=='h0': res.tiXAxisString = 'Time [years]'

if use_height_coord:
   res.tiYAxisString = 'Height [km]'
   res.nglYAxisType = "LinearAxis"
   res.trYMinF = 20
else:
   res.tiYAxisString = 'Pressure [hPa]'
   res.trYReverse = True
   # res.trYLog = True # doesn't work due to irregular spacing :(
   res.trYMinF = 5
   res.trYMaxF = 100

   tm_vals = [1,10,50,100,200]
   res.tmYLMode = 'Explicit'
   res.tmYLValues = tm_vals
   res.tmYLLabels = tm_vals
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
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

def write_file(case_in,var_in,data_avg):
   tmp_file = f'data/QBO.hov.v1.{case_in}.{var_in}.lat_{lat1}_{lat2}.nc'
   print('    writing to file: '+tmp_file)
   ds_out = xr.Dataset( coords=data_avg.coords )
   ds_out[var[v]] = data_avg
   ds_out.to_netcdf(path=tmp_file,mode='w')
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   tvar = var[v]
   area_name = 'area'
   #----------------------------------------------------------------------------
   # read the data
   #----------------------------------------------------------------------------
   data_list,time_list,lev_list = [],[],[]
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)

      data_dir_tmp,data_sub_tmp = None, None
      # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
      if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      if case_sub[c] is not None: data_sub_tmp = case_sub[c]

      case_obj = he.Case( name=case[c], atm_comp='eam', time_freq=None,
                          data_dir=data_dir_tmp, data_sub=data_sub_tmp )

      if 'lat1' in vars() : case_obj.lat1,case_obj.lat2 = lat1,lat2
      if 'lon1' in vars() : case_obj.lon1,case_obj.lon2 = lon1,lon2

      # avoid creating large chunks
      with dask.config.set(**{'array.slicing.split_large_chunks': True}):  

         data = case_obj.load_data(tvar,     htype=htype,num_files=num_files,lev=lev)
         area = case_obj.load_data(area_name,htype=htype,num_files=num_files).astype(np.double)

         # hc.print_stat(area,name='area',stat='naxsh',indent='    ',compact=True)
         
         # print(); print(data)
         # print(); print(data['lev'])
         # print(); print(area)
         # exit()

         if use_height_coord: Z = case_obj.load_data('Z3',htype=htype,num_files=num_files)

         hc.print_time_length(data.time,indent=' '*4)

         ### print stats before time averaging
         # if print_stats: hc.print_stat(data,name=var[v],stat='naxsh',indent='    ',compact=True)

         avg_dims= ('ncol')
         if 'lat' in data.dims: avg_dims= ('lat','lon')
         # data_avg = ( (data*area).sum(dim=avg_dims) / area.sum(dim=avg_dims) )
         data_avg = data.mean(dim=avg_dims)

         ### print stats after time averaging
         if print_stats: hc.print_stat(data_avg,name=var[v],stat='naxsh',indent='    ',compact=True)

         # print(); print(data)
         # print(); print(data_avg)

         if use_height_coord: Z_avg = ( (Z*area).sum(dim='ncol') / area.sum(dim='ncol') ).mean(dim='time') / 1e3
         #----------------------------------------------------------------------
         # adjust time to represent the middle of the month instead of the end
         time = data_avg.time.values
         time_orig = copy.deepcopy(time)
         for i,t in enumerate(time):
            if i==0:
               dt = pd.Timedelta('15 days')
            else:
               dt = ( time_orig[i] - time_orig[i-1] ) / 2
            time[i] = time_orig[i] - dt
         data_avg['time'] = time
         
         # Create simple time coordinate for plotting
         time = data_avg['time.year'].values + data_avg['time.month'].values/12.

         if print_time: print(); print(time); print()
         #----------------------------------------------------------------------
         # add data to lists
         data_list.append( data_avg.transpose().values )
         time_list.append( time )
         if     use_height_coord: lev_list.append( Z_avg.values )
         if not use_height_coord: lev_list.append( data_avg['lev'].values )
         #----------------------------------------------------------------------
         if write_to_file: write_file(case[c],var[v],data_avg)

   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   if var[v]=='U' and add_obs and not use_height_coord:
      print(f'    case: {obs_case}')
      # obs_scratch_path = '/lcrc/group/e3sm/diagnostics/observations/Atm/time-series'
      obs_scratch_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series'
      if obs_case=='ERAi': obs_data_file = f'{obs_scratch_path}/ERA-Interim/ua_197901_201612.nc'
      if obs_case=='ERA5': obs_data_file = f'{obs_scratch_path}/ERA5/ua_197901_201912.nc'
      ds = xr.open_dataset(obs_data_file)
      area = calculate_obs_area(ds['lon'].values,ds['lat'].values,ds['lon_bnds'].values,ds['lat_bnds'].values)
      area = xr.DataArray( area, coords=[ds['lat'],ds['lon']] )
      data = ds['ua']
      data = data.sel(lat=slice(lat1,lat2))
      area = area.sel(lat=slice(lat1,lat2))
      data_avg = (data*area).sum(dim=('lon','lat')) / area.sum(dim=('lon','lat'))


      ### truncate obs to match AMIP period
      # amip_yr1 = 1984
      # amip_yr2 = amip_yr1 + int(np.ceil(num_files/12))
      # yr = data_avg['time.year'].values
      # t_beg,t_end = None,None
      # for t,y in enumerate(yr):
      #    if t_beg is None and y==amip_yr1: t_beg = t
      #    if t_end is None and y> amip_yr2: t_end = t-1
      # if t_beg is not None and t_end is not None:
      #    # print(f'amip_yr1: {amip_yr1}  t_beg: {t_beg}')
      #    # print(f'amip_yr2: {amip_yr2}  t_end: {t_end}')
      #    data_avg = data_avg.isel(time=slice(t_beg,t_end))
      # else:
      #    exit(f'Something went wrong? t_beg: {t_beg}  t_end: {t_end}')

      time = data_avg['time.year'].values + data_avg['time.month'].values/12.
      sim_t1 = time_list[0][ 0]
      sim_t2 = time_list[0][-1]
      t_beg,t_end = None,None
      for i,t in enumerate(time):
         if t_beg is None and t==sim_t1: t_beg = i
         if t_end is None and t==sim_t2: t_end = i
      if t_beg is not None and t_end is not None:
         # print(f'\n  sim_t1/sim_t2: {sim_t1} / {sim_t2}  t_beg/t_end: {t_beg} / {t_end} \n')
         data_avg = data_avg.isel(time=slice(t_beg,t_end))
         time = data_avg['time.year'].values + data_avg['time.month'].values/12.
      else:
         exit(f'Something went wrong? t_beg: {t_beg}  t_end: {t_end}')

      if print_time: print(); print(time); print()

      k_list = []
      for k,p in enumerate(data_avg['plev'].values/1e2):
         if p in lev: k_list.append(k)
      data_avg = data_avg.isel(plev=k_list)

      
      lev_list.insert(0, data_avg['plev'].values/1e2 )
      time_list.insert(0, time )
      data_list.insert(0, data_avg.transpose().values )
      case_name.insert(0, obs_case )
   
      if write_to_file: write_file(case_name[0],var[v],data_avg)

   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)

   # tres.cnFillPalette = np.array( cmocean.cm.diff(np.linspace(0,1,256)) )
   # tres.cnFillPalette = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )
   tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )

   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])

   tres.cnLevelSelectionMode = 'ExplicitLevels'
   tres.cnLevels = np.linspace(-50,50,num=21)

   # nlev = 21
   # aboutZero = False
   # if var[v] in ['U','V'] : aboutZero = True
   # clev_tup = ngl.nice_cntr_levels(data_min, data_max, cint=None, max_steps=nlev, returnLevels=False, aboutZero=aboutZero )
   # if clev_tup==None: 
   #    tres.cnLevelSelectionMode = 'AutomaticLevels'   
   # else:
   #    cmin,cmax,cint = clev_tup
   #    tres.cnLevels = np.linspace(cmin,cmax,num=nlev)
   #    tres.cnLevelSelectionMode = 'ExplicitLevels'

   num_plot = (num_case+1) if add_obs else num_case
   for c in range(num_case):
      
      time = time_list[c]
      # time = ( time - time[0] ).astype('float') / 86400e9
      # time = ( time - time[0] ).astype('float') / (60*60*24*365)

      tres.sfYArray = lev_list[c]#.astype(int)
      tres.sfXArray = time
      # tres.sfXArray = np.linspace( 1./12., float(num_files)/12., num=len(time) )

      # print(); print(tres.sfXArray)

      if c==num_case-1: 
         tres.tmXBOn = True
         tres.tiXAxisOn = True


      plot.append( ngl.contour(wks, data_list[c], tres) )

      #-------------------------------------------------------------------------
      # Set strings at top of plot
      #-------------------------------------------------------------------------
      var_str = var[v]
      if var[v]=='PRECT' : var_str = 'Precipitation [mm/day]'
      if var[v]=='U'     : var_str = 'Zonal Wind [m/s]'

      ctr_str = ''
      # if var[v] in ['PRECT','PRECC','PRECL'] : ctr_str = 'Mean: '+'%.2f'%avg_X+' [mm/day]'
      hs.set_subtitles(wks, plot[len(plot)-1], case_name[c], ctr_str, var_str, font_height=0.015)

#-------------------------------------------------------------------------------
# Finalize plot
#-------------------------------------------------------------------------------

# layout = [len(plot),1]
# layout = [num_var,num_case] if var_x_case else [num_case,num_var]


layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent = 5
if use_common_label_bar: 
   pnl_res.nglPanelLabelBar   = True
   # pnl_res.lbTopMarginF       =  0.2
   # pnl_res.lbBottomMarginF    = -0.2
   pnl_res.lbLeftMarginF      = 0.5+0.5
   pnl_res.lbRightMarginF     = 0.5
   # pnl_res.lbTitleString      = "m/s"
   # pnl_res.lbTitlePosition    = "bottom"
   # pnl_res.lbLabelFontHeightF = 0.001
   # pnl_res.lbTitleFontHeightF = 0.01

### add panel labels
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.01

# if num_var==1  : layout = [num_case,num_var]
# if num_case==1 : layout = [num_var,num_case]

ngl.panel(wks,plot[0:len(plot)],layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
