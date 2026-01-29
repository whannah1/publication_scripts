import os, copy, glob, ngl, xarray as xr, numpy as np, warnings, string
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

var,lev_list = [],[]
def add_var(var_name,lev=-1,mask=None): 
   var.append(var_name); lev_list.append(lev)
#---------------------------------------------------------------------------------------------------
tmp_path_co2_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_v2  = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
add_case('BEST',                                                 n='BEST',  c='black')
# add_case('HadSST',                                                 n='HadSST',  c='black')
add_case('v2.LR.historical_0101',                                  n='E3SMv2',  p=tmp_path_hst_v2, s='archive/atm/hist')
add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',p=tmp_path_hst_mmf,s='archive/atm/hist')
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2',n='MMF PI 1xCO2',c='blue' ,p=tmp_path_co2_mmf,s='archive/atm/hist')
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.2xCO2',n='MMF PI 2xCO2',c='green',p=tmp_path_co2_mmf,s='archive/atm/hist')
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.4xCO2',n='MMF PI 4xCO2',c='red'  ,p=tmp_path_co2_mmf,s='archive/atm/hist')
#---------------------------------------------------------------------------------------------------

htype,yr1,yr2 = 'h0',1950,2014


# region,lat1,lat2,lon1,lon2 = 'Nino3'  ,-5,5,360-150,360-90
# region,lat1,lat2,lon1,lon2 = 'Nino4'  ,-5,5,160,360-150
region,lat1,lat2,lon1,lon2 = 'Nino3.4',-5,5,190,240

fig_file,fig_type = 'figs/F09-ENSO-timeseries','png'
tmp_file_head = 'data/ENSO-timeseries'


recalculate   = False

print_stats   = True

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
num_case = len(case)


if 'clr' not in vars() or clr==[]: clr = ['black']*num_case
if 'dsh' not in vars() or dsh==[]: dsh = np.zeros(num_case)

if 'lev' not in vars(): lev = np.array([0])

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
plot = [None]*(num_case)

wkres = ngl.Resources()
npix=1024; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
if 'legend_file' in locals(): lgd_wks = ngl.open_wks('png',legend_file,wkres)

res = hs.res_xy()
res.vpHeightF = 0.1
res.tmYLLabelFontHeightF   = 0.01
res.tmXBLabelFontHeightF   = 0.01
res.tiXAxisFontHeightF     = 0.01
res.tiYAxisFontHeightF     = 0.01
res.xyLineThicknessF       = 4

lres = hs.res_xy()
lres.xyDashPattern    = 0
lres.xyLineThicknessF = 1
lres.xyLineColor      = "black"

#---------------------------------------------------------------------------------------------------
scrip_file_path  = os.getenv('HOME')+'/E3SM/data_grid/ne30pg2_scrip.nc'
scrip_ds = xr.open_mfdataset(scrip_file_path).rename({'grid_size':'ncol'})
area = scrip_ds['grid_area']
lat  = scrip_ds['grid_center_lat']
lon  = scrip_ds['grid_center_lon']

mask = xr.DataArray( np.ones(len(area),dtype=bool), coords={'ncol':scrip_ds['ncol']} )
mask = mask & (lat>=lat1) & (lat<=lat2)
mask = mask & (lon>=lon1) & (lon<=lon2)
mask.load()

area = area.where( mask, drop=True)
#---------------------------------------------------------------------------------------------------

time_list,data_list = [],[]
for c in range(num_case):
   print('    case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)
   tmp_file = f'{tmp_file_head}.{case[c]}.reg_{region}.lat_{lat1}_{lat2}.lon_{lon1}_{lon2}.nc'
   #----------------------------------------------------------------------------
   if recalculate:
      if case[c]=='HadSST':
         file_name = '/global/cfs/cdirs/m3312/whannah/obs_data/HadSST/HadISST_sst.remap_ne30pg2.nc'
         ds = xr.open_dataset( file_name )
         data = ds['sst']
      elif case[c]=='BEST':
         obs_file = f'/global/cfs/cdirs/m3312/whannah/obs_data/BEST/Land_and_Ocean_LatLong1.remap_ne30pg2.nc'
         ds = xr.open_dataset(obs_file)
         data = ds['temperature']
         file_time_list = [None]*len(data['time'])
         # convert anomalies to absolute temperature and create time coord
         for t,tt in enumerate(data['time'].values):
            yr_val = int(np.floor(tt))
            mn_ind = int(np.floor((tt-yr_val)*12))
            file_time_list[t] = np.datetime64(f'{yr_val}-{(mn_ind+1):02d}')
            data[t,:] = data[t,:] + ds['climatology'].isel(month_number=mn_ind)
         with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            time = xr.DataArray(file_time_list,dims=('time'))
         data = data.assign_coords(time=time)
      else:
         file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/*.eam.{htype}.*'
         file_list = sorted(glob.glob(file_path))
         ds = xr.open_mfdataset( file_list )
         data = ds['TS']
      #-------------------------------------------------------------------------
      # subset in time
      data = data.where( data['time.year']>=yr1, drop=True)
      data = data.where( data['time.year']<=yr2, drop=True)
      #----------------------------------------------------------------------
      # regional subset
      data = data.where( mask, drop=True)
      #-------------------------------------------------------------------------
      # remove annual cycle
      for n in range(12):
         month_index = slice(n,len(data.time),12)
         data[month_index,:] = data.isel(time=month_index) - data.isel(time=month_index).mean(dim='time')
      #-------------------------------------------------------------------------
      # convert to anomalies
      data = data - data.mean(dim='time')
      # detrend in time
      data = data - xr.polyval(data['time'], data.polyfit(dim='time', deg=1).polyfit_coefficients)
      # calculate spatial mean
      data = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )
      #-------------------------------------------------------------------------
      # Write to file 
      #-------------------------------------------------------------------------
      if os.path.isfile(tmp_file) : os.remove(tmp_file)
      tmp_ds = xr.Dataset( coords=data.coords )
      tmp_ds['TS'] = data
      tmp_ds.to_netcdf(path=tmp_file,mode='w')
   else:
      tmp_ds = xr.open_dataset( tmp_file )
      data = tmp_ds['TS']
   #----------------------------------------------------------------------------
   if print_stats: hc.print_stat(data,name='',stat='naxs',indent=' '*6,compact=True)
   
   plot_time = data['time.year'].values + data['time.month'].values/12.

   data_list.append( data.values )
   time_list.append( plot_time )
   
#-------------------------------------------------------------------------------
# Create plot

tres = copy.deepcopy(res)
tres.tiXAxisString = 'Year'

# Make sure plot bounds are consistent
tres.trYMinF = -4 #np.min([np.nanmin(d) for d in data_list])
tres.trYMaxF =  4 #np.max([np.nanmax(d) for d in data_list])
tres.trXMinF = yr1#np.min([np.nanmin(d) for d in time_list])
tres.trXMaxF = yr2#np.max([np.nanmax(d) for d in time_list])

tmp_num_case = num_case

for c in range(tmp_num_case):
   
   ip = c
   
   tres.xyLineColor   = clr[c]
   tres.xyDashPattern = dsh[c]

   tplot = ngl.xy(wks, time_list[c], data_list[c], tres)
   
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

   #----------------------------------------------------------------------------

   hs.set_subtitles(wks, plot[ip], case_name[c], region, '', font_height=0.02)

#-------------------------------------------------------------------------------
# Add Legend
# lgres = ngl.Resources()
# lgres.vpWidthF, lgres.vpHeightF  = 0.08, 0.08  
# lgres.lgLabelFontHeightF = 0.01
# lgres.lgLineThicknessF   = 4
# lgres.lgMonoLineColor    = False
# # lgres.lgMonoDashIndex    = True
# lgres.lgLineColors       = clr
# lgres.lgDashIndexes      = dsh
# lgres.lgLabelJust    = 'CenterLeft'
# # pid = ngl.legend_ndc(wks, len(name), name, 0.5, 0.4, lgres)  # 3x2
# # pid = ngl.legend_ndc(wks, len(name), name, 0.5, 0.1, lgres)  # 3x2
# # pid = ngl.legend_ndc(wks, len(name), name, 0.3, 0.5, lgres)  # 1x2

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
hc.printline()

layout = [num_case,1]

pnl_res = hs.setres_panel()
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.018

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
if 'legend_file' in locals(): hc.trim_png(legend_file)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
