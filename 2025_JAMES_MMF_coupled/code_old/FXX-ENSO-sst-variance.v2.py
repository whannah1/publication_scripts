import os, copy, string, ngl, xarray as xr, numpy as np, cmocean, glob
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
var,lev_list = [],[]
def add_var(var_name,lev=-1): 
   var.append(var_name); lev_list.append(lev)
#---------------------------------------------------------------------------------------------------
tmp_path_hst_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_v2  = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
scrip_file_path  = os.getenv('HOME')+'/E3SM/data_grid/ne30pg2_scrip.nc'
add_case('HadSST',                                                 n='HadSST',  c='black')
add_case('v2.LR.historical_0101',                                  n='E3SMv2',  c='cyan',    p=tmp_path_hst_v2, s='archive/atm/hist')
add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',c='blue', p=tmp_path_hst_mmf,s='archive/atm/hist')
#---------------------------------------------------------------------------------------------------

htype,yr1,yr2 = 'h0',1950,2014

fig_file,fig_type = 'figs/FXX-ENSO-sst-variance','png'

write_file    = False
print_stats   = True
overlay_cases = False

recalculate = False

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
npix=1024*2; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
if 'legend_file' in locals(): lgd_wks = ngl.open_wks('png',legend_file,wkres)

res = hs.res_contour_fill_map()
res.tmYLLabelFontHeightF  = 0.008
res.tmXBLabelFontHeightF  = 0.008
# res.lbLabelFontHeightF    = 0.015
res.lbLabelBarOn          = False
res.tmXBOn                = False
res.tmYLOn                = False


# res.cnFillPalette = "MPL_viridis"
# res.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
res.cnFillPalette = np.array( cmocean.cm.amp(np.linspace(0,1,256)) )
res.cnLevelSelectionMode = "ExplicitLevels"

# res.cnLevels = np.linspace(0,3,11)
# res.cnLevels = np.arange(0.4,3.6+0.4,0.4)
res.cnLevels = np.arange(0.3,3.0+0.3,0.3)

res.mpMinLatF = -25
res.mpMaxLatF =  25
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
scrip_ds = xr.open_mfdataset(scrip_file_path).rename({'grid_size':'ncol'})
area = scrip_ds['grid_area']
lat  = scrip_ds['grid_center_lat']
lon  = scrip_ds['grid_center_lon']

landfrac_ds = xr.open_dataset('/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4pg2_16xdel2.c20200108.nc')
landfrac = landfrac_ds['LANDFRAC']
#---------------------------------------------------------------------------------------------------

time_list,data_list = [],[]
for c in range(num_case):
   print('    case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)
   tmp_file = os.getenv('HOME')+f'/Research/E3SM/data_temp/ENSO.sst-variance.{case[c]}.nc'
   #----------------------------------------------------------------------------
   if recalculate:
      if case[c]=='HadSST':
         file_name = '/global/cfs/cdirs/m3312/whannah/obs_data/HadSST/HadISST_sst.remap_ne30pg2.nc'
         ds = xr.open_dataset( file_name )
         data = ds['sst']
      else:
         file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/*.eam.{htype}.*'
         file_list = sorted(glob.glob(file_path))
         ds = xr.open_mfdataset( file_list )
         data = ds['TS']
      #-------------------------------------------------------------------------
      # subset in time
      ds = ds.where( ds['time.year']>=yr1, drop=True)
      ds = ds.where( ds['time.year']<=yr2, drop=True)
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
      # calculate variance
      # data = data.var(dim='time')
      data = data.std(dim='time')
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
   
   data = data.where(landfrac<0.5,np.nan)

   data_list.append( data.values )


#-------------------------------------------------------------------------------
# Create plot
#-------------------------------------------------------------------------------
# data_min = np.min([np.nanmin(d) for d in data_list])
# data_max = np.max([np.nanmax(d) for d in data_list])

tmp_num_case = num_case

pdum = [None]*tmp_num_case*3

for c in range(tmp_num_case):

   tres = copy.deepcopy(res)
   
   tres.cnFillMode       = "AreaFill"
   tres.sfXArray         = scrip_ds['grid_center_lon'].values
   tres.sfYArray         = scrip_ds['grid_center_lat'].values

   # tres.cnFillMode       = "CellFill"
   # tres.sfXArray         = scrip_ds['grid_center_lon'].values
   # tres.sfYArray         = scrip_ds['grid_center_lat'].values
   # tres.sfXCellBounds    = scrip_ds['grid_corner_lon'].values
   # tres.sfYCellBounds    = scrip_ds['grid_corner_lat'].values

   ip = c

   plot[ip] = ngl.contour_map(wks,np.ma.masked_invalid(data_list[c]),tres)

   #----------------------------------------------------------------------------
   ### draw boxes around ENSO regions
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

   #----------------------------------------------------------------------------

   # hs.set_subtitles(wks, plot[ip], case_name[c], '', f'SST Variance [K~S~2~N~]', font_height=0.02)
   hs.set_subtitles(wks, plot[ip], case_name[c], '', f'SST Std Deviation [K]', font_height=0.02)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
hc.printline()

layout = [num_case,1]

pnl_res = hs.setres_panel()
pnl_res.nglPanelLabelBar = True
# pnl_res.lbLeftMarginF      = -0.3
# pnl_res.lbRightMarginF     = -0.3
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.015

ngl.panel(wks,plot[0:len(plot)],layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
if 'legend_file' in locals(): hc.trim_png(legend_file)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
