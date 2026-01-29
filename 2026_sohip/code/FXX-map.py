import os, glob, ngl, subprocess as sp, numpy as np, xarray as xr
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
import copy, string
import cmocean
# export PYNGL_RANGS=~/.conda/envs/pyn_env/lib/ncarg/database/rangs
host = hc.get_host()
case_root = '/pscratch/sd/w/whannah/scream_scratch/pm-gpu'
grid_root = '/global/cfs/cdirs/m4842/whannah/files_grid'
#-------------------------------------------------------------------------------
case_opts_list = []
case_list = []
def add_case(case_in,**kwargs):
   case_list.append(case_in)
   tmp_opts = {}
   for k, val in kwargs.items(): tmp_opts[k] = val
   if '256x2-eq-ind-v1' in case_in: tmp_opts['g'] = f'{grid_root}/2025-sohip-256x2-eq-ind-v1-pg2_scrip.nc'
   if '256x2-eq-ind-v1' in case_in: tmp_opts['g'] = f'{grid_root}/2025-sohip-256x2-eq-ind-v1-pg2_scrip.nc'
   if '256x2-ptgnia-v1' in case_in: tmp_opts['g'] = f'{grid_root}/2025-sohip-256x2-ptgnia-v1-pg2_scrip.nc'
   if '256x2-sc-ind-v1' in case_in: tmp_opts['g'] = f'{grid_root}/2025-sohip-256x2-sc-ind-v1-pg2_scrip.nc'
   if '256x2-sc-pac-v1' in case_in: tmp_opts['g'] = f'{grid_root}/2025-sohip-256x2-sc-pac-v1-pg2_scrip.nc'
   if '256x2-se-pac-v1' in case_in: tmp_opts['g'] = f'{grid_root}/2025-sohip-256x2-se-pac-v1-pg2_scrip.nc'
   if '256x2-sw-ind-v1' in case_in: tmp_opts['g'] = f'{grid_root}/2025-sohip-256x2-sw-ind-v1-pg2_scrip.nc'
   if '256x3-eq-ind-v1' in case_in: tmp_opts['g'] = f'{grid_root}/2025-sohip-256x3-eq-ind-v1-pg2_scrip.nc'
   if '256x3-eq-ind-v1' in case_in: tmp_opts['g'] = f'{grid_root}/2025-sohip-256x3-eq-ind-v1-pg2_scrip.nc'
   if '256x3-ptgnia-v1' in case_in: tmp_opts['g'] = f'{grid_root}/2025-sohip-256x3-ptgnia-v1-pg2_scrip.nc'
   if '256x3-sc-ind-v1' in case_in: tmp_opts['g'] = f'{grid_root}/2025-sohip-256x3-sc-ind-v1-pg2_scrip.nc'
   if '256x3-sc-pac-v1' in case_in: tmp_opts['g'] = f'{grid_root}/2025-sohip-256x3-sc-pac-v1-pg2_scrip.nc'
   if '256x3-se-pac-v1' in case_in: tmp_opts['g'] = f'{grid_root}/2025-sohip-256x3-se-pac-v1-pg2_scrip.nc'
   if '256x3-sw-ind-v1' in case_in: tmp_opts['g'] = f'{grid_root}/2025-sohip-256x3-sw-ind-v1-pg2_scrip.nc'
   tmp_opts['p'] = case_root
   tmp_opts['s'] = 'run'
   case_opts_list.append(tmp_opts)
#-------------------------------------------------------------------------------
var,var_str,klev_list = [],[],[]
unit_fac_list = []
def add_var(var_name,vstr=None,klev=0,unit_fac=None): 
   var.append(var_name);
   var_str.append(var_name if vstr is None else vstr)
   klev_list.append(klev)
   unit_fac_list.append(unit_fac)
#-------------------------------------------------------------------------------

# add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-19.09.NN_420',n='256x2-eq-ind-v1',xlat=0, xlon=  90,tlat= -6.99,tlon=  84.74,slat= 10.24,slon=  94.16)
# add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-21.02.NN_420',n='256x2-eq-ind-v1',xlat=-5, xlon=  80,tlat= -3.05,tlon=  75.97,slat= 13.56,slon=  85.35)
# add_case('2025-SOHIP-RRM-00.256x2-ptgnia-v1.2023-06-13.19',       n='256x2-ptgnia-v1',xlat=-60,xlon= -50,tlat=-49.46,tlon= -60.24,slat=  None,slon=   None)
add_case('2025-SOHIP-RRM-00.256x2-sc-ind-v1.2023-06-21.09',       n='256x2-sc-ind-v1',xlat=-50,xlon=  80,tlat=-52.49,tlon=  67.04,slat=-51.03,slon=  98.64)
# add_case('2025-SOHIP-RRM-00.256x2-sc-pac-v1.2023-06-14.15',       n='256x2-sc-pac-v1',xlat=-35,xlon=-135,tlat=-34.73,tlon=-136.73,slat=-43.76,slon=-114.47)
# add_case('2025-SOHIP-RRM-00.256x2-se-pac-v1.2023-06-12.16',       n='256x2-se-pac-v1',xlat=-50,xlon= -95,tlat=-49.60,tlon= -94.45,slat=-51.80,slon= -63.70)
# add_case('2025-SOHIP-RRM-00.256x2-sw-ind-v1.2023-06-12.06',       n='256x2-sw-ind-v1',xlat=-50,xlon=  45,tlat=-49.61,tlon=  45.20,slat=-51.79,slon=  75.97)

# add_case('2025-SOHIP-RRM-00.256x3-eq-ind-v1.2023-06-19.09',       n='256x3-eq-ind-v1',xlat=-5, xlon=  80,tlat= -6.99,tlon=  84.74,slat= 10.24,slon=  94.16)
# add_case('2025-SOHIP-RRM-00.256x3-eq-ind-v1.2023-06-21.02',       n='256x3-eq-ind-v1',xlat=-5, xlon=  80,tlat= -3.05,tlon=  75.97,slat= 13.56,slon=  85.35)
# add_case('2025-SOHIP-RRM-00.256x3-ptgnia-v1.2023-06-13.19',       n='256x3-ptgnia-v1',xlat=-60,xlon= -50,tlat=-49.46,tlon= -60.24,slat=  None,slon=   None)
# add_case('2025-SOHIP-RRM-00.256x3-sc-ind-v1.2023-06-21.09',       n='256x3-sc-ind-v1',xlat=-50,xlon=  80,tlat=-52.49,tlon=  67.04,slat=-51.03,slon=  98.64)
# add_case('2025-SOHIP-RRM-00.256x3-sc-pac-v1.2023-06-14.15',       n='256x3-sc-pac-v1',xlat=-35,xlon=-135,tlat=-34.73,tlon=-136.73,slat=-43.76,slon=-114.47)
# add_case('2025-SOHIP-RRM-00.256x3-se-pac-v1.2023-06-12.16',       n='256x3-se-pac-v1',xlat=-50,xlon= -95,tlat=-49.60,tlon= -94.45,slat=-51.80,slon= -63.70)
# add_case('2025-SOHIP-RRM-00.256x3-sw-ind-v1.2023-06-12.06',       n='256x3-sw-ind-v1',xlat=-50,xlon=  45,tlat=-49.61,tlon=  45.20,slat=-51.79,slon=  75.97)


# view width/height in degrees
dx = 30
dy = 20

htype = 'output.scream.2D.10min.INSTANT.nmins_x10'
# first_file,num_files = -1,1

#-------------------------------------------------------------------------------

# add_var('ps')
add_var('precip_total_surf_mass_flux', vstr='precip',unit_fac=86400*1e3)
# add_var('VapWaterPath')
# add_var('LiqWaterPath')
# add_var('IceWaterPath')
# add_var('surf_sens_flux')
# add_var('surface_upward_latent_heat_flux')
# add_var('wind_speed_10m')
# add_var('SW_flux_up_at_model_top')
# add_var('SW_flux_dn_at_model_top')
# add_var('LW_flux_up_at_model_top')

fig_file,fig_type = os.getenv('HOME')+'/Research/E3SM/pub_figs/sohip/figs/FXX.map','png'

#-------------------------------------------------------------------------------
# lat1,lat2 = -40,40
# lat1,lat2 = -60,60

# lat1,lat2,lon1,lon2 = 0,30,360-65,360-15 # Atlantic wide view

# lat1,lat2,lon1,lon2 = 0,40,180-60,180+0 # Pacific wide view
# lat1,lat2,lon1,lon2 = -40,40,180-60,180+60 # Pacific wider view
# lat1,lat2,lon1,lon2 = -50,50,180-80,180+80 # Pacific wider-er view

# lat1,lat2,lon1,lon2 = -20,0,220,240
# lat1,lat2,lon1,lon2 = 10,30,360-100,360-55

# lat1,lat2,lon1,lon2 = -10,70,180,340             # N. America
# lat1,lat2,lon1,lon2 = 25, 50, 360-125, 360-75    # CONUS
# lat1,lat2,lon1,lon2 = -40,40,90,240              # MC + West Pac
# lat1,lat2,lon1,lon2 = -20,20,150,200              # MC + West Pac
# lat1,lat2,lon1,lon2 = -5,5,150,160              # MC + West Pac
# lat1,lat2,lon1,lon2 = 0,60,50,120                # India
# lat1,lat2,lon1,lon2 =  20,60,360-60,360-10              # Atlantic
# lat1,lat2,lon1,lon2 = -45,-45+60,360-120,360-40  # S. America
# lat1,lat2,lon1,lon2 = -0,25,360-115,360-75  # Panama

# lat1,lat2,lon1,lon2 = 15,35,360-100,360-75 # Caribbean
# lat1,lat2,lon1,lon2 = 27,30,360-92,360-88 # Louisana coast (Katrina landfall)
# lat1,lat2,lon1,lon2 = 25,30,360-92,360-86 # Louisana coast wider (Katrina landfall)

# lat1,lat2,lon1,lon2 =10,40,100,140  # Taiwan

# lat1,lat2,lon1,lon2 = 0,10,0,10
#-------------------------------------------------------------------------------
use_remap,remap_grid = False,'90x180' # 90x180 / 180x360

use_snapshot,ss_t = False,-1

print_stats = True

var_x_case = False

num_plot_col = int(np.sqrt(len(case_list)))

use_common_label_bar = True

if use_snapshot: print(); print(f'{hc.tcolor.RED}WARNING - snapshot mode enabled! (ss_t={ss_t}){hc.tcolor.ENDC}'); print()

#---------------------------------------------------------------------------------------------------
# Set up plot resources

if case_list==[]: raise ValueError('ERROR - case list is empty!')
num_var,num_case = len(var),len(case_list)

subtitle_font_height = 0.01

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix

wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*(num_var*num_case)
   
res = hs.res_contour_fill_map()
res.tmYLLabelFontHeightF         = 0.008
res.tmXBLabelFontHeightF         = 0.008
res.lbLabelFontHeightF           = 0.01
# res.tmXBOn                       = False
# res.tmYLOn                       = False
# res.mpGeophysicalLineColor       = 'white'
# res.mpCenterLonF                 = 180
# res.mpProjection                 = 'Robinson'
# res.pmTickMarkDisplayMode        = "Never"

res.mpDataBaseVersion = 'MediumRes'
# res.mpDataBaseVersion = 'HighRes'

mres = ngl.Resources()
mres.nglDraw,mres.nglFrame         = False,False
mres.xyMarkLineMode = 'Markers'
mres.xyMarkerColor = 'red'
mres.xyMarkerSizeF = 0.008

#---------------------------------------------------------------------------------------------------
def get_comp(case):
   comp = 'scream'
   return comp
#---------------------------------------------------------------------------------------------------
def get_ctr_str(glb_avg=None): return ''
# def get_ctr_str(glb_avg=None):
#    ctr_str = ''
#    if 'lat1' in globals():
#       lat1_str = f'{lat1}N' if lat1>=0 else f'{(lat1*-1)}S'
#       lat2_str = f'{lat2}N' if lat2>=0 else f'{(lat2*-1)}S'
#       ctr_str += f' {lat1_str}:{lat2_str} '
#    if 'lon1' in globals():
#       lon1_str = f'{lon1}E' #if lon1>=0 and lon1<=360 else f'{(lon1*-1)}S'
#       lon2_str = f'{lon2}E' #if lon2>=0 and lon2<=360 else f'{(lon2*-1)}S'
#       ctr_str += f' {lon1_str}:{lon2_str} '
#    # if glb_avg is not None:
#    #    # add logic here t display global average value?
#    return ctr_str
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   data_list,area_list,lat_list,lon_list = [],[],[],[]
   glb_avg_list = []
   std_list,cnt_list = [],[]
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case_list[c]+hc.tcolor.ENDC)
      case_opts = case_opts_list[c]
      case_name = case_opts['n']
      case_root = case_opts['p']
      case_sub  = case_opts['s']
      #-------------------------------------------------------------------------
      # read the data
      file_path = f'{case_root}/{case_list[c]}/{case_sub}/*{htype}*'
      
      file_list = sorted(glob.glob(file_path))
      if 'first_file' in locals(): file_list = file_list[first_file:]
      if 'num_files'  in locals(): file_list = file_list[:num_files]

      for f in file_list: print(f'    {hc.tcolor.YELLOW}{f}{hc.tcolor.ENDC}')

      if file_list==[]:
         print(); print('ERROR - file_list is empty!')
         print(); print(f'file_path: {file_path}')
         print()
      #-------------------------------------------------------------------------
      ds = xr.open_mfdataset( file_list )
      area = ds['area']
      #-------------------------------------------------------------------------
      tvar = var[v]
      
      # var_chk = var[v] in ds.variables
      # if var[v]=='precip_total_surf_mass_flux' and not var_chk: tvar = 'precip_liq_surf_mass_flux'
      # if var[v]=='horiz_winds_at_model_bot_u': tvar = 'horiz_winds_at_model_bot'
      # if var[v]=='horiz_winds_at_model_bot_v': tvar = 'horiz_winds_at_model_bot'
      # if var[v]=='RESTOM': tvar = 'SW_flux_dn_at_model_top'

      data = ds[tvar]
      
      # if var[v]=='precip_total_surf_mass_flux' and not var_chk: data = data + ds['precip_ice_surf_mass_flux']
      # if var[v]=='horiz_winds_at_model_bot_u': data = data.isel(dim2=0)
      # if var[v]=='horiz_winds_at_model_bot_v': data = data.isel(dim2=1)
      # if var[v]=='RESTOM':
      #    data = data - ds['SW_flux_up_at_model_top']
      #    data = data - ds['LW_flux_up_at_model_top']

      #-------------------------------------------------------------------------
      # # Get rid of lev dimension
      # if 'lev' in data.dims:
      #    if klev_list[v] is not None:
      #       data = data.isel(lev=klev_list[v])
      #-------------------------------------------------------------------------
      # adjust units
      if unit_fac_list[v] is not None: data = data * unit_fac_list[v]
      #-------------------------------------------------------------------------
      # print stats before time averaging
      if print_stats:
         hc.print_stat(data,name=var[v],stat='naxsh',indent='    ',compact=True)
      #-------------------------------------------------------------------------
      # average over time dimension
      if 'time' in data.dims : 
         if use_snapshot:
            data = data.isel(time=ss_t)
         else:
            data = data.mean(dim='time')
      #-------------------------------------------------------------------------
      # Calculate area weighted global mean
      if 'area' in locals() :
         gbl_mean = ( (data*area).sum() / area.sum() ).values 
         print(hc.tcolor.CYAN+f'      Area Weighted Global Mean : {gbl_mean:6.4}'+hc.tcolor.ENDC)
         glb_avg_list.append(gbl_mean)
      #-------------------------------------------------------------------------
      # append to data lists

      # if case_list[c]=='TRMM' and 'lon1' not in locals(): 
      #    data_list.append( ngl.add_cyclic(data.values) )
      # else:
      #    data_list.append( data.values )

      data_list.append( data.values )

      if 'area' in locals() : area_list.append( area.values )
   #------------------------------------------------------------------------------------------------
   # Plot averaged data
   #------------------------------------------------------------------------------------------------
   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])

   for c in range(num_case):
      #-------------------------------------------------------------------------
      case_opts = case_opts_list[c]
      ctr_lon   = case_opts['xlon']
      ctr_lat   = case_opts['xlat']
      grid_file = case_opts['g']
      grid_ds = xr.open_dataset(grid_file)
      grid_center_lon = grid_ds['grid_center_lon']
      grid_center_lat = grid_ds['grid_center_lat']
      grid_corner_lon = grid_ds['grid_corner_lon']
      grid_corner_lat = grid_ds['grid_corner_lat']
      #-------------------------------------------------------------------------
      # Set colors
      tres = copy.deepcopy(res)
      rain_clr_map = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
      amp_clr_map  = np.array( cmocean.cm.amp(np.linspace(0,1,256)) )
      bal_clr_map  = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
      
      tres.cnFillPalette = 'MPL_viridis'

      # if var[v]=='U'                          : tres.cnFillPalette = bal_clr_map
      # if 'horiz_winds' in var[v]              : tres.cnFillPalette = bal_clr_map
      # if var[v]=='T_2m'                       : tres.cnFillPalette = amp_clr_map
      # if var[v]=='precip_total_surf_mass_flux': tres.cnFillPalette = rain_clr_map
      # if var[v]=='precip_ice_surf_mass'       : tres.cnFillPalette = rain_clr_map
      # if var[v]=='precip_liq_surf_mass'       : tres.cnFillPalette = rain_clr_map
      # if var[v]=='LiqWaterPath'               : tres.cnFillPalette = rain_clr_map
      # if var[v]=='IceWaterPath'               : tres.cnFillPalette = rain_clr_map
      
      #-------------------------------------------------------------------------
      # Set explicit contour levels
      #-------------------------------------------------------------------------
      # if 'precip' in var[v]      : tres.cnLevels = np.arange(1,20+1,1)
      # if 'precip' in var[v]      : tres.cnLevels = np.arange(4,80+4,4)
      if 'precip' in var[v]      : tres.cnLevels = np.logspace( 0, 3, num=30).round(decimals=2)
      # if 'prec' in var[v]      : tres.cnLevels = np.logspace( -2, 2, num=20).round(decimals=2)

      if var[v]=='VapWaterPath'  : tres.cnLevels = np.arange(10,90+1,1)
      if var[v]=='LiqWaterPath'  : tres.cnLevels = np.arange(1,61+3,3)/1e2
      if var[v]=='IceWaterPath'  : tres.cnLevels = np.arange(1,61+3,3)/1e1

      # if 'precip' in var[v]      : tres.cnLevels = np.logspace(  -0.5, 2.6, num=30).round(decimals=2)
      # if 'precip' in var[v]      : tres.cnLevels = np.logspace(  -0.2, 2.6, num=40).round(decimals=2)
      if var[v]=='LiqWaterPath'  : tres.cnLevels = np.logspace( -2, 1.2, num=20).round(decimals=2)
      if var[v]=='IceWaterPath'  : tres.cnLevels = np.logspace( -2, 1.2, num=20).round(decimals=2)

      #-------------------------------------------------------------------------
      ### print color levels
      # if hasattr(tres,'cnLevels') : 
      #    print(f'\ntres.cnLevels:')
      #    msg = ''
      #    for cl in range(len(tres.cnLevels)):
      #       msg += f'{tres.cnLevels[cl]}, '
      #    print(f'[ {msg} ]\n')
      #-------------------------------------------------------------------------
      # set non-explicit contour levels
      if hasattr(tres,'cnLevels') : 
         tres.cnLevelSelectionMode = 'ExplicitLevels'
      else:
         nlev = 21
         aboutZero = False
         # if var[v] in []: aboutZero = True
         # if var[v]=='U': aboutZero = True
         clev_tup = ngl.nice_cntr_levels(data_min, data_max, cint=None, max_steps=nlev, \
                                         returnLevels=False, aboutZero=aboutZero )
         if clev_tup==None: 
            tres.cnLevelSelectionMode = 'AutomaticLevels'   
         else:
            cmin,cmax,cint = clev_tup
            tres.cnLevels = np.linspace(cmin,cmax,num=nlev)
            tres.cnLevelSelectionMode = 'ExplicitLevels'
      #-------------------------------------------------------------------------
      tres.mpMinLatF = max(-90,ctr_lat-dy)
      tres.mpMaxLatF = min( 90,ctr_lat+dy)
      tres.mpMinLonF = ctr_lon-dx
      tres.mpMaxLonF = ctr_lon+dx
      #-------------------------------------------------------------------------
      # Create plot
      #-------------------------------------------------------------------------
      tres.lbLabelBarOn = False if use_common_label_bar else True

      tres.cnFillMode    = 'CellFill'
      tres.sfXArray      = grid_center_lon.values
      tres.sfYArray      = grid_center_lat.values
      tres.sfXCellBounds = grid_corner_lon.values
      tres.sfYCellBounds = grid_corner_lat.values

      ip = v*num_case+c if var_x_case else c*num_var+v

      plot[ip] = ngl.contour_map(wks,np.ma.masked_invalid(data_list[c]),tres) 
      
      #-------------------------------------------------------------------------
      slat = case_opts['slat'] # ISS position
      slon = case_opts['slon']
      tlat = case_opts['tlat'] # tangent point of retreival
      tlon = case_opts['tlon']
      #-------------------------------------------------------------------------
      if slat is not None and slon is not None:
         mres.xyMarker = 5 # X
         xx = np.array([1,1])*slon
         yy = np.array([1,1])*slat
         ngl.overlay(plot[ip], ngl.xy(wks, xx, yy, mres) )
      #-------------------------------------------------------------------------
      if tlat is not None and tlon is not None:
         mres.xyMarker = 16 # round circle
         xx = np.array([1,1])*tlon
         yy = np.array([1,1])*tlat
         ngl.overlay(plot[ip], ngl.xy(wks, xx, yy, mres) )
      #-------------------------------------------------------------------------
      # set plot subtitles

      ctr_str = ''
      # if glb_avg_list != []: ctr_str = f'glb mean: {glb_avg_list[c]:6.4}'

      hs.set_subtitles( wks, plot[ip], \
                        left_string=case_name, \
                        # center_string=ctr_str, \ 
                        right_string=var_str[v], \
                        # center_sub_string=ctr_str, \
                        font_height=subtitle_font_height)
#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

# title_str = None
# title_str = 'ANN'
# if months==[1,2,12]: title_str = 'DJF'
# if months==[6,7,8]: title_str = 'JJA'
# if title_str is not None:
#    textres =  ngl.Resources()
#    textres.txFontHeightF =  0.025
#    ngl.text_ndc(wks,title_str,0.5,.7,textres)

layout = [num_var,num_case] if var_x_case else [num_case,num_var]



if num_case==1 or num_var==1:
   layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
   
pnl_res = hs.setres_panel()

### use common label bar
if use_common_label_bar:
   pnl_res.nglPanelLabelBar = True
   pnl_res.lbTitleString      = f'{var_str[v]}'
   pnl_res.lbTitlePosition    = 'bottom'
   pnl_res.lbTitleFontHeightF = 0.02
   pnl_res.nglPanelLabelBarLabelFontHeightF = 0.01

### add panel labels
# pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
# pnl_res.nglPanelFigureStringsJust        = "TopLeft"
# pnl_res.nglPanelFigureStringsFontHeightF = 0.01
# if layout==[3,2] : pnl_res.nglPanelFigureStringsFontHeightF = 0.015

pnl_res.nglPanelYWhiteSpacePercent = 10

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.printline()
if use_snapshot: print(); print(f'{hc.tcolor.RED}WARNING - snapshot mode enabled! (ss_t={ss_t}){hc.tcolor.ENDC}'); print()

hc.trim_png(fig_file,root=os.getenv('HOME')+'/Research/E3SM/pub_figs/sohip')
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
