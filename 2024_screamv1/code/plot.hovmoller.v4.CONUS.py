# v4 is similar to v3 except panels are placed 'manually' using ImageMagick
#-------------------------------------------------------------------------------
import os, ngl, copy, string, xarray as xr, numpy as np, glob, dask, numba, cmocean, subprocess as sp
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
data_dir,data_sub = None,None
#-------------------------------------------------------------------------------
case_name,case,case_dir,case_sub = [],[],[],[]
clr,dsh,mrk = [],[],[]
obs_flag = []
def add_case(case_in,n='',p=None,s='',g=None,d=0,c='black',m=0,r=False,obs=False):
   global case_name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); case_name.append(n)
   case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
   obs_flag.append(obs)
#-------------------------------------------------------------------------------
var = []
file_type_list = []
obs_var_list = []
obs_file_type_list = []
def add_var(var_name,file_type,obs_var=None,obs_file_type=None): 
   var.append(var_name)
   file_type_list.append(file_type)
   obs_var_list.append(obs_var)
   obs_file_type_list.append(obs_file_type)
#-------------------------------------------------------------------------------
scrip_file_sim = '/pscratch/sd/w/whannah/e3sm_scratch/files_grid/ne1024pg2_scrip.nc'
scrip_file_obs = '/pscratch/sd/w/whannah/e3sm_scratch/files_grid/IMERG_1800x3600_scrip.nc'
scrip_file_era = '~/HICCUP/files_grid/scrip_ERA5_721x1440.nc'
sim_data_root  = '/global/cfs/cdirs/e3sm/gsharing/EAMxx'
obs_data_root  = '/pscratch/sd/w/whannah/Obs'

sim_topo_file = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_ne1024np4pg2_16xconsistentSGH_20190528_converted.nc'

# add_case('DYAMOND2_SCREAMv1' ,n='SCREAMv1 Jan 2020',p=data_root)
# add_case('Apr1_2013_SCREAMv1',n='SCREAMv1 Apr 2013',p=data_root)
# add_case('DYAMOND1_SCREAMv1' ,n='SCREAMv1 Aug 2016',p=data_root)
# add_case('Oct1_2013_SCREAMv1',n='SCREAMv1 Oct 2013',p=data_root)

add_case('IMERG',             n='IMERG', p=obs_data_root,s='hourly_nc',obs=True)
add_case('Apr1_2013_SCREAMv1',n='SCREAM',p=sim_data_root)


#-------------------------------------------------------------------------------

# add_var('LW_flux_up@tom','output.scream.TOMVars.INSTANT')

# add_var('VapWaterPath','output.scream.VertIntegrals.INSTANT')
# add_var('IceWaterPath','output.scream.VertIntegrals.INSTANT')

add_var('precip','output.scream.SurfVars.INSTANT','precipitation','3B-HHR.MS.MRG.3IMERG')

# add_var('precip_liq_surf_mass','output.scream.SurfVars.INSTANT','precipAvg','3B-DAY.MS.MRG.3IMERG.2016')
# add_var('precip_ice_surf_mass','output.scream.SurfVars.INSTANT')
# add_var('surf_evap','output.scream.SurfVars.INSTANT')
# add_var('horiz_winds@bot','output.scream.SurfVars.INSTANT')

#-------------------------------------------------------------------------------
# tmp_data_path = os.getenv('HOME')+'/Research/E3SM/pub_figs/2023_screamv1_4season/data'
tmp_data_path = 'data'

# day_center,day_span = 8,1
day_center,day_span = 5,2; date_str_hov,date_str_map = 'April 3-7, 2013','April 5, 2013'

# fig_file,fig_type = f'figs/hovmoller.v4.CONUS.day_{day_center}.span_{day_span}','png'
fig_file,fig_type = f'figs/hovmoller-v4-CONUS-day-{day_center}-span-{day_span}','png'


# lat1,lat2 = 30,48 # Carbone & Tuttle 2008
# lat1,lat2 = 32,45 # Surcel et al. 2010
# lat1,lat2 = 38,49 # Clark et al. 2007 (Fig 7)
lat1,lat2 = 35,45
lon1,lon2 = 360-120,360-78
# lon1,lon2 = 360-114,360-78

# var_x_case = False

subtitle_font_height = 0.025

print_stats = False

use_common_labelbar = True

skip_plot_generation = False  # use this for testing imageMagick commands

# label_list = list(string.ascii_lowercase)
label_list = ['c','d','a','b']
label_cnt = 0
#---------------------------------------------------------------------------------------------------
def run_cmd(cmd,verbose=True,indent='  '):
   cmd_str = '\n' + indent + hc.tcolor.GREEN + cmd + hc.tcolor.ENDC + '\n'
   if verbose: print(cmd_str)
   proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, universal_newlines=True)
   (msg, err) = proc.communicate()
   if verbose and msg!='': print(f'  msg: {msg}')
   if err!='' and not verbose: print(cmd_str)
   if err!='': print(f'err: {err}'); exit()
   return msg

#---------------------------------------------------------------------------------------------------
# def trim_png(fig_file_in,verbose=False):
#    run_cmd(f'convert -trim {fig_file_in}.{fig_type} {fig_file_in}.trim.{fig_type}')
def trim_png_sq(fig_file_in,verbose=False):
   fig_type = 'png'
   cmd = f"convert -trim  {fig_file_in}.{fig_type} -format '%wx%h\n' info:"
   (msg, err) = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, universal_newlines=True).communicate()
   dim_list = np.array( [ int(d) for d in msg.replace('\n','').split('x') ] )
   max_dim = np.max(dim_list)
   run_cmd(f'convert -trim -gravity center -extent {max_dim}x{max_dim} {fig_file_in}.{fig_type} {fig_file_in}.{fig_type}')
#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_var = len(case),len(var)

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
# wks = ngl.open_wks(fig_type,fig_file,wkres)

# plot = [None]*num_case*num_var*2
res = hs.res_contour_fill()
res.vpHeightF = 0.38
# res.vpWidthF  = 0.4
res.tmYLLabelFontHeightF         = 0.01
res.tmXBLabelFontHeightF         = 0.01
res.tiYAxisString = 'Time [days]'
res.tiXAxisString = 'Longitude'
res.lbLabelFontHeightF    = 0.015
res.lbTitlePosition       = 'Bottom'
res.lbTopMarginF          = -0.2
res.lbBottomMarginF       =  0.2+0.01
# map_res.lbLeftMarginF     = 0.5
# map_res.lbRightMarginF    = 0.5
# res.cnFillPalette = 'MPL_viridis'
res.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
# res.cnFillPalette = np.array( cmocean.cm.amp(np.linspace(0,1,256)) )

map_res = hs.res_contour_fill()
# map_res.vpWidthF   = 0.2
# map_res.vpHeightF  = 0.2
map_res.mpLimitMode           = 'LatLon'
map_res.mpMinLatF             = lat1-8-5
map_res.mpMaxLatF             = lat2+8+5
map_res.mpMinLonF             = lon1-6
map_res.mpMaxLonF             = lon2+6
map_res.tiYAxisString = 'Latitude'
map_res.tiXAxisString = 'Longitude'
map_res.tmYLLabelFontHeightF  = 0.01
map_res.tmXBLabelFontHeightF  = 0.01
map_res.lbLabelFontHeightF    = 0.01
map_res.lbLeftMarginF         = 0.5
map_res.lbRightMarginF        = 0.5
map_res.lbTitlePosition       = 'Bottom'
# map_res.tmXBOn                = False
# map_res.tmYLOn                = False
map_res.mpGridAndLimbOn       = False
map_res.mpGeophysicalLineThicknessF = 2
map_res.cnFillPalette         = res.cnFillPalette

cnt_res = hs.res_contour()
cnt_res.tmXBOn                = False
cnt_res.tmYLOn                = False
cnt_res.cnLineThicknessF      = 2.
# cnt_res.mpLimitMode           = map_res.mpLimitMode
# cnt_res.mpMinLatF             = map_res.mpMinLatF
# cnt_res.mpMaxLatF             = map_res.mpMaxLatF
# cnt_res.mpMinLonF             = map_res.mpMinLonF
# cnt_res.mpMaxLonF             = map_res.mpMaxLonF

if use_common_labelbar:
   res.lbLabelBarOn = False
   map_res.lbLabelBarOn = False

lres = hs.res_xy()
lres.xyDashPattern    = 0
lres.xyLineThicknessF = 4
lres.xyLineColor      = 'red'

pgres = ngl.Resources()
pgres.nglDraw,pgres.nglFrame = False,False
pgres.gsLineColor = 'red'
pgres.gsLineThicknessF = 6

if 'lev' not in locals(): lev = np.array([0])

# pdum = [None]*len(plot) # for map box

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

sim_grid_ds = xr.open_dataset( scrip_file_sim ).rename({'grid_size':'ncol'})
# sim_grid_ds['grid_center_lon'] = xr.where(sim_grid_ds['grid_center_lon']<0,sim_grid_ds['grid_center_lon']+360,sim_grid_ds['grid_center_lon'])
# sim_grid_ds['grid_corner_lon'] = xr.where(sim_grid_ds['grid_corner_lon']<0,sim_grid_ds['grid_corner_lon']+360,sim_grid_ds['grid_corner_lon'])

obs_grid_ds = xr.open_dataset( scrip_file_obs ).rename({'grid_size':'ncol'})
obs_grid_ds['grid_center_lon'] = xr.where(obs_grid_ds['grid_center_lon']<0,obs_grid_ds['grid_center_lon']+360,obs_grid_ds['grid_center_lon'])
obs_grid_ds['grid_corner_lon'] = xr.where(obs_grid_ds['grid_corner_lon']<0,obs_grid_ds['grid_corner_lon']+360,obs_grid_ds['grid_corner_lon'])

era_grid_ds = xr.open_dataset( scrip_file_era ).rename({'grid_size':'ncol'})
era_grid_ds['grid_center_lon'] = xr.where(era_grid_ds['grid_center_lon']<0,era_grid_ds['grid_center_lon']+360,era_grid_ds['grid_center_lon'])
era_grid_ds['grid_corner_lon'] = xr.where(era_grid_ds['grid_corner_lon']<0,era_grid_ds['grid_corner_lon']+360,era_grid_ds['grid_corner_lon'])

#---------------------------------------------------------------------------------------------------
# plot separate colorbar
#---------------------------------------------------------------------------------------------------

precip_clev = np.logspace( np.log10(0.5), np.log10(150), num=10).round(decimals=4)
lb_clr = np.array( cmocean.cm.rain(np.linspace(0,1,len(precip_clev))) )

if use_common_labelbar:
   lb_file = 'figs/hovmoller.v4.CONUS.colorbar'
   lbwkres = ngl.Resources()
   wks = ngl.open_wks(fig_type,lb_file,lbwkres)
   lb_res = ngl.Resources()
   # lb_res.vpWidthF             = 0.3
   lb_res.vpHeightF            = 0.1
   lb_res.lbFillColors         = lb_clr
   lb_res.lbMonoFillPattern    = True
   lb_res.lbFillPattern        = 'SolidFill'
   lb_res.lbOrientation        = 'Horizontal'
   lb_res.lbLabelFontHeightF   = 0.012
   lb_res.lbTitleFontHeightF   = 0.015
   lb_res.lbTitlePosition      = 'Bottom'
   lb_res.lbTitleString        = '[mm/day]'
   # lb_res.lbTopMarginF         =  0.25
   # lb_res.lbBottomMarginF      = -0.1
   # lb_res.lbLeftMarginF        =  0.2
   # lb_res.lbRightMarginF       =  0.2
   clev_str = []
   for c in precip_clev: clev_str.append(f'{c:.1f}')
   plot = ngl.labelbar_ndc(wks, len(clev_str), clev_str, 0.3,.3,lb_res)
   ngl.frame(wks)
   run_cmd(f'convert -trim -gravity center -extent 2048x180 {lb_file}.{fig_type} {lb_file}.{fig_type}')
   ngl.destroy(wks)

#---------------------------------------------------------------------------------------------------
# combine images
def combine_images():
   fig_list_str = ''
   print('\n  Combining figures...')
   for v in range(num_var):
      for n in range (2):
         for c in range(num_case):
            fig_file_hov = f'{fig_file}.hov.v_{v}.c_{c}.{fig_type}'
            fig_file_map = f'{fig_file}.map.v_{v}.c_{c}.{fig_type}'
            # if n==0: 
            #    print(' '*2+f'fig_file_hov: {fig_file_hov}')
            #    fig_list_str += ' '+fig_file_hov
            #    print(' '*2+f'fig_file_map: {fig_file_map}')
            #    fig_list_str += ' '+fig_file_map
            # if n==0: 
            #    print(' '*2+f'fig_file_hov: {fig_file_hov}')
            #    fig_list_str += ' '+fig_file_hov
            # if n==1: 
            #    print(' '*2+f'fig_file_map: {fig_file_map}')
            #    fig_list_str += ' '+fig_file_map

            fig_file_hov_trim = fig_file_hov.replace(f'.{fig_type}',f'.trim.{fig_type}')
            fig_file_map_trim = fig_file_map.replace(f'.{fig_type}',f'.trim.{fig_type}')
            run_cmd(f'convert -trim {fig_file_hov} {fig_file_hov_trim}')
            run_cmd(f'convert -trim {fig_file_map} {fig_file_map_trim}')
            if n==0: 
               print(' '*2+f'fig_file_map: {fig_file_map_trim}')
               fig_list_str += ' '+fig_file_map_trim
            if n==1: 
               print(' '*2+f'fig_file_hov: {fig_file_hov_trim}')
               fig_list_str += ' '+fig_file_hov_trim

   # geo = '+0+0'
   # geo,layout = '1024x1024>+10+10',f'{num_var*2}x{num_case}'
   # run_cmd(f'montage {fig_list_str} -tile {layout} -geometry {geo} {fig_file}.{fig_type}')
   # run_cmd(f'montage {fig_list_str} -geometry 1024x1024+5 -tile {num_var*2}x{num_case} {fig_file}.{fig_type}')
   run_cmd(f'montage {fig_list_str} -geometry 1024x800+5 -tile {num_var*2}x{num_case} {fig_file}.{fig_type}')

   if use_common_labelbar:
      run_cmd(f'montage {fig_file}.{fig_type} {lb_file}.{fig_type} -geometry +10 -tile 1x2  -gravity center {fig_file}.{fig_type}')
      # run_cmd(f'montage {fig_file}.{fig_type} {lb_file}.{fig_type} -tile 1x2  -gravity center {fig_file}.{fig_type}')

   run_cmd(f'convert -trim {fig_file}.{fig_type} {fig_file}.{fig_type}')

   print(f'\n{fig_file}.{fig_type}\n')
#---------------------------------------------------------------------------------------------------

if skip_plot_generation: combine_images(); exit()
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   time_list = []
   lon_list = []
   hov_list = []
   map_list = []
   hgt_list = []
   for c in range(num_case):
      if obs_flag[c]:
         tfile_type = obs_file_type_list[v]
         tvar = obs_var_list[v]
      else:
         tfile_type = file_type_list[v]
         tvar = var[v]
      if c==0: print(' '*2+'var: '+hc.tcolor.GREEN+tvar+hc.tcolor.ENDC)
      print('\n'+' '*4+'case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)
      #-------------------------------------------------------------------------
      # Load precalculated hovmoller data from v2
      tmp_file = f'{tmp_data_path}/CONUS.hovmoller.horz.v1.tmp.{tvar}.{case[c]}'
      if 'lat1' in globals(): tmp_file += f'.lat1_{lat1}.lat2_{lat2}'
      if 'lon1' in globals(): tmp_file += f'.lon1_{lon1}.lon2_{lon2}'
      tmp_file += '.nc'
      print(' '*6+f'Reading hov data from file: {hc.tcolor.MAGENTA}{tmp_file}{hc.tcolor.ENDC}')
      tmp_ds = xr.open_dataset( tmp_file )
      hov    = tmp_ds[tvar]
      #-------------------------------------------------------------------------
      hov['time'] = hov['time'] - hov['time'][0]
      hov = hov.sel(time=slice(day_center-day_span, day_center+day_span),drop=True)

      lon_bins = hov['lon'].values
      time     = hov['time'].values

      time_list.append(time)
      lon_list.append(lon_bins)
      hov_list.append(hov.values)

      #-------------------------------------------------------------------------
      # idenfity the file to load for map plot
      file_list = sorted(glob.glob(f'{case_dir[c]}/{case[c]}/{case_sub[c]}/{tfile_type}*'))
      f_ind = int(day_center) if not obs_flag[c] else int(day_center*48)
      #-------------------------------------------------------------------------
      # load data for map
      print(' '*6+f'Reading map data from file: {hc.tcolor.LIGHTRED}{file_list[f_ind]}{hc.tcolor.ENDC}')
      group_name = 'Grid' if case[c]=='IMERG' else None
      data_ds = xr.open_dataset( file_list[f_ind], group=group_name)
      if tvar=='precip':
         data = data_ds['precip_liq_surf_mass'] + data_ds['precip_ice_surf_mass']
      else:
         data = data_ds[tvar]
      if obs_flag[c]:
         data = data.isel(time=0)
         data['lon'] = xr.where(data['lon']<0,data['lon']+360,data['lon'])
         data = data.stack(ncol=('lat','lon'))
      else:
         data = data.isel(time=slice(0,2)).mean(dim='time') # average 2x 15-min to compare with IMERG
      #-------------------------------------------------------------------
      # convert units
      if var[v]=='precip':
         if obs_flag[c]:
            data = data*24. # mm/hr to mm/day
         else:
            data = (data/100)*86400 # kg/m2 to mm/day using dtime=100 (ne1024)
      #-------------------------------------------------------------------
      map_list.append(data.values)
      #-------------------------------------------------------------------
      if print_stats: hc.print_stat(hov,name=tvar,compact=True,indent=' '*6)
      if print_stats: hc.print_stat(data,name=tvar,compact=True,indent=' '*6)
   
      #-------------------------------------------------------------------------
      # Load geopotential height for contours
      #-------------------------------------------------------------------------
      hgt_lev = 500
      if obs_flag[c]:
         hgt_case = 'ERA5'
         hgt_sub = 'daily'
         tfile_type = 'ERA5.daily.Z.atm.'
         tvar_hgt = 'z'
      else:
         hgt_case = case[c]
         hgt_sub = case_sub[c]
         tfile_type = 'output.scream.PresLevs.INSTANT.nmins_x15'
         tvar_hgt = f'VerticalLayerInterface@{hgt_lev}mb'
         # tvar_hgt = 'VerticalLayerInterface@200mb'
      # if tfile_type is None:
      #    hgt_list.append(None)
      # else:
      if True:
         #-------------------------------------------------------------------------
         # idenfity the file to load for map plot of geopotential heights
         file_list = sorted(glob.glob(f'{case_dir[c]}/{hgt_case}/{hgt_sub}/{tfile_type}*'))
         f_ind = int(day_center) #if not obs_flag[c] else int(day_center*48)
         #-------------------------------------------------------------------------
         # load data for map
         print(' '*6+f'Reading map data from file: {hc.tcolor.LIGHTRED}{file_list[f_ind]}{hc.tcolor.ENDC}')
         data_ds = xr.open_dataset( file_list[f_ind])
         data = data_ds[tvar_hgt]
         if obs_flag[c]:
            data = data.rename({'longitude':'lon','latitude':'lat'})
            data = data.isel(time=0).sel(level=hgt_lev)
            data['lon'] = xr.where(data['lon']<0,data['lon']+360,data['lon'])
            data = data.stack(ncol=('lat','lon'))
         else:
            data = data.isel(time=slice(0,2)).mean(dim='time') # average 2x 15-min to compare with IMERG
            mask = data.values<0
            topo_ds = xr.open_dataset( sim_topo_file, group=group_name)
            data = data + topo_ds['PHIS'].values / 9.81
            data = np.ma.masked_where( mask, data.values )
         #-------------------------------------------------------------------
         # # convert units
         if     obs_flag[c]: data = data/1e3/9.81 # m2/s2 to km
         if not obs_flag[c]: data = data/1e3      # m to km
         #-------------------------------------------------------------------
         if print_stats: hc.print_stat(data,name=tvar_hgt,compact=True,indent=' '*6)
         hgt_list.append(data)
      
   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   # unit_str = ''
   # if var in ['PRECT','PRECC','PRECL']   : unit_str = '[mm/day]'
   # res.tiXAxisString = unit_str
   # var_str = var
   # if var=='PRECT' : var_str = 'Precipitation'
   # res.trXMinF = bin_ds['bins'].min().values
   # res.trXMaxF = bin_ds['bins'].max().values

   #-------------------------------------------------------------------------
   # Create hovmoller plot
   #-------------------------------------------------------------------------
   if var[v]=='precip': 
      # res.cnLevels = np.logspace( -1, 1.35, num=10).round(decimals=4)
      # res.cnLevels = np.logspace( np.log10(0.5), np.log10(150), num=10).round(decimals=4)
      res.cnLevels = precip_clev
      res.cnLevelSelectionMode = 'ExplicitLevels'
   #-------------------------------------------------------------------------
   for c in range(num_case):
      fig_file_tmp = f'{fig_file}.hov.v_{v}.c_{c}'
      wks = ngl.open_wks(fig_type,fig_file_tmp,wkres)
      #-------------------------------------------------------------------------
      tres = copy.deepcopy(res)
      # tres.tmXBMode      = 'Explicit'
      # tres.tmXBValues    = np.arange( len(lon_bins) )
      # tres.tmXBLabels    = lon_bins
      # if var[v]=='precip': tres.lbTitleString = 'mm/day'
      tres.sfYArray = time_list[c]
      tres.sfXArray = lon_list[c]
      plot = ngl.contour(wks, hov_list[c] ,tres) 
      #-------------------------------------------------------------------------
      # add line to show time of map
      ngl.overlay( plot, ngl.xy(wks, np.array([-1e3,1e8]), day_center*np.array([1,1]), lres) )
      #-------------------------------------------------------------------------
      # hs.set_subtitles(wks, plot, case_name[c], '', tvar, font_height=subtitle_font_height)
      hs.set_subtitles(wks, plot, case_name[c], '', date_str_hov, font_height=subtitle_font_height)
      #-------------------------------------------------------------------------
      pnl_res = ngl.Resources()
      pnl_res.nglPanelFigureStringsJust        = 'TopLeft'
      pnl_res.nglPanelFigureStringsFontHeightF = 0.02
      pnl_res.nglPanelFigureStrings            = [label_list[label_cnt]] ; label_cnt+=1
      ngl.panel(wks,[plot],[1,1],pnl_res)
      # ngl.draw(plot); ngl.frame(wks)
      trim_png_sq(fig_file_tmp,verbose=False)
      ngl.destroy(wks)

   #----------------------------------------------------------------------------
   # Create map plot
   #----------------------------------------------------------------------------
   if var[v]=='precip': 
      # map_res.cnLevels = np.logspace( -3, 0.5, num=30).round(decimals=4)
      # map_res.cnLevels = np.logspace( -1, 2.2, num=10).round(decimals=4)
      # map_res.cnLevels = np.logspace( np.log10(0.5), np.log10(150), num=10).round(decimals=4)
      map_res.cnLevels = precip_clev
      # map_res.cnLevels = res.cnLevels
      map_res.cnLevelSelectionMode = 'ExplicitLevels'
   #----------------------------------------------------------------------------
   for c in range(num_case):
      fig_file_tmp = f'{fig_file}.map.v_{v}.c_{c}'
      wks = ngl.open_wks(fig_type,fig_file_tmp,wkres)
      #-------------------------------------------------------------------------
      tres = copy.deepcopy(map_res)
      grid_ds = sim_grid_ds if not obs_flag[c] else obs_grid_ds
      tres.cnFillMode    = 'CellFill'
      tres.sfXArray      = grid_ds['grid_center_lon'].values
      tres.sfYArray      = grid_ds['grid_center_lat'].values
      tres.sfXCellBounds = grid_ds['grid_corner_lon'].values
      tres.sfYCellBounds = grid_ds['grid_corner_lat'].values
      # if var[v]=='precip': tres.lbTitleString = 'mm/day'
      #-------------------------------------------------------------------------
      plot = ngl.contour_map(wks,map_list[c],tres)
      #-------------------------------------------------------------------------
      if hgt_list[c] is not None:
         tres = copy.deepcopy(cnt_res)
         grid_ds = sim_grid_ds if not obs_flag[c] else era_grid_ds
         tres.cnLevelSelectionMode = 'AutomaticLevels'
         tres.sfXArray      = grid_ds['grid_center_lon'].values
         tres.sfYArray      = grid_ds['grid_center_lat'].values
         tres.cnLevelSelectionMode = 'ExplicitLevels'
         tres.cnLevels             = np.arange(4.5,6.0+0.05,0.05)
         ngl.overlay( plot , ngl.contour(wks,hgt_list[c],tres) )
         # plot = ngl.contour_map(wks,hgt_list[c],tres)
      # else:
      #    plot = ngl.contour_map(wks,map_list[c],tres)
      #-------------------------------------------------------------------------
      # draw box around map region used for inset
      bx = np.array([lon1,lon2,lon2,lon1,lon1])
      by = np.array([lat1,lat1,lat2,lat2,lat1])
      pgres.gsLineColor = 'red'
      pdum = ngl.add_polyline(wks,plot,bx,by,pgres)
      #-------------------------------------------------------------------------
      # draw boxes around regions for next plot
      box_list = []
      box_list.append([37,40,360-108,360-105,'MTN'])
      box_list.append([37,40,360-104,360-101,'HP'])
      box_list.append([37,40,360-100,360- 97,'MP'])
      box_list.append([37,40,360- 96,360- 93,'LP'])
      box_list.append([31,34,360- 90,360- 82,'SOUTHEAST'])
      bdum = [None]*len(box_list)
      for b in range(len(box_list)):
         [by1,by2,bx1,bx2,bstr] = box_list[b] 
         bx = np.array([bx1,bx2,bx2,bx1,bx1])
         by = np.array([by1,by1,by2,by2,by1])
         pgres.gsLineColor = 'blue'
         pgres.gsLineThicknessF = 4
         pgres.gsLineDashPattern = 1 if bstr=='SOUTHEAST' else 0
         bdum[b] = ngl.add_polyline(wks,plot,bx,by,pgres)
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      # hs.set_subtitles(wks, plot, case_name[c], '', tvar, font_height=subtitle_font_height)
      hs.set_subtitles(wks, plot, case_name[c], '', date_str_map, font_height=subtitle_font_height)
      #-------------------------------------------------------------------------
      pnl_res = ngl.Resources()
      pnl_res.nglFrame = False
      pnl_res.nglPanelFigureStringsJust        = 'TopLeft'
      pnl_res.nglPanelFigureStringsFontHeightF = 0.02
      pnl_res.nglPanelFigureStrings            = [label_list[label_cnt]] ; label_cnt+=1
      ngl.panel(wks,[plot],[1,1],pnl_res)
      #-------------------------------------------------------------------------
      # Add text labels to regional boxes
      ttres               = ngl.Resources()
      ttres.nglDraw       = True
      ttres.txFontHeightF = 0.015
      ttres.txFontColor   = 'blue'
      tdum = [None]*len(box_list)
      for b in range(len(box_list)):
         [by1,by2,bx1,bx2,bstr] = box_list[b] 
         bx = np.array([bx1,bx2,bx2,bx1,bx1])
         by = np.array([by1,by1,by2,by2,by1])
         ndcx,ndcy = ngl.datatondc(plot, (bx1+bx2)/2, by1-(by2-by1)*0.8 )
         tdum[b] = ngl.text_ndc(wks, bstr, ndcx, ndcy, ttres)
      ngl.frame(wks)
      #-------------------------------------------------------------------------
      # ngl.draw(plot); ngl.frame(wks)
      trim_png_sq(fig_file_tmp,verbose=False)
      run_cmd(f'convert -trim -gravity center -extent 2048x2048 {fig_file_tmp}.{fig_type} {fig_file_tmp}.{fig_type}')
      ngl.destroy(wks)

      # exit('\n  Stopping to check annotations.')

#---------------------------------------------------------------------------------------------------
combine_images()
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
