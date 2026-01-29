import os, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
import cmocean
#-------------------------------------------------------------------------------
case_dir,case_sub = [],[]
case,name = [],[]
clr,dsh,mrk = [],[],[]
def add_case(case_in,n=None,p=None,s=None,d=0,c='black',m=0):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); name.append(n)
   case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
var,lev_list = [],[]
def add_var(var_name,lev=-1): var.append(var_name); lev_list.append(lev)
#-------------------------------------------------------------------------------
cscratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl'
gscratch = '/global/cfs/cdirs/m3312/whannah/2022-QBO-TEST'
pscratch = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu'

# add_case('MAC-PG',    n='MAC')
# add_case('GPM-PG',    n='IMERG')

add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01',       n='control',     c='red'  ,p=gscratch,s='run')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01', n='L72 smoothed',c='green',p=gscratch,s='run')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01', n='L80 refined', c='blue' ,p=gscratch,s='run')

add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72.GW-MOD-0', n='E3SM L72',      d=1,c='black', p=pscratch,s='run')
add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-nsu40',    n='E3SM L72-nsu40',d=1,c='red',   p=pscratch,s='run')
# add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-rscl',     n='E3SM L72-rscl', d=1,c='purple',p=pscratch,s='run')
# add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-rlim',     n='E3SM L72-rlim', d=1,c='pink',  p=pscratch,s='run')
#-------------------------------------------------------------------------------

# add_var('PRECT')
# add_var('TGCLDLWP')
# add_var('TGCLDIWP')
# add_var('FLNT')
add_var('U',lev=200)

fig_type = 'png'
fig_file = 'figs/FXX-clim-map'

# lat1,lat2 = -40,40
# lon1,lon2 = 90,260

htype,years,months,first_file,num_files = 'h0',[],[],0,10*12
# htype,years,months,first_file,num_files = 'h0',[],[],0,20*12 ; num_files_obs = 365*20
# htype,years,months,first_file,num_files = 'h0',[],[],0,1 ; num_files_obs = 30

use_remap,remap_grid = False,'180x360' # 90x180 / 180x360

plot_diff,add_diff,diff_base = True,False,0

chk_significance  = False
print_stats = True
var_x_case = True
num_plot_col = 2

common_colorbar = False

recalculate = True

tmp_file_head = 'data/climatology'

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)

subtitle_font_height = 0.02 - 0.0012*num_var - 0.0014*(num_case+int(add_diff))
subtitle_font_height = max(subtitle_font_height,0.005)

subtitle_font_height = 0.008

if 'diff_case' not in vars(): diff_case = [(i+1) for i in range(num_case-1)]
if 'lev' not in vars(): lev = np.array([0])

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix

wks = ngl.open_wks(fig_type,fig_file,wkres)
if plot_diff and add_diff: 
   plot = [None]*(num_var*(num_case*2-1))
else:
   plot = [None]*(num_var*num_case)
   
res = hs.res_contour_fill_map()
if 'lat1' in vars() : res.mpMinLatF = lat1; res.mpMaxLatF = lat2
if 'lon1' in vars() : res.mpMinLonF = lon1; res.mpMaxLonF = lon2

res.tmYLLabelFontHeightF         = 0.008
res.tmXBLabelFontHeightF         = 0.008
res.lbLabelFontHeightF           = 0.018
res.tmXBOn                       = False
res.tmYLOn                       = False
# res.mpGeophysicalLineColor       = 'white'

if common_colorbar:
   res.lbLabelBarOn = False
else:
   res.lbLabelBarOn = True

def get_comp(case):
   comp = 'eam'
   if 'OBS'  in case: comp = None
   if 'MAC'  in case: comp = None
   if 'GPM'  in case: comp = None
   if 'CESM' in case: comp = 'cam'
   return comp

def get_name(case):
   tmp_name = case
   if 'MAC'  in case: tmp_name = 'MAC'
   if 'GPM'  in case: tmp_name = 'GPM'
   return tmp_name

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   data_list,area_list,lat_list,lon_list = [],[],[],[]
   std_list,cnt_list = [],[]
   if 'lev_list' in locals(): lev = lev_list[v]
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)

      data_dir_tmp,data_sub_tmp = None, None
      # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
      if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      if case_sub[c] is not None: data_sub_tmp = case_sub[c]

      case_obj = he.Case( name=case[c], atm_comp=get_comp(case[c]), data_dir=data_dir_tmp, data_sub=data_sub_tmp  )
      # case_obj = he.Case( name=get_name(case[c]), time_freq='daily')

      case_obj.set_coord_names(var[v])
      
      remap_str = 'remap_ne30pg2'
      use_remap = False
      if case[c]=='MAC-PG' : use_remap = True; remap_str=f'remap_ne30pg2'
      if case[c]=='MAC-FV' : use_remap = False
      if case[c]=='GPM-PG' : use_remap = True; remap_str=f'remap_ne30pg2'
      if case[c]=='GPM-FV' : use_remap = True; remap_str=f'remap_180x360'

      if case[c] in ['MAC-FV']: case_obj.file_path = f'{case_obj.data_dir}/{case_obj.name}/daily_cwp/*'
      if case[c] in ['MAC-PG']: case_obj.file_path = f'{case_obj.data_dir}/{case_obj.name}/daily_cwp_ne30pg2/*'

      if case[c] in ['GPM-PG']: case_obj.file_path = f'{case_obj.data_dir}/{case_obj.name}/daily_ne30pg2/*'
      if case[c] in ['GPM-FV']: case_obj.file_path = f'{case_obj.data_dir}/{case_obj.name}/daily_180x360/*'

      #-------------------------------------------------------------------------
      # read the data
      #-------------------------------------------------------------------------   
      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}.nc'
      print(f'    tmp_file: {tmp_file}')
      if recalculate :

         num_files_tmp = num_files

         if 'MAC' in case[c]: num_files_tmp = num_files_obs
         if 'GPM' in case[c]: num_files_tmp = num_files_obs

         data = case_obj.load_data(var[v], component=get_comp(case[c]),htype=htype,ps_htype=htype,
                                         years=years,months=months,lev=lev,
                                         first_file=first_file,num_files=num_files_tmp,
                                         use_remap=use_remap,remap_str=remap_str)

         area = case_obj.load_data('area',htype=htype,num_files=1,
                                    use_remap=use_remap,remap_str=remap_str).astype(np.double)
         
         # Get rid of lev dimension
         if 'lev' in data.dims : data = data.isel(lev=0)

         # print stats before time averaging
         if print_stats: hc.print_stat(data,name=var[v],stat='naxsh',indent='    ')

         #----------------------------------------------------------------------
         #----------------------------------------------------------------------

         # average over time dimension
         if 'time' in data.dims :
            time = data.time
            hc.print_time_length(data.time,indent=' '*6)
            data = data.mean(dim='time')
      
         # combine lat/lon dimensions for plotting
         if '-FV' in case[c]: data = data.stack(ncol=("lat", "lon"))

         if 'GPM' in case[c]: data = data*24. # convert from mm/hr => mm/day
         
         #----------------------------------------------------------------------
         # Calculate area weighted global mean
         #----------------------------------------------------------------------
         if 'area' in locals() :
            gbl_mean = ( (data*area).sum() / area.sum() ).values 
            print(f'      Area Weighted Global Mean : {gbl_mean:6.4}')
         #----------------------------------------------------------------------
         # write to temporary file
         #----------------------------------------------------------------------
         print(f'    writing to file: {tmp_file}')
         data.name = var[v]
         tmp_ds = xr.Dataset()
         tmp_ds['data'] = data
         if 'gbl_mean' in locals(): tmp_ds['gbl_mean'] = gbl_mean
         if 'time'     in locals(): tmp_ds['time'] = time
         print(f'      writing to file: {tmp_file}')
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         tmp_ds = xr.open_dataset( tmp_file )
         data = tmp_ds['data']
         gbl_mean = tmp_ds['gbl_mean']

      if 'time'in locals(): del time
      #-------------------------------------------------------------------------
      # append to data lists
      #-------------------------------------------------------------------------
      if case[c]=='TRMM' and 'lon1' not in locals(): 
         data_list.append( ngl.add_cyclic(data.values) )
      else:
         data_list.append( data.values )

      if 'area' in locals() : area_list.append( area.values )
      #-------------------------------------------------------------------------
      # save baseline for diff map
      #-------------------------------------------------------------------------
      if plot_diff :
         if c==diff_base:
            data_baseline = data.copy()

   #------------------------------------------------------------------------------------------------
   # Plot averaged data
   #------------------------------------------------------------------------------------------------
   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])

   if plot_diff:
      tmp_data = data_list - data_list[diff_base]
      for c in range(num_case): tmp_data[c] = data_list[c] - data_list[diff_base]
      diff_data_min = np.min([np.nanmin(d) for d in tmp_data])
      diff_data_max = np.max([np.nanmax(d) for d in tmp_data])

   for c in range(num_case):
      # case_obj = he.Case( name=case[c], atm_comp=get_comp(case[c]), data_dir=data_dir, data_sub=data_sub )

      # case_obj = he.Case( name=get_name(case[c]) )
      # case_obj.set_coord_names(var[v])

      num_case_alt = num_case*2-1 if (plot_diff and add_diff) else num_case
      ip = v*num_case_alt+c if var_x_case else c*num_var+v
      #-------------------------------------------------------------------------
      # Set colors
      #-------------------------------------------------------------------------
      tres = copy.deepcopy(res)
      # tres.cnFillPalette = "MPL_viridis"
      tres.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )

      #-------------------------------------------------------------------------
      # Set explicit contour levels
      #-------------------------------------------------------------------------
      if var[v] in ['PRECT','PRECC']   : tres.cnLevels = np.arange(5,150+5,5)/1e1
      # if var[v] in ['PRECT','PRECC']   : tres.cnLevels = np.logspace( -2, 1.31, num=60).round(decimals=2)
      # if var[v]=="TGCLDIWP"            : tres.cnLevels = np.logspace( -2, 0.25, num=20).round(decimals=2)
      # if var[v]=="TGCLDLWP"            : tres.cnLevels = np.logspace( -2, 0.25, num=20).round(decimals=2)
      if var[v]=="TGCLDIWP"            : tres.cnLevels = np.arange(0.005,0.155,0.01)
      if var[v]=="TGCLDLWP"            : tres.cnLevels = np.arange(0.01,0.25,0.015)

      #-------------------------------------------------------------------------
      # set non-explicit contour levels
      #-------------------------------------------------------------------------
      if hasattr(tres,'cnLevels') : 
         tres.cnLevelSelectionMode = 'ExplicitLevels'
      else:
         nlev = 21
         aboutZero = False
         if var[v] in ['SPTLS','SPQTLS','U','V','VOR','DIV',
                       'U850','V850','U200','V200',
                       'MMF_CVT_TEND_T','MMF_CVT_TEND_Q',] : 
            aboutZero = True
         clev_tup = ngl.nice_cntr_levels(data_min, data_max, cint=None, max_steps=nlev, \
                                         returnLevels=False, aboutZero=aboutZero )
         if clev_tup==None: 
            tres.cnLevelSelectionMode = 'AutomaticLevels'   
         else:
            cmin,cmax,cint = clev_tup
            tres.cnLevels = np.linspace(cmin,cmax,num=nlev)
            tres.cnLevelSelectionMode = 'ExplicitLevels'

      #-------------------------------------------------------------------------
      # set alternate variable names
      #-------------------------------------------------------------------------
      var_str = var[v]
      if var[v]=='PRECT':     var_str = 'Precipitation'
      if var[v]=='TMQ':       var_str = 'CWV'
      if var[v]=='TGCLDLWP':  var_str = 'Liquid Water Path'
      if var[v]=='TGCLDIWP':  var_str = 'Ice Water Path'

      lev_str = None
      if lev>0: lev_str = f'{lev}mb'
      if lev<0: lev_str = f'k={(lev*-1)}'
      if lev_str is not None and var[v] in ['U','V','OMEGA','T','Q','Z3']:
         var_str = f'{lev_str} {var[v]}'
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      scrip_file_path = 'scrip_files/ne30pg2_scrip.nc'
      if case[c]=='MAC-FV': scrip_file_path = 'scrip_files/cmip6_180x360_scrip.20210722.nc'
      if case[c]=='GPM-FV': scrip_file_path = 'scrip_files/cmip6_180x360_scrip.20210722.nc'
      
      scrip_ds = xr.open_dataset(scrip_file_path)

      tres.cnFillMode    = "CellFill"
      tres.sfXArray      = scrip_ds['grid_center_lon'].values
      tres.sfYArray      = scrip_ds['grid_center_lat'].values
      tres.sfXCellBounds = scrip_ds['grid_corner_lon'].values
      tres.sfYCellBounds = scrip_ds['grid_corner_lat'].values
      #-------------------------------------------------------------------------
      # Create plot
      #-------------------------------------------------------------------------
      
      if plot_diff and c==diff_base : base_name = name[c]

      if not plot_diff  or (plot_diff and add_diff) or (plot_diff and c==diff_base) : 

         plot[ip] = ngl.contour_map(wks,data_list[c],tres) 
         #----------------------------------------------------------------------
         # set plot subtitles
         #----------------------------------------------------------------------
         ctr_str = ''

         hs.set_subtitles(wks, plot[ip], name[c], ctr_str, var_str, font_height=subtitle_font_height)

      #-------------------------------------------------------------------------
      # create difference plot
      #-------------------------------------------------------------------------
      if plot_diff and c in diff_case :
         
         data_list[c] = data_list[c] - data_baseline.values

         tres.cnFillPalette = 'BlueWhiteOrangeRed'
         tres.cnLevelSelectionMode = "ExplicitLevels"
         
         if hasattr(tres,'cnLevels') : del tres.cnLevels
         # if var[v] in ['PRECT','PRECC','PRECL'] : tres.cnLevels = np.arange(-5,5+1,1)
         if not hasattr(tres,'cnLevels') : 
            if np.min(data_list[c])==np.max(data_list[c]) : 
               print(hc.tcolor.RED+'WARNING: Difference is zero!'+hc.tcolor.ENDC)
            else:
               cmin,cmax,cint,clev = ngl.nice_cntr_levels(diff_data_min, diff_data_max,    \
                                                          cint=None, max_steps=21,      \
                                                          returnLevels=True, aboutZero=True )
               tres.cnLevels = np.linspace(cmin,cmax,num=21)
         
         tres.lbLabelBarOn = True

         ipd = ip
         if add_diff and     var_x_case: ipd = ip+1
         if add_diff and not var_x_case: ipd = ip+num_var*(num_case-1)

         plot[ipd] = ngl.contour_map(wks,data_list[c],tres)

         
         ctr_str = ''
         # case_name = name[c]+' - '+base_name
         if 'name' in vars():
            case_name = name[c]
         else:
            case_name = case_obj.short_name
         
         # hs.set_subtitles(wks, plot[ipd], case_name, '', var_str+' (Diff)', font_height=subtitle_font_height)
         hs.set_subtitles(wks, plot[ipd], case_name, 'Difference', var_str, font_height=subtitle_font_height)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

num_case_alt = num_case*2-1 if (plot_diff and add_diff) else num_case
layout = [num_var,num_case_alt] if var_x_case else [num_case_alt,num_var]

if not (add_diff):
   if num_case==1 or num_var==1:
      layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
   
pnl_res = hs.setres_panel()

### add panel labels
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.01
if layout==[3,2] : pnl_res.nglPanelFigureStringsFontHeightF = 0.015

if common_colorbar: 
   pnl_res.nglPanelLabelBar = True
   pnl_res.lbTitleFontHeightF = 0.01
   pnl_res.nglPanelLabelBarLabelFontHeightF = 0.01
   pnl_res.lbTitlePosition = 'Bottom'
   if var[v]=='PRECT':    pnl_res.lbTitleString = 'mm/day'
   if var[v]=='TGCLDLWP': pnl_res.lbTitleString = 'kg/m2'


pnl_res.nglPanelYWhiteSpacePercent = 5

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
