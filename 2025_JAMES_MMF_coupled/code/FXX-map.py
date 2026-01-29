# v2 is a cleaned up version of v1 - functionality is mostly the same
import os, glob, ngl, subprocess as sp, numpy as np, xarray as xr, dask, copy, string, cmocean
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
var,lev_list,var_str = [],[],[]
def add_var(var_in,lev=-1,name=None): 
   var.append(var_in); lev_list.append(lev)
   var_str.append(var_in if name is None else name)
#-------------------------------------------------------------------------------
scrip_file_path  = os.getenv('HOME')+'/E3SM/data_grid/ne30pg2_scrip.nc'
tmp_path_hst_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_v2  = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
tmp_path_co2_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_sub = 'archive/atm/hist'
add_case('v2.LR.historical_0101',                                  n='E3SMv2',  c='cyan',    p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0151',                                  n='E3SMv2',  c='orange', p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0201',                                  n='E3SMv2',  c='green',  p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0251',                                  n='E3SMv2',  c='purple',   p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0301',                                  n='E3SMv2',  c='pink', p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical',                                       n='E3SMv2 ens',  c='red',  p=None, s=None)
add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',c='blue', p=tmp_path_hst_mmf,s='archive/atm/hist')
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2',n='1xCO2',c='blue' ,p=tmp_path_co2_mmf,s=tmp_sub)
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.2xCO2',n='2xCO2',c='green',p=tmp_path_co2_mmf,s=tmp_sub)
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.4xCO2',n='4xCO2',c='red'  ,p=tmp_path_co2_mmf,s=tmp_sub)
#-------------------------------------------------------------------------------

fig_file,fig_type = 'figs/FXX-map','png'

# add_var('PRECT',name='Precipitation')
# add_var('TS',name='Tsfc')
# add_var('TMQ')
# add_var('LHFLX')
# add_var('SHFLX')
# add_var('TGCLDLWP')
# add_var('TGCLDIWP')
# add_var('NET_TOA_RAD')
# add_var('FSNT'); add_var('FLNT')
# add_var('FSNS'); add_var('FLNS')
add_var('FSNS')
add_var('FSNSC')

# add_var('CLDLOW'); add_var('CLDLOW_CAL')

# add_var('CLDTOT_ISCCP')

#-------------------------------------------------------------------------------
# lat1,lat2 = -60,60

htype,yr1,yr2 = 'ha',1979,2014
# htype,first_file,num_files = 'ha',0,10
# htype,first_file,num_files = 'ha',100,20

use_remap,remap_grid = False,'90x180'

plot_diff,add_diff = True,False

print_stats          = True
var_x_case           = False
num_plot_col         = 1 # len(case)
use_common_label_bar = False

#---------------------------------------------------------------------------------------------------
# Set up plot resources
if case==[]: raise ValueError('ERROR - case list is empty!')
num_var,num_case = len(var),len(case)

subtitle_font_height = 0.015

diff_base = 0

if 'scrip_file_path' not in locals(): scrip_file_path = None

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix

wks = ngl.open_wks(fig_type,fig_file,wkres)
if plot_diff and add_diff: 
   plot = [None]*(num_var*(num_case*2-1))
else:
   plot = [None]*(num_var*num_case)
   
res = hs.res_contour_fill_map()
if 'lat1' in locals() : res.mpMinLatF = lat1; res.mpMaxLatF = lat2
if 'lon1' in locals() : res.mpMinLonF = lon1; res.mpMaxLonF = lon2

res.tmYLLabelFontHeightF         = 0.008
res.tmXBLabelFontHeightF         = 0.008
res.lbLabelFontHeightF           = 0.01
res.tmXBOn                       = False
res.tmYLOn                       = False
# res.mpGeophysicalLineColor       = 'white'

# res.mpProjection = 'Mollweide'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   data_list = []
   glb_avg_list = []
   lat_list,lon_list = [],[]
   if 'lev_list' in locals(): lev = lev_list[v]
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)
      #-------------------------------------------------------------------------
      data_dir_tmp,data_sub_tmp = None, None
      if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
      if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      if case_sub[c] is not None: data_sub_tmp = case_sub[c]
      case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp )
      case_obj.set_coord_names(var[v])
      #-------------------------------------------------------------------------
      # read the data
      # if use_remap or ('CESM' in case[c] and 'ne30' not in case[c]):
      #    lat = case_obj.load_data('lat',htype=htype,
      #                            use_remap=use_remap,remap_str=f'remap_{remap_grid}')
      #    lon = case_obj.load_data('lon',htype=htype,
      #                            use_remap=use_remap,remap_str=f'remap_{remap_grid}')
      # else:
      #    aname = case_obj.area_name
      #    area = case_obj.load_data(aname,htype=htype)
      #-------------------------------------------------------------------------
      # read the data
      # with dask.config.set(**{'array.slicing.split_large_chunks': True}):
      #    data = case_obj.load_data(var[v],htype=htype,ps_htype=htype,lev=lev,
      #                                    first_file=first_file,num_files=num_files,
      #                                    use_remap=use_remap,remap_str=f'remap_{remap_grid}')
      #-------------------------------------------------------------------------
      # read the data
      file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/*.eam.{htype}.*'
      file_list_all = sorted(glob.glob(file_path))
      file_list1 = [] # data for anomalies
      file_list2 = [] # data for baseline
      for f in range(len(file_list_all)):
         yr = int(file_list_all[f][-7:-7+4])
         if yr>= yr1 and yr<= yr2: file_list1.append(file_list_all[f])
      if file_list1==[]: exit(f'\nERROR: no files found for file_path:\n{file_path}\n')   
      ds = xr.open_mfdataset( file_list1 )
      #-------------------------------------------------------------------------
      data = ds[var[v]]
      #-------------------------------------------------------------------------
      # Special handling of various specific circumstances
      if 'lev' in data.dims : data = data.isel(lev=0)
      if var[v]=='TS': data = data-273
      #-------------------------------------------------------------------------
      # print stats before time averaging
      if print_stats: hc.print_stat(data,name=var[v],stat='naxsh',indent='    ',compact=True)
      #-------------------------------------------------------------------------
      # average over time dimension
      if 'time' in data.dims : 
         hc.print_time_length(data.time,indent=' '*6)
         data = data.mean(dim='time')
      #-------------------------------------------------------------------------
      # Calculate area weighted global mean
      if 'area' in locals() :
         gbl_mean = ( (data*area).sum() / area.sum() ).values 
         print(hc.tcolor.CYAN+f'      Area Weighted Global Mean : {gbl_mean:6.4}'+hc.tcolor.ENDC)
         # glb_avg_list.append(gbl_mean)
      #-------------------------------------------------------------------------
      # append to data lists
      if case[c]=='TRMM' and 'lon1' not in locals(): data = ngl.add_cyclic(data.values)
      
      data_list.append( data.values )
      #-------------------------------------------------------------------------
      # save baseline for diff map
      if plot_diff :
         if c==diff_base:
            data_baseline = data.copy()
   #----------------------------------------------------------------------------
   # calculate common limits for consistent contour levels
   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])
   if plot_diff:
      tmp_data = data_list - data_list[diff_base]
      for c in range(num_case): tmp_data[c] = data_list[c] - data_list[diff_base]
      diff_data_min = np.min([np.nanmin(d) for d in tmp_data])
      diff_data_max = np.max([np.nanmax(d) for d in tmp_data])
   #----------------------------------------------------------------------------
   # Plot averaged data
   for c in range(num_case):

      data_dir_tmp,data_sub_tmp = None, None
      if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
      if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      if case_sub[c] is not None: data_sub_tmp = case_sub[c]
      
      case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp )
      case_obj.set_coord_names(var[v])
      #-------------------------------------------------------------------------
      # Set color palette
      tres = copy.deepcopy(res)
      # tres.cnFillPalette = "MPL_viridis"
      # tres.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
      # tres.cnFillPalette = np.array( cmocean.cm.amp(np.linspace(0,1,256)) )
      # tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
      tres.cnFillPalette = np.array( cmocean.cm.thermal(np.linspace(0,1,256)) )
      # if var[v] in ['P-E']                      : tres.cnFillPalette = "BlueWhiteOrangeRed"
      # if var[v] in ['CLDLOW','CLDMED','CLDHGH'] : tres.cnFillPalette = "CBR_wet"
      # if var[v] in ['TGCLDLWP','TGCLDIWP']      : tres.cnFillPalette = "MPL_viridis"
      # if var[v] in ['DYN_QLIQ']                 : tres.cnFillPalette = "MPL_viridis"
      # if var[v] in ['TS','PS']                  : tres.cnFillPalette = 'BlueWhiteOrangeRed'
      # # if var[v] in ['TS','PS']                  : tres.cnFillPalette = 'WhiteBlueGreenYellowRed'
      # # if var[v] in ['TS','PS'] : tres.cnFillPalette = np.array( cmocean.cm.amp(np.linspace(0,1,256)) )
      # if var[v] in ['U','V','UBOT','VBOT','U850','V850','U200','V200']: 
      #    # tres.cnFillPalette = "BlueWhiteOrangeRed"
      #    tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )

      #-------------------------------------------------------------------------
      # Set contour levels
      if var[v] == 'PRECT'       : tres.cnLevels = np.arange(2,20+2,2)
      if var[v]=='LHFLX'         : tres.cnLevels = np.arange(5,205+5,5)
      if var[v]=='TS'            : tres.cnLevels = np.arange(-54,54+8,8)

      # if var[v]=='TGCLDIWP'      : tres.cnLevels = np.arange(1,30+1,1)*1e-2
      # if var[v]=='TGCLDLWP'      : tres.cnLevels = np.logspace( -2, 0.25, num=60).round(decimals=2)
      

      if var[v] in ['U','V']: 
         if plot_diff and not add_diff and c!=diff_base :
            if lev==850: tres.cnLevels = np.arange( -8, 8+1,1)
            if lev==200: tres.cnLevels = np.arange(-16,16+2,2)
         else:
            if lev==850: tres.cnLevels = np.arange(-20,20+2,2)
            if lev==200: tres.cnLevels = np.arange(-60,60+6,6)
      #-------------------------------------------------------------------------
      # set non-explicit contour levels
      if hasattr(tres,'cnLevels') : 
         tres.cnLevelSelectionMode = 'ExplicitLevels'
      else:
         nlev = 21
         aboutZero = False
         if var[v] in ['U','V','VOR','DIV','U850','V850','U200','V200']: 
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
      # Create plot
      if use_remap \
      or case_obj.obs \
      or ('CESM' in case[c] and 'ne30' not in case[c]) :
         if case[c]=='ERAi': lat,lon = erai_lat,erai_lon
         hs.set_cell_fill(tres,case_obj=case_obj,lat=lat,lon=lon)
      else:
         if case[c]=='v2.LR.amip_0101' : case_obj.grid = 'ne30pg2'
         hs.set_cell_fill(tres,case_obj=case_obj,htype=htype,scrip_file_path=scrip_file_path)

      tres.lbLabelBarOn = False if use_common_label_bar else True      
         
      # if plot_diff and c==diff_base : base_name = name[c]

      num_case_alt = num_case*2-1 if (plot_diff and add_diff) else num_case
      ip = v*num_case_alt+c if var_x_case else c*num_var+v

      if not plot_diff  or (plot_diff and add_diff) or (plot_diff and c==diff_base) : 

         plot[ip] = ngl.contour_map(wks,np.ma.masked_invalid(data_list[c]),tres) 
         
         # set plot subtitles
         ctr_str = f'{glb_avg_list[c]:6.4}' if glb_avg_list!=[] else ''
         hs.set_subtitles(wks, plot[ip], name[c], ctr_str, var_str[v], font_height=subtitle_font_height)

      #-------------------------------------------------------------------------
      # create difference plot
      if plot_diff and c!=diff_base :
         
         data_list[c] = data_list[c] - data_baseline.values
         
         # tres.cnFillPalette = 'BlueWhiteOrangeRed'
         # tres.cnFillPalette = np.array( cmocean.cm.diff(np.linspace(0,1,256)) )
         # tres.cnFillPalette = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )
         tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )

         if hasattr(tres,'cnLevels') : del tres.cnLevels
         if var[v] in ['PRECT','PRECC','PRECL'] : tres.cnLevels = np.arange(-5,5+1,1)
         if var[v]=='TS'                        : tres.cnLevels = np.arange(-28,28+4,4)
         if not hasattr(tres,'cnLevels') : 
            if np.min(data_list[c])==np.max(data_list[c]) : 
               print(hc.tcolor.RED+'WARNING: Difference is zero!'+hc.tcolor.ENDC)
            else:
               cmin,cmax,cint,clev = ngl.nice_cntr_levels(diff_data_min, diff_data_max,    \
                                                          cint=None, max_steps=21,      \
                                                          returnLevels=True, aboutZero=True )
               tres.cnLevels = np.linspace(cmin,cmax,num=21)
         
         tres.cnLevelSelectionMode = "ExplicitLevels"
         # tres.cnLevelSelectionMode = "AutomaticLevels" # override the level settings and just use auto

         ipd = ip
         if add_diff and     var_x_case: ipd = ip+(num_case-1)
         if add_diff and not var_x_case: ipd = ip+num_var*(num_case-1)

         tres.lbLabelBarOn = True

         plot[ipd] = ngl.contour_map(wks,np.ma.masked_invalid(data_list[c]),tres)
         
         # ctr_str = 'Difference'
         # if glb_avg_list != []: 
         #    glb_diff = glb_avg_list[c] - glb_avg_list[diff_base]
         #    ctr_str += f' ({glb_diff:6.4})'
         # hs.set_subtitles(wks, plot[ipd], name[c], ctr_str, var_str[v], font_height=subtitle_font_height)

         lstr = f'{name[c]} - {name[0]}'
         hs.set_subtitles(wks, plot[ipd], lstr, '', var_str[v], font_height=subtitle_font_height)

#---------------------------------------------------------------------------------------------------
# Finalize plot

num_case_alt = num_case*2-1 if (plot_diff and add_diff) else num_case
layout = [num_var,num_case_alt] if var_x_case else [num_case_alt,num_var]


if not (plot_diff and add_diff):
   if num_case==1 or num_var==1:
      layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
   
pnl_res = hs.setres_panel()
if use_common_label_bar: pnl_res.nglPanelLabelBar = True
### add panel labels
# pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
# pnl_res.nglPanelFigureStringsJust        = "TopLeft"
# pnl_res.nglPanelFigureStringsFontHeightF = 0.01
# if layout==[3,2] : pnl_res.nglPanelFigureStringsFontHeightF = 0.015

pnl_res.nglPanelYWhiteSpacePercent = 5

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
