import os, ngl, subprocess as sp, numpy as np, xarray as xr
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
import copy, string
import cmocean
host = hc.get_host()
"""
commands to regrid GPCP data

ncremap -G ttl='GPCP grid 72x144'#latlon=72,144#lat_typ=uni#lon_typ=grn_wst \
-g $HOME/E3SM/data_grid/gpcp_72x144_scrip.nc

ncremap --alg_typ=aave \
--grd_src=$HOME/E3SM/data_grid/gpcp_72x144_scrip.nc \
--grd_dst=$HOME/E3SM/data_grid/ne30pg2_scrip.nc \
--map=$HOME/maps/map_72x144_to_ne30pg2_aave.nc

ncremap -m $HOME/maps/map_72x144_to_ne30pg2_aave.nc \
-i ~/Data/Obs/GPCP/precip.mon.ltm.1981-2010.nc \
-o ~/Data/Obs/GPCP/precip.mon.ltm.1981-2010.ne30pg2.nc

"""
#-------------------------------------------------------------------------------
name,case,case_dir,case_sub,case_grid = [],[],[],[],[]
def add_case(case_in,n=None,p=None,s=None,g=None,c=None):
   global name,case,case_dir,case_sub
   tmp_name = case_in if n is None else n
   case.append(case_in); name.append(tmp_name); case_dir.append(p); case_sub.append(s); case_grid.append(g)
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
# add_case('v2.LR.historical_0101',                                  n='E3SMv2',  c='cyan',    p=tmp_path_hst_v2, s='archive/atm/hist')
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

add_var('PRECT')

fig_file,fig_type = 'figs/FXX-map-precip-hist','png'
tmp_file_head = 'data/precip-hist'

# lat1,lat2 = -60,60

# htype,first_file,num_files = 'ha',31,2 # 1981-1986
htype,first_file,num_files = 'ha',31,30 # 1981-2010 (30 years)

use_remap,remap_grid = False,'90x180' # 90x180 / 180x360

plot_diff,add_diff = True,True

print_stats = True
var_x_case  = False
add_obs     = True

num_plot_col = 2#len(case)

use_common_label_bar = False

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
if case==[]: raise ValueError('ERROR - case list is empty!')
num_var,num_case = len(var),len(case)

subtitle_font_height = 0.015


if 'scrip_file_path' not in locals(): scrip_file_path = None

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

tmp_num_case = num_case
if add_obs: tmp_num_case += 1
if plot_diff and add_diff: tmp_num_case = tmp_num_case*2-1
plot = [None]*(num_var*tmp_num_case)
   
res = hs.res_contour_fill_map()
if 'lat1' in vars() : res.mpMinLatF = lat1; res.mpMaxLatF = lat2
if 'lon1' in vars() : res.mpMinLonF = lon1; res.mpMaxLonF = lon2

res.tmYLLabelFontHeightF         = 0.008
res.tmXBLabelFontHeightF         = 0.008
res.lbLabelFontHeightF           = 0.01
res.tmXBOn                       = False
res.tmYLOn                       = False
# res.mpGeophysicalLineColor       = 'white'

res.mpProjection = 'Mollweide'
res.pmTickMarkDisplayMode = "Never"
res.mpPerimOn = True

scripfile = xr.open_dataset(scrip_file_path)
# res.cnFillMode    = 'CellFill'
res.sfXArray      = scripfile['grid_center_lon'].rename({'grid_size':'ncol'}).values#.where( mask,drop=True).values
res.sfYArray      = scripfile['grid_center_lat'].rename({'grid_size':'ncol'}).values#.where( mask,drop=True).values
# res.sfXCellBounds = scripfile['grid_corner_lon'].rename({'grid_size':'ncol'}).values#.where( mask,drop=True).values 
# res.sfYCellBounds = scripfile['grid_corner_lat'].rename({'grid_size':'ncol'}).values#.where( mask,drop=True).values

#---------------------------------------------------------------------------------------------------
def get_tmp_file(case,var):
   return f'{tmp_file_head}.{case[c]}.nc'
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   data_list = []
   glb_avg_list = []
   if 'lev_list' in locals(): lev = lev_list[v]
   for c in range(num_case):
      tmp_file = get_tmp_file(case[c],var[v])
      print(' '*4+f'case: {hc.tclr.GREEN}{case[c]}{hc.tclr.END}  =>  {tmp_file}')

      if case[c]=='v2.LR.historical':
         v2_amip_ens_list = []
         v2_amip_ens_list.append('v2.LR.historical_0101')
         v2_amip_ens_list.append('v2.LR.historical_0151')
         v2_amip_ens_list.append('v2.LR.historical_0201')
         v2_amip_ens_list.append('v2.LR.historical_0251')
         v2_amip_ens_list.append('v2.LR.historical_0301')
         ens_cnt = 0
         for e,ens_member in enumerate(v2_amip_ens_list):
            tmp_file = get_tmp_file(ens_member,var[v])
            tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
            ens_member_data = tmp_ds[var[v]]
            if ens_cnt==0: 
               data = xr.zeros_like(ens_member_data)
            #    ens_min = ens_member_data.copy()
            #    ens_max = ens_member_data.copy()
            # else:
            #    for t in range(len(ens_member_data)):
            #       ens_min[t] = np.min([ens_min[t].values,ens_member_data[t].values])
            #       ens_max[t] = np.max([ens_max[t].values,ens_member_data[t].values])
            data = ( data*ens_cnt + ens_member_data ) / (ens_cnt+1)
            ens_cnt += 1
      else:
         #----------------------------------------------------------------------
         data_dir_tmp,data_sub_tmp = None, None
         if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
         if case_dir[c] is not None: data_dir_tmp = case_dir[c]
         if case_sub[c] is not None: data_sub_tmp = case_sub[c]
         
         case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp )
         case_obj.set_coord_names(var[v])
      
         # read the data
         tmp_first_file = first_file
         if 'v2.LR.historical' in case[c]: tmp_first_file = 50 + first_file

         area = case_obj.load_data('area',htype=htype)
         data = case_obj.load_data(var[v],htype=htype,ps_htype=htype,lev=lev,
                                         first_file=tmp_first_file,num_files=num_files,
                                         use_remap=use_remap,remap_str=f'remap_{remap_grid}')

         data['units'],data = 'mm/day',data*86400.*1e3

         # average over time dimension
         if 'time' in data.dims : 
            # hc.print_time_length(data.time,indent=' '*6)
            data = data.mean(dim='time')

         # print stats after time averaging
         if print_stats: hc.print_stat(data,name=var[v],stat='naxsh',indent='    ',compact=True)

         #----------------------------------------------------------------------
         # Calculate area weighted global mean
         if 'area' in locals() :
            gbl_mean = ( (data*area).sum() / area.sum() ).values 
            print(hc.tcolor.CYAN+f'      Area Weighted Global Mean : {gbl_mean:6.4}'+hc.tcolor.ENDC)
         #----------------------------------------------------------------------
         # write to file
         data.load()
         tmp_ds = xr.Dataset( coords=data.coords )
         tmp_ds[var[v]] = data
         tmp_ds['first_file'] = tmp_first_file
         tmp_ds['num_files']  = num_files
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      #-------------------------------------------------------------------------
      # append to data lists
      data_list.append( data.values )

   #------------------------------------------------------------------------------------------------
   # Add Observations
   if add_obs: 
      # obs_file = os.getenv('HOME')+'/Data/Obs/GPCP/precip.mon.mean.nc'
      obs_file = os.getenv('HOME')+'/Data/Obs/GPCP/precip.mon.ltm.1981-2010.ne30pg2.nc'
      ds = xr.open_dataset(obs_file)
      data_list.insert(0, ds['precip'].mean(dim='time').values )
      name.insert(0,'GPCP')
      case.insert(0,'GPCP')
      num_case = num_case+1

   #------------------------------------------------------------------------------------------------
   if plot_diff : data_baseline = data_list[0]

   #------------------------------------------------------------------------------------------------
   # Plot averaged data
   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])

   if plot_diff:
      tmp_data = [];
      for c in range(num_case): tmp_data.append( data_list[c] - data_baseline )
      diff_data_min = np.min([np.nanmin(d) for d in tmp_data])
      diff_data_max = np.max([np.nanmax(d) for d in tmp_data])

   for c in range(num_case):
      #-------------------------------------------------------------------------
      # Set colors
      #-------------------------------------------------------------------------
      tres = copy.deepcopy(res)
      # tres.cnFillPalette = "MPL_viridis"
      tres.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
      # tres.cnFillPalette = np.array( cmocean.cm.amp(np.linspace(0,1,256)) )

      #-------------------------------------------------------------------------
      # Set explicit contour levels
      #-------------------------------------------------------------------------
      if var[v] in ['PRECT','PRECC','PRECL']   : tres.cnLevels = np.arange(1,15+2,2)
      # if var[v] in ['PRECT','PRECC','PRECL']   : tres.cnLevels = np.arange(2,20+2,2)
      # if var[v] in ['PRECT','PRECC','PRECL']   : tres.cnLevels = np.arange(5,100+5,5) # for std dev
      # if var[v] in ['PRECT','PRECC','PRECL']   : tres.cnLevels = np.logspace( -2, 1.4, num=60).round(decimals=2)

      #-------------------------------------------------------------------------
      # set non-explicit contour levels
      #-------------------------------------------------------------------------
      if hasattr(tres,'cnLevels') : 
         tres.cnLevelSelectionMode = 'ExplicitLevels'
      else:
         nlev = 41
         aboutZero = False
         if var[v] in ['U','V'] : 
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
      # set unit string for colorbar
      unit_str = None
      if var[v]=='PRECT': unit_str = 'mm/day'
      if unit_str is not None: tres.lbTitleString = f'[{unit_str}]'
      #-------------------------------------------------------------------------
      # Create plot
      #-------------------------------------------------------------------------
      if use_common_label_bar: 
         tres.lbLabelBarOn = False
      else:
         tres.lbLabelBarOn = True
         tres.lbTitlePosition = 'Bottom'
         tres.lbTitleFontHeightF = 0.015
         tres.lbLabelFontHeightF = 0.015
         
         
      if plot_diff and c==0 : base_name = name[c]

      num_case_alt = num_case*2-1 if (plot_diff and add_diff) else num_case
      ip = v*num_case_alt+c if var_x_case else c*num_var+v

      if not plot_diff  or (plot_diff and add_diff) or (plot_diff and c==0) : 

         plot[ip] = ngl.contour_map(wks,np.ma.masked_invalid(data_list[c]),tres) 
         
         #----------------------------------------------------------------------
         # set plot subtitles
         #----------------------------------------------------------------------
         
         # ctr_str = 'Mean: '+'%.2f'%gbl_mean+' [mm/day]'

         ctr_str = f'{glb_avg_list[c]:6.4}' if glb_avg_list != [] else ''

         hs.set_subtitles(wks, plot[ip], name[c], ctr_str, var_str[v], font_height=subtitle_font_height)

      #-------------------------------------------------------------------------
      # create difference plot
      #-------------------------------------------------------------------------
      if plot_diff and c>0 :
         
         data_list[c] = data_list[c] - data_baseline

         # tres.cnFillPalette = "MPL_viridis"
         # tres.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
         # tres.cnFillPalette = np.array( cmocean.cm.amp(np.linspace(0,1,256)) )
         tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )

         tres.cnLevelSelectionMode = "ExplicitLevels"
         
         if hasattr(tres,'cnLevels') : del tres.cnLevels
         if var[v] in ['PRECT','PRECC','PRECL'] : tres.cnLevels = np.arange(-5,5+1,1)
         if not hasattr(tres,'cnLevels') : 
            if np.min(data_list[c])==np.max(data_list[c]) : 
               print(hc.tcolor.RED+'WARNING: Difference is zero!'+hc.tcolor.ENDC)
            else:
               cmin,cmax,cint,clev = ngl.nice_cntr_levels(diff_data_min, diff_data_max,    \
                                                          cint=None, max_steps=21,      \
                                                          returnLevels=True, aboutZero=True )
               tres.cnLevels = np.linspace(cmin,cmax,num=21)
         

         # if use_common_label_bar: 
         #    tres.lbLabelBarOn = False
         # else:
         #    tres.lbLabelBarOn = True
         tres.lbLabelBarOn = True

         ipd = ip
         if add_diff and     var_x_case: ipd = ip+(num_case-1)
         if add_diff and not var_x_case: ipd = ip+num_var*(num_case-1)

         # plot[ipd] = ngl.contour_map(wks,data_list[c],tres)
         plot[ipd] = ngl.contour_map(wks,np.ma.masked_invalid(data_list[c]),tres)
         
         ctr_str = ''
         # case_name = name[c]+' - '+base_name
         if 'name' in vars():
            case_name = name[c]
         else:
            case_name = case_obj.short_name

         # ctr_str = 'Diff'
         ctr_str = 'Bias wrt GPCP'
         # if glb_avg_list != []: 
         #    glb_diff = glb_avg_list[c] - glb_avg_list[0]
         #    ctr_str += f' ({glb_diff:6.4})'
         
         # hs.set_subtitles(wks, plot[ipd], case_name, '', var_str+' (Diff)', font_height=subtitle_font_height)
         hs.set_subtitles(wks, plot[ipd], case_name, ctr_str, var_str[v], font_height=subtitle_font_height)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

num_plot_col = 1
layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

pnl_res = hs.setres_panel()

### use common label bar
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
