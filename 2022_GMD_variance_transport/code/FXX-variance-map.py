import os, ngl, xarray as xr, numpy as np, copy
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
import cmocean
data_dir,data_sub = None,None
case,name,clr,dsh = [],[],[],[]
var,lev_list = [],[]
def add_case(case_in,n='',c='black',d=0):
   global name,case
   case.append(case_in); name.append(n); clr.append(c); dsh.append(d)

def add_var(var_name,lev=-1): var.append(var_name); lev_list.append(lev)
#-------------------------------------------------------------------------------
# add_case('MAC-PG',    n='MAC')
# add_case('GPM-PG',    n='GPM')
add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.00',           n='E3SM-MMF')
add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.00',n='E3SM-MMF+BVT')
add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_1.00',n='E3SM-MMF+FVT')
#-------------------------------------------------------------------------------

add_var('PRECT')
# add_var('TMQ')
# add_var('U850')
# add_var('TGCLDLWP')
# add_var('TGCLDIWP')

fig_type = 'png'
fig_file = 'figs/FXX-variance-map'

lat1,lat2 = -60,60

num_yr = 10
# htype,first_file,num_files = 'h1',0,365*num_yr ; num_files_obs = 12*num_yr
# use_remap,remap_grid = False,'90x180' # 90x180 / 180x360

plot_diff = True

num_plot_col = 2

tmp_file_head = f'data/variance'
# tmp_file_head = f'data/variance.nyr_{num_yr}'

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)

diff_base = 0

# if 'lev' not in vars(): lev = np.array([0])
if 'lev' not in vars(): lev = np.array([0]*num_var)
# if len(lev) != num_var : lev = [lev]*num_var

wks = ngl.open_wks(fig_type,fig_file)
plot = [None]*(num_var*num_case)
res = hs.res_contour_fill_map()
if 'lat1' in locals(): res.mpMinLatF = lat1
if 'lat2' in locals(): res.mpMaxLatF = lat2
if 'lon1' in locals(): res.mpMinLonF = lon1
if 'lon2' in locals(): res.mpMaxLonF = lon2

### use this for aquaplanets
# res.mpOutlineBoundarySets = "NoBoundaries"
# res.mpCenterLonF = 0.

res.tmYLLabelFontHeightF         = 0.008
res.tmXBLabelFontHeightF         = 0.008
res.lbLabelFontHeightF           = 0.012

if 'TRMM' in case :
   res.mpMinLatF = -50
   res.mpMaxLatF =  50

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
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
   data_list = []
   if 'lev_list' in locals(): lev = lev_list[v]
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)

      #-------------------------------------------------------------------------
      # read the pre-calculated data
      #-------------------------------------------------------------------------   
      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}.nc'
      print('    tmp_file: '+tmp_file+'')
   
      ds = xr.open_dataset( tmp_file )
      data = ds[var[v]]

      data_list.append( data.values )

      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      # data_dir_tmp,data_sub_tmp = None, None
      # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
      # if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      # if case_sub[c] is not None: data_sub_tmp = case_sub[c]
      
      # case_obj = he.Case( name=case[c], atm_comp=get_comp(case[c]), data_dir=data_dir_tmp, data_sub=data_sub_tmp )
      # case_obj.set_coord_names(var[v])

      case_obj = he.Case( name=get_name(case[c]), time_freq='daily' )
      case_obj.set_coord_names(var[v])
      
      remap_str = 'remap_ne30pg2'
      use_remap = False
      if case[c]=='MAC-PG' : use_remap = True; remap_str=f'remap_ne30pg2'
      if case[c]=='MAC-FV' : use_remap = False
      if case[c]=='GPM-PG' : use_remap = True; remap_str=f'remap_ne30pg2'
      if case[c]=='GPM-FV' : use_remap = True; remap_str=f'remap_180x360'

      if case[c] in ['MAC-FV']: case_obj.file_path_daily = f'{case_obj.data_dir}/{case_obj.name}/daily_cwp/*'
      if case[c] in ['MAC-PG']: case_obj.file_path_daily = f'{case_obj.data_dir}/{case_obj.name}/daily_cwp_ne30pg2/*'

      if case[c] in ['GPM-PG']: case_obj.file_path_daily = f'{case_obj.data_dir}/{case_obj.name}/daily_ne30pg2/*'
      if case[c] in ['GPM-FV']: case_obj.file_path_daily = f'{case_obj.data_dir}/{case_obj.name}/daily_180x360/*'

      #-------------------------------------------------------------------------
      # print global mean
      #-------------------------------------------------------------------------
      # print('    Loading area for global mean...')
      area = case_obj.load_data('area', component=get_comp(case[c]),
                                 htype='h1',num_files=1,
                                 use_remap=use_remap,remap_str=remap_str)
      gbl_mean = ( (data*area).sum() / area.sum() ).values 
      print('\n      Area Weighted Global Mean : '+'%f'%gbl_mean+'\n')

      #-------------------------------------------------------------------------
      # save baseline for diff map
      #-------------------------------------------------------------------------
      if plot_diff :
         if c==diff_base:
            data_baseline = data.copy()
   
   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])

   if plot_diff:
      tmp_data = data_list - data_list[diff_base]
      for c in range(num_case): tmp_data[c] = data_list[c] - data_list[diff_base]
      diff_data_min = np.min([np.nanmin(d) for d in tmp_data])
      diff_data_max = np.max([np.nanmax(d) for d in tmp_data])

   for c in range(num_case):
      # ip = v*num_case+c
      ip = c*num_var+v
      #-------------------------------------------------------------------------
      # Set colors and contour levels
      #-------------------------------------------------------------------------
      tres = copy.deepcopy(res)
      tres.cnFillPalette = "MPL_viridis"
      # if var[v] in ['PRECT','PRECC','PRECL']    : tres.cnFillPalette = "WhiteBlueGreenYellowRed"

      if var[v] in ['PRECT']        : tres.cnLevels = np.arange(10,1000+10,10)
      if var[v] in ['CRM_PREC']        : tres.cnLevels = np.arange(50,5000+50,50)
      if var[v]=="TGCLDLWP"         : tres.cnLevels = np.arange(0.005,0.1+0.005,0.005)
      
      if hasattr(tres,'cnLevels') : 
         tres.cnLevelSelectionMode = "ExplicitLevels"
      else:
         aboutZero = False
         clev_tup = ngl.nice_cntr_levels(data.min().values, data.max().values,       \
                                         cint=None, max_steps=21,              \
                                         returnLevels=False, aboutZero=aboutZero )
         if clev_tup==None: 
            tres.cnLevelSelectionMode = "AutomaticLevels"   
         else:
            cmin,cmax,cint = clev_tup
            tres.cnLevels = np.linspace(cmin,cmax,num=21)
            tres.cnLevelSelectionMode = "ExplicitLevels"

      var_str = var[v]
      if var[v]=="PRECT"      : var_str = "Precipitation Variance"
      if var[v]=="TGCLDLWP"   : var_str = "LWP Variance"

      # different color options for difference plot
      if plot_diff and c!=diff_base :
         # tres.cnFillPalette = 'BlueWhiteOrangeRed'
         # tres.cnFillPalette = np.array( cmocean.cm.diff(np.linspace(0,1,256)) )
         # tres.cnFillPalette = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )
         tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
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
      #-------------------------------------------------------------------------
      # Create plot
      #-------------------------------------------------------------------------
      hs.set_cell_fill(tres,case_obj)
   
      if 'name' in vars():
         case_name = name[c]
      else:
         case_name = case_obj.short_name

      # if case[c]=='TRMM' : 
      #    data = ngl.add_cyclic(data)
      # else:
      #    data = data.values
      
      if plot_diff and c!=diff_base: data_list[c] = data_list[c] - data_baseline.values

      plot[ip] = ngl.contour_map(wks,data_list[c],tres)

      ctr_str = ''
      if plot_diff: ctr_str = 'diff'
      # if var[v] in ['PRECT','PRECC','PRECL'] and 'gbl_mean' in vars() : 
      #    ctr_str = 'Mean: '+'%.2f'%gbl_mean+' [mm/day]'
      # hs.set_subtitles(wks, plot[len(plot)-1], case_name, ctr_str, var_str, font_height=0.01)
      hs.set_subtitles(wks, plot[ip], case_name, ctr_str, var_str, font_height=0.015)

      
#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

# layout = [len(plot),1]
# layout = [num_var,num_case]
layout = [num_case,num_var]

# if num_var==1  : layout = [int(np.ceil(num_case/2.)),2]
if num_var==1  : layout = [num_case,num_var]
if num_case==1 : layout = [num_var,num_case]

# if num_case==1 or num_var==1:
#    layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

# if num_var==1 and num_case==4 : layout = [2,2]
# if num_var==1 and num_case==6 : layout = [3,2]
# if num_case==1 and num_var==4 : layout = [2,2]

ngl.panel(wks,plot[0:len(plot)],layout,hs.setres_panel())
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------