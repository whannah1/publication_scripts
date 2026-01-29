#---------------------------------------------------------------------------------------------------
# Plot the zonal mean of the specified variables
#---------------------------------------------------------------------------------------------------
import os, ngl, xarray as xr, numpy as np, copy
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import cmocean
np.seterr(divide='ignore', invalid='ignore')
np.errstate(divide='ignore', invalid="ignore")
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
add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.00',           n='E3SM-MMF'    ,c='black')
# add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.VT_0.00',      n='E3SM-MMF+BVT',c='green')
# add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.VT_1.00',      n='E3SM-MMF+FVT',c='blue')
add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.00',n='E3SM-MMF+BVT',c='green')
add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_1.00',n='E3SM-MMF+FVT',c='blue')
#-------------------------------------------------------------------------------

add_var('PRECT')
# add_var('TMQ')
add_var('U850')
# add_var('TGCLDLWP')
# add_var('TGCLDIWP')
# add_var('LHFLX')
# add_var('SHFLX')
# add_var('FLNT')
# add_var('FSNT')


fig_type = 'png'
fig_file = 'figs/FXX-zonal-mean-variance'

dlat = 2

plot_diff = True

num_plot_col = 2

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var = len(var)
num_case = len(case)

wks = ngl.open_wks('png',fig_file)
plot = [None]*num_var*2
# plot = [None]*num_var
# if plot_diff: plot = [None]*num_var*2
res = hs.res_xy()
res.vpHeightF = 0.3
res.xyLineThicknessF = 8

if 'clr' not in vars(): 
   if num_case>1 : clr = np.linspace(2,len( ngl.retrieve_colormap(wks) )-1,num_case,dtype=int)
   else : clr = ['black']

# if num_case>1 and 'dsh' not in vars(): dsh = np.arange(0,num_case,1)
if 'dsh' not in vars(): 
   if num_case>1 : dsh = np.zeros(num_case)
   else : dsh = [0]
res.xyLineColors   = clr
res.xyDashPatterns = dsh

# res.tiXAxisString = 'Latitude'
res.tiXAxisString = 'sin( Latitude )'

class tcolor:
   ENDC,RED,GREEN,YELLOW,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[33m','\033[35m','\033[36m'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
msg_list = []
for v in range(num_var):
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   data_list_cli,data_list_var = [],[]
   if 'lev_list' in locals(): lev = lev_list[v]
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)
      case_obj = he.Case( name=case[c] )
      #-------------------------------------------------------------------------
      # read the data
      #-------------------------------------------------------------------------
      area = case_obj.load_data('area',htype='h0',num_files=1).astype(np.double)
      lat = case_obj.load_data('lat',  htype='h0',num_files=1)

      tmp_file_head_var = 'data/variance'
      tmp_file_head_cli = 'data/climatology'

      tmp_file_var = f'{tmp_file_head_var}.{case[c]}.{var[v]}.nc'
      tmp_file_cli = f'{tmp_file_head_cli}.{case[c]}.{var[v]}.nc'

      print(f'    tmp_files:')
      print(f'      {tmp_file_var}')
      print(f'      {tmp_file_cli}')

      ds_var = xr.open_dataset( tmp_file_var )
      ds_cli = xr.open_dataset( tmp_file_cli )
      
      data_var = ds_var[var[v]]
      data_cli = ds_cli[var[v]]

      msg = hc.print_stat(data_cli,name=var[v]+' mean    ',stat='naxsh',indent=(' '*6),compact=True)
      msg = hc.print_stat(data_var,name=var[v]+' variance',stat='naxsh',indent=(' '*6),compact=True)
      # msg_list.append('  case: '+case[c]+'\n'+msg)
      if 'area' in vars() :
         gbl_mean_cli = ( (data_cli*area).sum() / area.sum() ).values 
         gbl_mean_var = ( (data_var*area).sum() / area.sum() ).values 
         print(f'      Area Weighted Global Mean / Variance : {gbl_mean_cli:6.4} / {gbl_mean_var:6.4}')

      #-------------------------------------------------------------------------
      # Calculate time and zonal mean
      #-------------------------------------------------------------------------
      bin_ds_var = hc.bin_YbyX( data_var, lat, bin_min=-90, bin_max=90, bin_spc=dlat, wgt=area )
      bin_ds_cli = hc.bin_YbyX( data_cli, lat, bin_min=-90, bin_max=90, bin_spc=dlat, wgt=area )

      data_list_var.append( bin_ds_var['bin_val'].values )
      data_list_cli.append( bin_ds_cli['bin_val'].values )

      sin_bins = np.sin( bin_ds_cli['bins'].values *np.pi/180.)

   #----------------------------------------------------------------------------
   # Take difference from first case
   #----------------------------------------------------------------------------
   if plot_diff :
      diff_list_var = copy.deepcopy(data_list_var)
      # diff_list_cli = copy.deepcopy(data_list_cli)
      for c in range(num_case): 
         diff_list_var[c] = data_list_var[c] - data_list_var[0]
         # diff_list_cli[c] = data_list_cli[c] - data_list_cli[0]
      # data_list_var = diff_list_var
      # data_list_cli = diff_list_cli
   #----------------------------------------------------------------------------
   # set up plot stuff
   #----------------------------------------------------------------------------
   subtitle_font_height = 0.015
   unit_str,var_str = '',var[v]
   if var[v] =='PRECT'                      : unit_str = '[mm/day]'; var_str = 'Precipitation'
   if var[v] in ['PRECT','PRECC','PRECL']   : unit_str = '[mm/day]'
   if var[v] in ['LHFLX','SHFLX']           : unit_str = '[W/m2]'
   if var[v]=='TMQ': unit_str = '[mm]'

   res.tiYAxisString = unit_str

   res.trXMinF = -1. #np.min( sin_bins )
   res.trXMaxF =  1. #np.max( sin_bins )

   lat_tick = np.array([-90,-60,-30,0,30,60,90])
   res.tmXBMode = "Explicit"
   res.tmXBValues = np.sin( lat_tick*3.14159/180. )
   res.tmXBLabels = lat_tick

   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   if plot_diff:
      plot[v*2  ] = ngl.xy(wks, sin_bins, np.ma.masked_invalid(  np.stack(data_list_var) ), res)
      plot[v*2+1] = ngl.xy(wks, sin_bins, np.ma.masked_invalid(  np.stack(diff_list_var) ), res)
      hs.set_subtitles(wks, plot[v*2  ], "Variance","", var_str, font_height=subtitle_font_height)
      hs.set_subtitles(wks, plot[v*2+1], "Variance Difference", "", var_str, font_height=subtitle_font_height)
   else:
      plot[v*2  ] = ngl.xy(wks, sin_bins, np.ma.masked_invalid(  np.stack(data_list_cli) ), res)
      plot[v*2+1] = ngl.xy(wks, sin_bins, np.ma.masked_invalid(  np.stack(data_list_var) ), res)
      hs.set_subtitles(wks, plot[v*2  ], "Time Mean","", var_str, font_height=subtitle_font_height)
      hs.set_subtitles(wks, plot[v*2+1], "Variance of Daily Means", "", var_str, font_height=subtitle_font_height)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

ngl.panel(wks,plot,layout,hs.setres_panel())
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------