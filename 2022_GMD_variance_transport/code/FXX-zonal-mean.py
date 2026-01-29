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
# add_case('GPM-PG',    n='IMERG')
# add_case('MAC-PG',    n='MAC')
add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.00',           n='E3SM-MMF'    ,c='black')
add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.00',n='E3SM-MMF+BVT',c='green')
add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_1.00',n='E3SM-MMF+FVT',c='blue')
#-------------------------------------------------------------------------------

# add_var('PRECT')
# add_var('TMQ')
# add_var('LHFLX')
# add_var('SHFLX')
# add_var('TGCLDLWP')
# add_var('TGCLDIWP')
# add_var('FLNT')
# add_var('FSNT')
add_var('NET_TOA_RAD')
# add_var('FLNS'); add_var('FSNS')


# htype,first_file,num_files = 'h2',5,5
htype,first_file,num_files = 'h0',0,10*12
# htype,first_file,num_files = 'h0',0,2*12
# htype,first_file,num_files = 'h0',0,2 # for debugging


fig_type = 'png'
fig_file = 'figs/F02-zonal-mean'

dlat = 2

plot_diff = False

print_stats = True

num_plot_col = 2

if plot_diff : num_plot_col = 2
#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var = len(var)
num_case = len(case)

wks = ngl.open_wks('png',fig_file)
plot = [None]*num_var
if plot_diff: plot = [None]*num_var*2
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
   data_list,std_list,cnt_list = [],[],[]
   if 'lev_list' in locals(): lev = lev_list[v]
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)

      tname = case[c]
      if 'MAC'  in case[c]: tname = 'MAC'
      if 'GPM'  in case[c]: tname = 'GPM'
      case_obj = he.Case( name=tname, time_freq='daily' )

      use_remap = False
      remap_str=f'remap_ne30pg2'
      if case[c]=='MAC-PG' : use_remap = True
      if case[c]=='GPM-PG' : use_remap = True

      if case[c] in ['MAC-PG']: case_obj.file_path = f'{case_obj.data_dir}/{case_obj.name}/daily_cwp_ne30pg2/*'
      if case[c] in ['GPM-PG']: case_obj.file_path = f'{case_obj.data_dir}/{case_obj.name}/daily_ne30pg2/*'
      #-------------------------------------------------------------------------
      # read the data
      #-------------------------------------------------------------------------
      if 'lon1' in vars() : case_obj.lon1 = lon1
      if 'lon2' in vars() : case_obj.lon2 = lon2
      if 'lev'  in vars() : case_obj.lev  = lev

      lat  = case_obj.load_data('lat', htype=htype,num_files=1,use_remap=use_remap,remap_str=remap_str)
      area = case_obj.load_data('area',htype=htype,num_files=1,use_remap=use_remap,remap_str=remap_str).astype(np.double)
      data = case_obj.load_data(var[v],htype=htype,num_files=num_files,use_remap=use_remap,remap_str=remap_str)

      # print()
      # print(data)
      # exit()

      hc.print_time_length(data.time,indent=(' '*6))

      if print_stats: 
         msg = hc.print_stat(data,name=var[v],stat='naxsh',indent=(' '*6),compact=True)
         msg_list.append('  case: '+case[c]+'\n'+msg)
         if 'area' in vars() :
            gbl_mean = ( (data*area).sum() / area.sum() ).values 
            print(f'      Area Weighted Global Mean : {gbl_mean:6.4}')

      #-------------------------------------------------------------------------
      # Calculate time and zonal mean
      #-------------------------------------------------------------------------
      # bin_ds = hc.bin_YbyX( data.mean(dim='time'), lat, lat_bins, bin_mode="explicit", wgt=area )
      # bin_ds = hc.bin_YbyX( data.mean(dim='time'), lat, bin_min=-88, bin_max=88, bin_spc=dlat, wgt=area )
      bin_ds = hc.bin_YbyX( data.mean(dim='time'), lat, bin_min=-90, bin_max=90, bin_spc=dlat, wgt=area )

      data_list.append( bin_ds['bin_val'].values )
      
      std_list.append( bin_ds['bin_std'].values )
      cnt_list.append( bin_ds['bin_cnt'].values )

      lat_bins = bin_ds['bins'].values

      sin_lat_bins = np.sin(lat_bins*np.pi/180.)

   #----------------------------------------------------------------------------
   # Take difference from first case
   #----------------------------------------------------------------------------
   if plot_diff :
      diff_list = copy.deepcopy(data_list)
      for c in range(num_case): diff_list[c] = data_list[c] - data_list[0]
   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   subtitle_font_height = 0.015
   unit_str,var_str = '',var[v]
   if var[v] =='PRECT'                      : unit_str = '[mm/day]'; var_str = 'Precipitation'
   if var[v] in ['PRECT','PRECC','PRECL']   : unit_str = '[mm/day]'
   if var[v] in ['LHFLX','SHFLX']           : unit_str = '[W/m2]'
   if var[v]=='TMQ'                         : unit_str = '[mm]'

   if var[v]=='TGCLDLWP' : unit_str = '[kg/m2]'; var_str = 'Liq Water Path'
   if var[v]=='TGCLDIWP' : unit_str = '[kg/m2]'; var_str = 'Ice Water Path'
   if var[v]=='FLNT'     : unit_str = '[W/m2]' ; var_str = 'TOA Net Longwave Radiation'
   if var[v]=='FSNT'     : unit_str = '[W/m2]' ; var_str = 'TOA Net Shortwave Radiation'

   res.tiYAxisString = unit_str

   res.trXMinF = -1. #np.min( sin_lat_bins )
   res.trXMaxF =  1. #np.max( sin_lat_bins )

   lat_tick = np.array([-90,-60,-30,0,30,60,90])
   res.tmXBMode = "Explicit"
   res.tmXBValues = np.sin( lat_tick*3.14159/180. )
   res.tmXBLabels = lat_tick

   tplot = ngl.xy(wks, sin_lat_bins, np.ma.masked_invalid(  np.stack(data_list) ), res)
   if plot_diff :
      plot[v*2  ] = tplot
      plot[v*2+1] = ngl.xy(wks, sin_lat_bins, np.ma.masked_invalid(  np.stack(diff_list) ), res)
      hs.set_subtitles(wks, plot[v*2  ], "", "", var_str, font_height=subtitle_font_height)
      hs.set_subtitles(wks, plot[v*2+1], "", "Difference", var_str, font_height=subtitle_font_height)
   else:
      plot[v] = tplot
      hs.set_subtitles(wks, plot[v], "", "", var_str, font_height=subtitle_font_height)


#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
# layout = [num_var,1]
# layout = [np.ceil(len(plot)/2),2]
# layout = [1,num_var]
# if num_var==4 : layout = [2,2]
# if num_var==6 : layout = [3,2]

layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

ngl.panel(wks,plot,layout,hs.setres_panel())
ngl.end()

# print()
# for msg in msg_list: print(msg)

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------