#---------------------------------------------------------------------------------------------------
# Plot the zonal mean of the specified variables
#---------------------------------------------------------------------------------------------------
import os, ngl, xarray as xr, numpy as np, copy, warnings, string
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

def add_var(var_name,lev=-1): var.append(var_name); #lev_list.append(lev)
#-------------------------------------------------------------------------------
add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.01',n='')
#-------------------------------------------------------------------------------

lev = np.array([10,30,50,75,100,125,150,200,250,300,350,400,450,500,
               550,600,650,700,750,800,825,850,875,900,925,950,975,1000])
# lev = np.array([5,10,30,50,100,150,200,300,400,500,600,700,800,850,925,975,1000])


num_plot_col = 3

# add_var('MMF_VT_T')
# add_var('MMF_VT_Q')
# add_var('MMF_VT_U')

add_var('MMF_VT_T')
add_var('MMF_VT_TLS')
add_var('MMF_VT_TEND_T')
if num_plot_col==4: add_var('MMF_VT_TNET')

add_var('MMF_VT_Q')
add_var('MMF_VT_QLS')
add_var('MMF_VT_TEND_Q')
if num_plot_col==4: add_var('MMF_VT_QNET')

add_var('MMF_VT_U')
add_var('MMF_VT_ULS')
add_var('MMF_VT_TEND_U')
if num_plot_col==4: add_var('MMF_VT_UNET')


htype,first_file,num_files = 'h0',0,10*12

fig_type = 'png'
fig_file = 'figs/F07-zonal-mean-vt-tracer'

temp_file_head = 'data/zonal-mean-VT'

print_stats = False

recalculate = False

lat1, lat2, dlat = -88., 88., 2

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var = len(var)
num_case = len(case)

wks = ngl.open_wks('png',fig_file)
plot = [None]*num_var

res = hs.res_contour_fill()

res.vpHeightF = 0.3

res.lbBottomMarginF = 0.3
res.lbTopMarginF    = 0.05
res.lbRightMarginF  = -0.2
res.lbLeftMarginF   = -0.2

res.trYReverse = True

res.tiXAxisString = 'sin( Latitude )'
res.tiYAxisString = 'Pressure [hPa]'

res.tiXAxisFontHeightF     = 0.02
res.tiYAxisFontHeightF     = 0.02
res.tmYLLabelFontHeightF   = 0.02
res.tmXBLabelFontHeightF   = 0.02
res.lbLabelFontHeightF     = 0.02
res.lbTitleFontHeightF     = 0.02
res.lbTitlePosition        = 'Bottom'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
msg_list = []
for v in range(num_var):
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   data_list,std_list,cnt_list = [],[],[]
   # if 'lev_list' in locals(): lev = lev_list[v]
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

      tmp_file = f'{temp_file_head}.{case[c]}.{var[v]}.nc'

      print('\n    tmp_file: '+tmp_file+'\n')

      if recalculate :
         #----------------------------------------------------------------------
         # read the data
         #----------------------------------------------------------------------
         if 'lon1' in vars() : case_obj.lon1 = lon1
         if 'lon2' in vars() : case_obj.lon2 = lon2
         if 'lev'  in vars() : case_obj.lev  = lev

         lat  = case_obj.load_data('lat', htype=htype,num_files=1,use_remap=use_remap,remap_str=remap_str)
         area = case_obj.load_data('area',htype=htype,num_files=1,use_remap=use_remap,remap_str=remap_str).astype(np.double)
         data = case_obj.load_data(var[v],htype=htype,first_file=first_file,num_files=num_files,
                                    lev=lev,use_remap=use_remap,remap_str=remap_str)

         # if 'NET' in var[v]:
         #    if 'VT_T' in var[v]:tvar1,tvar2 = 'MMF_VT_TEND_T','MMF_VT_TLS'
         #    if 'VT_Q' in var[v]:tvar1,tvar2 = 'MMF_VT_TEND_Q','MMF_VT_QLS'
         #    if 'VT_U' in var[v]:tvar1,tvar2 = 'MMF_VT_TEND_U','MMF_VT_ULS'
         #    data1 = case_obj.load_data(tvar1,htype=htype,first_file=first_file,num_files=num_files,
         #                               lev=lev,use_remap=use_remap,remap_str=remap_str)
         #    data2 = case_obj.load_data(tvar2,htype=htype,first_file=first_file,num_files=num_files,
         #                               lev=lev,use_remap=use_remap,remap_str=remap_str)
         #    data = data1 + data2
         # else:
         #    data = case_obj.load_data(var[v],htype=htype,first_file=first_file,num_files=num_files,
         #                               lev=lev,use_remap=use_remap,remap_str=remap_str)
         
         if print_stats:
            hc.print_time_length(data.time,indent=(' '*6))
            hc.print_stat(data,name=var[v],stat='naxsh',indent=(' '*6),compact=True)

            if 'area' in vars() :
               gbl_mean = ( (data*area).sum() / area.sum() ).values 
               print(f'      Area Weighted Global Mean : {gbl_mean:6.4}')

         #-------------------------------------------------------------------------
         # Calculate time and zonal mean
         #-------------------------------------------------------------------------
         with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)

            bin_ds = hc.bin_YbyX( data.mean(dim='time', skipna=True), lat, \
                                  bin_min=lat1, bin_max=lat2, \
                                  bin_spc=dlat, wgt=area, keep_lev=True )

         print('writing to file: '+tmp_file)
         bin_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         bin_ds = xr.open_dataset( tmp_file )

      sin_bins = np.sin(bin_ds['bins'].values*np.pi/180.)
      data_binned = np.ma.masked_invalid( bin_ds['bin_val'].transpose().values )
      data_list.append( data_binned )

   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)

   # set units for color bar
   unit_str = ''
   if var[v]=='MMF_VT_T'      : unit_str = 'K~S~2~N~'
   if var[v]=='MMF_VT_TLS'    : unit_str = 'K~S~2~N~/s'
   if var[v]=='MMF_VT_TEND_T' : unit_str = 'K~S~2~N~/s'
   if var[v]=='MMF_VT_Q'      : unit_str = 'kg~S~2~N~/kg~S~2~N~'
   if var[v]=='MMF_VT_QLS'    : unit_str = 'kg~S~2~N~/kg~S~2~N~/s'
   if var[v]=='MMF_VT_TEND_Q' : unit_str = 'kg~S~2~N~/kg~S~2~N~/s'
   if var[v]=='MMF_VT_U'      : unit_str = 'm~S~2~N~/s~S~2~N~'
   if var[v]=='MMF_VT_ULS'    : unit_str = 'm~S~2~N~/s~S~2~N~/s'
   if var[v]=='MMF_VT_TEND_U' : unit_str = 'm~S~2~N~/s~S~2~N~/s'
   tres.lbTitleString = unit_str

   tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
   if not( 'TEND' in var[v] or 'LS' in var[v] or 'NET' in var[v] ) : 
      tres.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )

   lat_tick = np.array([-90,-60,-30,0,30,60,90])
   tres.sfXArray = sin_bins
   tres.tmXBMode = "Explicit"
   tres.tmXBValues = np.sin( lat_tick*3.14159/180. )
   tres.tmXBLabels = lat_tick

   if v%num_plot_col!=0:
      tres.tiYAxisString = ''

   tres.sfYCStartV = np.min( lev )
   tres.sfYCEndV   = np.max( lev )

   data_min = np.min([np.min(d) for d in data_list])
   data_max = np.max([np.max(d) for d in data_list])

   aboutZero = False
   if 'TEND' in var[v]: aboutZero = True
   if 'LS'   in var[v]: aboutZero = True
   if 'NET'  in var[v]: aboutZero = True
   cmin,cmax,cint,clev = ngl.nice_cntr_levels(data_min, data_max, cint=None, max_steps=21, \
                                              returnLevels=True,aboutZero=aboutZero )
   tres.cnLevels = np.linspace(cmin,cmax,num=21)
   tres.cnLevelSelectionMode = "ExplicitLevels"

   ip = v # v*num_case+c if var_x_case else c*num_var+v

   plot[ip] = ngl.contour(wks, data_list[c], tres)

   #----------------------------------------------------------------------------
   # Set plot subtitles
   #----------------------------------------------------------------------------
   var_str = var[v]
   if 'MMF_VT' in var[v]: var_str = var_str.replace('MMF_VT','VT')
   if 'TEND' in var_str: var_str = var_str.replace('TEND_','')+'_TEND_CRM'
   if 'LS'   in var_str: var_str = var_str.replace('LS'   ,'')+'_TEND_GCM'
   if 'NET'  in var_str: var_str = var_str.replace('NET'  ,'')+'_TEND_NET'

   name_str = name[c]

   hs.set_subtitles(wks, plot[ip], name_str, '', var_str, font_height=0.010)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

pnl_res = hs.setres_panel()

pnl_res.nglPanelYWhiteSpacePercent = 0
pnl_res.nglPanelXWhiteSpacePercent = 2

pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.01

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

# print()
# for msg in msg_list: print(msg)

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
