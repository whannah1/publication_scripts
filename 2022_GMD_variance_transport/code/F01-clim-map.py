import os, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string
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
# add_case('GPM-PG',    n='IMERG')
# add_case('MAC-PG',    n='MAC')
add_case('OBS-PG', n='OBS')
add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.00',           n='E3SM-MMF')
add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.00',n='E3SM-MMF+BVT')
add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_1.00',n='E3SM-MMF+FVT')
#-------------------------------------------------------------------------------

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-v',dest='var',default=None,help='variable to plot')
(opts, args) = parser.parse_args()


if opts.var is None:
   add_var('PRECT')
   # add_var('TGCLDLWP')
   # add_var('TGCLDIWP')
   # add_var('TMQ')
   # add_var('U850')
else:
   add_var(opts.var)


fig_type = 'png'
fig_file = 'figs/F01-clim-map'

lat1,lat2 = -60,60
# lat1,lat2,lon1,lon2 = -30,40,90,260

htype,years,months,first_file,num_files = 'h0',[],[],0,10*12 ; num_files_obs = 365*10

# num_files_obs = 2 # test

use_remap,remap_grid = False,'180x360' # 90x180 / 180x360

plot_diff,add_diff,diff_base = False,False,0

recalculate = False

print_stats = True
var_x_case = False
num_plot_col = 2

common_colorbar = True

tmp_file_head = 'data/climatology'

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)

subtitle_font_height = 0.02 - 0.0012*num_var - 0.0014*(num_case+int(add_diff))
subtitle_font_height = max(subtitle_font_height,0.005)

if 'diff_case' not in vars(): diff_case = [(i+1) for i in range(num_case-1)]
if 'lev' not in vars(): lev = np.array([0])

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix

wks = ngl.open_wks(fig_type,fig_file,wkres)
if plot_diff and add_diff: 
   plot = [None]*(num_var*(num_case*2-1))
else:
   plot = [None]*(num_var*num_case)

anno = [None]*len(plot)
pdum = [None]*len(plot)
   
res = hs.res_contour_fill_map()
res.nglMaximize      = False
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

### resources for box around inset region
pgres = ngl.Resources()
pgres.nglDraw,pgres.nglFrame = False,False
pgres.gsLineColor = 'red'
pgres.gsLineThicknessF = 6

def get_comp(case):
   comp = 'eam'
   if 'OBS'  in case: comp = None
   if 'MAC'  in case: comp = None
   if 'GPM'  in case: comp = None
   if 'CESM' in case: comp = 'cam'
   return comp

def get_name(case,var):
   tmp_name = case
   if 'MAC'  in case: tmp_name = 'MAC'
   if 'GPM'  in case: tmp_name = 'GPM'
   if 'OBS' in case and var=='TGCLDLWP': tmp_name = 'MAC'
   if 'OBS' in case and var=='PRECT'   : tmp_name = 'GPM'
   return tmp_name

def is_obs(case):
   if 'MAC' in case: return True
   if 'GPM' in case: return True
   if 'OBS' in case: return True
   return False

def is_MAC_FV(case,var):
   return True if ('FV' in case and is_MAC(case,var)) else False
def is_MAC_PG(case,var):
   return True if ('PG' in case and is_MAC(case,var)) else False

def is_GPM_FV(case,var):
   return True if ('FV' in case and is_GPM(case,var)) else False
def is_GPM_PG(case,var): 
   return True if ('PG' in case and is_GPM(case,var)) else False

def is_MAC(case,var):
   return True if ('MAC' in case or ('OBS' in case and var=='TGCLDLWP')) else False

def is_GPM(case,var):
   return True if ('GPM' in case or ('OBS' in case and var=='PRECT'))else False
      
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   data_list,lat_list,lon_list = [],[],[]
   std_list,cnt_list = [],[]
   if 'lev_list' in locals(): lev = lev_list[v]
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)

      # data_sub_tmp = data_sub
      # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'

      if recalculate :
         # case_obj = he.Case( name=tname, atm_comp=get_comp(case[c]), data_dir=data_dir, data_sub=data_sub_tmp )
         case_obj = he.Case( name=get_name(case[c],var[v]), time_freq='daily' )
         case_obj.set_coord_names(var[v])
         remap_str = 'remap_ne30pg2'
         use_remap = False
         ddir = f'{case_obj.data_dir}/{case_obj.name}'
         if is_MAC_PG(case[c],var[v]): use_remap = True;remap_str=f'remap_ne30pg2';case_obj.file_path=f'{ddir}/daily_cwp_ne30pg2/*'
         if is_MAC_FV(case[c],var[v]): use_remap = False;                          case_obj.file_path=f'{ddir}/daily_cwp/*'
         if is_GPM_PG(case[c],var[v]): use_remap = True;remap_str=f'remap_ne30pg2';case_obj.file_path=f'{ddir}/daily_ne30pg2/*'
         if is_GPM_FV(case[c],var[v]): use_remap = True;remap_str=f'remap_180x360';case_obj.file_path=f'{ddir}/daily_180x360/*'
      
      #-------------------------------------------------------------------------
      # read the data
      #-------------------------------------------------------------------------   
      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}.nc'
      print('    tmp_file: '+tmp_file+'')

      if recalculate :

         num_files_tmp = num_files

         if is_obs(case[c]): num_files_tmp = num_files_obs
         
         tvar = var[v]
         if var[v]=='U850': tvar,lev = 'U',850
         if var[v]=='U200': tvar,lev = 'U',200

         data = case_obj.load_data(tvar, component=get_comp(case[c]),htype=htype,ps_htype=htype,
                                         years=years,months=months,lev=lev,
                                         first_file=first_file,num_files=num_files_tmp,
                                         use_remap=use_remap,remap_str=remap_str)

         if is_GPM(case[c],var[v]): data = data*24.   # for GPM convert from mm/hr => mm/day
         if is_MAC(case[c],var[v]): data = data/1e3   # For MAC convert from g/m2 => kg/m2

         # Get rid of lev dimension
         if 'lev' in data.dims : data = data.isel(lev=0)

         # combine lat/lon dimensions for plotting FV data
         if '-FV' in case[c]: 
            data = data.stack(ncol=("lat", "lon"))

         # print stats before time averaging
         if print_stats: hc.print_stat(data,name=var[v],stat='naxsh',indent='    ')

         # average over time dimension
         if 'time' in data.dims : 
            hc.print_time_length(data.time,indent=' '*6)
            data = data.mean(dim='time')
   
         # Calculate area weighted global mean
         if 'area' in locals() :
            gbl_mean = ( (data*area).sum() / area.sum() ).values 
            print(f'      Area Weighted Global Mean : {gbl_mean:6.4}')

         print('    writing to file: '+tmp_file)
         data.name = var[v]
         data.to_netcdf(path=tmp_file,mode='w')

         # load MAC total water path for later use to indicate suspicious data
         if is_MAC(case[c],var[v]):
            case_obj.file_path_daily = case_obj.file_path_daily.replace('cwp','twp')
            case_obj.file_path = case_obj.file_path_daily
            MAC_TWP = case_obj.load_data('TWP', component=get_comp(case[c]),htype=htype,ps_htype=htype,
                                            years=years,months=months,lev=lev,
                                            first_file=first_file,num_files=num_files_tmp,
                                            use_remap=use_remap,remap_str=remap_str)
            MAC_TWP = MAC_TWP.mean(dim='time')/1e3
            # MAC_TWP = xr.Dataset(MAC_TWP)
            MAC_TWP.name = var[v]
            tmp_file_twp = f'{tmp_file_head}.{case[c]}.TWP.nc'
            MAC_TWP.to_netcdf(path=tmp_file_twp,mode='w')

      else:
         ds = xr.open_dataset( tmp_file )
         data = ds[var[v]]

         if is_MAC(case[c],var[v]): 
            tmp_file_twp = f'{tmp_file_head}.{case[c]}.TWP.nc'
            ds_twp = xr.open_dataset( tmp_file_twp )
            MAC_TWP = ds_twp[var[v]]

      if is_MAC(case[c],var[v]): MAC_TWP = MAC_TWP.values
      #-------------------------------------------------------------------------
      # append to data lists
      #-------------------------------------------------------------------------
      if case[c]=='TRMM' and 'lon1' not in locals(): 
         data_list.append( ngl.add_cyclic(data.values) )
      else:
         data_list.append( data.values )

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
      if var[v] in ['PRECT','PRECC']   : tres.cnLevels = np.arange(1,20+1,1)
      # if var[v] in ['PRECT','PRECC']   : tres.cnLevels = np.logspace( -2, 1.31, num=60).round(decimals=2)
      # if var[v]=="TGCLDIWP"            : tres.cnLevels = np.logspace( -2, 0.25, num=20).round(decimals=2)
      # if var[v]=="TGCLDLWP"            : tres.cnLevels = np.logspace( -2, 0.25, num=20).round(decimals=2)
      if var[v]=="TGCLDIWP"            : tres.cnLevels = np.arange(0.005,0.155,0.01)
      if var[v]=="TGCLDLWP"            : tres.cnLevels = np.arange(0.005,0.25,0.005)
      if var[v]=="DYN_QLIQ"            : tres.cnLevels = np.logspace( -6, -4, num=40)

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
      if 'FV' in case[c]: scrip_file_path = 'scrip_files/cmip6_180x360_scrip.20210722.nc'

      scrip_ds = xr.open_dataset(scrip_file_path)

      tres.cnFillMode    = "CellFill"
      tres.sfXArray      = scrip_ds['grid_center_lon'].values
      tres.sfYArray      = scrip_ds['grid_center_lat'].values
      tres.sfXCellBounds = scrip_ds['grid_corner_lon'].values
      tres.sfYCellBounds = scrip_ds['grid_corner_lat'].values
      #-------------------------------------------------------------------------
      # Create plot
      #-------------------------------------------------------------------------
      
      tmp_data = data_list[c]

      # special considerations for MAC data
      if is_MAC(case[c],var[v]):
         # determine areas to add hatching - mean of CWP is small compared to total water path
         wp_ratio = tmp_data / MAC_TWP
         ones  = np.ones(len(tmp_data))
         zeros = np.zeros(len(tmp_data))
         hatch = np.where( wp_ratio<0.6, ones, zeros)
         # mask out land areas
         landfrac_ds = xr.open_dataset('data/E3SM_landfrac_ne30pg2.nc').isel(time=0)
         landfrac = landfrac_ds['LANDFRAC']
         ocn_mask = xr.DataArray( np.ones(len(tmp_data),dtype=bool), coords=[tmp_data] )
         ocn_mask = ocn_mask & (landfrac.values<0.5)
         tmp_data = np.ma.masked_array(tmp_data, mask=ocn_mask==0)
         # tmp_data = np.ma.masked_array(tmp_data, mask=tmp_data==np.nan) # also mask NaN
         
         
         
      if plot_diff and c==diff_base : base_name = name[c]

      if not plot_diff  or (plot_diff and add_diff) or (plot_diff and c==diff_base) : 

         plot[ip] = ngl.contour_map(wks,tmp_data,tres) 

         #----------------------------------------------------------------------
         # Hatching for MAC data
         #----------------------------------------------------------------------
         if is_MAC(case[c],var[v]):
            sres = hs.res_stippling()
            sres.sfXArray = tres.sfXArray
            sres.sfYArray = tres.sfYArray
            sres.cnFillScaleF = 1.
            sres.cnFillColor = 'white'
            # sres.cnFillPattern = np.array([-1,17]) # stippling
            sres.cnFillPatterns = np.array([-1,6])  # hatching
            ngl.overlay( plot[ip], ngl.contour(wks,hatch,sres) )
         #----------------------------------------------------------------------
         # Add inset to highlight region with checkerboard
         #----------------------------------------------------------------------
         # ires = hs.res_contour_fill_map()
         # ires.nglMaximize     = False
         # ires.tmXBOn          = False
         # ires.tmYLOn          = False
         # ires.vpHeightF       = 0.15
         # ires.vpWidthF        = 0.15
         # ires.cnFillPalette         = tres.cnFillPalette
         # ires.cnLevelSelectionMode  = tres.cnLevelSelectionMode
         # ires.cnLevels              = tres.cnLevels
         # ires.cnFillMode    = "CellFill"
         # ires.sfXArray      = scrip_ds['grid_center_lon'].values
         # ires.sfYArray      = scrip_ds['grid_center_lat'].values
         # ires.sfXCellBounds = scrip_ds['grid_corner_lon'].values
         # ires.sfYCellBounds = scrip_ds['grid_corner_lat'].values
         # ires.lbLabelBarOn = False
         # ires.mpLimitMode = "LatLon" 
         # # ires.mpMinLatF = -30
         # # ires.mpMaxLatF = 30
         # # ires.mpMinLonF = 150
         # # ires.mpMaxLonF = 180
         
         
         # iplot = ngl.contour_map(wks, tmp_data, ires)

         # # attach plot via annotation
         # ares = ngl.Resources()
         # ares.amZone = 1
         # ares.amOrthogonalPosF = 0.35
         # ares.amParallelPosF   = 0.75
         # anno = ngl.add_annotation(plot[ip], iplot ,ares)


         ires = copy.deepcopy(tres)
         
         ires.vpHeightF        = 0.06
         ires.vpWidthF         = 0.06
         
         ilat,ilon,idx = 15,162,6
         ires.mpMinLatF,ires.mpMaxLatF = ilat-idx,ilat+idx
         ires.mpMinLonF,ires.mpMaxLonF = ilon-idx,ilon+idx
         
         iplot = ngl.contour_map(wks, tmp_data, ires) 

         ### get NDC coordinates of where want to put the inset plot
         alat = -60
         alon = 360-20
         ndcx, ndcy = ngl.datatondc(plot[ip], alon, alat )

         ndcy = ndcy-0.2

         ### attach plot via annotation
         ares = ngl.Resources()
         ares.amZone = 1
         ares.amOrthogonalPosF,ares.amParallelPosF  = ndcy, ndcx
         anno[ip] = ngl.add_annotation(plot[ip], iplot ,ares)

         ### draw box around map region used for inset
         xbar = np.array([ilon-idx,ilon+idx,ilon+idx,ilon-idx,ilon-idx])
         ybar = np.array([ilat-idx,ilat-idx,ilat+idx,ilat+idx,ilat-idx])
         pdum[ip] = ngl.add_polyline(wks,plot[ip],xbar,ybar,pgres)

         #----------------------------------------------------------------------
         # set plot subtitles
         #----------------------------------------------------------------------
         if 'ctr_str' not in locals(): ctr_str = ''
         
         # ctr_str = 'Mean: '+'%.2f'%gbl_mean+' [mm/day]'

         # if 'lev' in locals() :
         #    if type(lev) not in [list,tuple]:
         #       if lev>0: ctr_str = f'{lev} mb'

         case_name = name[c]
         if is_MAC(case[c],var[v]): case_name = 'MAC'
         if is_GPM(case[c],var[v]): case_name = 'IMERG'

         hs.set_subtitles(wks, plot[ip], case_name, ctr_str, var_str, font_height=subtitle_font_height)

      #-------------------------------------------------------------------------
      # create difference plot
      #-------------------------------------------------------------------------
      if plot_diff and c in diff_case :
         
         tmp_data = tmp_data - data_baseline.values

         tres.cnFillPalette = 'BlueWhiteOrangeRed'
         tres.cnLevelSelectionMode = "ExplicitLevels"
         
         if hasattr(tres,'cnLevels') : del tres.cnLevels
         if var[v] in ['PRECT','PRECC','PRECL'] : tres.cnLevels = np.arange(-5,5+1,1)
         if not hasattr(tres,'cnLevels') : 
            if np.min(tmp_data)==np.max(tmp_data) : 
               print(hc.tcolor.RED+'WARNING: Difference is zero!'+hc.tcolor.ENDC)
            else:
               cmin,cmax,cint,clev = ngl.nice_cntr_levels(diff_data_min, diff_data_max,    \
                                                          cint=None, max_steps=21,      \
                                                          returnLevels=True, aboutZero=True )
               tres.cnLevels = np.linspace(cmin,cmax,num=21)
         
         ### override the level settings and just use auto
         # tres.cnLevelSelectionMode = "AutomaticLevels"

         tres.lbLabelBarOn = True

         ipd = ip
         if add_diff and     var_x_case: ipd = ip+1
         if add_diff and not var_x_case: ipd = ip+num_var*(num_case-1)

         plot[ipd] = ngl.contour_map(wks,tmp_data,tres)

         ctr_str = ''
         # case_name = name[c]+' - '+base_name
         case_name = name[c]
         if is_MAC(case[c],var[v]): case_name = 'MAC'
         if is_GPM(case[c],var[v]): case_name = 'IMERG'
         
         # hs.set_subtitles(wks, plot[ipd], case_name, '', var_str+' (Diff)', font_height=subtitle_font_height)
         hs.set_subtitles(wks, plot[ipd], case_name, 'Difference', var_str, font_height=subtitle_font_height)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
# if plot_diff : num_case = num_case+len(diff_case)   # use this to plot both before and after diff

num_case_alt = num_case*2-1 if (plot_diff and add_diff) else num_case
layout = [num_var,num_case_alt] if var_x_case else [num_case_alt,num_var]

if (not plot_diff) or (plot_diff and not add_diff):
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
