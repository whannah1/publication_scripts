import os, ngl, subprocess as sp
import numpy as np, xarray as xr
import hapy_common as hc
import hapy_E3SM   as he
import hapy_setres as hs
import copy, string
data_dir,data_sub = None,None
# data_dir,data_sub = '/project/projectdirs/m3312/jlee1046/E3SM/',''
#-------------------------------------------------------------------------------

### physgrid validation
name,case,git_hash = [],[],'cbe53b'
# case.append(f'E3SM.PGVAL.ne30_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}')
# case.append(f'E3SM.PGVAL.ne30pg2_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}')
# case.append(f'E3SM.PGVAL.ne30pg3_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}')
# case.append(f'E3SM.PGVAL.ne30pg4_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}')
case.append(f'E3SM.PGVAL.conusx4v1_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}')
case.append(f'E3SM.PGVAL.conusx4v1pg2_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}')


for c in case:
   if 'E3SM.PGVAL.ne30_r05_oECv3'         in c: name.append('ne30np4')
   if 'E3SM.PGVAL.ne30pg2_r05_oECv3'      in c: name.append('ne30pg2')
   if 'E3SM.PGVAL.ne30pg3_r05_oECv3'      in c: name.append('ne30pg3')
   if 'E3SM.PGVAL.ne30pg4_r05_oECv3'      in c: name.append('ne30pg4')
   if 'E3SM.PGVAL.ne30pg4_r05_oECv3'      in c: name.append('ne30pg4')
   if 'E3SM.PGVAL.conusx4v1_r05_oECv3'    in c: name.append('RRM np4')
   if 'E3SM.PGVAL.conusx4v1pg2_r05_oECv3' in c: name.append('RRM pg2')


#-------------------------------------------------------------------------------

var,lev_list = [],[]
def add_var(var_name,lev=-1): var.append(var_name); lev_list.append(lev)

add_var('PRECT')
add_var('TGCLDLWP')
add_var('TMQ')
add_var('U',lev=500)


fig_type = 'png'
fig_file = 'figs/F09-clim-map-RRM'

htype,years,months,first_file,num_files = 'h0',[],[],0,5*12

use_remap,remap_grid = True,'180x360' # 90x180 / 180x360

plot_diff,add_diff = True,False
diff_base = 0

use_snapshot,ss_t = False,28
chk_significance  = True
# t_crit=1.96  # 2-tail test w/ inf dof & P=0.05
t_crit=2.571  # 2-tail test w/ 5 dof & P=0.05

print_stats = True

var_x_case = True

scrip_dir = None # './scrip_files'

write_file = False
temp_dir = 'data/climatology'

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
class tcolor:
   ENDC,RED,GREEN,YELLOW,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[33m','\033[35m','\033[36m'

num_var = len(var)
num_case = len(case)

subtitle_font_height = 0.02 - 0.0015*num_var - 0.0014*(num_case+int(add_diff)*2.)

if 'diff_case' not in vars(): diff_case = [(i+1) for i in range(num_case-1)]
if 'lev' not in vars(): lev = np.array([0])

wkres = ngl.Resources()
npix=2048
# npix=4096
wkres.wkWidth  = npix
wkres.wkHeight = npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
if plot_diff and add_diff: 
   plot = [None]*(num_var*(num_case*2-1))
else:
   plot = [None]*(num_var*num_case)
   
res = hs.res_contour_fill_map()
if 'lat1' in vars() : res.mpMinLatF = lat1
if 'lat2' in vars() : res.mpMaxLatF = lat2
if 'lon1' in vars() : res.mpMinLonF = lon1
if 'lon2' in vars() : res.mpMaxLonF = lon2


res.tmYLLabelFontHeightF         = 0.008
res.tmXBLabelFontHeightF         = 0.008
res.lbLabelFontHeightF           = 0.018

# res.mpGeophysicalLineColor = 'white'

res.tmXBOn = False
res.tmYLOn = False

# res.mpLimitMode = "LatLon" 
# res.mpMinLatF   = 45 -15
# res.mpMaxLatF   = 45 +15
# res.mpMinLonF   = 180-15
# res.mpMaxLonF   = 180+15

if 'PGVAL' in case[0]:
   res.lbLabelFontHeightF  = 0.022
   subtitle_font_height    = 0.02

def get_comp(case):
   comp = 'cam'
   if 'CESM' in case: comp = 'cam'
   return comp

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+tcolor.MAGENTA+var[v]+tcolor.ENDC)
   data_list,area_list,lat_list,lon_list = [],[],[],[]
   std_list,cnt_list = [],[]
   if 'lev_list' in locals(): lev = lev_list[v]
   for c in range(num_case):
      print('    case: '+tcolor.GREEN+case[c]+tcolor.ENDC)

      tmp_file = f'{temp_dir}/clim_map_{case[c]}_{var[v]}.nc'

      if write_file:

         data_sub_tmp = data_sub
         if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
            
         case_obj = he.Case( name=case[c], atm_comp=get_comp(case[c]), data_dir=data_dir, data_sub=data_sub_tmp )

         # case_obj.time_freq = 'daily'
         # case_obj.time_freq = 'monthly'
         # case_obj = he.Case( name=case[c], data_dir='/global/cscratch1/sd/hewang/acme_scratch/cori-knl/' )

         case_obj.set_coord_names(var[v])
         
         #-------------------------------------------------------------------------
         # read the data
         #-------------------------------------------------------------------------   

         if use_remap or ('CESM' in case[c] and 'ne30' not in case[c]):
            # lat = case_obj.load_data('lat',component='cam',htype=htype,use_remap=use_remap,remap_str=f'remap_{remap_grid}')
            # lon = case_obj.load_data('lon',component='cam',htype=htype,use_remap=use_remap,remap_str=f'remap_{remap_grid}')
            lat = case_obj.load_data('lat',component=get_comp(case[c]),htype=htype,num_files=1,use_remap=use_remap,remap_str=f'remap_{remap_grid}')
            lon = case_obj.load_data('lon',component=get_comp(case[c]),htype=htype,num_files=1,use_remap=use_remap,remap_str=f'remap_{remap_grid}')
         else:
            aname = case_obj.area_name
            area = case_obj.load_data(aname,component=get_comp(case[c]),htype=htype,num_files=1)

         data = case_obj.load_data(var[v],component=get_comp(case[c]),
                                   htype=htype,ps_htype=htype,
                                   years=years,months=months,lev=lev,
                                   first_file=first_file,num_files=num_files,
                                   use_remap=use_remap,remap_str=f'remap_{remap_grid}')

         # Get rid of lev dimension
         if 'lev' in data.dims : data = data.isel(lev=0)

         if print_stats: hc.print_stat(data,name=var[v],stat='naxsh',indent='    ')

         # average over time dimension
         if 'time' in data.dims : 
            hc.print_time_length(data.time,indent=' '*6)
            if use_snapshot:
               data = data.isel(time=ss_t)
               print(hc.tcolor.RED+'WARNING - plotting snapshot!!!'+hc.tcolor.ENDC)
            else:
               if chk_significance :
                  std_list.append( data.std(dim='time').values )
                  cnt_list.append( float(len(data.time)) )
               data = data.mean(dim='time')

         #-------------------------------------------------------------------------
         #-------------------------------------------------------------------------

         # Calculate area weighted global mean
         if 'area' in locals() :
            gbl_mean = ( (data*area).sum() / area.sum() ).values 
            print(f'      Area Weighted Global Mean : {gbl_mean:6.4}')

         # # Calculate RMSE
         # if c==0:baseline = data
         # if c>0:
         #    rmse = np.sqrt( np.mean( np.square( data.to_masked_array() - baseline.to_masked_array() )))
         #    print(f'      Root Mean Square Error    : {rmse:6.4}')

         if case[c]=='TRMM' and 'lon1' not in locals(): 
            data_list.append( ngl.add_cyclic(data.values) )
         else:
            data_list.append( data.values )

         if 'area' in locals() : area_list.append( area.values )

         lat_list.append(lat.values)
         lon_list.append(lon.values)

         #-------------------------------------------------------------------------
         # Write to file
         #-------------------------------------------------------------------------
         tmp_ds = xr.Dataset()
         if use_remap:
            tmp_ds['data'] = xr.DataArray(data_list[c], dims=('dim_0','dim_1'))
            tmp_ds[ 'std'] = xr.DataArray( std_list[c], dims=('dim_0','dim_1'))
            tmp_ds[ 'cnt'] = xr.DataArray( cnt_list[c])
            tmp_ds[ 'lat'] = xr.DataArray( lat_list[c], dims=('dim_0'))
            tmp_ds[ 'lon'] = xr.DataArray( lon_list[c], dims=('dim_1'))
         else:
            tmp_ds['data'] = xr.DataArray(data_list[c])
            tmp_ds[ 'std'] = xr.DataArray( std_list[c])
            tmp_ds[ 'cnt'] = xr.DataArray( cnt_list[c])
         tmp_ds.to_netcdf(path=tmp_file,mode='w')

      else:

         tmp_ds = xr.open_dataset( tmp_file )

         data_list.append(tmp_ds['data'].values)
         std_list.append( tmp_ds[ 'std'].values)
         cnt_list.append( tmp_ds[ 'cnt'].values)
         lat_list.append( tmp_ds[ 'lat'].values)
         lon_list.append( tmp_ds[ 'lon'].values)

   #------------------------------------------------------------------------------------------------
   # Check significance using t-test - https://stattrek.com/hypothesis-test/difference-in-means.aspx
   #------------------------------------------------------------------------------------------------
   if plot_diff and chk_significance :
      print('\n  Checking significance of differences...')
      sig_list = []
      for c in range(1,num_case):
         t_stat = hc.calc_t_stat(D0=data_list[0],D1=data_list[c],
                                 S0=std_list[0], S1=std_list[c],
                                 N0=cnt_list[0], N1=cnt_list[c],t_crit=t_crit)

         sig_list.append( np.where( np.absolute(t_stat) > t_crit, 1, 0 ) )

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
      case_obj = he.Case( name=case[c], atm_comp=get_comp(case[c]), populate_files=False,
                          data_dir=data_dir, data_sub=data_sub )
      case_obj.set_coord_names(var[v])

      num_case_alt = num_case*2-1 if (plot_diff and add_diff) else num_case
      ip = v*num_case_alt+c if var_x_case else c*num_var+v
      #-------------------------------------------------------------------------
      # Set colors and contour levels
      #-------------------------------------------------------------------------
      tres = copy.deepcopy(res)
      tres.cnFillPalette = "MPL_viridis"
      if var[v] in ['P-E']                      : tres.cnFillPalette = "BlueWhiteOrangeRed"
      if var[v] in ['CLDLOW','CLDMED','CLDHGH'] : tres.cnFillPalette = "CBR_wet"
      if var[v] in ['TGCLDLWP','TGCLDIWP']      : tres.cnFillPalette = "MPL_viridis"
      if var[v] in ['DYN_QLIQ']                 : tres.cnFillPalette = "MPL_viridis"
      if var[v] in ['TS']                       : tres.cnFillPalette = "BlueRed"
      if var[v] in ['U','V','UBOT','VBOT','U850','V850','U200','V200',
                    'MMF_CVT_TEND_T','MMF_CVT_TEND_Q']: 
         tres.cnFillPalette = "BlueWhiteOrangeRed"
         # tres.cnFillPalette = 'MPL_RdYlBu'


      if var[v] in ['PRECT','PRECC']   : tres.cnLevels = np.arange(2,20+2,2)
      # if var[v] in ['PRECT','PRECC']   : tres.cnLevels = np.logspace( -2, 1.31, num=60).round(decimals=2)
      # if var[v]=='LHFLX'               : tres.cnLevels = np.arange(5,205+5,5)
      if var[v]=='P-E'                 : tres.cnLevels = np.linspace(-10,10,21)
      # if var[v]=="TS"                  : tres.cnLevels = np.arange(0,40+2,2)
      if var[v]=="RH"                  : tres.cnLevels = np.arange(10,100+1,1)
      if var[v]=="TGCLDIWP"            : tres.cnLevels = np.arange(1,30+1,1)*1e-2
      if var[v]=="TGCLDLWP"            : tres.cnLevels = np.logspace( -2, 0.25, num=60).round(decimals=2)
      if var[v]=="DYN_QLIQ"            : tres.cnLevels = np.logspace( -6, -4, num=40)
      if var[v]=="CLDLIQ"              : tres.cnLevels = np.logspace( -7, -3, num=60)
      if var[v]=='SPTLS'               : tres.cnLevels = np.linspace(-0.001,0.001,81)
      if var[v]=="PHIS"                : tres.cnLevels = np.arange(-0.01,0.01+0.001,0.001)
      if 'MMF_MC' in var[v]            : tres.cnLevels = np.linspace(-1,1,21)*0.01
      if 'OMEGA' in var[v]             : tres.cnLevels = np.linspace(-1,1,21)*0.6

      if var[v] in ['MMF_DU','MMF_DV','ZMMTU','ZMMTV','uten_Cu','vten_Cu']:
         tres.cnLevels = np.arange(-20,20+1,1)

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

      var_str = var[v]
      if var[v]=='PRECT':     var_str = 'Precipitation'
      if var[v]=='TMQ':       var_str = 'CWV'
      if var[v]=='TGCLDLWP':  var_str = 'LWP'
      if var[v]=='TGCLDIWP':  var_str = 'IWP'
      # if 'lev' in locals() and  var_str = f'{lev} mb {var[v]}'
      
      if 'PGVAL' in case[c]:
         # if var[v]=='TMQ':       var_str = 'Col. Water Vapor'
         # if var[v]=='TGCLDLWP':  var_str = 'Liq. Water Path'
         # if var[v]=='TGCLDIWP':  var_str = 'Ice Water Path'
         if var[v]=='PRECT':     var_str = 'Precip'
         if var[v]=='TMQ':       var_str = 'CWV'
         if var[v]=='TGCLDLWP':  var_str = 'LWP'
         if var[v]=='TGCLDIWP':  var_str = 'IWP'
      

      lev_str = None
      if lev>0: lev_str = f'{lev}mb'
      if lev<0: lev_str = f'k={(lev*-1)}'
      if lev_str is not None and var[v] in ['U','V','OMEGA','T','Q','Z3']:
         # var_str = f'{lev_str} {var[v]}'
         var_str = f'{var[v]}{lev}'
      # if var[v]=='U':         var_str = f'{lev_str} Zonal Wind'
      # if var[v]=='V':         var_str = f'{lev_str} Meridional Wind'
      # if var[v]=='OMEGA':     var_str = f'{lev_str} Omega'

      if 'RCE' in case[c] or 'AQP' in case[c] :
         tres.mpOutlineBoundarySets = 'NoBoundaries'
         tres.mpCenterLonF = 0.
      #-------------------------------------------------------------------------
      # Create plot
      #-------------------------------------------------------------------------
      if use_remap :
         hs.set_cell_fill(tres,case_obj=case_obj,lat=lat_list[c],lon=lon_list[c])
      else:
         hs.set_cell_fill(tres,case_obj=case_obj,htype=htype,scrip_dir=scrip_dir)

      if plot_diff : 
         if c==diff_base : base_name = name[c]

      if not plot_diff  or (plot_diff and add_diff) or (plot_diff and c==diff_base) : 

         plot[ip] = ngl.contour_map(wks,data_list[c],tres) 

         # scripfile = case_obj.get_scrip()
         
         #----------------------------------------------------------------------
         #----------------------------------------------------------------------
         ctr_str = ''
         if 'name' not in vars(): case_name = case_obj.short_name
         if 'name'     in vars(): case_name = name[c]
            
         # ctr_str = 'Mean: '+'%.2f'%gbl_mean+' [mm/day]'

         # if 'lev' in locals() :
         #    if type(lev) not in [list,tuple]:
         #       if lev>0: ctr_str = f'{lev} mb'
         
         hs.set_subtitles(wks, plot[ip], case_name, ctr_str, var_str, font_height=subtitle_font_height)

      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      if plot_diff and c in diff_case :
         
         data_list[c] = data_list[c] - data_list[diff_base]

         # tres.cnFillPalette = 'ncl_default'
         tres.cnFillPalette = 'BlueWhiteOrangeRed'
         # tres.cnFillPalette = 'MPL_bwr'
         tres.cnLevelSelectionMode = "ExplicitLevels"
         
         if hasattr(tres,'cnLevels') : del tres.cnLevels
         if var[v] in ['PRECT','PRECC','PRECL'] : tres.cnLevels = np.arange(-5,5+1,1)
         # if var[v] in ['TGCLDLWP'] : tres.cnLevels = np.arange(-28,28+4,4)*1e-3
         # if var[v] in ['TGCLDIWP'] : tres.cnLevels = np.arange(-10,10+2,2)*1e-3
         if not hasattr(tres,'cnLevels') : 
            # if np.min(data_list[c]).values==np.max(data_list[c]).values : 
            if np.min(data_list[c])==np.max(data_list[c]) : 
               print(hc.tcolor.RED+'WARNING: Difference is zero!'+hc.tcolor.ENDC)
            else:
               cmin,cmax,cint,clev = ngl.nice_cntr_levels(diff_data_min, diff_data_max,    \
                                                          cint=None, max_steps=21,      \
                                                          returnLevels=True, aboutZero=True )
               tres.cnLevels = np.linspace(cmin,cmax,num=21)
         
         ### override the level settings and just use auto
         # tres.cnLevelSelectionMode = "AutomaticLevels"

         if use_remap:
            hs.set_cell_fill(tres,case_obj=case_obj,lat=lat_list[c],lon=lon_list[c])
         else:
            hs.set_cell_fill(tres,case_obj=case_obj)

         tres.lbLabelBarOn = True

         ipd = ip
         if add_diff and     var_x_case: ipd = ip+1
         if add_diff and not var_x_case: ipd = ip+num_var*(num_case-1)

         plot[ipd] = ngl.contour_map(wks,data_list[c],tres)

         #-----------------------------------
         ####################################
         # overlay stippling for significance
         ####################################
         #-----------------------------------
         if chk_significance:
            sres = ngl.Resources()
            sres.nglDraw                = False
            sres.nglFrame               = False
            sres.cnFillOn               = True
            sres.cnLinesOn              = False
            sres.cnLineLabelsOn         = False
            sres.cnInfoLabelOn          = False
            sres.lbLabelBarOn           = False
            sres.cnMonoFillPattern      = False
            sres.cnMonoFillColor        = True
            sres.cnFillColor            = 'black'
            sres.cnLevelSelectionMode   = "ExplicitLevels"
            sres.cnLevels               = np.array([0.5])
            sres.cnFillPatterns         = np.array([-1,17])
            # sres.cnFillScaleF, sres.cnFillDotSizeF = 0.6, 0.0015
            sres.cnFillScaleF, sres.cnFillDotSizeF = 0.3, 0.0010
            if use_remap:
               lat2D =               np.repeat( lat_list[c][...,None],len(lon_list[c]),axis=1)
               lon2D = np.transpose( np.repeat( lon_list[c][...,None],len(lat_list[c]),axis=1) )
               sres.sfXArray      = lon2D
               sres.sfYArray      = lat2D
            else:
               sres.sfXArray      = lon
               sres.sfYArray      = lat
            if np.sum(sig_list[c-1]) != 0:
               ngl.overlay( plot[ip], ngl.contour(wks,sig_list[c-1],sres) )
               # plot[ip] = ngl.contour(wks,sig_list[c-1],sres)
         #-----------------------------------
         ####################################
         #-----------------------------------

         ctr_str = ''
         # case_name = name[c]+' - '+base_name
         if 'name' in vars():
            case_name = name[c]
         else:
            case_name = case_obj.short_name
         
         hs.set_subtitles(wks, plot[ipd], case_name, '', var_str+' (Diff)', font_height=subtitle_font_height)
         # hs.set_subtitles(wks, plot[ipd], case_name, '', var_str, font_height=subtitle_font_height)
   

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
# if plot_diff : num_case = num_case+len(diff_case)   # use this to plot both before and after diff

num_case_alt = num_case*2-1 if (plot_diff and add_diff) else num_case
layout = [num_var,num_case_alt] if var_x_case else [num_case_alt,num_var]

if not (plot_diff and add_diff):
   if num_case==1 :
      if num_var==4 : layout = [2,2]
      if num_var==6 : layout = [3,2]
   if num_var==1 :
      if num_case==4: layout = [2,2]

pnl_res = hs.setres_panel()

### add panel labels
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.01
pnl_res.nglPanelFigureStringsFontHeightF = 0.015

pnl_res.nglPanelYWhiteSpacePercent = 8

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
