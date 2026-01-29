# v2 is a cleaned up version of v1 - functionality is mostly the same
import os, ngl, subprocess as sp, numpy as np, xarray as xr, dask, copy, string, cmocean, glob
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
host = hc.get_host()
#-------------------------------------------------------------------------------
name,case,case_root,case_sub = [],[],[],[]
scrip_file_list = []
def add_case(case_in,n=None,p=None,s=None,g=None,c=None,scrip_file=None):
   global name,case,case_root,case_sub
   tmp_name = case_in if n is None else n
   case.append(case_in); name.append(tmp_name); 
   case_root.append(p); case_sub.append(s);
   scrip_file_list.append(scrip_file)
#-------------------------------------------------------------------------------
var,lev_list,var_str = [],[],[]
def add_var(var_name,lev=0,s=None): 
   var.append(var_name); lev_list.append(lev); 
   if s is None:
      var_str.append(var_name)
   else:
      var_str.append(s)
#-------------------------------------------------------------------------------
fig_file,fig_type = 'figs/FXX-map','png'
#-------------------------------------------------------------------------------

### 2024 AQP/aqua CESS runs with grid sensitivity
tmp_scratch,tmp_sub = '/pscratch/sd/w/whannah/2024-AQP-CESS','archive/atm/hist'
scrip_file_root = '/global/cfs/projectdirs/m3312/whannah/HICCUP/files_grid/'

add_case('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne30pg2_ne30pg2.NN_32.SSTP_0K', n='AQP EAM ne30 +0K',p=tmp_scratch,s=tmp_sub,scrip_file=f'{scrip_file_root}/scrip_ne30pg2.nc')
# add_case('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne30pg2_ne30pg2.NN_32.SSTP_4K',   n='AQP EAM ne30 +4K',p=tmp_scratch,s=tmp_sub,scrip_file=f'{scrip_file_root}/scrip_ne30pg2.nc')
# add_case('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne45pg2_ne45pg2.NN_64.SSTP_0K',   n='AQP EAM ne45 +0K',p=tmp_scratch,s=tmp_sub,scrip_file=f'{scrip_file_root}/scrip_ne45pg2.nc')
# add_case('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne60pg2_ne60pg2.NN_128.SSTP_0K',  n='AQP EAM ne60 +0K',p=tmp_scratch,s=tmp_sub,scrip_file=f'{scrip_file_root}/scrip_ne60pg2.nc')
# add_case('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne90pg2_ne90pg2.NN_256.SSTP_0K',  n='AQP EAM ne90 +0K',p=tmp_scratch,s=tmp_sub,scrip_file=f'{scrip_file_root}/scrip_ne90pg2.nc')
add_case('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne120pg2_ne120pg2.NN_512.SSTP_0K',n='AQP EAM ne120 +0K',p=tmp_scratch,s=tmp_sub,scrip_file=f'{scrip_file_root}/scrip_ne90pg2.nc')

htype,first_file,last_file = 'h0',int(300/10),int(310/10)
# htype,first_file,last_file = 'h0',int(300/10),int(500/10)
# htype,first_file,last_file = 'h0',int(300/10),int(320/10)

#-------------------------------------------------------------------------------
# add_var('FLNTC')
# add_var('FSNTC')

add_var('PS')

# add_var('LWCF')
# add_var('SWCF')

# add_var('PRECC',s='Precipitation')
# add_var('TGCLDLWP',s='Liq Water Path')
# add_var('TGCLDIWP',s='Ice Water Path')
# add_var('Z3',lev=-63)
# add_var('TS')
# add_var('PSL')
# add_var('PRECSC')
# add_var('PRECC'); add_var('PRECL')
# add_var('P-E')
# add_var('TMQ')
# add_var('LHFLX')
# add_var('SHFLX')
# add_var('TS')
# add_var('NET_TOA_RAD')
# add_var('RESTOM')
# add_var('FLUTOA')
# add_var('FSNS'); add_var('FLNS')

# add_var('U10')
# add_var('TBOT')
# add_var('QBOT')
# add_var('TAUX')
# add_var('WSPD_BOT')

### use for h0
# add_var('U',lev=975)
# add_var('V',lev=975)
# add_var('U',lev=850,s='U850')
# add_var('U',lev=500)
# add_var('U',lev=200,s='U200')

# add_var('OMEGA',lev=500)


#-------------------------------------------------------------------------------


use_remap,remap_grid = False,'90x180'

plot_diff,add_diff = False,False

use_snapshot,ss_t    = False,-1
print_stats          = True
var_x_case           = True

num_plot_col         = len(case)

use_common_label_bar = False

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
if case==[]: raise ValueError('ERROR - case list is empty!')
num_var,num_case = len(var),len(case)

subtitle_font_height = 0.01

diff_base = 0

wkres = ngl.Resources()
npix = 2048*2; wkres.wkWidth,wkres.wkHeight=npix,npix

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

res.mpOutlineBoundarySets = 'NoBoundaries'
res.mpCenterLonF = 0.

res.mpProjection = 'Mollweide'
# res.mpProjection = "Orthographic"

#---------------------------------------------------------------------------------------------------
def load_data(var,case,case_root,case_sub):
   global htype,first_file,num_files
   file_list = sorted(glob.glob(f'{case_root}/{case}/{case_sub}/*eam.{htype}*'))
   # file_list = file_list[first_file:first_file+num_files]
   file_list = file_list[first_file:last_file+1]
   ds = xr.open_mfdataset(file_list)
   ds = ds.mean(dim='time')
   da = ds[var]
   if var=='PRECC': da = da*86400*1e3
   if var=='PRECT': da = da*86400*1e3
   return da
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
      #-------------------------------------------------------------------------
      # data_dir_tmp,data_sub_tmp = None, None
      # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
      # if case_root[c] is not None: data_dir_tmp = case_root[c]
      # if case_sub[c] is not None: data_sub_tmp = case_sub[c]
      # case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp,
      #                     populate_files=True )
      # case_obj.set_coord_names(var[v])
      #-------------------------------------------------------------------------
      # read the data
      #-------------------------------------------------------------------------   
      # with dask.config.set(**{'array.slicing.split_large_chunks': True}):
      #    data = case_obj.load_data(var[v],htype=htype,ps_htype=htype,lev=lev,
      #                                    first_file=first_file,num_files=num_files,
      #                                    use_remap=use_remap,remap_str=f'remap_{remap_grid}')
      #-------------------------------------------------------------------------
      data = load_data(var[v],case[c],case_root[c],case_sub[c])
      #-------------------------------------------------------------------------
      # Special handling of various specific circumstances
      #-------------------------------------------------------------------------
      if case[c]=='ERAi': erai_lat,erai_lon = data['lat'],data['lon']
      if 'crm_nx' in data.dims : data = data.mean(dim=('crm_nx','crm_ny')).isel(crm_nz=15)
      if 'lev' in data.dims : data = data.isel(lev=0)
      #-------------------------------------------------------------------------
      # print stats before time averaging
      #-------------------------------------------------------------------------
      if print_stats: hc.print_stat(data,name=var[v],stat='naxsh',indent='    ',compact=True)
      #-------------------------------------------------------------------------
      # average over time dimension
      #-------------------------------------------------------------------------
      if 'time' in data.dims : 
         hc.print_time_length(data.time,indent=' '*6,print_span=True, print_length=False)
         if use_snapshot:
            data = data.isel(time=ss_t); print(hc.tcolor.RED+'WARNING - snapshot mode enabled'+hc.tcolor.ENDC)
         else:
            data = data.mean(dim='time')
      #-------------------------------------------------------------------------
      # Calculate area weighted global mean
      #-------------------------------------------------------------------------
      if 'area' in locals() :
         gbl_mean = ( (data*area).sum() / area.sum() ).values 
         print(hc.tcolor.CYAN+f'      Area Weighted Global Mean : {gbl_mean:6.4}'+hc.tcolor.ENDC)
         # glb_avg_list.append(gbl_mean)
      #-------------------------------------------------------------------------
      # append to data lists
      #-------------------------------------------------------------------------
      if case[c]=='TRMM' and 'lon1' not in locals(): data = ngl.add_cyclic(data.values)
      
      data_list.append( data.values )
      #-------------------------------------------------------------------------
      # save baseline for diff map
      #-------------------------------------------------------------------------
      if plot_diff :
         if c==diff_base:
            data_baseline = data.copy()
   #----------------------------------------------------------------------------
   # calculate common limits for consistent contour levels
   #----------------------------------------------------------------------------
   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])
   if plot_diff:
      tmp_data = data_list - data_list[diff_base]
      for c in range(num_case): tmp_data[c] = data_list[c] - data_list[diff_base]
      diff_data_min = np.min([np.nanmin(d) for d in tmp_data])
      diff_data_max = np.max([np.nanmax(d) for d in tmp_data])
   #----------------------------------------------------------------------------
   # Plot averaged data
   #----------------------------------------------------------------------------
   for c in range(num_case):
      #-------------------------------------------------------------------------
      # Set color palette
      #-------------------------------------------------------------------------
      tres = copy.deepcopy(res)
      tres.cnFillPalette = "MPL_viridis"
      # tres.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
      # tres.cnFillPalette = np.array( cmocean.cm.amp(np.linspace(0,1,256)) )
      if var[v] in ['P-E']                      : tres.cnFillPalette = "BlueWhiteOrangeRed"
      if var[v] in ['CLDLOW','CLDMED','CLDHGH'] : tres.cnFillPalette = "CBR_wet"
      if var[v] in ['TGCLDLWP','TGCLDIWP']      : tres.cnFillPalette = "MPL_viridis"
      if var[v] in ['DYN_QLIQ']                 : tres.cnFillPalette = "MPL_viridis"
      # if var[v] in ['TS','PS']                  : tres.cnFillPalette = 'BlueWhiteOrangeRed'
      # if var[v] in ['TS','PS']                  : tres.cnFillPalette = 'WhiteBlueGreenYellowRed'
      # if var[v] in ['TS','PS'] : tres.cnFillPalette = np.array( cmocean.cm.amp(np.linspace(0,1,256)) )
      if var[v] in ['U','V','UBOT','VBOT','U850','V850','U200','V200','MMF_DU','MMF_DV']: 
         # tres.cnFillPalette = "BlueWhiteOrangeRed"
         tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )

      #-------------------------------------------------------------------------
      # Set contour levels
      #-------------------------------------------------------------------------
      # if var[v] in ['PRECT','PRECC','PRECL']   : tres.cnLevels = np.arange(4,80+4,4)
      if var[v] in ['PRECT','PRECC','PRECL']   : tres.cnLevels = np.arange(2,20+2,2)
      # if var[v] in ['PRECT','PRECC','PRECL']   : tres.cnLevels = np.arange(5,100+5,5) # for std dev
      # if var[v] in ['PRECT','PRECC','PRECL']   : tres.cnLevels = np.logspace( -3, 2.5, num=60).round(decimals=2)
      # if var[v]=='LHFLX'               : tres.cnLevels = np.arange(5,205+5,5)
      if var[v]=='P-E'                 : tres.cnLevels = np.linspace(-10,10,21)
      if var[v]=='RH'                  : tres.cnLevels = np.arange(10,100+1,1)
   
      if var[v]=='PS'                  : tres.cnLevels = np.arange(940,1040+1,1)*1e2
      if var[v] in ['TGCLDIWP','TGPRCIWP']: tres.cnLevels = np.arange(1,30+1,1)*1e-2
      if var[v] in ['TGCLDLWP','TGPRCLWP']: tres.cnLevels = np.logspace( -2, 0.25, num=60).round(decimals=2)
      #-------------------------------------------------------------------------
      # set non-explicit contour levels
      #-------------------------------------------------------------------------
      if hasattr(tres,'cnLevels') : 
         tres.cnLevelSelectionMode = 'ExplicitLevels'
      else:
         nlev = 41
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
      #-------------------------------------------------------------------------

      ds_grid = xr.open_dataset(scrip_file_list[c])
      tres.cnFillMode    = 'CellFill'
      tres.sfXArray      = ds_grid['grid_center_lon'].values
      tres.sfYArray      = ds_grid['grid_center_lat'].values
      tres.sfXCellBounds = ds_grid['grid_corner_lon'].values
      tres.sfYCellBounds = ds_grid['grid_corner_lat'].values

      tres.lbLabelBarOn = False if use_common_label_bar else True      
         
      # if plot_diff and c==diff_base : base_name = name[c]

      num_case_alt = num_case*2-1 if (plot_diff and add_diff) else num_case
      ip = v*num_case_alt+c if var_x_case else c*num_var+v

      if not plot_diff  or (plot_diff and add_diff) or (plot_diff and c==diff_base) : 

         plot[ip] = ngl.contour_map(wks,np.ma.masked_invalid(data_list[c]),tres) 
         
         # set plot subtitles
         ctr_str = ''
         if glb_avg_list != []: ctr_str = f'{glb_avg_list[c]:6.4}'
         hs.set_subtitles(wks, plot[ip], name[c], ctr_str, var_str[v], font_height=subtitle_font_height)

      #-------------------------------------------------------------------------
      # create difference plot
      #-------------------------------------------------------------------------
      if plot_diff and c!=diff_base :
         
         data_list[c] = data_list[c] - data_baseline.values
         
         tres.cnFillPalette = 'BlueWhiteOrangeRed'
         # tres.cnFillPalette = np.array( cmocean.cm.diff(np.linspace(0,1,256)) )
         # tres.cnFillPalette = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )
         # tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )

         if hasattr(tres,'cnLevels') : del tres.cnLevels
         if var[v] in ['PRECT','PRECC','PRECL'] : tres.cnLevels = np.arange(-5,5+1,1)
         if var[v] in ['MMF_DU','MMF_DV']       : tres.cnLevels = np.linspace(-2,2,11)
         if var[v]=='TS'                        : tres.cnLevels = np.arange(-40,40+5,5)/1e1
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
         
         ctr_str = 'Diff'
         if glb_avg_list != []: 
            glb_diff = glb_avg_list[c] - glb_avg_list[diff_base]
            ctr_str += f' ({glb_diff:6.4})'
         
         hs.set_subtitles(wks, plot[ipd], name[c], ctr_str, var_str[v], font_height=subtitle_font_height)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

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
