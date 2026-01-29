# plot the fractional occurence of all possible patterns
import os, ngl, sys, numba, copy, xarray as xr, numpy as np, re, string
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import pg_checkerboard_utilities as pg
import cmocean

case,name,clr,dsh,pat = [],[],[],[],[]
var,var_name,lev_list = [],[],[]
def add_case(case_in,n='',c='black',d=0,p=0):
   global name,case
   case.append(case_in); name.append(n); clr.append(c); dsh.append(d); pat.append(p)

def add_var(v,n=None,lev=-1): var.append(v); var_name.append(n); lev_list.append(lev)
#-------------------------------------------------------------------------------
# add_case('MAC-FV',  n='MAC 1x1 deg'                                          ,c='black',p=0)
# add_case('MAC-PG',  n='MAC ne30pg2'                                          ,c='gray',p=0)
# add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.no_dcape.00',n='E3SM (no DCAPE)',c='green',p=0)
add_case('E3SM.CHX.ne30pg2_ne30pg2.F-MMF1.00',            n='E3SM-MMF'       ,c='blue' ,p=0)
add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.00',         n='E3SM'           ,c='red'  ,p=0)
#-------------------------------------------------------------------------------

var = []
# add_var('PRECT')
add_var('TGCLDLWP',n='Liq Water Path')
# add_var('TGCLDIWP',n='Ice Water Path')
# add_var('TMQ')
# add_var('LHFLX')
# add_var('FLNT')
# add_var('FSNT')
# add_var('Q850')
# add_var('T850')
# add_var('U850')


# use_regional_subset = False
lat1,lat2 = -50,50
lon1,lon2 = 110,360


fig_type = 'png'
fig_file = 'figs/F09-chx-phase-map'
tmp_file = 'data/occurrence-chx-only'


common_colorbar   = True   # use common set of color levels for all plots
ocean_only        = False
add_inset         = True # only works with single case

var_x_case           = True    # controls plot panel arrangement
print_stats          = True    #
verbose              = True    #

#---------------------------------------------------------------------------------------------------
# generate sets of possible neighbor states
#---------------------------------------------------------------------------------------------------

### only pure checkboard
tmp_file='data/occurrence-chx-only'; rotate_sets=False; sets=pg.chx_only_sets

### find index of pure checkerboard
# chx_idx = []
# for s,tset in enumerate(sets):
#    for chx in chx_sets:
#       if all(tset==chx): chx_idx.append(s)


(num_set,set_len) = sets.shape
set_coord,nn_coord = np.arange(num_set),np.arange(set_len)
sets.assign_coords(coords={'set':set_coord,'neighbors':nn_coord})

set_labels = pg.get_set_labels(sets)

for s,label in enumerate(set_labels): set_labels[s] = label.replace(' ','')
# for label in set_labels: print(label)
# exit()

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_comp(case):
   comp = 'eam'
   if case=='MAC': comp = None
   if 'CESM' in case: comp = 'cam'
   if case=='E3SM.RGMA.ne30pg2_r05_oECv3.FC5AV1C-L.00': comp = 'cam'
   return comp

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_var,num_set = len(case),len(var),len(sets)

wkres = ngl.Resources()
npix=4096; wkres.wkWidth,wkres.wkHeight=npix,npix 
wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*(num_var)

res = hs.res_contour_fill_map()
res.tmXBOn = False
res.tmYLOn = False
# res.lbTitlePosition = 'Bottom'
# res.lbTitleFontHeightF = 0.01
# res.lbTitleString = 'Count'
if common_colorbar:
   res.lbLabelBarOn = False
else:
   res.lbLabelBarOn = True

if 'lat1' in vars() : res.mpMinLatF = lat1
if 'lat2' in vars() : res.mpMaxLatF = lat2
if 'lon1' in vars() : res.mpMinLonF = lon1
if 'lon2' in vars() : res.mpMaxLonF = lon2

### resources for box around inset region
pgres = ngl.Resources()
pgres.nglDraw,pgres.nglFrame = False,False
pgres.gsLineColor = 'red'
pgres.gsLineThicknessF = 6

### resources for visual aid lines between inset and main map
plres = ngl.Resources()
# plres.nglDraw,pgres.nglFrame = False,False
plres.gsLineColor = 'black'
plres.gsLineThicknessF = 2

### resources for attaching inset plot as annotation
ares                = ngl.Resources()
ares.amZone         = 1

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
cnt_ds_list = []
for c in range(num_case):
   print('\n  case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)

   tname = case[c]
   if 'MAC'  in case[c]: tname = 'MAC'
   case_obj = he.Case( name=tname, time_freq='daily' )
   if 'lev' not in vars() : lev = np.array([-1])

   comp = 'eam'
   if case[c]=='EAR5': comp = None
   if case[c]=='E3SM.RGMA.ne30pg2_r05_oECv3.FC5AV1C-L.00': comp = 'cam'

   use_remap = False
   remap_str=f'remap_ne30pg2'
   if case[c]=='MAC-FV' : use_remap = False
   if case[c]=='MAC-PG' : use_remap = True

   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   for v in range(num_var) :
      case_tmp_file = f'{tmp_file}.daily.{case[c]}.{var[v]}.nc'
      
      print('    var: '+hc.tcolor.GREEN+f'{var[v]:10}'+hc.tcolor.ENDC+'    '+case_tmp_file)
      
      cnt_ds = xr.open_dataset( case_tmp_file )

      #-------------------------------------------------------------------------
      # print the time length
      #-------------------------------------------------------------------------
      msg = '    num_time: '
      msg += str(cnt_ds['num_time'].values)
      msg += ' days'
      print(msg)
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      
      cnt_ds['cnt'] = cnt_ds['cnt'] / cnt_ds['num_valid'] 

      # print(); print(cnt_ds); print()
      ### replace NaNs with -1
      # cnt_ds['cnt'] = np.where( cnt_ds['cnt'].values==np.nan, np.zeros(cnt_ds['cnt'].shape), cnt_ds['cnt'].values )
      # cnt_ds = cnt_ds.where( cnt_ds==np.nan, 0. )

      if print_stats: hc.print_stat(cnt_ds['cnt'],name='final count dataset',stat='naxs',indent='    ')

      cnt_ds_list.append(cnt_ds)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

# res.cnFillPalette = "MPL_viridis"
res.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
# res.cnFillPalette = np.array( cmocean.cm.amp(np.linspace(0,1,256)) )

if common_colorbar:
   data_min = np.zeros(num_set)
   data_max = np.zeros(num_set)
   for cnt_ds in cnt_ds_list:
      for s in range(num_set):
         data_min[s] = np.min([ data_min[s], cnt_ds['cnt'].isel(set=s).min().values ])
         data_max[s] = np.max([ data_max[s], cnt_ds['cnt'].isel(set=s).max().values ])

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
list_cnt = 0
plot = [None]*(num_var*num_case*num_set)
anno = [None]*(num_var*num_case*num_set)
pdum = [None]*(num_var*num_case*num_set)
ldum1 = [None]*(num_var*num_case*num_set)
ldum2 = [None]*(num_var*num_case*num_set)
ldum3 = [None]*(num_var*num_case*num_set)
ldum4 = [None]*(num_var*num_case*num_set)

for c in range(num_case):
   
   tname = case[c]
   if 'MAC'  in case[c]: tname = 'MAC'
   case_obj = he.Case( name=tname, time_freq='daily' )

   scrip_file_path = 'scrip_files/ne30pg2_scrip.nc'
   if case[c]=='MAC-FV' : scrip_file_path = 'scrip_files/cmip6_180x360_scrip.20210722.nc'

   scripfile = xr.open_dataset(scrip_file_path)

   if ocean_only: 
      if case[c]=='MAC-FV' : 
         print('ocean_only not implemented for MAC-FV!')
      else:
         landfrac_ds = xr.open_dataset('data/E3SM_landfrac_ne30pg2.nc').isel(time=0)
         mask = xr.DataArray( np.ones(len(cnt_ds['ncol']),dtype=bool), coords=[cnt_ds['ncol']] )
         mask = mask & (landfrac_ds['LANDFRAC'].values<0.5)
         scripfile = scripfile.where(mask.rename({'ncol':'grid_size'}),drop=True)

   res.cnFillMode    = "CellFill"
   res.sfXArray      = scripfile['grid_center_lon'].values
   res.sfYArray      = scripfile['grid_center_lat'].values
   res.sfXCellBounds = scripfile['grid_corner_lon'].values
   res.sfYCellBounds = scripfile['grid_corner_lat'].values

   for v in range(num_var) :

      tres = copy.deepcopy(res)
   
      cnt_ds = cnt_ds_list[list_cnt]
      list_cnt += 1

      for s in range(num_set):

         # print(); print(data_min[s])
         # print(); print(data_max[s])
         # print(); exit()

         if common_colorbar:
            if 'nlev' not in locals(): nlev = 21
            aboutZero = False
            # (cmin,cmax,cint) = ngl.nice_cntr_levels(data_min[s], data_max[s], 
            #                                        cint=None, max_steps=nlev, 
            #                                        returnLevels=False, aboutZero=aboutZero )
            # tres.cnLevels = np.linspace(cmin,cmax,num=nlev)
            tres.cnLevels = np.arange(0.004,0.06+0.004,0.004)
            tres.cnLevelSelectionMode = 'ExplicitLevels'

         if num_var>1:
            ip = s*num_var+v if var_x_case else v*num_set+s
         else:
            ip = s*num_case+c if var_x_case else c*num_set+s

         tmp_data = cnt_ds['cnt'].isel(set=s)

         if ocean_only: tmp_data = tmp_data.where(mask,drop=True)

         ### replace NaNs with zero
         # tmp_data = np.where( tmp_data.values==np.nan, 0., tmp_data.values )
         tmp_data = np.where( np.isfinite(tmp_data.values), tmp_data.values, 0. )
         

         # plot[ip] = ngl.contour_map(wks, tmp_data.values, tres) 
         plot[ip] = ngl.contour_map(wks, tmp_data, tres) 

         set_str = set_labels[s]
         subtitle_font_height = 0.015
         hs.set_subtitles(wks, plot[ip], name[c], set_str, var_name[v], font_height=subtitle_font_height)

         #----------------------------------------------------------------------
         #----------------------------------------------------------------------
         if add_inset:
            ires = copy.deepcopy(tres)
            if num_case==1:
               ires.vpHeightF        = 0.1
               ires.vpWidthF         = 0.1
            if num_case==2:
               ires.vpHeightF        = 0.05
               ires.vpWidthF         = 0.05
            ires.nglMaximize         = False

            if case[c]=='E3SM.CHX.ne30pg2_ne30pg2.F-MMF1.00'   :ilat,ilon,idx = -20,360-115,3
            if case[c]=='E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.00':ilat,ilon,idx = 15,360-40,3
            ires.mpMinLatF,ires.mpMaxLatF = ilat-idx,ilat+idx
            ires.mpMinLonF,ires.mpMaxLonF = ilon-idx,ilon+idx
            
            iplot = ngl.contour_map(wks, tmp_data, ires) 

            ### get NDC coordinates of where want to put the inset plot
            alat = -50
            alon = 360-65
            ndcx, ndcy = ngl.datatondc(plot[ip], alon, alat )

            ### attach plot via annotation
            ares.amOrthogonalPosF,ares.amParallelPosF  = ndcy, ndcx
            anno[ip] = ngl.add_annotation(plot[ip], iplot ,ares)

            ### draw box around map region used for inset
            xbar = np.array([ilon-idx,ilon+idx,ilon+idx,ilon-idx,ilon-idx])
            ybar = np.array([ilat-idx,ilat-idx,ilat+idx,ilat+idx,ilat-idx])
            pdum[ip] = ngl.add_polyline(wks,plot[ip],xbar,ybar,pgres)

            ### draw lines between map region and annotation
            # if case[c]=='E3SM.CHX.ne30pg2_ne30pg2.F-MMF1.00':
            #    lx = np.array([ilon-idx,360-90])
            #    ly = np.array([ilat-idx,-39])
            #    ldum1[ip] = ngl.add_polyline(wks,plot[ip],lx,ly,plres)

            #    lx = np.array([ilon-idx,360-90])
            #    ly = np.array([ilat+idx,0])
            #    ldum4[ip] = ngl.add_polyline(wks,plot[ip],lx,ly,plres)

            # if case[c]=='E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.00':
            #    lx = np.array([ilon-idx,360-90])
            #    ly = np.array([ilat-idx,-39])
            #    ldum1[ip] = ngl.add_polyline(wks,plot[ip],lx,ly,plres)

            #    lx = np.array([ilon-idx,360-90])
            #    ly = np.array([ilat+idx,0])
            #    ldum4[ip] = ngl.add_polyline(wks,plot[ip],lx,ly,plres)


         #----------------------------------------------------------------------
         #----------------------------------------------------------------------

         


#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
# hs.set_plot_labels(wks, plot, font_height=0.01, justify='left')

if num_var>1:
   layout = [num_set,num_var] if var_x_case else [num_var,num_set]
else:
   layout = [num_set,num_case] if var_x_case else [num_case,num_set]

# if num_case==1 and num_set>3: 
#    num_plot_col = 6
#    layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent = 5
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.015   

if common_colorbar: 
   pnl_res.nglPanelLabelBar = True
   pnl_res.lbTitleFontHeightF = 0.015
   # pnl_res.lbLabelFontHeightF = 0.01
   pnl_res.nglPanelLabelBarLabelFontHeightF = 0.015
   pnl_res.lbTitlePosition = 'Bottom'
   pnl_res.lbTitleString = 'Fractional Occurrence'

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
