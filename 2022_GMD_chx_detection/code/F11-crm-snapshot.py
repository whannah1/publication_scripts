# plot a snapshot of CRM states along with the corresponding map of GCM data
import os, ngl, sys, numba, copy, xarray as xr, numpy as np, re, string
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import pg_checkerboard_utilities as pg
import cmocean

case,name,clr,dsh,pat = [],[],[],[],[]
def add_case(case_in,n='',c='black',d=0,p=0):
   global name,case
   case.append(case_in); name.append(n); clr.append(c); dsh.append(d); pat.append(p)
#-------------------------------------------------------------------------------
add_case('E3SM.CHX.ne30pg2_ne30pg2.F-MMF1.00', n='E3SM-MMF')
#-------------------------------------------------------------------------------

gcm_var = 'TGCLDLWP'
crm_var = 'CRM_QV'

ilat,ilon = 10,143 ; dx=30; lat1,lat2 , lon1,lon2 = ilat-dx,ilat+dx , ilon-dx,ilon+dx

# file_num = 365*4+180

fig_type = 'png'
fig_file = 'figs/F11-crm-snapshot'

time_stamp = '0005-06-30'
scratch_path = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl'
data_file_h1 = f'{scratch_path}/{case[0]}/run/{case[0]}.eam.h1.{time_stamp}-00000.nc'
data_file_h2 = f'{scratch_path}/{case[0]}/run/{case[0]}.eam.h2.{time_stamp}-00000.nc'

#---------------------------------------------------------------------------------------------------
# Set up plot
#---------------------------------------------------------------------------------------------------
num_case = len(case)

wkres = ngl.Resources()
npix=4096; wkres.wkWidth,wkres.wkHeight=npix,npix # use this for plotting all patterns w/ rotation
wks = ngl.open_wks(fig_type,fig_file,wkres)
# plot = [None]

#---------------------------------------------------------------------------------------------------
# Set up plot resources for GCM data
#---------------------------------------------------------------------------------------------------
gcm_res = hs.res_contour_fill_map()
gcm_res.vpHeightF = 0.5
gcm_res.vpWidthF  = 0.5
gcm_res.vpXF = 0.1
# gcm_res.vpYF = 0.9
gcm_res.nglMaximize = False
gcm_res.tmXBOn,gcm_res.tmYLOn = False,False
gcm_res.lbTitlePosition       = 'Bottom'
gcm_res.lbTitleFontHeightF    = 0.015
gcm_res.lbLabelFontHeightF    = 0.015
gcm_res.lbTitleString         = 'kg/m2'
# gcm_res.lbLabelBarOn = False
if 'lat1' in vars() : gcm_res.mpMinLatF,gcm_res.mpMaxLatF = lat1,lat2
if 'lon1' in vars() : gcm_res.mpMinLonF,gcm_res.mpMaxLonF = lon1,lon2
scrip_file_path = 'scrip_files/ne30pg2_scrip.nc'
scripfile = xr.open_dataset(scrip_file_path)
gcm_res.cnFillMode    = "CellFill"
gcm_res.sfXArray      = scripfile['grid_center_lon'].values
gcm_res.sfYArray      = scripfile['grid_center_lat'].values
gcm_res.sfXCellBounds = scripfile['grid_corner_lon'].values
gcm_res.sfYCellBounds = scripfile['grid_corner_lat'].values
gcm_res.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )

if gcm_var=="TGCLDLWP": gcm_res.cnLevels = np.arange(0.1,1.,0.1)

if hasattr(gcm_res,'cnLevels') : gcm_res.cnLevelSelectionMode = 'ExplicitLevels'

### resources for box around inset region
pgres = ngl.Resources()
pgres.nglDraw,pgres.nglFrame = False,False
pgres.gsLineColor = 'red'
pgres.gsLineThicknessF = 12

#---------------------------------------------------------------------------------------------------
# Set up plot resources for CRM data
#---------------------------------------------------------------------------------------------------

crm_res = hs.res_contour_fill()
crm_res.vpHeightF = 0.4
crm_res.tmXBOn= False ; # crm_res.tmYLOn = False
crm_res.lbLabelBarOn = False
crm_res.nglYAxisType = "LinearAxis"

crm_res.tmYLMode   = 'Explicit'
crm_res.tmYLValues = np.arange(0,30,2)
crm_res.tmYLLabels = np.arange(0,30,2).astype(int)

crm_res.cnFillPalette = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )

#---------------------------------------------------------------------------------------------------
# Load SCRIP data
#---------------------------------------------------------------------------------------------------
center_lat = scripfile['grid_center_lat'].values
center_lon = scripfile['grid_center_lon'].values
corner_lon = scripfile['grid_corner_lon'].values
corner_lat = scripfile['grid_corner_lat'].values
ncol = len(center_lat)

#---------------------------------------------------------------------------------------------------
# find adjacent neighbors around main point
#---------------------------------------------------------------------------------------------------

if ilat==6 and ilon==139: 
   mcol = 9246
   neighborhood = [9246, 9241, 9243, 9244, 9245, 9247, 9364, 9365, 9361]
   neighborbear = [  -1, -129,  -95,  179,  124,   84,    0,   48,  -56]

if ilat==10 and ilon==143:
   mcol = 9489
   neighborhood =  [9489, 9370, 9371, 9488, 9374, 9492, 9491, 9494, 9490]
   neighborbear =  [  -1, -133,  179,  -97,  125,   82,    0,   45,  -55]

if 'neighborhood' not in locals():
   dist = pg.calc_great_circle_distance( ilat, center_lat[:], ilon, center_lon[:] )
   mcol = np.argsort(dist,kind='mergesort')[0]
   num_neighbors = 8
   neighbors = np.full([ncol,num_neighbors+1],-1)
   bearings  = np.full([ncol,num_neighbors+1],-1)
   edge_flag = np.full([ncol,num_neighbors+1],-1)
   # Find adjacent neighbors
   pg.find_neighbors(center_lat, center_lon, corner_lat, corner_lon, neighbors, bearings, edge_flag)
   neighborhood = neighbors[mcol,:]
   neighborbear = bearings[mcol,:]
   print(f' mcol:         {mcol}')
   print(f' neighborhood: {neighborhood}')
   print(f' neighborbear: {neighborbear}')

neighborbear = np.array(neighborbear)
neighborhood = np.array(neighborhood)

tmp_neighborhood = neighborhood[1:]
neighborhood[1:] = tmp_neighborhood[ np.argsort(neighborbear[1:]) ]

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
cnt_ds_list = []
for c in range(num_case):
   print('\n  case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)

   # case_obj = he.Case( name=case[c], time_freq='daily' )
   # if 'lev' not in vars() : lev = np.array([-1])
   # comp = 'eam'
   # use_remap = False

   #----------------------------------------------------------------------------
   # plot GCM data on map
   #----------------------------------------------------------------------------
   print('    var: '+hc.tcolor.GREEN+f'{gcm_var:10}'+hc.tcolor.ENDC)

   # gcm_data = case_obj.load_data(gcm_var,htype='h1',num_files=1,first_file=file_num).isel(time=0)
   
   ds = xr.open_dataset(data_file_h1)
   gcm_data = ds[gcm_var].isel(time=0)

   gcm_plot = ngl.contour_map(wks, gcm_data, gcm_res) 

   ### draw a box that outlines the points used for the CRM snapshot
   neighb_ind = [1,3,5,7,1]
   corner_ind = [0,3,2,1,0]
   clat = corner_lat[neighborhood[neighb_ind],corner_ind]
   clon = corner_lon[neighborhood[neighb_ind],corner_ind]
   pdum = ngl.add_polyline(wks,gcm_plot,clon,clat,pgres)

   var_name = gcm_var
   if var_name=='TGCLDLWP': var_name = 'Liq Water Path'
   if var_name=='PRECT':    var_name = 'Precipitation'
   hs.set_subtitles(wks, gcm_plot, var_name, '', '', font_height=0.015)

   ngl.draw(gcm_plot)
   # ngl.frame(wks)

   #----------------------------------------------------------------------------
   # Add some boxes and lines to clean up the presentation
   #----------------------------------------------------------------------------
   poly_res = ngl.Resources()
   poly_res.gsFillColor = 'white'
   poly_res.gsLineColor = 'red'

   # x1,x2 = 0.5,1.0

   x1,x2 = 0.45,0.95
   y1,y2 = 0.32,0.68

   # adjustment after adding colorbar
   y1 = y1 - 0.01
   y2 = y2 + 0.01


   bx = [x1,x1,x2,x2,x1]
   by = [y1,y2,y2,y1,y1]


   poly_res.gsLineThicknessF = 6
   poly_res.gsLineDashPattern = 1
   for i in range(len(clat)):
   # for i in [0,1]:
      ndcx, ndcy = ngl.datatondc(gcm_plot, clon[i], clat[i] )
      lx = [ ndcx, bx[i] ]
      ly = [ ndcy, by[i] ]
      ngl.polyline_ndc(wks,lx,ly,poly_res)


   poly_res.gsLineThicknessF = 20
   poly_res.gsLineDashPattern = 0
   ngl.polygon_ndc(wks,bx,by,poly_res)
   ngl.polyline_ndc(wks,bx,by,poly_res)

   

   #----------------------------------------------------------------------------
   # Make grid of CRM snapshots
   #----------------------------------------------------------------------------
   print('    var: '+hc.tcolor.GREEN+f'{crm_var:10}'+hc.tcolor.ENDC)
   
   # crm_data = case_obj.load_data(crm_var,htype='h2',num_files=1,first_file=file_num).isel(ncol=neighborhood,time=0,crm_ny=0)
   # hgt_data = case_obj.load_data('Z3',   htype='h2',num_files=1,first_file=file_num).isel(ncol=neighborhood,time=0)

   ds = xr.open_dataset(data_file_h2)
   crm_data = ds[crm_var].isel(ncol=neighborhood,time=0,crm_ny=0)
   hgt_data = ds['Z3'   ].isel(ncol=neighborhood,time=0)

   # adjust units
   if crm_var=='CRM_QV': crm_data = crm_data*1e3

   # Remove horz mean
   crm_data = crm_data - crm_data.mean(dim='crm_nx')

   data_min = np.min(crm_data.values)
   data_max = np.max(crm_data.values)

   nclev = 21 ; aboutZero = False
   (cmin,cmax,cint) = ngl.nice_cntr_levels(data_min, data_max,  cint=None, max_steps=nclev, returnLevels=False, aboutZero=aboutZero )
   crm_res.cnLevels = np.linspace(cmin,cmax,num=nclev)
   crm_res.cnLevelSelectionMode = 'ExplicitLevels'

   plot_pos = [4,6,3,0,1,2,5,8,7] # use this to arrange plots according to their geographic position

   crm_plots = [None]*9
   for i in range(len(crm_plots)):
      ip = plot_pos[i]
      # tmp_hgt = hgt_data[-48:,i].values
      nlev_gcm = len(hgt_data['lev'])
      nlev = 35 # omit top levels in stratosphere for cleaner plot
      tmp_crm = crm_data[:nlev,:,i]
      tmp_hgt = hgt_data[nlev_gcm-nlev:nlev_gcm,i].sortby('lev', ascending=False) /1e3
      crm_res.sfYArray = tmp_hgt.values
      crm_plots[ip] = ngl.contour( wks, tmp_crm.values, crm_res ) 
      tlat,tlon = center_lat[neighborhood[i]],center_lon[neighborhood[i]]
      hs.set_subtitles(wks, crm_plots[ip],'',f'{tlat:5.1f}N {tlon:5.1f}E','', font_height=0.01)
   
   pnl_res = hs.setres_panel()

   pnl_res.nglPanelLabelBar = True
   pnl_res.lbTitleFontHeightF = 0.01
   pnl_res.nglPanelLabelBarLabelFontHeightF = 0.01
   pnl_res.lbTitlePosition = 'Bottom'
   pnl_res.lbTitleString = 'g/kg'
   pnl_res.lbTopMarginF       =  0.2
   pnl_res.lbBottomMarginF    = -0.2

   pnl_res.nglPanelLeft  = 0.45
   pnl_res.nglPanelRight = 0.95
   ngl.panel(wks,crm_plots,[3,3],pnl_res)

   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------

   
      

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

# pnl_res = hs.setres_panel()
# pnl_res.nglPanelYWhiteSpacePercent = 5
# pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
# pnl_res.nglPanelFigureStringsJust        = "TopLeft"
# pnl_res.nglPanelFigureStringsFontHeightF = 0.015   


# pnl_res.nglPanelLabelBar = True
# pnl_res.lbTitleFontHeightF = 0.015
# # pnl_res.lbLabelFontHeightF = 0.01
# pnl_res.nglPanelLabelBarLabelFontHeightF = 0.015
# pnl_res.lbTitlePosition = 'Bottom'
# pnl_res.lbTitleString = 'Fractional Occurrence'

# ngl.panel(wks,[plot],[1,1],pnl_res)

# ngl.panel(wks,crm_plots,[3,3],pnl_res)

ngl.end()

hc.trim_png(fig_file)
# print(f'\n{fig_file}.png\n')
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------