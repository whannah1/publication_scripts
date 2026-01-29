# This script provides a visual check of the method for finding neighbors
# in all methods the great circle bearing is used to sort the neighbors
# v1 - use nearest neighbors based on distance
# v2 - use the SCRIP file corner information to determine edge and corner neighbors
import os, ngl, sys, numba, copy, xarray as xr, numpy as np, string
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import pg_checkerboard_utilities as pg
case,name = [],[]
def add_case(case_in,n=''):
   global name,case
   case.append(case_in); name.append(n)
#-------------------------------------------------------------------------------

add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.00',        n='E3SM')


fig_type = 'png'
fig_file = 'figs/F03-neighbor-examples'

htype,first_file,num_files = 'h1',0,1


# lat1,lat2,lon1,lon2 = 30,40,140,150
num_pts = 16
rng_seed = 1239852 # 2500

# pt_list = [5197,0,14531,21152] # good
pt_list = [3060,21152] # better

# pt_list = [16495] # polar point


print_neighbor_coords = False

if 'pt_list' in locals(): num_pts = len(pt_list)

num_plot_cols = 3

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case = len(case)

wkres = ngl.Resources()
npix=1024
wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*(num_case*num_pts)

res = hs.res_default()
res.mpGridAndLimbOn = False
res.mpLimitMode = 'LatLon'
res.tmYLLabelFontHeightF = 0.02
res.tmXBLabelFontHeightF = 0.02

### for drawing grid cells
pgres = ngl.Resources()
pgres.nglDraw,pgres.nglFrame = False,False
pgres.gsLineColor = 'red'
pgres.gsLineThicknessF = 5.0

### for labelling nearest neighbors
txres               = ngl.Resources()
# txres.txFont        = "helvetica-bold"
# txres.txFontColor   = "coral4"
txres.txFontHeightF = 0.035
# txres.txJust        = "BottomCenter"

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for c in range(num_case):
   print('\n  case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)
   case_obj = he.Case( name=case[c], time_freq='daily' )

   # specify plot width in degrees 
   if 'ne4'  in case[c]: dx = 30
   if 'ne30' in case[c]: dx = 7
   if 'ne45' in case[c]: dx = 5

   #----------------------------------------------------------------------------
   # load the coordinate data
   #----------------------------------------------------------------------------
   comp = 'eam'

   use_remap = False
   remap_str=f'remap_ne30pg2'
   #----------------------------------------------------------------------------
   # get lat and lon coordinates
   #----------------------------------------------------------------------------

   lat = case_obj.load_data('lat',htype=htype,num_files=1,component=comp,use_remap=use_remap,remap_str=remap_str)
   lon = case_obj.load_data('lon',htype=htype,num_files=1,component=comp,use_remap=use_remap,remap_str=remap_str)

   #*************************************************
   ### use this for finding columns in a certain area
   # print('WARNING! - looping through coordinate data (are you sure you want to do that?)')
   # for n in range(len(lat)):
   #    # if lon[n]>-1 and lon[n]<1 and lat[n]>30:
   #    if lat[n] < -80 and lat[n] > -85 : 
   #       print(f'  n: {n}   lat: {lat[n].values:6.2f}   lon: {lon[n].values:6.2f}')
   #*************************************************

   if case[c] in ['MAC'] :
      nlat,nlon = len(lat), len(lon)
      ncol = len(lat) * len(lon)
      ilat,ilon = lat,lon
      lat,lon = np.zeros(ncol), np.zeros(ncol)
      for j in range(nlat):
         for i in range(nlon):
            n = j*nlon + i
            lat[n],lon[n] = ilat[j],ilon[i]
      lat = xr.DataArray(lat,coords=[np.arange(ncol)],dims='ncol')
      lon = xr.DataArray(lon,coords=[np.arange(ncol)],dims='ncol')

   ncol = len(lat)
   #----------------------------------------------------------------------------
   # find nearest neighbors
   #----------------------------------------------------------------------------

   if case[c] in ['MAC'] :
      scripfile = xr.open_dataset(os.getenv('HOME')+'/E3SM/data_grid/fv0.9x1.25_070727.nc')
   else:
      scripfile = case_obj.get_scrip()
   
   # print(case_obj)
   # print()
   # print(scripfile)
   # exit()


   # center_lon = scripfile['grid_center_lon'][:].values
   # center_lat = scripfile['grid_center_lat'][:].values
   corner_lon = scripfile['grid_corner_lon'][:,:].values
   corner_lat = scripfile['grid_corner_lat'][:,:].values
     
   ### check for logitude wrapping issues    
   # for n in range(ncol):
   #    for c in range(4):
   #       corner_lon[n,c] = np.round(corner_lon[n,c],8)
   #       if corner_lon[n,c]<0.01 or corner_lon[n,c]>359.09 :
   #          print(corner_lon[n,c])
   # exit()

   ### code for investigating precision issue
   ### there was a problem where checking exact equality was missing neighbors
   # cnt = 0
   # for n in range(ncol):
   #    for c in range(4):
   #       lat_tmp = corner_lat[n,c]
   #       lon_tmp = corner_lon[n,c]
   #       if lat_tmp>30 and lat_tmp<33 and lon_tmp>340 and lon_tmp<343 :
   #          msg = f'  n: {n}  c: {c}   corner_lat: {lat_tmp}   corner_lon: {lon_tmp}'
   #          if cnt==0:
   #             lat_base = corner_lat[n,c]
   #             lon_base = corner_lon[n,c]
   #          if lat_tmp == lat_base and lon_tmp==lon_base:
   #             msg += '  !!!!!!!!'
   #          print(msg)
   #          cnt += 1
   # exit()


   neighbors = np.full([ncol,8+1],-1)
   bearings  = np.full([ncol,8+1],-1)
   edge_flag = np.full([ncol,8+1],-1)

   pg.find_neighbors(lat.values, lon.values, corner_lat, corner_lon, neighbors, bearings, edge_flag)

   neighbors = neighbors[:,1:8+1]
   bearings  = bearings [:,1:8+1]
   edge_flag = edge_flag[:,1:8+1]
   
   # sort neighbors by bearing
   bear = bearings[:,:]
   bear = np.where(bear<0,bear+360,bear)
   neighbors_sorted = pg.sort_neighbors( neighbors[:,:], bear )
   bearings_sorted  = pg.sort_neighbors( bearings[:,:],  bear )
   edge_flag_sorted = pg.sort_neighbors( edge_flag[:,:], bear )

   # if 'pt_list' in locals(): 
   #    print(f' unsorted:')
   #    print(f'  n: {neighbors[pt_list[0],:]}')
   #    print(f'  b: {bearings[pt_list[0],:]}')
   #    print(f'  e: {edge_flag[pt_list[0],:]}')
   #    print(f' sorted:')
   #    print(f'  n: {neighbors_sorted[pt_list[0],:]}')
   #    print(f'  b: {bearings_sorted[pt_list[0],:]}')
   #    print(f'  e: {edge_flag_sorted[pt_list[0],:]}')

   pg.adjust_sorted_neighbors( ncol, neighbors_sorted, bearings_sorted, edge_flag_sorted )

   # if 'pt_list' in locals(): 
   #    print(f' post-adjustment:')
   #    print(f' n: {neighbors_sorted[pt_list[0],:]}')
   #    print()

   ### use unsorted neighbors
   # neighbors_sorted = neighbors[:,1:8+1]

   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   np.random.seed(rng_seed*1000 + c^2)

   if 'pt_list' not in locals():
      pts_ind = np.random.choice( np.arange(ncol), size=num_pts, replace=False )
      print(f'pts_ind: {pts_ind}')
   else:
      pts_ind = pt_list

   # if num_pts == 1: pts_ind = [0]

   for n in range(num_pts):
      xn = pts_ind[n]

      ip = c*num_pts + n

      tres = copy.deepcopy(res)

      xlat = lat[xn].values
      xlon = lon[xn].values

      # print the index, lat, and lon for each point
      print(f' n/idx: {n:4} / {xn:6}      lat/lon: {xlat:8.2f}  / {xlon:8.2f}')

      tres.mpMinLatF = np.max([xlat-dx,-90])
      tres.mpMaxLatF = np.min([xlat+dx, 90])
      tres.mpMinLonF = xlon-dx
      tres.mpMaxLonF = xlon+dx
      
      # tres.mpCenterLonF = xlon
      if tres.mpMinLonF<180 and tres.mpMaxLonF>180: tres.mpCenterLonF = xlon


      # print(f' plot bounds:   {tres.mpMinLatF} / {tres.mpMaxLatF}    {tres.mpMinLonF} / {tres.mpMaxLonF} ')

      plot[ip] = ngl.map(wks,tres)

      for nn in neighbors[xn,:] :
         nn = int(nn)
         cx,cy = corner_lon[nn,:], corner_lat[nn,:]
         
         if print_neighbor_coords: print(f'   n: {n:3} nn: {nn:4}  cx: {cx}  cy: {cy} ')
         
         xx = np.array([ cx[0], cx[1], cx[2], cx[3], cx[0]])
         yy = np.array([ cy[0], cy[1], cy[2], cy[3], cy[0]])
         dummy = ngl.add_polyline(wks,plot[ip],xx,yy,pgres)

      for nn,nnid in enumerate(neighbors_sorted[xn,:]) :
         i = int(nnid)
         dummy = ngl.add_text(wks, plot[ip], f'{nn+1}', lon[i], lat[i], txres)
         # if lon[i]==xlon and lat[i]==xlat: text = 'X'

      # dummy = ngl.add_text(wks, plot[ip], f'{0}', xlon, xlat, txres)


      subtitle_font_height = 0.01
      # hs.set_subtitles(wks, plot[ip], name[c], '', '', font_height=subtitle_font_height)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
# hs.set_plot_labels(wks, plot, font_height=0.01, justify='left')

layout = [num_case,num_pts]

# if num_case==1: layout = [int(np.ceil(num_pts/num_plot_cols)),num_plot_cols]


pnl_res = hs.setres_panel()
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.01
pnl_res.nglPanelYWhiteSpacePercent = 5

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------