import os, numpy as np, xarray as xr, ngl, copy
import hapy_setres as hs
home = os.getenv('HOME')

grid_file_list,topo_file_list,name = [],[],[]
scratch_dir = '/global/cfs/cdirs/e3sm/inputdata'
# scratch_dir = '/gpfs/alpine/cli115/scratch/hannah6/inputdata'
topo_dir = f'{scratch_dir}/atm/cam/topo'


# grid_file_list.append(f'{home}/E3SM/data_grid/ne30np4_scrip.nc');      name.append('ne30np4')
# grid_file_list.append(f'{home}/E3SM/data_grid/ne30pg2_scrip.nc');      name.append('ne30pg2')
# grid_file_list.append(f'{home}/E3SM/data_grid/ne30pg3_scrip.nc');      name.append('ne30pg3')
# grid_file_list.append(f'{home}/E3SM/data_grid/ne30pg4_scrip.nc');      name.append('ne30pg4')

grid_file_list.append(f'{home}/E3SM/data_grid/ne30pg2_scrip.nc');      name.append('ne30pg2')
grid_file_list.append(f'{home}/E3SM/data_grid/conusx4v1pg2_scrip.nc'); name.append('RRM pg2')

for grid in grid_file_list:
    if 'ne30np4'  in grid: topo_file_list.append(f'{topo_dir}/USGS-gtopo30_ne30np4_16xdel2-PFC-consistentSGH.nc')
    if 'ne30pg2'  in grid: topo_file_list.append(f'{topo_dir}/USGS-gtopo30_ne30np4pg2_16xdel2.c20200108.nc')
    if 'ne30pg3'  in grid: topo_file_list.append(f'{topo_dir}/USGS-gtopo30_ne30np4pg3_16xdel2.c20200504.nc')
    if 'ne30pg4'  in grid: topo_file_list.append(f'{topo_dir}/USGS-gtopo30_ne30np4pg4_16xdel2.c20200504.nc')
    if 'x4v1_'    in grid: topo_file_list.append(f'{topo_dir}/USGS_conusx4v1-tensor12x_consistentSGH_c150924.nc')
    if 'x4v1pg2_' in grid: topo_file_list.append(f'{topo_dir}/USGS_conusx4v1pg2_12x_consistentSGH_20200609.nc')

fig_type = 'png'
fig_file = 'figs/F07-grid-comparison-2'

#-------------------------------------------------------------------------------
# Set up workstation
#-------------------------------------------------------------------------------
wks = ngl.open_wks(fig_type,fig_file)
plot = []
res = ngl.Resources()
res.nglDraw,res.nglFrame         = False,False
res.tmXTOn                       = False
res.tmXBMajorOutwardLengthF      = 0.
res.tmXBMinorOutwardLengthF      = 0.
res.tmYLMajorOutwardLengthF      = 0.
res.tmYLMinorOutwardLengthF      = 0.
res.tmYLLabelFontHeightF         = 0.015
res.tmXBLabelFontHeightF         = 0.015
res.tiXAxisFontHeightF           = 0.015
res.tiYAxisFontHeightF           = 0.015
res.tmXBMinorOn,res.tmYLMinorOn  = False,False
res.tmXBOn,res.tmYLOn  = False,False

res.cnFillOn                     = True
res.cnLinesOn                    = False
res.cnLineLabelsOn               = False
res.cnInfoLabelOn                = False
res.lbOrientation                = "Horizontal"
res.lbLabelFontHeightF           = 0.008
res.mpGridAndLimbOn              = False


res.cnFillMode      = 'CellFill'
res.cnCellFillEdgeColor = "black"

res.lbLabelBarOn = False
res.cnFillPalette   = "OceanLakeLandSnow"
res.cnLevelSelectionMode = "ExplicitLevels"
res.cnLevels = np.arange(0,4800+100,100)

# res.cnFillPalette   = "OceanLakeLandSnow"


### turn off color fill
# res.cnFillOpacityF = 0.0

res.mpLimitMode     = 'LatLon'
# CONUS
# res.mpMinLatF,res.mpMaxLatF      =  12, 70
# res.mpMinLonF,res.mpMaxLonF      =  360-150, 360-50

# zoomed-out CONUS
res.mpMinLatF,res.mpMaxLatF      =  10, 72
res.mpMinLonF,res.mpMaxLonF      =  360-160, 360-42
res.cnLevels = np.arange(0,3000+100,100)

# # zoomed-in CONUS
# res.mpMinLatF,res.mpMaxLatF      =  25, 50
# res.mpMinLonF,res.mpMaxLonF      =  360-125, 360-75
# res.cnLevels = np.arange(0,2500+100,100)

# # West CONUS
# res.mpMinLatF,res.mpMaxLatF      =  30, 60
# res.mpMinLonF,res.mpMaxLonF      =  360-150, 360-100

#-------------------------------------------------------------------------------
# Load data and create plot
#-------------------------------------------------------------------------------
for f,(grid_file,topo_file) in enumerate( zip(grid_file_list,topo_file_list) ):

  ds_grid = xr.open_dataset(grid_file)
  ds_topo = xr.open_dataset(topo_file)

  topo = ds_topo['PHIS'].values
  topo = topo / 9.81
  land = ds_topo['LANDFRAC'].values
  topo = np.where(land>0.5,topo,-1e3)

  tres = copy.deepcopy(res)
  tres.sfXArray      = ds_grid['grid_center_lon'].values
  tres.sfYArray      = ds_grid['grid_center_lat'].values
  tres.sfXCellBounds = ds_grid['grid_corner_lon'].values
  tres.sfYCellBounds = ds_grid['grid_corner_lat'].values

  # res.tiXAxisString = 'normalized level index'
  # res.tiYAxisString = 'lev [mb]'

  plot.append( ngl.contour_map(wks, topo, tres) )

  hs.set_subtitles(wks, plot[len(plot)-1], '', name[f], '', font_height=0.02)
  
#-------------------------------------------------------------------------------
# Finalize plot
#-------------------------------------------------------------------------------
pres = ngl.Resources()
pres.nglPanelLabelBar = True
pres.nglPanelFigureStrings            = ['a','b','c','d','e','f','g','h']
pres.nglPanelFigureStringsJust        = "TopLeft"
pres.nglPanelFigureStringsFontHeightF = 0.015
pres.nglPanelYWhiteSpacePercent = 10    
pres.nglPanelXWhiteSpacePercent = 5

# layout = [len(plot),1]
layout = layout = [2,np.ceil(len(plot)/2)]
ngl.panel(wks,plot,layout,pres); ngl.end()

# trim white space
fig_file = f'{fig_file}.{fig_type}'
os.system(f'convert -trim +repage {fig_file}   {fig_file}')
print(f'\n{fig_file}\n')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------