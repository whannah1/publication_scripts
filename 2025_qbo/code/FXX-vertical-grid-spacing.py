import os, numpy as np, xarray as xr, ngl, copy
import hapy_common as hc
home = os.getenv('HOME')

vert_file_list,name,clr,dsh = [],[],[],[]


def add_grid(grid_in,n=None,d=0,c='black'):
   global vert_file_list,name,clr,dsh
   vert_file_list.append(grid_in); name.append(n); dsh.append(d); clr.append(c)


grid_path = f'{home}/E3SM/vert_grid_files'

add_grid(f'{grid_path}/L72_E3SM.nc',                                         n='L72',           c='gray')
# add_grid(f'{grid_path}/L72_E3SM_smooth.nsmooth_40.upper_lev_only.nc',        n='L72 smoothed',  c='red')
# add_grid(f'{grid_path}/L72_E3SM_refine+smooth.nsmooth_40.upper_lev_only.nc', n='L72 refined',   c='blue')
# add_grid(f'{grid_path}/L72_E3SM_refine-limit.nc',                            n='L72 limit dzdz',c='pink')
# add_grid(f'{grid_path}/L72_E3SM_refine-scale.nc',                            n='L72 scaled',    c='purple')

add_grid(f'{grid_path}/L72_E3SM_refine+smooth.nsmooth_40.upper_lev_only.nc', n='L80',   c='cyan', d=1)
# add_grid(f'{grid_path}/2023_refined_grid_test_L72+00.nsmooth_00.nc',n='L72+00',c='black')
# add_grid(f'{grid_path}/2023_refined_grid_test_L72+05.nsmooth_00.nc',n='L72+05',c='red')
# add_grid(f'{grid_path}/2023_refined_grid_test_L72+10.nsmooth_00.nc',n='L72+10',c='green')
# add_grid(f'{grid_path}/2023_refined_grid_test_L72+15.nsmooth_00.nc',n='L72+10',c='blue')
# add_grid(f'{grid_path}/2023_refined_grid_test_L72+20.nsmooth_00.nc',n='L72+20',c='purple')
# add_grid(f'{grid_path}/2023_refined_grid_test_L72+25.nsmooth_00.nc',n='L72+25',c='pink')

add_grid(f'{grid_path}/2023_refined_grid_test_L72+08.nsmooth_40.nc',n='L72+05',c='magenta')

# add_grid(f'{grid_path}/2023_refined_grid_test_L72+00.nsmooth_40.nc',n='L72+00',c='black')
# add_grid(f'{grid_path}/2023_refined_grid_test_L72+05.nsmooth_40.nc',n='L72+05',c='red')
# add_grid(f'{grid_path}/2023_refined_grid_test_L72+10.nsmooth_40.nc',n='L72+10',c='green')
# add_grid(f'{grid_path}/2023_refined_grid_test_L72+15.nsmooth_40.nc',n='L72+10',c='blue')
# add_grid(f'{grid_path}/2023_refined_grid_test_L72+20.nsmooth_40.nc',n='L72+20',c='purple')
# add_grid(f'{grid_path}/2023_refined_grid_test_L72+25.nsmooth_40.nc',n='L72+25',c='pink')

# cscratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl'
# gscratch = '/global/cfs/cdirs/m3312/whannah/2022-QBO-TEST'
# pscratch = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu'
# tcase='E3SM.QBO-TEST-02.F2010.ne30pg2.L72-nsu40';   add_grid(f'{pscratch}/{tcase}/run/{tcase}.eam.h0.0001-01.nc',n='L72',     c='black')
# tcase='E3SM.QBO-TEST-02.F2010.ne30pg2.L72.GW-MOD-0';add_grid(f'{pscratch}/{tcase}/run/{tcase}.eam.h0.0001-01.nc',n='L72 smth',c='red')

fig_file,fig_type = 'figs/FXX-vertical-grid-spacing','png'

print_table = False
use_height = False # for Y-axis, or else use pressure
add_zoomed_plot = False
add_refine_box = False

# set limits for second plot zoomed in on lower levels
# zoom_top_idx = -10
zoom_top_idx = -30

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

res.xyLineThicknessF             = 6.

res.xyMarkLineMode = 'MarkLines'
res.xyMarkerSizeF = 0.005
res.xyMarker = 16

if use_height:
  res.trYReverse = False
else:
  res.trYReverse = True

# res.xyYStyle = "Log"

#-------------------------------------------------------------------------------
# print stuff
#-------------------------------------------------------------------------------
if print_table:

  mlev_list, zlev_list = [],[]
  for f,vert_file in enumerate(vert_file_list):
    ds = xr.open_dataset(vert_file)
    # mlev = ds['hyam'].values*1000 + ds['hybm'].values*1000
    mlev = ds['hyai'].values*1000 + ds['hybi'].values*1000
    # zlev = ???
    mlev_list.append(mlev)
    # zlev_list.append(zlev)

  max_len = 0
  for mlev in mlev_list: 
    if len(mlev)>max_len: max_len = len(mlev)

  for k in range(max_len):
    msg = f'{k:3}  '
    # k2 = max_len-k-1
    # msg = f'{k:3}  ({k2:3})'
    for mlev in mlev_list: 
      if k < len(mlev):
        k2 = len(mlev)-k#-1
        msg += f'    ({k2:3})    {mlev[k]:8.1f}'
      else:
        msg += ' '*(12+4+5)
    print(msg)

  # exit()

#-------------------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------------------
mlev_list = []
dlev_list = []
for f,vert_file in enumerate(vert_file_list):

  ds = xr.open_dataset(vert_file)

  mlev = ds['hyam'].values*1000 + ds['hybm'].values*1000
  ilev = ds['hyai'].values*1000 + ds['hybi'].values*1000

  # rough estimate of height from pressure
  ilevz = np.log(ilev/1e3) * -6740.
  mlevz = np.log(mlev/1e3) * -6740.

  print()
  hc.print_stat(mlev, name=name[f]+' mlev',stat='nxh',indent='    ',compact=True)
  hc.print_stat(mlevz,name=name[f]+' zlev',stat='nxh',indent='    ',compact=True)

  
  dlevz = mlevz*0.
  for k in range(len(mlev)): dlevz[k] = ilevz[k] - ilevz[k+1]

  if use_height:
    mlev_list.append(mlevz)
    dlev_list.append(dlevz)
  else:
    mlev_list.append(mlev)
    dlev_list.append(dlevz)

#-------------------------------------------------------------------------------
# Create plot
#-------------------------------------------------------------------------------
dlev_min = np.min([np.nanmin(d) for d in dlev_list])
dlev_max = np.max([np.nanmax(d) for d in dlev_list])

mlev_min = np.min([np.nanmin(d) for d in mlev_list])
mlev_max = np.max([np.nanmax(d) for d in mlev_list])

# set limits for second plot w/ linear scale
dlev_min_2 = np.min([np.nanmin(d[zoom_top_idx:]) for d in dlev_list])
dlev_max_2 = np.max([np.nanmax(d[zoom_top_idx:]) for d in dlev_list])
mlev_min_2 = np.min([np.nanmin(d[zoom_top_idx:]) for d in mlev_list])
mlev_max_2 = np.max([np.nanmax(d[zoom_top_idx:]) for d in mlev_list])


for f,vert_file in enumerate(vert_file_list):

  mlev = mlev_list[f]
  dlev = dlev_list[f]

  if use_height:
    res.tiXAxisString = 'Grid Spacing [m]'
    res.tiYAxisString = 'Approx. Height [m]'
  else:
    res.tiXAxisString = 'Grid Spacing [m]'
    res.tiYAxisString = 'Approx. Pressure [hPa]'

  res.xyDashPattern = dsh[f]
  res.xyLineColor = clr[f]
  res.xyMarkerColor = clr[f]

  tres1 = copy.deepcopy(res)
  tres2 = copy.deepcopy(res)

  tres1.trXMinF = dlev_min
  tres1.trXMaxF = dlev_max + (dlev_max-dlev_min)*0.05
  tres1.trYMinF = mlev_min - mlev_min/2
  tres1.trYMaxF = 1e3 #mlev_max

  tres2.trXMinF = dlev_min_2
  tres2.trXMaxF = dlev_max_2 + (dlev_max_2-dlev_min_2)*0.05
  tres2.trYMinF = mlev_min_2 #- mlev_min_2/2
  tres2.trYMaxF = 1e3 #mlev_max_2

  if use_height: 
    tres1.xyYStyle = "Linear"
    tres2.xyYStyle = "Linear"
  else:
    tres1.xyYStyle = "Log"
    tres2.xyYStyle = "Linear"

  tplot1 = ngl.xy(wks, dlev, mlev, tres1)
  
  ### add plot zoomed in on lowest levels
  if add_zoomed_plot: 
    tplot2 = ngl.xy(wks, dlev, mlev, tres2)

  if f==0:
    plot.append(tplot1)
    if add_zoomed_plot: plot.append(tplot2)
  else:
    ngl.overlay(plot[0],tplot1)
    # if len(plot)>=2: ngl.overlay(plot[1],tplot2)
    if add_zoomed_plot: ngl.overlay(plot[1],tplot2)

#-------------------------------------------------------------------------------
# add lines to visually see how grids line up with first case (i.e. control/default)
#-------------------------------------------------------------------------------
# if len(plot)>1:
#   tres3 = copy.deepcopy(tres2)
#   tres3.xyDashPattern = 1
#   tres3.xyLineColor = 'black'
#   tres3.xyLineThicknessF = 1.
#   mlev = mlev_list[0]
#   for k in range(10):
#     kk = len(mlev)-1-k
#     xx = np.array([ -1e5, 1e5 ])
#     yy = np.array([ mlev[kk], mlev[kk] ])
#     ngl.overlay(plot[1], ngl.xy(wks, xx, yy, tres3) )

#-------------------------------------------------------------------------------
# indicate refinement levels
#-------------------------------------------------------------------------------
# if add_refine_box:
#   pgres = ngl.Resources()
#   pgres.nglDraw                = False
#   pgres.nglFrame               = False
#   pgres.gsLineColor            = 'black'
#   pgres.gsLineThicknessF       = 10.0
#   pgres.gsFillIndex            = 0
#   pgres.gsFillColor            = 'red'
#   pgres.gsFillOpacityF         = 0.3

#   rx1,rx2 = 0,1e6
#   rx = [rx1,rx2,rx2,rx1,rx1]

#   rz1,rz2 = 10e3,45e3
#   if use_height: 
#     rk1,rk2 = rz1,rz2
#   else:
#     rk1,rk2 = rp1,rp2
#   ry = [rk1,rk1,rk2,rk2,rk1]

#   pdum = ngl.add_polygon(wks, plot[0], rx, ry, pgres)

#-------------------------------------------------------------------------------
# Add legend
#-------------------------------------------------------------------------------
lgres = ngl.Resources()
lgres.vpWidthF           = 0.1
lgres.vpHeightF          = 0.2
lgres.lgLabelFontHeightF = 0.012
lgres.lgMonoDashIndex    = True
lgres.lgLineLabelsOn     = False
lgres.lgLineThicknessF   = 20
lgres.lgLabelJust        = 'CenterLeft'
lgres.lgLineColors       = clr[::-1]

for n in range(len(name)): name[n] = ' '+name[n]

if add_zoomed_plot:
  lpx, lpy = 0.3, 0.45
  # lpx, lpy = 0.8, 0.45
else:
  lpx, lpy = 0.6, 0.3

pid = ngl.legend_ndc(wks, len(name), name[::-1], lpx, lpy, lgres)

#-------------------------------------------------------------------------------
# Finalize plot
#-------------------------------------------------------------------------------
ngl.panel(wks,plot[0:len(plot)],[1,len(plot)],ngl.Resources()); ngl.end()

# trim white space
fig_file = f'{fig_file}.{fig_type}'
os.system(f'convert -trim +repage {fig_file}   {fig_file}')
fig_file = fig_file.replace(f'{home}/E3SM/','')
print(f'\n{fig_file}\n')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
