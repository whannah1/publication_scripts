# plot the spectra of spatial variance for specific CRM columns
import os, ngl, sys, numba, copy, xarray as xr, numpy as np, re, string, cmocean
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import pg_checkerboard_utilities as pg
import xrft

case,name = [],[]
def add_case(case_in,n=''):
   global name,case
   case.append(case_in); name.append(n)
var,var_str,clr,dsh = [],[],[],[]
def add_var(var_name,n=None,c='black',d=0): 
   if n is None: n = var_name
   var.append(var_name); var_str.append(n); clr.append(c); dsh.append(d)
#-------------------------------------------------------------------------------
add_case('E3SM.CHX.ne30pg2_ne30pg2.F-MMF1.00', n='E3SM-MMF')
#-------------------------------------------------------------------------------

add_var('MMF_QTLS',c='red')
# add_var('MMF_DQ',  c='blue')

# ilat,ilon = 6,139
ilat,ilon = 10,143

num_files = 1
first_file = 365*4+180-int(num_files/2)

fig_type = 'png'
fig_file = 'figs/FXX-crm-forcing'

#---------------------------------------------------------------------------------------------------
# Set up plot
#---------------------------------------------------------------------------------------------------
num_case = len(case)
num_var  = len(var)

wkres = ngl.Resources()
npix=2048; wkres.wkWidth,wkres.wkHeight=npix,npix # use this for plotting all patterns w/ rotation
wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*9

plot_pos = [4,6,3,0,1,2,5,8,7] # use this to arrange plots according to their geographic position

res = hs.res_xy()
res.vpHeightF = 0.4
res.xyLineThicknessF = 8
# res.tmYLLabelFontHeightF = 0.008
# res.tmXBLabelFontHeightF = 0.008

lres = hs.res_xy()
lres.xyLineColor = 'red'
lres.xyLineThicknessF = 1

#---------------------------------------------------------------------------------------------------
# Load SCRIP data
#---------------------------------------------------------------------------------------------------
scrip_file_path = 'scrip_files/ne30pg2_scrip.nc'
scripfile = xr.open_dataset(scrip_file_path)
center_lat = scripfile['grid_center_lat'].values
center_lon = scripfile['grid_center_lon'].values
corner_lon = scripfile['grid_corner_lon'].values
corner_lat = scripfile['grid_corner_lat'].values
ncol = len(center_lat)

#---------------------------------------------------------------------------------------------------
# find adjacent neighbors around main point - see code/F09-crm-snapshot.py
#---------------------------------------------------------------------------------------------------

if ilat==6 and ilon==139: 
   mcol = 9246
   neighborhood = [9246, 9241, 9243, 9244, 9245, 9247, 9364, 9365, 9361]
   neighborbear = [  -1, -129,  -95,  179,  124,   84,    0,   48,  -56]

if ilat==10 and ilon==143:
   mcol = 9489
   neighborhood = [9489, 9370, 9371, 9488, 9374, 9492, 9491, 9494, 9490]
   neighborbear = [  -1, -133,  179,  -97,  125,   82,    0,   45,  -55]

neighborbear = np.array(neighborbear)
neighborhood = np.array(neighborhood)

tmp_neighborhood = neighborhood[1:]
neighborhood[1:] = tmp_neighborhood[ np.argsort(neighborbear[1:]) ]

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# for c in range(num_case):
c = 0

print('\n  case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)

case_obj = he.Case( name=case[c], time_freq='daily' )

#-------------------------------------------------------------------------------
# Load height data [ find height for subset
#-------------------------------------------------------------------------------

hgt_data = case_obj.load_data('Z3',   htype='h2',num_files=num_files,first_file=first_file).isel(ncol=neighborhood,time=0)

hgt_data = hgt_data.mean(dim='ncol')
hgt_data = hgt_data.sortby('lev', ascending=False).values / 1e3
nlev = -1
for k in range(len(hgt_data)):
   if hgt_data[k]>12: nlev = k

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

for v in range(num_var) :

   #-------------------------------------------------------------------------------
   # Load and subset CRM data - throw out data above 10km
   #-------------------------------------------------------------------------------

   print('    var: '+hc.tcolor.GREEN+f'{var[v]:10}'+hc.tcolor.ENDC)

   crm_data = case_obj.load_data(var[v],htype='h2',num_files=num_files,first_file=first_file).isel(ncol=neighborhood,crm_ny=0)

   crm_data = crm_data[:,:nlev,:,:]

   print(crm_data)

   #-------------------------------------------------------------------------------
   # Create plots
   #-------------------------------------------------------------------------------
   res.xyXStyle = "Log"
   res.trXMinF = np.min(x_coord)
   res.trXMaxF = np.max(x_coord)
   res.trYMaxF = np.max(data_fft.values)

   res.xyLineColor = clr[v]

   for i in range(len(plot)):
      ip = plot_pos[i]

      tplot = ngl.xy( wks, x_coord, data_fft[:,i].values, res ) 
      if v==0: plot[ip] = tplot
      if v>0: ngl.overlay(plot[ip], tplot)

      # Ocerlay Markov red noise spectra
      ngl.overlay(plot[ip], ngl.xy( wks, x_coord, markov, lres ) )

      # add subtitles
      tlat,tlon = center_lat[neighborhood[i]],center_lon[neighborhood[i]]
      hs.set_subtitles(wks, plot[ip],'',f'{tlat:5.1f}N {tlon:5.1f}E','', font_height=0.01)

#-------------------------------------------------------------------------------
# Finalize plot
#-------------------------------------------------------------------------------
pnl_res = hs.setres_panel()

ngl.panel(wks,plot,[3,3],pnl_res)

ngl.end()

hc.trim_png(fig_file)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
