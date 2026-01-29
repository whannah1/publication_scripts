import os, ngl, copy, xarray as xr, numpy as np
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import cmocean
cscratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl'
gscratch = '/global/cfs/cdirs/m3312/whannah/2022-QBO-TEST'
pscratch = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu'
#-------------------------------------------------------------------------------
case_dir,case_sub = [],[]
case,name = [],[]
clr,dsh,mrk = [],[],[]
def add_case(case_in,n=None,p=None,s=None,d=0,c='black',m=0):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); name.append(n)
   case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
var,lev_list = [],[]
# def add_var(var_name,lev=-1): var.append(var_name); lev_list.append(lev)
##------------------------------------------------------------------------------
add_case('E3SM.BTAU-TEST-00.F2010.ne30pg2',n='E3SM',c='red',p=pscratch,s='run')
#-------------------------------------------------------------------------------

### define pressure levels for plotting - initial calculations stay on native grid
plev = [100,50,30,10]
# plev = [100]

var = ['BTAUXS']

num_plot_col = 2

fig_file,fig_type = 'figs/FXX-BTAUX-output-chk','png'


lat1,lat2 = -10,10


htype,first_file,num_files = 'h0',0,1

recalculate = True

plot_diff = True
common_colorbar = False

print_stats      = False

if num_files==12:
   tmp_file_head = 'data/BTAUX-spectrum-1yr'
else:
   tmp_file_head = 'data/BTAUX-spectrum'


tvar_list1,tau_list = [],[]
for n in reversed(range(32)): tvar_list1.append(f'BTAUXSn{n+1:02d}'); tau_list.append(2.5*(n+1)*-1)
for n in range(0,32+1):       tvar_list1.append(f'BTAUXSp{n:02d}');   tau_list.append(2.5*n)

data_shape  = (len(tau_list))
data_coords = [('bins', np.array(tau_list))]

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_plev = len(case),len(plev)

if 'diff_mode' not in locals(): diff_mode = 0
if 'plev' not in locals(): plev = np.array([0])

wks = ngl.open_wks(fig_type,fig_file)
plot = [None]*(num_plev)
res = hs.res_xy()
# res.vpWidthF = 0.4
# res.xyMarkLineMode = "MarkLines"
# res.xyMarkerSizeF = 0.008
# res.xyMarker = 16
# res.xyLineThicknessF = 8
# res.tmYLLabelFontHeightF         = 0.008
# res.tmXBLabelFontHeightF         = 0.008
# res.tmXBAutoPrecision = False
# res.tmXBPrecision = 2
res.trXMinF = -40
res.trXMaxF =  40
res.tiXAxisFontHeightF  = 0.02
res.tiYAxisFontHeightF  = 0.02
res.tiYAxisString       = 'Convective GW Tau [Pa]'
res.tiXAxisString       = 'Gravity Wave Phase Speed [m/s/day]'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for c in range(num_case):
   hc.printline()
   print(hc.tcolor.CYAN+'    case: '+case[c]+hc.tcolor.ENDC)
   #-------------------------------------------------------------------------
   data_dir_tmp,data_sub_tmp = None, None
   if case_dir[c] is not None: data_dir_tmp = case_dir[c]
   if case_sub[c] is not None: data_sub_tmp = case_sub[c]
   case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp  )
   if 'lat1' in locals(): case_obj.lat1 = lat1; case_obj.lat2 = lat2
   if 'lon1' in locals(): case_obj.lon1 = lon1; case_obj.lon2 = lon2
   #-------------------------------------------------------------------------
   area = case_obj.load_data('area',htype=htype,first_file=first_file,num_files=1,).astype(np.double)   
   #-------------------------------------------------------------------------
   for l in range(num_plev):
      print(hc.tcolor.GREEN+f'  plev: {plev[l]}'+hc.tcolor.ENDC)
      # data_list,bin_list = [],[]
      data1 = xr.DataArray(np.zeros(data_shape), coords=data_coords )
      data2 = xr.DataArray(np.zeros(data_shape), coords=data_coords )
      #-------------------------------------------------------------------------
      # read the post-interpolated spectral GW data
      for t,tvar in enumerate(tvar_list1):
         data_tmp = case_obj.load_data(tvar,lev=plev[l],htype=htype,first_file=first_file,num_files=num_files)
         data1[t] = ( (data_tmp.isel(lev=0).mean(dim='time')*area).sum(dim='ncol') / area.sum(dim='ncol') ).values
      #-------------------------------------------------------------------------
      # read the online-interpolated spectral GW data
      tvar_list2 = []
      for n in reversed(range(32)): tvar_list2.append(f'BTAUXSn{n+1:02d}_{plev[l]}mb')
      for n in range(0,32+1):       tvar_list2.append(f'BTAUXSp{n:02d}_{plev[l]}mb')
      for t,tvar in enumerate(tvar_list2):
         data_tmp = case_obj.load_data(tvar,htype=htype,first_file=first_file,num_files=num_files)
         data2[t] = ( (data_tmp.mean(dim='time')*area).sum(dim='ncol') / area.sum(dim='ncol') )
      #-------------------------------------------------------------------------
      # convert from m/s/s to m/s/day
      data1 = data1*86400.
      data2 = data2*86400.
      #-------------------------------------------------------------------------
      # print stuff
      if print_stats: hc.print_stat(data1,name=f'{var[v]} after averaging',indent='    ',compact=True)
      if print_stats: hc.print_stat(data2,name=f'{var[v]} after averaging',indent='    ',compact=True)
      #-------------------------------------------------------------------------
      # Create plot - overlay all cases for each var
      tres = copy.deepcopy(res)

      # if plot_diff:
      #    for c in range(1,num_case): 
      #       data_list[c] = data_list[c] - data_list[0]
      # diff_mag_max = np.max(np.absolute(data_list[1:]))

      # if (plot_diff and c==0) or not plot_diff:
      #    tres.cnLevelSelectionMode = 'ExplicitLevels'
      #    # tres.cnLevels = np.arange(-20,20+2,2)*1e-5 # use for units = m/s/s
      #    tres.cnLevels = np.arange(-18,18+2,2) # use for units = m/s/day
      #    # tres.cnFillPalette = "MPL_viridis"
      #    tres.cnFillPalette = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )
      
      # if plot_diff and c>0:
      #    tres.cnLevelSelectionMode = 'ExplicitLevels'
      #    tres.cnLevels = np.linspace(-1*diff_mag_max,diff_mag_max,21)
      #    tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )

      # tres.sfXArray = bin_list[c]
      # tres.sfYArray = lev_list[c]

      # ip = v*num_case+c
      ip = c*num_plev+l

      tres.xyLineColors = ['blue','red']
      tres.xyDashPatterns = [0,1]
      

      plot[ip] = ngl.xy(wks, np.array(tau_list), np.stack([data1,data2]), tres)

      hs.set_subtitles(wks, plot[ip], '', '', f'plev={plev[l]}', font_height=0.01)

      # add vertical line
      lres = hs.res_xy()
      lres.xyLineThicknessF = 1
      lres.xyDashPattern = 0
      lres.xyLineColor = 'black'
      ngl.overlay(plot[ip],ngl.xy(wks, np.array([0,0]), np.array([-1e3,1e8]), lres))


      

   # reg_str = ''
   # var_str = var[v]

   # if 'lat1' in locals(): 
   #    lat1_str = f'{lat1}N' if lat1>=0 else f'{(lat1*-1)}S'
   #    lat2_str = f'{lat2}N' if lat2>=0 else f'{(lat2*-1)}S'
   #    reg_str += f' {lat1_str}:{lat2_str} '
   # if 'lon1' in locals(): 
   #    lon1_str = f'{lon1}E' #if lon1>=0 and lon1<=360 else f'{(lon1*-1)}S'
   #    lon2_str = f'{lon2}E' #if lon2>=0 and lon2<=360 else f'{(lon2*-1)}S'
   #    reg_str += f' {lon1_str}:{lon2_str} '

   # if plot_diff: var_str += ' (diff)'

   # hs.set_subtitles(wks, plot[ip], var_str, '', reg_str, font_height=0.01)


#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

pres = hs.setres_panel()
pres.nglPanelTop      =  0.93

if common_colorbar: 
   pres.nglPanelLabelBar = True
   pres.lbTitleFontHeightF = 0.01
   pres.nglPanelLabelBarLabelFontHeightF = 0.01
   pres.lbTitlePosition = 'Bottom'
   pres.lbTitleString = '' # what are the units of "Beres Tau"?

ngl.panel(wks,plot,layout,pres)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
