import os, ngl, string, copy, xarray as xr, numpy as np
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
data_dir,data_sub = None,None
#---------------------------------------------------------------------------------------------------
'''
file1=figs/throttle.profile.v1.32xM.png
file2=figs/throttle.profile.v1.fixed_NxM.png
montage ${file1} ${file2} -geometry 2048x500+5 -tile 1x2 figs/throttle.profile.v1.png
'''
#---------------------------------------------------------------------------------------------------
case_name,case,case_dir,case_sub,case_grid,clr,dsh,mrk = [],[],[],[],[],[],[],[]
remap_flag = []
def add_case(case_in,n=None,p=None,s=None,g=None,d=0,c='black',m=0,r=False):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   if n is None:
      tmp_name = ''
   else:
      tmp_name = n
   case.append(case_in); case_name.append(tmp_name)
   case_dir.append(p); case_sub.append(s); case_grid.append(g)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
   remap_flag.append(r)

var,varstr,vclr,vdsh = [],[],[],[]
def add_var(var_name,s=None,clr='black',dsh=0): 
   var.append(var_name); vclr.append(clr); vdsh.append(dsh)
   varstr.append(s)

lev = np.array([30,50,75,100,125,150,200,250,300,350,400,450,500,550,600,650,700,750,800,825,850,875,900,925,950,975,1000])
# lev = np.array([1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200])
# lev = np.array([50,100,150,200,300,400,500,600,700,750,800,850,875,900,925,975])
# lev = np.array([0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200])
# lev = np.array([1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,825,850,875,900,925,950,975,1000])

#---------------------------------------------------------------------------------------------------
### CRM throttling tests



# suffix = '32xM'; fig_lbl = list(string.ascii_lowercase)[0:]; subtitle_prefix = 'Subset 1: '
# add_case('E3SM.2023-THR-TEST-02.ne4pg2_oQU480.F2010-MMF1.NXY_32_1'   ,n='32x1', c='firebrick', d=0)
# add_case('E3SM.2023-THR-TEST-02.ne4pg2_oQU480.F2010-MMF1.NXY_32_4'   ,n='32x4', c='orangered', d=0)
# add_case('E3SM.2023-THR-TEST-02.ne4pg2_oQU480.F2010-MMF1.NXY_32_8'   ,n='32x8', c='orange', d=0)

suffix = 'fixed_NxM'; fig_lbl = list(string.ascii_lowercase)[4:]; subtitle_prefix = 'Subset 2: '
add_case('E3SM.2023-THR-TEST-02.ne4pg2_oQU480.F2010-MMF1.NXY_32_32'  ,n='32x32', c='deeppink', d=0)
add_case('E3SM.2023-THR-TEST-02.ne4pg2_oQU480.F2010-MMF1.NXY_128_8'  ,n='128x8', c='darkorchid', d=0)
add_case('E3SM.2023-THR-TEST-02.ne4pg2_oQU480.F2010-MMF1.NXY_256_4'  ,n='256x4', c='dodgerblue', d=0)

### suffix = 'Nx8'
### add_case('E3SM.2023-THR-TEST-02.ne4pg2_oQU480.F2010-MMF1.NXY_32_8'   ,n='32x8',  c='magenta', d=0)
### add_case('E3SM.2023-THR-TEST-02.ne4pg2_oQU480.F2010-MMF1.NXY_64_8'   ,n='64x8',  c='red', d=0)
### add_case('E3SM.2023-THR-TEST-02.ne4pg2_oQU480.F2010-MMF1.NXY_128_8'  ,n='128x8', c='green', d=0)
### add_case('E3SM.2023-THR-TEST-02.ne4pg2_oQU480.F2010-MMF1.NXY_256_8'  ,n='256x8', c='blue', d=0)
### add_case('E3SM.2023-THR-TEST-02.ne4pg2_oQU480.F2010-MMF1.NXY_512_8'  ,n='512x8', c='purple', d=0)

#---------------------------------------------------------------------------------------------------
add_var('CLDLIQ'  ,s='Cloud Liquid Amount')
add_var('CLDICE'  ,s='Cloud Ice Amount')
add_var('RH'      ,s='Relative Humidity')
add_var('MMF_MCUP',s='CRM Upward Mass Flux')

#---------------------------------------------------------------------------------------------------

if 'suffix' not in locals(): raise NameError('The variable "suffix" must be defined')

num_plot_col = len(var)

fig_file,fig_type = f'figs/throttle.profile.v1.{suffix}','png'

htype,years,months,first_file,num_files = 'h0',[],[],0,12

plot_diff = False

print_stats        = False

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var = len(var)
num_case = len(case)

if 'dsh' not in locals(): 
   if num_case>1 : dsh = np.zeros(num_case)
   else : dsh = [0]

if 'lev' not in vars(): lev = np.array([0])

wks = ngl.open_wks(fig_type,fig_file)
plot = [None]*(num_var)
res = hs.res_xy()
# res.vpWidthF = 0.4
# res.xyMarkLineMode = "MarkLines"
res.xyMarkerSizeF = 0.008
res.xyMarker = 16
res.xyLineThicknessF = 8
res.tmYLLabelFontHeightF = 0.02
res.tmXBLabelFontHeightF = 0.02
res.tiXAxisFontHeightF   = 0.02
res.tiYAxisFontHeightF   = 0.02

res.tmXBAutoPrecision = False
res.tmXBPrecision = 2

res.tiYAxisString = 'Pressure [hPa]'
res.trYReverse = True
# res.xyYStyle = 'Log'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
data_list_list,lev_list_list = [],[]
for v in range(num_var):
   hc.printline()
   print(hc.tcolor.GREEN+'  var: '+var[v]+hc.tcolor.ENDC)
   data_list,lev_list = [],[]
   for c in range(num_case):
      print(hc.tcolor.CYAN+'    case: '+case[c]+hc.tcolor.ENDC)
      #-------------------------------------------------------------------------
      data_dir_tmp,data_sub_tmp = None, None
      if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      if case_sub[c] is not None: data_sub_tmp = case_sub[c]
      case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp )
      #-------------------------------------------------------------------------
      tvar = var[v]
      # if tvar=='OMEGA' and 'pg2' in case[c] : tvar = 'DYN_OMEGA'
      # if tvar=='MMF_MCU': tvar = 'MMF_MCUUP'
      #-------------------------------------------------------------------------
      # read the data
      area = case_obj.load_data('area',htype=htype,first_file=first_file,num_files=num_files).astype(np.double)
      data = case_obj.load_data(tvar,  htype=htype,first_file=first_file,num_files=num_files,lev=lev)
      #-------------------------------------------------------------------------
      # # Handle derived variables
      # if tvar=='MMF_MCU': 
      #    data = data + case_obj.load_data('MMF_MCUDN',htype=htype,first_file=first_file,num_files=tnum_files,lev=lev)
      #-------------------------------------------------------------------------
      if var[v]=='CLDLIQ'  : data = data*1e3
      if var[v]=='CLDICE'  : data = data*1e3
      # if var[v]=='RH'      : 
      if var[v]=='MMF_MCUP': data = data*86400. # kg/m2/s => kg/m2/day
      #-------------------------------------------------------------------------
      # hc.print_time_length(data.time,indent=' '*6)
      # if print_stats: hc.print_stat(data,name=f'{var[v]} before averaging',indent=' '*4,compact=True)
      #-------------------------------------------------------------------------
      data = data.mean(dim='time')
      data = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )
      #-------------------------------------------------------------------------
      if print_stats: hc.print_stat(data,name=f'{var[v]} after averaging',indent=' '*4,compact=True)
      #-------------------------------------------------------------------------
      data_list.append( data.values )
      lev_list.append( data['lev'].values )

   data_list_list.append(data_list)
   lev_list_list.append(lev_list)


#-------------------------------------------------------------------------------
# Create plot - overlay all cases for each var
#-------------------------------------------------------------------------------
for v in range(num_var):
   data_list = data_list_list[v]
   lev_list = lev_list_list[v]
   
   tres = copy.deepcopy(res)

   # ip = v*num_case+c
   ip = c*num_var+v
   
   baseline = data_list[0]
   if plot_diff:
      for c in range(num_case): 
         data_list[c] = data_list[c] - baseline

   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])
   tres.trXMinF = data_min
   tres.trXMaxF = data_max

   # if var[v]=='CLDLIQ'  : tres.trXMinF, tres.trXMaxF = 0, 0.035
   # if var[v]=='CLDICE'  : tres.trXMinF, tres.trXMaxF = 0, 0.012
   # if var[v]=='RH'      : tres.trXMinF, tres.trXMaxF = 0, 84
   # if var[v]=='MMF_MCUP': tres.trXMinF, tres.trXMaxF = 0, 0.004

   if var[v]=='CLDLIQ'  : tres.trXMinF, tres.trXMaxF = 0, 35
   if var[v]=='CLDICE'  : tres.trXMinF, tres.trXMaxF = 0, 12
   if var[v]=='RH'      : tres.trXMinF, tres.trXMaxF = 0, 84
   if var[v]=='MMF_MCUP': tres.trXMinF, tres.trXMaxF = 0, 320

   # adjust the panels that were inconsistent

   # if var[v] in ['CLDICE','MMF_MCUP']:
   #    tres.vpWidthF  = 0.6+0.1
   #    tres.vpHeightF = 0.6+0.1
   # else:
   #    tres.vpWidthF  = 0.6-0.1
   #    tres.vpHeightF = 0.6-0.1

   # if var[v] in ['CLDICE','MMF_MCUP']:
   #    tres.tmXBLabelFontHeightF = tres.tmXBLabelFontHeightF*0.5

   if var[v]=='CLDLIQ'  : tres.tiXAxisString = '[g/kg]'
   if var[v]=='CLDICE'  : tres.tiXAxisString = '[g/kg]'
   if var[v]=='RH'      : tres.tiXAxisString = '[%]'
   if var[v]=='MMF_MCUP': tres.tiXAxisString = '[kg/m2/day]'

   ip = v

   for c in range(num_case):
      tres.xyLineColor   = clr[c]
      tres.xyMarkerColor = clr[c]
      tres.xyDashPattern = dsh[c]

      tplot = ngl.xy(wks, data_list[c], lev_list[c], tres)
      
      if (c==1 and plot_diff) or (c==0 and not plot_diff) :
         plot[ip] = tplot
      elif (plot_diff and c>0) or not plot_diff:
         ngl.overlay(plot[ip],tplot)

   ### add vertical line
   lres = hs.res_xy()
   lres.xyLineThicknessF = 1
   lres.xyDashPattern = 0
   lres.xyLineColor = 'black'
   ngl.overlay(plot[ip],ngl.xy(wks, np.array([0,0]), np.array([-1e3,1e8]), lres))
   
   # hs.set_subtitles(wks, plot[ip], '', ctr_str, var_str, font_height=0.008)
   hs.set_subtitles(wks, plot[ip], subtitle_prefix+varstr[v], '', '', font_height=0.008)


#-------------------------------------------------------------------------------
# Add legend
#-------------------------------------------------------------------------------
if num_case>1:
   lgres = ngl.Resources()
   lgres.vpWidthF           = 0.05
   lgres.vpHeightF          = 0.08
   lgres.lgLabelFontHeightF = 0.008
   lgres.lgMonoDashIndex    = True
   lgres.lgLineLabelsOn     = False
   lgres.lgLineThicknessF   = 10
   lgres.lgLabelJust        = 'CenterLeft'
   lgres.lgLineColors       = clr
   lgres.lgDashIndexes      = dsh

   labels = case_name
   for l,lbl in enumerate(labels): labels[l] = f'  {lbl}'

   pid = ngl.legend_ndc(wks, len(labels), labels, 0.13, 0.55, lgres)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

#-- draw a common title string on top of the panel
textres               =  ngl.Resources()
# textres.txFontHeightF =  0.01                  #-- title string size
# ngl.text_ndc(wks,f'time step = {ss_t}',0.5,.97,textres)  #-- add title to plot
# textres.txFontHeightF =  0.02                  #-- title string size
# if layout[0]==1: y_pos = 0.7
# if layout[0]>=2: y_pos = 0.9
# ngl.text_ndc(wks,f'time step = {ss_t}',0.5,y_pos,textres)  #-- add title to plot

pres = hs.setres_panel()
pres.nglPanelTop      =  0.93

pres.nglPanelFigureStrings            = fig_lbl
pres.nglPanelFigureStringsJust        = "TopRight"
pres.nglPanelFigureStringsFontHeightF = 0.01

pres.nglPanelXWhiteSpacePercent = 5
pres.nglPanelYWhiteSpacePercent = 5

ngl.panel(wks,plot,layout,pres)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
