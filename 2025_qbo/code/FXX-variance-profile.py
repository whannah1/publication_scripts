import os, ngl, copy, xarray as xr, numpy as np
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
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
def add_var(var_name,lev=-1): var.append(var_name); lev_list.append(lev)
##------------------------------------------------------------------------------
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01',       n='E3SM control',     c='red')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01', n='E3SM L72 smoothed',c='green')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01', n='E3SM L80 refined', c='blue')
#-------------------------------------------------------------------------------

# lev = np.array([30,50,75,100,125,150,200,250,300,350,400,450,500,550,600,650,700,750,800,825,850,875,900,925,950,975,1000])
# lev = np.array([50,100,150,200,300,400,500,600,700,750,800,850,875,900,925,975])
# lev = np.array([1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200])

# var = []
add_var('T')
# add_var('Q')
# add_var('RH')
add_var('U')
# add_var('V')
# add_var('OMEGA')
# add_var('THETA')
# add_var('CLOUD')
# add_var('CLDLIQ')
# add_var('CLDICE')
# add_var('QRL'); add_var('QRS')

# fields sugested by Yaga
add_var('UTGWSPEC')  # U tendency - gravity wave spectrum
add_var('BUTGWSPEC') # Beres U tendency - gravity wave spectrum
add_var('UTGWORO')   # U tendency - orographic gravity wave drag
add_var('BTAUE')     # Reynolds stress from Beres scheme for waves propagating East
add_var('BTAUW')     # Reynolds stress from Beres scheme for waves propagating West

# Beres U-tendency across 5 spectral bins
add_var('BUTEND1') # U tendency   c < -40
add_var('BUTEND2') # U tendency  -40 < c < -15
add_var('BUTEND3') # U tendency  -15 < c <  15
add_var('BUTEND4') # U tendency   15 < c <  40
add_var('BUTEND5') # U tendency   40 < c

# add_var('O3')


num_plot_col = len(var)
num_plot_col = 4


fig_type = "png"
fig_file = 'figs/FXX-variance-profile'


lat1,lat2 = -30,30
# lat1,lat2,lon1,lon2 = 15,25,360-60,360-50

# xlat,xlon,dy,dx = 60,120,10,10;
# if 'xlat' in locals(): lat1,lat2,lon1,lon2 = xlat-dy,xlat+dy,xlon-dx,xlon+dx


# htype,first_file,num_files = 'h0',0,2
htype,first_file,num_files = 'h0',0,20*12
# htype,first_file,num_files = 'h1',0,1
# htype,first_file,num_files = 'h2',0,0

plot_diff = False

use_height_coord   = True

omit_bot,bot_k     = False,-2
omit_top,top_k     = False,30

print_stats        = False
print_profile      = False

recalculate = True

tmp_file_head = 'data/variance.profile'

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)

if 'diff_mode' not in locals(): diff_mode = 0

if 'lev' not in vars(): lev = np.array([0])

wks = ngl.open_wks(fig_type,fig_file)
plot = [None]*(num_var)
res = hs.res_xy()
# res.vpWidthF = 0.4
# res.xyMarkLineMode = "MarkLines"
res.xyMarkerSizeF = 0.008
res.xyMarker = 16
res.xyLineThicknessF = 8
# res.tmYLLabelFontHeightF         = 0.008
# res.tmXBLabelFontHeightF         = 0.008

res.tmXBAutoPrecision = False
res.tmXBPrecision = 2

if not use_height_coord: 
   res.trYReverse = True
   # res.xyYStyle = 'Log'


def get_comp(case):
   comp = 'eam'
   if 'CESM' in case: comp = 'cam'
   return comp
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
data_list_list,lev_list_list = [],[]
for v in range(num_var):
   hc.printline()
   print(hc.tcolor.GREEN+'  var: '+var[v]+hc.tcolor.ENDC)
   data_list,lev_list = [],[]
   for c in range(num_case):
      print(hc.tcolor.CYAN+'    case: '+case[c]+hc.tcolor.ENDC)

      data_dir_tmp,data_sub_tmp = None, None
      # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
      if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      if case_sub[c] is not None: data_sub_tmp = case_sub[c]
      
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}'
      if 'lat1' in vars(): tmp_file = tmp_file + f'.lat1_{lat1}.lat2_{lat2}'
      if 'lon1' in vars(): tmp_file = tmp_file + f'.lon1_{lon1}.lon2_{lon2}'
      tmp_file = f'{tmp_file}.nc'
      print(f'      tmp_file: {tmp_file}')

      if recalculate :

         case_obj = he.Case( name=case[c], atm_comp=get_comp(case[c]), data_dir=data_dir_tmp, data_sub=data_sub_tmp  )

         #-------------------------------------------------------------------------
         # read the data   
         if 'lat1' in vars(): case_obj.lat1 = lat1; case_obj.lat2 = lat2
         if 'lon1' in vars(): case_obj.lon1 = lon1; case_obj.lon2 = lon2
         
         data = case_obj.load_data(var[v],htype=htype,first_file=first_file,num_files=num_files,lev=lev)
         area = case_obj.load_data('area',htype=htype,first_file=first_file,num_files=1,).astype(np.double)   
         Z = case_obj.load_data('Z3',htype=htype,first_file=first_file,num_files=num_files)

         hc.print_time_length(data.time,indent=' '*6)

         #-------------------------------------------------------------------------
         if 'time' in data.dims : 
            time = data.time
            hc.print_time_length(data.time,indent=' '*6)

         data = data.var(dim='time')
         Z = Z.mean(dim='time')

         #-------------------------------------------------------------------------
         # area weighted spatial average 
         if case[c]=='ERA5':
            data = ( (data*area).sum(dim=['lat','lon']) / area.sum(dim=['lat','lon']) )#.values 
         else:
            Z    = ( (Z   *area).sum(dim='ncol') / area.sum(dim='ncol') )#.values 
            data = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )#.values 

         #-------------------------------------------------------------------------
         # options for omitting top or bottom levels
         if omit_bot:
            data = data[:bot_k]
            Z = Z[:bot_k]

         if omit_top:
            data = data[top_k:]
            Z = Z[top_k:]

         #----------------------------------------------------------------------
         # write to temporary file
         #----------------------------------------------------------------------
         tmp_ds = xr.Dataset()
         tmp_ds['data'] = data
         tmp_ds['time'] = time
         if 'lat1' in vars(): tmp_ds['lat1']=lat1; tmp_ds['lat2']=lat2
         if 'lon1' in vars(): tmp_ds['lon1']=lon1; tmp_ds['lon2']=lon2
         tmp_ds['Z']    = Z
         print(f'      writing to file: {tmp_file}')
         tmp_ds.to_netcdf(path=tmp_file,mode='w')

      else:
         tmp_ds = xr.open_dataset( tmp_file )
         data = tmp_ds['data']
         time = tmp_ds['time']
         Z    = tmp_ds['Z']
      
      #-------------------------------------------------------------------------
      # print stuff
      if print_stats:
         hc.print_stat(data,name=f'{var[v]} after averaging',indent='    ',compact=True)
         # hc.print_stat(Z,name=f'Z after averaging',indent='    ')

      if print_profile:
         print()
         for xx in data.values: print(f'    {xx}')
         print()

      #-------------------------------------------------------------------------
      # append final data to list
      data_list.append( data.values )
      if use_height_coord:
         lev_list.append( Z.values )
      else:
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
   if diff_mode==1 :
      baseline1 = data_list[0]
      baseline2 = data_list[1]
   if diff_mode==2 :
      baseline1 = data_list[0]
      baseline2 = data_list[2]
   if plot_diff:
      for c in range(num_case): 
         if diff_mode==1 :
            if c==0 or c==2 : baseline = baseline1
            if c==1 or c==3 : baseline = baseline2
         if diff_mode==2 :
            if c==0 or c==1 : baseline = baseline1
            if c==2 or c==3 : baseline = baseline2
         data_list[c] = data_list[c] - baseline

   
   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])
   tres.trXMinF = data_min
   tres.trXMaxF = data_max
   ip = v

   for c in range(num_case):
      tres.xyLineColor   = clr[c]
      tres.xyMarkerColor = clr[c]
      tres.xyDashPattern = dsh[c]

      tplot = ngl.xy(wks, data_list[c], lev_list[c], tres)

      if diff_mode==0 :
         if (c==1 and plot_diff) or (c==0 and not plot_diff) :
            plot[ip] = tplot
         elif (plot_diff and c>0) or not plot_diff:
            ngl.overlay(plot[ip],tplot)

      if diff_mode==1 :
         if (c==2 and plot_diff) or (c==0 and not plot_diff) :
            plot[ip] = tplot
         elif (plot_diff and c!=0 and c!=1) or not plot_diff:
            ngl.overlay(plot[ip],tplot)

      if diff_mode==2 :
         if (c==1 and plot_diff) or (c==0 and not plot_diff) :
            plot[ip] = tplot
         elif (plot_diff and c!=0 and c!=2) or not plot_diff:
            ngl.overlay(plot[ip],tplot)

   ### add vertical line
   lres = hs.res_xy()
   lres.xyLineThicknessF = 1
   lres.xyDashPattern = 0
   lres.xyLineColor = 'black'
   ngl.overlay(plot[ip],ngl.xy(wks, np.array([0,0]), np.array([-1e3,1e8]), lres))

   reg_str = ''
   var_str = var[v]


   if 'lat1' in locals(): 
      lat1_str = f'{lat1}N' if lat1>=0 else f'{(lat1*-1)}S'
      lat2_str = f'{lat2}N' if lat2>=0 else f'{(lat2*-1)}S'
      reg_str += f' {lat1_str}:{lat2_str} '
   if 'lon1' in locals(): 
      lon1_str = f'{lon1}E' #if lon1>=0 and lon1<=360 else f'{(lon1*-1)}S'
      lon2_str = f'{lon2}E' #if lon2>=0 and lon2<=360 else f'{(lon2*-1)}S'
      reg_str += f' {lon1_str}:{lon2_str} '

   if plot_diff: var_str += ' (diff)'

   hs.set_subtitles(wks, plot[ip], var_str, '', reg_str, font_height=0.01)


#-------------------------------------------------------------------------------
# Add legend
#-------------------------------------------------------------------------------
if num_case>1:
   lgres = ngl.Resources()
   lgres.vpWidthF           = 0.05
   lgres.vpHeightF          = 0.08
   lgres.lgLabelFontHeightF = 0.012
   lgres.lgMonoDashIndex    = True
   lgres.lgLineLabelsOn     = False
   lgres.lgLineThicknessF   = 8
   lgres.lgLabelJust        = 'CenterLeft'
   lgres.lgLineColors       = clr
   lgres.lgDashIndexes      = dsh

   lx,ly = 0.5,0.45
   if num_var==2: lx,ly = 0.3,0.45
   if num_var==4: lx,ly = 0.05,0.5

   # pid = ngl.legend_ndc(wks, len(case_name), case_name, lx, ly, lgres)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

#-- draw a common title string on top of the panel
textres               =  ngl.Resources()
# textres.txFontHeightF =  0.01                  #-- title string size
# ngl.text_ndc(wks,f'time step = {ss_t}',0.5,.97,textres)  #-- add title to plot
textres.txFontHeightF =  0.02                  #-- title string size
if layout[0]==1: y_pos = 0.7
if layout[0]>=2: y_pos = 0.9
# ngl.text_ndc(wks,f'time step = {ss_t}',0.5,y_pos,textres)  #-- add title to plot

pres = hs.setres_panel()
pres.nglPanelTop      =  0.93

ngl.panel(wks,plot,layout,pres)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
