import os, ngl, xarray as xr, numpy as np
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import copy, cftime, warnings
import cmocean
#---------------------------------------------------------------------------------------------------
case_name,case,case_dir,case_sub,case_grid,clr,dsh,mrk = [],[],[],[],[],[],[],[]
def add_case(case_in,n=None,p=None,s=None,g=None,d=0,c='black',m=0):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   if n is None:
      tmp_name = ''
   else:
      tmp_name = n
   case.append(case_in); case_name.append(tmp_name)
   case_dir.append(p); case_sub.append(s); case_grid.append(g)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
#-------------------------------------------------------------------------------
var_x,var_y,lev_list = [],[],[]
def add_vars(var_x_name,var_y_name,lev=None): 
   # if lev==-1: lev = np.array([0])
   var_x.append(var_x_name)
   var_y.append(var_y_name)
   lev_list.append(lev)
#-------------------------------------------------------------------------------
# add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.00',           n='E3SM-MMF',     c='red')
# add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.00',n='E3SM-MMF+BVT', c='green')
# add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_1.00',n='E3SM-MMF+FVT', c='blue')
add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.02',n='E3SM-MMF+BVT', c='green')
#-------------------------------------------------------------------------------

lev = np.array([10,30,50,75,100,125,150,200,250,300,350,400,450,500,
               550,600,650,700,750,800,825,850,875,900,925,950,975,1000])
# lev = np.array([5,10,30,50,100,150,200,300,400,500,600,700,800,850,925,975,1000])



from optparse import OptionParser
parser = OptionParser()
parser.add_option('-v',dest='var',default=None,help='variable to plot')
(opts, args) = parser.parse_args()


if opts.var is None:
   # add_vars('CWV','PRECT')
   
   # tmp_xvar = 'CWV'
   tmp_xvar = 'PRECT'

   add_vars(tmp_xvar,'MMF_VT_T')
   add_vars(tmp_xvar,'MMF_VT_TLS')
   add_vars(tmp_xvar,'MMF_VT_TEND_T')

   add_vars(tmp_xvar,'MMF_VT_Q')
   add_vars(tmp_xvar,'MMF_VT_QLS')
   add_vars(tmp_xvar,'MMF_VT_TEND_Q')

   # add_vars(tmp_xvar,'MMF_VT_U')
   # add_vars(tmp_xvar,'MMF_VT_TEND_U')
   # add_vars(tmp_xvar,'MMF_VT_ULS')

else:
   tvar = opts.var
   add_vars('CWV',f'MMF_VT_{tvar}')
   add_vars('CWV',f'MMF_VT_TEND_{tvar}')
   add_vars('CWV',f'MMF_VT_{tvar}LS')
   # add_vars('PRECT',opts.var)


htype,first_file,num_files = 'h1',1,365
# htype,first_file,num_files = 'h1',5,10

htype_x,htype_y = htype,htype

# lat1,lat2,lon1,lon2 = -5,5,-5,5
# lat1,lat2 = -5,5
# lat1,lat2 = -60,60

recalculate = False

temp_file_head = 'data/bin-VT'

fig_type,fig_file = 'png','figs/F08-bin-vt-tracer'

num_plot_col = 3

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var = len(var_x)
num_case = len(case)

# if 'clr' not in vars(): 
#    if num_case>1 : clr = np.linspace(2,len( ngl.retrieve_colormap(wks) )-1,num_case,dtype=int)
#    else : clr = ['black']

# if 'dsh' not in vars(): 
#    if num_case>1 : dsh = np.zeros(num_case)
#    else : dsh = [0]

wks = ngl.open_wks(fig_type,fig_file)
plot = [None]*(num_var*2)
res = hs.res_xy()
res.vpHeightF = 0.3
# res.tmYLLabelFontHeightF         = 0.008
# res.tmXBLabelFontHeightF         = 0.008
res.xyLineColors   = clr
res.xyDashPatterns = dsh

cres = hs.res_contour_fill()
cres.trYReverse = True
cres.tiXAxisFontHeightF     = 0.02
cres.tiYAxisFontHeightF     = 0.02
cres.tmYLLabelFontHeightF   = 0.02
cres.tmXBLabelFontHeightF   = 0.02
cres.lbLabelFontHeightF     = 0.02
cres.lbTitleFontHeightF     = 0.02
cres.lbTitlePosition        = 'Bottom'
cres.tiYAxisString = 'Pressure [hPa]'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_data_dir(c):
   global case_sub
   data_dir_tmp = None
   # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
   if case_dir[c] is not None: data_dir_tmp = case_dir[c]
   return data_dir_tmp

def get_data_sub(c):
   global case_dir
   data_sub_tmp = None
   if case_sub[c] is not None: data_sub_tmp = case_sub[c]
   return data_sub_tmp

def get_comp(case):
   comp = 'eam'
   if 'INCITE2019' in case: comp = 'cam'
   if 'RGMA' in case: comp = 'cam'
   if 'CESM' in case: comp = 'cam'
   if 'MAML' in case: comp = 'eam_0001'
   if 'E3SM.PI-CPL.v1.' in case: comp = 'cam'
   return comp
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   print('  X / Y : '+hc.tcolor.MAGENTA+var_x[v]+hc.tcolor.ENDC+' / '+hc.tcolor.MAGENTA+var_y[v]+hc.tcolor.ENDC)
   # print('  var_x: '+hc.tcolor.MAGENTA+var_x[v]+hc.tcolor.ENDC)
   # print('  var_y: '+hc.tcolor.MAGENTA+var_y[v]+hc.tcolor.ENDC)
   bin_list = []
   cnt_list = []
   val_list = []
   lev_list = []
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)
      
      tcase_x, tcase_y = case[c], case[c]   
      if case[c]=="OBS": tcase_x, tcase_y = obs_case_X, obs_case_Y

      case_obj_x = he.Case( name=tcase_x, atm_comp=get_comp(case[c]), data_dir=get_data_dir(c), data_sub=get_data_sub(c)  )
      case_obj_y = he.Case( name=tcase_y, atm_comp=get_comp(case[c]), data_dir=get_data_dir(c), data_sub=get_data_sub(c)  )

      tmp_file = f'{temp_file_head}.{case[c]}.{var_x[v]}.{var_y[v]}'
      if 'lat1' in vars(): tmp_file = tmp_file + f'.lat1_{lat1}.lat2_{lat2}'
      if 'lon1' in vars(): tmp_file = tmp_file + f'.lon1_{lat1}.lon2_{lat2}'
      tmp_file = tmp_file + f'.nc'

      print('\n    tmp_file: '+tmp_file+'\n')

      if recalculate :
         #-------------------------------------------------------------------------
         # read the data
         #-------------------------------------------------------------------------
         # if 'lev_list' in locals(): lev = lev_list[v]
         if 'lat1' in vars(): case_obj_x.lat1, case_obj_y.lat1 = lat1, lat1
         if 'lat2' in vars(): case_obj_x.lat2, case_obj_y.lat2 = lat2, lat2
         if 'lon1' in vars(): case_obj_x.lon1, case_obj_y.lon1 = lon1, lon1
         if 'lon2' in vars(): case_obj_x.lon2, case_obj_y.lon2 = lon2, lon2

         X = case_obj_x.load_data(var_x[v],htype=htype_x,first_file=first_file,num_files=num_files,lev=lev)
         Y = case_obj_y.load_data(var_y[v],htype=htype_y,first_file=first_file,num_files=num_files,lev=lev)

         if 'lev' in X.dims: exit('X variable needs to be 2D!')
         
         #-------------------------------------------------------------------------
         # Convert to daily mean
         #-------------------------------------------------------------------------
         X = X.resample(time='D').mean(dim='time')
         Y = Y.resample(time='D').mean(dim='time')

         # print stuff
         hc.print_stat(X,name=var_x[v],stat='naxsh',indent='    ',compact=True)
         hc.print_stat(Y,name=var_y[v],stat='naxsh',indent='    ',compact=True)

         #-------------------------------------------------------------------------
         #-------------------------------------------------------------------------
         if var_x[v]=='PRECT'    : bin_min, bin_max, bin_spc = 0, 120, 4
         if var_x[v]=='CWV'      : bin_min, bin_max, bin_spc = 2, 68, 2
         if var_x[v]=='FLNT'     : bin_min, bin_max, bin_spc = 120, 300, 5

         bin_ds = hc.bin_YbyX( Y, X, bin_min=bin_min, bin_max=bin_max, bin_spc=bin_spc, keep_lev=True )

         print('writing to file: '+tmp_file)
         bin_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         bin_ds = xr.open_dataset( tmp_file )

      # print(); print(bin_ds); print(); 
      # exit()

      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      if c==0:
         use_lev = False
         if 'lev' in bin_ds['bin_val'].dims: use_lev = True
         # # if use_lev and v==0: plot = [None]*(num_var*num_case*2)
         # if use_lev and v==0: plot = [None]*(num_var*num_case)

      if use_lev:
         lev_list.append(bin_ds['bin_val']['lev'].values)
         val_list.append(bin_ds['bin_val'].transpose().values)
         cnt_list.append( bin_ds['bin_cnt'].sum(dim='lev').values / bin_ds['bin_cnt'].sum().values )
      else:
         val_list.append(bin_ds['bin_val'].values)

      bin_list.append(bin_ds['bins'].values)
   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   if use_lev:
      cres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
      if not( 'TEND' in var_y[v] or 'LS' in var_y[v] or 'NET' in var_y[v] ) : 
         cres.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )

      cres.sfXArray   = bin_list[c]
      cres.sfYCStartV = np.min( lev_list[c] )
      cres.sfYCEndV   = np.max( lev_list[c] )

      # set units for color bar
      unit_str = ''
      if var_y[v]=='MMF_VT_T'      : unit_str = '[K~S~2~N~]'
      if var_y[v]=='MMF_VT_TLS'    : unit_str = '[K~S~2~N~/s]'
      if var_y[v]=='MMF_VT_TEND_T' : unit_str = '[K~S~2~N~/s]'
      if var_y[v]=='MMF_VT_Q'      : unit_str = '[kg~S~2~N~/kg~S~2~N~]'
      if var_y[v]=='MMF_VT_QLS'    : unit_str = '[kg~S~2~N~/kg~S~2~N~/s]'
      if var_y[v]=='MMF_VT_TEND_Q' : unit_str = '[kg~S~2~N~/kg~S~2~N~/s]'
      if var_y[v]=='MMF_VT_U'      : unit_str = '[m~S~2~N~/s~S~2~N~]'
      if var_y[v]=='MMF_VT_ULS'    : unit_str = '[m~S~2~N~/s~S~2~N~/s]'
      if var_y[v]=='MMF_VT_TEND_U' : unit_str = '[m~S~2~N~/s~S~2~N~/s]'
      cres.lbTitleString = unit_str

      if var_x[v]=='PRECT' : cres.tiXAxisString = 'Precipitation [mm/day]'


      # if any([s in var_y[v]] for s in ['TEND','LS']): 
      if 'TEND' in var_y[v] or 'LS' in var_y[v] or 'NET' in var_y[v] : 
         data_min = np.min([np.nanmin(d) for d in val_list])
         data_max = np.max([np.nanmax(d) for d in val_list])
         if var_y[v] in ['MMF_VT_TLS','MMF_VT_TEND_T']: data_min,data_max = -10,10
         if var_y[v] in ['MMF_VT_QLS','MMF_VT_TEND_Q']: data_min,data_max = -0.002,0.002
         num_clev,aboutZero = 21,True
         (cmin,cmax,cint) = ngl.nice_cntr_levels(data_min, data_max, 
                                    cint=None, max_steps=num_clev, 
                                    returnLevels=False, aboutZero=aboutZero )
         cres.cnLevels = np.linspace(cmin,cmax,num=num_clev)
         cres.cnLevelSelectionMode = "ExplicitLevels"
      else:
         cres.cnLevelSelectionMode = "AutomaticLevels"

      for c in range(num_case):
         var_str = var_y[v]
         if 'MMF_VT' in var_str: var_str = var_str.replace('MMF_VT','VT')
         if 'TEND' in var_str: var_str = var_str.replace('TEND_','')+'_TEND_CRM'
         if 'LS'   in var_str: var_str = var_str.replace('LS'   ,'')+'_TEND_GCM'
         if 'NET'  in var_str: var_str = var_str.replace('NET'  ,'')+'_TEND_NET'

         if use_lev and c==0 and v==0: plot = [None]*(num_var*num_case)
         ip1 = v*num_case+c
         plot[ip1] = ngl.contour(wks,np.ma.masked_invalid(val_list[c]),cres) 
         res.xyLineColors = 'black'
         hs.set_subtitles(wks, plot[ip1], '', '', var_str, font_height=0.015)

         ### alt version - add second plot for bin count
         # if use_lev and c==0 and v==0: plot = [None]*(num_var*num_case*2)
         # ip1, ip2 = v*num_case*2+c, v*num_case*2+c+num_var
         # plot[ip1] = ngl.contour(wks,np.ma.masked_invalid(val_list[c]),cres) 
         # res.xyLineColors = 'black'
         # plot[ip2] = ngl.xy(wks,bin_list[c],cnt_list[c],res) 
         # hs.set_subtitles(wks, plot[ip1], '', '', var_y[v], font_height=0.015)
         # hs.set_subtitles(wks, plot[ip2], '', '', var_y[v], font_height=0.015)
         
   else:
      if num_case==1:
         plot_data_x  = bin_list[0]
         plot_data_y1 = val_list[0]
         plot_data_y2 = cnt_list
      else:
         plot_data_x  = np.stack(bin_list)
         plot_data_y1 = np.stack(val_list)
         plot_data_y2 = np.stack(cnt_list)
      plot_data_y1 = np.ma.masked_invalid(plot_data_y1) # mask NaN values
      ip1, ip2 = v, v+num_var
      plot[ip1] = ngl.xy(wks,plot_data_x,plot_data_y1,res) 
      plot[ip2] = ngl.xy(wks,plot_data_x,plot_data_y2,res) 
      # hs.set_subtitles(wks, plot[ip1], '', '', var_y[v], font_height=0.015)
      # hs.set_subtitles(wks, plot[ip2], '', '', var_y[v], font_height=0.015)

#-------------------------------------------------------------------------------
# Add legend
#-------------------------------------------------------------------------------
# if num_case>1:
#    lgres = ngl.Resources()
#    lgres.vpWidthF           = 0.05
#    lgres.vpHeightF          = 0.08
#    lgres.lgLabelFontHeightF = 0.012
#    lgres.lgMonoDashIndex    = True
#    lgres.lgLineLabelsOn     = False
#    lgres.lgLineThicknessF   = 8
#    lgres.lgLabelJust        = 'CenterLeft'
#    lgres.lgLineColors       = clr
#    lgres.lgDashIndexes      = dsh

#    lx,ly = 0.2,0.9
#    # if num_var==2: lx,ly = 0.3,0.45
#    # if num_var==4: lx,ly = 0.05,0.5

#    pid = ngl.legend_ndc(wks, len(case_name), case_name, lx, ly, lgres)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

if use_lev:
   # layout = [num_case,num_var]
   # layout = [num_case*2,num_var]
   layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
else:
   layout = [2,num_var]

ngl.panel(wks,plot,layout,hs.setres_panel())
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------