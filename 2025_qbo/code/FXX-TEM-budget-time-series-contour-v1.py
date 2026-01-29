import os, ngl, glob, subprocess as sp, numpy as np, xarray as xr, copy, string, dask
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
import cmocean
scratch = '/global/cfs/cdirs/m3312/whannah/2022-QBO-TEST'
#-------------------------------------------------------------------------------
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
var = []
def add_var(var_name,lev=-1): var.append(var_name)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

remap_str,search_str = 'remap_90x180','h0.tem.'; first_file,num_files = 12*0,12*20

### Runs for QBO paper
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01',       n='E3SM control',     p=scratch,s=f'data_{remap_str}_tem')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01', n='E3SM L72 smoothed',p=scratch,s=f'data_{remap_str}_tem')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01', n='E3SM L80 refined', p=scratch,s=f'data_{remap_str}_tem')



lev = np.array([ 1.,  2.,  3.,  5.,  7., 10., 20., 30., 50., 70., 100.,125.,150.,])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# add_var('u'        )
# add_var('dudt'     )
# add_var('utendepfd') 
# add_var('utendvtem')
# add_var('utendwtem')
add_var('BUTGWSPEC')
# add_var('UTGWSPEC' )

# lat1,lat2 = -5,5
lat1,lat2 = -10,10

fig_file,fig_type = f'figs/FXX-TEM-budget-time-series-contour.{var[0]}','png'

recalculate = True

add_u       = True
use_z_coord = False
plot_diff   = False
print_stats = True

tmp_file_head = 'data/TEM-budget.zonal_mean'

var_x_case = False
use_common_label_bar = True
num_plot_col = 1

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)

# if 'lev' not in vars(): lev = np.array([0])

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*(num_var*num_case)
res = hs.res_contour_fill()
res.vpHeightF = 0.1
res.tmYLLabelFontHeightF         = 0.01
res.tmXBLabelFontHeightF         = 0.01
res.tiXAxisFontHeightF           = 0.01
res.tiYAxisFontHeightF           = 0.01
# res.lbLabelBarOn = False

# disable these by default - turn back on for bottom panel(s)
res.tmXBOn = False
res.tiXAxisOn = False

res.tiXAxisString = 'Time'
res.tiXAxisString = 'Time [years]'

if use_common_label_bar: res.lbLabelBarOn = False

if use_z_coord:
   res.tiYAxisString = 'Height [km]'
   res.nglYAxisType = 'LinearAxis'
   res.trYMinF = 20
else:
   res.tiYAxisString = 'Pressure [hPa]'
   res.trYReverse = True
   # res.trYLog = True # doesn't work due to irregular spacing :(
   # res.trYMinF = 5
   # res.trYMaxF = 100e2

   tm_vals = np.array([1,10,50,100,200])
   res.tmYLMode = 'Explicit'
   res.tmYLValues = tm_vals*1e2
   res.tmYLLabels = tm_vals

ures = hs.res_contour()
ures.cnLineThicknessF = 4.

def neg_dash_contours(res_in, clev):
   dsh = np.zeros((len(clev)),'i')
   for k in range(len(clev)):
      # if (clev[k] < 0.): dsh[k] = 6
      if (clev[k] < 0.): dsh[k] = 2
   res_in.cnLineDashPatterns    = dsh
   res_in.cnMonoLineDashPattern = False
   res_in.cnFillOn  = False
   res_in.cnLinesOn = True

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print(hc.tcolor.GREEN+'  var: '+var[v]+hc.tcolor.ENDC)
   tvar = var[v]
   area_name = 'area'
   #----------------------------------------------------------------------------
   # read the data
   #----------------------------------------------------------------------------
   data_list,time_list,lev_list = [],[],[]
   u_data_list = []
   for c in range(num_case):
      print(hc.tcolor.CYAN+'    case: '+case[c]+hc.tcolor.ENDC)

      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}'
      if 'lat1' in locals(): tmp_file = tmp_file + f'.lat1_{lat1}.lat2_{lat2}'
      if 'lon1' in locals(): tmp_file = tmp_file + f'.lon1_{lon1}.lon2_{lon2}'
      tmp_file = f'{tmp_file}.nc'
      
      if recalculate :

         data_dir_tmp,data_sub_tmp = None, None
         # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
         if case_dir[c] is not None: data_dir_tmp = case_dir[c]
         if case_sub[c] is not None: data_sub_tmp = case_sub[c]

         case_obj = he.Case( name=case[c], atm_comp='eam', time_freq=None,
                             data_dir=data_dir_tmp, data_sub=data_sub_tmp )

         # avoid creating large chunks
         with dask.config.set(**{'array.slicing.split_large_chunks': True}):  

            file_path = f'{scratch}/{case[c]}/data_{remap_str}_tem/*.eam.{search_str}*'
            file_list = sorted(glob.glob(file_path))

            if 'first_file' in locals(): file_list = file_list[first_file:]
            if 'num_files'  in locals(): file_list = file_list[:num_files]

            ds = xr.open_mfdataset(file_list)

            data = ds[tvar]
            if add_u: u_data = ds['u'].load()
            if use_z_coord: Z = ds['?']

            # hc.print_time_length(data.time,indent=' '*4)

            data = data.sel(plev=slice(100e2,1e2))
            if add_u: u_data = u_data.sel(plev=slice(100e2,1e2))

            unit_str = ''
            if var[v]=='dudt'        : data = data*86400.; unit_str = 'm/s/day'
            if var[v]=='RES'         : data = data*86400.; unit_str = 'm/s/day'
            if var[v]=='TEND'        : data = data*86400.; unit_str = 'm/s/day'
            if var[v]=='utendepfd'   : data = data*86400.; unit_str = 'm/s/day'
            if var[v]=='utendepfd_lp': data = data*86400.; unit_str = 'm/s/day'
            if var[v]=='utendvtem'   : data = data*86400.; unit_str = 'm/s/day'
            if var[v]=='utendwtem'   : data = data*86400.; unit_str = 'm/s/day'
            if var[v]=='epfy'        : data = data*86400.; unit_str = 'm3/s/day'
            if var[v]=='epfz'        : data = data*86400.; unit_str = 'm3/s/day'
            if var[v]=='depfydy'     : data = data*86400.; unit_str = ''
            if var[v]=='depfzdp'     : data = data*86400.; unit_str = ''
            if var[v]=='psitem'      : data = data*86400.; unit_str = 'kg/day'
            if var[v]=='vtem'        :                     unit_str = 'm/s'
            if var[v]=='wtem'        :                     unit_str = 'm/s'
            if var[v]=='BUTGWSPEC'   : data = data*86400.; unit_str = 'm/s/day'
            if var[v]=='UTGWSPEC'    : data = data*86400.; unit_str = 'm/s/day'
            if var[v]=='UTGWORO'     : data = data*86400.; unit_str = 'm/s/day'

            avg_dims= ('lat')
            data_avg = data.sel(lat=slice(lat1,lat2)).mean(dim=avg_dims)
            if add_u: u_data = u_data.sel(lat=slice(lat1,lat2)).mean(dim=avg_dims)

            ### print stats after time averaging
            if print_stats: hc.print_stat(data_avg,name=var[v],stat='naxsh',indent='    ',compact=True)

            # if use_z_coord: Z_avg = ( (Z*area).sum(dim='ncol') / area.sum(dim='ncol') ).mean(dim='time') / 1e3

            yr = data_avg['time.year'].values
            mn = data_avg['time.month'].values
            time = yr + mn/12.

         #    #----------------------------------------------------------------------
         #    # write to temporary file
         #    tmp_ds = xr.Dataset()
         #    tmp_ds['lat']  = lat
         #    tmp_ds['plev'] = plev
         #    # tmp_ds['time'] = time
         #    tmp_ds['data'] = data
         #    # print(f'      writing to file: {tmp_file}')
         #    tmp_ds.to_netcdf(path=tmp_file,mode='w')
         # else:
         #    tmp_ds = xr.open_dataset( tmp_file )
         #    lat  = tmp_ds['lat']
         #    plev = tmp_ds['plev']
         #    # time = tmp_ds['time']
         #    data = tmp_ds['data']

         data_list.append( data_avg.transpose().values )
         if add_u: u_data_list.append(u_data.transpose().values)
         # time_list.append( data_avg['time'].values )
         time_list.append( time )
         if use_z_coord:
            lev_list.append( Z_avg.values )
         else:
            # lev_list.append( data_avg['lev'].values )
            lev_list.append( data_avg['plev'].values )

   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)

   # tres.cnFillPalette = np.array( cmocean.cm.diff(np.linspace(0,1,256)) )
   # tres.cnFillPalette = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )
   tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )

   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])

   tres.cnLevelSelectionMode = 'ExplicitLevels'

   tres.cnLevels = np.linspace(-50,50,num=21)/1e2

   # nlev = 21
   # aboutZero = False
   # if var[v] in ['U','V'] : aboutZero = True
   # clev_tup = ngl.nice_cntr_levels(data_min, data_max, cint=None, max_steps=nlev, returnLevels=False, aboutZero=aboutZero )
   # if clev_tup==None: 
   #    tres.cnLevelSelectionMode = 'AutomaticLevels'   
   # else:
   #    cmin,cmax,cint = clev_tup
   #    tres.cnLevels = np.linspace(cmin,cmax,num=nlev)
   #    tres.cnLevelSelectionMode = 'ExplicitLevels'

   for c in range(num_case):
      
      time = time_list[c]
      # time = ( time - time[0] ).astype('float') / 86400e9
      # time = ( time - time[0] ).astype('float') / (60*60*24*365)

      tres.sfYArray = lev_list[c]#.astype(int)
      tres.sfXArray = time
      # tres.sfXArray = np.linspace( 1./12., float(num_files)/12., num=len(time) )

      if c==num_case-1: 
         tres.tmXBOn = True
         tres.tiXAxisOn = True

      ip = v*num_case+c if var_x_case else c*num_var+v

      plot[ip] = ngl.contour(wks, data_list[c], tres)


      if add_u:
         tmp_ures = copy.deepcopy(ures)
         tmp_ures.sfXArray = tres.sfXArray
         tmp_ures.sfYArray = tres.sfYArray
         tmp_ures.cnLevelSelectionMode = "ExplicitLevels"
         tmp_ures.cnLevels = np.arange(-60,60+10,10)
         # if not hasattr(tmp_ures,'cnLevels'): 
         #    nlev_u = 11
         #    u_data_min = np.min(np.ma.masked_invalid(u_data_list[c]))
         #    u_data_max = np.max(np.ma.masked_invalid(u_data_list[c]))
         #    (cmin,cmax,cint) = ngl.nice_cntr_levels(u_data_min, u_data_max,
         #                       cint=None, max_steps=nlev_u,
         #                       returnLevels=False, aboutZero=False )
         #    tmp_ures.cnLevels = np.linspace(cmin,cmax,num=nlev_u)
         #    print(f'contour levels U - min: {cmin}  max: {cmax}  nlev: {nlev_u}')
         neg_dash_contours(tmp_ures, tmp_ures.cnLevels)
         ngl.overlay(plot[ip], ngl.contour(wks, np.ma.masked_invalid(u_data_list[c]), tmp_ures) )

      #-------------------------------------------------------------------------
      # Set strings at top of plot
      #-------------------------------------------------------------------------
      var_str = var[v]
      if var[v]=='PRECT' : var_str = 'Precipitation [mm/day]'
      if var[v]=='U'     : var_str = 'Zonal Wind [m/s]'

      ctr_str = ''
      # if var[v] in ['PRECT','PRECC','PRECL'] : ctr_str = 'Mean: '+'%.2f'%avg_X+' [mm/day]'
      hs.set_subtitles(wks, plot[ip], case_name[c], ctr_str, var_str, font_height=0.015)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

# layout = [len(plot),1]
# layout = [num_var,num_case] if var_x_case else [num_case,num_var]


layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent = 5
if use_common_label_bar: 
   pnl_res.nglPanelLabelBar   = True
   # pnl_res.lbTopMarginF       =  0.2
   # pnl_res.lbBottomMarginF    = -0.2
   pnl_res.lbLeftMarginF      = 0.5+0.5
   pnl_res.lbRightMarginF     = 0.5
   # pnl_res.lbTitleString      = "m/s"
   # pnl_res.lbTitlePosition    = "bottom"
   # pnl_res.lbLabelFontHeightF = 0.001
   # pnl_res.lbTitleFontHeightF = 0.01

### add panel labels
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.01

# if num_var==1  : layout = [num_case,num_var]
# if num_case==1 : layout = [num_var,num_case]

ngl.panel(wks,plot[0:len(plot)],layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
