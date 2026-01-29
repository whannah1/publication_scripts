import os, ngl, copy, xarray as xr, numpy as np
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import cmocean
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
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01',       n='E3SM control',     c='red')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01', n='E3SM L72 smoothed',c='green')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01', n='E3SM L80 refined', c='blue')
#-------------------------------------------------------------------------------

### define pressure levels for plotting - initial calculations stay on native grid
# lev = np.array([30,50,75,100,125,150,200,250,300,350,400,450,500,550,600,650,700,750,800,825,850,875,900,925,950,975,1000])
# lev = np.array([50,100,150,200,300,400,500,600,700,750,800,850,875,900,925,975])
plev = np.array([0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200,250,300,400,500,600,700,800,850,900])

var = ['BTAUXS']

tvar_list,tau_list = [],[]
for n in reversed(range(32)): tvar_list.append(f'BTAUXSn{n+1:02d}'); tau_list.append((n+1)*-1)
for n in range(0,32+1):       tvar_list.append(f'BTAUXSp{n:02d}');   tau_list.append(n)

num_plot_col = len(case)

fig_file,fig_type = 'figs/FXX-BTAUX-ratio','png'


# lat1,lat2 = -10,10
lat1,lat2 = -30,30
# lat1,lat2,lon1,lon2 = 15,25,360-60,360-50

# xlat,xlon,dy,dx = 60,120,10,10;
# if 'xlat' in locals(): lat1,lat2,lon1,lon2 = xlat-dy,xlat+dy,xlon-dx,xlon+dx


# htype,first_file,num_files = 'h0',0,2
htype,first_file,num_files = 'h0',0,20*12
# htype,first_file,num_files = 'h1',0,1
# htype,first_file,num_files = 'h2',0,0

# recalculate = False

plot_diff = False

common_colorbar = False

# Y-axis options
use_height_coord = False
# omit_bot,bot_k   = False,-2
# omit_top,top_k   = False,30
print_stats      = False
print_profile    = False

tmp_file_head = 'data/BTAUX-spectrum'

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_var = len(case),len(var)

if 'diff_mode' not in locals(): diff_mode = 0
if 'plev' not in locals(): plev = np.array([0])

wks = ngl.open_wks(fig_type,fig_file)
plot = [None]*(num_case)
res = hs.res_contour_fill()
# res.vpWidthF = 0.4
# res.xyMarkLineMode = "MarkLines"
# res.xyMarkerSizeF = 0.008
# res.xyMarker = 16
# res.xyLineThicknessF = 8
# res.tmYLLabelFontHeightF         = 0.008
# res.tmXBLabelFontHeightF         = 0.008
# res.tmXBAutoPrecision = False
# res.tmXBPrecision = 2

res.tiXAxisFontHeightF  = 0.02
res.tiYAxisFontHeightF  = 0.02
# res.tiYAxisString       = 'Pressure [hPa]'
res.tiXAxisString       = ''

if common_colorbar:
   res.lbLabelBarOn = False
else:
   res.lbLabelBarOn = True
   res.lbTitleFontHeightF = 0.015
   res.lbLabelFontHeightF = 0.015

if not use_height_coord: 
   res.trYReverse = True
   # res.xyYStyle = 'Log'

def get_comp(case):
   comp = 'eam'
   if 'CESM' in case: comp = 'cam'
   return comp
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# data_list_list,lev_list_list = [],[]
for v in range(num_var):
   hc.printline()
   print(hc.tcolor.GREEN+'  var: '+var[v]+hc.tcolor.ENDC)
   data_list,lev_list,bin_list = [],[],[]
   for c in range(num_case):
      print(hc.tcolor.CYAN+'    case: '+case[c]+hc.tcolor.ENDC)

      data_dir_tmp,data_sub_tmp = None, None
      # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
      if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      if case_sub[c] is not None: data_sub_tmp = case_sub[c]
      
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}'
      if 'lat1' in locals(): tmp_file = tmp_file + f'.lat1_{lat1}.lat2_{lat2}'
      if 'lon1' in locals(): tmp_file = tmp_file + f'.lon1_{lon1}.lon2_{lon2}'
      tmp_file = f'{tmp_file}.nc'
      print(f'      tmp_file: {tmp_file}')

      # if recalculate :

      #    case_obj = he.Case( name=case[c], atm_comp=get_comp(case[c]), data_dir=data_dir_tmp, data_sub=data_sub_tmp  )

      #    if 'lat1' in locals(): case_obj.lat1 = lat1; case_obj.lat2 = lat2
      #    if 'lon1' in locals(): case_obj.lon1 = lon1; case_obj.lon2 = lon2

      #    #-------------------------------------------------------------------------
      #    # read the height and area data
      #    area = case_obj.load_data('area',htype=htype,first_file=first_file,num_files=1,).astype(np.double)   
      #    Z = case_obj.load_data('Z3',lev=plev,htype=htype,first_file=first_file,num_files=num_files)
         
      #    ncol = Z.ncol

      #    if 'time' in Z.dims : 
      #       time = Z.time
      #       hc.print_time_length(Z.time,indent=' '*6)

      #    lev = Z.lev

      #    # Time mean
      #    Z = Z.mean(dim='time')

      #    # area weighted spatial mean
      #    Z = ( (Z   *area).sum(dim='ncol') / area.sum(dim='ncol') )

      #    #-------------------------------------------------------------------------
      #    # read the spectral GW data
      #    for t,tvar in enumerate(tvar_list):
      #       data_tmp = case_obj.load_data(tvar,lev=plev,htype=htype,first_file=first_file,num_files=num_files)
      #       data_tmp = data_tmp.mean(dim='time')
      #       data_tmp = ( (data_tmp*area).sum(dim='ncol') / area.sum(dim='ncol') )

      #       if t==0:
      #          shape  = (len(lev),len(tau_list))
      #          coords = [('lev', lev),('bins', np.array(tau_list))]
      #          data = xr.DataArray(np.zeros(shape), coords=coords )

      #       data[:,t] = data_tmp.values

      #    #----------------------------------------------------------------------
      #    # write to temporary file
      #    tmp_ds = xr.Dataset()
      #    tmp_ds['data'] = data
      #    tmp_ds['Z']    = Z
      #    if 'time' in locals(): tmp_ds['time'] = time
      #    if 'lat1' in locals(): tmp_ds['lat1'] = lat1; tmp_ds['lat2']=lat2
      #    if 'lon1' in locals(): tmp_ds['lon1'] = lon1; tmp_ds['lon2']=lon2
      #    print(f'      writing to file: {tmp_file}')
      #    tmp_ds.to_netcdf(path=tmp_file,mode='w')
      # else:
      if True:
         tmp_ds = xr.open_dataset( tmp_file )
         data = tmp_ds['data']
         bins = tmp_ds['bins']
         time = tmp_ds['time']
         Z    = tmp_ds['Z']
      
      #-------------------------------------------------------------------------
      # # options for omitting top or bottom levels
      # if omit_bot: data=data[     :bot_k]; Z=Z[     :bot_k]
      # if omit_top: data=data[top_k:     ]; Z=Z[top_k:     ]
      #-------------------------------------------------------------------------
      # print stuff
      if print_stats: hc.print_stat(data,name=f'{var[v]} after averaging',indent='    ',compact=True)
      if print_profile: 
         print()
         for xx in data.values: print(f'    {xx}')
         print()
      #-------------------------------------------------------------------------
      # calculate p/n ratio at each level (omit zero phase speed?)
      nbin_old = len(bins)
      nbin_new = int((nbin_old-1)/2)
      ratio_data = xr.full_like(data[:,:nbin_new],-999)
      ratio_bins = xr.full_like(bins[:nbin_new],-999)
      for b in range(nbin_new):
         b1 = b
         b2 = nbin_old-1-b
         b3 = nbin_new-1-b
         # print(f'  b: {b}  b1: {bins[b1].values}  b2: {bins[b2].values}')
         ratio_data[:,b3] = data[:,b2] / data[:,b1]
         ratio_bins[b3]   = bins[b2]
         # print(f'  b: {b}  b1: {bins[b1].values}  b2: {bins[b2].values}  b3: {ratio_bins[b3].values}')
      data = ratio_data
      bins = ratio_bins
      # print(); print(data)
      # print(); print(bins)
      # exit()
      #-------------------------------------------------------------------------
      # append final data to list
      data_list.append( data.values )
      bin_list.append( bins.values )
      if     use_height_coord: lev_list.append( Z.values )
      if not use_height_coord: lev_list.append( data['lev'].values )

   # data_list_list.append(data_list)
   # lev_list_list.append(lev_list)

#-------------------------------------------------------------------------------
# Create plot - overlay all cases for each var
#-------------------------------------------------------------------------------
   
   tres = copy.deepcopy(res)

   if plot_diff:
      for c in range(1,num_case): 
         data_list[c] = data_list[c] - data_list[0]

   diff_mag_max = np.max(np.absolute(data_list[1:]))

   tres.tmYLMode = "Explicit"
   tres.tmYLValues = plev[::2]
   tres.tmYLLabels = plev[::2]

   # data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])
   # tres.trXMinF = data_min
   # tres.trXMaxF = data_max
   # ip = v

   tres.trXMaxF = 10

   # data_list_interp = []
   # data_tmp = ngl.vinth2p( data_list[c], hya, hyb, lev, PS_dum.values, 
   #                            interp_type, P0, 1, extrap_flag)


   for c in range(num_case):

      if (plot_diff and c==0) or not plot_diff:
         tres.cnLevelSelectionMode = 'ExplicitLevels'
         # np.linspace(0.,data_max,21)
         tres.cnLevels = np.arange(4,36+4,4)*1e-1
         # tres.cnLevels = np.arange(-20,20+2,2)*1e-5
         tres.cnFillPalette = "MPL_viridis"
         # tres.cnFillPalette = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )
      
      if plot_diff and c>0:
         tres.cnLevelSelectionMode = 'ExplicitLevels'
         tres.cnLevels = np.linspace(-1*diff_mag_max,diff_mag_max,21)
         tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )

      tres.sfXArray = bin_list[c]
      tres.sfYArray = lev_list[c]

      # tres.sfXCStartV = min( bin_list[c] )
      # tres.sfXCEndV   = max( bin_list[c] )
      # tres.sfYCStartV = min( lev_list[c] )
      # tres.sfYCEndV   = max( lev_list[c] )

      # ip = v*num_case+c
      ip = c*num_var+v
      
      plot[ip] = ngl.contour(wks, data_list[c], tres)

      rstr = ''
      if plot_diff and c>0: rstr = 'diff'

      hs.set_subtitles(wks, plot[ip], name[c], '', rstr, font_height=0.01)

   ### add vertical line
   # lres = hs.res_xy()
   # lres.xyLineThicknessF = 1
   # lres.xyDashPattern = 0
   # lres.xyLineColor = 'black'
   # ngl.overlay(plot[ip],ngl.xy(wks, np.array([0,0]), np.array([-1e3,1e8]), lres))

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
