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
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01',       n='E3SM control',     c='red'  ,p=gscratch,s='run')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01', n='E3SM L72 smoothed',c='green',p=gscratch,s='run')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01', n='E3SM L80 refined', c='blue' ,p=gscratch,s='run')

# add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72.GW-MOD-0', n='E3SM L72'   ,d=1,c='red',   p=pscratch,s='run')
# add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-nsu40',    n='E3SM L72-nsu40',d=1,c='red',   p=pscratch,s='run')
#add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-rscl',     n='E3SM L72-rscl', d=1,c='purple',p=pscratch,s='run')
#add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-rlim',     n='E3SM L72-rlim', d=1,c='pink',  p=pscratch,s='run')


# scratch = '/global/cfs/cdirs/m4310/whannah/E3SM'
# add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L72' ,n='E3SM L72', d=0,c='red',     p=scratch,s='archive/atm/hist')
# add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L80' ,n='E3SM L80', d=0,c='blue',    p=scratch,s='archive/atm/hist')
# add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L128',n='E3SM L128',d=0,c='purple',  p=scratch,s='run')

# scratch = '/global/cfs/cdirs/m4310/wandiyu/E3SM/'
# add_case('E3SM.2023-SCIDAC.ne30pg2_EC30to60E2r2.F20TR-MMF1.L64',n='MMF L64',d=1,c='red',     p=scratch,s='archive/atm/hist')
# add_case('E3SM.2023-SCIDAC.ne30pg2_EC30to60E2r2.F20TR-MMF1.L72',n='MMF L72',d=1,c='blue',    p=scratch,s='archive/atm/hist')

### 2024 MMF GWD test
pscratch,psub = '/pscratch/sd/w/whannah/e3sm_scratch/pm-gpu','run'
# pscratch,psub = '/global/cfs/cdirs/m4310/whannah/E3SM',''
add_case('E3SM.2024-MMF-GW-test.ne30pg2_EC30to60E2r2.F20TR-MMF1.L72.CTL',n='MMF CTL',p=pscratch,s=psub,c='red')
add_case('E3SM.2024-MMF-GW-test.ne30pg2_EC30to60E2r2.F20TR-MMF1.L72.EXP',n='MMF EXP',p=pscratch,s=psub,c='blue')

#-------------------------------------------------------------------------------

### define pressure levels for plotting - initial calculations stay on native grid
# lev = np.array([30,50,75,100,125,150,200,250,300,350,400,450,500,550,600,650,700,750,800,825,850,875,900,925,950,975,1000])
# lev = np.array([50,100,150,200,300,400,500,600,700,750,800,850,875,900,925,975])
# plev = np.array([0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200,250,300,400,500,600,700,800,850,900])
plev = np.array([1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200,250,300,400,500,600,700,800,850,900])

var = ['BTAUXS']

tvar_list,tau_list = [],[]
for n in reversed(range(32)): tvar_list.append(f'BTAUXSn{n+1:02d}'); tau_list.append(2.5*(n+1)*-1)
for n in range(0,32+1):       tvar_list.append(f'BTAUXSp{n:02d}');   tau_list.append(2.5*n)

num_plot_col = len(case)

fig_file,fig_type = 'figs/FXX-BTAUX-spectrum','png'


lat1,lat2 = -5,5
# lat1,lat2 = -10,10
# lat1,lat2 = -30,30
# lat1,lat2,lon1,lon2 = 15,25,360-60,360-50

# xlat,xlon,dy,dx = 60,120,10,10;
# if 'xlat' in locals(): lat1,lat2,lon1,lon2 = xlat-dy,xlat+dy,xlon-dx,xlon+dx


htype,first_file,num_files = 'h0',0,1
# htype,first_file,num_files = 'h0',0,12*5
# htype,first_file,num_files = 'h0',0,20*12
# htype,first_file,num_files = 'h1',0,10
# htype,first_file,num_files = 'h2',0,0

recalculate = True

plot_diff = False

common_colorbar = False

# Y-axis options
use_height_coord = False
omit_bot,bot_k   = False,-2
omit_top,top_k   = False,30
print_stats      = False
print_profile    = False

if num_files==12:
   tmp_file_head = 'data/BTAUX-spectrum-1yr'
else:
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

res.trXMinF = -40
res.trXMaxF =  40

res.tiXAxisFontHeightF  = 0.02
res.tiYAxisFontHeightF  = 0.02
res.tiYAxisString       = 'Pressure [hPa]'
res.tiXAxisString       = 'Gravity Wave Phase Speed [m/s/day]'


ures = hs.res_xy()
ures.xyLineThicknessF = 4
ures.xyDashPattern = 0
ures.tiYAxisString = ''
ures.tiXAxisString = ''

# ures.tmXTOn = True
# ures.tmXTLabelsOn = True


if common_colorbar:
   res.lbLabelBarOn = False
else:
   res.lbLabelBarOn = True
   res.lbTitleString = 'Convective GW Tau [Pa]'
   res.lbTitlePosition = 'Bottom'
   res.lbTitleFontHeightF = 0.015
   res.lbLabelFontHeightF = 0.015

if not use_height_coord: 
   res.trYReverse = True
   ures.trYReverse = True
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
   wind_list = []
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

      if recalculate :

         case_obj = he.Case( name=case[c], atm_comp=get_comp(case[c]), data_dir=data_dir_tmp, data_sub=data_sub_tmp  )

         if 'lat1' in locals(): case_obj.lat1 = lat1; case_obj.lat2 = lat2
         if 'lon1' in locals(): case_obj.lon1 = lon1; case_obj.lon2 = lon2

         #-------------------------------------------------------------------------
         # read the height and area data
         area = case_obj.load_data('area',htype=htype,first_file=first_file,num_files=1,).astype(np.double)   
         Z = case_obj.load_data('Z3',lev=plev,htype=htype,first_file=first_file,num_files=num_files)
         U = case_obj.load_data('U', lev=plev,htype=htype,first_file=first_file,num_files=num_files)
         V = case_obj.load_data('V', lev=plev,htype=htype,first_file=first_file,num_files=num_files)
         
         ncol = Z.ncol

         if 'time' in Z.dims : 
            time = Z.time
            hc.print_time_length(Z.time,indent=' '*6)

         lev = Z.lev.values

         # Time mean
         Z = Z.mean(dim='time')
         U = U.mean(dim='time')
         V = V.mean(dim='time')


         # area weighted spatial mean
         Z = ( (Z*area).sum(dim='ncol') / area.sum(dim='ncol') )
         U = ( (U*area).sum(dim='ncol') / area.sum(dim='ncol') )
         V = ( (V*area).sum(dim='ncol') / area.sum(dim='ncol') )

         #-------------------------------------------------------------------------
         # read the spectral GW data
         for t,tvar in enumerate(tvar_list):
            # data_tmp = case_obj.load_data(tvar,lev=plev,htype=htype,first_file=first_file,num_files=num_files)
            print('WARNING - disabled pressure level interpolation')
            data_tmp = case_obj.load_data(tvar,htype=htype,first_file=first_file,num_files=num_files)

            for k in range(len(data_tmp.lev.values)):
            # if True:
            #    k = 11
               msg = f'{tvar}  t: {t:3d}  k: {k:3d}'
               msg+=f'  min: {data_tmp.isel(lev=k).min().values: 20.12f}'
               msg+=f'  avg: {data_tmp.isel(lev=k).mean().values:20.12f}'
               msg+=f'  max: {data_tmp.isel(lev=k).max().values: 20.12f}'
               print(msg)

            data_tmp = data_tmp.mean(dim='time')
            data_tmp = ( (data_tmp*area).sum(dim='ncol') / area.sum(dim='ncol') )

            # # for k in range(len(data_tmp.lev.values)):
            # if True:
            #    k = 11
            #    msg = f'{tvar}  t: {t:3d}  k: {k:3d}'
            #    msg+=f'  min: {data_tmp.isel(lev=k).min().values:12.2f}'
            #    msg+=f'  avg: {data_tmp.isel(lev=k).mean().values:12.2f}'
            #    msg+=f'  max: {data_tmp.isel(lev=k).max().values:12.2f}'
            #    print(msg)

            if t==0:
               shape  = (len(lev),len(tau_list))
               coords = [('lev', lev),('bins', np.array(tau_list))]
               data = xr.DataArray(np.zeros(shape), coords=coords )

            data[:,t] = data_tmp.values

         # exit('EXITING')
         #----------------------------------------------------------------------
         # write to temporary file
         tmp_ds = xr.Dataset()
         tmp_ds['data'] = data
         tmp_ds['Z']    = Z
         tmp_ds['U']    = U
         tmp_ds['V']    = V
         if 'time' in locals(): tmp_ds['time'] = time
         if 'lat1' in locals(): tmp_ds['lat1'] = lat1; tmp_ds['lat2']=lat2
         if 'lon1' in locals(): tmp_ds['lon1'] = lon1; tmp_ds['lon2']=lon2
         print(f'      writing to file: {tmp_file}')
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         tmp_ds = xr.open_dataset( tmp_file )
         data = tmp_ds['data']
         time = tmp_ds['time']
         Z    = tmp_ds['Z']
         U    = tmp_ds['U']
         V    = tmp_ds['V']

      data = data*86400. # convert from m/s/s to m/s/day

      # print();print(data)
      # print();print(data.values)
      # exit()

      #-------------------------------------------------------------------------
      # # try removing the mean at each level
      for k in range(len(data.lev.values)):
         # tmp = np.mean(data[k,:].values)
         # print(f'  k: {k}    {tmp}')
         # hc.print_stat(data[k,:],name=f'{var[v]} k: {k}',indent='    ',compact=True)
         msg = f'{var[v]} k: {k}'
         msg+=f'  min: {data[k,:].min().values:12.2f}'
         msg+=f'  avg: {data[k,:].mean().values:12.2f}'
         msg+=f'  max: {data[k,:].max().values:12.2f}'
         print(msg)
         # data[k,:] = data[k,:] - np.mean(data[k,:].values)
      exit()
      
      #-------------------------------------------------------------------------
      # options for omitting top or bottom levels
      if omit_bot: data=data[     :bot_k]; Z=Z[     :bot_k]; U=U[     :bot_k]; V=V[     :bot_k]
      if omit_top: data=data[top_k:     ]; Z=Z[top_k:     ]; U=U[top_k:     ]; V=V[top_k:     ]
      #-------------------------------------------------------------------------
      # print stuff
      if print_stats: hc.print_stat(data,name=f'{var[v]} after averaging',indent='    ',compact=True)
      if print_profile: 
         print()
         for xx in data.values: print(f'    {xx}')
         print()
      #-------------------------------------------------------------------------
      # append final data to list
      data_list.append( data.values )
      # bin_list.append(tmp_ds['bins'].values )
      bin_list.append( np.array(tau_list) )
      if     use_height_coord: lev_list.append( Z.values )
      if not use_height_coord: lev_list.append( data['lev'].values )

      wind_list.append(U.values)

   # data_list_list.append(data_list)
   # lev_list_list.append(lev_list)

#-------------------------------------------------------------------------------
# Create plot - overlay all cases for each var
#-------------------------------------------------------------------------------
   
   tres = copy.deepcopy(res)

   if plot_diff:
      for c in range(1,num_case): 
         data_list[c] = data_list[c] - data_list[0]

   if plot_diff and num_case>1: diff_mag_max = np.max(np.absolute(data_list[1:]))

   tres.tmYLMode = "Explicit"
   tres.tmYLValues = plev[::2]
   tres.tmYLLabels = plev[::2]

   # data_min = np.min([np.nanmin(d) for d in data_list])
   # data_max = np.max([np.nanmax(d) for d in data_list])
   # tres.trXMinF = data_min
   # tres.trXMaxF = data_max
   # ip = v

   # data_list_interp = []
   # data_tmp = ngl.vinth2p( data_list[c], hya, hyb, lev, PS_dum.values, 
   #                            interp_type, P0, 1, extrap_flag)


   for c in range(num_case):

      if (plot_diff and c==0) or not plot_diff:
         tres.cnLevelSelectionMode = 'ExplicitLevels'
         # tres.cnLevels = np.arange(-20,20+2,2)*1e-5 # use for units = m/s/s
         # tres.cnLevels = np.arange(-18,18+2,2) # use for units = m/s/day
         # tres.cnFillPalette = "MPL_viridis"
         tres.cnFillPalette = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )
      
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

      # Overlay wind profile
      for c in range(num_case):
         ures.xyLineColor = clr[c]
         ures.xyDashPattern = dsh[c]
         ngl.overlay(plot[ip],ngl.xy(wks, wind_list[c], lev_list[c], ures))

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
