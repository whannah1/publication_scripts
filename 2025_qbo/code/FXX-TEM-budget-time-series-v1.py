import os, glob, ngl, copy, xarray as xr, numpy as np
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
var,lev_list,vclr,vdsh = [],[],[],[]
vstr_list = []
def add_var(var_name,lev=-1,c='black',d=0,vstr=''): 
   var.append(var_name); lev_list.append(lev)
   vclr.append(c); vdsh.append(d)
   vstr_list.append(vstr)
##------------------------------------------------------------------------------
scratch = '/global/cfs/cdirs/m4310/whannah/E3SM'
# scratch = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu'
add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L72' ,n='E3SM L72', p=scratch,s='archive/atm/hist')
add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L80' ,n='E3SM L80', p=scratch,s='archive/atm/hist')
add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L128',n='E3SM L128',p=scratch,s='archive/atm/hist')

remap_str,search_str = 'remap_90x180','h0.tem.'; first_file,num_files = 12*0,12*20

#-------------------------------------------------------------------------------
# add_var('dudt'     ,c='black',d=1)
# add_var('RES'      ,c='red')
# add_var('TEND'     ,c='blue')

# add_var('dudt'     ,c='gray',d=0)
add_var('utendepfd',vstr='EP-Flux Div',c='blue') # u-wind tendency due to TEM EP flux divergence
# add_var('utendvtem',c='orange')
# add_var('utendwtem',c='pink')
add_var('BUTGWSPEC',vstr='Convective GWD',c='magenta')
# add_var('UTGWSPEC' ,c='blue')


# add_var('UTGWORO'  ,c='purple') # this is basically zero?

#-------------------------------------------------------------------------------

fig_file,fig_type = 'figs/FXX-TEM-time-series-v1','png'

recalculate = True

add_u = True
plot_diff = False

tmp_file_head = 'data/TEM-budget.zonal_mean'

print_stats = True

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_var = len(case),len(var)

if 'plev' not in locals(): plev = np.array([0])

wkres = ngl.Resources()
npix=2048; wkres.wkWidth,wkres.wkHeight=npix,npix
# wks = ngl.open_wks(fig_type,fig_file,wkres)

# plot = [None]*(num_case)
res = hs.res_xy()
res.vpHeightF = 0.1
res.nglMaximize  = False
res.tmYROn       = False
res.tmYRMinorOn  = False
res.tmYLLabelFontHeightF   = 0.01
res.tmXBLabelFontHeightF   = 0.01
res.tiXAxisFontHeightF     = 0.01
res.tiYAxisFontHeightF     = 0.01
# res.tiYAxisString          = 'Pressure [hPa]'
# res.tiXAxisString          = 'Latitude'

res.trYMinF,res.trYMaxF =-0.4,0.4

#-------------------------------------------------------------------------------
# create legend in separate file
#-------------------------------------------------------------------------------
legend_file = fig_file+'.legend'
wkres = ngl.Resources()
npix = 2048 ; wkres.wkWidth,wkres.wkHeight=npix,npix
lgd_wks = ngl.open_wks(fig_type,legend_file,wkres)

lgres = ngl.Resources()
lgres.pmTickMarkDisplayMode = 'NoCreate'
lgd_plot = ngl.blank_plot(lgd_wks, lgres)
gsres                  = ngl.Resources()
gsres.gsLineThicknessF = 30
txres                  = ngl.Resources()
txres.txFontHeightF    = 0.01
txres.txJust           = 'CenterLeft'
lgd_var = copy.deepcopy(vstr_list)
lgd_clr = copy.deepcopy(vclr)
if add_u: lgd_var.insert(0,'U'); lgd_clr.insert(0,'black')
num_lgd = len(lgd_var)
anno_id,tx_id,pgon_id = [None]*num_lgd,[None]*num_lgd,[None]*num_lgd
for v in range(num_lgd):
   txy = 0.3
   txx = 0.5*(v+1)/num_var
   tx_id[v] = ngl.text_ndc(lgd_wks, lgd_var[v], txx, txy, txres)

   gsres.gsLineColor = lgd_clr[v]
   px = txx - 0.03; pw = 0.01
   plx = np.array([px-pw,px+pw])
   ply = np.array([txy,txy])
   pgon_id[v] = ngl.polyline_ndc(lgd_wks, plx, ply, gsres)

ngl.draw(lgd_plot)
ngl.frame(lgd_wks)
hc.trim_png(legend_file)

# exit()

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# data_list_list,lev_list_list = [],[]
data_list = []
u_data_list = []
lev_list  = []
lat_list  = []
data_baseline_list = []
panel_list = []

for c in range(num_case):
   if 'plot'    in locals(): del plot
   if 'uplot'   in locals(): del uplot
   if 'wks_tmp' in locals(): del wks_tmp

   fig_file_tmp = f'{fig_file}.{case[c]}'
   wks_tmp = ngl.open_wks(fig_type,fig_file_tmp,wkres)

   print('  '+hc.tcolor.CYAN+'case: '+case[c]+hc.tcolor.ENDC)
   
   plot = None
   for v in range(num_var):
      print('    '+hc.tcolor.GREEN+'var: '+var[v]+hc.tcolor.ENDC)

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
      
      if recalculate :

         file_path = f'{case_dir[c]}/{case[c]}/data_{remap_str}_tem/*.eam.{search_str}*'
         file_list = sorted(glob.glob(file_path))

         if 'first_file' in locals(): file_list = file_list[first_file:]
         if 'num_files'  in locals(): file_list = file_list[:num_files]

         # for ff in file_list: print(ff)

         # ds = xr.open_mfdataset(file_list)
         #----------------------------------------------------------------------
         # manually combine files due to issues with open_mfdataset on Perlmutter
         ds = xr.open_dataset(file_list[0])
         for f in file_list[1:]:
           ds_tmp = xr.open_dataset(f)
           ds = xr.concat([ds, ds_tmp], dim='time')
         #----------------------------------------------------------------------

         hc.print_time_length(ds.time,indent=' '*6)

         # time = ( time - time[0] ).astype('float').values / 86400e9
         time = np.arange(0,len(ds.time),1) / 12

         tvar = var[v]
         if var[v]=='dudt': tvar = 'u'
         if var[v]=='RES':  tvar = 'u'
         if var[v]=='TEND': tvar = 'utendepfd'

         # Load data 
         data = ds[tvar].load()
         data = data.sel(lat=slice(-5,5)).mean(dim=('lat'))
         data = data.sel(plev=20e2)

         if var[v]=='dudt': data = data.differentiate('time')

         if var[v]=='RES': 
            dudt = data.differentiate('time')
            tend = ds['utendepfd'].sel(lat=slice(-5,5)).mean(dim=('lat')).sel(plev=20e2) \
                  +ds['utendvtem'].sel(lat=slice(-5,5)).mean(dim=('lat')).sel(plev=20e2) \
                  +ds['utendwtem'].sel(lat=slice(-5,5)).mean(dim=('lat')).sel(plev=20e2) \
                  +ds['BUTGWSPEC'].sel(lat=slice(-5,5)).mean(dim=('lat')).sel(plev=20e2) \
                  +ds['UTGWSPEC' ].sel(lat=slice(-5,5)).mean(dim=('lat')).sel(plev=20e2) \
                  +ds['UTGWORO'  ].sel(lat=slice(-5,5)).mean(dim=('lat')).sel(plev=20e2)
            data = dudt - tend

         if var[v]=='TEND': 
            data = ds['utendepfd'].sel(lat=slice(-5,5)).mean(dim=('lat')).sel(plev=20e2) \
                  +ds['utendvtem'].sel(lat=slice(-5,5)).mean(dim=('lat')).sel(plev=20e2) \
                  +ds['utendwtem'].sel(lat=slice(-5,5)).mean(dim=('lat')).sel(plev=20e2) \
                  +ds['BUTGWSPEC'].sel(lat=slice(-5,5)).mean(dim=('lat')).sel(plev=20e2) \
                  +ds['UTGWSPEC' ].sel(lat=slice(-5,5)).mean(dim=('lat')).sel(plev=20e2) \
                  +ds['UTGWORO'  ].sel(lat=slice(-5,5)).mean(dim=('lat')).sel(plev=20e2)

         if add_u and v==num_var-1:
            u_data = ds['u'].load()
            u_data = u_data.sel(lat=slice(-5,5)).mean(dim=('lat'))
            u_data = u_data.sel(plev=20e2)
            u_data_list.append(u_data.values)

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

      #-------------------------------------------------------------------------
      # Unit conversions

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

      #-------------------------------------------------------------------------
      # if 'time' in locals(): hc.print_time_length(time,indent=' '*6)

      if print_stats:
         msg = hc.print_stat(data,name=var[v],stat='naxsh',indent=(' '*6),compact=True,fmt='f')
         
      #-------------------------------------------------------------------------
      # Create plot
      #-------------------------------------------------------------------------
   
      tres = copy.deepcopy(res)

      # tres.tmYLMode = "Explicit"
      # tres.tmYLValues = plev[::4].values
      # tres.tmYLLabels = plev[::4].values / 1e2

      # tres.sfXArray = lat.values
      # tres.sfYArray = plev.values

      # tres.lbTitleString   = f'[{unit_str}]'

      # if plot_diff:
      #    if c==0:
      #       data_baseline_list.append(data)
      #    else:
      #       data = data - data_baseline_list[v]

      tres.xyLineColor   = vclr[v]
      tres.xyDashPattern = vdsh[v]
      
      # if var[v]=='u':
      #    tres.tiYAxisString      = "U component of wind (m/s)"
      #    tres.tiYAxisSide        = "Right"
      #    tres.tiYAxisFontColor   = "purple"
      #    # res2.tiXAxisFontHeightF = ngl.get_float(plot1,"tiXAxisFontHeightF")
      #    uplot = ngl.xy(wks, time, data.values, tres)
      # else:

      # tplot = ngl.xy(wks, time, np.ma.masked_invalid(data.values), tres)
      tplot = ngl.xy(wks_tmp, time, data.values, tres)

      if plot is None:
         plot = tplot
      else:
         ngl.overlay(plot,tplot)

      if v==(num_var-1):
         cstr = ''
         # if plot_diff and c>0: cstr = f'diff'
         hs.set_subtitles(wks_tmp, plot, name[c], '', '', center_sub_string=cstr, font_height=0.01)

   # add reference line
   lres = hs.res_xy()
   lres.xyLineThicknessF = 4
   lres.xyDashPattern = 0
   lres.xyLineColor = 'gray'
   ngl.overlay(plot,ngl.xy(wks_tmp, np.array([-1e3,1e8]), np.array([0,0]), lres))

   if add_u:
      ngl.draw(plot)
      ures = copy.deepcopy(res)
      ures.vpHeightF        = ngl.get_float(plot,'vpHeightF')
      ures.vpWidthF         = ngl.get_float(plot,'vpWidthF')
      ures.vpXF             = ngl.get_float(plot,'vpXF')
      ures.vpYF             = ngl.get_float(plot,'vpYF')
      ures.tmYUseLeft       = False
      ures.tmXBOn           = False
      ures.tmXBLabelsOn     = False
      ures.tmXBMinorOn      = False
      ures.tmYLOn           = False
      ures.tmYLLabelsOn     = False
      ures.tmYLMinorOn      = False
      ures.tmYRLabelsOn     = True
      ures.tmYROn           = True
      ures.tiYAxisString    = 'Zonal Wind'
      ures.tiYAxisSide      = 'Right'
      ures.xyLineThicknessF = res.xyLineThicknessF*2
      ures.trYMinF          = -35
      ures.trYMaxF          = 35
      
      uplot = ngl.xy(wks_tmp, time, u_data.values, ures)   
      ngl.draw(uplot)
   #    ngl.draw(plot)
   # else:
   #    ngl.draw(plot)

   ngl.draw(plot)
   #------------------------------------------------------------------------------------------------
   # finalize the individual panel into separate file
   #------------------------------------------------------------------------------------------------
   ngl.frame(wks_tmp)
   hc.trim_png(fig_file_tmp)

   panel_list.append(f'{fig_file_tmp}.{fig_type}')

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
ngl.end()

if fig_type=='png':
   cmd  = f'montage'
   cmd += f' -density 300'
   cmd += f' -tile 1x{(len(panel_list)+1)}'
   cmd += f' -geometry +20+20'
   # cmd += f' -composite -gravity south '
   # cmd += f' -gravity south '
   cmd += f' {" ".join(panel_list)}'
   cmd += f' {legend_file}.{fig_type}'
   cmd += f' {fig_file}.{fig_type}'

   print(hc.tcolor.CYAN+cmd+hc.tcolor.ENDC+'\n')
   os.system(cmd)

   hc.trim_png(fig_file)


#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
