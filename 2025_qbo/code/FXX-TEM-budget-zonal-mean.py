import os, glob, ngl, copy, xarray as xr, numpy as np
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import cmocean
# scratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl'
scratch = '/global/cfs/cdirs/m3312/whannah/2022-QBO-TEST'
#-------------------------------------------------------------------------------
case_dir,case_sub = [],[]
case,name = [],[]
clr,dsh,mrk = [],[],[]
def add_case(case_in,n=None,p=None,s=None,d=0,c='black',m=0):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); name.append(n)
   case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
var,vx_list,vy_list,cnt_list,lev_list = [],[],[],[],[]
def add_var(var_name,lev=-1,vx=None,vy=None,c=None): 
   var.append(var_name)
   lev_list.append(lev)
   cnt_list.append(c)
   vx_list.append(vx)
   vy_list.append(vy)
##------------------------------------------------------------------------------
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01',       n='E3SM')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01',       n='E3SM control',     c='red')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01', n='E3SM L72 smoothed',c='green')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01', n='E3SM L80 refined', c='blue')

remap_str,search_str = 'remap_90x180','h0.tem.'; num_files = 12*20
# remap_str,search_str = 'remap_90x180','h2.tem.'; num_files = 1

#-------------------------------------------------------------------------------
# add_var('utendepfd') # u-wind tendency due to TEM EP flux divergence
# add_var('utendvtem') # u-wind tendency due to TEM northward wind advection and coriolis
# add_var('utendwtem') # u-wind tendency due to TEM upward wind advection

# add_var('utendepfd',   c='u') # u-wind tendency due to TEM EP flux divergence
# add_var('utendepfd',   c='BUTGWSPEC') # u-wind tendency due to TEM EP flux divergence

add_var('utendepfd',c='u') # u-wind tendency due to TEM EP flux divergence
add_var('utendvtem',c='u') # u-wind tendency due to TEM northward wind advection and coriolis
add_var('utendwtem',c='u') # u-wind tendency due to TEM upward wind advection

# add_var('u',    c='u')
# add_var('vtem', c='u')
# add_var('wtem', c='u')

# add_var('BUTGWSPEC',c='u')
# add_var('UTGWSPEC', c='u')
# add_var('UTGWORO',  c='u')


# add_var('epfy')     # Northward component of the Eliassen-Palm flux
# add_var('epfz')     # Upward component of the Eliassen-Palm flux
# add_var('vtem')     # Transformed Eulerian mean northward wind
# add_var('wtem')     # Transformed Eulerian mean upward wind
# add_var('psitem')   # Transformed Eulerian mean mass stream function
# add_var('depfydy')  # Meridional derivative of northward component of the Eliassen-Palm flux
# add_var('depfzdz')  # Vertical derivative of upward component of the Eliassen-Palm flux

# add_var('utendepfd',vx='epfy',vy='epfz') # u-wind tendency due to TEM EP flux divergence

#-------------------------------------------------------------------------------

# fig_file,fig_type = 'figs/FXX-TEM-budget-zonal-mean','png'

if var[0]=='utendepfd': fig_file,fig_type = f'figs/TEM-budget-20yr-epflx','png'
if var[0]=='u'        : fig_file,fig_type = f'figs/TEM-budget-20yr-winds','png'
if var[0]=='BUTGWSPEC': fig_file,fig_type = f'figs/TEM-budget-20yr-gwdrg','png'

recalculate = True

plot_diff = True
# common_colorbar = False

tmp_file_head = 'data/TEM-budget.zonal_mean'

print_stats = False

num_plot_col = 3
var_x_case = True

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_var = len(case),len(var)

if 'plev' not in locals(): plev = np.array([0])

wkres = ngl.Resources()
npix=4096; wkres.wkWidth,wkres.wkHeight=npix,npix

wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*(num_case*num_var)
res = hs.res_contour_fill()
res.vpHeightF = 0.4
res.trYReverse = True
res.lbTitleFontHeightF     = 0.015
res.lbLabelFontHeightF     = 0.015
# res.tmYLLabelFontHeightF   = 0.01
# res.tmXBLabelFontHeightF   = 0.01
# res.tiXAxisFontHeightF     = 0.02
# res.tiYAxisFontHeightF     = 0.02
res.tiYAxisString          = 'Pressure [hPa]'
res.tiXAxisString          = 'Latitude'
res.lbTitleString   = ''
res.lbTitlePosition = 'Bottom'

# res.lbLabelBarOn = False if common_colorbar else True

res.trYMaxF = 500e2

vres = hs.res_default()

cres = hs.res_contour()
cres.cnLineThicknessF = 4.

def neg_dash_contours(res_in, clev):
   dsh = np.zeros((len(clev)),'i')
   for k in range(len(clev)):
      if (clev[k] < 0.): dsh[k] = 6
   res_in.cnLineDashPatterns    = dsh
   res_in.cnMonoLineDashPattern = False
   res_in.cnFillOn  = False
   res_in.cnLinesOn = True

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# data_list_list,lev_list_list = [],[]
data_list = []
lev_list  = []
lat_list  = []
data_vx_list = []
data_vy_list = []
data_baseline_list = []

for c in range(num_case):
   print(hc.tcolor.CYAN+'    case: '+case[c]+hc.tcolor.ENDC)

   # print(hc.tcolor.RED+'WARNING: alternate config for method comparison! '+case[c]+hc.tcolor.ENDC)
   # if c==0:
   #    if 'h0' in search_str: remap_str,search_str = 'remap_90x180','h0.tem.'
   #    if 'h2' in search_str: remap_str,search_str = 'remap_90x180','h2.tem.'
   # if c==1:
   #    if 'h0' in search_str: remap_str,search_str = 'remap_90x180','h0.tem_y.'
   #    if 'h2' in search_str: remap_str,search_str = 'remap_90x180','h2.tem_y.'

   for v in range(num_var):
      # print(hc.tcolor.GREEN+'  var: '+var[v]+hc.tcolor.ENDC)

      data_dir_tmp,data_sub_tmp = None, None
      # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
      if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      if case_sub[c] is not None: data_sub_tmp = case_sub[c]

      add_vectors = False
      if vx_list[v] is not None and vy_list[v] is not None: add_vectors = True
      
      add_contour = False
      if cnt_list[v] is not None: add_contour = True
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}'
      if 'lat1' in locals(): tmp_file = tmp_file + f'.lat1_{lat1}.lat2_{lat2}'
      if 'lon1' in locals(): tmp_file = tmp_file + f'.lon1_{lon1}.lon2_{lon2}'
      tmp_file = f'{tmp_file}.nc'
      
      # print(f'      tmp_file: {tmp_file}')

      if recalculate :

         # file_path = f'{scratch}/{case[c]}/data_{remap_str}_prs/*.eam.{search_str}*'
         file_path = f'{scratch}/{case[c]}/data_{remap_str}_tem/*.eam.{search_str}*'
         file_list = sorted(glob.glob(file_path))

         # exit(file_path)

         if 'num_files' in locals(): file_list = file_list[:num_files]

         # for ff in file_list: print(ff)

         ds = xr.open_mfdataset(file_list)

         lat  = ds['lat']
         plev = ds['plev']
         # time = ds['time']

         # Load data for contours/shading
         data = ds[var[v]].load()

         if 'time' in ds.dims: 
            hc.print_time_length(ds.time,indent=' '*6)
            data = data.mean(dim='time')

         # Load data for vectors
         if add_vectors: 
            data_vx = ds[vx_list[v]].load()#.mean(dim='time')
            data_vy = ds[vy_list[v]].load()#.mean(dim='time')
            if 'time' in ds.dims: 
               data_vx = data_vx.mean(dim='time')
               data_vy = data_vy.mean(dim='time')

         if add_contour: 
            data_cnt = ds[cnt_list[v]].load()#.mean(dim='time')
            if 'time' in ds.dims: 
               data_cnt = data_cnt.mean(dim='time')

         #----------------------------------------------------------------------
         # write to temporary file
         tmp_ds = xr.Dataset()
         tmp_ds['lat']  = lat
         tmp_ds['plev'] = plev
         # tmp_ds['time'] = time
         tmp_ds['data'] = data
         if add_vectors:
            tmp_ds['data_vx'] = data_vx
            tmp_ds['data_vy'] = data_vy
         if add_contour: 
            tmp_ds['data_cnt'] = data_cnt
         # print(f'      writing to file: {tmp_file}')
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         tmp_ds = xr.open_dataset( tmp_file )
         lat  = tmp_ds['lat']
         plev = tmp_ds['plev']
         # time = tmp_ds['time']
         data = tmp_ds['data']
         if add_vectors:
            data_vx = tmp_ds['data_vx']
            data_vy = tmp_ds['data_vy']
         if add_contour: 
            data_cnt = tmp_ds['data_cnt']

      #-------------------------------------------------------------------------
      # Unit conversions

      unit_str = ''
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
      
      
      if add_vectors:
         if vx_list[v]=='epfy': data_vx = data_vx/86400.
         if vy_list[v]=='epfz': data_vy = data_vy/86400.

      if add_contour: 
         if cnt_list[v]=='psitem': data_cnt = data_cnt/86400.

      #-------------------------------------------------------------------------
      # if 'time' in locals(): hc.print_time_length(time,indent=' '*6)

      if print_stats:
         msg = hc.print_stat(data,name=var[v],stat='naxsh',indent=(' '*6),compact=True,fmt='f')
         if add_vectors: 
            msg = hc.print_stat(data_vx,name=vx_list[v],stat='naxsh',indent=(' '*6),compact=True)
            msg = hc.print_stat(data_vy,name=vy_list[v],stat='naxsh',indent=(' '*6),compact=True)
         if add_contour:
            msg = hc.print_stat(data_cnt,name=cnt_list[v],stat='naxsh',indent=(' '*6),compact=True)
         
      #-------------------------------------------------------------------------
      # Create plot
      #-------------------------------------------------------------------------
   
      tres = copy.deepcopy(res)

      tres.tmYLMode = "Explicit"
      tres.tmYLValues = plev[::4].values
      tres.tmYLLabels = plev[::4].values / 1e2

      tres.sfXArray = lat.values
      tres.sfYArray = plev.values

      tres.lbTitleString   = f'[{unit_str}]'

      if plot_diff:
         if c==0:
            data_baseline_list.append(data)
         else:
            data = data - data_baseline_list[v]

      if plot_diff and c>0: 
         dummy = False
         if var[v]=='utendepfd': tres.cnLevels = np.arange(-10,10+2,2)/1e1
         if var[v]=='utendvtem': tres.cnLevels = np.arange(-10,10+2,2)/1e1
         if var[v]=='utendwtem': tres.cnLevels = np.arange(-10,10+2,2)/1e1
         if var[v]=='u'        : tres.cnLevels = np.arange(-8,8+2,2)
         if var[v]=='vtem'     : tres.cnLevels = np.arange(-10,10+2,2)*1e-1
         if var[v]=='wtem'     : tres.cnLevels = np.arange(-10,10+2,2)*1e-4
         if var[v]=='BUTGWSPEC': tres.cnLevels = np.arange(-0.1,0.1+0.02,0.02)
         if var[v]=='UTGWSPEC' : tres.cnLevels = np.arange(-0.1,0.1+0.02,0.02)
         if var[v]=='UTGWORO'  : tres.cnLevels = np.arange(-0.1,0.1+0.02,0.02)
      else:
         if var[v]=='utendepfd': tres.cnLevels = np.arange(-10,10+2,2)
         if var[v]=='utendvtem': tres.cnLevels = np.arange(-10,10+2,2)
         if var[v]=='utendwtem': tres.cnLevels = np.arange(-10,10+2,2)
         # if var[v]=='u'        : tres.cnLevels = np.arange(-8,8+2,2)
         if var[v]=='vtem'     : tres.cnLevels = np.arange(-10,10+2,2)/4
         if var[v]=='wtem'     : tres.cnLevels = np.arange(-10,10+2,2)*1e-2 /2
         if var[v]=='BUTGWSPEC': tres.cnLevels = np.arange(-0.3,0.3+0.03,0.03)
         if var[v]=='UTGWSPEC' : tres.cnLevels = np.arange(-0.3,0.3+0.03,0.03)
         if var[v]=='UTGWORO'  : tres.cnLevels = np.arange(-0.3,0.3+0.03,0.03)
         # if var[v]=='epfy'   : tres.cnLevels = np.arange(-100,100+5,5)*1e-6
         # if var[v]=='epfz'   : tres.cnLevels = np.arange(-100,100+5,5)*1e-6
         # if var[v]=='vtem'   : tres.cnLevels = np.arange(-100,100+5,5)*1e-6
         # if var[v]=='wtem'   : tres.cnLevels = np.arange(-100,100+5,5)*1e-6
         if var[v]=='psitem' : tres.cnLevels = np.arange(-20,20+2,2)/1e1
         # if var[v]=='depfydy': tres.cnLevels = np.arange(-100,100+5,5)*1e-6
         # if var[v]=='depfzdp': tres.cnLevels = np.arange(-100,100+5,5)*1e-6

      # nlev = 31
      # data_min = np.min(np.ma.masked_invalid(data.values))
      # data_max = np.max(np.ma.masked_invalid(data.values))
      # (cmin,cmax,cint) = ngl.nice_cntr_levels(data_min, data_max,
      #                    cint=None, max_steps=nlev,returnLevels=False, aboutZero=True )
      # tres.cnLevels = np.linspace(cmin,cmax,num=nlev)
      
      if hasattr(tres,'cnLevels'): tres.cnLevelSelectionMode = "ExplicitLevels"

      ip = v*num_case+c if var_x_case else c*num_var+v
      
      plot[ip] = ngl.contour(wks, np.ma.masked_invalid(data.values), tres)

      if add_contour:
         tmp_cres = copy.deepcopy(cres)
         tmp_cres.sfXArray = tres.sfXArray
         tmp_cres.sfYArray = tres.sfYArray
         tmp_cres.cnLevelSelectionMode = "ExplicitLevels"
         if cnt_list[v]=='u':          tmp_cres.cnLevels = np.arange(-80,80+10,10)
         if cnt_list[v]=='psitem':     tmp_cres.cnLevels = np.arange(-20,20+2,2)/1e1
         if cnt_list[v]=='BUTGWSPEC':  tmp_cres.cnLevels = np.arange(-55,55+10,10)*1e-7
         if not hasattr(tmp_cres,'cnLevels'): 
            nlev = 21
            data_min = np.min(np.ma.masked_invalid(data_cnt.values))
            data_max = np.max(np.ma.masked_invalid(data_cnt.values))
            (cmin,cmax,cint) = ngl.nice_cntr_levels(data_min, data_max,
                               cint=None, max_steps=nlev,returnLevels=False, aboutZero=False )
            tmp_cres.cnLevels = np.linspace(cmin,cmax,num=nlev)
            print(f'contour levels {cnt_list[v]:16} - min: {cmin}  max: {cmax}  nlev: {nlev}')
         neg_dash_contours(tmp_cres, tmp_cres.cnLevels)
         ngl.overlay(plot[ip], ngl.contour(wks, np.ma.masked_invalid(data_cnt.values), tmp_cres) )

      if add_vectors:
         tmp_vres = copy.deepcopy(vres)
         tmp_vres.vfXArray = tres.sfXArray
         tmp_vres.vfYArray = tres.sfYArray
         tmp_vres.vcMinFracLengthF = 0.5
         tmp_vres.vcMinDistanceF   = 0.025
         # tmp_vres.vcRefMagnitudeF  = 5.0
         # tmp_vres.vcRefLengthF     = 0.01
         if vx_list[v]=='epfy': tmp_vres.vcRefMagnitudeF  = 500e3
         vx_tmp = np.ma.masked_invalid(data_vx.values)
         vy_tmp = np.ma.masked_invalid(data_vy.values) *1e3
         ngl.overlay(plot[ip], ngl.vector(wks, vx_tmp, vy_tmp, tmp_vres) )
         # plot[ip] = ngl.vector(wks, vx_tmp, vy_tmp, tmp_vres)

      cstr = ''
      if plot_diff and c>0: cstr = f'diff'
      hs.set_subtitles(wks, plot[ip], name[c], '', var[v], center_sub_string=cstr, font_height=0.01)

      #-------------------------------------------------------------------------
      # # add reference line
      # lres = hs.res_xy()
      # lres.xyLineThicknessF = 1
      # lres.xyDashPattern = 0
      # lres.xyLineColor = 'black'
      # ngl.overlay(plot[ip],ngl.xy(wks, np.array([0,0]), np.array([-1e3,1e8]), lres))

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

if num_case > 1:
   layout = [num_var,num_case] if var_x_case else [num_case,num_var]
else:
   layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

pres = hs.setres_panel()
pres.nglPanelTop      =  0.93

# if common_colorbar: 
#    pres.nglPanelLabelBar = True
#    pres.lbTitleFontHeightF = 0.01
#    pres.nglPanelLabelBarLabelFontHeightF = 0.01
#    pres.lbTitlePosition = 'Bottom'
#    pres.lbTitleString = '' # what are the units of "Beres Tau"?

ngl.panel(wks,plot,layout,pres)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
