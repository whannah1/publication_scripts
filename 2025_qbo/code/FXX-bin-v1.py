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
def add_var(var_name,lev=-1): var.append(var_name); lev_list.append(lev)
var_x,var_y,lev_list = [],[],[]
def add_vars(var_x_name,var_y_name,lev=None): 
   # if lev==-1: lev = np.array([0])
   var_x.append(var_x_name)
   var_y.append(var_y_name)
   lev_list.append(lev)
##------------------------------------------------------------------------------
cscratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl'
gscratch = '/global/cfs/cdirs/m3312/whannah/2022-QBO-TEST'
pscratch = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu'
##------------------------------------------------------------------------------
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01',       n='E3SM control',     c='red'  ,p=gscratch,s='run')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01', n='E3SM L72 smoothed',c='green',p=gscratch,s='run')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01', n='E3SM L80 refined', c='blue' ,p=gscratch,s='run')

add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72.GW-MOD-0', n='E3SM L72 M0',c='red'  ,d=1,p=pscratch,s='run')
add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72.GW-MOD-1', n='E3SM L72 M1',c='green',d=1,p=pscratch,s='run')
#-------------------------------------------------------------------------------

lev = np.array([1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200])

# fields sugested by Yaga
# add_var('UTGWSPEC')  # U tendency - gravity wave spectrum
# add_var('BUTGWSPEC') # Beres U tendency - gravity wave spectrum
# add_var('UTGWORO')   # U tendency - orographic gravity wave drag
# add_var('BTAUE')     # Reynolds stress from Beres scheme for waves propagating East
# add_var('BTAUW')     # Reynolds stress from Beres scheme for waves propagating West

add_vars('U','BUTGWSPEC') # Beres U tendency - gravity wave spectrum


num_plot_col = 2


fig_type = "png"
fig_file = 'figs/FXX-bin-v1'

tmp_file_head = 'data/bin.v1'

lat1,lat2 = -15,15

htype,first_file,num_files = 'h2',0,365*1

recalculate = True

plot_diff = True

# Y-axis options
omit_bot,bot_k   = False,-2
omit_top,top_k   = False,30
print_stats      = False
print_profile    = False



#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var_y),len(case)

if 'lev' not in vars(): lev = np.array([0])

wks = ngl.open_wks(fig_type,fig_file)
plot = [None]*(num_case*num_var)

res = hs.res_contour_fill()
res.tiXAxisFontHeightF  = 0.02
res.tiYAxisFontHeightF  = 0.02
res.tiYAxisString       = 'Pressure [hPa]'

res.trYReverse = True
# res.xyYStyle = 'Log'  # doesn't work due to irregular spacing
# res.trYLog = True     # doesn't work due to irregular spacing

res.tmYLMode = 'Explicit'
res.tmYLValues = lev
res.tmYLLabels = lev


def get_comp(case):
   comp = 'eam'
   if 'CESM' in case: comp = 'cam'
   return comp
#---------------------------------------------------------------------------------------------------
# @numba.jit(nopython=True)
# def calculate_equiangular_area(lon,lat,lon_bnds,lat_bnds):
#    re = 6.37122e06  # radius of earth
#    nlat,nlon = len(lat),len(lon)
#    area = np.empty((nlat,nlon),np.float64)
#    for j in range(nlat):
#       for i in range(nlon):
#          dlon = np.absolute( lon_bnds[j,1] - lon_bnds[j,0] )
#          dlat = np.absolute( lat_bnds[j,1] - lat_bnds[j,0] )
#          dx = re*dlon*np.pi/180.
#          dy = re*dlat*np.pi/180.
#          area[j,i] = dx*dy
#    return area
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   # print(hc.tcolor.GREEN+'  var: '+var_y[v]+hc.tcolor.ENDC)
   print('  X / Y : '+hc.tcolor.MAGENTA+var_x[v]+hc.tcolor.ENDC+' / '+hc.tcolor.MAGENTA+var_y[v]+hc.tcolor.ENDC)
   bin_list = []
   lev_list = []
   cnt_list = []
   val_list = []
   for c in range(num_case):
      print(hc.tcolor.CYAN+'    case: '+case[c]+hc.tcolor.ENDC)

      data_dir_tmp,data_sub_tmp = None, None
      # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
      if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      if case_sub[c] is not None: data_sub_tmp = case_sub[c]
      
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      tmp_file = f'{tmp_file_head}.{case[c]}.{var_x[v]}.{var_y[v]}'
      if 'lat1' in locals(): tmp_file = tmp_file + f'.lat1_{lat1}.lat2_{lat2}'
      if 'lon1' in locals(): tmp_file = tmp_file + f'.lon1_{lon1}.lon2_{lon2}'
      tmp_file = f'{tmp_file}.nc'
      print(f'      tmp_file: {tmp_file}')

      if recalculate :
         case_obj = he.Case( name=case[c], 
                             data_dir=data_dir_tmp, 
                             data_sub=data_sub_tmp  )
         if 'lat1' in locals(): case_obj.lat1 = lat1; case_obj.lat2 = lat2
         if 'lon1' in locals(): case_obj.lon1 = lon1; case_obj.lon2 = lon2
         #-------------------------------------------------------------------------
         # read the data   

         area = case_obj.load_data('area',htype=htype,num_files=1,).astype(np.double)   

         Y = case_obj.load_data(var_y[v],htype=htype,
                                 first_file=first_file,
                                 num_files=num_files,lev=lev)
         X = case_obj.load_data(var_x[v],htype=htype,
                                 first_file=first_file,
                                 num_files=num_files,lev=lev)

         # U = case_obj_x.load_data('U',htype=htype,
         #                         first_file=first_file,
         #                         num_files=num_files,lev=lev)
         # V = case_obj_x.load_data(var_x[v],htype=htype,
         #                         first_file=first_file,
         #                         num_files=num_files,lev=lev)
         # X = np.sqrt( np.square(U) + np.square(V) )

         #-------------------------------------------------------------------------
         # if 'time' in data.dims : 
         #    time = data.time
         #    hc.print_time_length(data.time,indent=' '*6,print_span=True, print_length=False)

         X = X.resample(time='D').mean(dim='time')
         Y = Y.resample(time='D').mean(dim='time')

         hc.print_stat(X,name=f'{var_x[v]}',indent='    ',compact=True)
         hc.print_stat(Y,name=f'{var_y[v]}',indent='    ',compact=True)

         #-------------------------------------------------------------------------
         # bin the data
         #-------------------------------------------------------------------------
         if var_x[v]=='U': bin_min, bin_max, bin_spc = -40, 40, 5

         # Create the output variables for binning
         nbin    = np.round( ( bin_max - bin_min + bin_spc )/bin_spc ).astype(np.int)
         bins    = np.linspace(bin_min,bin_max,nbin)
         bin_coord = xr.DataArray( bins )
         nlev,ntime  = len(Y['lev']),len(Y['time'])
         bin_val = xr.DataArray( np.full((nbin,nlev),np.nan,dtype=Y.dtype), 
                                 coords=[ ('bin', bin_coord), ('lev', Y['lev']) ] )
         bin_cnt = xr.DataArray( np.zeros((nbin,nlev),dtype=Y.dtype), 
                                 coords=bin_val.coords, dims=bin_val.dims )
         # do some preliminary adjustments to the input data
         wgt, *__ = xr.broadcast(area, Y[:,0,:]) ; wgt = wgt.transpose('time','ncol')
         X_tmp = np.where(np.isfinite(X.values),X.values,bin_min-bin_spc*1e3)
         # Loop through bins - perform area weighted average and average in time
         for b in range(nbin):
            bin_bot = bin_min - bin_spc/2. + bin_spc*(b  )
            bin_top = bin_min - bin_spc/2. + bin_spc*(b+1)
            for k in range(nlev):
               condition = xr.DataArray( np.full(X[:,k,:].shape,False,dtype=bool), coords=X[:,k,:].coords )
               condition.values = ( X_tmp[:,k,:] >=bin_bot ) & ( X_tmp[:,k,:]  <bin_top )
               if np.sum(condition.values)>0 :
                  Ywgt_sum = (Y[:,k,:]*wgt).where(condition,drop=True).sum(dim='ncol',skipna=True)
                  wgt_sum  =           wgt .where(condition,drop=True).sum(dim='ncol',skipna=True)
                  bin_val[b,k] = ( Ywgt_sum / wgt_sum ).mean(dim='time',skipna=True).values
                  bin_cnt[b,k] = np.sum(condition.values)
         # Store the result in a dataset
         dims = ('bins','lev')
         bin_ds = xr.Dataset()
         bin_ds['bin_val'] = (dims, bin_val )
         bin_ds['bin_cnt'] = (dims, bin_cnt )
         bin_ds['bin_pct'] = (dims, bin_cnt/bin_cnt.sum()*1e2 )
         bin_ds.coords['bins'] = ('bins',bin_coord)
         bin_ds.coords['lev'] = ( 'lev', xr.DataArray(Y['lev']) )
         # write the result to a file
         print('writing to file: '+tmp_file)
         bin_ds['time'] = Y.time
         bin_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         bin_ds = xr.open_dataset( tmp_file )

      #-------------------------------------------------------------------------
      val_list.append(bin_ds['bin_val'].transpose().values)
      cnt_list.append(bin_ds['bin_pct'].values)
      lev_list.append(bin_ds['lev'].values)
      bin_list.append(bin_ds['bins'].values)

   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------

   if plot_diff:
      data_min,data_max = np.nanmin(val_list[0]),np.nanmax(val_list[0])
      for c in range(1,num_case): val_list[c] = val_list[c] - val_list[0]
      diff_mag_max = np.nanmax(np.absolute(val_list[1:]))
   else:
      data_min = np.min([np.nanmin(d) for d in val_list])
      data_max = np.max([np.nanmax(d) for d in val_list])
   
   for c in range(num_case):
      tres = copy.deepcopy(res)
      tres.cnFillPalette = "MPL_viridis"
      tres.sfYArray = lev_list[c]
      tres.sfXArray = bin_list[c]

      tres.cnLevelSelectionMode = 'ExplicitLevels'
      # tres.cnLevels = np.linspace(-50,50,num=21)
      nlev = 21
      (cmin,cmax,cint) = ngl.nice_cntr_levels(data_min, data_max, cint=None, 
                         max_steps=nlev, returnLevels=False, aboutZero=True )
      tres.cnLevels = np.linspace(cmin,cmax,num=nlev)

      if plot_diff and c>0:
         tres.cnLevels = np.linspace(-1*diff_mag_max,diff_mag_max,21)
         tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )

      tres.tiXAxisString = var_x[v]
      if var_x[v]=='U': tres.tiXAxisString = 'Zonal Wind [m/s]'

      # ip = v*num_case+c
      ip = c*num_var+v

      plot[ip] = ngl.contour(wks, np.ma.masked_invalid(val_list[c]) ,tres)

      hs.set_subtitles(wks, plot[ip], name[c], '', var_y[v], font_height=0.015)

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
