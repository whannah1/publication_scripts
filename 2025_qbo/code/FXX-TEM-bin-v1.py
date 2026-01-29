import os, ngl, copy, xarray as xr, numpy as np, numba, glob
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import cmocean
# from dask.distributed import Client
# c = Client(n_workers=os.cpu_count()-2, threads_per_worker=1)
# c = Client(n_workers=2, threads_per_worker=1)
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
var_x,var_y,var_str = [],[],[]
def add_vars(var_x_name,var_y_name,var_str_in=''): 
   # if lev==-1: lev = np.array([0])
   var_x.append(var_x_name)
   var_y.append(var_y_name)
   var_str.append(var_str_in)
##------------------------------------------------------------------------------

scratch = '/global/cfs/cdirs/m4310/whannah/E3SM'
# scratch = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu'
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01'      ,n='L72',p=scratch,s='archive/atm/hist')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01',n='L80',p=scratch,s='archive/atm/hist')

# add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L72' ,n='E3SMv2 L72', p=scratch,s='archive/atm/hist')
# add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L80' ,n='E3SMv2 L80', p=scratch,s='archive/atm/hist')
# add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L128',n='E3SM L128',p=scratch,s='archive/atm/hist')
#-------------------------------------------------------------------------------

# lev = np.array([1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200])

add_vars('u','BUTGWSPEC','Convective GWD') # Beres U tendency - gravity wave spectrum
add_vars('u','utendepfd','EP Flux Divergence') # u-wind tendency due to TEM EP flux divergence

num_plot_col = len(var_y)

fig_type,fig_file = 'png','figs/FXX-TEM-bin-v1'
tmp_file_head = 'data/TEM-bin.v1'

lat1,lat2 = -5,5
# lat1,lat2 = -15,15

# htype,first_file,num_files = 'h2',0,365*1
# remap_str,search_str = 'remap_90x180','h0.tem.'; first_file,num_files = 0,12*20
remap_str,search_str = 'remap_90x180','h2.tem.'; first_file,num_files = 0,365*35

recalculate = False

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

wks = ngl.open_wks(fig_type,fig_file)
plot = [None]*(num_case*num_var)

res = hs.res_contour_fill()
res.vpHeightF = 0.5
res.tiXAxisFontHeightF  = 0.02
res.tiYAxisFontHeightF  = 0.02
res.tiYAxisString       = 'Pressure [hPa]'

res.lbLabelFontHeightF = 0.02

res.trYMaxF = 150

res.trYReverse = True
# res.xyYStyle = 'Log'  # doesn't work due to irregular spacing
# res.trYLog = True     # doesn't work due to irregular spacing

# res.tmYLMode = 'Explicit'
# res.tmYLValues = lev
# res.tmYLLabels = lev


def get_comp(case):
   comp = 'eam'
   if 'CESM' in case: comp = 'cam'
   return comp
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

      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      tmp_file = f'{tmp_file_head}.{case[c]}.{var_x[v]}.{var_y[v]}'
      if 'lat1' in locals(): tmp_file = tmp_file + f'.lat1_{lat1}.lat2_{lat2}'
      if 'lon1' in locals(): tmp_file = tmp_file + f'.lon1_{lon1}.lon2_{lon2}'
      tmp_file = f'{tmp_file}.nc'
      print(f'      tmp_file: {tmp_file}')

      if recalculate :
         file_path = f'{case_dir[c]}/{case[c]}/data_{remap_str}_tem/*.eam.{search_str}*'
         file_list = sorted(glob.glob(file_path))
         if 'first_file' in locals(): file_list = file_list[first_file:]
         if 'num_files'  in locals(): file_list = file_list[:num_files]
         # for ff in file_list: print(ff)
         # exit()
         #-------------------------------------------------------------------------
         # This doesn't work - throws many "H5" errors about not finding attributes
         # ds = xr.open_mfdataset(file_list,parallel=False)
         # X = ds[var_x[v]].sel(lat=slice(lat1,lat2),plev=slice(200e2,0))
         # Y = ds[var_y[v]].sel(lat=slice(lat1,lat2),plev=slice(200e2,0))
         # print(ds)
         # # ds.load()
         # exit()
         #-------------------------------------------------------------------------
         plev1,plev2 = 200e2,0
         #-------------------------------------------------------------------------
         ds = xr.open_dataset(file_list[0])
         X = ds[var_x[v]].sel(lat=slice(lat1,lat2),plev=slice(plev1,plev2))
         Y = ds[var_y[v]].sel(lat=slice(lat1,lat2),plev=slice(plev1,plev2))
         for f in file_list[1:]:
            ds = xr.open_dataset(f)
            X_tmp = ds[var_x[v]].sel(lat=slice(lat1,lat2),plev=slice(plev1,plev2))
            Y_tmp = ds[var_y[v]].sel(lat=slice(lat1,lat2),plev=slice(plev1,plev2))
            X = xr.concat([X, X_tmp], dim='time')
            Y = xr.concat([Y, Y_tmp], dim='time')
         #-------------------------------------------------------------------------

         hc.print_stat(X,name=f'{var_x[v]}',indent=' '*6,compact=True)
         hc.print_stat(Y,name=f'{var_y[v]}',indent=' '*6,compact=True)

         #-------------------------------------------------------------------------
         # bin the data
         #-------------------------------------------------------------------------
         if var_x[v]=='u': bin_min, bin_max, bin_spc = -20, 20, 1

         # Create the output variables for binning
         nbin    = np.round( ( bin_max - bin_min + bin_spc )/bin_spc ).astype(int)
         bins    = np.linspace(bin_min,bin_max,nbin)
         bin_coord = xr.DataArray( bins )
         nlev,ntime  = len(Y['plev']),len(Y['time'])
         bin_val = xr.DataArray( np.full((nbin,nlev),np.nan,dtype=Y.dtype), 
                                 coords=[ ('bin', bin_coord.values), ('plev', Y['plev'].values) ] )
         bin_cnt = xr.DataArray( np.zeros((nbin,nlev),dtype=Y.dtype), 
                                 coords=bin_val.coords, dims=bin_val.dims )
         # do some preliminary adjustments to the input data
         X_tmp = np.where(np.isfinite(X.values),X.values,bin_min-bin_spc*1e3)
         # Loop through bins - perform area weighted average and average in time
         for b in range(nbin):
            bin_bot = bin_min - bin_spc/2. + bin_spc*(b  )
            bin_top = bin_min - bin_spc/2. + bin_spc*(b+1)
            for k in range(nlev):
               condition = xr.DataArray( np.full(X[:,k,:].shape,False,dtype=bool), coords=X[:,k,:].coords )
               condition.values = ( X_tmp[:,k,:] >=bin_bot ) & ( X_tmp[:,k,:]  <bin_top )
               if np.sum(condition.values)>0 :
                  bin_val[b,k] = Y[:,k,:].where(condition,drop=True).sum(skipna=True).values
                  bin_cnt[b,k] = np.sum(condition.values)
         # Store the result in a dataset
         dims = ('bins','plev')
         bin_ds = xr.Dataset()
         bin_ds['bin_val'] = bin_val
         bin_ds['bin_cnt'] = bin_cnt
         bin_ds['bin_pct'] = bin_cnt/bin_cnt.sum()*1e2
         bin_ds.coords['bins'] = bin_coord
         bin_ds.coords['plev'] = xr.DataArray(Y['plev'])
         # write the result to a file
         print(' '*6+'writing to file: '+tmp_file)
         bin_ds['time'] = Y.time
         bin_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         bin_ds = xr.open_dataset( tmp_file )

      #-------------------------------------------------------------------------
      val_list.append(bin_ds['bin_val'].transpose().values)
      cnt_list.append(bin_ds['bin_pct'].values)
      lev_list.append(bin_ds['plev'].values)
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
      # tres.cnFillPalette = "MPL_viridis"
      tres.cnFillPalette = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )
      tres.sfYArray = lev_list[c]/1e2
      tres.sfXArray = bin_list[c]

      tres.tmYLMode = 'Explicit'
      tres.tmYLValues = lev_list[c]/1e2
      tres.tmYLLabels = lev_list[c]/1e2

      tres.cnLevelSelectionMode = 'ExplicitLevels'
      # tres.cnLevels = np.linspace(-50,50,num=21)
      data_min, data_max = -0.3,0.4
      nlev = 21
      (cmin,cmax,cint) = ngl.nice_cntr_levels(data_min, data_max, cint=None, 
                         max_steps=nlev, returnLevels=False, aboutZero=True )
      tres.cnLevels = np.linspace(cmin,cmax,num=nlev)

      if plot_diff and c>0:
         diff_mag_max = 0.08
         tres.cnLevels = np.linspace(-1*diff_mag_max,diff_mag_max,21)
         tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )

      tres.tiXAxisString = var_x[v]
      if var_x[v]=='u': tres.tiXAxisString = 'Zonal Wind [m/s]'

      # ip = v*num_case+c
      ip = c*num_var+v

      plot[ip] = ngl.contour(wks, np.ma.masked_invalid(val_list[c]) ,tres)

      hs.set_subtitles(wks, plot[ip], name[c], '', var_str[v], font_height=0.015)

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
