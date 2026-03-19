from sohip_methods import *
case_root = '/pscratch/sd/w/whannah/scream_scratch/pm-gpu'
#---------------------------------------------------------------------------------------------------
'''
salloc --nodes 1 --qos interactive --time 04:00:00 --constraint cpu --account=e3sm
source activate ux_env
time python code/calc.curtain.v1.py
'''
#---------------------------------------------------------------------------------------------------
case_opts_list = []
case_list = []
def add_case(case_in,**kwargs):
    case_list.append(case_in)
    tmp_opts = {}
    for k, val in kwargs.items(): tmp_opts[k] = val
    tmp_opts['g'] = get_grid_file(case_in)
    tmp_opts['p'] = case_root
    tmp_opts['s'] = 'data_reduced'
    case_opts_list.append(tmp_opts)
#---------------------------------------------------------------------------------------------------
# add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-19.09.NN_420',n='256x2-eq-ind-v1',xtime='2023-06-19 21:00',xlat=0,  xlon=  90,tlat= -6.99,tlon=  84.74,slat= 10.24,slon=  94.16)
# add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-21.02.NN_420',n='256x2-eq-ind-v1',xtime='2023-06-21 21:00',xlat=-5, xlon=  80,tlat= -3.05,tlon=  75.97,slat= 13.56,slon=  85.35)
add_case('2025-SOHIP-RRM-00.256x2-ptgnia-v1.2023-06-13.19',       n='256x2-ptgnia-v1',xtime='2023-06-14 02:00',xlat=-50,xlon= -60,tlat=-49.46,tlon= -60.24,slat=-51.66,slon= -28.91)
# add_case('2025-SOHIP-RRM-00.256x2-sc-ind-v1.2023-06-21.09',       n='256x2-sc-ind-v1',xtime='2023-06-21 15:00',xlat=-50,xlon=  80,tlat=-52.49,tlon=  67.04,slat=-51.03,slon=  98.64)
# add_case('2025-SOHIP-RRM-00.256x2-sc-pac-v1.2023-06-14.15',       n='256x2-sc-pac-v1',xtime='2023-06-15 04:00',xlat=-35,xlon=-135,tlat=-34.73,tlon=-136.73,slat=-43.76,slon=-114.47)
# add_case('2025-SOHIP-RRM-00.256x2-se-pac-v1.2023-06-12.16',       n='256x2-se-pac-v1',xtime='2023-06-13 04:00',xlat=-50,xlon= -95,tlat=-49.60,tlon= -94.45,slat=-51.80,slon= -63.70)
# add_case('2025-SOHIP-RRM-00.256x2-sw-ind-v1.2023-06-12.06',       n='256x2-sw-ind-v1',xtime='2023-06-12 19:00',xlat=-50,xlon=  45,tlat=-49.61,tlon=  45.20,slat=-51.79,slon=  75.97)

# add_case('2025-SOHIP-RRM-00.256x3-ptgnia-v1.2023-06-13.19',       n='256x3-ptgnia-v1',xlat=-50,xlon= -60,tlat=-49.46,tlon= -60.24,slat=  None,slon=   None)
#---------------------------------------------------------------------------------------------------
var_list,var_opts_list = [],[]
def add_var(var_name,**kwargs): 
    var_list.append(var_name);
    tmp_opts = {}
    for k, val in kwargs.items(): tmp_opts[k] = val
    var_opts_list.append(tmp_opts)
#---------------------------------------------------------------------------------------------------

# htype = 'output.scream.2D.10min.INSTANT.nmins_x10'
htype = 'output.scream.3D.10min.INSTANT.nmins_x10'


add_var('T_mid')
add_var('p_mid')
add_var('qv')
add_var('horiz_winds_u')
add_var('horiz_winds_v')
add_var('omega')


# output_data_root = f'{case_root}/curtain_data'
output_data_root = '/global/cfs/cdirs/m4842/whannah/curtain_data'

#---------------------------------------------------------------------------------------------------

path_len_km = 1200   # total path distance [km]
path_spc_km = 2      # spacing between interpolated path points [km]
path_ncells = 2      # number of cells to consider nearest to each point (ncll)

#---------------------------------------------------------------------------------------------------

print_stats = False
var_x_case  = False

#---------------------------------------------------------------------------------------------------
if case_list==[]: raise ValueError('ERROR - case list is empty!')
num_var,num_case = len(var_list),len(case_list)
#---------------------------------------------------------------------------------------------------
for c in range(num_case):
    print(); print('case: '+hapy.tclr.GREEN+case_list[c]+hapy.tclr.END)
    case_opts = case_opts_list[c]
    case_root,case_sub = case_opts['p'],case_opts['s']
    #---------------------------------------------------------------------------
    # create the path between ISS and tangent positions
    print('\n'+' '*2+'creating path data...')
    # path_lat, path_lon, path_npts, path_coord = [None]*num_case, [None]*num_case, [None]*num_case, [None]*num_case
    slat,slon = float(case_opts['slat']),float(case_opts['slon'])
    tlat,tlon = float(case_opts['tlat']),float(case_opts['tlon'])
    # define path outward from location of tangent point
    (path_npts,path_coord,path_lat,path_lon) = calculate_path( tlat, tlon, slat, slon, path_len_km, path_spc_km )
    print(f'    path_len_km  : {path_len_km}')
    print(f'    path_spc_km  : {path_spc_km}')
    print(f'    path_npts    : {path_npts}')
    print(f'    path_lat     : {np.min(path_lat):6.2f}  -  {np.max(path_lat):6.2f}')
    print(f'    path_lon     : {np.min(path_lon):6.2f}  -  {np.max(path_lon):6.2f}')
    #---------------------------------------------------------------------------
    # read the data
    file_path = f'{case_root}/{case_list[c]}/{case_sub}/*{htype}*'
    
    file_list_all = sorted(glob.glob(file_path))
    if 'first_file' in locals(): file_list_all = file_list[first_file:]
    if 'num_files'  in locals(): file_list_all = file_list[:num_files]

    # for f in file_list_all: print(' '*6+f'{hapy.tclr.YELLOW}{f}{hapy.tclr.END}')

    if file_list_all==[]:
        print(); print('ERROR - file_list is empty!')
        print(); print(f'file_path: {file_path}')
        print()
    #-----------------------------------------------------------------------
    dt = datetime.datetime.strptime(case_opts['xtime'], '%Y-%m-%d %H:%M')
    target_time = cftime.DatetimeNoLeap(dt.year, dt.month, dt.day, dt.hour, dt.minute, 0)
    file_list = reduce_file_list_to_target_time(file_list_all,target_time)
    #-----------------------------------------------------------------------
    # expand list of files
    for i,f in enumerate(file_list_all):
        if f==file_list[0]: break
    file_list = file_list_all[i-1:i+1+1]
    #---------------------------------------------------------------------------
    print()
    for f in file_list: print(' '*2+f'{hapy.tclr.YELLOW}{f}{hapy.tclr.END}')
    #---------------------------------------------------------------------------
    ds = xr.open_mfdataset(file_list, data_vars='all')
    # ds = ds.sel(time=target_time, method='nearest') # select time nearest SOHIP observation
    #---------------------------------------------------------------------------
    # find ncol indives and distance weighting for path interpolation
    print('\n'+' '*2+'finding column indices for interpolation...')
    nlev = len(ds['lev'])
    (ncol_idx, dist_wgt) = find_path_ncol_wgt( path_npts, nlev, path_ncells, path_lat, path_lon,
                                               ds['lat'].values, ds['lon'].values )
    ncol_idx = ncol_idx.astype(int)
    # tmp_coords = {''}
    # ncol_idx = xr.DataArray(ncol_idx.astype(int),coords=tmp_coords)
    # dist_wgt = xr.DataArray(dist_wgt,            coords=tmp_coords)
    #---------------------------------------------------------------------------
    zmid = None
    ds_tmp = None
    #---------------------------------------------------------------------------
    for v in range(num_var):
        var_opts = var_opts_list[v]
        print('\n'+' '*2+'var: '+hapy.tclr.MAGENTA+var_list[v]+hapy.tclr.END)
        #-----------------------------------------------------------------------
        os.makedirs(output_data_root, exist_ok=True)
        f_name = case_opts['n']
        f_time = case_opts['xtime'].replace(' ','_').replace(':','_')
        tmp_file = f'{output_data_root}/{f_name}.{f_time}.curtain.ncells_{path_ncells}.len_{int(path_len_km)}.spc_{int(path_spc_km)}.nc'
        #-----------------------------------------------------------------------
        # print(); print(tmp_file); exit()
        #-----------------------------------------------------------------------
        tvar = var_list[v]
        if 'horiz_winds' in var_list[v]: tvar = 'horiz_winds'
        if var_list[v]=='p_mid': tvar = 'ps'
        data = ds[tvar]
        if v==0: zmid = ds['z_mid']
        #-----------------------------------------------------------------------
        if 'horiz_winds' in var_list[v]:
            if var_list[v]=='horiz_winds_u': data = data.isel(dim2=0)
            if var_list[v]=='horiz_winds_v': data = data.isel(dim2=1)
        if var_list[v]=='p_mid':
            p0 = 1e5
            data = ds['hyam']*p0 + ds['hybm']*ds['ps']
            data = data.transpose('time','ncol','lev')
            data.load()
        #-----------------------------------------------------------------------
        # adjust units
        if 'unit_fac' in var_opts: data = data * var_opts['unit_fac']
        #-----------------------------------------------------------------------
        # print('\n'+' '*4+'loading reduced data...')
        # data.load()
        # zmid.load()
        #-----------------------------------------------------------------------
        # print()
        # print(data)
        # print()
        # exit()
        #-----------------------------------------------------------------------
        # interp_coords = {'path_coord':path_coord,'lev':ds['lev']}
        # data_interp = xr.DataArray(np.zeros([path_npts,nlev]), coords=interp_coords)
        # zmid_interp = xr.DataArray(np.zeros([path_npts,nlev]), coords=interp_coords)
        # for n in range(path_npts):
        #     tmp_wgt = dist_wgt[n,:,np.newaxis]
        #     data_interp[n,:] = np.sum(data.isel(ncol=ncol_idx[n,:])*tmp_wgt) / np.sum(tmp_wgt)
        #     zmid_interp[n,:] = np.sum(zmid.isel(ncol=ncol_idx[n,:])*tmp_wgt) / np.sum(tmp_wgt)
        #-----------------------------------------------------------------------
        @numba.njit()
        def interpolate_to_path_loc( data, zmid, path_npts, ncol_idx, dist_wgt ):
            ntime = data.shape[0]
            nlev  = data.shape[2]
            data_interp = np.zeros(( ntime, path_npts, nlev ))
            zmid_interp = np.zeros(( ntime, path_npts, nlev ))
            for n in range(path_npts):
                tmp_wgt = dist_wgt[n,:,np.newaxis]
                data_interp[:,n,:] = np.sum( data[:,ncol_idx[n,:],:]*tmp_wgt, axis=1 ) / np.sum(tmp_wgt)
                zmid_interp[:,n,:] = np.sum( zmid[:,ncol_idx[n,:],:]*tmp_wgt, axis=1 ) / np.sum(tmp_wgt)
            return ( data_interp, zmid_interp )
        #-----------------------------------------------------------------------
        print(''+' '*4+'interpolating to path...')
        (data_interp,zmid_interp) = interpolate_to_path_loc( data.values, zmid.values, path_npts, ncol_idx, dist_wgt )
        interp_coords = {'time':data['time'],'path_coord':path_coord,'lev':ds['lev']}
        data_interp = xr.DataArray(data_interp, coords=interp_coords)
        zmid_interp = xr.DataArray(zmid_interp, coords=interp_coords)
        #-----------------------------------------------------------------------
        # @numba.njit()
        # def interp_to_path_loc( data, path_npts, ncol_idx, dist_wgt ): # data shape: (ntime, ncol, nlev)
        #     ntime = data.shape[0]
        #     nlev  = data.shape[2]
        #     data_interp = np.zeros((ntime, path_npts, nlev))
        #     for t in range(ntime):
        #         for n in range(path_npts):
        #             tmp_wgt = dist_wgt[np.newaxis,n,:,np.newaxis]
        #             # tmp_wgt = dist_wgt[n,:,np.newaxis]
        #             data_interp[t,n,:] = np.sum(data[t,ncol_idx[n,:],:]*tmp_wgt) / np.sum(tmp_wgt)
        #     return data_interp
        # #-----------------------------------------------------------------------
        # print('\n'+' '*4+'interpolating to path...')
        # interp_coords = {'time':data['time'],'path_coord':path_coord,'lev':ds['lev']}
        # data_interp = xr.DataArray(interp_to_path_loc( data.values, path_npts, ncol_idx, dist_wgt ), coords=interp_coords)
        # zmid_interp = xr.DataArray(interp_to_path_loc( zmid.values, path_npts, ncol_idx, dist_wgt ), coords=interp_coords)
        #-----------------------------------------------------------------------
        # print('\n'+' '*4+'interpolating to path...')
        # interp_coord = {'time':data['time'],'path_coord':path_coord,'lev':ds['lev']}
        # interp_shape = ( len(data['time']), path_npts, len(data['lev']) )
        # data_interp = xr.DataArray(np.zeros(interp_shape), coords=interp_coord)
        # zmid_interp = xr.DataArray(np.zeros(interp_shape), coords=interp_coord)
        # for n in range(path_npts):
        #     print(f'  n: {n}')
        #     tmp_wgt = dist_wgt[np.newaxis,n,:,np.newaxis]
        #     print(f'    !!!')
        #     data_interp[:,n,:] = (data.isel(ncol=ncol_idx[n,:])*tmp_wgt).sum(dim='ncol') / np.sum(tmp_wgt)
        #     print(f'    !!!')
        #     zmid_interp[:,n,:] = (zmid.isel(ncol=ncol_idx[n,:])*tmp_wgt).sum(dim='ncol') / np.sum(tmp_wgt)
        #     print(f'    !!!')
        #-----------------------------------------------------------------------
        # print()
        # print(data_interp)
        # print()
        # exit()
        #-----------------------------------------------------------------------
        # interpolate to height coordinate
        print(''+' '*4+'interpolating to height coordinate...')
        target_heights = np.arange(10e3,55e3+250,200)
        data_interp_hgt = hapy.interp_to_height( data_interp, zmid_interp, target_heights,
                                                 lev_dim='lev', height_dim='height',extrapolate=False)
        #-----------------------------------------------------------------------
        if print_stats: hapy.print_stat(data_interp_hgt,name=var_list[v],stat='naxsh',indent=' '*4,compact=True)
        #-----------------------------------------------------------------------
        # print()
        # print(data_interp_hgt)
        # print()
        # exit()
        #-----------------------------------------------------------------------
        # build output dataset
        if v==0:
            ds_tmp = xr.Dataset()
            ds_tmp['path_lat']  = path_lat
            ds_tmp['path_lon']  = path_lon
            ds_tmp.attrs['path_len_km'] = path_len_km
            ds_tmp.attrs['path_spc_km'] = path_spc_km
            ds_tmp.attrs['path_ncells'] = path_ncells
            ds_tmp.attrs['tlat']        = tlat
            ds_tmp.attrs['tlon']        = tlon
            ds_tmp.attrs['slat']        = slat
            ds_tmp.attrs['slon']        = slon
        #-----------------------------------------------------------------------
        # add current variable
        ds_tmp[var_list[v]] = data_interp_hgt
    #---------------------------------------------------------------------------
    print('\n'+' '*2+f'writing file...  => {tmp_file}')
    ds_tmp.to_netcdf(path=tmp_file,mode='w')
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
