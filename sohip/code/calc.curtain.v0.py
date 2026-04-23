from sohip_methods import *
case_root = '/pscratch/sd/w/whannah/scream_scratch/pm-gpu'
#---------------------------------------------------------------------------------------------------
# This script is just for exploring the range of path_ncells
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
add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-19.09.NN_420',n='256x2-eq-ind-v1',xtime='2023-06-19 21:00',xlat=0,  xlon=  90,tlat= -6.99,tlon=  84.74,slat= 10.24,slon=  94.16)
# add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-21.02.NN_420',n='256x2-eq-ind-v1',xtime='2023-06-21 21:00',xlat=-5, xlon=  80,tlat= -3.05,tlon=  75.97,slat= 13.56,slon=  85.35)
# add_case('2025-SOHIP-RRM-00.256x2-ptgnia-v1.2023-06-13.19',       n='256x2-ptgnia-v1',xtime='2023-06-14 02:00',xlat=-50,xlon= -60,tlat=-49.46,tlon= -60.24,slat=-51.66,slon= -28.91)
# add_case('2025-SOHIP-RRM-00.256x2-sc-ind-v1.2023-06-21.09',       n='256x2-sc-ind-v1',xtime='2023-06-21 15:00',xlat=-50,xlon=  80,tlat=-52.49,tlon=  67.04,slat=-51.03,slon=  98.64)
# add_case('2025-SOHIP-RRM-00.256x2-sc-pac-v1.2023-06-14.15',       n='256x2-sc-pac-v1',xtime='2023-06-15 04:00',xlat=-35,xlon=-135,tlat=-34.73,tlon=-136.73,slat=-43.76,slon=-114.47)
# add_case('2025-SOHIP-RRM-00.256x2-se-pac-v1.2023-06-12.16',       n='256x2-se-pac-v1',xtime='2023-06-13 04:00',xlat=-50,xlon= -95,tlat=-49.60,tlon= -94.45,slat=-51.80,slon= -63.70)
# add_case('2025-SOHIP-RRM-00.256x2-sw-ind-v1.2023-06-12.06',       n='256x2-sw-ind-v1',xtime='2023-06-12 19:00',xlat=-50,xlon=  45,tlat=-49.61,tlon=  45.20,slat=-51.79,slon=  75.97)

num_files = 1
#---------------------------------------------------------------------------------------------------
# htype = 'output.scream.2D.10min.INSTANT.nmins_x10'
htype = 'output.scream.3D.10min.INSTANT.nmins_x10'
#---------------------------------------------------------------------------------------------------

# wgt_method = 'inv_dist'
# path_ncells = 2      # number of cells to consider nearest to each point (ncll)
path_len_km = 1200   # total path distance [km]
path_spc_km = 2      # spacing between interpolated path points [km]

# use larger cell count for "background" profile
wgt_method = 'area'
path_ncells_list = []
# path_ncells_list.append(  5) # ~ 7km radius
# path_ncells_list.append(  6) # ~ 7.6 km radius
# path_ncells_list.append(  7) # ~ 7.6 km radius
# path_ncells_list.append(  8) # ~ 8.0 km radius
# path_ncells_list.append(  9) # ~ 8.7 km radius
# path_ncells_list.append( 10) # ~ 9.6 km radius
# path_ncells_list.append( 16) # ~11.4 km radius
# path_ncells_list.append( 18) # ~12.3 km radius
# path_ncells_list.append( 20) # ~12.5 km radius
# path_ncells_list.append( 22) # ~13.6 km radius
# path_ncells_list.append( 24) # ~14.0 km radius
# path_ncells_list.append( 25) # ~14.4 km radius
# path_ncells_list.append( 26) # ~14.4 km radius
# path_ncells_list.append( 28) # ~14.8 km radius
# path_ncells_list.append( 30) # ~15.4 km radius
# path_ncells_list.append( 35) # ~17.10 km radius
# path_ncells_list.append( 40) # ~17.63 km radius
# path_ncells_list.append( 45) # ~18.85 km radius
# path_ncells_list.append( 50) # ~19.78 km radius
# path_ncells_list.append( 55) # ~20.94 km radius
# path_ncells_list.append( 60) # ~21.64 km radius
# path_ncells_list.append( 65) # ~22.38 km radius
# path_ncells_list.append( 70) # ~24.03 km radius
# path_ncells_list.append( 75) # ~24.27 km radius
# path_ncells_list.append( 80) # ~24.87 km radius
# path_ncells_list.append( 85) # ~?? km radius
# path_ncells_list.append( 90) # ~26.74 km radius
# path_ncells_list.append( 95) # ~27.24 km radius
# path_ncells_list.append(100) # ~27.76 km radius


#---------------------------------------------------------------------------------------------------
if case_list==[]: raise ValueError('ERROR - case list is empty!')
num_case = len(case_list)
#---------------------------------------------------------------------------------------------------
# for c in range(num_case):
c = 0
#-------------------------------------------------------------------------------
print(); print('case: '+hapy.tclr.GREEN+case_list[c]+hapy.tclr.END)
case_opts = case_opts_list[c]
case_root,case_sub = case_opts['p'],case_opts['s']
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
# read the data
file_path = f'{case_root}/{case_list[c]}/{case_sub}/*{htype}*'
file_list_all = sorted(glob.glob(file_path))
if 'first_file' in locals(): file_list_all = file_list_all[first_file:]
if 'num_files'  in locals(): file_list_all = file_list_all[:num_files]
#-------------------------------------------------------------------------------
# for f in file_list_all: print(' '*6+f'{hapy.tclr.YELLOW}{f}{hapy.tclr.END}')
if file_list_all==[]:
    print(); print('ERROR - file_list is empty!')
    print(); print(f'file_path: {file_path}')
    print()
#-------------------------------------------------------------------------------
dt = datetime.datetime.strptime(case_opts['xtime'], '%Y-%m-%d %H:%M')
target_time = cftime.DatetimeNoLeap(dt.year, dt.month, dt.day, dt.hour, dt.minute, 0)
file_list = reduce_file_list_to_target_time(file_list_all,target_time)
#-------------------------------------------------------------------------------
# expand list of files
for i,f in enumerate(file_list_all):
    if f==file_list[0]: break
file_list = file_list_all[i-1:i+1+1]
#-------------------------------------------------------------------------------
print()
for f in file_list: print(' '*2+f'{hapy.tclr.YELLOW}{f}{hapy.tclr.END}')
#-------------------------------------------------------------------------------
ds = xr.open_mfdataset(file_list, data_vars='all')
# ds = ds.sel(time=target_time, method='nearest') # select time nearest SOHIP observation
#---------------------------------------------------------------------------------------------------
for path_ncells in path_ncells_list:
    #---------------------------------------------------------------------------
    # find ncol indives and distance weighting for path interpolation
    # print('\n'+' '*2+'finding column indices for interpolation...')
    nlev = len(ds['lev'])
    (ncol_idx, dist_wgt, area_wgt) = find_path_ncol_wgt( path_npts, nlev, path_ncells, 
                                                         path_lat, path_lon,
                                                         ds['lat'].values, ds['lon'].values, 
                                                         ds['area'].isel(time=0,missing_dims='ignore').values )
    # ncol_idx = ncol_idx.astype(int)
    #---------------------------------------------------------------------------
    # print()
    hapy.print_stat(1./dist_wgt,name=f'path_ncells: {path_ncells:3}',stat='ax')
    # print()
    
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
