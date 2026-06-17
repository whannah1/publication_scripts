from sohip_methods import *
import pymsis
case_root = '/pscratch/sd/w/whannah/scream_scratch/pm-gpu'
#---------------------------------------------------------------------------------------------------
''' NOTES
June 2026 update:
preliminary "SOHIP simulator" results show that we need both a longer path length
and we need to blend MSIS data into the curtain.
'''
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
add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-19.09.NN_420',n='256x2-eq-ind-v1',xtime='2023-06-19 21:00',xlat=0,  xlon=  90,tlat= -6.99,tlon=  84.74,slat= 10.24,slon=  94.16)
add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-21.02.NN_420',n='256x2-eq-ind-v1',xtime='2023-06-21 21:00',xlat=-5, xlon=  80,tlat= -3.05,tlon=  75.97,slat= 13.56,slon=  85.35)
add_case('2025-SOHIP-RRM-00.256x2-ptgnia-v1.2023-06-13.19',       n='256x2-ptgnia-v1',xtime='2023-06-14 02:00',xlat=-50,xlon= -60,tlat=-49.46,tlon= -60.24,slat=-51.66,slon= -28.91)
add_case('2025-SOHIP-RRM-00.256x2-sc-ind-v1.2023-06-21.09',       n='256x2-sc-ind-v1',xtime='2023-06-21 15:00',xlat=-50,xlon=  80,tlat=-52.49,tlon=  67.04,slat=-51.03,slon=  98.64)
add_case('2025-SOHIP-RRM-00.256x2-sc-pac-v1.2023-06-14.15',       n='256x2-sc-pac-v1',xtime='2023-06-15 04:00',xlat=-35,xlon=-135,tlat=-34.73,tlon=-136.73,slat=-43.76,slon=-114.47)
add_case('2025-SOHIP-RRM-00.256x2-se-pac-v1.2023-06-12.16',       n='256x2-se-pac-v1',xtime='2023-06-13 04:00',xlat=-50,xlon= -95,tlat=-49.60,tlon= -94.45,slat=-51.80,slon= -63.70)
add_case('2025-SOHIP-RRM-00.256x2-sw-ind-v1.2023-06-12.06',       n='256x2-sw-ind-v1',xtime='2023-06-12 19:00',xlat=-50,xlon=  45,tlat=-49.61,tlon=  45.20,slat=-51.79,slon=  75.97)

# add_case('2025-SOHIP-RRM-00.256x3-ptgnia-v1.2023-06-13.19',       n='256x3-ptgnia-v1',xlat=-50,xlon= -60,tlat=-49.46,tlon= -60.24,slat=  None,slon=   None)
#---------------------------------------------------------------------------------------------------
var_list,var_opts_list = [],[]
def add_var(var_name,**kwargs): 
    var_list.append(var_name);
    tmp_opts = {}
    for k, val in kwargs.items(): tmp_opts[k] = val
    var_opts_list.append(tmp_opts)
#---------------------------------------------------------------------------------------------------

curtain_data_root = '/global/cfs/cdirs/m4842/whannah/curtain_data'

#---------------------------------------------------------------------------------------------------

wgt_method = 'inv_dist'

path_len_km = 1200   # total path distance [km]
# path_len_km = 2000   # total path distance [km]

path_spc_km = 2      # spacing between interpolated path points [km]
path_ncells = 2      # number of cells to consider nearest to each point (ncll)


### use larger cell count for "background" profile
# path_ncells = 10 # ~ 9.6km radius
# path_ncells = 20 # ~12.5km radius
# path_ncells = 30 # ~15.4km radius
# path_ncells = 40 # ~17.6km radius
# path_ncells = 50 # ~19.78km radius
# path_ncells = 60 # ~21.64km radius
path_ncells = 80 # ~24.87km radius
wgt_method = 'area'

# target_heights = np.arange(10e3,55e3+250,200) # this is used for the initial curtain calculation

new_target_heights = np.arange(10e3,60e3+200,200) # small extension for testing
# new_target_heights = np.arange(10e3,200e3+200,200)

blend_altitude_beg = 50e3

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

    dt = datetime.datetime.strptime(case_opts['xtime'], '%Y-%m-%d %H:%M')
    target_time = cftime.DatetimeNoLeap(dt.year, dt.month, dt.day, dt.hour, dt.minute, 0)
    #---------------------------------------------------------------------------
    f_name = case_opts['n']
    f_time = case_opts['xtime'].replace(' ','_').replace(':','_')
    tmp_file = f'{curtain_data_root}/{f_name}.{f_time}.curtain.ncells_{path_ncells}.len_{int(path_len_km)}.spc_{int(path_spc_km)}'
    # if wgt_method == 'inv_dist': tmp_file += '.wgt_inv_dist'
    if wgt_method == 'area': tmp_file += '.wgt_area'
    tmp_file += '.nc'
    #---------------------------------------------------------------------------
    ds = xr.open_dataset(tmp_file)
    curtain_height_km = ds['height']/1e3
    #---------------------------------------------------------------------------
    print(); print(ds)
    # print(); print(ds['height'].values/1e3)
    # print()
    exit()
    #---------------------------------------------------------------------------
    msis_Tmid_idx = pymsis.Variable["TEMPERATURE"]
    msis_rhod_idx = pymsis.Variable["MASS_DENSITY"]
    #---------------------------------------------------------------------------
    # calculate MSIS profile
    f107 = 150; f107a = 150; ap = 7; aps = [[ap] * 7]
    
    msis_data = np.squeeze( pymsis.calculate( np.datetime64(case_opts['xtime']), 
                                              case_opts['tlon'],
                                              case_opts['tlat'], 
                                              new_target_heights, f107, f107a, aps) )
    msis_data_Tmid = msis_data[:,msis_Tmid_idx]
    msis_data_rhod = msis_data[:,msis_rhod_idx]

    #---------------------------------------------------------------------------
    msis_data2 = np.squeeze( pymsis.calculate( np.datetime64(case_opts['xtime']), 
                                              case_opts['slon']-1,
                                              case_opts['slat']+1, 
                                              new_target_heights, f107, f107a, aps) )
    msis_data_Tmid2 = msis_data2[:,msis_Tmid_idx]
    #---------------------------------------------------------------------------

    for k in range(len(msis_data_Tmid)):
        tmp_height = ds['height'].isel(height=k).values
        tmp1 = msis_data_Tmid[k]
        # tmp2 = ds['T_mid'].isel(height=k).sel(path_coord=0,time=target_time, method='nearest').values
        tmp2 = msis_data_Tmid2[k]
        diff = tmp2-tmp1
        print(f'  {tmp_height}  {tmp1:6.2f}  {tmp2:6.2f}  {diff}')
        # print(f'  {tmp1}')

    exit()
    #---------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
