from sohip_methods import *
case_root = '/pscratch/sd/w/whannah/scream_scratch/pm-gpu'
#---------------------------------------------------------------------------------------------------
case_opts_list = []
case_list = []
def add_case(case_in,**kwargs):
    case_list.append(case_in)
    tmp_opts = {}
    for k, val in kwargs.items(): tmp_opts[k] = val
    tmp_opts['g'] = get_grid_file(case_in)
    tmp_opts['p'] = case_root; tmp_opts['s'] = 'run'
    case_opts_list.append(tmp_opts)
#---------------------------------------------------------------------------------------------------
# add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-19.09.NN_420',n='256x2-eq-ind-v1',xtime='2023-06-19 21:00',xlat=0,  xlon=  90,tlat= -6.99,tlon=  84.74,slat= 10.24,slon=  94.16)
# add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-21.02.NN_420',n='256x2-eq-ind-v1',xtime='2023-06-21 21:00',xlat=-5, xlon=  80,tlat= -3.05,tlon=  75.97,slat= 13.56,slon=  85.35)
add_case('2025-SOHIP-RRM-00.256x2-ptgnia-v1.2023-06-13.19',       n='256x2-ptgnia-v1',xtime='2023-06-14 02:00',xlat=-50,xlon= -60,tlat=-49.46,tlon= -60.24,slat=  None,slon=   None)
add_case('2025-SOHIP-RRM-00.256x2-sc-ind-v1.2023-06-21.09',       n='256x2-sc-ind-v1',xtime='2023-06-21 15:00',xlat=-50,xlon=  80,tlat=-52.49,tlon=  67.04,slat=-51.03,slon=  98.64)
add_case('2025-SOHIP-RRM-00.256x2-sc-pac-v1.2023-06-14.15',       n='256x2-sc-pac-v1',xtime='2023-06-15 04:00',xlat=-35,xlon=-135,tlat=-34.73,tlon=-136.73,slat=-43.76,slon=-114.47)
add_case('2025-SOHIP-RRM-00.256x2-se-pac-v1.2023-06-12.16',       n='256x2-se-pac-v1',xtime='2023-06-13 04:00',xlat=-50,xlon= -95,tlat=-49.60,tlon= -94.45,slat=-51.80,slon= -63.70)
add_case('2025-SOHIP-RRM-00.256x2-sw-ind-v1.2023-06-12.06',       n='256x2-sw-ind-v1',xtime='2023-06-12 19:00',xlat=-50,xlon=  45,tlat=-49.61,tlon=  45.20,slat=-51.79,slon=  75.97)

# add_case('2025-SOHIP-RRM-00.256x3-ptgnia-v1.2023-06-13.19',       n='256x3-ptgnia-v1',xlat=-50,xlon= -60,tlat=-49.46,tlon= -60.24,slat=  None,slon=   None)

# htype = 'output.scream.2D.10min.INSTANT.nmins_x10'
htype = 'output.scream.3D.10min.INSTANT.nmins_x10'

#---------------------------------------------------------------------------------------------------
if case_list==[]: raise ValueError('ERROR - case list is empty!')
num_case = len(case_list)
#---------------------------------------------------------------------------------------------------
for c in range(num_case):
    print(); print('case: '+hapy.tclr.GREEN+case_list[c]+hapy.tclr.END)
    case_opts = case_opts_list[c]
    case_root,case_sub = case_opts['p'],case_opts['s']
    #---------------------------------------------------------------------------
    # read the data
    file_path = f'{case_root}/{case_list[c]}/{case_sub}/*{htype}*'
    file_list_all = sorted(glob.glob(file_path))
    if 'first_file' in locals(): file_list_all = file_list[first_file:]
    if 'num_files'  in locals(): file_list_all = file_list[:num_files]
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
    file_list = file_list_all[i-2:i+2+1]
    #---------------------------------------------------------------------------
    for src_file in file_list:
        # print(' '*2+f'{hapy.tclr.YELLOW}{src_file}{hapy.tclr.END}')
        #-----------------------------------------------------------------------
        src_file_name =  src_file.replace(f'{case_root}/{case_list[c]}/{case_sub}/','')
        out_file_root = f'{case_root}/{case_list[c]}/data_reduced'
        os.makedirs(out_file_root, exist_ok=True)
        dst_file_name = src_file_name.replace('.nc','.reduced.nc')
        out_file = f'{out_file_root}/{dst_file_name}'
        #-----------------------------------------------------------------------
        print()
        print(f'  src_file: {hapy.tclr.YELLOW}{src_file}{hapy.tclr.END}')
        print(f'  out_file: {hapy.tclr.CYAN}{out_file}{hapy.tclr.END}')
        # print()
        # exit()
        #-----------------------------------------------------------------------
        ds = xr.open_dataset(src_file)
        #-----------------------------------------------------------------------
        if src_file == file_list[0]:
            print('\n'+' '*2+'creating mask...')
            mdx = 20 # width in degrees of area centered on tanget point
            mlat1,mlat2 = case_opts['tlat']-mdx/2., case_opts['tlat']+mdx/2.
            mlon1,mlon2 = case_opts['tlon']-mdx/2., case_opts['tlon']+mdx/2.
            if mlon1<0: mlon1 = (mlon1+360)%360
            if mlon2<0: mlon2 = (mlon2+360)%360
            lat_full,lon_full = ds['lat'],ds['lon']
            mask = xr.ones_like(lat_full,dtype=bool)
            mask = mask & (lat_full>=mlat1) & (lat_full<=mlat2)
            mask = mask & (lon_full>=mlon1) & (lon_full<=mlon2)
            mask.load()
        #-----------------------------------------------------------------------
        # print(); print(ds)
        # print(); print(mask)
        # print()
        #-----------------------------------------------------------------------
        ds = ds.where(mask,drop=True)
        ds.compute()
        ds.to_netcdf(path=out_file,mode='w')
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
