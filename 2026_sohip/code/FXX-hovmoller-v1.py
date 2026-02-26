from sohip_methods import *
case_root = '/pscratch/sd/w/whannah/scream_scratch/pm-gpu'
#-------------------------------------------------------------------------------
'''
python -i
exec( open( 'code/FXX-hovmoller-v1.py' ).read() )
'''
#-------------------------------------------------------------------------------
case_opts_list = []
case_list = []
def add_case(case_in,**kwargs):
    case_list.append(case_in)
    tmp_opts = {}
    for k, val in kwargs.items(): tmp_opts[k] = val
    tmp_opts['g'] = get_grid_file(case_in)
    tmp_opts['p'] = case_root
    tmp_opts['s'] = 'run'
    case_opts_list.append(tmp_opts)
#-------------------------------------------------------------------------------
# add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-19.09.NN_420',n='256x2-eq-ind-v1',xtime='2023-06-19 21:00',xlat=0,  xlon=  90,tlat= -6.99,tlon=  84.74,slat= 10.24,slon=  94.16)
# add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-21.02.NN_420',n='256x2-eq-ind-v1',xtime='2023-06-21 21:00',xlat=-5, xlon=  80,tlat= -3.05,tlon=  75.97,slat= 13.56,slon=  85.35)
# add_case('2025-SOHIP-RRM-00.256x2-ptgnia-v1.2023-06-13.19',       n='256x2-ptgnia-v1',xtime='2023-06-14 02:00',xlat=-50,xlon= -60,tlat=-49.46,tlon= -60.24,slat=  None,slon=   None)
add_case('2025-SOHIP-RRM-00.256x2-sc-ind-v1.2023-06-21.09',       n='256x2-sc-ind-v1',xtime='2023-06-21 15:00',xlat=-50,xlon=  80,tlat=-52.49,tlon=  67.04,slat=-51.03,slon=  98.64)
add_case('2025-SOHIP-RRM-00.256x2-sc-pac-v1.2023-06-14.15',       n='256x2-sc-pac-v1',xtime='2023-06-15 04:00',xlat=-35,xlon=-135,tlat=-34.73,tlon=-136.73,slat=-43.76,slon=-114.47)
add_case('2025-SOHIP-RRM-00.256x2-se-pac-v1.2023-06-12.16',       n='256x2-se-pac-v1',xtime='2023-06-13 04:00',xlat=-50,xlon= -95,tlat=-49.60,tlon= -94.45,slat=-51.80,slon= -63.70)
add_case('2025-SOHIP-RRM-00.256x2-sw-ind-v1.2023-06-12.06',       n='256x2-sw-ind-v1',xtime='2023-06-12 19:00',xlat=-50,xlon=  45,tlat=-49.61,tlon=  45.20,slat=-51.79,slon=  75.97)

# add_case('2025-SOHIP-RRM-00.256x3-ptgnia-v1.2023-06-13.19',       n='256x3-ptgnia-v1',xlat=-50,xlon= -60,tlat=-49.46,tlon= -60.24,slat=  None,slon=   None)

#-------------------------------------------------------------------------------
var_list,var_opts_list = [],[]
def add_var(var_name,**kwargs): 
    var_list.append(var_name);
    tmp_opts = {}
    for k, val in kwargs.items(): tmp_opts[k] = val
    var_opts_list.append(tmp_opts)
#-------------------------------------------------------------------------------

htype_2D = 'output.scream.2D.10min.INSTANT.nmins_x10'
htype_3D = 'output.scream.3D.10min.INSTANT.nmins_x10'

# add_var('T_mid',obs_var='T',vstr='Temperature',htype=htype_3D)

add_var('T_mid',      vstr='Temperature', htype=htype_3D)
# add_var('horiz_winds',vstr='wind speed', htype=htype_3D)

# add_var('ps')
# add_var('omega')
# add_var('qv')
# add_var('T_mid')
# add_var('z_mid')

#-------------------------------------------------------------------------------

fig_file = 'figs/FXX-hovmoller-v1.png'

# first_file,num_files = -1,1

#-------------------------------------------------------------------------------

print_stats = False
var_x_case  = False
# plot_diff   = False

# num_plot_col = int(np.sqrt(len(case_list)))
# use_common_label_bar = True

target_heights = np.arange(20e3,55e3+250,200)

#---------------------------------------------------------------------------------------------------
if case_list==[]: raise ValueError('ERROR - case list is empty!')
num_var,num_case = len(var_list),len(case_list)
#---------------------------------------------------------------------------------------------------
# set up figure object
(d1,d2) = (num_var,num_case) if var_x_case else (num_case,num_var)
fdx,fdy=10,15;figsize = (fdx*num_case,fdy*num_var) if var_x_case else (fdx*num_var,fdy*num_case)
fig = plt.figure(figsize=figsize)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
    var_opts = var_opts_list[v]
    hapy.print_line()
    print('  var: '+hapy.tclr.MAGENTA+var_list[v]+hapy.tclr.END)
    data_list,lat_list,lon_list = [],[],[]
    std_list,cnt_list = [],[]
    for c in range(num_case):
        print('\n'+' '*4+'case: '+hapy.tclr.GREEN+case_list[c]+hapy.tclr.END)
        case_opts = case_opts_list[c]
        case_root = case_opts['p']
        case_sub  = case_opts['s']
        #-----------------------------------------------------------------------
        # read the data
        htype = var_opts['htype']
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
        #-----------------------------------------------------------------------
        print()
        for f in file_list: print(' '*6+f'{hapy.tclr.YELLOW}{f}{hapy.tclr.END}')
        #-----------------------------------------------------------------------
        # ds = ux.open_mfdataset(case_opts['g'], file_list, data_vars='all')
        ds = xr.open_mfdataset(file_list, data_vars='all')
        # ds = ds.sel(time=target_time, method='nearest') # select time nearest SOHIP observation
        #-----------------------------------------------------------------------
        data = ds[var_list[v]]
        zmid = ds['z_mid']
        #-----------------------------------------------------------------------
        if var_list[v]=='horiz_winds':
            data = np.sqrt( np.square(data.isel(dim2=0)) + np.square(data.isel(dim2=1)) )
        #-----------------------------------------------------------------------
        # adjust units
        if 'unit_fac' in var_opts:
            data = data * var_opts['unit_fac']
        #-----------------------------------------------------------------------
        tlat,tlon = float(case_opts['tlat']),float(case_opts['tlon'])
        # data = data.sel(lat=tlat, lon=tlon, method="nearest")
        # zmid = zmid.sel(lat=tlat, lon=tlon, method="nearest")
        min_dist_ncol = find_closest_cells(tlat,tlon,data['lat'].values,data['lon'].values,num_cells=1)
        #-----------------------------------------------------------------------
        data = data.isel(ncol=min_dist_ncol[0],drop=True)
        zmid = zmid.isel(ncol=min_dist_ncol[0],drop=True)
        #-----------------------------------------------------------------------
        # print(); print(data)#; exit()
        #-----------------------------------------------------------------------
        # interpolate to height coordinate
        data_interp_hgt = hapy.interp_to_height( data, zmid, target_heights,
                                                 lev_dim='lev', height_dim='height',extrapolate=False)
        #-----------------------------------------------------------------------
        # print(); print(data_interp_hgt); exit()
        #-----------------------------------------------------------------------
        if print_stats: hapy.print_stat(data_interp_hgt,name=var_list[v],stat='naxsh',indent=' '*4,compact=True)
        #-----------------------------------------------------------------------
        data_interp_hgt = data_interp_hgt - data_interp_hgt.mean(dim='time')
        #-----------------------------------------------------------------------
        data_list.append( data_interp_hgt )
    #-----------------------------------------------------------------------------------------------
    data_min = np.min([np.nanmin(d) for d in data_list])
    data_max = np.max([np.nanmax(d) for d in data_list])
    #---------------------------------------------------------------------------
    # if plot_diff:
    #     tmp_data = copy.deepcopy(data_list)
    #     for c in range(num_case): tmp_data[c] = data_list[c] - data_list[diff_base]
    #     # diff_data_max = np.max([np.nanmax(d) for d in tmp_data])
    #     # diff_data_min = np.min([np.nanmin(d) for d in tmp_data])
    #     diff_data_max = np.max([np.nanmax(np.absolute(d)) for d in tmp_data])
    #     diff_data_min = -1*diff_data_max
    #     for c in range(num_case):
    #         if c!=diff_base:
    #             data_list[c] = data_list[c] - data_list[diff_base]
    #---------------------------------------------------------------------------
    # set color map
    cmap = 'viridis'
    # cmap = cmocean.cm.rain
    # cmap = cmocean.cm.balance
    # cmocean.cm.amp
    #---------------------------------------------------------------------------
    # set color bar levels
    clev = None
    # if 'precip' in var_list[v]      : clev = np.logspace( 0, 3, num=30).round(decimals=2)
    # if var_list[v]=='VapWaterPath'  : clev = np.arange(10,90+1,1)
    #---------------------------------------------------------------------------
    # Create plot
    for c in range(num_case):
        case_opts = case_opts_list[c]
        #-----------------------------------------------------------------------
        img_kwargs = {}
        # img_kwargs['origin']    = 'lower'
        img_kwargs['cmap']      = cmap
        img_kwargs['levels']    = 60

        # if plot_diff and c!=0:
        #     img_kwargs['cmap']   = cmocean.cm.balance
        #     img_kwargs['vmin']   = diff_data_min
        #     img_kwargs['vmax']   = diff_data_max
        #     clev = None

        if clev is not None:
            img_kwargs['norm'] = mcolors.BoundaryNorm(clev, ncolors=256)

        (d1,d2,ip) = (num_var,num_case,v*num_case+c+1) if var_x_case else (num_case,num_var,c*num_var+v+1)
        ax = fig.add_subplot(d1,d2,ip)        

        ax.set_title(case_opts['n'],    fontsize=20, loc='left')
        # ax.set_title(case_opts['xtime'],fontsize=20, loc='center')
        ax.set_title(var_opts['vstr'],  fontsize=20, loc='right')

        data = data_list[c].values
        hgth = data_list[c]['height'].values
        
        time_vals = data_list[c]['time'].values
        time_idx    = np.arange(len(time_vals))
        time_labels = [t.strftime('%m/%d %H:%M') for t in time_vals]


        contour = ax.contourf(time_idx, hgth, data.T, **img_kwargs)
        ax.set_box_aspect(0.6)

        tsec_diff_list = []
        for i,t in enumerate(time_vals):
            tsec_diff = abs((t - target_time).total_seconds())
            tsec_diff_list.append(tsec_diff)
        tsec_diff_min_index = min(range(len(tsec_diff_list)), key=lambda i: tsec_diff_list[i])

        line_x = np.array([tsec_diff_min_index,tsec_diff_min_index])
        line_y = np.array([min(hgth),max(hgth)])
        ax.plot(line_x, line_y, color='black', linewidth=1)

        tick_stride = max(1, len(time_idx) // 20)
        ax.set_xticks(time_idx[::tick_stride])
        ax.set_xticklabels(time_labels[::tick_stride], rotation=45, ha='right')

        cbar = fig.colorbar(contour, ax=ax, fraction=0.02, orientation='vertical')
        cbar.ax.tick_params(labelsize=15)

#---------------------------------------------------------------------------------------------------
# Finalize plot
fig.savefig(fig_file, dpi=200, bbox_inches='tight')
plt.close(fig)

print(f'\n{fig_file}\n')

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
