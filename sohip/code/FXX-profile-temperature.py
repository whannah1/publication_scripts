from sohip_methods import *
case_root = '/pscratch/sd/w/whannah/scream_scratch/pm-gpu'
#-------------------------------------------------------------------------------
'''
python -i
exec( open( 'code/FXX-curtain-v2.py' ).read() )
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
# add_case('2025-SOHIP-RRM-00.256x2-sc-ind-v1.2023-06-21.09',       n='256x2-sc-ind-v1',xtime='2023-06-21 15:00',xlat=-50,xlon=  80,tlat=-52.49,tlon=  67.04,slat=-51.03,slon=  98.64)
add_case('2025-SOHIP-RRM-00.256x2-sc-pac-v1.2023-06-14.15',       n='256x2-sc-pac-v1',xtime='2023-06-15 04:00',xlat=-35,xlon=-135,tlat=-34.73,tlon=-136.73,slat=-43.76,slon=-114.47)
# add_case('2025-SOHIP-RRM-00.256x2-se-pac-v1.2023-06-12.16',       n='256x2-se-pac-v1',xtime='2023-06-13 04:00',xlat=-50,xlon= -95,tlat=-49.60,tlon= -94.45,slat=-51.80,slon= -63.70)
# add_case('2025-SOHIP-RRM-00.256x2-sw-ind-v1.2023-06-12.06',       n='256x2-sw-ind-v1',xtime='2023-06-12 19:00',xlat=-50,xlon=  45,tlat=-49.61,tlon=  45.20,slat=-51.79,slon=  75.97)

# add_case('2025-SOHIP-RRM-00.256x3-ptgnia-v1.2023-06-13.19',       n='256x3-ptgnia-v1',xlat=-50,xlon= -60,tlat=-49.46,tlon= -60.24,slat=  None,slon=   None)

#-------------------------------------------------------------------------------
var_list,var_opts_list = [],[]
def add_var(var_name,**kwargs): 
    var_list.append(var_name);
    tmp_opts = {}
    for k, val in kwargs.items(): tmp_opts[k] = val
    var_opts_list.append(tmp_opts)
#-------------------------------------------------------------------------------
add_var('T_mid',         vstr='T_mid [K]')
# add_var('p_mid',         vstr='Pressure [mb]')

#-------------------------------------------------------------------------------

fig_file = 'figs/FXX-profile-temperature.png'

curtain_data_root = '/global/cfs/cdirs/m4842/whannah/curtain_data'

path_len_km = 1200   # total path distance [km]
path_spc_km = 2      # path point spacing [km]
path_ncells = 2      # number of cells to average for each path point

# first_file,num_files = -1,1

polyfit_order = 7

max_lz = 2000

#-------------------------------------------------------------------------------

print_stats = False
# var_x_case  = False
# plot_diff   = False

num_plot_col = 1
# use_common_label_bar = True

#---------------------------------------------------------------------------------------------------
if case_list==[]: raise ValueError('ERROR - case list is empty!')
num_var,num_case = len(var_list),len(case_list)
#---------------------------------------------------------------------------------------------------
# set up figure object

# num_panels = num_case * num_var

# if 'num_plot_col' in locals() and num_plot_col is not None:
#     nrows = int(np.ceil(num_panels / float(num_plot_col)))
#     ncols = num_plot_col
# else:
#     (nrows,ncols) = (num_var,num_case) if var_x_case else (num_case,num_var)

# # fdx,fdy=10,15;figsize = (fdx*num_case,fdy*num_var) if var_x_case else (fdx*num_var,fdy*num_case)
# # fig = plt.figure(figsize=figsize)

# fig, axs = plt.subplots(nrows, ncols, figsize=(10*ncols, 5*nrows), squeeze=False)
# axs_flat = axs.flatten()

# Build figure
nrows = int(np.ceil(num_var / float(num_plot_col)))
ncols = num_plot_col
fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 5*nrows), squeeze=False)
axes_flat = axes.flatten()


#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
    var_opts = var_opts_list[v]
    hapy.print_line()
    print('  var: '+hapy.tclr.MAGENTA+var_list[v]+hapy.tclr.END)
    # SCREAM data lists
    scream_plev_list = []
    scream_zlev_list = []
    scream_data_list = []
    # COSMIC data list
    cosmic_lat_list = []
    cosmic_lon_list = []
    cosmic_plev_list = []
    cosmic_zlev_list = []
    cosmic_data_list = []
    # data_list,lat_list,lon_list = [],[],[]
    # std_list,cnt_list = [],[]
    for c in range(num_case):
        print('\n'+' '*4+'case: '+hapy.tclr.GREEN+case_list[c]+hapy.tclr.END)
        case_opts = case_opts_list[c]
        case_root = case_opts['p']
        case_sub  = case_opts['s']
        #-----------------------------------------------------------------------
        if var_opts.get('bkgd_path_ncells'):
            bkgd_path_ncells = var_opts['bkgd_path_ncells']
        else:
            bkgd_path_ncells = 50
        #-----------------------------------------------------------------------
        # read the data
        f_name = case_opts['n']
        f_time = case_opts['xtime'].replace(' ','_').replace(':','_')

        main_file_path = f'{curtain_data_root}/{f_name}.{f_time}.curtain.ncells_{path_ncells}.len_{int(path_len_km)}.spc_{int(path_spc_km)}.nc'
        bkgd_file_path = f'{curtain_data_root}/{f_name}.{f_time}.curtain.ncells_{bkgd_path_ncells}.len_{int(path_len_km)}.spc_{int(path_spc_km)}.wgt_area.nc'
        #-----------------------------------------------------------------------
        # print(f'\nfile_path: {file_path}\n')
        #-----------------------------------------------------------------------
        ds_main = xr.open_dataset(main_file_path)
        # ds_bkgd = xr.open_dataset(bkgd_file_path)
        #-----------------------------------------------------------------------
        dt = datetime.datetime.strptime(case_opts['xtime'], '%Y-%m-%d %H:%M')
        target_time = cftime.DatetimeNoLeap(dt.year, dt.month, dt.day, dt.hour, dt.minute, 0)
        #-----------------------------------------------------------------------
        ds_main = ds_main.sel(time=target_time, method='nearest') # select time nearest SOHIP observation
        # ds_bkgd = ds_bkgd.sel(time=target_time, method='nearest') # select time nearest SOHIP observation
        #-----------------------------------------------------------------------
        if var_list[v]=='wind_speed':
            data = np.sqrt( np.square(ds_main['horiz_winds_u']) + np.square(ds_main['horiz_winds_u']) )
            # bkgd = np.sqrt( np.square(ds_bkgd['horiz_winds_u']) + np.square(ds_bkgd['horiz_winds_u']) )
        else:
            data = ds_main[var_list[v]]
            # bkgd = ds_bkgd[var_list[v]]
        plev = ds_main['p_mid']
        #-----------------------------------------------------------------------
        # select center of curtain
        plev = plev.sel(path_coord=0,method='nearest')
        data = data.sel(path_coord=0,method='nearest')
        #-----------------------------------------------------------------------
        # print(); print(data)
        # print(); print(bkgd)
        # print()
        # exit()
        #-----------------------------------------------------------------------
        # if var_opts.get('method'):
        #     if var_opts['method']=='anomaly_path':
        #         data = data - data.mean(dim='path_coord')
        #     if var_opts['method']=='anomaly_poly':
        #         for i in range(len(data['path_coord'])):
        #             tmp_data = data[i,:]
        #             nan_mask = ~np.isnan(tmp_data)
        #             tmp_data = tmp_data[nan_mask]
        #             hght = tmp_data['height']
        #             bkgd_poly = np.polyval(np.polyfit(hght, tmp_data, polyfit_order), hght)
        #             data[i,nan_mask] = data[i,nan_mask] - bkgd_poly
        #     if var_opts['method']=='anomaly_disc':
        #         data = data - bkgd
        #     if var_opts['method']=='Ep':
        #         if var_opts.get('max_lz'):
        #             max_lz_loc = var_opts['max_lz']
        #         else:
        #             max_lz_loc = max_lz
        #         data = calc_gw_ep(data, poly_order=polyfit_order, min_lz=None, max_lz=max_lz_loc, bkgd=bkgd)
        #         print()
        #         hapy.print_stat(data)
        #-----------------------------------------------------------------------
        # # adjust units
        # if 'unit_fac' in var_opts:
        #     data = data * var_opts['unit_fac']
        #-----------------------------------------------------------------------
        # if print_stats: hapy.print_stat(data_interp_hgt,name=var_list[v],stat='naxsh',indent=' '*4,compact=True)
        #-----------------------------------------------------------------------
        scream_zlev_list.append( data['height'].values/1e3 )
        scream_plev_list.append( plev.values )
        scream_data_list.append( data.values )
        #-----------------------------------------------------------------------
        # load COSMIC data
        cosmic_file_list = cosmic_get_file_list(target_time, 
                           xlat=case_opts['tlat'], xlon=case_opts['tlon'], 
                           max_dist_km=1000, window_minutes=60,)
        tmp_cosmic_lat_list = []
        tmp_cosmic_lon_list = []
        tmp_cosmic_zlev_list = []
        tmp_cosmic_plev_list = []
        tmp_cosmic_data_list = []
        for f in cosmic_file_list:
            ds_cosmic = xr.open_dataset(f)
            tmp_cosmic_lat_list.append(ds_cosmic.attrs['lat'])
            tmp_cosmic_lon_list.append(ds_cosmic.attrs['lon'])
            tmp_cosmic_zlev_list.append(ds_cosmic['MSL_alt'].values)       # km
            tmp_cosmic_plev_list.append(ds_cosmic['Pres'].values)          # mb
            tmp_cosmic_data_list.append(ds_cosmic['Temp'].values+273.15)   # C => K
            
        cosmic_lat_list.append( tmp_cosmic_lat_list)
        cosmic_lon_list.append( tmp_cosmic_lon_list)
        cosmic_zlev_list.append( tmp_cosmic_zlev_list )
        cosmic_plev_list.append( tmp_cosmic_plev_list )
        cosmic_data_list.append( tmp_cosmic_data_list )
        #-----------------------------------------------------------------------
        # print(); print(cosmic_zlev)
        # print(); print(cosmic_plev)
        # print(); print(cosmic_data)
        # print()
        # hapy.print_stat(cosmic_zlev,name='zlev')
        # hapy.print_stat(cosmic_plev,name='plev')
        # hapy.print_stat(cosmic_data,name='data')
        # print()
        # exit()
    #-----------------------------------------------------------------------------------------------
    # data_min = np.min([np.nanmin(d) for d in scream_data_list])
    # data_max = np.max([np.nanmax(d) for d in scream_data_list])

    # for d in cosmic_data_list:
    #     print()
    #     print(d)
    #     print(type(d))
    #     for dd in d:
    #         print()
    #         print(dd)
    #         print(type(dd))
    #         print(dd.shape)
    #         # print(np.nanmin(dd))
    #         exit()
    #---------------------------------------------------------------------------
    scream_data_min = np.min([np.nanmin(d) for d in scream_data_list])
    scream_data_max = np.max([np.nanmax(d) for d in scream_data_list])
    # cosmic_data_min = np.min([np.nanmin(np.array(d)) for d in cosmic_data_list])
    # cosmic_data_min = np.max([np.nanmax(np.array(d)) for d in cosmic_data_list])
    # cosmic_data_min = np.nanmin(np.concatenate([np.atleast_1d(np.array(d)).ravel() for d in cosmic_data_list]))
    # cosmic_data_max = np.nanmax(np.concatenate([np.atleast_1d(np.array(d)).ravel() for d in cosmic_data_list]))
    cosmic_data_min = np.nanmin(np.concatenate([np.concatenate([np.atleast_1d(x) for x in d]) for d in cosmic_data_list]))
    cosmic_data_max = np.nanmax(np.concatenate([np.concatenate([np.atleast_1d(x) for x in d]) for d in cosmic_data_list]))
    # cosmic_data_min = np.min([np.nanmin(d) for d in cosmic_data_list])
    # cosmic_data_max = np.max([np.nanmax(d) for d in cosmic_data_list])
    data_min = np.min([scream_data_min,cosmic_data_min])
    data_max = np.max([scream_data_max,cosmic_data_max])
    #---------------------------------------------------------------------------
    # Create plot
    ax = axes_flat[v]
    for c in range(num_case):
        case_opts = case_opts_list[c]
        ls = case_opts.get('ls', None)
        ax.plot(scream_data_list[c],
                scream_zlev_list[c],
                color = 'black',
                linestyle = 'solid',
                linewidth = 1.5,)
                # label = case_opts['n'])
        cosmic_clr = ['red','green','blue','purple']
        for i in range(len(cosmic_data_list[c])):
            ax.plot(cosmic_data_list[c][i],
                    cosmic_zlev_list[c][i],
                    color = cosmic_clr[i],
                    linestyle = 'solid',
                    linewidth = 1.5,)

    # # Vertical zero line
    # ylims = (hght_list[0].min(), hght_list[0].max())
    # ax.axvline(x=0, color='black', linewidth=0.8, linestyle='solid')

    ax.set_xlim( data_min,
                 data_max )
    ax.set_ylim( np.min(scream_zlev_list),
                 np.max(scream_zlev_list) )

#---------------------------------------------------------------------------------------------------
# Finalize plot
plt.tight_layout(w_pad=1, h_pad=2)
fig.savefig(fig_file, dpi=150, bbox_inches='tight')
plt.close(fig)

print(f'\n{fig_file}\n')

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
