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
# add_case('2025-SOHIP-RRM-00.256x2-sc-pac-v1.2023-06-14.15',       n='256x2-sc-pac-v1',xtime='2023-06-15 04:00',xlat=-35,xlon=-135,tlat=-34.73,tlon=-136.73,slat=-43.76,slon=-114.47)
# add_case('2025-SOHIP-RRM-00.256x2-se-pac-v1.2023-06-12.16',       n='256x2-se-pac-v1',xtime='2023-06-13 04:00',xlat=-50,xlon= -95,tlat=-49.60,tlon= -94.45,slat=-51.80,slon= -63.70)
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
add_var('T_mid',         vstr='Temperature [K]')
# add_var('T_mid',         vstr='Temperature [K]')
# add_var('T_mid',         vstr='Temperature [K]')
# add_var('T_mid',         vstr='Temperature [K]')

# add_var('density',       vstr='Density [kg/m3]')
# add_var('wind_speed',    vstr='Wind Speed [m/s]')
# add_var('p_mid',         vstr='Pressure [mb]')

#-------------------------------------------------------------------------------

fig_file = 'figs/FXX-stockwell-v1.png'

curtain_data_root = '/global/cfs/cdirs/m4842/whannah/curtain_data'

path_len_km = 1200   # total path distance [km]
path_spc_km = 2      # path point spacing [km]
path_ncells = 2      # number of cells to average for each path point

# vert_avg_scale = 10e3 # scale for rolling vertical average used for background removal [m]

#-------------------------------------------------------------------------------

print_stats = False
var_x_case  = True
# plot_diff   = False

# num_plot_col = 2#int(np.sqrt(len(case_list)))
# use_common_label_bar = True

#---------------------------------------------------------------------------------------------------
if case_list==[]: raise ValueError('ERROR - case list is empty!')
num_var,num_case = len(var_list),len(case_list)
#---------------------------------------------------------------------------------------------------
# set up figure object

num_panels = num_case * num_var

if 'num_plot_col' in locals() and num_plot_col is not None:
    nrows = int(np.ceil(num_panels / float(num_plot_col)))
    ncols = num_plot_col
else:
    (nrows,ncols) = (num_var,num_case) if var_x_case else (num_case,num_var)

# fdx,fdy=10,15;figsize = (fdx*num_case,fdy*num_var) if var_x_case else (fdx*num_var,fdy*num_case)
# fig = plt.figure(figsize=figsize)

fig, axs = plt.subplots(nrows, ncols*3, figsize=(16*ncols, 5*nrows), squeeze=False,
                        gridspec_kw={'width_ratios': [1,3,1]*ncols})
axs_flat = axs.flatten()

#---------------------------------------------------------------------------------------------------
for v in range(num_var):
    var_opts = var_opts_list[v]
    hapy.print_line()
    print('  var: '+hapy.tclr.MAGENTA+var_list[v]+hapy.tclr.END)
    data_list = []
    bkgd_list = []
    S_tx_list = []
    S_ep_list = []
    for c in range(num_case):
        print('\n'+' '*4+'case: '+hapy.tclr.GREEN+case_list[c]+hapy.tclr.END)
        case_opts = case_opts_list[c]
        case_root = case_opts['p']
        case_sub  = case_opts['s']
        #-----------------------------------------------------------------------
        # read the data
        f_name = case_opts['n']
        f_time = case_opts['xtime'].replace(' ','_').replace(':','_')
        file_path = f'{curtain_data_root}/{f_name}.{f_time}.curtain.ncells_{path_ncells}.len_{int(path_len_km)}.spc_{int(path_spc_km)}.nc'
        #-----------------------------------------------------------------------
        # print(f'\nfile_path: {file_path}\n')
        #-----------------------------------------------------------------------
        ds = xr.open_dataset(file_path)
        #-----------------------------------------------------------------------
        # print(); print(ds); print()
        # exit()
        #-----------------------------------------------------------------------
        dt = datetime.datetime.strptime(case_opts['xtime'], '%Y-%m-%d %H:%M')
        target_time = cftime.DatetimeNoLeap(dt.year, dt.month, dt.day, dt.hour, dt.minute, 0)
        ds = ds.sel(time=target_time, method='nearest') # select time nearest SOHIP observation
        num_path_pts = len(ds.path_coord)
        #-----------------------------------------------------------------------
        if var_list[v]=='wind_speed':
            data = np.sqrt( np.square(ds['horiz_winds_u']) + np.square(ds['horiz_winds_v']) )
        elif var_list[v]== 'density':
            Rd = 287.05 # J/kg/K
            data = ds['p_mid'] / ( Rd * ds['T_mid'] )
        else:
            data = ds[var_list[v]]
        #-----------------------------------------------------------------------
        # print(); print(data); print()
        # exit()
        #-----------------------------------------------------------------------
        # # calculate anomaly
        # if 'anomaly' in var_opts and var_opts['anomaly']:
        #     data = data - data.mean(dim='path_coord')
        #     bg = data.rolling(height=n_bg, center=True, min_periods=1).mean()
        #     bg = data.rolling(height=21, path_coord=31, center=True, min_periods=1).mean()
        #     data = data - bg
        #-----------------------------------------------------------------------
        # # adjust units
        # if 'unit_fac' in var_opts:
        #     data = data * var_opts['unit_fac']
        #-----------------------------------------------------------------------
        if print_stats: hapy.print_stat(data,name=var_list[v],stat='naxsh',indent=' '*4,compact=True)
        #-----------------------------------------------------------------------
        # reduce height dimension based on missing/NaN values
        min_hgt_pts = len(ds.height)
        for p in range(num_path_pts):
            tmp_data = data.isel(path_coord=p).values
            # Create a boolean mask for non-NaN values
            mask = ~np.isnan(tmp_data)
            tmp_data = tmp_data[mask]
            min_hgt_pts = min(min_hgt_pts,len(tmp_data))
        data = data.isel(height=slice(0,min_hgt_pts))
        #-----------------------------------------------------------------------
        # calculate vertical wavelength
        height = data.height.values
        dz = height[1] - height[0] # altitude spacing [m]
        vert_frq = np.fft.rfftfreq(len(height), d=dz) # Frequency axis (cycles per unit height), first value is 0
        vwav = 1.0 / vert_frq[1:]
        #-----------------------------------------------------------------------
        # obtain anomaly w/ polynomial fit
        dtmp = data.sel(path_coord=0)
        anom = dtmp.copy()*0.
        order = 7
        # if var_list[v]== 'density': order = 5

        # if v==0: order = 3
        # if v==1: order = 5
        # if v==2: order = 7
        # if v==3: order = 9

        bkgd = np.polyval(np.polyfit(height, dtmp.values, order), height)
        anom = ( dtmp - bkgd ) / bkgd # Remove background and normalize: x'/x_bar
        #-----------------------------------------------------------------------
        # # obtain anomaly w/ Savitzky-Golay filter
        # dtmp = data.sel(path_coord=0)
        # anom = dtmp.copy()*0.
        # # window of ~15 km, polynomial order 3 within the window
        # window_pts = int(15e3 / dz)
        # if window_pts % 2 == 0: window_pts += 1  # must be odd
        # bkgd = savgol_filter(dtmp.values, window_length=window_pts, polyorder=3)
        # anom = ( dtmp - bkgd ) / bkgd # Remove background and normalize: x'/x_bar
        #-----------------------------------------------------------------------
        # # Calculate Stockwell transform - average over many columns
        # S_tx_abs_path = np.empty([len(vwav)+1,min_hgt_pts,num_path_pts]) 
        # for p in range(num_path_pts):
        #     tmp_data = data.isel(path_coord=p).values
        #     S_tx_abs_path[:,:,p] = abs( st.st(tmp_data) ) # abs => magnitude of the transform for plotting
        # S_tx_abs = np.mean(S_tx_abs_path[:,:,:],axis=-1)
        #-----------------------------------------------------------------------
        # Calculate Stockwell transform for single column
        S_tx = st.st(anom)
        # if v==0:
        #     S_tx = st.st(anom)
        # if v==1:
        #     window = scipy.signal.windows.tukey(len(anom), alpha=0.1)  # 10% taper at each end
        #     anom_tapered = anom * window
        #     S_tx = st.st(anom_tapered)
        S_tx = S_tx[1:,:]
        #-----------------------------------------------------------------------
        # apply wave length mask
        min_lz, max_lz = dz*2, 5000
        vwav_mask = (vwav >= min_lz) & (vwav <= max_lz)
        S_tx = S_tx[vwav_mask,:]
        vwav = vwav[vwav_mask]
        #-----------------------------------------------------------------------
        S_sq = np.square(np.abs(S_tx))   # x^2 => proportional to power
        #-----------------------------------------------------------------------
        # if v==0:
        if True:
            # define cone of influence (COI) for S-transform: at distance d from the edge,
            # require N_cycles full cycles to fit between level and edge
            N_cycles   = 2        # increase this to be more conservative
            coi_top    = (height.max() - height) / N_cycles
            coi_bottom = (height - height.min()) / N_cycles
            coi        = np.minimum(coi_top, coi_bottom)
            coi_mask   = vwav[:, np.newaxis] <= coi[np.newaxis, :]
        #-----------------------------------------------------------------------
        # if v==1:
        #     # define cone of influence (COI) for S-transform: at distance d from the edge,
        #     # wavelengths longer than 2*d are unreliable        
        #     if v==0: vwav_crit = 2
        #     if v==1: vwav_crit = 4
        #     coi_top    = vwav_crit * (height.max() - height)   # max reliable wavelength near top
        #     coi_bottom = vwav_crit * (height - height.min())   # max reliable wavelength near bottom
        #     coi        = np.minimum(coi_top, coi_bottom)  # most restrictive at each level
        #     coi_mask   = vwav[:, np.newaxis] <= coi[np.newaxis, :]  # shape (n_wav, n_hgt)
        #-----------------------------------------------------------------------
        S_ep = 0.5 * np.sum(S_sq * coi_mask, axis=0) # GW potential energy
        #-----------------------------------------------------------------------
        # convert to xarray DataArray
        S_sq = xr.DataArray(S_sq,dims=('vwav','hght'))
        S_sq['vwav'] = vwav
        S_sq['hght'] = height
        S_ep = xr.DataArray(S_ep,dims=('hght'))
        S_ep['hght'] = height
        #-----------------------------------------------------------------------
        S_tx_list.append( S_sq )
        S_ep_list.append( S_ep )
        data_list.append( dtmp )
        bkgd_list.append( bkgd )
    #-----------------------------------------------------------------------------------------------
    data_min = np.min([np.nanmin(d) for d in S_tx_list])
    data_max = np.max([np.nanmax(d) for d in S_tx_list])
    #---------------------------------------------------------------------------
    print()
    print(f'  data_min: {data_min}')
    print(f'  data_max: {data_max}')
    print()
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

    # clev = np.logspace( -4, 1, num=60)
    clev = np.logspace( np.log10(data_min), np.log10(data_max), num=60)
    #---------------------------------------------------------------------------
    # Create plot
    for c in range(num_case):
        case_opts = case_opts_list[c]
        #-----------------------------------------------------------------------
        img_kwargs = {}
        # img_kwargs['origin']    = 'lower'
        img_kwargs['cmap']      = cmap
        img_kwargs['levels']    = clev
        # img_kwargs['norm'] = mcolors.LogNorm(vmin=1e-5, vmax=10)
        img_kwargs['norm'] = mcolors.LogNorm(vmin=np.min(clev), vmax=np.max(clev))
        # img_kwargs['norm'] = mcolors.LogNorm(vmin=np.log10(data_min), vmax=np.log10(data_max))

        # if plot_diff and c!=0:
        #     img_kwargs['cmap']   = cmocean.cm.balance
        #     img_kwargs['vmin']   = diff_data_min
        #     img_kwargs['vmax']   = diff_data_max
        #     clev = None

        # if clev is not None:
        #     img_kwargs['norm'] = mcolors.BoundaryNorm(clev, ncolors=256)
        

        # (d1,d2,ip) = (num_var,num_case,v*num_case+c+1) if var_x_case else (num_case,num_var,c*num_var+v+1)
        # ax = fig.add_subplot(d1,d2,ip)    

        fontsize = 15

        S_tx = S_tx_list[c].values

        ip = v*num_case+c if var_x_case else c*num_var+v
        ax = axs_flat[ip*3+1]

        ax.set_title(case_opts['n'],    fontsize=fontsize, loc='left')
        ax.set_title(case_opts['xtime'],fontsize=fontsize, loc='center')
        ax.set_title(var_opts['vstr'],  fontsize=fontsize, loc='right')

        vwav = S_tx_list[c]['vwav'].values
        hgth = S_tx_list[c]['hght'].values / 1e3
        contour = ax.contourf(vwav, hgth, S_tx.T, **img_kwargs)

        ax.set_box_aspect(0.6)
        ax.set_xscale('log')

        # ax.set_xlim([min_lz, max_lz])
        ax.set_xlim(right=max_lz)

        cbar = fig.colorbar(contour, ax=ax, fraction=0.02, orientation='vertical')
        cbar.ax.tick_params(labelsize=fontsize)

        # ax.set_xlabel('Distance from Tanget Point [km]',fontsize=fontsize)
        ax.set_xlabel('vertical wavelength [m]',fontsize=fontsize)
        ax.set_ylabel('Altitude [km]',          fontsize=fontsize)

        ax_ab = axs_flat[ip*3]
        ax_ab.plot(data_list[c], hgth, color='C0', label='anom')
        ax_ab.set_ylim(ax.get_ylim())
        ax_ab.set_xlabel('anom',          fontsize=fontsize)
        ax_ab.set_ylabel('Altitude [km]', fontsize=fontsize)
        ax_ab.tick_params(labelsize=fontsize)
        ax_ab_bkgd = ax_ab.twiny()
        ax_ab_bkgd.plot(bkgd_list[c], hgth, color='C1', label='bkgd')
        ax_ab_bkgd.set_xlabel(var_opts['vstr'], fontsize=fontsize)
        ax_ab_bkgd.tick_params(labelsize=fontsize)
        ax_ab.legend(handles=ax_ab.get_lines()+ax_ab_bkgd.get_lines(),
                     labels=['total field','bkgd'], fontsize=fontsize-2)

        ax_ep = axs_flat[ip*3+2]
        ax_ep.plot(S_ep_list[c].values, hgth)
        ax_ep.set_ylim(ax.get_ylim())
        ax_ep.set_xlabel('GW $E_p$',   fontsize=fontsize)
        ax_ep.set_ylabel('Altitude [km]', fontsize=fontsize)
        ax_ep.tick_params(labelsize=fontsize)

#---------------------------------------------------------------------------------------------------
# Finalize plot
plt.tight_layout(w_pad=1, h_pad=2)
fig.savefig(fig_file, dpi=200, bbox_inches='tight')
plt.close(fig)

print(f'\n{fig_file}\n')

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
