#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''
Calculate diagnostic variables for QBOi data submission using Gerber and Manzini 2016 and idl code from Yaga Richter
'''
import sys, glob, numpy as np, xarray as xr
import scipy.integrate as integrate

scratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl'
#-------------------------------------------------------------------------------
class tclr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
#-------------------------------------------------------------------------------
case,name,case_dir,case_sub,clr,dsh,mrk = [],[],[],[],[],[],[]
def add_case(case_in,n=None,p=None,s=None,d=0,c='black',m=0):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); name.append(n); case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
var,lev_list = [],[]
def add_var(var_name,lev=-1): var.append(var_name); lev_list.append(lev)
#-------------------------------------------------------------------------------



add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01')

num_files = 32
# num_files = 1

lev_name,T_name,U_name,V_name,W_name = 'plev','T','U','V','OMEGA'

t_axis,z_axis,y_axis,x_axis = 0,1,2,3
#-------------------------------------------------------------------------------
# define contsants:
p0 = 101325 # Pa surface pressure;
R = 287.058 # J K-1 kg-1, gas constant for dry air
Cp = 1004.64 # J K-1 kg-1, specific heat for dry air, at constant pressure
g0 = 9.80665 # m s-1, global average of gravity at mean sea level
a = 6371230 # m, Earth's radius
Omega = 7.29212e-5 # s-1, Earth's rotation rate (2*np.pi/86400)
# f = 2*Omega*np.sin(phi) # Coriolis parameter
pi = 3.14159 # mathematical constant
k = R/Cp
# H = 7000.0 # Assume scale height

#-------------------------------------------------------------------------------

# fname = str(sys.argv[1]) # '/work/kit/imk-asf/px5501/qboi_exp1/IMK_MESSy______20141101_0000_6h02_pl.nc'

num_case = len(case)
for c in range(num_case):
    print(f'  case: {tclr.CYAN}{case[c]}{tclr.END}')

    idir = f'{scratch}/{case[c]}/data_remap_90x180_prs'
    odir = f'{scratch}/{case[c]}/data_remap_90x180_tem'

    file_path = f'{idir}/*eam.h2.0*'
    file_list = sorted(glob.glob(file_path))
    file_list.remove(file_list[-1]) # last file is empty

    if 'num_files' in locals(): file_list = file_list[:num_files]

    for fname in file_list: 
        f_out = fname.replace('.h2.','.h2.tem_y.').replace(idir,odir)

        print()
        print(f'  {fname}')

        #-------------------------------------------------------------------
        ds = xr.open_dataset(fname)

        levels = ds[lev_name]
        lat = ds.lat
        lon = ds.lon
        T = ds[T_name]
        u = ds[U_name]
        v = ds[V_name]
        w = ds[W_name]

        # Calculate Altitude in meters
        # z = -np.log(levels/levels[-1])*H

        # Calculate latitudes in radians
        rlats = lat * pi/180.

        # Convert pressure from hPa to Pa
        pr0 = levels.data
        # Define pressure at each lat and lon
        ilon = len(lon)
        ilat = len(lat)
        # this is really slow!!
        # print('Here it is very slow!')
        p = np.array([([levels,]*ilat),]*ilon).T

        # Calculate Potential temperature
        th = T*(p0/p)**k

        # Calculate f
        f = 2. * Omega * np.sin(rlats)

        # Calculate Zonally Averaged quantities and zonal anomalies
        thbar, ubar, vbar, wbar = th.mean('lon'), u.mean('lon'), v.mean('lon'), w.mean('lon')

        # calculate vertical derivatives
        uzbar = np.gradient(ubar,pr0,axis=z_axis)
        thzbar = np.gradient(thbar,pr0,axis=z_axis)

        # Calculate zonal anomalies
        thp, up, vp, wp = th-thbar, u-ubar, v-vbar, w-wbar

        # Calculate Zonally Averaged quantities
        vpthpbar = (vp*thp).mean('lon')
        vpupbar  = (vp*up).mean('lon')
        wpupbar  = (wp*up).mean('lon')   # this is really omega'u'

        ac = a*np.cos(rlats)

        ac2d = np.array([ac,] * (len(levels.data)) )

        Fphi = ac2d*(uzbar*vpthpbar/thzbar-vpupbar)
        # Fphihat = p/p0*Fphi

        temp = np.gradient((ubar*np.cos(rlats)),rlats,axis=y_axis)
        Fz = ac.data*((f.data-ac.data**(-1.)*temp)*vpthpbar/thzbar-wpupbar)
        # Fzhat = -H/p0*Fz

        temp = np.gradient(Fphi*np.cos(rlats),rlats,axis=y_axis)
        Fphiphi = ac.data**(-1.)*temp

        # derivative for all latitudes:
        Fzz = np.gradient(Fz,pr0, axis=z_axis)

        DELF = ((Fphiphi+Fzz)/ac2d) # *24. *3600.  #in m/s/d

        # # new output for comparing calculation methods
        # dFydy = Fphiphi
        # dFzdp = Fzz    

        # Calculare Residual Velocities for all latitudes
        psi = vpthpbar/thzbar
        Psi = np.zeros(np.shape(psi))
        Psi = integrate.cumtrapz(vbar,pr0,initial=0,axis=z_axis)-psi
        Psi = 2*pi*ac.data/g0*Psi
        temp = np.gradient(psi,pr0,axis=z_axis)
        vres = vbar-temp

        # derivative for all levels:
        temp = np.gradient((vpthpbar*np.cos(rlats)/thzbar),rlats,axis=y_axis)
        wres = wbar + (1./ac.data)*temp
        # derivative for all levels:
        temp = np.gradient((ubar*np.cos(rlats)),rlats,axis=y_axis)
        VADV = -vres*((1./ac.data)*temp-f.data)


        WADV = -wres*uzbar # *24.*3600.
        # VADV = VADV*24.*3600.

        fy = Fphi # northward component of Eliassen-Palm flux
        fz = Fz # upward component of Eliassen-Palm flux
        vstar = vres # residual mean northward wind
        wstar = wres # residual mean upward wind
        psistar = Psi # mean mass stream function
        utenddivf = DELF # tendency of eastward wind due to Eliassen-Palm flux divergence
        utendw = WADV # tendency of eastward wind due to upward wind direction
        utendv = VADV # tendency of eastward wind due to northward wind direction and the Coriolis term

        #-------------------------------------------------------------------
        # print(); print(fy); print()
        # utenddivf = xr.DataArray(utenddivf, coords=[fy[lev_name], fy.lat], dims=[lev_name, 'lat'])
        # utenddivf = xr.DataArray(utenddivf, coords=[fy.time, fy[lev_name], fy.lat], dims=['time',lev_name, 'lat'])
        utenddivf = xr.DataArray(utenddivf, coords=fy.coords, dims=fy.dims)

        fy.name,               fy.attrs['long_name'],        fy.attrs['units'] = 'fy',        'northward EP-flux',                'N * m^-1'
        fz.name,               fz.attrs['long_name'],        fz.attrs['units'] = 'fz',        'upward EP-flux',                   'N * m^-1'
        vstar.name,         vstar.attrs['long_name'],     vstar.attrs['units'] = 'vstar',     'residual northward wind',          'm * s^-1'
        wstar.name,         wstar.attrs['long_name'],     wstar.attrs['units'] = 'wstar',     'residual upward wind',             'm * s^-1'
        psistar.name,     psistar.attrs['long_name'],   psistar.attrs['units'] = 'psistar',   'residual stream function',         'kg * s^-1'
        utenddivf.name, utenddivf.attrs['long_name'], utenddivf.attrs['units'] = 'utenddivf', 'u-tendency by EP-flux divergence', 'm * s^-2'
        
        # dFydy = xr.DataArray(dFydy, coords=fy.coords, dims=fy.dims)
        # dFzdp = xr.DataArray(dFzdp, coords=fy.coords, dims=fy.dims)
        # dFydy.name, dFydy.attrs['long_name'], dFydy.attrs['units'] = 'dFydy', 'EP-flux div Y-component', 'm * s^-2'
        # dFzdp.name, dFzdp.attrs['long_name'], dFzdp.attrs['units'] = 'dFzdp', 'EP-flux div Z-component', 'm * s^-2'

        #-------------------------------------------------------------------
        ds_out = xr.Dataset()

        ds_out['fy']        = fy
        ds_out['fz']        = fz
        ds_out['vstar']     = vstar
        ds_out['wstar']     = wstar
        ds_out['psistar']   = psistar
        # ds_out['utenddivf'] = utenddivf
        ds_out['utendepfd'] = utenddivf
        
        # ds_out['epfy']      = Fphi
        # ds_out['epfz']      = Fz
        # ds_out['depfydy']   = dFydy 
        # ds_out['depfzdz']   = dFzdp

        print(f'  {f_out}')
        ds_out.to_netcdf(f_out)

        #-------------------------------------------------------------------
        # fy.to_netcdf('fz_EMAC_QBOiExp4_r1i1p1_'+fname[50:58]+'.nc')
        # fz.to_netcdf('fy_EMAC_QBOiExp4_r1i1p1_'+fname[50:58]+'.nc')
        # vstar.to_netcdf('vstar_EMAC_QBOiExp4_r1i1p1_'+fname[50:58]+'.nc')
        # wstar.to_netcdf('wstar_EMAC_QBOiExp4_r1i1p1_'+fname[50:58]+'.nc')
        # psistar.to_netcdf('psistar_EMAC_QBOiExp4_r1i1p1_'+fname[50:58]+'.nc')
        # utenddivf.to_netcdf('utenddivf_EMAC_QBOiExp4_r1i1p1_'+fname[50:58]+'.nc')
        #-------------------------------------------------------------------

