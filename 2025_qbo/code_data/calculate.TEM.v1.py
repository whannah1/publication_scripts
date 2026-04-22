import os, ngl, copy, glob, xarray as xr, numpy as np
class tclr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
#-------------------------------------------------------------------------------
case,name,case_dir,case_sub,clr,dsh,mrk = [],[],[],[],[],[],[]
def add_case(case_in,n=None,p=None,s=None,d=0,c='black',m=0):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); name.append(n); case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
var,lev_list = [],[]
def add_var(var_name,lev=-1): var.append(var_name); lev_list.append(lev)
##------------------------------------------------------------------------------
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01')

# num_files = 32
num_files = 1

#-------------------------------------------------------------------------------

scratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl'

debug = False

#---------------------------------------------------------------------------------------------------
# constants
H     = 7e3          # m         assumed mean scale heightof the atmosphere
P0    = 101325       # Pa        surface pressure
Rd    = 287.058      # J/kg/K    gas constant for dry air
cp    = 1004.64      # J/kg/K    specific heat for dry air
g     = 9.80665      # m/s       global average of gravity at MSLP
a     = 6.37123e6    # m         Earth's radius
omega = 7.29212e-5   # 1/s       Earth's rotation rate
pi    = 3.14159

#---------------------------------------------------------------------------------------------------
num_case = len(case)
for c in range(num_case):
   print(f'    case: {tclr.CYAN}{case[c]}{tclr.END}')

   idir = f'{scratch}/{case[c]}/data_remap_90x180_prs'
   odir = f'{scratch}/{case[c]}/data_remap_90x180_tem'

   file_path = f'{idir}/*eam.h2.0*'
   file_list = sorted(glob.glob(file_path))
   file_list.remove(file_list[-1]) # last file is empty

   if 'num_files' in locals(): file_list = file_list[:num_files]

   if debug:
      print(hc.tcolor.RED+'WARNING! - only producing a single file for debugging!'+hc.tcolor.ENDC)

   #----------------------------------------------------------------------------
   # loop through files to calculate TEM terms 
   for f in file_list: 

      # exit(f)

      ds = xr.open_dataset(f)
      nlat = len(ds[ 'lat'].values)
      nlev = len(ds['plev'].values)

      lat = np.deg2rad(ds['lat'])
      cos_lat = np.cos(lat)

      dlat = xr.full_like(ds['lat'],np.nan)
      for j in range(1,nlat-1):
         dlat[j] = a*cos_lat[j] * ( lat[j+1] - lat[j-1] )

      dp = xr.full_like(ds['plev'],np.nan)
      for k in range(1,nlev-1):
         dp[k] = ds['plev'][k-1] - ds['plev'][k+1]

      z = -H*np.log(ds['plev']/P0)
      dz = xr.full_like(ds['plev'],np.nan)
      for k in range(1,nlev-1):
         dz[k] = z[k-1] - z[k+1]

      #-------------------------------------------------------------------------
      # basic zonal means and anomalies

      TH = ds['T'] * np.power( P0/ds['plev'], Rd/cp )
      
      TH_b = TH.mean(dim='lon')
      U_b  = ds['U'].mean(dim='lon')
      V_b  = ds['V'].mean(dim='lon')
      W_b  = ds['OMEGA'].mean(dim='lon')

      TH_p = TH          - TH_b
      U_p  = ds['U']     - U_b
      V_p  = ds['V']     - V_b
      W_p  = ds['OMEGA'] - W_b

      #-------------------------------------------------------------------------
      # EP flux vectors

      dTHdp = xr.full_like(TH_b,np.nan)
      dUdp  = xr.full_like( U_b,np.nan)
      for k in range(1,nlev-1):
         dTHdp[:,k,:] = ( TH_b[:,k-1,:] - TH_b[:,k+1,:] ) / dp[k]
         dUdp[:,k,:]  = (  U_b[:,k-1,:] -  U_b[:,k+1,:] ) / dp[k]

      # eddy stream function
      gamma = (V_p*TH_p).mean(dim='lon') / dTHdp

      fcor = 2*omega*np.sin(lat)

      dUdy = xr.full_like(U_b,np.nan)
      for j in range(1,nlat-1):
         dUdy[:,:,j] =  (  U_b[:,:,j+1]*cos_lat[j+1] \
                         - U_b[:,:,j-1]*cos_lat[j-1] ) / dlat[j]

      ### original version based on Gerber and Manzini
      F_y = a*cos_lat * (       dUdp *gamma - (U_p*V_p).mean(dim='lon') )
      F_z = a*cos_lat * ( (fcor-dUdy)*gamma - (U_p*W_p).mean(dim='lon') )

      ### alt definition from NOAA document
      # F_y = a*cos_lat * -(U_p*V_p).mean(dim='lon')
      # F_z = a*cos_lat * fcor * ( (V_p*TH_p).mean(dim='lon') / dTHdp )

      ### alt version based on Yaga IDL code
      ### for k=0,lev-1 do  rac(*,k)=rhobar(*,k)*6.37e6*cos(rlats)
      ### Fphi=rac*(uzbar*vpthpbar/thzbar-vpupbar)
      ### ac=6.37e6*cos(rlats)
      ### for k=0,lev-1 do begin
      ###     temp=DERIV(rlats,(ubar(*,k)*cos(rlats)))
      ###     Fz(*,k)=rac(*,k)*((f-ac^(-1.)*temp)*vpthpbar(*,k)/thzbar(*,k)-wpupbar(*,k))
      ###     temp=DERIV(rlats,Fphi(*,k)*cos(rlats))
      ###     Fphiphi(*,k)=ac^(-1.)*temp
      ### endfor
      ### for j=0,lat-1 do Fzz(j,*)=DERIV(z,Fz(j,*))
      ### DELF=((Fphiphi+Fzz)/rac)*24. *3600.  ;in m/s/d
      # dTHdz = xr.full_like(TH_b,np.nan)
      # dUdz  = xr.full_like( U_b,np.nan)
      # for k in range(1,nlev-1):
      #    dTHdz[:,k,:] = ( TH_b[:,k-1,:] - TH_b[:,k+1,:] ) / dz[k]
      #    dUdz[:,k,:]  = (  U_b[:,k-1,:] -  U_b[:,k+1,:] ) / dz[k]
      # rho = ds['plev']/(9.81*H)
      # wz = -ds['OMEGA']*H/ds['plev']
      # wz_b = wz.mean(dim='lon')
      # wz_p = wz - wz_b
      # gamma_dz = (V_p*TH_p).mean(dim='lon') / dTHdz
      # F_y = rho*a*cos_lat * (       dUdz *gamma_dz - (U_p*V_p).mean(dim='lon') )
      # F_z = rho*a*cos_lat * ( (fcor-dUdy)*gamma_dz - (U_p*wz_p).mean(dim='lon') )

      F_y = F_y.transpose('time','plev','lat')
      F_z = F_z.transpose('time','plev','lat')

      #-------------------------------------------------------------------------
      # debug print statements

      # print(); print(F_y)
      # print(); print(F_z)
      # exit()

      # hc.print_stat(TH_b,name='TH_b',stat='naxsh',indent=(' '*6),compact=True)
      # hc.print_stat( U_b,name='U_b',stat='naxsh',indent=(' '*6),compact=True)
      # hc.print_stat(,name='',stat='naxsh',indent=(' '*6),compact=True)
      # exit()

      #-------------------------------------------------------------------------
      # EP flux divergence
      
      dFydy = xr.full_like(F_y,np.nan)
      for j in range(1,nlat-1):
         dFydy[:,:,j] =  (  F_y[:,:,j+1]*cos_lat[j+1] \
                          - F_y[:,:,j-1]*cos_lat[j-1] ) / dlat[j]

      ### original version based on Gerber and Manzini
      dFzdp = xr.full_like(F_z,np.nan)
      for k in range(1,nlev-1):
         dFzdp[:,k,:] = ( F_z[:,k-1,:] - F_z[:,k+1,:] ) / dp[k]

      EP_div = dFydy + dFzdp

      ### alt version based on Yaga IDL code
      # dFzdz = xr.full_like(F_z,np.nan)
      # for k in range(1,nlev-1):
      #    dFzdz[:,k,:] = ( F_z[:,k-1,:] - F_z[:,k+1,:] ) / dz[k]

      # EP_div = dFydy + dFzdz

      #-------------------------------------------------------------------------
      # EP flux divergence transformed to log-pressure

      F_y_lp = F_y * ds['plev']/P0
      F_z_lp = F_z * -H/P0

      dFzdz_lp = xr.full_like(F_z,np.nan)
      for k in range(1,nlev-1):
         dFzdz_lp[:,k,:] = ( F_z_lp[:,k-1,:] - F_z_lp[:,k+1,:] ) / dz[k]

      dFydy_lp = xr.full_like(F_y,np.nan)
      for j in range(1,nlat-1):
         dFydy_lp[:,:,j] =  (  F_y_lp[:,:,j+1]*cos_lat[j+1] \
                             - F_y_lp[:,:,j-1]*cos_lat[j-1] ) / dlat[j]

      EP_div_lp = dFydy_lp + dFzdz_lp
      # EP_div_lp = EP_div * ds['plev']/P0

      #-------------------------------------------------------------------------
      # TEM meridional and vertical  velocities

      dgamma_dp = xr.full_like(U_b,np.nan)
      for k in range(1,nlev-1):
         dgamma_dp[:,k,:] = ( gamma[:,k-1,:] - gamma[:,k+1,:] ) / dp[k]

      dgamma_dy = xr.full_like(U_b,np.nan)
      for j in range(1,nlat-1):
         dgamma_dy[:,:,j] =  (  gamma[:,:,j+1]*cos_lat[j+1] \
                              - gamma[:,:,j-1]*cos_lat[j-1] ) / dlat[j]

      V_star = V_b - dgamma_dp
      W_star = W_b - dgamma_dy

      #-------------------------------------------------------------------------
      # TEM mass stream function

      dp_int = xr.full_like(ds['plev'],np.nan)
      for k in range(1,nlev-1):
         pint1 = ( ds['plev'][k-1] + ds['plev'][k-0] ) / 2.
         pint2 = ( ds['plev'][k-0] + ds['plev'][k+1] ) / 2.
         dp_int[k] =  pint1 - pint2

      gamma_mass = xr.full_like(U_b,np.nan)
      for k in range(1,nlev-1):
         tmp_integral = xr.full_like(U_b.mean(dim='plev'),0)
         for kk in range(1,k):
            tmp_integral[:,:] = tmp_integral[:,:] + ( V_b[:,kk,:]*dp_int[kk] - gamma[:,kk,] )
         gamma_mass[:,k,:] = ( (2*pi*a*cos_lat[:]/g) * tmp_integral[:,:] ).transpose('time','lat')

      #-------------------------------------------------------------------------
      # TEM northward and upward advection

      dUdt_y = V_star * ( fcor - dUdy )
      dUdt_z = -1 * W_star * dUdp

      #-------------------------------------------------------------------------
      # gravity wave tendencies

      # BUTGWSPEC = ds['BUTGWSPEC'].mean(dim='lon') # Beres U tendency - gravity wave spectrum  (convection)
      # UTGWSPEC  = ds['UTGWSPEC' ].mean(dim='lon') # C&M U tendency - gravity wave spectrum    (frontogenesis)
      # UTGWORO   = ds['UTGWORO'  ].mean(dim='lon') # U tendency - orographic gravity wave drag

      #-------------------------------------------------------------------------
      f_out = f.replace('.h2.','.h2.tem.').replace(idir,odir)
      ds_out = xr.Dataset()

      ds_out['vtem']         = V_star
      ds_out['wtem']         = W_star
      ds_out['psitem']       = gamma_mass 
      # ds_out['epfy']         = F_y
      # ds_out['epfz']         = F_z
      # ds_out['depfydy']      = dFydy 
      # ds_out['depfzdp']      = dFzdp 
      ds_out['depfydy']      = dFydy_lp 
      ds_out['depfzdz']      = dFzdz_lp 
      # ds_out['utendepfd']    = EP_div
      # ds_out['utendepfd_lp'] = EP_div_lp
      # ds_out['utendvtem']    = dUdt_y
      # ds_out['utendwtem']    = dUdt_z
      # ds_out['u']            = ds['U' ].mean(dim='lon')
      # ds_out['z']            = ds['Z3'].mean(dim='lon')

      ds_out['vtem']        .attrs['long_name'] = 'Transformed Eulerian mean northward wind'
      ds_out['wtem']        .attrs['long_name'] = 'Transformed Eulerian mean upward wind'
      ds_out['psitem']      .attrs['long_name'] = 'Transformed Eulerian mean mass stream function'
      # ds_out['epfy']        .attrs['long_name'] = 'Northward component of the Eliassen-Palm flux'
      # ds_out['epfz']        .attrs['long_name'] = 'Upward component of the Eliassen-Palm flux'
      # ds_out['depfydy']     .attrs['long_name'] = 'Meridional derivative of northward component of the Eliassen-Palm flux'
      # ds_out['depfzdp']     .attrs['long_name'] = 'Vertical derivative of upward component of the Eliassen-Palm flux'
      # ds_out['utendepfd']   .attrs['long_name'] = 'Tendency of eastward wind due to TEM Eliassen-Palm flux divergence'
      # ds_out['utendepfd_lp'].attrs['long_name'] = 'Tendency of eastward wind due to TEM Eliassen-Palm flux divergence (log-pressure)'
      # ds_out['utendvtem']   .attrs['long_name'] = 'Tendency of eastward wind due to TEM northward wind advection and coriolis'
      # ds_out['utendwtem']   .attrs['long_name'] = 'Tendency of eastward wind due to TEM upward wind advection'
      # ds_out['u']           .attrs['long_name'] = 'Zonal mean eastward wind'
      # ds_out['z']           .attrs['long_name'] = 'Geopotential Height'
      
      ds_out['vtem']        .attrs['units'] = 'm/s'
      ds_out['wtem']        .attrs['units'] = 'm/2'
      ds_out['psitem']      .attrs['units'] = 'kg/s'
      # ds_out['epfy']        .attrs['units'] = 'm3/s2'
      # ds_out['epfz']        .attrs['units'] = 'm3/s2'
      # ds_out['depfydy']     .attrs['units'] = 'm/s2'
      # ds_out['depfzdp']     .attrs['units'] = 'm/s2'
      # ds_out['utendepfd']   .attrs['units'] = 'm/s2'
      # ds_out['utendepfd_lp'].attrs['units'] = 'm/s2'
      # ds_out['utendvtem']   .attrs['units'] = 'm/s2'
      # ds_out['utendwtem']   .attrs['units'] = 'm/s2'
      # ds_out['u']           .attrs['units'] = 'm/s'

      # ds_out['BUTGWSPEC']  = BUTGWSPEC
      # ds_out['UTGWSPEC' ]  = UTGWSPEC
      # ds_out['UTGWORO'  ]  = UTGWORO

      print(' '*4+f'writing to file: {f_out}')
      ds_out.to_netcdf(path=f_out,mode='w')

      # exit()

      if debug:
         print(hc.tcolor.RED+'WARNING! - only producing a single file for debugging!'+hc.tcolor.ENDC)
         exit()

      #-------------------------------------------------------------------------

print(); print('done.'); print()
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------