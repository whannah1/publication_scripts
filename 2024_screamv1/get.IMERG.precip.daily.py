import os, subprocess as sp
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
#-------------------------------------------------------------------------------
# see https://gpm.nasa.gov/data/directory
#-------------------------------------------------------------------------------

path_top = 'https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3'
# product = 'GPM_3IMERGDF.06'
product = 'GPM_3IMERGDF.07'

# dst_path = '/global/cscratch1/sd/whannah/Obs/IMERG/daily'
# dst_path = '/global/cscratch1/sd/whannah/Obs/GPM/daily'
dst_path = '/pscratch/sd/w/whannah/Obs/IMERG/daily'

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
class tcolor:
   ENDC,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
def run_cmd(cmd,suppress_output=False,execute=True):
  # if suppress_output : cmd = cmd + ' > /dev/null'
  msg = tcolor.GREEN + cmd + tcolor.ENDC
  print(f'\n{msg}')
  os.system(cmd)
  return
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# yr1,yr2 = 2005,2010
# yr1,yr2 = 2011,2015
# for yr in range(yr1,yr2+1):

# yrs = [2000,2001,2002,2003,2004, 2016,2017,2018,2019]
# for yr in yrs:

# yr1,mn1,yr2,mn2 = 2000,1,2000,12 # single year
# SCREAMv1 paper
# yr1,mn1 = 2016, 8; yr2,mn2=yr1,mn1+1
# yr1,mn1 = 2013, 4; yr2,mn2=yr1,mn1+1
# yr1,mn1 = 2013,10; yr2,mn2=yr1,mn1+1
# yr1,mn1 = 2020, 1; yr2,mn2=yr1,mn1+1


#-------------------------------------------------------------------------------
def add_yrmn(yrmn_list,yr1,mn1,yr2,mn2):
  all_done = False
  cnt = 0
  yrmn_beg = int(yr1*1e2+mn1)
  yrmn_end = int(yr2*1e2+mn2)
  while not all_done:
    yrmn = int( yrmn_beg + int(cnt/12)*100 + cnt )
    yrmn_list.append(yrmn)
    if yrmn==yrmn_end: all_done = True
    cnt+=1
#-------------------------------------------------------------------------------
yrmn_list = []

# yr1,mn1 = 2016, 8; yr2,mn2=yr1,mn1+1 ; add_yrmn(yrmn_list,yr1,mn1,yr2,mn2)
# yr1,mn1 = 2013, 4; yr2,mn2=yr1,mn1+1 ; add_yrmn(yrmn_list,yr1,mn1,yr2,mn2)
yr1,mn1 = 2013,10; yr2,mn2=yr1,mn1+1 ; add_yrmn(yrmn_list,yr1,mn1,yr2,mn2)
yr1,mn1 = 2020, 1; yr2,mn2=yr1,mn1+1 ; add_yrmn(yrmn_list,yr1,mn1,yr2,mn2)



#-------------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# for yr in range(yr1,yr2+1):
#   for mn in range(mn1,mn2+1):
for yrmn in yrmn_list:
  yr = int(yrmn/100)
  mn = yrmn-yr*100
    
  file_path = f'{path_top}/{product}/{yr}/{mn:02d}/'

  #-----------------------------------------------------------------------------
  ### debugging stuff
  # print(f'{yrmn}  {file_path}')
  # continue
  #-----------------------------------------------------------------------------

  ### list all files
  # run_cmd(f'wget -q -nH -nd "{file_path}" -O - | grep IMERG | cut -f4 -d\\"')

  ### skip a specific month?
  # if yr==2005 and mn==1: continue

  ### download files for each month
  cmd  = f'wget'
  cmd += f' --load-cookies ~/.urs_cookies'
  cmd += f' --save-cookies ~/.urs_cookies'
  cmd += f' --keep-session-cookies'
  # cmd += f' -O' # overwrite existing files
  cmd += f' -r -c -nH -nd -np -A nc4 --content-disposition'
  cmd += f' "{file_path}"'
  cmd += f' --directory-prefix={dst_path}'
  
  print('\n'+clr.GREEN+cmd+clr.END)

  run_cmd(cmd)

  # exit()

print('\ndone.')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
