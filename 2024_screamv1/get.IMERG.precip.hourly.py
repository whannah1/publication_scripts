import os, subprocess as sp, datetime
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
#-------------------------------------------------------------------------------
# see https://gpm.nasa.gov/data/directory
#-------------------------------------------------------------------------------

path_top = 'https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3'
product = 'GPM_3IMERGHH.07'

# dst_path = '/global/cscratch1/sd/whannah/Obs/IMERG/daily'
# dst_path = '/global/cscratch1/sd/whannah/Obs/GPM/daily'
dst_path = '/pscratch/sd/w/whannah/Obs/IMERG/hourly'

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
def add_yrmndy(yrmndy_list,yr1,mn1,dy1,yr2,mn2,dy2):
  done = False
  cnt = 0
  # yrmndy_beg = int(yr1*1e4+mn1*1e2+dy1)
  # yrmndy_end = int(yr2*1e4+mn2*1e2+dy2)
  beg_date = datetime.date(yr1, mn1, dy1)
  end_date = datetime.date(yr2, mn2, dy2)
  while not done:
    cur_date = beg_date + datetime.timedelta(days=cnt)
    # day_of_year = cur_date.timetuple().tm_yday
    yrmndy = int( cur_date.year*1e4 + cur_date.month*1e2 + cur_date.day )
    yrmndy_list.append(yrmndy)
    if cur_date==end_date: done = True
    cnt+=1
#-------------------------------------------------------------------------------
yrmndy_list = []

# yr1,mn1,dy1=2016, 8, 1; yr2,mn2,dy2=2016, 9, 9 ; add_yrmndy(yrmndy_list,yr1,mn1,dy1,yr2,mn2,dy2)
yr1,mn1,dy1=2013, 4, 1; yr2,mn2,dy2=2013, 5,10 ; add_yrmndy(yrmndy_list,yr1,mn1,dy1,yr2,mn2,dy2)
# yr1,mn1,dy1=2013,10, 1; yr2,mn2,dy2=2013,11, 9 ; add_yrmndy(yrmndy_list,yr1,mn1,dy1,yr2,mn2,dy2)
# yr1,mn1,dy1=2020, 1,20; yr2,mn2,dy2=2020, 2,28 ; add_yrmndy(yrmndy_list,yr1,mn1,dy1,yr2,mn2,dy2)


# exit()

#-------------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# for yr in range(yr1,yr2+1):
#   for mn in range(mn1,mn2+1):
for yrmndy in yrmndy_list:
  yr = int(yrmndy/1e4)
  mn = int((yrmndy-yr*1e4)/1e2)
  dy = int(yrmndy-yr*1e4-mn*1e2)
  cur_date = datetime.date(yr, mn, dy)
  day_of_year = cur_date.timetuple().tm_yday
    
  file_path = f'{path_top}/{product}/{yr}/{day_of_year:03d}/'

  #-----------------------------------------------------------------------------
  ### debugging stuff
  # print(f'{yrmndy}  {file_path}')
  # continue
  #-----------------------------------------------------------------------------

  ### list all files
  # run_cmd(f'wget -q -nH -nd "{file_path}" -O - | grep IMERG | cut -f4 -d\\"')
  # exit()

  ### skip a specific month?
  # if yr==2005 and mn==1: continue

  ### download files for each month
  cmd  = f'wget'
  cmd += f' --load-cookies ~/.urs_cookies'
  cmd += f' --save-cookies ~/.urs_cookies'
  cmd += f' --keep-session-cookies'
  # cmd += f' -O' # overwrite existing files
  # cmd += f' -r -c -nH -nd -np -A nc4 --content-disposition'
  cmd += f' -r -c -nH -nd -np -A HDF5 --content-disposition'
  cmd += f' "{file_path}"'
  cmd += f' --directory-prefix={dst_path}'
  
  print('\n'+clr.GREEN+cmd+clr.END)

  run_cmd(cmd)

  # exit()

print('\ndone.')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
