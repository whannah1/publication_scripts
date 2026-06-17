from sohip_methods import *
#-------------------------------------------------------------------------------
'''
This script identifies how many COSMIC retreivals are avilable for a given set
of max_dist_km and window_minutes
'''
#-------------------------------------------------------------------------------

opts_list = []

opts_list.append({'xtime':'2023-06-19 21:00','xlat':  0,'xlon':  90})
opts_list.append({'xtime':'2023-06-21 21:00','xlat': -5,'xlon':  80})
opts_list.append({'xtime':'2023-06-14 02:00','xlat':-50,'xlon': -60})
opts_list.append({'xtime':'2023-06-21 15:00','xlat':-50,'xlon':  80})
opts_list.append({'xtime':'2023-06-15 04:00','xlat':-35,'xlon':-135})
opts_list.append({'xtime':'2023-06-13 04:00','xlat':-50,'xlon': -95})
opts_list.append({'xtime':'2023-06-12 19:00','xlat':-50,'xlon':  45})

print()
for win_length in [20,60]:
  print(f'window_minutes = {win_length}')
  for opts in opts_list:
    xtime, xlat, xlon = opts['xtime'], opts['xlat'], opts['xlon']

    dt = datetime.datetime.strptime(xtime, '%Y-%m-%d %H:%M')

    # target_time = cftime.DatetimeJulian(dt.year, dt.month, dt.day, dt.hour, dt.minute, 0) # for obs/IMERG
    target_time = cftime.DatetimeNoLeap(dt.year, dt.month, dt.day, dt.hour, dt.minute, 0)

    cosmic_file_list_1 = cosmic_get_file_list(target_time,                               window_minutes=win_length)
    cosmic_file_list_2 = cosmic_get_file_list(target_time, xlat, xlon, max_dist_km=2000, window_minutes=win_length)
    cosmic_file_list_3 = cosmic_get_file_list(target_time, xlat, xlon, max_dist_km=1000, window_minutes=win_length)
    cosmic_file_list_4 = cosmic_get_file_list(target_time, xlat, xlon, max_dist_km= 500, window_minutes=win_length)
      
    msg = f'  {xtime}   num files = '
    msg +=       f'{len(cosmic_file_list_1):3}'
    msg += ' / 2000km '+f'{len(cosmic_file_list_2):3}'
    msg += ' / 1000km '+f'{len(cosmic_file_list_3):3}'
    msg += ' / 500km ' +f'{len(cosmic_file_list_4):3}'
    print(msg)

  print()

#-------------------------------------------------------------------------------