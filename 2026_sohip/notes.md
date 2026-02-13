
# SOHIP Observation date/time

```shell
2023-06-12 18:30 - South West Indian ocean
2023-06-13 04:00 - South East Pacific
2023-06-15 03:30 - South Central Pacific
2023-06-19 21:00 - Equatorial Indian Ocean
2023-06-21 14:30 - South Central Indian Ocean
2023-06-21 21:00 - Equatorial Indian Ocean
```

# ERA5 validation data

```shell
# Dates for RRM grids - init ~12 hour prior to observation
2023-06-14 01:30 - init 2023-06-13 19Z - Patagonia / South Atlantic
2023-06-12 18:30 - init 2023-06-12 06Z - South West Indian ocean
2023-06-13 04:00 - init 2023-06-12 16Z - South East Pacific
2023-06-15 03:30 - init 2023-06-14 15Z - South Central Pacific
2023-06-19 21:00 - init 2023-06-19 09Z - Equatorial Indian Ocean
2023-06-21 14:30 - init 2023-06-21 02Z - South Central Indian Ocean
2023-06-21 21:00 - init 2023-06-21 09Z - Equatorial Indian Ocean

# commands to get initialization data - copied from HICCUP notes
DATE=20230613; HR=19; nohup python -u ~/HICCUP/get_hindcast_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --output-root=/global/cfs/projectdirs/m4842/whannah/HICCUP > ${HOME}/E3SM_grid_support/2025-SOHIP-RRM/logs_hiccup/log.2025-SOHIP-get-ERA5.${DATE}.${HR} &
DATE=20230612; HR=06; nohup python -u ~/HICCUP/get_hindcast_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --output-root=/global/cfs/projectdirs/m4842/whannah/HICCUP > ${HOME}/E3SM_grid_support/2025-SOHIP-RRM/logs_hiccup/log.2025-SOHIP-get-ERA5.${DATE}.${HR} &
DATE=20230612; HR=16; nohup python -u ~/HICCUP/get_hindcast_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --output-root=/global/cfs/projectdirs/m4842/whannah/HICCUP > ${HOME}/E3SM_grid_support/2025-SOHIP-RRM/logs_hiccup/log.2025-SOHIP-get-ERA5.${DATE}.${HR} &
DATE=20230614; HR=15; nohup python -u ~/HICCUP/get_hindcast_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --output-root=/global/cfs/projectdirs/m4842/whannah/HICCUP > ${HOME}/E3SM_grid_support/2025-SOHIP-RRM/logs_hiccup/log.2025-SOHIP-get-ERA5.${DATE}.${HR} &
DATE=20230619; HR=09; nohup python -u ~/HICCUP/get_hindcast_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --output-root=/global/cfs/projectdirs/m4842/whannah/HICCUP > ${HOME}/E3SM_grid_support/2025-SOHIP-RRM/logs_hiccup/log.2025-SOHIP-get-ERA5.${DATE}.${HR} &
DATE=20230621; HR=02; nohup python -u ~/HICCUP/get_hindcast_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --output-root=/global/cfs/projectdirs/m4842/whannah/HICCUP > ${HOME}/E3SM_grid_support/2025-SOHIP-RRM/logs_hiccup/log.2025-SOHIP-get-ERA5.${DATE}.${HR} &
DATE=20230621; HR=09; nohup python -u ~/HICCUP/get_hindcast_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --output-root=/global/cfs/projectdirs/m4842/whannah/HICCUP > ${HOME}/E3SM_grid_support/2025-SOHIP-RRM/logs_hiccup/log.2025-SOHIP-get-ERA5.${DATE}.${HR} &

# same commands uising get_validation_data
DATE=20230613; HR=19; python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
DATE=20230612; HR=06; python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
DATE=20230612; HR=16; python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
DATE=20230614; HR=15; python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
DATE=20230619; HR=09; python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
DATE=20230621; HR=02; python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
DATE=20230621; HR=09; python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5

# get all dates at once:
python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230611 --final-date=20230622 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5


# get all dates - one at a time:
python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230611 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230612 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230613 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5

python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230614 --final-date=20230616 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230614 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230615 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230616 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5

python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230617 --final-date=20230622 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230617 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230618 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230619 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230620 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230621 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230622 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5


```