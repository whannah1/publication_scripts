import os, subprocess as sp

#-------------------------------------------------------------------------------
# case_list = []
# case_list.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne30pg2_ne30pg2.NN_32.SSTP_0K')
# case_list.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne45pg2_ne45pg2.NN_64.SSTP_0K')
# case_list.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne60pg2_ne60pg2.NN_128.SSTP_0K')
# case_list.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne90pg2_ne90pg2.NN_256.SSTP_0K')
# case_list.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne120pg2_ne120pg2.NN_512.SSTP_0K')
# case_list.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne30pg2_ne30pg2.NN_32.SSTP_0K.ALT-NCPL_72')
# case_list.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne45pg2_ne45pg2.NN_64.SSTP_0K.ALT-NCPL_72')
# case_list.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne60pg2_ne60pg2.NN_128.SSTP_0K.ALT-NCPL_72')
# case_list.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne90pg2_ne90pg2.NN_256.SSTP_0K.ALT-NCPL_72')
# case_list.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne120pg2_ne120pg2.NN_512.SSTP_0K.ALT-NCPL_72')

# print()
# for case in case_list:
#   res = case.split('.')[4].split('_')[0]
  
#   case_root = f'/pscratch/sd/w/whannah/2024-AQP-CESS/{case}'
#   os.chdir(f'{case_root}/case_scripts')

#   # os.system(f'grep \'id=\"ATM_NCPL\' ./*xml')
  
#   output1 = sp.check_output( f'grep \'id=\"ATM_NCPL\' ./*xml',shell=True,
#                             universal_newlines=True).split('\n')[0]
#   output2 = sp.check_output( f'grep se_tstep ./user_nl_eam',shell=True,
#                             universal_newlines=True).split('\n')[0]

#   print(f'{res:10}  {output1:60}  {output2:40}')

#   if res=='ne120pg2': print()

# print()
# exit()

#-------------------------------------------------------------------------------
case2_list = []
case1 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne30pg2_ne30pg2.NN_32.SSTP_0K'
case2_list.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne45pg2_ne45pg2.NN_64.SSTP_0K')
case2_list.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne60pg2_ne60pg2.NN_128.SSTP_0K')
case2_list.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne90pg2_ne90pg2.NN_256.SSTP_0K')
case2_list.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne120pg2_ne120pg2.NN_512.SSTP_0K')
# case1 = 'E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne30pg2_ne30pg2.NN_32.SSTP_0K'
# case2_list.append('E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne120pg2_ne120pg2.NN_512.SSTP_0K')
for case2 in case2_list:
  case1_root = f'/pscratch/sd/w/whannah/2024-AQP-CESS/{case1}'
  case2_root = f'/pscratch/sd/w/whannah/2024-AQP-CESS/{case2}'
  file_subpath = 'run/atm_in'
  # file_subpath = 'run/drv_in'

  case1_file = f'{case1_root}/{file_subpath}'
  case2_file = f'{case2_root}/{file_subpath}'

  cmd = f'diff {case1_file} {case2_file}'

  # print(cmd)
  # output = sp.check_output(cmd,shell=True,universal_newlines=True)
  # print(output)
  # print()

  print()
  os.system(cmd)
  print()
#-------------------------------------------------------------------------------
  
  
  

