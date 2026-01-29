
def get_lutable(version,exp):
    '''
    ***** Note: please follow this example to add your case info. ******

    'case_shortname': [
    CTL_casename, CTL_startyear, CTL_endyear,
    P4K_casename, P4K_startyear, P4K_endyear,
    ]

    '''

    mmf_yr1 =  60
    mmf_yr2 = 120

    lu_table = {
    'mmf_2023cpl_2x1x':[
        'E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2',mmf_yr1,mmf_yr2,
        'E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.2xCO2',mmf_yr1,mmf_yr2,
    ],
    'mmf_2023cpl_4x1x':[
        'E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2',mmf_yr1,mmf_yr2,
        'E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.4xCO2',mmf_yr1,mmf_yr2,
    ],
    'mmf_2023cpl_4x2x':[
        'E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.2xCO2',mmf_yr1,mmf_yr2,
        'E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.4xCO2',mmf_yr1,mmf_yr2,
    ],
    'v1_coupled':[
    'None',None,None,
    'None',None,None,
    ],
    'v2_coupled':[
    'None',None,None,
    'None',None,None,
    ], 
    'v1_amip4K':[
    'None',None,None,
    'None',None,None,
    ], 
    'v1':[
    '20211208.F2010C5-CMIP6-LR.IC.ne30_oECv3.compy.1080',2,3,
    '20211210.F2010C5-CMIP6-LR.IC.p4Ka.ne30_oECv3.compy.1080',2,3,
    ],
    'v2test':[
    '20221013.v2.F2010-CICE.IC.back.clubb.g1.1.ne30pg2_EC.compy',2,3,
    '20221013.v2.F2010-CICE.IC.back.clubb.g1.1.p4Ka.ne30pg2_EC.compy',2,3,
    ],
    'v2test2':[
    '20221013.v2.F2010-CICE.IC.back.clubb.g1.1.ne30pg2_EC.compy2',3,4,
    '20221013.v2.F2010-CICE.IC.back.clubb.g1.1.p4Ka.ne30pg2_EC.compy2',3,4,
    ],
    }

    if exp == 'amip':
        return lu_table[version][0],lu_table[version][1],lu_table[version][2]
    else:
        return lu_table[version][3],lu_table[version][4],lu_table[version][5]


