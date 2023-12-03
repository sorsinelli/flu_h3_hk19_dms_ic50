cohorts = [
	'hk19_2-5_years',
	'hk19_15-20_years',
	'hk19_40-45_years',
	'perth09_2-4_years',
	'perth09_30-34_years'
]

escape_color = [0.0, 0.447, 0.698]

for cohort in cohorts:
    cmd.reinitialize()
    infile=f'{cohort}_normalized_max.pdb'
    cmd.load(infile)
    cmd.hide('everything', f'{cohort}_normalized')
    cmd.show('surface', f'{cohort}_normalized')
    cmd.set_color('escape_color', escape_color)
    cmd.spectrum('b', 'white escape_color', minimum=0, maximum=1)
    
    # site B view:
    
    cmd.set_view(         [0.979569793,   -0.153052837,    0.130378887,
    -0.201045781,   -0.737559319,    0.644656718,
    -0.002502405,   -0.657701254,   -0.753268123,
     0.000000000,    0.000000000, -338.260589600,
    50.343509674,  -29.065820694,  -53.074180603,
  -1555.376098633, 2231.897705078,  -20.000000000] )
    
    # full protein view:    
#     cmd.set_view([0.916666031, 0.107226826, -0.385001540,
#                 0.399538547,   -0.222682491,   0.889258325,
#                 0.009619311,   -0.968975425,   -0.246967122,
#                 0.000000000,    0.000000000, -417.736297607,
#                 50.343509674,  -29.065820694,  -53.074180603,
#                 329.346496582,  506.126098633,  -20.000000000]) 

	# top down view:
	
#     cmd.set_view([0.966774940, -0.252064854, -0.042563122,
#                 -0.252220124,   -0.967667818,   0.001727089,
#                 -0.041621406,   0.009066927,   -0.999092460,
#                 0.000000000,    0.000000000, -417.736297607,
#                 50.343509674,  -29.065820694,  -53.074180603,
#                 329.346496582,  506.126098633,  -20.000000000]) 



        
    outfile = cohort + '_max_siteB.png'
    
    cmd.set('ray_shadow', 'off')
    cmd.set('light_count', 1)
    cmd.set('depth_cue', 0)
    cmd.set('ambient', 0.5)
    cmd.png(outfile, ray=1)


# sera_max = {
# #             'AUSAB-11': 18.2,
# #             '1C04-5G04_cocktail': 4.24,
#             'AUSAB-07': 41.8,
#             'AUSAB-13': 17,
#            }


# +
# #     'AUSAB-07_sum_epitope_1': 52.1,
# #     'AUSAB-07_sum_epitope_2': 5.3,
#     'AUSAB-11_sum_epitope_1': 17.3,
#     'AUSAB-11_sum_epitope_2': 7.4,
#     'AUSAB-13_sum_epitope_1': 18.9,
#     'AUSAB-13_sum_epitope_2': 3,
#     'AUSAB-05_sum_epitope_1': 1.23
# -

# epitope_colors = {'1': [0.0, 0.447, 0.698],
# #                   '2': [0.8, 0.475, 0.655],
# #                   '3': [0.298, 0.208, 0.286]
#                  }
# 
# negative_color = [0.902, 0.624, 0.0]

# +
# cmd.set_color('negative_color', negative_color)

# for serum in sera_max:
#     for epitope, max_color in epitope_colors.items():
#         cmd.reinitialize()
#         infile=f'{serum}_sum_epitope_{epitope}.pdb'
#         cmd.load(infile)
# #         cmd.bg_color('white')
#         cmd.hide('everything', f'{serum}_sum_epitope_{epitope}')
#         cmd.show('surface', f'{serum}_sum_epitope_{epitope}')
#         cmd.set_color('max_color', epitope_colors[epitope])
#         cmd.set_color('negative_color', negative_color)
#         cmd.spectrum('b', 'negative_color white max_color', minimum=-(sera_max[serum]), maximum=sera_max[serum])
#         
#         cmd.fetch('5thf')
#         
#         cmd.align(f'{serum}_sum_epitope_{epitope}',  '/5thf//A+B')
#         cmd.copy('obj01', f'{serum}_sum_epitope_{epitope}')
#         cmd.align('obj01',  '/5thf//C+D')
#         cmd.copy('obj02', f'{serum}_sum_epitope_{epitope}')
#         cmd.align('obj02',  '/5thf//E+F')
# #         cmd.show('surface', 'obj01')
# #         cmd.show('surface', 'obj02')
#         
#         cmd.hide('everything', '5thf')
#         
#         cmd.set_view([0.837208807, 0.353551418, 0.417231888,
#                       -0.483371109,    0.835233569,    0.262171030,
#                       -0.255793542,   -0.421170801,    0.870161235,
#                       0.000000000,    0.000000000, -451.576568604,
#                       70.882568359,   23.363365173,   18.654100418,
#                       -71602.406250000, 72505.570312500,  -20.000000000])
#         
#         outfile = os.path.splitext(infile)[0] + '.png'
#         cmd.png(outfile, ray=1)


# +
# set_view (\
#      0.837208807,    0.353551418,    0.417231888,\
#     -0.483371109,    0.835233569,    0.262171030,\
#     -0.255793542,   -0.421170801,    0.870161235,\
#      0.000000000,    0.000000000, -451.576568604,\
#     70.882568359,   23.363365173,   18.654100418,\
#   -71602.406250000, 72505.570312500,  -20.000000000 )