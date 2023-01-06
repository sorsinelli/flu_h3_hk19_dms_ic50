sera_max = {
#     'AUSAB-07': 52.1,
            'AUSAB-11': 17.3,
            'AUSAB-13': 18.9,
            'AUSAB-05': 1.23,
           }


# +
# #     'AUSAB-07_sum_epitope_1': 52.1,
# #     'AUSAB-07_sum_epitope_2': 5.3,
#     'AUSAB-11_sum_epitope_1': 17.3,
#     'AUSAB-11_sum_epitope_2': 7.4,
#     'AUSAB-13_sum_epitope_1': 18.9,
#     'AUSAB-13_sum_epitope_2': 3,
#     'AUSAB-05_sum_epitope_1': 1.23
# -

epitope_colors = {'1': [0.0, 0.447, 0.698],
                  '2': [0.8, 0.475, 0.655]}

negative_color = [0.902, 0.624, 0.0]

# +
# cmd.set_color('negative_color', negative_color)

for serum in sera_max:
    for epitope, max_color in epitope_colors.items():
        cmd.reinitialize()
        infile=f'{serum}_sum_epitope_{epitope}.pdb'
        cmd.load(infile)
#         cmd.bg_color('white')
        cmd.hide('everything', f'{serum}_sum_epitope_{epitope}')
        cmd.show('surface', f'{serum}_sum_epitope_{epitope}')
        cmd.set_color('max_color', epitope_colors[epitope])
        cmd.set_color('negative_color', negative_color)
        cmd.spectrum('b', 'negative_color white max_color', minimum=-(sera_max[serum]), maximum=sera_max[serum])
        
        cmd.fetch('5thf')
        
        cmd.align(f'{serum}_sum_epitope_{epitope}',  '/5thf//A+B')
        cmd.copy('obj01', f'{serum}_sum_epitope_{epitope}')
        cmd.align('obj01',  '/5thf//C+D')
        cmd.copy('obj02', f'{serum}_sum_epitope_{epitope}')
        cmd.align('obj02',  '/5thf//E+F')
#         cmd.show('surface', 'obj01')
#         cmd.show('surface', 'obj02')
        
        cmd.hide('everything', '5thf')
        
        cmd.set_view([0.837208807, 0.353551418, 0.417231888,
                      -0.483371109,    0.835233569,    0.262171030,
                      -0.255793542,   -0.421170801,    0.870161235,
                      0.000000000,    0.000000000, -451.576568604,
                      70.882568359,   23.363365173,   18.654100418,
                      -71602.406250000, 72505.570312500,  -20.000000000])
        
        outfile = os.path.splitext(infile)[0] + '.png'
        cmd.png(outfile, ray=1)


# +
# set_view (\
#      0.837208807,    0.353551418,    0.417231888,\
#     -0.483371109,    0.835233569,    0.262171030,\
#     -0.255793542,   -0.421170801,    0.870161235,\
#      0.000000000,    0.000000000, -451.576568604,\
#     70.882568359,   23.363365173,   18.654100418,\
#   -71602.406250000, 72505.570312500,  -20.000000000 )
