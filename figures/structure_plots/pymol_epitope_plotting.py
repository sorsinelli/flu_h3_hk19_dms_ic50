# Define your site lists
A_sites=['122', '124', '126', '130-133', '135', '137', '138', '140', '142-146', '150', '152', '168']
B_sites=['128-129', '155-160', '163', '165', '186-190', '192-194', '196-198']
C_sites=['44-48', '50', '51', '53-54', '273', '275', '276', '278-280', '294', '297', '299', '300', '304-305', '307-312']
D_sites=['96', '102', '103', '117', '121', '167', '170-177', '179', '182', '201', '203', '207-209', '212-219', '226-230', '238', '240', '242', '244', '246-248']
E_sites=['57', '59', '62-63', '67', '75', '78', '80-83', '86-88', '91', '92', '94', '109', '260-262', '265']
# A_sites=['122-146']
# B_sites=['155-169', '186-198']
# C_sites=['44-54', '273-280']
# D_sites=['201-219']
# E_sites=['62-65', '78-94', '260-265']

# Define the colors
cmd.set_color('color_B', [51, 51, 51])
cmd.set_color('color_A', [116, 116, 116])
cmd.set_color('color_E', [82, 82, 82])
cmd.set_color('color_D', [220, 220, 220])
cmd.set_color('color_C', [188, 188, 188])
# cmd.set_color('color_A', [0, 0, 0])
# cmd.set_color('color_B', [86, 86, 86])
# cmd.set_color('color_C', [32, 32, 32])
# cmd.set_color('color_D', [58, 58, 58])
# cmd.set_color('color_E', [116, 116, 116])
# cmd.set_color('color_C', [252, 205, 229])
# cmd.set_color('color_B', [188, 128, 189])
# cmd.set_color('color_A', [204, 235, 197])
# cmd.set_color('color_D', [153, 153, 153])
# cmd.set_color('color_E', [255, 237, 111])
# cmd.set_color('color_C', [255, 157, 166])
# cmd.set_color('color_B', [178, 121, 162])
# cmd.set_color('color_A', [157, 117, 93])
# cmd.set_color('color_D', [153, 153, 153])
# cmd.set_color('color_E', [238, 202, 59])
cmd.set_color('RBD', [213, 94, 0])

structure='4o5n'
# cmd.load(f'{infile}.pdb')

cmd.hide('everything', structure)
cmd.show('surface', structure)
cmd.set('surface_quality', 1)
cmd.set('ambient', 0.5)
cmd.color('white', 'all')

# Loop through chains A, C, and E
for chain in ['A', 'C', 'E']:
    for resi in A_sites:
    	cmd.color('color_A', f'resi {resi} and chain {chain}')
    for resi in B_sites:
    	cmd.color('color_B', f'resi {resi} and chain {chain}')
    for resi in C_sites:
    	cmd.color('color_C', f'resi {resi} and chain {chain}')
    for resi in D_sites:
    	cmd.color('color_D', f'resi {resi} and chain {chain}')
    for resi in E_sites:
    	cmd.color('color_E', f'resi {resi} and chain {chain}')
#     cmd.color('color_A', f'resi {A_sites[0]} and chain {chain}')
#     cmd.color('color_B', f'resi {B_sites[0]} and chain {chain} or resi {B_sites[1]} and chain {chain}')
#     cmd.color('color_C', f'resi {C_sites[0]} and chain {chain} or resi {C_sites[1]} and chain {chain}')
# # 	cmd.color('color_D', f'resi {D_sites[0]} and chain {chain}')
#     cmd.color('color_E', f'resi {E_sites[0]} and chain {chain} or resi {E_sites[1]} and chain {chain} or resi {E_sites[2]} and chain {chain}')


cmd.set_view(     [0.916666031,    0.107226826,   -0.385001540,
     0.399538547,   -0.222682491,    0.889258325,
     0.009619311,   -0.968975425,   -0.246967122,
     0.000000000,    0.000000000, -417.736297607,
    50.343509674,  -29.065820694,  -53.074180603,
   329.346496582,  506.126098633,  -20.000000000] )

cmd.set('ray_shadow', 'off')
cmd.set('light_count', 1)
cmd.set('depth_cue', 0)

cmd.png('base.png', ray=1)

cmd.hide('all')
cmd.select('RBS', 'resi 190 and chain A or resi 194 and chain A or resi 155 and chain A or resi 153 and chain A or resi 98 and chain A or resi 134-137 and chain A or resi 226 and chain A or resi 190 and chain E or resi 194 and chain E or resi 155 and chain E or resi 153 and chain E or resi 98 and chain E or resi 134-137 and chain E or resi 226 and chain E')
# or resi 190 and chain C or resi 194 and chain C or resi 155 and chain C or resi 153 and chain C or resi 98 and chain C or resi 134-137 and chain C or resi 226 and chain C 
# cmd.select('RBS', 'resi 190 and chain A or resi 194 and chain A or resi 195 and chain A or resi 183 and chain A or resi 153 and chain A or resi 228 and chain A or resi 225-226 and chain A or resi 190 and chain C or resi 194 and chain C or resi 195 and chain C or resi 183 and chain C or resi 153 and chain C or resi 228 and chain C or resi 225-226 and chain C or resi 190 and chain E or resi 194 and chain E or resi 195 and chain E or resi 183 and chain E or resi 153 and chain E or resi 228 and chain E or resi 225-226 and chain E')
# cmd.select('RBS', 'resi 201-202 and chain A or resi 193 and chain A or resi 240 and chain A or resi 191 and chain A or resi 158 and chain A or resi 237-238 and chain A or resi 95 and chain A or resi 139-141 and chain A or resi 201-202 and chain C or resi 193 and chain C or resi 240 and chain C or resi 191 and chain C or resi 158 and chain C or resi 237-238 and chain C or resi 95 and chain C or resi 139-141 and chain C or resi 201-202 and chain E or resi 193 and chain E or resi 240 and chain E or resi 191 and chain E or resi 158 and chain E or resi 237-238 and chain E or resi 95 and chain E or resi 139-141 and chain E')

cmd.show('surface', 'RBS')
cmd.set('ray_trace_mode', 2)
cmd.set('ray_trace_color', 'RBD')
cmd.set('ray_opaque_background', 0)

cmd.png('overlay.png', ray=1)