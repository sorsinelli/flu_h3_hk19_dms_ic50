# Define sites in each antigenic region
A_sites=['122', '124', '126', '130-133', '135', '137', '138', '140', '142-146', '150', '152', '168']
B_sites=['128-129', '155-160', '163', '165', '186-190', '192-194', '196-198']
C_sites=['44-48', '50', '51', '53-54', '273', '275', '276', '278-280', '294', '297', '299', '300', '304-305', '307-312']
D_sites=['96', '102', '103', '117', '121', '167', '170-177', '179', '182', '201', '203', '207-209', '212-219', '226-230', '238', '240', '242', '244', '246-248']
E_sites=['57', '59', '62-63', '67', '75', '78', '80-83', '86-88', '91', '92', '94', '109', '260-262', '265']

# Start a fresh PyMol session
cmd.reinitialize()

# Define the colors
cmd.set_color('color_B', [51, 51, 51])
cmd.set_color('color_A', [116, 116, 116])
cmd.set_color('color_E', [82, 82, 82])
cmd.set_color('color_D', [220, 220, 220])
cmd.set_color('color_C', [188, 188, 188])
cmd.set_color('RBD', [213, 94, 0])

# get H3 HA structure (4O5N)
cmd.load('../../data/PDBs/4o5n.pdb')

# color surface white
cmd.hide('everything', '4o5n')
cmd.show('surface', '4o5n')
cmd.set('surface_quality', 1)
cmd.set('ambient', 0.5)
cmd.color('white', 'all')

# Loop through chains A, C, and E and color by antigenic region
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

# set visualization parameters
cmd.set('ray_shadow', 'off')
cmd.set('light_count', 1)
cmd.set('depth_cue', 0)

# Set view for full protein
cmd.set_view([0.916666031, 0.107226826, -0.385001540,
            0.399538547,   -0.222682491,   0.889258325,
            0.009619311,   -0.968975425,   -0.246967122,
            0.000000000,    0.000000000, -417.736297607,
            50.343509674,  -29.065820694,  -53.074180603,
            329.346496582,  506.126098633,  -20.000000000]) 

# save full protein PNG
cmd.png('figure_2/antigenic_regions.png', ray=1)

# Repeat to save version zoomed in on site B
cmd.set_view(         [0.979569793,   -0.153052837,    0.130378887,
    -0.201045781,   -0.737559319,    0.644656718,
    -0.002502405,   -0.657701254,   -0.753268123,
     0.000000000,    0.000000000, -338.260589600,
    50.343509674,  -29.065820694,  -53.074180603,
  -1555.376098633, 2231.897705078,  -20.000000000] )

cmd.png('figure_4/antigenic_regions.png', ray=1)

# select RBD residues
cmd.hide('all')
cmd.select('RBS', 'resi 190 and chain A or resi 194 and chain A or resi 155 and chain A or resi 153 and chain A or resi 98 and chain A or resi 134-137 and chain A or resi 226 and chain A or resi 190 and chain E or resi 194 and chain E or resi 155 and chain E or resi 153 and chain E or resi 98 and chain E or resi 134-137 and chain E or resi 226 and chain E')

# save just the RBD, zoomed on site B
cmd.show('surface', 'RBS')
cmd.set('ray_opaque_background', 0)
cmd.png('figure_4/RBD_overlay.png', ray=1)

# save just the RBD outline, zoomed on site B
cmd.set('ray_trace_mode', 2)
cmd.set('ray_trace_color', 'RBD')
cmd.png('figure_4/RBD_outline.png', ray=1)

# get RBD overlay + outline for full protein view
cmd.set_view([0.916666031, 0.107226826, -0.385001540,
            0.399538547,   -0.222682491,   0.889258325,
            0.009619311,   -0.968975425,   -0.246967122,
            0.000000000,    0.000000000, -417.736297607,
            50.343509674,  -29.065820694,  -53.074180603,
            329.346496582,  506.126098633,  -20.000000000]) 

cmd.png('figure_2/RBD_outline.png', ray=1)

# set back to normal ray trace mode and save RBD surface
cmd.set('ray_trace_mode', 0)
cmd.png('figure_2/RBD_overlay.png', ray=1)