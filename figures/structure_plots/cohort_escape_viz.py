# define cohort naming
cohorts = [
	'hk19_2-5_years',
	'hk19_15-20_years',
	'hk19_40-45_years',
	'perth09_2-4_years',
	'perth09_30-34_years'
]

# define color corresponding to escape
escape_color = [0.0, 0.447, 0.698]

# Loop through cohorts and save PNGs of the resulting structures
# Saves both a view of the full protein and a view focused on antigenic region B
for cohort in cohorts:
    cmd.reinitialize()
    infile=f'escape_score_pdbs/{cohort}_normalized_max.pdb'
    cmd.load(infile)
    cmd.hide('everything', f'{cohort}_normalized')
    cmd.show('surface', f'{cohort}_normalized')
    cmd.set_color('escape_color', escape_color)
    
    # recolor by b-factors, with white being neutral
    cmd.spectrum('b', 'white escape_color', minimum=0, maximum=1)
    
    # adjust visualization of the protein
    cmd.set('ray_shadow', 'off')
    cmd.set('light_count', 1)
    cmd.set('depth_cue', 0)
    cmd.set('ambient', 0.5)
    
    
    # set full protein view:    
    cmd.set_view([0.916666031, 0.107226826, -0.385001540,
                0.399538547,   -0.222682491,   0.889258325,
                0.009619311,   -0.968975425,   -0.246967122,
                0.000000000,    0.000000000, -417.736297607,
                50.343509674,  -29.065820694,  -53.074180603,
                329.346496582,  506.126098633,  -20.000000000]) 
    
    # save full protein figure
    outfile = 'figure_5/' + cohort + '_max.png'
    cmd.png(outfile, ray=1)
    
    
    # set site B view:
    cmd.set_view(         [0.979569793,   -0.153052837,    0.130378887,
    -0.201045781,   -0.737559319,    0.644656718,
    -0.002502405,   -0.657701254,   -0.753268123,
     0.000000000,    0.000000000, -338.260589600,
    50.343509674,  -29.065820694,  -53.074180603,
  -1555.376098633, 2231.897705078,  -20.000000000] )
     
    # save site B figure
    outfile = 'figure_5/' + cohort + '_max_siteB.png'
    cmd.png(outfile, ray=1)