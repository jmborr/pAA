work_dir = '/home/jbq/research/pAA/solvated_acrylate/flag_bound_water'
for sc in 'sf sfNN sfNY sfYN sfYY'.split():
    LoadSassena(Filename='{}/{}.h5'.format(work_dir,sc), OutputWorkspace=sc)