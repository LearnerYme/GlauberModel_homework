from functions import Woods_Saxon, thickness, overlap, Nb

#initialize
Au = {'name':'Au', 'radius':6.38, 'd': 0.535, 'A':197, 'sigma':4.2}
Cu = {'name':'Cu', 'radius':4.20641, 'A':7, 'd': 0.5977, 'sigma':'Unknown'}
Pb = {'name':'Pb', 'radius':6.62, 'A':7, 'd': 0.546, 'sigma':'Unknown'}

#draw Woods Saxon and thickness
def job1(nuclei):
    #problem 1
    #Woods Saxon
    inst_woods_saxon = Woods_Saxon(nuclei)
    args_ws = inst_woods_saxon.args
    args_ws['title'] = nuclei['name']
    args_ws['save'] = True
    args_ws['path'] = './%s_Woods_Saxon.png'%nuclei['name']
    inst_woods_saxon.plot_func(args_ws)
    #problem 2
    #Thickness
    inst_thickness = thickness(inst_woods_saxon)
    args_th = inst_thickness.args
    args_th['title'] = nuclei['name']
    args_th['save'] = True
    args_th['path'] = './%s_thickness.png'%nuclei['name']
    inst_thickness.plot_func(args_th)
    return inst_thickness

#for Au, draw overlap and N collision
def job2(nuclei):
    #problem 1 and 2
    inst_thickness = job1(nuclei)
    #problem 3
    #overlap
    inst_overlap = overlap(inst_thickness)
    args_ov = inst_overlap.args
    args_ov['title'] = nuclei['name']
    args_ov['save'] = True
    args_ov['path'] = './%s_overlap.png'%nuclei['name']
    args_ov['second_tget_bins'] = 15
    inst_overlap.plot_func(args_ov)
    #N collision
    inst_Nb = Nb(inst_overlap)
    args_nb = inst_Nb.args
    args_nb['title'] = nuclei['name']
    args_nb['save'] = True
    args_nb['path'] = './%s_Nb.png'%nuclei['name']
    args_nb['sigma'] = nuclei['sigma']
    inst_Nb.plot_func(args_nb)
    return

#apply above functions on exact nucleus
for nuclei in [Cu, Pb]:
    job1(nuclei)
job2(Au)