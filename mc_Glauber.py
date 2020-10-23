from functions import Nbp_mc, probability

#initialize
Au = {'name':'Au', 'radius':6.38, 'd':0.535, 'A':197, 'sigma':4.2}
'''
#draw Nb and Np
inst_Nbp_mc = Nbp_mc(Au, Au, 1000)
args_nb = inst_Nbp_mc.args
args_nb['title'] = 'Au + Au at 200GeV'
args_nb['save'] = True
args_nb['path'] = './AuAu_Nb_mc.png'
args_nb['b_bins'] = 15
args_nb['rho_bins'] = 20
args_nb['theta_bins'] = 20
inst_Nbp_mc.plot_func(args_nb)
inst_Nbp_mc.save_data()
'''
#draw Pb PNb and PNp
inst_prob_mc = probability(15)
args_prob = inst_prob_mc.args
args_prob['save'] = True
args_prob['path'] = './AuAu_P_mc.png'
inst_prob_mc.event_generator(100000)
inst_prob_mc.plot_func(args_prob)