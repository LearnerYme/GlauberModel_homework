from utils import rect_vector, polar_vector, find_bin, integ, rdg, nucleon_generator, bsq
import numpy as np
import matplotlib.pyplot as plt

class custom_function():
    def __init__(self):
        self.args = {
            'title':None,
            'save':False,
            'path':'./res.png'
        }
        return

    def function(self, *args, **kwargs):
        pass

    def get_value(self, *args):
        self.value = self.function(*args)
        return

    def plot_func(self, *args, **kwargs):
        pass

class Woods_Saxon(custom_function):
    def __init__(self, nuclei):
        super(Woods_Saxon, self).__init__()
        self.radius = nuclei['radius']
        self.d = nuclei['d']
        self.A = nuclei['A']
        self.rho_0 = 0.17
        self.args['x'] = np.linspace(0, 15, 200)
        self.args['label'] = 'Woods Saxon'
        return

    def function(self, r):
        return self.rho_0/(1 + np.exp((r - self.radius)/self.d))

    def plot_func(self, args=None):
        if args is None:
            args = self.args
        self.x = args['x']
        self.get_value(self.x)
        with plt.style.context(['ieee', 'science', 'no-latex', 'grid']):
            fig, ax = plt.subplots()
            ax.plot(self.x, self.value, label=args['label'])
            ax.set_xlabel(r'$r\ \mathrm{fm}$')
            ax.set_ylabel(r'$\rho_A(r)$')
            ax.legend()
            if type(args['title']) is str:
                ax.set_title('%s Radius=%.4f, d=%.4f'%(args['title'], self.radius, self.d))
            if args['save']:
                fig.savefig(args['path'])
        return

class thickness(custom_function):
    def __init__(self, func:custom_function):
        super(thickness, self).__init__()
        self.woods_saxon = func
        self.args['var'] = np.linspace(0, 15, 200)
        self.args['func_value'] = self.woods_saxon.value
        self.args['tget_var_range'] = [0, 15]
        self.args['integ_var_range'] = [0, 15]
        self.args['tget_bins'] = 200
        self.args['integ_bins'] = 200
        self.args['label'] = 'Thickness'
        return

    @staticmethod
    def param_func(s, z, fargs=None):
        return np.sqrt(s**2 + z**2)

    def function(self, args=None):
        if args is None:
            args = self.args
        var = args['var']
        func_value = args['func_value']
        s_var_range = args['tget_var_range']
        z_var_range = args['integ_var_range']
        s_bins = args['tget_bins']
        z_bins = args['integ_bins']
        return integ.integ1d(var, func_value, s_var_range, z_var_range, s_bins, z_bins, thickness.param_func)

    def plot_func(self, args=None):
        if args is None:
            args = self.args
        self.get_value(args)
        self.x = self.value[0]
        self.value = self.value[1]
        with plt.style.context(['ieee', 'science', 'no-latex', 'grid']):
            fig, ax = plt.subplots()
            ax.plot(self.x, self.value, label=args['label'])
            ax.set_xlabel(r'$\vec{r}_\perp\ \mathrm{fm}$')
            ax.set_ylabel(r'$T_A(\vec{r}_\perp)\ \mathrm{fm}^{-2}$')
            ax.legend()
            if type(args['title']) is str:
                ax.set_title('%s'%args['title'])
            if args['save']:
                fig.savefig(args['path'])     

class overlap(custom_function):
    def __init__(self, func:custom_function):
        super(overlap, self).__init__()
        self.thickness = func
        self.args['var'] = self.thickness.x
        self.args['func_value'] = self.thickness.value
        self.args['first_tget_var_range'] = [0, 30]
        self.args['second_tget_var_range'] = [0, 15]
        self.args['first_integ_var_range'] = [0, 2*np.pi]
        self.args['second_integ_var_range'] = [0, 30]
        self.args['second_tget_bins'] = 200
        self.args['first_integ_bins'] = 100
        self.args['second_integ_bins'] = 100
        self.args['label'] = 'Overlap'
        return

    def function(self, args=None):
        if args is None:
            args = self.args
        var = args['var']
        func_value = args['func_value']
        b_var_range = args['second_tget_var_range']
        theta_var_range = args['first_integ_var_range']
        rho_var_range = args['second_integ_var_range']
        b_bins = args['second_tget_bins']
        theta_bins = args['first_integ_bins']
        rho_bins = args['second_integ_bins']
        #hard to use utils.integ for this vector integral...
        #directly integral in this method
        b_arr = np.linspace(b_var_range[0], b_var_range[1], b_bins)
        s_arr = np.linspace(rho_var_range[0], rho_var_range[1], rho_bins)
        theta_arr = np.linspace(theta_var_range[0], theta_var_range[1], theta_bins)
        s_stride = (rho_var_range[1] - rho_var_range[0])/rho_bins
        theta_stride = (theta_var_range[1] - theta_var_range[0])/theta_bins
        O_res = np.zeros_like(b_arr)
        for bidx, bitem in enumerate(b_arr, 0):
            b_vec = rect_vector(bitem, 0)
            res = 0
            print('Computing b = %.4f...'%bitem)
            for sidx, sitem in enumerate(s_arr, 0):
                print('\r Computing s = %.4f...'%sitem, end='')
                for thetaidx, thetaitem in enumerate(theta_arr, 0):
                    s_vec = polar_vector(sitem, thetaitem).rect_vec
                    sA_scalar = sitem
                    sB_scalar = (s_vec + b_vec).get_length
                    tA = func_value[find_bin.find_bin1d(sA_scalar, var)]
                    tB = func_value[find_bin.find_bin1d(sB_scalar, var)]
                    phase_interval = sitem*s_stride*theta_stride
                    res += tA*tB*phase_interval
            print('')
            O_res[bidx] = res
        return b_arr, O_res
        
    def plot_func(self, args=None):
        if args is None:
            args = self.args
        self.get_value(args)
        self.x = self.value[0]
        self.value = self.value[1]
        with plt.style.context(['ieee', 'science', 'no-latex', 'grid']):
            fig, ax = plt.subplots()
            ax.plot(self.x, self.value, label=args['label'])
            ax.set_xlabel(r'$b\ \mathrm{fm}$')
            ax.set_ylabel(r'$T_{AA}(\vec{r}_\perp)\ \mathrm{fm}^{-2}$')
            ax.legend()
            if type(args['title']) is str:
                ax.set_title('%s'%args['title'])
            if args['save']:
                fig.savefig(args['path'])

class Nb(custom_function):
    def __init__(self, func:custom_function):
        super(Nb, self).__init__()
        self.overlap = func
        self.args['label'] = r'$<N_{\mathrm{coll}}>$'
        self.args['sigma'] = 0
        return

    def function(self, args=None):
        if args is None:
            args = self.args
        return self.overlap.value*args['sigma']

    def plot_func(self, args=None):
        if args is None:
            args = self.args
        self.get_value(args)
        self.x = self.overlap.x
        with plt.style.context(['ieee', 'science', 'no-latex', 'grid']):
            fig, ax = plt.subplots()
            ax.plot(self.x, self.value, label=args['label'])
            ax.set_xlabel(r'$b\ \mathrm{fm}$')
            ax.set_ylabel(r'$<N_{\mathrm{coll}}>$')
            ax.legend()
            if type(args['title']) is str:
                ax.set_title('%s'%args['title'])
            if args['save']:
                fig.savefig(args['path'])

class thickness_mc(custom_function):
    def __init__(self, nuclei):
        super(thickness_mc, self).__init__()
        self.ws = Woods_Saxon(nuclei)
        self.nucleon = nucleon_generator(self.ws)
        self.s, self.theta = self.nucleon.get_particle()#note: here I rename phi -> theta
        self.args['rho_range'] = [0, 15]
        self.args['theta_range'] = [0, 2*np.pi]
        self.args['rho_bins'] = 200
        self.args['theta_bins'] = 100
        return

    def function(self, args=None):
        if args is None:
            args = self.args
        rho_bins = args['rho_bins']
        theta_bins = args['theta_bins']
        rho_arr = np.linspace(args['rho_range'][0], args['rho_range'][1], rho_bins)
        theta_arr = np.linspace(args['theta_range'][0], args['theta_range'][1], theta_bins)
        T_res = np.zeros((rho_bins, theta_bins))
        for p_idx in range(self.ws.A):
            rho_idx, theta_idx = find_bin.find_bin2d(self.s[p_idx], self.theta[p_idx], rho_arr, theta_arr)
            T_res[rho_idx, theta_idx] += 1
        return T_res, rho_arr, theta_arr

class Nbp_mc(custom_function):
    def __init__(self, nuclei1, nuclei2, entries=1000):
        super(Nbp_mc, self).__init__()
        self.nuclei1 = nuclei1
        self.nuclei2 = nuclei2
        self.entries = entries
        self.nucleon_tk_1 = None
        self.nucleon_tk_2 = None
        self.args['sigma'] = 4.2#indeed, the sigma of nuclei1 and 2 is not really exact
        self.args['b_range'] = [0, 15]
        self.args['b_bins'] = 200
        self.args['rho_range'] = [0, 15]
        self.args['theta_range'] = [0, 2*np.pi]
        self.args['rho_bins'] = 200
        self.args['theta_bins'] = 100
        self.args['label1'] = r'$<N_{\mathrm{coll}}>$'
        self.args['label2'] = r'$<N_{\mathrm{part}}>$'
        return

    def nucleon_thickness_init(self, args=None):
        if args is None:
            args = self.args
        self.nucleon_tk_1 = thickness_mc(self.nuclei1)
        self.nucleon_tk_2 = thickness_mc(self.nuclei2)
        self.T_res_1, self.rho_arr, self.theta_arr = self.nucleon_tk_1.function(args)
        self.T_res_2, _, _ = self.nucleon_tk_2.function(args)
        self.s1 = self.nucleon_tk_1.s
        self.s2 = self.nucleon_tk_2.s
        self.theta1 = self.nucleon_tk_1.theta
        self.theta2 = self.nucleon_tk_2.theta
        return

    def function(self, args=None):
        r'''
        Note: as O = sum of i: {Ti rho_i delta_rho delta_theta}
        Here Ti = thickness i / rho_i, so we apply:
        O = sum of i: {thickness i delta_tho delta_theta}
        '''
        if args is None:
            args = self.args
        self.nucleon_thickness_init(args)
        b_arr = np.linspace(args['b_range'][0], args['b_range'][1], args['b_bins'])
        Nb_arr = np.zeros_like(b_arr)
        Np_arr = np.zeros_like(b_arr)
        s1 = self.nucleon_tk_1.s
        theta1 = self.nucleon_tk_1.theta
        s2 = self.nucleon_tk_2.s
        theta2 = self.nucleon_tk_2.theta
        rho_stride = self.rho_arr[1] - self.rho_arr[0]
        theta_stride = self.theta_arr[1] - self.theta_arr[0]
        for b_idx, b_item in enumerate(b_arr, 0):
            b_vec = rect_vector(b_item, 0)
            for idx in range(self.nuclei1['A']):
                sA = polar_vector(s1[idx], theta1[idx])
                sB = polar_vector(s2[idx], theta2[idx])
                sA_rho_bin = find_bin.find_bin1d(sA.rho, self.rho_arr)
                sA_theta_bin = find_bin.find_bin1d(sA.theta, self.theta_arr)
                sB_rho_bin = find_bin.find_bin1d(sB.rho, self.rho_arr)
                sB_theta_bin = find_bin.find_bin1d(sB.theta, self.theta_arr)
                T1 = self.T_res_1[sA_rho_bin, sA_theta_bin]
                T2_inverse = self.T_res_2[sB_rho_bin, sB_theta_bin]
                sApb = (sA.rect_vec + b_vec).polar_vec
                sApb_rho_bin = find_bin.find_bin1d(sApb.rho, self.rho_arr)
                sApb_theta_bin = find_bin.find_bin1d(sApb.theta, self.theta_arr)
                sBpb = (sB.rect_vec + b_vec).polar_vec
                sBpb_rho_bin = find_bin.find_bin1d(sBpb.rho, self.rho_arr)
                sBpb_theta_bin = find_bin.find_bin1d(sBpb.theta, self.theta_arr)
                T1_inverse = self.T_res_1[sBpb_rho_bin, sBpb_theta_bin]
                T2 = self.T_res_2[sApb_rho_bin, sApb_theta_bin]
                Nb_arr[b_idx] += T1*T2*self.args['sigma']*rho_stride*theta_stride
                Np_A = (1 - (1 - T2/self.nuclei2['A']*self.args['sigma'])**self.nuclei2['A'])*T1*rho_stride*theta_stride
                Np_B = (1 - (1 - T1_inverse/self.nuclei1['A']*self.args['sigma'])**self.nuclei1['A'])*T2_inverse*rho_stride*theta_stride
                Np_arr[b_idx] += Np_A + Np_B
        return Nb_arr, Np_arr

    def plot_func(self, args=None):
        if args is None:
            args = self.args
        self.x = np.linspace(args['b_range'][0], args['b_range'][1], args['b_bins'])
        Nb_res = np.zeros_like(self.x)
        Np_res = np.zeros_like(self.x)
        for i in range(self.entries):
            print('Entry: %d of %d...'%(i+1, self.entries))
            self.get_value(args)
            Nb_tmp, Np_tmp = self.value
            Nb_res += Nb_tmp
            Np_res += Np_tmp
        self.Nb_res = Nb_res / self.entries
        self.Np_res = Np_res / self.entries
        with plt.style.context(['ieee', 'science', 'no-latex', 'grid']):
            fig, ax = plt.subplots()
            ax.plot(self.x, self.Nb_res, 'r--', label=args['label1'])
            ax.plot(self.x, self.Np_res, 'b:', label=args['label2'])
            ax.set_xlabel(r'$b\ \mathrm{fm}$')
            ax.set_ylabel(r'$<N_{\mathrm{coll}}>, <N_{\mathrm{part}}>$')
            ax.legend()
            if type(args['title']) is str:
                ax.set_title('%s'%args['title'])
            if args['save']:
                fig.savefig(args['path'])
        
    def save_data(self, path='./NbNp_data.csv'):
        b = self.x.reshape(-1, 1)
        Nb = self.Nb_res.reshape(-1, 1)
        Np = self.Np_res.reshape(-1, 1)
        res = np.concatenate((b, Nb, Np), 1)
        np.savetxt(path, res, '%.3f', ',')
        print('Data saved in %s!'%path)
        return
        
class probability(custom_function):
    def __init__(self, bmax, datapath='./NbNp_data.csv'):
        super(probability, self).__init__()
        data = np.loadtxt(datapath, delimiter=',')
        self.b_arr = data[:, 0]
        self.Nb_arr = data[:, 1]
        self.Np_arr = data[:, 2]
        self.bsq = bsq(bmax)
        self.rdg = rdg(self.bsq.function, [0, bmax], [0, bmax**2])
        return

    def event_generator(self, entries=1000):
        b = []
        Nb = []
        Np = []
        for i in range(entries):
            print('Entry: %d of %d'%(i+1, entries))
            b_tmp = self.rdg.get_random()
            b_tmp_bin = find_bin.find_bin1d(b_tmp, self.b_arr)
            b.append(b_tmp)
            Nb.append(self.Nb_arr[b_tmp_bin])
            Np.append(self.Np_arr[b_tmp_bin])
        self.b = np.array(b)
        self.Nb = np.array(Nb)
        self.Np = np.array(Np)
        self.args['label1'] = 'b'
        self.args['label2'] = r'$<N_{\mathrm{coll}}>$'
        self.args['label3'] = r'$<N_{\mathrm{part}}>$'
        return

    def plot_func(self, args=None):
        if args is None:
            args = self.args
        with plt.style.context(['ieee', 'science', 'no-latex', 'grid']):
            fig, ax = plt.subplots(2, 2)
            fig.set_size_inches(8, 8)
            ax[0, 0].hist(self.b, label=args['label1'], density=True)
            ax[0, 0].set_xlabel(r'$b\ \mathrm{fm}$')
            ax[0, 0].set_ylabel(r'$P(b)$')
            ax[0, 0].legend()
            ax[0, 1].hist(self.Nb, label=args['label2'])
            ax[0, 1].set_xlabel(r'$<N_{\mathrm{coll}}>$')
            ax[0, 1].set_ylabel(r'$P(<N_{\mathrm{coll}}>)$')
            ax[0, 1].legend()
            ax[1, 0].hist(self.Np, label=args['label3'])
            ax[1, 0].set_xlabel(r'$<N_{\mathrm{part}}>$')
            ax[1, 0].set_ylabel(r'$P(<N_{\mathrm{part}}>)$')
            ax[1, 0].legend()
            fig.tight_layout()
            if args['save']:
                fig.savefig(args['path'])


