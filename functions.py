from utils import *
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
    def __init__(self, radius, d, rho_0=0.17):
        super(Woods_Saxon, self).__init__()
        self.radius = radius
        self.d = d
        self.rho_0 = rho_0
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
        self.args['label'] = r'$<N_{\mathrm{part}}>$'
        self.args['sigma'] = 7
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
            ax.set_ylabel(r'$<N_{\mathrm{part}}>$')
            ax.legend()
            if type(args['title']) is str:
                ax.set_title('%s'%args['title'])
            if args['save']:
                fig.savefig(args['path'])

        