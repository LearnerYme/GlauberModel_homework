import numpy as np
import numpy.random as rd

class rect_vector():
    def __init__(self, x, y):
        self.x = x
        self.y = y
        return

    def __add__(self, other):
        return rect_vector(self.x + other.x, self.y + other.y)

    def __mul__(self, factor):
        return rect_vector(self.x*factor, self.y*factor)

    @property
    def get_length(self):
        return (self.x**2 + self.y**2)**0.5

    @property
    def polar_vec(self):
        return polar_vector(self.get_length, np.arctan2(self.y, self.x))

class polar_vector():
    def __init__(self, rho, theta):
        self.rho = rho
        self.theta = theta
        return

    @property
    def rect_vec(self):
        return rect_vector(self.rho*np.cos(self.theta), self.rho*np.sin(self.theta))

class find_bin():
    def __init__(self):
        r'''
        This class find location of a variable in a array-like sequence, or a location of a pixel in a matrix-like object.
        '''
        pass

    @staticmethod
    def find_bin1d(x, seq):
        if x <= seq[0]:
            return 0
        else:
            return (x>seq).sum() - 1

    @staticmethod
    def find_bin2d(x, y, seqx, seqy):
        bx = find_bin.find_bin1d(x, seqx)
        by = find_bin.find_bin1d(y, seqy)
        return [bx, by]

class integ():
    def __init__(self):
        r'''
        This class requires a list ( or array) of TBD x, and function value y,
        '''
        pass

    @staticmethod
    def integ1d_mc(function, x_range, f_range, entries=1000):
        r'''
        Not finished yet, do not use this method.
        '''
        x = rd.uniform(x_range[0], x_range[1], entries)
        f = rd.uniform(f_range[0], f_range[1], entries)
        fx = function(x)
        res = ((fx - f) >= 0).sum() * 1.0 / entries
        return res

    #old: def integ1d(y_range, y_func, x_range, s_range, xbins, sbins, relation):
    @staticmethod
    def integ1d(var, func_value, tget_var_range, integ_var_range, tget_bins, integ_bins, param_func, fargs=None):
        r'''
        This method will do a 1d integral over a pre-known-valued function, i.e. the discrete function_value-variable phase point.
        Here shows the meaning of arguments:
        func_value = func(var)
        var = param_func(tget_var, integ_var)
        func(tget_var) = \integral {integ_var_low, integ_var_high, func(param_func(tget_var, integ_var))d(integ_var)}
        '''
        tget_arr = np.linspace(tget_var_range[0], tget_var_range[1], tget_bins)
        integ_arr = np.linspace(integ_var_range[0], integ_var_range[1], integ_bins)
        phase_interval = np.ones_like(tget_arr) * (integ_var_range[1] - integ_var_range[0]) / integ_bins
        res_arr = np.zeros_like(tget_arr)
        for idx, item in enumerate(tget_arr, 0):
            res = np.expand_dims(param_func(item, integ_arr, fargs), 1)
            resbin = np.apply_along_axis(find_bin.find_bin1d, 1, res, seq=var)
            res_arr[idx] = (func_value[resbin]*phase_interval).sum()
        return tget_arr, res_arr

    #old: def integ2d(y_range, y_func, x_range, s_range, an_range, xbins, sbins, anbins, relation1, relation2):
    @staticmethod
    def integ2d(var, func_value, first_tget_var_range, second_tget_var_range, \
                first_integ_var_range, second_integ_var_range, first_tget_bins, second_tget_bins, \
                first_integ_bins, second_integ_bins, first_param_func, second_param_func, \
                f1args=None, f2args=None):
        r'''
        This method will do a 2d integral over a pre-known-valued function, i.e. the discrete function_value-variable phase point.
        Here shows the meaning of arguments:
        func_value = func(var)
        var = param_func(tget_var, integ_var)
        func(tget_var) = \integral {integ_var_low, integ_var_high, func(param_func(tget_var, integ_var))d(integ_var)}
        Note that for a multiple process, the second_integ_var_range and corresponding bins might be quite large. Whatever, it
            depends on the size of factor, i.e. if you do a integral over the angle, or sin\cos like function, this factor is small
            so it is okay.
        '''
        tget_arr, res_arr = integ.integ1d(var, func_value, first_tget_var_range, first_integ_var_range, first_tget_bins, first_integ_bins, first_param_func, f1args)
        tget_arr, res_arr = integ.integ1d(tget_arr, res_arr, second_tget_var_range, second_integ_var_range, second_tget_bins, second_integ_bins, second_param_func, f2args)
        return tget_arr, res_arr

class rdg():
    def __init__(self, function, xlim, ylim):
        self.function = function
        self.xlim = xlim
        self.ylim = ylim
        return
    
    def get_random(self):
        while True:
            x = rd.uniform(self.xlim[0], self.xlim[1])
            f = rd.uniform(self.xlim[0], self.xlim[1])
            fx = self.function(x)
            if f <= fx:
                return x
            
class nucleon_generator():
    def __init__(self, woods_saxon):
        self.A = woods_saxon.A
        self.r = []
        ymax = woods_saxon.function(0)
        ws_rdg = rdg(woods_saxon.function, [0, 15], [0, ymax])
        for _ in range(self.A):
            self.r.append(ws_rdg.get_random())
        self.r = np.array(self.r)
        self.theta = rd.uniform(0, np.pi, self.A)
        self.phi = rd.uniform(0, 2*np.pi, self.A)
        self.s = self.r*np.sin(self.theta)
        return

    def get_particle(self):
        r'''
        Return two array: rho and phi of nuclei matter of this particle.
        '''
        return self.s, self.phi

class bsq():
    def __init__(self, bmax):
        self.bmax = bmax
        self.a = 1.0 / bmax**3 * 3
        return

    def function(self, x):
        return self.a*x**2