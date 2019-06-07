#!/usr/bin/env python

import os
from datetime import datetime
from math import  sqrt, exp, log

from scipy.special import erf, erfinv
import numpy as np
import netCDF4

from optics_ import mie
import gads


psd_number  = 1      # number  size distribution 
psd_surface = 2      # surface size distrubution

scheme_cs    = 'Chou-Suarez'
scheme_rrtmg = 'RRTMG'

bands_cs     = {'shortwave': ( ( 0.175,  0.225), # 0
                               ( 0.225,  0.285), # 1
                               ( 0.285,  0.300), # 2 := 2 & 0
                               ( 0.300,  0.325), # 3
                               ( 0.325,  0.400), # ...
                               ( 0.400,  0.690),
                               ( 0.690,  1.220),
                               ( 1.220,  2.270),
                               ( 2.270,  3.850) ),
                'longwave':  ( (29.412, 40.000),
                               (18.519, 29.412),
                               (12.500, 18.519),
                               (10.204, 12.500),
                               ( 9.091, 10.204),
                               ( 8.230,  9.091),
                               ( 7.246,  8.230),
                               ( 5.263,  7.246),
                               ( 3.333,  5.263),
                               (16.129, 18.519) ),
                'units': '1e-6 m'}


bands_rrtmg  = {'shortwave':  ( ( 2600,  3250),
                                ( 3250,  4000),
                                ( 4000,  4650),
                                ( 4650,  5150),
                                ( 5150,  6150), 
                                ( 6150,  7700), 
                                ( 7700,  8050), 
                                ( 8050, 12850), 
                                (12850, 16000), 
                                (16000, 22650), 
                                (22650, 29000), 
                                (29000, 38000),
                                (38000, 50000), 
                                (  820,  2600) ),
                 'longwave':  ( (   10,   350), 
                                (  350,   500),
                                (  500,   630),
                                (  630,   700),
                                (  700,   820),
                                (  820,   980),
                                (  980,  1080),
                                ( 1080,  1180),
                                ( 1180,  1390),
                                ( 1390,  1480),
                                ( 1480,  1800),
                                ( 1800,  2080),
                                ( 2080,  2250),
                                ( 2250,  2380),
                                ( 2380,  2600),
                                ( 2600,  3250) ),
                 'units': 'cm-1'}


def geos5_bands(scheme=scheme_cs, units='m', spectrum=None):

    '''
    Returns bands of GEOS-5 radiation scheme.
    '''

    assert scheme == scheme_cs

    assert units == 'm'


    bands = None

    if scheme == scheme_cs:
        bands = bands_cs

    if scheme == scheme_rrtmg:
        bands = bands_rrtmg


    if spectrum != None:
        if spectrum.lower() == 'shortwave' or spectrum.lower() == 'solar':
            _spectrum = 'shortwave'

        if spectrum.lower() == 'longwave'  or spectrum.lower() == 'infrared':
            _spectrum = 'longwave'
    else:
        _spectrum = None
         
    
    if _spectrum == None:
       result = 1e-6 * np.array(bands['shortwave'] + bands['longwave'])
    else:
       result = 1e-6 * np.array(bands[_spectrum])

    return result


class LUT:
    '''
    Class to compute and write aerosol optics lookup tables 
    following the approach of Ghan and Zaveri (JGR, 2007).
    '''

    def __init__(self, wavelengths=None, bands=None, 
                       N_re=10, N_im=20, N_size=100, N_cheb=20,
                       components=('OC', 'BC', 'SU', 'SS', 'DU', 'AMM', 'WATER'), 
                       sigma=2.0, surface_mode_diameter=(2*0.01e-6, 2*25.0e-6), 
                       deliq=0.80, cryst=0.35,
                       scheme=scheme_cs, dir='/home/adarmeno/sandbox/colarco/radiation/gads/optdat/',
                       N_integration_bins=10000, verbose=False):


        assert (wavelengths != None) != (bands != None)

        # LUT's parameters
        self.wavelengths = wavelengths                       # monochromatic wavelengths, 'm'
        self.bands = bands                                   # bands,                     'm'

        self.N_re = N_re                                     # number of discrete refractive indexes
        self.N_im = N_im                                     # for real and imaginary parts

        self.N_size = N_size                                 # number of discrete sizes
        
        self.N_cheb = N_cheb                                 # number of Chebyshev's expansion terms  

        self.N_integration_bins = N_integration_bins         # number of bins used to integrate optics over the size distribution

        # aerosol mode parameters
        self.components = components                         # aerosol components in aerosol mode
        self.sigma      = sigma                              # geometric standard deviation of aerosol mode size distribution
        self.deliq      = deliq                              # deliq. point
        self.cryst      = cryst                              # cryst. point
        self.Dgs_min    = np.min(surface_mode_diameter[:2])  # range of surface mode sizes, 'm'
        self.Dgs_max    = np.max(surface_mode_diameter[:2])  # ...

        # misc
        self.verbose = verbose                               # verbosity
        self.scheme  = scheme                                # flag to handle the odd/compound bands in Chow-Suarez scheme
        self.gads_optics_dir = dir                           # full path to GADS aerosol files


        # quire GADS optics
        n_re_min, n_re_max, n_im_min, n_im_max, n = gads.refractive_index(components  = self.components, 
                                                                          bands       = self.bands, 
                                                                          wavelengths = self.wavelengths, 
                                                                          dir         = self.gads_optics_dir, 
                                                                          verbose     = verbose)

        print 'Re(n)', n_re_min, n_re_max
        print 'Im(n)', n_im_min, n_im_max
        print 'n    =', n 

        if self.bands != None and self.scheme == scheme_cs:
            # widths of band 0 and 2
            d0 = self.bands[0][1] - self.bands[0][0]
            d2 = self.bands[2][1] - self.bands[2][0]

            # weights for band 0 and 2
            w0 = d0 / (d0 + d2)
            w2 = d2 / (d0 + d2)

            for s in self.components:
                n[s]['re'][2] = w0*n[s]['re'][0] + w2*n[s]['re'][2]
                n[s]['im'][2] = w0*n[s]['im'][0] + w2*n[s]['im'][2] 

                n_re_min[2] = np.min((n_re_min[2], n[s]['re'][2]))
                n_re_max[2] = np.max((n_re_min[2], n[s]['re'][2]))
                n_im_min[2] = np.min((n_im_min[2], n[s]['im'][2]))
                n_im_max[2] = np.max((n_im_min[2], n[s]['im'][2]))
            

        print 'Re(n)', n_re_min, n_re_max
        print 'Im(n)', n_im_min, n_im_max
        print 'n    =', n

        # set refractive indexes for sampling
        self.n_component = n      # refractive indexes of the aerosol components at wavelengths/bands

        self.n_re = None          # sampled refractive_indexes - real part
        self.n_im = None          # sampled refractive_indexes - imaginary part

        if self.bands != None:
            self.n_re = [self.__sample_refractive_index_real(min=n_re_min[b], max=n_re_max[b]) for b in range(len(self.bands))]
            self.n_im = [self.__sample_refractive_index_imaginary(min=n_im_min[b], max=n_im_max[b]) for b in range(len(self.bands))]

        if self.wavelengths != None:
            self.n_re = [self.__sample_refractive_index_real(min=n_re_min[w], max=n_re_max[w]) for w in range(len(self.wavelengths))]
            self.n_im = [self.__sample_refractive_index_imaginary(min=n_im_min[w], max=n_im_max[w]) for w in range(len(self.wavelengths))]


        # set sizes for sampling
        self.Dgs = self.__sample_size()


    def __sample_refractive_index_real(self, min, max):
        '''
        Returns the values of the real component of the refractive index
        used to create the Mie LUT.

        Ghan and Zaveri (JGR, 2007) recommend values that are equally 
        spaced.
        '''
       
        result = np.linspace(min, max, num=self.N_re, endpoint=True)
        return result


    def __sample_refractive_index_imaginary(self, min, max):
        '''
        Returns the values of the imaginary component of the refractive index
        used to create the Mie LUT.

        Ghan and Zaveri (JGR, 2007) recommend values that are equally 
        spaced in log-space.
        '''
 
        # note that 0.0 is added as the first value
        result = [0.0, ] + [max * pow(max/min, float(n - self.N_im + 1)/(self.N_im-2)) for n in range(1, self.N_im)]
        return result


    def __sample_size(self, logspace=True):
        '''
        Returns discretized values of aerosol mode size.
        '''

        if logspace:
            # x = (2log(D) - log(D_min) - log(D_max)) / (log(D_max) - log(D_min))
            x = np.linspace(-1.0, 1.0, num=self.N_size, endpoint=True) 
            result = np.exp(0.5 * (x * (log(self.Dgs_max) - log(self.Dgs_min)) + log(self.Dgs_max) + log(self.Dgs_min)))
        else:
            result = np.linspace(self.Dgs_min, self.Dgs_max, num=self.N_size, endpoint=True)

        return result


    def __psd_lognormal_bounds(self, Dg, sigma, eps=1e-5):
        '''
        Computes lower and upper bounds such that
            CDF(x_low) = eps/2 = (1 - integral)/2
            1 - CDF(x_up) = eps/2 = (1 - integral)/2       

        The PSD is assumed to be normalized.
        '''

        def CDF(x, Dg, sigma):
            return 0.5 + 0.5*erf((log(x) - log(Dg))/sqrt(2*log(sigma)*log(sigma)))

        def CDFinv(cdf, Dg, sigma):
            return Dg * exp(sqrt(2*log(sigma)*log(sigma)) * erfinv(2*cdf - 1))

        D_low = CDFinv(0.5*eps, Dg, sigma)
        D_up  = CDFinv(1.0 - 0.5*eps, Dg, sigma)

        return (D_low, D_up)
    

    def _monochromatic_optics(self, refractive_index=complex(1.75, -4.5e-1), wavelength=0.550e-6, psd=psd_surface):
        '''
        Monochromatic optics calculations for a list of surface sizes and 
        fixed wavelength and refractive index:

            wavelength              wavelength for Mie calculations, 'm'
            integration_bins        number of bins used to integrate the size distribution
            psd                     particle size distribution, psd_number=1 or psd_surface=2

        Returns mass specific extinction, scattering and asymmetry parameter.
        '''

        # integration methods
        nim_midpoint = 1
        nim_simpson  = 2

        ext = np.zeros(self.N_size)
        sca = np.zeros_like(ext)
        g   = np.zeros_like(ext)

        for i in range(self.N_size):

            D = self.Dgs[i]
            D_min, D_max = self.__psd_lognormal_bounds(D, self.sigma)

            ext[i], sca[i], g[i] = mie.scattering_lognormal(psd,
                                                            D,
                                                            self.sigma,
                                                            refractive_index,
                                                            wavelength,
                                                            D_min,
                                                            D_max,
                                                            intervals=self.N_integration_bins,
                                                            specific=True,
                                                            method=nim_midpoint)

        return ext, sca, g


    def __compute_monochromatic_optics(self):
        '''
        Computes monochromatic LUT. 
        '''

        assert self.wavelengths != None

        # x in [-1, 1] 
        x = (2*np.log(self.Dgs) - log(self.Dgs_min) - log(self.Dgs_max)) / (log(self.Dgs_max) - log(self.Dgs_min))


        fit_dims = (self.N_cheb, self.N_re, self.N_im, len(self.wavelengths))

        c_sca = np.zeros(fit_dims)
        c_ext = np.zeros(fit_dims)
        c_g   = np.zeros(fit_dims)

        for i in range(self.N_re):
            for j in range(self.N_im):

                for l in range(len(self.wavelengths)):

                    print i, j, l

                    n = complex(self.n_re[l][i], self.n_im[l][j])
                    w = self.wavelengths[l]

                    ext, sca, g = self._monochromatic_optics(refractive_index=n, wavelength=w)

                    c_ext[:,i,j,l] = np.polynomial.chebyshev.chebfit(x, ext, self.N_cheb-1)
                    c_sca[:,i,j,l] = np.polynomial.chebyshev.chebfit(x, sca, self.N_cheb-1)
                    c_g  [:,i,j,l] = np.polynomial.chebyshev.chebfit(x, g,   self.N_cheb-1) 

        return ext, sca, g, c_ext, c_sca, c_g


    def __compute_bandaveraged_optics(self, glq_order=7):
        '''
        Computes band-averaged LUT. 
        '''

        assert self.bands != None

        # Gauss-Legendre quadrature (n=3)
        glq_w_n3 = (0.8888888888888888, 0.5555555555555556, 0.5555555555555556)
        glq_x_n3 = (0.0000000000000000,-0.7745966692414834, 0.7745966692414834)

        # Gauss-Legendre quadrature (n=5)
        glq_w_n5 = (0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891)
        glq_x_n5 = (0.0000000000000000,-0.5384693101056831, 0.5384693101056831,-0.9061798459386640, 0.9061798459386640)

        # Gauss-Legendre quadrature (n=7)
        glq_w_n7 = (0.4179591836734694, 0.3818300505051189, 0.3818300505051189, 0.2797053914892766, 0.2797053914892766, 0.1294849661688697, 0.1294849661688697)
        glq_x_n7 = (0.0000000000000000, 0.4058451513773972,-0.4058451513773972,-0.7415311855993945, 0.7415311855993945,-0.9491079123427585, 0.9491079123427585)

        # Gauss-Legendre quadrature (n=9)
        glq_w_n9 = (0.3302393550012598, 0.1806481606948574, 0.1806481606948574, 0.0812743883615744, 0.0812743883615744, 0.3123470770400029, 0.3123470770400029, 0.2606106964029354, 0.2606106964029354)
        glq_x_n9 = (0.0000000000000000,-0.8360311073266358, 0.8360311073266358,-0.9681602395076261, 0.9681602395076261,-0.3242534234038089, 0.3242534234038089,-0.6133714327005904, 0.6133714327005904)


        if glq_order == 3:
            glq_w = glq_w_n3
            glq_x = glq_x_n3
        elif glq_order == 5:
            glq_w = glq_w_n5
            glq_x = glq_x_n5
        elif glq_order == 7:
            glq_w = glq_w_n7
            glq_x = glq_x_n7
        else:
            glq_w = glq_w_n9
            glq_x = glq_x_n9


        # x in [-1, 1] 
        x = (2*np.log(self.Dgs) - log(self.Dgs_min) - log(self.Dgs_max)) / (log(self.Dgs_max) - log(self.Dgs_min))

        dims = (self.N_size, self.N_re, self.N_im, len(self.bands))

        ext = np.zeros(dims)
        sca = np.zeros(dims)
        g   = np.zeros(dims)

        for i in range(self.N_re):
            for j in range(self.N_im):

                for l in range(len(self.bands)):
                    
                    print i, j, l
                    n = complex(self.n_re[l][i], self.n_im[l][j])

                    # integrate
                    w_l, w_u = self.bands[l] # band range
                    b_w = 0.5*(w_u - w_l)    # band half width
                    b_c = 0.5*(w_u + w_l)    # band center

                    for q in range(len(glq_w)):
                        w = b_w * glq_x[q] + b_c

                        _ext, _sca, _g = self._monochromatic_optics(refractive_index=n, wavelength=w)

                        ext[:,i,j,l] = ext[:,i,j,l] + 0.5 * glq_w[q] * _ext
                        sca[:,i,j,l] = sca[:,i,j,l] + 0.5 * glq_w[q] * _sca
                        g[:,i,j,l]   = g[:,i,j,l]   + 0.5 * glq_w[q] * _g


        # compound bands
        # code is not safe because it assumes that band[0] and band[2] are in the list of bands
        if self.scheme == scheme_cs:
 
            # widths of band 0 and 2
            d0 = self.bands[0][1] - self.bands[0][0]
            d2 = self.bands[2][1] - self.bands[2][0]

            # weights for band 0 and 2
            w0 = d0 / (d0 + d2)
            w2 = d2 / (d0 + d2)

            ext[:,:,:,2] = w0*ext[:,:,:,0] + w2*ext[:,:,:,2]
            sca[:,:,:,2] = w0*sca[:,:,:,0] + w2*sca[:,:,:,2]
            g  [:,:,:,2] = w0*  g[:,:,:,0] + w2*g  [:,:,:,2]

       
        # ... let's assume we are working with Chow-Suarez scheme for now
        fit_dims = (self.N_cheb, self.N_re, self.N_im, len(self.bands))

        c_sca = np.zeros(fit_dims)
        c_ext = np.zeros(fit_dims)
        c_g   = np.zeros(fit_dims)

        for i in range(self.N_re):
            for j in range(self.N_im):

                for l in range(len(self.bands)):

                   c_ext[:,i,j,l] = np.polynomial.chebyshev.chebfit(x, ext[:,i,j,l], self.N_cheb-1)
                   c_sca[:,i,j,l] = np.polynomial.chebyshev.chebfit(x, sca[:,i,j,l], self.N_cheb-1)
                   c_g  [:,i,j,l] = np.polynomial.chebyshev.chebfit(x,   g[:,i,j,l], self.N_cheb-1) 

        return ext, sca, g, c_ext, c_sca, c_g


    def compute(self):
        '''
        Computes and returns ext, sca, asymmetry paramter and corresponding 
        coefficients of chebyshev expansion.
        '''
        
        if self.wavelengths != None:
            result = self.__compute_monochromatic_optics()

        if  self.bands != None:
            result = self.__compute_bandaveraged_optics()

        return result

    def save(self, file, c_ext=None, c_sca=None, c_g=None, title='', comment='', history='', n_chars=80):

        if self.wavelengths != None:
            self.__save_monochromatic(file, c_ext, c_sca, c_g, title, comment, history, n_chars)

        if  self.bands != None:
            self.__save_band_averaged(file, c_ext, c_sca, c_g, title, comment, history, n_chars)


    def __save_monochromatic(self, file, c_ext, c_sca, c_g, title, comment, history, n_chars):

        f = netCDF4.Dataset(file, 'w', format='NETCDF4')

        # global attributes
        f.Conventions   = 'CF'
        f.institution   = 'NASA '
        f.source        = 'NASA/GSFC/GMAO GEOS-5 Aerosol Group'
        f.contact       = 'Anton Darmenov, anton.s.darmenov@nasa.gov'
        f.creation_date = datetime.isoformat(datetime.now())
        f.title         = title
        f.comment       = comment
        f.references    = 'Ghan, S. J., and R. A. Zaveri (2007, doi:10.1029/2006JD007927); Global Aerosol Data Set (GADS).'
        f.history       = history
        f.aerosol_method= 'modal'
        f.optics        = 'monochromatic'
        f.mode_width    = self.sigma
        f.mode_deliq    = self.deliq
        f.mode_cryst    = self.cryst
        f.Dgs_min       = self.Dgs_min
        f.Dgs_max       = self.Dgs_max
        f.x             = '(2*np.log(Dgs) - log(Dgs_min) - log(Dgs_max)) / (log(Dgs_max) - log(Dgs_min))'

        # dimensions
        f.createDimension('band', len(self.wavelengths))

        f.createDimension('n_re', self.N_re)
        f.createDimension('n_im', self.N_im)
        f.createDimension('k', self.N_cheb)
        f.createDimension('component', len(self.components))
        f.createDimension('nchars', n_chars)
        f.createDimension('range', 2)

        # variables
        var_k          = f.createVariable('k', 'i4', ('k'))
        var_band       = f.createVariable('band', 'i4', ('band'))
        var_wavelength = f.createVariable('wavelength', 'f8', ('range', 'band'))
        var_component  = f.createVariable('component', 'S1', ('component', 'nchars'))

        var_component_n_re = f.createVariable('component_n_re', 'f8', ('component', 'band'))
        var_component_n_im = f.createVariable('component_n_im', 'f8', ('component', 'band'))

        var_n_re = f.createVariable('n_re', 'f8', ('n_re', 'band'))
        var_n_im = f.createVariable('n_im', 'f8', ('n_im', 'band'))

        assert c_ext.shape == (self.N_cheb, self.N_re, self.N_im, len(self.wavelengths))
        assert c_sca.shape == c_ext.shape
        assert c_g.shape   == c_ext.shape

        var_c_ext = f.createVariable('c_ext', 'f8', ('k', 'n_re', 'n_im', 'band'))
        var_c_sca = f.createVariable('c_sca', 'f8', ('k', 'n_re', 'n_im', 'band'))
        var_c_asy = f.createVariable('c_asy', 'f8', ('k', 'n_re', 'n_im', 'band'))

        # variables attributes
        var_band.long_name     = 'radiation_band'
        var_band.standard_name = 'radiation_band'
        var_band.units         = '1'

        var_wavelength.long_name     = 'wavelength_range'
        var_wavelength.standard_name = 'wavelength_range_of_radiation_band'
        var_wavelength.units         = 'm'

        var_k.long_name     = 'degree_of_chebyshev_polynomial'
        var_k.standard_name = 'degree_of_chebyshev_polynomial'
        var_k.units         = '1'

        var_component.long_name     = 'aerosol_component'
        var_component.standard_name = 'name_of_aerosol_component'

        var_component_n_re.long_name     = 'component_refractive_index_real_part'
        var_component_n_re.standard_name = 'real_part_of_refractive_index_of_aerosol_component'
        var_component_n_re.units         = '1'

        var_component_n_im.long_name     = 'component_refractive_index_real_part'
        var_component_n_im.standard_name = 'imaginary_part_of_refractive_index_of_aerosol_component'
        var_component_n_im.units         = '1'

        var_n_re.long_name     = 'refractive_index_real_part'
        var_n_re.standard_name = 'real_part_of_refractive_index_of_ambient_aerosol'
        var_n_re.units         = '1'

        var_n_im.long_name     = 'refractive_index_imaginary_part'
        var_n_im.standard_name = 'imaginary_part_of_refractive_index_of_ambient_aerosol'
        var_n_im.units         = '1'
    
        var_c_ext.long_name     = 'c_ext'
        var_c_ext.standard_name = 'coefficients_of_chebyshev_expansion_of_extinction_coefficient'
        var_c_ext.units         = 'm-2 kg'

        var_c_sca.long_name     = 'c_sca'
        var_c_sca.standard_name = 'coefficients_of_chebyshev_expansion_of_scattering_coefficient'
        var_c_sca.units         = 'm-2 kg'

        var_c_asy.long_name     = 'c_asy'
        var_c_asy.standard_name = 'coefficients_of_chebyshev_expansion_of_asymmetry_parameter'
        var_c_asy.units         = '1'

        # data
        var_band[:] = 1 + np.arange(len(self.wavelengths))

        var_wavelength[0,:] = self.wavelengths[:]
        var_wavelength[1,:] = self.wavelengths[:]

        var_k[:] = range(self.N_cheb)

        for c in range(len(self.components)):
            var_component[c] = netCDF4.stringtoarr(self.components[c], n_chars)

            var_component_n_re[c,:] = self.n_component[self.components[c]]['re']
            var_component_n_im[c,:] = self.n_component[self.components[c]]['im']

        for w in range(len(self.wavelengths)):
            var_n_re[:,w] = self.n_re[w]
            var_n_im[:,w] = self.n_im[w]

        for w in range(len(self.wavelengths)):
            var_c_ext[:,:,:,w] = c_ext[:,:,:,w]
            var_c_sca[:,:,:,w] = c_sca[:,:,:,w]
            var_c_asy[:,:,:,w] = c_g[:,:,:,w]

        f.close()


    def __save_band_averaged(self, file, c_ext, c_sca, c_g, title, comment, history, n_chars):

        f = netCDF4.Dataset(file, 'w', format='NETCDF4')

        # global attributes
        f.Conventions   = 'CF'
        f.institution   = 'NASA '
        f.source        = 'NASA/GSFC/GMAO GEOS-5 Aerosol Group'
        f.contact       = 'Anton Darmenov, anton.s.darmenov@nasa.gov'
        f.creation_date = datetime.isoformat(datetime.now())
        f.title         = title
        f.comment       = comment
        f.references    = 'Ghan, S. J., and R. A. Zaveri (2007, doi:10.1029/2006JD007927); Global Aerosol Data Set (GADS).'
        f.history       = history
        f.aerosol_method= 'modal'
        f.optics        = 'band averaged'
        f.mode_width    = self.sigma
        f.mode_deliq    = self.deliq
        f.mode_cryst    = self.cryst
        f.Dgs_min       = self.Dgs_min
        f.Dgs_max       = self.Dgs_max
        f.x             = '(2*np.log(Dgs) - log(Dgs_min) - log(Dgs_max)) / (log(Dgs_max) - log(Dgs_min))'


        # dimensions
        f.createDimension('band', len(self.bands[1:]))

        f.createDimension('n_re', self.N_re)
        f.createDimension('n_im', self.N_im)
        f.createDimension('k', self.N_cheb)
        f.createDimension('component', len(self.components))
        f.createDimension('nchars', n_chars)
        f.createDimension('range', 2)

        # variables
        var_k          = f.createVariable('k', 'i4', ('k'))
        var_band       = f.createVariable('band', 'i4', ('band'))
        var_wavelength = f.createVariable('wavelength', 'f8', ('range', 'band'))
        var_component  = f.createVariable('component', 'S1', ('component', 'nchars'))

        var_component_n_re = f.createVariable('component_n_re', 'f8', ('component', 'band'))
        var_component_n_im = f.createVariable('component_n_im', 'f8', ('component', 'band'))

        var_n_re = f.createVariable('n_re', 'f8', ('n_re', 'band'))
        var_n_im = f.createVariable('n_im', 'f8', ('n_im', 'band'))

        assert c_ext.shape == (self.N_cheb, self.N_re, self.N_im, len(self.bands))
        assert c_sca.shape == c_ext.shape
        assert c_g.shape   == c_ext.shape

        var_c_ext = f.createVariable('c_ext', 'f8', ('k', 'n_re', 'n_im', 'band'))
        var_c_sca = f.createVariable('c_sca', 'f8', ('k', 'n_re', 'n_im', 'band'))
        var_c_asy = f.createVariable('c_asy', 'f8', ('k', 'n_re', 'n_im', 'band'))

        # variables attributes
        var_band.long_name     = 'radiation_band'
        var_band.standard_name = 'radiation_band'
        var_band.units         = '1'

        var_wavelength.long_name     = 'wavelength_range'
        var_wavelength.standard_name = 'wavelength_range_of_radiation_band'
        var_wavelength.units         = ('m')

        var_k.long_name     = 'degree_of_chebyshev_polynomial'
        var_k.standard_name = 'degree_of_chebyshev_polynomial'
        var_k.units         = '1'

        var_component.long_name     = 'aerosol_component'
        var_component.standard_name = 'name_of_aerosol_component'

        var_component_n_re.long_name     = 'component_refractive_index_real_part'
        var_component_n_re.standard_name = 'real_part_of_refractive_index_of_aerosol_component'
        var_component_n_re.units         = '1'

        var_component_n_im.long_name     = 'component_refractive_index_real_part'
        var_component_n_im.standard_name = 'imaginary_part_of_refractive_index_of_aerosol_component'
        var_component_n_im.units         = '1'

        var_n_re.long_name     = 'refractive_index_real_part'
        var_n_re.standard_name = 'real_part_of_refractive_index_of_ambient_aerosol'
        var_n_re.units         = '1'

        var_n_im.long_name     = 'refractive_index_imaginary_part'
        var_n_im.standard_name = 'imaginary_part_of_refractive_index_of_ambient_aerosol'
        var_n_im.units         = '1'
    
        var_c_ext.long_name     = 'c_ext'
        var_c_ext.standard_name = 'coefficients_of_chebyshev_expansion_of_extinction_coefficient'
        var_c_ext.units         = 'm-2 kg'

        var_c_sca.long_name     = 'c_sca'
        var_c_sca.standard_name = 'coefficients_of_chebyshev_expansion_of_scattering_coefficient'
        var_c_sca.units         = 'm-2 kg'

        var_c_asy.long_name     = 'c_asy'
        var_c_asy.standard_name = 'coefficients_of_chebyshev_expansion_of_asymmetry_parameter'
        var_c_asy.units         = '1'
   
        # data
        var_band[:] = 1 + np.arange(len(self.bands[1:]))
        
        for b in 1 + np.arange(len(self.bands[1:])):
            var_wavelength[:, b-1] = self.bands[b]

        var_k[:] = range(self.N_cheb)

        for c in range(len(self.components)):
            var_component[c] = netCDF4.stringtoarr(self.components[c], n_chars)

            var_component_n_re[c,:] = self.n_component[self.components[c]]['re'][1:]
            var_component_n_im[c,:] = self.n_component[self.components[c]]['im'][1:]

        for b in 1 + np.arange(len(self.bands[1:])):
            var_n_re[:,b-1] = self.n_re[b]
            var_n_im[:,b-1] = self.n_im[b]

        for b in 1 + np.arange(len(self.bands[1:])):
            var_c_ext[:,:,:,b-1] = c_ext[:,:,:,b]
            var_c_sca[:,:,:,b-1] = c_sca[:,:,:,b]
            var_c_asy[:,:,:,b-1] = c_g[:,:,:,b]

        f.close()


def test_specific_monochromatic_optics(figure='optics-specific.png'):

    def num2surf(Dgn, sigma):
        return Dgn * np.exp(2.0*np.log(sigma)*np.log(sigma))
    
    def surf2num(Dgs, sigma):
        return Dgs / np.exp(2.0*np.log(sigma)*np.log(sigma))

    from matplotlib import pyplot as plt

    wavelengths = (0.550e-6,) 

    lut = LUT(wavelengths=wavelengths, 
              surface_mode_diameter=(2*0.01e-6, 2*25.0e-6), 
              sigma=2.0, 
              components=('su', 'water'), 
              N_size=100, 
              N_integration_bins=10000,
              verbose=False)

    ext, sca, g = lut._monochromatic_optics(wavelength=0.550e-6, refractive_index=complex(1.9, -0.6), psd=psd_surface)

    plt.clf()
    
    fig, ax1 = plt.subplots()
    l1 = ax1.plot(1e6*0.5*lut.Dgs, 1e-3*ext      , color='gray',  label=r'specific extinction, $m^2 g^{-1}$')
    l2 = ax1.plot(1e6*0.5*lut.Dgs, 1e-3*(ext-sca), color='red',   label=r'specific absorption, $m^2 g^{-1}$')
    l3 = ax1.plot(1e6*0.5*lut.Dgs, 1e-3*(    sca), color='blue',  label=r'specific scattering, $m^2 g^{-1}$')
    ax1.set_xscale('log', basex=10)

    ax2 = ax1.twinx()
    ax2.set_ylim((0.0, 1.0))
    l4 = ax2.plot(1e6*0.5*lut.Dgs, g,       linestyle='-' , color='black', label=r'asymmetry parameter')
    l5 = ax2.plot(1e6*0.5*lut.Dgs, sca/ext, linestyle='--', color='black', label=r'single scattering albedo')

    lines = l1 + l2 + l3 + l4 + l5
    lbls  = [l.get_label() for l in lines]
    ax1.legend(lines, lbls, loc='lower right', frameon=False, prop={'size': 18})
    
    ax1.set_xlabel(r'Surface mode radius, ${\mu}m$', fontsize=18)

    plt.title(r'Black Carbon, $\lambda = 0.55 {\mu}m$, $n=1.9 -i0.6$, $\sigma = 2$', fontsize=18)
    plt.savefig('optics-specific.psd-surface.png', bbox_inches='tight')


    lut = LUT(wavelengths=wavelengths, 
              surface_mode_diameter=(surf2num(2*0.01e-6,2), surf2num(2*25.0e-6,2)), 
              sigma=2.0, 
              components=('su', 'water'), 
              N_size=100, 
              N_integration_bins=10000,
              verbose=False)

    ext, sca, g = lut._monochromatic_optics(wavelength=0.550e-6, refractive_index=complex(1.9, -0.6), psd=psd_number)
    
    plt.clf()
    plt.plot(1e6*0.5*num2surf(lut.Dgs,2), 1e-3*ext      , color='gray',  label=r'specific extinction, $m^2 g^{-1}$')
    plt.plot(1e6*0.5*num2surf(lut.Dgs,2), 1e-3*(ext-sca), color='red',   label=r'specific absorption, $m^2 g^{-1}$')
    plt.plot(1e6*0.5*num2surf(lut.Dgs,2), 1e-3*(    sca), color='blue',  label=r'specific scattering, $m^2 g^{-1}$')
    plt.plot(1e6*0.5*num2surf(lut.Dgs,2),              g, color='black', label=r'asymmetry parameter')
    plt.xscale('log', basex=10)
    plt.legend(loc='upper right', frameon=False, prop={'size': 18})
    plt.title(r'Black Carbon, $\lambda = 0.55 {\mu}m$, $n=1.9 -i0.6$, $\sigma = 2$', fontsize=18)
    plt.xlabel(r'Surface mode radius, ${\mu}m$', fontsize=18)
    plt.savefig('optics-specific.psd-number.png', bbox_inches='tight')




if __name__ == '__main__':

    tests = True
    
    gads_dir = '/home/adarmeno/sandbox/colarco/radiation/gads/optdat/'

    # greate GADS object to read the GADS/OPAC wavelegths
    insol = gads.GADS(os.path.join(gads_dir, 'inso00'))


    wavelengths_gads = insol.wavelengths()
    bands_geos5 = geos5_bands(scheme=scheme_cs)

    ''' 
    # monochromatic LUT
    wavelengths = (0.50e-6, 0.550e-6, 0.600e-6)

    lut = LUT(wavelengths=wavelengths, 
              surface_mode_diameter=(2*0.01e-6, 2*25.0e-6), 
              sigma=2.0, 
              components=('su', 'water', 'oc', 'bc', 'du'), 
              N_re=3, 
              N_im=4,
              N_cheb=5,
              N_size=10,
              N_integration_bins=20,
              verbose=False)

    ext, sca, g, c_ext, c_sca, c_g = lut.compute()

    lut.save('foo-mono.nc4', c_ext, c_sca, c_g, title='', comment='', history='')
    '''

    
    # band-averaged LUT
    bands = bands_geos5[:]

    lut = LUT(bands=bands, 
              surface_mode_diameter=(2*0.01e-6, 2*25.0e-6), 
              sigma=2.0, 
              components=('su', 'water', 'oc', 'bc', 'du'), 
              N_re=3, 
              N_im=4,
              N_cheb=5,
              N_size=10,
              N_integration_bins=20,
              verbose=False)

    ext, sca, g, c_ext, c_sca, c_g = lut.compute()

    lut.save('foo-band.nc4', c_ext, c_sca, c_g, title='', comment='', history='')
    
    
    if tests:
        test_specific_monochromatic_optics()
