#!/usr/bin/env python

import os
import numpy as np

import scipy.interpolate

class GADS:

    def __init__(self, file):
        
        self.file = file

        self.parse_optical_parameters()


    def parse_optical_parameters(self):

        _optical_parameters = []

        with open(self.file) as f:
            for line in f:
                data = line.split()
       
                if data[0] == '#' and len(data[1:]) == 9:
	    
                    try: 
	                values = [float(v) for v in data[1:]]
                    except:
                        continue

                    _optical_parameters.append(values)

        optical_parameters = np.array(_optical_parameters)

        self._wavelength = 1e-6*optical_parameters[:, 0]
        self._refractive_index_real =  np.array(optical_parameters[:,-2])
        self._refractive_index_imag = -np.array(optical_parameters[:,-1])    # note that GADS Im(n) <= 0, we make them positive 
   

    def refractive_index(self, wavelengths=None, extend=True, bands=None, k=1, s=0):
        '''
        Returns the rafractive indexes at list of wavelengths or
        bands.
        '''

        result = {'re': [], 'im': []}

        if wavelengths == None and bands == None:
            result['re'] = self._refractive_index_real
            result['im'] = self._refractive_index_imag
        else:
            if extend:
                N = 100

                _w    = np.zeros(len(self._wavelength) + 2*N)
                _n_re = np.zeros(len(self._refractive_index_real) + 2*N)
                _n_im = np.zeros(len(self._refractive_index_imag) + 2*N)

                _w[:N]      = np.linspace(0.0, np.min(self._wavelength), N, endpoint=False)
                _n_re[:N]   = self._refractive_index_real[0]
                _n_im[:N]   = self._refractive_index_imag[0]

                _w[N:-N] = self._wavelength[:]
                _n_re[N:-N] = self._refractive_index_real[:]
                _n_im[N:-N] = self._refractive_index_imag[:]
                             
                _w[-N:]  = np.linspace(10*np.max(self._wavelength), np.max(self._wavelength), N, endpoint=False)[::-1]
                _n_re[-N:]  = self._refractive_index_real[-1]
                _n_im[-N:]  = self._refractive_index_imag[-1]
            else:
                _w    = self._wavelength
                _n_re = self._refractive_index_imag
                _n_im = self._refractive_index_imag

            spline_n_re = scipy.interpolate.UnivariateSpline(_w, _n_re, k=k, s=s)
            spline_n_im = scipy.interpolate.UnivariateSpline(_w, _n_im, k=k, s=s)

            if wavelengths != None:
                 result['re'] = spline_n_re(wavelengths)
                 result['im'] = spline_n_im(wavelengths)
            else:
                 for b in bands:
                     assert b[1] > b[0]
                     _n_re = spline_n_re.integral(b[0], b[1]) / (b[1] - b[0])
                     _n_im = spline_n_im.integral(b[0], b[1]) / (b[1] - b[0])

                     result['re'].append(_n_re)
                     result['im'].append(_n_im)

                 result['re'] = np.array(result['re'])
                 result['im'] = np.array(result['im'])

        return result

    def wavelengths(self, bands=None):
        '''
        Returns an array of all wavelenghts or the waveleghts in a list of 
        bands.
        '''
        
        if bands == None:
            result = self._wavelength
        else:
            result = []
            for b in band:
                i = np.logical_and(self._wavelengt >= b[0], self.wavelength <= b[1])

                result.append(self._wavelengt[i])
        
        return np.array(result)



def refractive_index(components=('OC', 'BC', 'SU', 'SS', 'DU', 'AMM', 'SOA', 'POM', 'WATER'), bands=None, wavelengths=None, 
                     dir='/home/adarmeno/sandbox/colarco/radiation/gads/optdat/', verbose=False):

    '''
    Returns min and max of real and imaginary part of refractive indexes as well as 
    a dictionary with refractive indexes of aerosol components interpolated at the 
    specified wavelengths or their mean values in bands.
    '''
    
    _components = [c.lower() for c in components]


    # Organic carbon
    waso = GADS(os.path.join(dir, 'waso00'))
    oc    = waso

    # Black carbon
    soot  = GADS(os.path.join(dir, 'soot00'))
    bc    = soot

    # Sulfate
    suso  = GADS(os.path.join(dir, 'suso00'))
    su    = suso

    # Dust: same spectral refractive indexes
    miam  = GADS(os.path.join(dir, 'miam00'))
    micm  = GADS(os.path.join(dir, 'micm00'))
    minm  = GADS(os.path.join(dir, 'minm00'))
    mitr  = GADS(os.path.join(dir, 'mitr00'))
    du    = miam

    # Sea salt: same spectral refractive indexes
    ssam  = GADS(os.path.join(dir, 'ssam00'))
    sscm  = GADS(os.path.join(dir, 'sscm00'))
    ss    = ssam

    # water (clouds)
    cucc = GADS(os.path.join(dir, 'cucc00'))
    cucp = GADS(os.path.join(dir, 'cucp00'))
    cuma = GADS(os.path.join(dir, 'cuma00'))
    stco = GADS(os.path.join(dir, 'stco00'))
    stma = GADS(os.path.join(dir, 'stma00'))
    water= cucc


    species = {}
    if 'oc'    in _components: species['oc' ]   = oc
    if 'bc'    in _components: species['bc' ]   = bc
    if 'su'    in _components: species['su' ]   = su
    if 'ss'    in _components: species['ss' ]   = ss
    if 'du'    in _components: species['du' ]   = du
    if 'amm'   in _components: species['amm']   = su      # ammonium == sulfate
    if 'pom'   in _components: species['pom']   = oc
    if 'soa'   in _components: species['soa']   = oc
    if 'water' in _components: species['water'] = water
    
    #_s = species.keys()[-1]
    #for s in species.keys():
    #    assert np.array_equal(species[s].wavelengths(), species[_s].wavelengths())

    n = None
    if wavelengths != None:
        n = {}
        for s in species.keys():
            n[s] = species[s].refractive_index(wavelengths=wavelengths)
    
    if bands != None:
        n = {}
        for s in species.keys():
            n[s] = species[s].refractive_index(bands=bands)

   
    _s = species.keys()[-1]
    _l = len(n[_s]['re']) # number of bands or wavelengths
    n_re_min = np.zeros(_l)
    n_re_max = np.zeros(_l)
    n_im_min = np.zeros(_l)
    n_im_max = np.zeros(_l)

    for i in range(_l):
        # refractive indexes in this band or at this wavelength
        _n_re = [n[s]['re'][i] for s in species.keys()]
        _n_im = [n[s]['im'][i] for s in species.keys()]
    
        n_re_min[i] = np.min(_n_re)
        n_re_max[i] = np.max(_n_re)

        n_im_min[i] = np.min(_n_im)
        n_im_max[i] = np.max(_n_im)

    if verbose:
        for i in range(_l):
            print 'band %2i    n_re = [%.3e, %.3e]    n_im = [%.3e, %.3e]' % (i+1, n_re_min[i], n_re_max[i], n_im_min[i], n_im_max[i])
         
        print
        print

    result = (n_re_min, n_re_max, n_im_min, n_im_max, n)

    return result


