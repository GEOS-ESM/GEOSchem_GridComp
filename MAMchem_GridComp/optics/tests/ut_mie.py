#! /usr/bin/env python

import optics_

PSD_NUMBER  = 1
PSD_SURFACE = 2

def scattering(psd, Dg, sigma, refractive_index, wavelength, D_min, D_max, N=10000, specific=True, method=0):
    
    ext, sca, asy = optics_.mie.scattering_lognormal(psd, Dg, sigma, refractive_index, wavelength, D_min, D_max, N, specific, method)

    return (ext, sca, asy)


if __name__ == '__main__':

   # test parameters ----------------------------
   verbose            = False

   bins               = 10000
   wavelength         = 0.50e-6
   component          = 'soot'
 
   # OPAC parameters of the *number* size distribution
   #
   # Source(s): 
   #     1. GEISA-2003 OPAC Database Content
   #     2. Van de Hulst, 1981
   #     3. OPAC software package (BAMS98)
   #         
   # URL(s): 
   #     1. http://ether.ipsl.jussieu.fr/etherTypo/index.php?id=1058&L=0
   #     2. http://ether.ipsl.jussieu.fr/etherTypo/fileadmin/files/GEISA/hess-etal2.pdf
   #

   OPAC = {}
   OPAC['insoluble'           ] = dict(m = complex(1.53, -8.0e-3), sigma = 2.51, median = 0.471,  r_min = 0.005, r_max = 20.0)
   OPAC['water-soluble'       ] = dict(m = complex(0.00, -0.00  ), sigma = 2.24, median = 0.0212, r_min = 0.005, r_max = 20.0) 
   OPAC['soot'                ] = dict(m = complex(1.75, -4.5e-1), sigma = 2.00, median = 0.0118, r_min = 0.005, r_max = 20.0) 
   OPAC['sea salt (acc. mode)'] = dict(m = complex(0.00, -0.00  ), sigma = 2.03, median = 0.209 , r_min = 0.005, r_max = 20.0) 
   OPAC['sea salt (coa. mode)'] = dict(m = complex(0.00, -0.00  ), sigma = 2.03, median = 1.75  , r_min = 0.005, r_max = 60.0) 
   OPAC['mineral (nuc. mode)' ] = dict(m = complex(1.53, -7.8e-3), sigma = 1.95, median = 0.07  , r_min = 0.005, r_max = 20.0) 
   OPAC['mineral (acc. mode)' ] = dict(m = complex(1.53, -7.8e-3), sigma = 2.00, median = 0.39  , r_min = 0.005, r_max = 20.0) 
   OPAC['mineral (coa. mode)' ] = dict(m = complex(1.53, -7.8e-3), sigma = 2.15, median = 1.90  , r_min = 0.005, r_max = 60.0) 
   OPAC['mineral transported' ] = dict(m = complex(1.53, -7.8e-3), sigma = 2.20, median = 0.50  , r_min = 0.020, r_max =  5.0) 
   OPAC['sulfate droplets'    ] = dict(m = complex(0.00, -0.00  ), sigma = 2.03, median = 0.0695, r_min = 0.005, r_max = 20.0) 
   # --------------------------------------------

   if component in OPAC.keys():
       refractive_index = OPAC[component]['m']
       r_min            = OPAC[component]['r_min']  * 1e-6
       r_max            = OPAC[component]['r_max']  * 1e-6
       median           = OPAC[component]['median'] * 1e-6
       sigma            = OPAC[component]['sigma']
   else:
       refractive_index = complex(1.75, -0.45) 
       r_min            = 0.00001  * 1e-6
       r_max            = 1000.00  * 1e-6
       median           = 20.0e-6 #0.0118 * 1e-6
       sigma            = 1.5


   print 'Results using prescribed range of sizes:'
   print 'Size range = (%.3e, %.3e) microns' % (1e6*r_min, 1e6*r_max)

   b_ext, b_sca, g = scattering(PSD_NUMBER,
                                2*median, 
                                sigma, 
                                refractive_index, 
                                wavelength, 
                                2*r_min, 
                                2*r_max, 
                                N=bins, 
                                specific=False) 
   
   print 'Effective extinction factor and scattering factor units are m-1, assuming number density of 1 particle m-3'
   print 'b_ext = %.3e %s' % (b_ext, 'm-1, [N] = 1 particle m-3')
   print 'b_sca = %.3e %s' % (b_sca, 'm-1, [N] = 1 particle m-3')
   print 'g     = %.3e %s' % (g,     '1  , [N] = 1 particle m-3')

   print

   print 'Effective extinction factor and scattering factor units are km-1, assuming number density of 1 particle cm-3'
   print 'b_ext = %.3e %s' % (b_ext * 1e9, 'km-1, [N] = 1 particle cm-3')
   print 'b_sca = %.3e %s' % (b_sca * 1e9, 'km-1, [N] = 1 particle cm-3')
   print 'g     = %.3e %s' % (g,           '1  ,  [N] = 1 particle cm-3')

   print
   print


   print 'Results using automatically set range of sizes:'
   r_low, r_up  = psd_lognormal_bounds(median * exp(2*log(sigma)*log(sigma)), sigma, eps=1e-5)
   print 'Size range = (%.3e, %.3e)' % (1e6*r_low, 1e6*r_up)
   b_ext, b_sca, g = population(psd_lognormal, refractive_index, 
                                               wavelength, 
                                               r_low, 
                                               r_up, 
                                               N=bins, 
                                               median=median,
                                               sigma=sigma)
   
   print 'Effective extinction factor and scattering factor units are m-1, assuming number density of 1 particle m-3'
   print 'b_ext = %.3e %s' % (b_ext, 'm-1, [N] = 1 particle m-3')
   print 'b_sca = %.3e %s' % (b_sca, 'm-1, [N] = 1 particle m-3')
   print 'g     = %.3e %s' % (g,     '1  , [N] = 1 particle m-3')

   print

   print 'Effective extinction factor and scattering factor units are km-1, assuming number density of 1 particle cm-3'
   print 'b_ext = %.3e %s' % (b_ext * 1e9, 'km-1, [N] = 1 particle cm-3')
   print 'b_sca = %.3e %s' % (b_sca * 1e9, 'km-1, [N] = 1 particle cm-3')
   print 'g     = %.3e %s' % (g,           '1  ,  [N] = 1 particle cm-3')
