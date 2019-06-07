#!/usr/bin/env python

import numpy as np

import gads


if __name__ == '__main__':

    from matplotlib import pyplot as plt

    OPT_DATA = '/home/adarmeno/sandbox/colarco/radiation/gads/optdat/'

    BANDS_CS     = {'shortwave': ( ( 0.225,  0.285),
                                   ( 0.175,  0.225), (0.285,  0.300),
                                   ( 0.300,  0.325),
                                   ( 0.325,  0.400),
                                   ( 0.400,  0.690),
                                   ( 0.690,  1.220),
                                   ( 1.220,  2.270),
                                   ( 2.270,  3.850) ),
                    'longwave':  ( (18.519, 29.412),
                                   (12.500, 18.519),
                                   (10.204, 12.500),
                                   ( 9.091, 10.204),
                                   ( 8.230,  9.091),
                                   ( 7.246,  8.230),
                                   ( 5.263,  7.246),
                                   ( 3.333,  5.263),
                                   (16.129, 18.519) ),
                    'units': '1e-6 m'}


    BANDS_RRTMG  = {'shortwave': ( ( 2600,  3250),
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


    # test min/max refractive insex of aerosol components
    bands = 1e-6 * np.array(BANDS_CS['shortwave'] + BANDS_CS['longwave'])
    _ = refractive_index(components=('OC', 'BC', 'SU', 'SS', 'DU', 'AMM', 'WATER'), bands=bands, verbose=True )


    # Organic carbon
    insol = GADS(os.path.join(OPT_DATA, 'inso00'))
    oc    = insol

    # Black carbon
    soot  = GADS(os.path.join(OPT_DATA, 'soot00'))
    bc    = soot


    # Sulfate
    waso  = GADS(os.path.join(OPT_DATA, 'waso00'))
    su    = waso

    # Dust: same spectral refractive indexes
    miam  = GADS(os.path.join(OPT_DATA, 'miam00'))
    micm  = GADS(os.path.join(OPT_DATA, 'micm00'))
    minm  = GADS(os.path.join(OPT_DATA, 'minm00'))
    mitr  = GADS(os.path.join(OPT_DATA, 'mitr00'))
    du    = miam

    # Sea salt: same spectral refractive indexes
    ssam  = GADS(os.path.join(OPT_DATA, 'ssam00'))
    sscm  = GADS(os.path.join(OPT_DATA, 'sscm00'))
    ss    = ssam

    # water (clouds)
    cucc = GADS(os.path.join(OPT_DATA, 'cucc00'))
    cucp = GADS(os.path.join(OPT_DATA, 'cucp00'))
    cuma = GADS(os.path.join(OPT_DATA, 'cuma00'))
    stco = GADS(os.path.join(OPT_DATA, 'stco00'))
    stma = GADS(os.path.join(OPT_DATA, 'stma00'))
    water= cucc


    species = {'OC'   : (oc,    'gray' ),
               'BC'   : (bc,    'black'),
               'SU'   : (su,    'green'),
               'SS'   : (ss,    'lightblue' ),
               'DU'   : (du,    'red'  ),
               'WATER': (water, 'blue' )}
    

    assert np.array_equal(oc.wavelengths(), bc.wavelengths())
    assert np.array_equal(oc.wavelengths(), du.wavelengths())
    assert np.array_equal(oc.wavelengths(), su.wavelengths())
    assert np.array_equal(oc.wavelengths(), ss.wavelengths())


    # GADS/OPAC wavelength in 'um'
    w = 1e6*oc.wavelengths()
    print 'GADS/OPAC wavelengths [um] = ', w
    
    BANDS = BANDS_CS
    

    band        = 1e-6*np.array(BANDS['shortwave'] + BANDS['longwave'])
    band_center = [0.5*(b[0] + b[1]) for b in band]
    band_lbound = [b[0] for b in band]
    band_ubound = [b[1] for b in band]

    n_oc_center = oc.refractive_index(wavelengths=band_center)
    n_oc_lbound = oc.refractive_index(wavelengths=band_lbound)
    n_oc_ubound = oc.refractive_index(wavelengths=band_ubound)
    n_oc_ave    = oc.refractive_index(bands=band)

    n_bc_center = bc.refractive_index(wavelengths=band_center)
    n_bc_lbound = bc.refractive_index(wavelengths=band_lbound)
    n_bc_ubound = bc.refractive_index(wavelengths=band_ubound)
    n_bc_ave    = bc.refractive_index(bands=band)

    n_su_center = su.refractive_index(wavelengths=band_center)
    n_su_lbound = su.refractive_index(wavelengths=band_lbound)
    n_su_ubound = su.refractive_index(wavelengths=band_ubound)
    n_su_ave    = su.refractive_index(bands=band)

    n_du_center = du.refractive_index(wavelengths=band_center)
    n_du_lbound = du.refractive_index(wavelengths=band_lbound)
    n_du_ubound = du.refractive_index(wavelengths=band_ubound)
    n_du_ave    = du.refractive_index(bands=band)

    n_ss_center = ss.refractive_index(wavelengths=band_center)
    n_ss_lbound = ss.refractive_index(wavelengths=band_lbound)
    n_ss_ubound = ss.refractive_index(wavelengths=band_ubound)
    n_ss_ave    = ss.refractive_index(bands=band)

    n_wt_center = water.refractive_index(wavelengths=band_center)
    n_wt_lbound = water.refractive_index(wavelengths=band_lbound)
    n_wt_ubound = water.refractive_index(wavelengths=band_ubound)
    n_wt_ave    = water.refractive_index(bands=band)


    print '-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    print '                 Chou-Suarez                      OC                   BC                 SU                  DU                  SS                Water           Effective      '
    print '-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    print 'band    center[um]       bounds[um]         center  averaged    center  averaged    center  averaged    center  averaged    center  averaged    center  averaged    min     max    '
    print '-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    for n in range(len(band)): 
        print '%2d      %7.4f       %7.4f, %7.4f      %.4f  %.4f      %.4f  %.4f      %.4f  %.4f      %.4f  %.4f      %.4f  %.4f      %.4f  %.4f      %.4f  %.4f' % (n+1, 1e6*band_center[n], 1e6*band[n][0], 1e6*band[n][1],
                                                                            n_oc_center['re'][n], n_oc_ave['re'][n],
                                                                            n_bc_center['re'][n], n_bc_ave['re'][n],
                                                                            n_su_center['re'][n], n_su_ave['re'][n],
                                                                            n_du_center['re'][n], n_du_ave['re'][n],
                                                                            n_ss_center['re'][n], n_ss_ave['re'][n],
                                                                            n_wt_center['re'][n], n_wt_ave['re'][n],
                                                                            np.min((n_oc_center['re'][n], n_oc_ave['re'][n],
                                                                                    n_bc_center['re'][n], n_bc_ave['re'][n],
                                                                                    n_su_center['re'][n], n_su_ave['re'][n],
                                                                                    n_du_center['re'][n], n_du_ave['re'][n],
                                                                                    n_ss_center['re'][n], n_ss_ave['re'][n],
                                                                                    n_wt_center['re'][n], n_wt_ave['re'][n])),
                                                                            np.max((n_oc_center['re'][n], n_oc_ave['re'][n],
                                                                                    n_bc_center['re'][n], n_bc_ave['re'][n],
                                                                                    n_su_center['re'][n], n_su_ave['re'][n],
                                                                                    n_du_center['re'][n], n_du_ave['re'][n],
                                                                                    n_ss_center['re'][n], n_ss_ave['re'][n],
                                                                                    n_wt_center['re'][n], n_wt_ave['re'][n])))
    print '-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    print

    print '-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    print '                 Chou-Suarez                      OC                   BC                 SU                  DU                  SS                Water           Effective    '
    print '-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    print 'band    center[um]       bounds[um]         center  averaged    center  averaged    center  averaged    center  averaged    center  averaged    center  averaged    min     max'
    print '-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    for n in range(len(band)): 
        print '%2d        %7.4f       %7.4f, %7.4f    %.4f %.4f     %.4f %.4f     %.4f %.4f     %.4f %.4f     %.4f %.4f     %.4f %.4f     %.4f %.4f' % (n+1, 1e6*band_center[n], 1e6*band[n][0], 1e6*band[n][1],
                                                                            n_oc_center['im'][n], n_oc_ave['im'][n],
                                                                            n_bc_center['im'][n], n_bc_ave['im'][n],
                                                                            n_su_center['im'][n], n_su_ave['im'][n],
                                                                            n_du_center['im'][n], n_du_ave['im'][n],
                                                                            n_ss_center['im'][n], n_ss_ave['im'][n],
                                                                            n_wt_center['im'][n], n_wt_ave['im'][n],
                                                                            np.min((n_oc_center['im'][n], n_oc_ave['im'][n],
                                                                                    n_bc_center['im'][n], n_bc_ave['im'][n],
                                                                                    n_su_center['im'][n], n_su_ave['im'][n],
                                                                                    n_du_center['im'][n], n_du_ave['im'][n],
                                                                                    n_ss_center['im'][n], n_ss_ave['im'][n],
                                                                                    n_wt_center['im'][n], n_wt_ave['im'][n])),
                                                                            np.max((n_oc_center['im'][n], n_oc_ave['im'][n],
                                                                                    n_bc_center['im'][n], n_bc_ave['im'][n],
                                                                                    n_su_center['im'][n], n_su_ave['im'][n],
                                                                                    n_du_center['im'][n], n_du_ave['im'][n],
                                                                                    n_ss_center['im'][n], n_ss_ave['im'][n],
                                                                                    n_wt_center['im'][n], n_wt_ave['im'][n])))

    print '-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    print


    print 
    print band_center[0],                    n_ss_center['im'][0], n_ss_center['re'][0]
    print band_lbound[0],                    n_ss_lbound['im'][0], n_ss_lbound['re'][0]
    print band_ubound[0],                    n_ss_ubound['im'][0], n_ss_ubound['re'][0]
    print band_lbound[0], band_ubound[0],       n_ss_ave['im'][0],    n_ss_ave['re'][0]



    plt.clf()
    plt.title('GADS/OPAC refractive indexes - real part')
    for name, (s, color) in species.items():
        n = s.refractive_index(wavelengths=band_center)
        print name, n['re']
        print name, 'min/max =', np.min(n['re']), np.max(n['re'])

        n = s.refractive_index(bands=band)
        print name, n['re']
        print name, 'min/max =', np.min(n['re']), np.max(n['re'])
        print
        
        n_gads_opac = s.refractive_index()
        plt.plot(1e6*s.wavelengths(), n_gads_opac['re'], #label=name, 
                                                         linestyle='', 
                                                         marker='o', 
                                                         mec=color, 
                                                         mfc=color, 
                                                         ms=3, 
                                                         mew=1)

        w = 1e-6*np.linspace(0.01, 100.0, 10000)
        n = s.refractive_index(wavelengths=w)
        plt.plot(1e6*w, n['re'], label=name, 
                                 linestyle='-',
                                 color=color)
        plt.xlim(0, 40.0)
        plt.ylim(0, 4)                                                      
        plt.legend(loc='upper left', frameon=False)
        plt.xlabel('wavelength, um')
        plt.savefig('plt-refr_index.real.opac.png')



    plt.clf()
    plt.title('GADS/OPAC refractive indexes - imaginary part')
    for name, (s, color) in species.items():
        n = s.refractive_index(wavelengths=band_center)
        print name, n['im']
        print name, 'min/max =', np.min(n['im']), np.max(n['im'])

        n = s.refractive_index(bands=band)
        print name, n['im']
        print name, 'min/max =', np.min(n['im']), np.max(n['im'])
        print
        

        n_gads_opac = s.refractive_index() 
        plt.plot(1e6*s.wavelengths(), n_gads_opac['im'], #label=name, 
                                                         linestyle='', 
                                                         marker='o', 
                                                         mec=color, 
                                                         mfc=color, 
                                                         ms=3, 
                                                         mew=1)

      
        w = 1e-6*np.linspace(0.01, 100, 10000)
        n = s.refractive_index(wavelengths=w, extend=True)
        plt.plot(1e6*w, n['im'], label=name, 
                                 linestyle='-',
                                 color=color)
        plt.xlim(0.0, 40.0)
        plt.ylim(-2.0,0.2) 
        plt.legend(loc='lower left', frameon=False)
        plt.xlabel('wavelength, um')
        plt.savefig('plt-refr_index.imaginary.opac.png')

         
