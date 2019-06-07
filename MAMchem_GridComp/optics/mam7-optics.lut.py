#!/usr/bin/env python

'''
Create and save MAM7 mie lookup tables.
'''

import os

import argparse

import gads
import radiation

LUT_VERSION = '0.3-rc2'
LUT_FILE_TEMPLATE = 'opticsBands_MAM7_{mode}.v{version}.nc4'

MAM7 = {}

MAM7['AIT'] = dict(name = 'Aitken', 
                   components = ('water', 'su', 'amm', 'soa', 'ss'),
                   sigma = 1.6, 
                   dileq = 0.80, 
                   cryst = 0.35)
      
MAM7['ACC'] = dict(name = 'Accumulation', 
                   components = ('water', 'su', 'amm', 'soa', 'pom', 'bc', 'ss'),
                   sigma = 1.8, 
                   dileq = 0.80, 
                   cryst = 0.35)

MAM7['PCM'] = dict(name = 'Primary Carbon', 
                   components = ('water', 'pom', 'bc'),
                   sigma = 1.6, 
                   dileq = 0.80, 
                   cryst = 0.35)

MAM7['FSS'] = dict(name = 'Fine Seasalt', 
                   components = ('water', 'su', 'amm', 'ss'),   
                   sigma = 2.0, 
                   dileq = 0.80, 
                   cryst = 0.35) 

MAM7['CSS'] = dict(name = 'Coarse Seasalt', 
                   components = ('water', 'su', 'amm', 'ss'),
                   sigma = 2.0, 
                   dileq = 0.80, 
                   cryst = 0.35)

MAM7['FDU'] = dict(name = 'Fine Dust',
                   components=('water', 'su', 'amm', 'du'),
                   sigma = 1.8, 
                   dileq = 0.80, 
                   cryst = 0.35)

MAM7['CDU'] = dict(name = 'Coarse Dust',
                   components = ('water', 'su', 'amm', 'du'),
                   sigma = 1.8, 
                   dileq = 0.80, 
                   cryst = 0.35)



if __name__ == '__main__':

    '''
    Example usage:

    source $ESMADIR/src/g5_modules

    # create broadband optical properties LUTs for Chou-Suarez scheme
    nohup python ./mam7-optics.lut.py --chou-suarez --mode=ait ./_opticsBands_MAM7_AIT.v0.3-rc2.nc4 >& logs/log.lut-ait.v0.3-rc2 &
    nohup python ./mam7-optics.lut.py --chou-suarez --mode=acc ./_opticsBands_MAM7_ACC.v0.3-rc2.nc4 >& logs/log.lut-acc.v0.3-rc2 &
    nohup python ./mam7-optics.lut.py --chou-suarez --mode=pcm ./_opticsBands_MAM7_PCM.v0.3-rc2.nc4 >& logs/log.lut-pcm.v0.3-rc2 &
    nohup python ./mam7-optics.lut.py --chou-suarez --mode=fss ./_opticsBands_MAM7_FSS.v0.3-rc2.nc4 >& logs/log.lut-fss.v0.3-rc2 &
    nohup python ./mam7-optics.lut.py --chou-suarez --mode=css ./_opticsBands_MAM7_CSS.v0.3-rc2.nc4 >& logs/log.lut-css.v0.3-rc2 &
    nohup python ./mam7-optics.lut.py --chou-suarez --mode=cdu ./_opticsBands_MAM7_CDU.v0.3-rc2.nc4 >& logs/log.lut-cdu.v0.3-rc2 &
    nohup python ./mam7-optics.lut.py --chou-suarez --mode=fdu ./_opticsBands_MAM7_FDU.v0.3-rc2.nc4 >& logs/log.lut-fdu.v0.3-rc2 &

    # create monochromatic (at GADS wavelengths) optical properties LUTs for Chou-Suarez scheme
    nohup python ./mam7-optics.lut.py --gads --mode=ait ./_optics_MAM7_AIT.v0.3-rc2.nc4 >& logs/log.lut-ait &
    nohup python ./mam7-optics.lut.py --gads --mode=acc ./_optics_MAM7_ACC.v0.3-rc2.nc4 >& logs/log.lut-acc &
    nohup python ./mam7-optics.lut.py --gads --mode=pcm ./_optics_MAM7_PCM.v0.3-rc2.nc4 >& logs/log.lut-pcm &
    nohup python ./mam7-optics.lut.py --gads --mode=fss ./_optics_MAM7_FSS.v0.3-rc2.nc4 >& logs/log.lut-fss &
    nohup python ./mam7-optics.lut.py --gads --mode=css ./_optics_MAM7_CSS.v0.3-rc2.nc4 >& logs/log.lut-css &
    nohup python ./mam7-optics.lut.py --gads --mode=cdu ./_optics_MAM7_CDU.v0.3-rc2.nc4 >& logs/log.lut-cdu &
    nohup python ./mam7-optics.lut.py --gads --mode=fdu ./_optics_MAM7_FDU.v0.3-rc2.nc4 >& logs/log.lut-fdu &
    '''


    # parse arguments
    parser = argparse.ArgumentParser(description='MAM7 mie lookup tables generator.')
  
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--gads',        action='store_true', help='Monochromatic optics at the 61 GADS wavelengths.')
    group.add_argument('--chou-suarez', action='store_true', help='Band-averaged optics for the Chow-Suarez radiation scheme.')
    group.add_argument('--rrtmg',       action='store_true', help='Band-averaged optics for the RRTMG radiation scheme.')

    parser.add_argument('--mode', choices=('ait', 'acc', 'pcm', 'fdu', 'cdu', 'fss', 'css'), required=True, help='MAM7 mode: Aitken, Accumulation, Primary Carbon, Fine dust, Coarse dust, Fine seasalt, Coarse seasalt.')
    parser.add_argument('file', help='Output file.')

    args = parser.parse_args()


    # LUT setup
    mode = args.mode.upper()

    # next three choices are exclusive
    if args.gads:
        gads_dir = '/home/adarmeno/sandbox/colarco/radiation/gads/optdat/'

        # greate GADS object to read the GADS/OPAC wavelegths
        insol = gads.GADS(os.path.join(gads_dir, 'inso00'))
        wavelengths_gads = insol.wavelengths()
        
        wavelengths = wavelengths_gads     # ... or reduce to (470e-9, 550e-9, 660e-9, 870e-9)
        bands       = None
    else:
        wavelengths = None
    
        if args.chou_suarez:
            bands = radiation.geos5_bands(scheme=radiation.scheme_cs)

        if args.rrtmg:
            assert False # not fully implemented yet
            bands = radiation.geos5_bands(scheme=radiation.scheme_rrtmg)

    
    description = 'Aerosol optical properties: MAM7 - {name} Mode'.format(name=MAM7[mode]['name'])
    print description

    lut = radiation.LUT(wavelengths=wavelengths,
                        bands=bands,
                        surface_mode_diameter=(2*0.01e-6, 2*25.0e-6), 
                        sigma=MAM7[mode]['sigma'], 
                        components=MAM7[mode]['components'], 
                        N_re=20,
                        N_im=30,
                        N_cheb=20,
                        N_size=200,
                        N_integration_bins=10000,
                        verbose=False)

    ext, sca, g, c_ext, c_sca, c_g = lut.compute()

    lut.save(args.file, c_ext, c_sca, c_g, title=description, comment='', history='')


    

