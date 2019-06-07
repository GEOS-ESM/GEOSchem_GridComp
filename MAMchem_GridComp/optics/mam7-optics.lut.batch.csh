#!/bin/csh


ssh dali02 << EOF
    source /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/g5_modules
    cd /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/GEOSgcs_GridComp/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/MAMchem_GridComp/optics/
    nohup python ./mam7-optics.lut.py --chou-suarez --mode=ait ./_opticsBands_MAM7_AIT.v0.3-rc2.nc4 >& logs/log.lut-ait.v0.3-rc2 &
EOF

ssh dali03 << EOF
    source /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/g5_modules
    cd /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/GEOSgcs_GridComp/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/MAMchem_GridComp/optics/
    nohup python ./mam7-optics.lut.py --chou-suarez --mode=acc ./_opticsBands_MAM7_ACC.v0.3-rc2.nc4 >& logs/log.lut-acc.v0.3-rc2 &
EOF

ssh dali04 << EOF
    source /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/g5_modules
    cd /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/GEOSgcs_GridComp/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/MAMchem_GridComp/optics/
    nohup python ./mam7-optics.lut.py --chou-suarez --mode=pcm ./_opticsBands_MAM7_PCM.v0.3-rc2.nc4 >& logs/log.lut-pcm.v0.3-rc2 &
EOF

ssh dali05 << EOF
    source /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/g5_modules
    cd /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/GEOSgcs_GridComp/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/MAMchem_GridComp/optics/
    nohup python ./mam7-optics.lut.py --chou-suarez --mode=fss ./_opticsBands_MAM7_FSS.v0.3-rc2.nc4 >& logs/log.lut-fss.v0.3-rc2 &
EOF

ssh dali06 << EOF
    source /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/g5_modules
    cd /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/GEOSgcs_GridComp/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/MAMchem_GridComp/optics/
    nohup python ./mam7-optics.lut.py --chou-suarez --mode=css ./_opticsBands_MAM7_CSS.v0.3-rc2.nc4 >& logs/log.lut-css.v0.3-rc2 &
EOF

ssh dali07 << EOF
    source /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/g5_modules
    cd /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/GEOSgcs_GridComp/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/MAMchem_GridComp/optics/
    nohup python ./mam7-optics.lut.py --chou-suarez --mode=cdu ./_opticsBands_MAM7_CDU.v0.3-rc2.nc4 >& logs/log.lut-cdu.v0.3-rc2 &
EOF

ssh dali08 << EOF
    source /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/g5_modules
    cd /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/GEOSgcs_GridComp/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/MAMchem_GridComp/optics/
    nohup python ./mam7-optics.lut.py --chou-suarez --mode=fdu ./_opticsBands_MAM7_FDU.v0.3-rc2.nc4 >& logs/log.lut-fdu.v0.3-rc2 &
EOF


ssh dali09 << EOF
    source /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/g5_modules
    cd /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/GEOSgcs_GridComp/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/MAMchem_GridComp/optics/
    nohup python ./mam7-optics.lut.py --gads --mode=ait ./_optics_MAM7_AIT.v0.3-rc2.nc4 >& logs/log.lut-ait &
EOF

ssh dali10 << EOF
    source /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/g5_modules
    cd /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/GEOSgcs_GridComp/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/MAMchem_GridComp/optics/
    nohup python ./mam7-optics.lut.py --gads --mode=acc ./_optics_MAM7_ACC.v0.3-rc2.nc4 >& logs/log.lut-acc &
EOF

ssh dali11 << EOF
    source /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/g5_modules
    cd /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/GEOSgcs_GridComp/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/MAMchem_GridComp/optics/
    nohup python ./mam7-optics.lut.py --gads --mode=pcm ./_optics_MAM7_PCM.v0.3-rc2.nc4 >& logs/log.lut-pcm &
EOF

ssh dali12 << EOF
    source /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/g5_modules
    cd /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/GEOSgcs_GridComp/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/MAMchem_GridComp/optics/
    nohup python ./mam7-optics.lut.py --gads --mode=fss ./_optics_MAM7_FSS.v0.3-rc2.nc4 >& logs/log.lut-fss &
EOF

ssh dali13 << EOF
    source /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/g5_modules
    cd /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/GEOSgcs_GridComp/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/MAMchem_GridComp/optics/
    nohup python ./mam7-optics.lut.py --gads --mode=css ./_optics_MAM7_CSS.v0.3-rc2.nc4 >& logs/log.lut-css &
EOF

ssh dali14 << EOF
    source /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/g5_modules
    cd /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/GEOSgcs_GridComp/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/MAMchem_GridComp/optics/
    nohup python ./mam7-optics.lut.py --gads --mode=cdu ./_optics_MAM7_CDU.v0.3-rc2.nc4 >& logs/log.lut-cdu &
EOF

ssh dali15 << EOF
    source /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/g5_modules
    cd /home/adarmeno/models/geos-5/ganymed-4.0-radiation/src/GEOSgcs_GridComp/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/MAMchem_GridComp/optics/
    nohup python ./mam7-optics.lut.py --gads --mode=fdu ./_optics_MAM7_FDU.v0.3-rc2.nc4 >& logs/log.lut-fdu &
EOF

