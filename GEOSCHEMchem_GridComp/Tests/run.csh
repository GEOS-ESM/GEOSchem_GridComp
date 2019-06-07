#!/bin/csh -f
#
#  Update test driver and run it
#

# Make sure test binaries are up to date
# --------------------------------------
  (cd ..; make ut_GEOSCHEM.x)

# Make sure RC files are here
# ---------------------------
  if ( ! (-e Chem_Registry.rc) ) then
    ( cd ..; make test_rc )
  endif

# Make sure there is a symlink to ExtData
# ---------------------------------------
#  if ( ! (-e ExtData) ) then
#    echo "Please add symlink to ExtData directory"
#    exit 1
#  endif

# Run the poor test
# -----------------
  mpirun -np 2 ../ut_GEOSCHEM.x
