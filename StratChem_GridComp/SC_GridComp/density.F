        FUNCTION density(cs,cn,t)

        IMPLICIT NONE

! Author: David Considine

        REAL :: cs,cn,t,denss,densn,density

! Liquid solution density g/cm3

        denss = 1000.0+123.64*cs-5.6e-4*cs*t**2 
     &       -29.54*cs**1.5 + 1.814e-4*cs**1.5*t**2
     &       + 2.343*cs**2  -1.487e-3*cs**2*t 
     &       -1.324e-5*cs**2*t**2

        densn = 1000.0+85.107*cn-5.043e-4*cn*t**2 
     &       -18.96*cn**1.5 + 1.427e-4*cn**1.5*t**2
     &       + 1.458*cn**2  -1.198e-3*cn**2*t 
     &       -9.703e-6*cn**2*t**2

        density = 0.001/((1.0/denss*cs/(cs+cn)
     &  +1.0/densn*cn/(cs+cn)))

        RETURN
        END FUNCTION density
