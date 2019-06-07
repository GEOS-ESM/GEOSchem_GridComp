      SUBROUTINE o3colum(levels,prslay,press,m,o3,no,o3column,o2column,
     &                   nocolumn) 
 
! Calculate the overlying o3 and o2 column.

! 1. Working down from top of model, multiply o3 mixing ratio by pressure 
!    depth of grid box to calculate o3 column in molecules per cm**2 from 
!    bottom of current box to top of model.
! 2. Subtract portion of o3 column that is below center of current box.
! 3. Define the ozone values above 50 km from the gcm ozone climatology.
! 4. Use o3abv80km for upper boundary condition.

! Input: 
!  o3	 = ozone number density in molecules/cm**3
!  m	 = atmospheric number density in molecules/cm**3 

! Output: 
!  o3column = overlying o3 column in molecules/cm**2

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: levels
      REAL, INTENT(IN) :: prslay(0:levels),press(levels)
      REAL, INTENT(IN) :: m(levels),o3(levels),no(levels)

      REAL, INTENT(OUT) :: o3column(levels),o2column(levels)
      REAL, INTENT(OUT) :: nocolumn(levels)

      INTEGER :: k

      REAL :: avgd,o3vmr1,o3vmr2,conv,o3c,noc,o3abv80km,o2c,o2abv80km
      REAL :: pboxbot,pboxtop,pboxcen,noabv80km,o2vmr
      REAL :: novmr1,novmr2

      REAL :: dp(levels),dphalf(levels)

! Column amount above 80 km (Brasseur-Solomon,1986)

      o3abv80km = 1.10E+11 
      noabv80km = 1.90E+14

! Column amount above 80 km (molecules/cm^2)

      o2abv80km = 7.072926E+19

! Constant O2 mixing ratio

      o2vmr = 0.20946

      avgd = 6.022E+23 
      conv = 3.5223E-02
  
      do k=1,levels
       dp(k) = prslay(k-1) - prslay(k)  
       dphalf(k) = prslay(k-1) - press(k) 
      enddo

      o3c = o3abv80km
      o2c = o2abv80km
      noc = noabv80km

      do k=levels,2,-1

       o3vmr1 = o3(k)/m(k)
       o3vmr2 = o3(k-1)/m(k-1)
       novmr1 = no(k)/m(k)
       novmr2 = no(k-1)/m(k-1)

       o3c = o3c+avgd*conv*dp(k)*o3vmr1
       noc = noc+avgd*conv*dp(k)*novmr1
       o2c = o2c+avgd*conv*dp(k)*o2vmr
 
       o3column(k) = o3c-avgd*conv*dphalf(k)*
     &  	     	 (0.75*o3vmr1+0.25*o3vmr2)
       nocolumn(k) = noc-avgd*conv*dphalf(k)*
     &  	     	 (0.75*novmr1+0.25*novmr2)
       o2column(k) = o2c-avgd*conv*dphalf(k)*o2vmr

      enddo

      o3c = o3c+avgd*conv*dp(1)*o3(1)/m(1)
      noc = noc+avgd*conv*dp(1)*no(1)/m(1)
      o2c = o2c+avgd*conv*dp(1)*o2vmr
      o3column(1) = o3c-avgd*conv*dphalf(1)*o3(1)/m(1)
      nocolumn(1) = noc-avgd*conv*dphalf(1)*no(1)/m(1)
      o2column(1) = o2c-avgd*conv*dphalf(1)*o2vmr

      RETURN
      END SUBROUTINE o3colum
