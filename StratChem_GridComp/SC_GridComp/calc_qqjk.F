      SUBROUTINE calc_qqjk(k,k1Strat,levels,dt,daytime,aj,ak,qqj,qqk,
     &              oh,o1d,cl,ch4,n2o,f11,f12,f113,f114,f115,ch3ccl3,
     &         ccl4,ch3cl,ch3br,chbr3,ch2br2,ch2brcl,chbrcl2,chbr2cl,
     &              h1301,h1211,h1202,h2402,hcfc22,hcfc141b,hcfc142b,
     &                    hfc23,hfc32,hfc125,hfc134a,hfc143a,hfc152a,
     &                                 levels_cal,numreacs,numphoto )

      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      INTEGER, INTENT(IN) :: k,k1Strat,levels,levels_cal
      INTEGER, INTENT(IN) :: numreacs,numphoto
      LOGICAL, INTENT(IN) :: daytime
      REAL, INTENT(IN)    :: dt

      REAL(KIND=DBL), INTENT(IN) :: aj(numphoto),ak(levels,numreacs)
      REAL(KIND=DBL), INTENT(IN) :: oh,o1d,cl,ch4,n2o
      REAL(KIND=DBL), INTENT(IN) :: f11,f12,f113,f114,f115
      REAL(KIND=DBL), INTENT(IN) :: ch3ccl3,ccl4,ch3cl,ch3br
      REAL(KIND=DBL), INTENT(IN) :: chbr3,ch2br2,ch2brcl,chbrcl2,chbr2cl
      REAL(KIND=DBL), INTENT(IN) :: hcfc22,hcfc141b,hcfc142b
      REAL(KIND=DBL), INTENT(IN) :: h1301,h1211,h1202,h2402
      REAL(KIND=DBL), INTENT(IN) :: hfc23,hfc32,hfc125,hfc134a,hfc143a,hfc152a 

      REAL(KIND=DBL), INTENT(INOUT) :: qqk(35),qqj(20)

      IF(k .ge. k1Strat .and. k .le. levels_cal) THEN 
        IF(daytime) THEN
           qqk( 1) = qqk( 1) + ak(k,14)*oh*ch4*dt
           qqk( 2) = qqk( 2) + ak(k,75)*oh*ch3ccl3*dt
           qqk( 3) = qqk( 3) + ak(k,98)*oh*hcfc22*dt

#ifndef REDUCED
           qqk( 4) = qqk( 4) + ak(k,129)*oh*hcfc141b*dt
           qqk( 5) = qqk( 5) + ak(k,130)*oh*hcfc142b*dt
#endif

           qqk( 6) = qqk( 6) + ak(k,16)*oh*ch3cl*dt
           qqk( 7) = qqk( 7) + ak(k,97)*oh*ch3br*dt

#ifndef REDUCED
           qqk( 8) = qqk( 8) + ak(k,131)*oh*chbr3*dt
           qqk( 9) = qqk( 9) + ak(k,132)*oh*ch2br2*dt 
           qqk(10) = qqk(10) + ak(k,133)*oh*ch2brcl*dt
           qqk(11) = qqk(11) + ak(k,134)*oh*chbrcl2*dt
           qqk(12) = qqk(12) + ak(k,135)*oh*chbr2cl*dt
           qqk(13) = qqk(13) + ak(k,136)*oh*hfc23*dt
           qqk(14) = qqk(14) + ak(k,137)*oh*hfc32*dt
           qqk(15) = qqk(15) + ak(k,138)*oh*hfc125*dt
           qqk(16) = qqk(16) + ak(k,139)*oh*hfc134a*dt
           qqk(17) = qqk(17) + ak(k,140)*oh*hfc143a*dt
           qqk(18) = qqk(18) + ak(k,141)*oh*hfc152a*dt
#endif

           qqk(19) = qqk(19) + ak(k,83)*o1d*f11*dt
           qqk(20) = qqk(20) + ak(k,80)*o1d*f12*dt
           qqk(21) = qqk(21) + ak(k,99)*o1d*f113*dt
           qqk(22) = qqk(22) + ak(k,126)*o1d*hcfc22*dt

#ifndef REDUCED
           qqk(23) = qqk(23) + ak(k,127)*o1d*hcfc141b*dt
           qqk(24) = qqk(24) + ak(k,128)*o1d*hcfc142b*dt
#endif

           qqk(25) = qqk(25) + (ak(k,49)+ak(k,57)+ak(k,119))*o1d*ch4*dt
           qqk(26) = qqk(26) + (ak(k,45)+ak(k,78))*o1d*n2o*dt

#ifndef REDUCED
           qqk(27) = qqk(27) + ak(k,142)*o1d*hfc23*dt
           qqk(28) = qqk(28) + ak(k,143)*o1d*hfc32*dt
           qqk(29) = qqk(29) + ak(k,144)*o1d*hfc125*dt
           qqk(30) = qqk(30) + ak(k,145)*o1d*hfc134a*dt
           qqk(31) = qqk(31) + ak(k,146)*o1d*hfc143a*dt
           qqk(32) = qqk(32) + ak(k,147)*o1d*hfc152a*dt
           qqk(33) = qqk(33) + ak(k,148)*o1d*f114*dt
           qqk(34) = qqk(34) + ak(k,149)*o1d*f115*dt
#endif

           qqj( 1) = qqj( 1) + aj(29)*f11*dt
           qqj( 2) = qqj( 2) + aj(30)*f12*dt
           qqj( 3) = qqj( 3) + aj(34)*f113*dt
           qqj( 4) = qqj( 4) + aj(31)*ccl4*dt
           qqj( 5) = qqj( 5) + aj(14)*n2o*dt
           qqj( 6) = qqj( 6) + aj(32)*ch3ccl3*dt
           qqj( 7) = qqj( 7) + aj(21)*h1301*dt
           qqj( 8) = qqj( 8) + aj(22)*h1211*dt

#ifndef REDUCED
           qqj( 9) = qqj( 9) + aj(49)*h1202*dt
           qqj(10) = qqj(10) + aj(50)*h2402*dt
           qqj(11) = qqj(11) + aj(51)*hcfc141b*dt
           qqj(12) = qqj(12) + aj(52)*hcfc142b*dt
           qqj(13) = qqj(13) + aj(44)*chbr3*dt
           qqj(14) = qqj(14) + aj(45)*ch2br2*dt
           qqj(15) = qqj(15) + aj(46)*ch2brcl*dt
           qqj(16) = qqj(16) + aj(47)*chbrcl2*dt
           qqj(17) = qqj(17) + aj(48)*chbr2cl*dt
           qqj(18) = qqj(18) + aj(53)*f114*dt
           qqj(19) = qqj(19) + aj(54)*f115*dt
#endif
 
        ENDIF
      ENDIF 

      IF(k .lt. k1Strat) THEN
        IF(daytime) THEN

           qqk( 1) = qqk( 1) + ak(k,14)*2*oh*ch4*dt
           qqk( 2) = qqk( 2) + ak(k,75)*2*oh*ch3ccl3*dt
           qqk( 3) = qqk( 3) + ak(k,98)*2*oh*hcfc22*dt

#ifndef REDUCED
           qqk( 4) = qqk( 4) + ak(k,129)*2*oh*hcfc141b*dt
           qqk( 5) = qqk( 5) + ak(k,130)*2*oh*hcfc142b*dt
#endif

           qqk( 6) = qqk( 6) + ak(k,16)*2*oh*ch3cl*dt
           qqk( 7) = qqk( 7) + ak(k,97)*2*oh*ch3br*dt

#ifndef REDUCED
           qqk( 8) = qqk( 8) + ak(k,131)*2*oh*chbr3*dt
           qqk( 9) = qqk( 9) + ak(k,132)*2*oh*ch2br2*dt
           qqk(10) = qqk(10) + ak(k,133)*2*oh*ch2brcl*dt
           qqk(11) = qqk(11) + ak(k,134)*2*oh*chbrcl2*dt
           qqk(12) = qqk(12) + ak(k,135)*2*oh*chbr2cl*dt
           qqk(13) = qqk(13) + ak(k,136)*2*oh*hfc23*dt
           qqk(14) = qqk(14) + ak(k,137)*2*oh*hfc32*dt
           qqk(15) = qqk(15) + ak(k,138)*2*oh*hfc125*dt
           qqk(16) = qqk(16) + ak(k,139)*2*oh*hfc134a*dt
           qqk(17) = qqk(17) + ak(k,140)*2*oh*hfc143a*dt
           qqk(18) = qqk(18) + ak(k,141)*2*oh*hfc152a*dt
#endif

           qqj( 6) = qqj( 6) + aj(32)*ch3ccl3*dt
           qqj( 7) = qqj( 7) + aj(21)*h1301*dt
           qqj( 8) = qqj( 8) + aj(22)*h1211*dt

#ifndef REDUCED
           qqj( 9) = qqj( 9) + aj(49)*h1202*dt
           qqj(10) = qqj(10) + aj(50)*h2402*dt
           qqj(11) = qqj(11) + aj(51)*hcfc141b*dt
           qqj(12) = qqj(12) + aj(52)*hcfc142b*dt
           qqj(13) = qqj(13) + aj(44)*chbr3*dt
           qqj(14) = qqj(14) + aj(45)*ch2br2*dt
           qqj(15) = qqj(15) + aj(46)*ch2brcl*dt
           qqj(16) = qqj(16) + aj(47)*chbrcl2*dt
           qqj(17) = qqj(17) + aj(48)*chbr2cl*dt
#endif

        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE calc_qqjk
 

