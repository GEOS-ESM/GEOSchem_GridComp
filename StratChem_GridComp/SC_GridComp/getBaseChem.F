      SUBROUTINE getBaseChem(ox,nox,hno3,n2o5,ho2no2,clono2,clx,hcl,hocl,h2o2,brx,n2o,cl2,
     &                oclo,brcl,hbr,brono2,ch4,hobr,ch3ooh,co,h2o,hno3c,h2oc,f11,f12,f113,
     &               f114,f115,ccl4,hcfc22,hcfc141b,hcfc142b,chbr3,ch2br2,ch2brcl,chbrcl2,
     &                       chbr2cl,hfc23,hfc32,hfc125,hfc134a,hfc143a,hfc152a,co2b,sf6b,
     &              ch3ccl3,ch3cl,ch3br,h1301,h1211,h1202,h2402,o3,o3p,o1d,n,no,no2,no3,h,
     &                       oh,ho2,cl,clo,bro,br,cl2o2,ch2o,ch3o2,ro3ox,sBase,numSpecies)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: numSpecies

      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      REAL(KIND=DBL), INTENT(OUT) :: ox,nox,hno3,n2o5,ho2no2,clono2,clx,hcl
      REAL(KIND=DBL), INTENT(OUT) :: hocl,h2o2,brx,n2o,cl2,oclo,brcl,hbr
      REAL(KIND=DBL), INTENT(OUT) :: brono2,ch4,hobr,ch3ooh,co,h2o,hno3c
      REAL(KIND=DBL), INTENT(OUT) :: h2oc,f11,f12,f113,f114,f115,ccl4,ch3ccl3
      REAL(KIND=DBL), INTENT(OUT) :: hcfc22,hcfc141b,hcfc142b,co2b,sf6b
      REAL(KIND=DBL), INTENT(OUT) :: hfc23,hfc32,hfc125,hfc134a,hfc143a,hfc152a
      REAL(KIND=DBL), INTENT(OUT) :: chbr3,ch2br2,ch2brcl,chbrcl2,chbr2cl
      REAL(KIND=DBL), INTENT(OUT) :: ch3cl,ch3br,h1301,h1211,h1202,h2402
      REAL(KIND=DBL), INTENT(OUT) :: o3,o3p,o1d,n,no,no2,no3,h,oh,ho2,cl,clo
      REAL(KIND=DBL), INTENT(OUT) :: bro,br,cl2o2,ch2o,ch3o2,ro3ox

      REAL(KIND=DBL), INTENT(IN) :: sBase(numSpecies+1)

      INTEGER :: count

      count = 1
      ox = sBase(count)
      count = count+1
      nox = sBase(count)
      count = count+1
      hno3 = sBase(count)
      count = count+1
      n2o5 = sBase(count)
      count = count+1
      ho2no2 = sBase(count)
      count = count+1
      clono2 = sBase(count)
      count = count+1
      clx = sBase(count)
      count = count+1
      hcl = sBase(count)
      count = count+1
      hocl = sBase(count)
      count = count+1
      h2o2 = sBase(count)
      count = count+1
      brx = sBase(count)
      count = count+1
      n2o = sBase(count)
      count = count+1
      cl2 = sBase(count)
      count = count+1
      oclo = sBase(count)
      count = count+1
      brcl = sBase(count)
      count = count+1
      hbr = sBase(count)
      count = count+1
      brono2 = sBase(count)
      count = count+1
      ch4 = sBase(count)
      count = count+1
      hobr = sBase(count)
      count = count+1
      ch3ooh = sBase(count)
      count = count+1
      co = sBase(count)
      count = count+1
      h2o = sBase(count)
      count = count+1
      hno3c = sBase(count)
      count = count+1
      h2oc = sBase(count)
      count = count+1
      f11 = sBase(count)
      count = count+1
      f12 = sBase(count)
      count = count+1
      f113 = sBase(count)
      count = count+1
      f114 = sBase(count)
      count = count+1
      f115 = sBase(count)
      count = count+1
      hcfc22 = sBase(count)
      count = count+1
      hcfc141b = sBase(count)
      count = count+1
      hcfc142b = sBase(count)
      count = count+1
      ccl4 = sBase(count)
      count = count+1
      ch3ccl3 = sBase(count)
      count = count+1
      ch3cl = sBase(count)
      count = count+1
      ch3br = sBase(count)
      count = count+1
      h1301 = sBase(count)
      count = count+1
      h1211 = sBase(count)
      count = count+1
      h1202 = sBase(count)
      count = count+1
      h2402 = sBase(count)
      count = count+1
      chbr3 = sBase(count)
      count = count+1
      ch2br2 = sBase(count)
      count = count+1
      ch2brcl = sBase(count)
      count = count+1
      chbrcl2 = sBase(count)
      count = count+1
      chbr2cl = sBase(count)
      count = count+1
      hfc23 = sBase(count)
      count = count+1
      hfc32 = sBase(count)
      count = count+1
      hfc125 = sBase(count)
      count = count+1
      hfc134a = sBase(count)
      count = count+1
      hfc143a = sBase(count)
      count = count+1
      hfc152a = sBase(count)
      count = count+1
      co2b = sBase(count)
      count = count+1
      sf6b = sBase(count)
      count = count+1
      o3 = sBase(count)
      count = count+1
      o3p = sBase(count)
      count = count+1
      o1d = sBase(count)
      count = count+1
      n = sBase(count)
      count = count+1
      no = sBase(count)
      count = count+1
      no2 = sBase(count)
      count = count+1
      no3 = sBase(count)
      count = count+1
      h = sBase(count)
      count = count+1
      oh = sBase(count)
      count = count+1
      ho2 = sBase(count)
      count = count+1
      cl = sBase(count)
      count = count+1
      clo = sBase(count)
      count = count+1
      bro = sBase(count)
      count = count+1
      br = sBase(count)
      count = count+1
      cl2o2 = sBase(count)
      count = count+1
      ch2o = sBase(count)
      count = count+1
      ch3o2 = sBase(count)
      count = count+1
      ro3ox = sBase(count)

      RETURN
      END SUBROUTINE getBaseChem
