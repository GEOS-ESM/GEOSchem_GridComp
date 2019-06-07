      SUBROUTINE storeBaseChem(ox,nox,hno3,n2o5,ho2no2,clono2,clx,hcl,hocl,h2o2,brx,n2o,cl2,
     &                  oclo,brcl,hbr,brono2,ch4,hobr,ch3ooh,co,h2o,hno3c,h2oc,f11,f12,f113,
     &                 f114,f115,ccl4,hcfc22,hcfc141b,hcfc142b,chbr3,ch2br2,ch2brcl,chbrcl2,
     &                         chbr2cl,hfc23,hfc32,hfc125,hfc134a,hfc143a,hfc152a,co2b,sf6b,
     &                ch3ccl3,ch3cl,ch3br,h1301,h1211,h1202,h2402,o3,o3p,o1d,n,no,no2,no3,h,
     &                         oh,ho2,cl,clo,bro,br,cl2o2,ch2o,ch3o2,ro3ox,sBase,numSpecies)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: numSpecies

      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      REAL(KIND=DBL), INTENT(IN) :: ox,nox,hno3,n2o5,ho2no2,clono2,clx,hcl
      REAL(KIND=DBL), INTENT(IN) :: hocl,h2o2,brx,n2o,cl2,oclo,brcl,hbr
      REAL(KIND=DBL), INTENT(IN) :: brono2,ch4,hobr,ch3ooh,co,h2o,hno3c
      REAL(KIND=DBL), INTENT(IN) :: h2oc,f11,f12,f113,f114,f115,ccl4,ch3ccl3
      REAL(KIND=DBL), INTENT(IN) :: hcfc22,hcfc141b,hcfc142b,co2b,sf6b
      REAL(KIND=DBL), INTENT(IN) :: hfc23,hfc32,hfc125,hfc134a,hfc143a,hfc152a
      REAL(KIND=DBL), INTENT(IN) :: chbr3,ch2br2,ch2brcl,chbrcl2,chbr2cl
      REAL(KIND=DBL), INTENT(IN) :: ch3cl,ch3br,h1301,h1211,h1202,h2402
      REAL(KIND=DBL), INTENT(IN) :: o3,o3p,o1d,n,no,no2,no3,h,oh,ho2,cl,clo
      REAL(KIND=DBL), INTENT(IN) :: bro,br,cl2o2,ch2o,ch3o2,ro3ox

      REAL(KIND=DBL), INTENT(OUT) :: sBase(numSpecies+1)

      INTEGER :: count

      count = 1
      sBase(count) = ox
      count = count+1
      sBase(count) = nox
      count = count+1
      sBase(count) = hno3
      count = count+1
      sBase(count) = n2o5
      count = count+1
      sBase(count) = ho2no2
      count = count+1
      sBase(count) = clono2
      count = count+1
      sBase(count) = clx
      count = count+1
      sBase(count) = hcl
      count = count+1
      sBase(count) = hocl
      count = count+1
      sBase(count) = h2o2
      count = count+1
      sBase(count) = brx
      count = count+1
      sBase(count) = n2o
      count = count+1
      sBase(count) = cl2
      count = count+1
      sBase(count) = oclo
      count = count+1
      sBase(count) = brcl
      count = count+1
      sBase(count) = hbr
      count = count+1
      sBase(count) = brono2
      count = count+1
      sBase(count) = ch4
      count = count+1
      sBase(count) = hobr
      count = count+1
      sBase(count) = ch3ooh
      count = count+1
      sBase(count) = co
      count = count+1
      sBase(count) = h2o
      count = count+1
      sBase(count) = hno3c
      count = count+1
      sBase(count) = h2oc
      count = count+1
      sBase(count) = f11
      count = count+1
      sBase(count) = f12
      count = count+1
      sBase(count) = f113
      count = count+1
      sBase(count) = f114
      count = count+1
      sBase(count) = f115
      count = count+1
      sBase(count) = hcfc22
      count = count+1
      sBase(count) = hcfc141b
      count = count+1
      sBase(count) = hcfc142b
      count = count+1
      sBase(count) = ccl4
      count = count+1
      sBase(count) = ch3ccl3
      count = count+1
      sBase(count) = ch3cl
      count = count+1
      sBase(count) = ch3br
      count = count+1
      sBase(count) = h1301
      count = count+1
      sBase(count) = h1211
      count = count+1
      sBase(count) = h1202
      count = count+1
      sBase(count) = h2402
      count = count+1
      sBase(count) = chbr3
      count = count+1
      sBase(count) = ch2br2
      count = count+1
      sBase(count) = ch2brcl
      count = count+1
      sBase(count) = chbrcl2
      count = count+1
      sBase(count) = chbr2cl
      count = count+1
      sBase(count) = hfc23
      count = count+1
      sBase(count) = hfc32
      count = count+1
      sBase(count) = hfc125
      count = count+1
      sBase(count) = hfc134a
      count = count+1
      sBase(count) = hfc143a
      count = count+1
      sBase(count) = hfc152a
      count = count+1
      sBase(count) = co2b
      count = count+1
      sBase(count) = sf6b
      count = count+1
      sBase(count) = o3
      count = count+1
      sBase(count) = o3p
      count = count+1
      sBase(count) = o1d
      count = count+1
      sBase(count) = n
      count = count+1
      sBase(count) = no
      count = count+1
      sBase(count) = no2
      count = count+1
      sBase(count) = no3
      count = count+1
      sBase(count) = h
      count = count+1
      sBase(count) = oh
      count = count+1
      sBase(count) = ho2
      count = count+1
      sBase(count) = cl
      count = count+1
      sBase(count) = clo
      count = count+1
      sBase(count) = bro
      count = count+1
      sBase(count) = br
      count = count+1
      sBase(count) = cl2o2
      count = count+1
      sBase(count) = ch2o
      count = count+1
      sBase(count) = ch3o2
      count = count+1
      sBase(count) = ro3ox

      RETURN
      END SUBROUTINE storeBaseChem
