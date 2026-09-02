SUBROUTINE localZeroSet(ro3po3,ro1do3,ro3ox,rnno,rnono2,rno3no2,rno2nox,rclclo, &
	  rcloclx,rocloclo,rcl2o2clo,rbrobrx,rbrnbrx,po3,lo3,ph2o2,lh2o2,phno3, &
	  lhno3,phno3het,pn2o5,ln2o5,ln2o5het,pho2no2,lho2no2,pnox,lnox,lnoxsq, &
	    pnoxcl,lnoxcl,pclx,lclx,pclono2,lclono2,lclono2het,phcl,lhcl,phocl, &
	    lhocl,pbrono2,lbrono2het,pbrx,lbrx,phbr,lhbr,phobr,lhobr,pno2,lno2, &
	    pno3,lno3,poclo,pnoxa,lnoxa,pcl2,pbrcl,lo3hox,lo3nox,lo3oxsq,lo3cl, &
	   lo3cly,lo3oh,lo3brx,lo3chx,o3pe,o3e,oxe,ne,noe,no2e,no3e,noxe,n2o5e, &
	     h2o2e,hno3e,ho2no2e,hcle,hocle,cle,cloe,ocloe,cl2o2e,clono2e,clxe, &
					      bre,broe,brono2e,brxe,hbre,hobre, &
					 rmedsts,rmednat,rmedice,denssts,vfall)
 IMPLICIT NONE
 INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

 REAL(KIND=DBL), INTENT(OUT) :: ro3po3,ro1do3
 REAL(KIND=DBL), INTENT(OUT) :: ro3ox,rnno
 REAL(KIND=DBL), INTENT(OUT) :: rnono2,rno3no2
 REAL(KIND=DBL), INTENT(OUT) :: rno2nox,rclclo
 REAL(KIND=DBL), INTENT(OUT) :: rcloclx,rocloclo
 REAL(KIND=DBL), INTENT(OUT) :: rcl2o2clo,rbrobrx
 REAL(KIND=DBL), INTENT(OUT) :: rbrnbrx,po3
 REAL(KIND=DBL), INTENT(OUT) :: lo3,ph2o2
 REAL(KIND=DBL), INTENT(OUT) :: lh2o2,phno3
 REAL(KIND=DBL), INTENT(OUT) :: lhno3,phno3het
 REAL(KIND=DBL), INTENT(OUT) :: pn2o5,ln2o5
 REAL(KIND=DBL), INTENT(OUT) :: ln2o5het,pho2no2
 REAL(KIND=DBL), INTENT(OUT) :: lho2no2,pnox
 REAL(KIND=DBL), INTENT(OUT) :: lnox,lnoxsq
 REAL(KIND=DBL), INTENT(OUT) :: pnoxcl,lnoxcl
 REAL(KIND=DBL), INTENT(OUT) :: pclx,lclx
 REAL(KIND=DBL), INTENT(OUT) :: pclono2,lclono2
 REAL(KIND=DBL), INTENT(OUT) :: lclono2het,phcl
 REAL(KIND=DBL), INTENT(OUT) :: lhcl,phocl
 REAL(KIND=DBL), INTENT(OUT) :: lhocl,pbrono2
 REAL(KIND=DBL), INTENT(OUT) :: lbrono2het
 REAL(KIND=DBL), INTENT(OUT) :: pbrx,lbrx
 REAL(KIND=DBL), INTENT(OUT) :: phbr,lhbr
 REAL(KIND=DBL), INTENT(OUT) :: phobr,lhobr
 REAL(KIND=DBL), INTENT(OUT) :: pno2,lno2
 REAL(KIND=DBL), INTENT(OUT) :: pno3,lno3
 REAL(KIND=DBL), INTENT(OUT) :: poclo,pnoxa
 REAL(KIND=DBL), INTENT(OUT) :: lnoxa,pcl2
 REAL(KIND=DBL), INTENT(OUT) :: pbrcl,lo3hox
 REAL(KIND=DBL), INTENT(OUT) :: lo3nox,lo3oxsq
 REAL(KIND=DBL), INTENT(OUT) :: lo3cl,lo3cly
 REAL(KIND=DBL), INTENT(OUT) :: lo3oh,lo3brx
 REAL(KIND=DBL), INTENT(OUT) :: lo3chx
 REAL(KIND=DBL), INTENT(OUT) :: o3pe,o3e,oxe,ne,noe,no2e,no3e
 REAL(KIND=DBL), INTENT(OUT) :: noxe,n2o5e,h2o2e,hno3e,ho2no2e,hcle,hocle
 REAL(KIND=DBL), INTENT(OUT) :: cle,cloe,ocloe,cl2o2e,clono2e,clxe,bre,broe
 REAL(KIND=DBL), INTENT(OUT) :: brono2e,brxe,hbre,hobre

 REAL, INTENT(OUT) :: rmedsts,rmednat,rmedice,denssts,vfall

! Ratios
! ------
     ro3po3 = 0.0
     ro1do3 = 0.0
      ro3ox = 0.0
       rnno = 0.0
     rnono2 = 0.0
    rno3no2 = 0.0
    rno2nox = 0.0
     rclclo = 0.0
    rcloclx = 0.0
   rocloclo = 0.0
  rcl2o2clo = 0.0
    rbrobrx = 0.0
    rbrnbrx = 0.0

! Production and loss
! -------------------
	po3 = 0.0
	lo3 = 0.0
      ph2o2 = 0.0
      lh2o2 = 0.0
      phno3 = 0.0
      lhno3 = 0.0
   phno3het = 0.0
      pn2o5 = 0.0
      ln2o5 = 0.0
   ln2o5het = 0.0
    pho2no2 = 0.0
    lho2no2 = 0.0
       pnox = 0.0
       lnox = 0.0
     lnoxsq = 0.0
     pnoxcl = 0.0
     lnoxcl = 0.0
       pclx = 0.0
       lclx = 0.0
    pclono2 = 0.0
    lclono2 = 0.0
 lclono2het = 0.0
       phcl = 0.0
       lhcl = 0.0
      phocl = 0.0
      lhocl = 0.0
    pbrono2 = 0.0
 lbrono2het = 0.0
       pbrx = 0.0
       lbrx = 0.0
       phbr = 0.0
       lhbr = 0.0
      phobr = 0.0
      lhobr = 0.0
       pno2 = 0.0
       lno2 = 0.0
       pno3 = 0.0
       lno3 = 0.0
      poclo = 0.0
      pnoxa = 0.0
      lnoxa = 0.0
       pcl2 = 0.0
      pbrcl = 0.0
     lo3hox = 0.0
     lo3nox = 0.0
    lo3oxsq = 0.0
      lo3cl = 0.0
     lo3cly = 0.0
      lo3oh = 0.0
     lo3brx = 0.0
     lo3chx = 0.0

! Estimated species
! -----------------
       o3pe = 0.0
        o3e = 0.0
        oxe = 0.0
  	 ne = 0.0
        noe = 0.0
       no2e = 0.0
       no3e = 0.0
       noxe = 0.0
      n2o5e = 0.0
      h2o2e = 0.0
      hno3e = 0.0
    ho2no2e = 0.0
       hcle = 0.0
      hocle = 0.0
        cle = 0.0
       cloe = 0.0
      ocloe = 0.0
     cl2o2e = 0.0
    clono2e = 0.0
       clxe = 0.0
        bre = 0.0
       broe = 0.0
    brono2e = 0.0
       brxe = 0.0
       hbre = 0.0
      hobre = 0.0

! PSC profile data
! ----------------

    rmedsts = 0.0
    rmednat = 0.0
    rmedice = 0.0
    denssts = 0.0
      vfall = 0.0

 RETURN
END SUBROUTINE localZeroSet

