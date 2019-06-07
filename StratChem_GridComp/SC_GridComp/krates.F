      SUBROUTINE krates(gcSC,nl,m,kel,ak)
! =======================================================================

! ROUTINE
!   krates - new version for JPL 2006 - online calculation of rates 

! DESCRIPTION
!   Calculates and returns rate constants for thermal reactions

! ARGUMENTS
!  INPUT
!   kel : Temperature (K)
!     m : Number density

!  FROM the grid component
!   gcSC%cnsttab      1st column - k
!                2nd column - activation energy
!		 3rd column - k1
!		 4th column - e1
!		 5th - k2
!		 6th - e2
!  gcSC%indxs1-3     : 1,2,3 body indices of gcSC%cnsttab table

!  OUTPUT
!   ak	    : cm3 molecule-1 s-1 - rate constants

! -----------------------------------------------------------------
! JPL 2006: Numbers 1, 8, 31, 37, 70, and 108 are special cases.
! -----------------------------------------------------------------

! Revision history:
!             Randy Kawa   /home/kawa/tchem/fort/krate_calc_code.f
!  2 Jun 2009 Eric Nielsen Mods for including in GEOS-5, include Mohan
!                          Gupta's three reactions (123-125).
! =======================================================================
      USE SC_GridCompMod
      IMPLICIT NONE

      TYPE(SC_GridComp), INTENT(INOUT) :: gcSC  ! SC Grid Component

      INTEGER, INTENT(IN) :: nl
      REAL, INTENT(IN) :: kel(nl),m(nl)

      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      REAL(KIND=DBL), INTENT(OUT) :: ak(nl,gcSC%numreacs)

      REAL :: dum0,dum1,dum2,dum3,dum4

      INTEGER :: i,k

      OnLayer: DO k = 1,nl

! Set rates for 1st order reactions
! ---------------------------------
       DO i = 1,gcSC%num1d
        ak(k,gcSC%indxs1(i)) = gcSC%cnsttab(1,gcSC%indxs1(i))
       END DO

! 2-body reactions
! ----------------
       DO i = 1,gcSC%num2d
        ak(k,gcSC%indxs2(i)) = gcSC%cnsttab(1,gcSC%indxs2(i))*EXP(-gcSC%cnsttab(2,gcSC%indxs2(i))/kel(k))
       END DO

! Termolecular reaction rates
! ---------------------------
       DO i = 1,gcSC%num3d
        dum0 = 0.0d0
        dum1 = 0.0d0
        dum2 = 0.0d0
        dum3 = 0.0d0
        dum4 = 0.0d0
        dum0 = gcSC%cnsttab(1,gcSC%indxs3(i))*(kel(k)/300.0d0)**(-gcSC%cnsttab(2,gcSC%indxs3(i))) 
        dum1 = gcSC%cnsttab(3,gcSC%indxs3(i))*(kel(k)/300.0d0)**(-gcSC%cnsttab(4,gcSC%indxs3(i)))
        dum2 = dum0*m(k)/dum1
        dum4 = LOG10(dum2)
        dum3 = 0.6d0**(1./(1.+dum4*dum4))
        ak(k,gcSC%indxs3(i))=dum0/(1.+dum2)*dum3
       END DO	    

! ---------------------
! Two-body SPECIAL CASE
! ---------------------

!   8 O(1D) + M -> O(3P) + M  M= N2 O2
! ------------------------------------
       dum0 = gcSC%cnsttab(1,8)*EXP(-gcSC%cnsttab(2,8)/kel(k))
       ak(k,8) = 0.8d0*dum0+0.2d0*gcSC%cnsttab(3,8)*EXP(-gcSC%cnsttab(4,8)/kel(k))

! --------------------------
! Termolecular special cases
! --------------------------

!   1 O + O2 + M = O3 + M    
! -----------------------
       dum0=gcSC%cnsttab(1,1)*(kel(k)/300.0d0)**(-gcSC%cnsttab(2,1)) 
       dum1=1.0d0
       dum2=dum0*m(k)/dum1
       dum4=LOG10(dum2)
       dum3=0.6d0**(1./(1.+dum4*dum4))
       ak(k,1)=dum0/(1.+dum2)*dum3

!  36 CO + OH = CO2 + H became special 3-body chemical activation with JPL 06.
!  -- Old rate
!     No need to divide by M here since still used as 2-body in solver
! ----------------------------------------------------------------------------
!SRK   dum0 = (28.97d0*2870000.d0*kel(k)*1.66d-24)/1.013d+6
!SRK   ak(k,36) = gcSC%cnsttab(1,36)*(1.d0 + 0.6d0 *(m(k) * dum0)) 
!  -- New rate - Qing Liang 2016/02/04

!.... JPL 06-2 ; now CO + OH is composed of two separate reactions!
!.... Now it's density and temperature dependent.
!         M
! OH + CO -> HOCO , but HOCO + O2 -> HO2 + CO2 quickly ; termolecular
!         M
! OH + CO -> H + CO2, but H + O2 -> HO2 quickly ; chemical activation reaction
!       dum0=1.5d-13*(300.0d0/kel(k))**(-0.6)
!       dum1=2.1d9*(300.0d0/kel(k))**(-6.1)
!       dum2=dum0*m(k)/dum1
!       dum3=(dum0/(1.0d0+dum2))*0.6d0**(1.0d0/(1.0d0+(log10(dum2))**2))
!       dum4=dum3

!       dum0=m(k)*5.9d-33*(300.0d0/kel(k))**(1.4)
!       dum1=1.1d-12*(300.0d0/kel(k))**(-1.3)
!       dum2=dum0/dum1
!       dum3=(dum0/(1.0d0+dum2))*0.6d0**(1.0d0/(1.0d0+(log10(dum2))**2)) 
!       ak(k,36)=dum4+dum3       


!  37 HNO3 + OH = NO3 + H2O
! -------------------------
       dum0 = gcSC%cnsttab(1,37)*EXP(gcSC%cnsttab(2,37)/kel(k))
       dum1 = gcSC%cnsttab(3,37)*EXP(gcSC%cnsttab(4,37)/kel(k))
       dum2 = gcSC%cnsttab(5,37)*EXP(gcSC%cnsttab(6,37)/kel(k))*m(k)
       ak(k,37) = dum0+dum1/(1.+(dum1/dum2))

!  31 N2O5 + M = NO2 + NO3 + M
! ----------------------------
       ak(k,31)= ak(k,46)/((gcSC%cnsttab(1,31))*EXP(gcSC%cnsttab(2,31)/kel(k)))

!  70 HO2NO2 + M = HO2 + NO2 + M
! ------------------------------
       ak(k,70) = ak(k,34)/((gcSC%cnsttab(1,70))*EXP(gcSC%cnsttab(2,70)/kel(k)))

! 108 CL2O2 + M = CLO + CLO + M 
! -----------------------------
       ak(k,108) = ak(k,107)/((gcSC%cnsttab(1,108))*EXP(gcSC%cnsttab(2,108)/kel(k)))

      END DO OnLayer

      RETURN
      END SUBROUTINE krates
