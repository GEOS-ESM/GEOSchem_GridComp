      SUBROUTINE rdkrate(gcSC,fileName,verbose,rc)
! =======================================================================
! Read the ASCII file that contains the thermal reaction rate information.
! The calling routine is responsible for opening and closing the file.

! The file is read twice, the first scan is for detemining the number of 
! equations and the size of the allocatable arrays.  The second scan, after
! rewinding, of course, is for reading and storing the tabular data.

! FROM MODULES
!   gcSC%cnsttab - Table to fill in when reading in constants
!   gcSC%numreacs - Number of reactions in the list

! REVISION HISTORY:
!   Eric Nielsen June 2009
!   Eric Nielsen April 2012: Stand-alone version
! =======================================================================

      USE SC_GridCompMod
      IMPLICIT NONE

      TYPE(SC_GridComp), INTENT(INOUT) :: gcSC  ! SC Grid Component

      CHARACTER(LEN=*), INTENT(IN) :: fileName
      LOGICAL, INTENT(IN) :: verbose
      INTEGER, INTENT(OUT) :: rc

      INTEGER :: i, iostat, m
      INTEGER :: n, n1, n2, n3, nSCan
      INTEGER :: unitNumber
      LOGICAL :: exists, open, found, s
      CHARACTER(LEN=38) :: string
      CHARACTER(LEN=3) :: chn
      REAL :: r(6)

      rc = 0

! Find an available logical unit 
! ------------------------------
      found = .FALSE.
      i = 11

      DO WHILE (.NOT. found .AND. i <= 99)
       INQUIRE(UNIT=i, EXIST=exists, OPENED=open)
       IF(exists .AND. .NOT. open) THEN
        found = .TRUE.
        unitNumber = i
       END IF
       i = i+1
      END DO

      IF(.NOT. found) THEN
       WRITE(*,FMT="(/,'SC::rdKRates: No available logical units')")
       rc = 1
       RETURN
      END IF

! Open the ASCII file
! -------------------
      OPEN(UNIT=unitNumber, FILE=TRIM(fileName), FORM="formatted", ACTION="read", STATUS="old")

! First scan is for determining size of allocations
! -------------------------------------------------
      Scan: DO nScan = 1,2

       n  = 0
       n1 = 0
       n2 = 0
       n3 = 0

! For each reaction until end-of-file ...
! ---------------------------------------
       Rxns: DO

! Read next line of file into character string
! --------------------------------------------
        READ(unitNumber, FMT="(I3,1X,I1,1X,L1,1X,A38,3(ES10.3,F9.2))",
     &        IOSTAT=iostat) i,m,s,string,r

! Error exit
! ----------
        IF(iostat /= 0) EXIT

! Equation number
! ---------------
        n = n+1

! Save reaction data on second scan
! ---------------------------------
        IF(nScan == 2) THEN
         gcSC%KRxnName(n) =  TRIM(string)
         gcSC%KRxnOrder(n) = m
         gcSC%standardKRxn(n) = s
	 gcSC%cnsttab(1:6,n) = r(1:6)
	END IF

! Indicies tables
! ---------------
        SELECT CASE (m)
         CASE (1)
	  IF(s) THEN
	   n1 = n1 + 1
           IF(nScan == 2) gcSC%indxs1(n1) =  n
	  END IF
         CASE (2)
	  IF(s) THEN
	   n2 = n2 + 1
           IF(nScan == 2) gcSC%indxs2(n2) =  n
	  END IF
         CASE (3)
	  IF(s) THEN
	   n3 = n3 + 1
           IF(nScan == 2) gcSC%indxs3(n3) =  n
	  END IF
        END SELECT

       END DO Rxns

! After scanning the file once, do the allocations and rewind.
! ------------------------------------------------------------
       SetSizes: IF(nScan == 1) THEN

        gcSC%numreacs = n
        gcSC%num1d = n1
        gcSC%num2d = n2
        gcSC%num3d = n3

        ALLOCATE(gcSC%KRxnOrder(gcSC%numreacs))
        ALLOCATE(gcSC%standardKRxn(gcSC%numreacs))
        ALLOCATE(gcSC%KRxnName(gcSC%numreacs))
        ALLOCATE(gcSC%cnsttab(6,gcSC%numreacs))

        ALLOCATE(gcSC%indxs1(n1))
        ALLOCATE(gcSC%indxs2(n2))
        ALLOCATE(gcSC%indxs3(n3))

        gcSC%didKRxnsAlloc = .TRUE.

        REWIND(UNIT=unitNumber)

       END IF SetSizes

      END DO Scan

! Close the file
! --------------
      CLOSE(UNIT=unitNumber)

! Provide summary
! ---------------
      KInfo: IF(verbose) THEN

       WRITE(*,FMT="(' ')")
       WRITE(*,FMT="(' Thermal reactions: ',I3)") gcSC%numreacs
       WRITE(*,FMT="(' ')")
       WRITE(*,FMT="(' Rxn#  Standard  Order  Index  Reaction')")
       WRITE(*,FMT="(' ----  --------  -----  -----  ----------------')")

       n1 = 0
       n2 = 0
       n3 = 0

       DO n = 1,gcSC%numreacs

        chn="   "
        IF(gcSC%standardKRxn(n)) THEN

         SELECT CASE (gcSC%KRxnOrder(n))
          CASE (1)
           n1 = n1 + 1
           WRITE(chn,FMT="(I3)") n1
          CASE (2)
           n2 = n2 + 1
           WRITE(chn,FMT="(I3)") n2
          CASE (3)
           n3 = n3 + 1
           WRITE(chn,FMT="(I3)") n3
         END SELECT

        END IF
        WRITE(*,FMT="(I5,2X,L4,8X,I1,5X,A3,3X,A)") n,gcSC%standardKRxn(n),
     &  		      gcSC%KRxnOrder(n),chn,TRIM(gcSC%KRxnName(n))
       END DO
 
       WRITE(*,FMT="(' ')")

      END IF KInfo

      RETURN
      END SUBROUTINE rdkrate
