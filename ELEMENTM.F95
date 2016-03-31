MODULE ELEMENTM
    USE WAREHOUSE
    USE INTERPRETER
    IMPLICIT NONE
    PUBLIC  :: ELEMENT_IN
    PRIVATE :: TRUSSD
!    
    CONTAINS
!**************************************************
!
!   SUBROUTINES FOR ELEMENT INPUT
!
!**************************************************
    SUBROUTINE ELEMENT_IN(IIN,IOUT,ITEM1)
!{{{
!--------------------------------------------------
!   SELECT ELEMENT TYPE
!   NPAR()
!   TRUSS
!       NPAR(1) = 1
!       NPAR(2) = MTYPE
!       NPAR(3) = NUMELE
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT,ITEM1
!
    INTEGER :: I,J,L,M
    INTEGER :: NumGroup
    INTEGER,PARAMETER :: KwEle=1
    CHARACTER(20),DIMENSION(KwEle) :: KwElement
!       
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KwElement(1)='\TRUSS'
!
    WRITE(IOUT,2000)
!
    DO J=1,NEGROUP
       I = KeyWord(KwElement,KwEle)
       WRITE(IOUT,1000) KwElement(I)
       SELECT CASE(I)
       CASE(1)
          NPAR(1)=1
          CALL TRUSSD(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,7)
          WRITE(ITEM1) ((XYZ(L,M),L=1,6),M=1,NPAR(3)),((LM(L,M),L=1,6),M=1,NPAR(3))
          ! Deallocate variables
          DEALLOCATE (XYZ)
          DEALLOCATE (LM)
       CASE DEFAULT
          WRITE(*,*) "ERROR, INPUT ERROR IN ELEMENTD"
          STOP
       END SELECT 
    END DO
    RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 1000   FORMAT(/,A20)
 2000 FORMAT(//,' E L E M E N T   I N F O R M A T I O N',///)
!}}}  
    END SUBROUTINE ELEMENT_IN
!
    SUBROUTINE TRUSSD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE TRUSSD
! *      READ AND STORE TRUSS DATA
! *      TRUSS
! *         NPAR(1) = 1
! *         NPAR(2) = MTYPE
! *         NPAR(3) = NUMELE
! *      INPUT:
! *        MTYPE    -  MATERIAL SET NUMBER
! *    	   NUMELE   -  THE NUMBER OF TRUSS ELEMENTS  
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J
    INTEGER :: NumEle
    INTEGER :: N,II,JJ
    INTEGER :: NCOUNT
    INTEGER :: AllocateStatus
!
    CHARACTER(20),DIMENSION(3) :: KWTRUSS
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWTRUSS(1) = 'MTYPE:'
    KWTRUSS(2) = 'NUMELE:'
    KWTRUSS(3) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=0
!   
    I = KeyWord(KWTRUSS,3)
    DO WHILE(I<=3)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           NPAR(3)=J
         CASE(3)
           EXIT
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN TRUSSD"
           STOP
       END SELECT
       I = KeyWord(KWTRUSS,4)
    END DO
!
    WRITE(IOUT,1000) NPAR(2),NPAR(3)
!
    IF(NPAR(3).EQ.0) THEN
            WRITE(*,*) "ERROR, NUMBER OF TRUSS ELEMENT .EQ.0"
            STOP
    END IF
    NumEle=NPAR(3)
!
!-------- Allocate memory
!
        ALLOCATE (XYZ(6,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
        ALLOCATE (LM(6,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!
        WRITE(IOUT,2000)
        NCOUNT=0
        READ(IIN,*) N, II, JJ
        DO WHILE(N.NE.0)
           IF(N.GT.NumEle) THEN
              WRITE(*,*) "ERROR N .GT. NumEle"
              STOP
           END IF
           NCOUNT=NCOUNT+1
!
           XYZ(1,N)=X(II)
           XYZ(2,N)=Y(II)
           XYZ(3,N)=Z(II)
!
           XYZ(4,N)=X(JJ)
           XYZ(5,N)=Y(JJ)
           XYZ(6,N)=Z(JJ)
!           
           DO I=1,3
             LM(I,N)=ID(I,II)
             LM(I+3,N)=ID(I,JJ)
           END DO
           WRITE(IOUT,2050) N,II,JJ
!
           READ(IIN,*) N,II,JJ
        END DO
!
        IF(NCOUNT.NE.NumEle) THEN
              WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumEle"
              STOP
        END IF
    RETURN
!--------------------------------------------------
 1000 Format(/,&
              '    THE ELELMENT TYPE: TRUSS',/,&
              '    THE MATERIAL SET NUMBER          =',I5,/,&
              '    THE NUMBER OF TRUSS ELEMENTS     =',I5,/)
 2000 FORMAT(/,&
              ' ELEMENT     NODE     NODE       ',/,&
              ' NUMBER-N      I        J       ',/)
 2050 FORMAT(I5,6X,I5,4X,I5)
!}}}
    END SUBROUTINE TRUSSD
!**************************************************
!
!   SUBROUTINES FOR ELEMENT INPUT
!
!**************************************************
END MODULE ELEMENTM
