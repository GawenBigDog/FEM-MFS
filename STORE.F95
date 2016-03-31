program STORE
    implicit none



!--------------------------------------------------
!   INTERFACE FOR SUBROUTINE LOADD
!--------------------------------------------------
    TYPE LoadType
          INTEGER,DIMENSION(NLOAD) :: IDIRN=0
          REAL(8),DIMENSION(NLOAD) :: FLOAD=0
          INTEGER,DIMENSION(NLOAD) :: NFUNC=0
          REAL(8),DIMENSION(NEQ)   :: R=0
          INTEGER,DIMENSION(NEQ)   :: NRF=0
    END TYPE LoadType
    TYPE(LoadType) :: LOAD
!--------------------------------------------------
!   INTERFACE FOR SUBROUTINE INITLD
!--------------------------------------------------
    TYPE IncType
          REAL(8),DIMENSION(3,NUMNP) :: UI=0
          REAL(8),DIMENSION(3,NUMNP) :: VI=0
          REAL(8),DIMENSION(NEQ)     :: U=0
          REAL(8),DIMENSION(NEQ)     :: V=0
    END TYPE IncType
    TYPE(IncType) :: INCON
!--------------------------------------------------
!   INTERFACE FOR SUBROUTINE MATERIALD
!--------------------------------------------------
    REAL(8),DIMENSION(7,NMAT) :: MATP=0 
!--------------------------------------------------
!   INTERFACE FOR SUBROUTINE MATERIALD
!--------------------------------------------------


end program STORE
!
    SUBROUTINE LOADD
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
! *     SUBROUTINE LOADD 
! *      READ AND STORE LOAD DATA                            
! *      INPUT: 
! *       LOAD CASE
! *         IDIRN(NLOAD): LOAD DIRECTION CODES (1:X DIRECTION,2:Y DIRECTION; 3: Z DIRECTION)
! *         FLOAD(NLOAD): MAGNITUDE OF LOAD
! *         NFUNC(NLOAD): LOAD FUNCTION CODES 
! *               1: CONSTANT
! *               2: LINEAR
! *               3: TRIGONOMETRIC
! *         R(NEQ)         : NODAL LOAD DATA
! *         NRFUNC(NEQ)    : LOAD FUNCTION CODES                               
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
       IMPLICIT NONE
       INTEGER :: I,J
!--------------------------------------------------
!      Internal variables
!--------------------------------------------------
       INTEGER :: N
       INTEGER :: IDIRNN,NFUNCC
       REAL(8) :: FLOADD
!--------------------------------------------------
!      Assistant variables
!--------------------------------------------------
       INTEGER :: NumNEQ=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       WRITE (IOUT,2015) PROJECT%NLOAD
       WRITE (IOUT,2020)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       N = 0
       I=0
       READ (IIN,*) N,IDIRNN,FLOADD,NFUNCC
       DO WHILE(N.NE.0)
         IF(N.GT.PROJECT%NUMNP) THEN
           WRITE(*,*) "ERROR,NODE NUMBER .GT. NUMNP IN LOAD_DATA!"
           STOP
         END IF
         I=I+1        
         IF(I.GT.PROJECT%NLOAD) THEN
            WRITE(*,*) "ERROR, NUMBER OF INPUTTED LOADS .GT. NLOAD!"
            STOP
         END IF
         LOAD%IDIRN(I)=IDIRNN
         LOAD%FLOAD(I)=FLOADD
         LOAD%NFUNC(I)=NFUNCC        
         WRITE (IOUT,2030) N,IDIRNN,FLOADD,NFUNCC
!
         NumNEQ=NODE%ID(LOAD%IDIRN(I),N)
         IF (NumNEQ.NE.0) THEN
            LOAD%R(NumNEQ)=LOAD%FLOAD(I)
            LOAD%NRF(NumNEQ)=LOAD%NFUNC(I)
         END IF
!
         READ (IIN,*) N,IDIRNN,FLOADD,NFUNCC
       END DO
       IF(I.NE.PROJECT%NLOAD) THEN
          WRITE(*,*) "ERROR, NUMBER OF INPUTTED LOADS .LT. NLOAD!"
          STOP
       END IF
!
       WRITE(IOUT,2040)
       DO I=1,PROJECT%NEQ
          WRITE(IOUT,2050) I,LOAD%R(I),LOAD%NRF(I)  
       END DO
!
       RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 2015  FORMAT(/,'NUMBER OF LOADS =',I5,/)
 2020  FORMAT(' NODE         IDRIN         LOAD           NFUNC',/,&
              'NUMBER        NUMBER      MAGNITUDE        NUMBER')
 2030  FORMAT(I5,9X,I5,3X,F13.3,5X,I6)
 2040  FORMAT(/,'EQUATION        LOAD           FUNCTION',/&
              ' NUMBER       MAGNITUDE          CODES ')
 2050  FORMAT(I5,5X,F13.3,7X,I6)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!}}}
    END SUBROUTINE LOADD
!
    SUBROUTINE INITLD
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
! *     SUBROUTINE INITLDT
! *      READ AND STORE INITIAL CONDITIONS                           
! *      INPUT: 
! *        UI(3,NUMNP)    : X Y Z DISPLACEMENT
! *	   VI(3,NUMNP)    : X Y Z VELOCITY
! *        U(NEQ)         : DISPLACEMENT VARIABLES
! *        V(NEQ)         : VELOCITY VARIABLES
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
       IMPLICIT NONE
       INTEGER :: I,N,K
       REAL(8) :: UXX,UYY,UZZ,VXX,VYY,VZZ
!--------------------------------------------------
!       Assistant variables
!--------------------------------------------------
       INTEGER :: NCOUNT
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       WRITE (IOUT,2000)
       WRITE (IOUT,2010)
       WRITE (IOUT,2020)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
       N = 0
       NCOUNT=0
       READ (IIN,*) N,UXX,UYY,UZZ,VXX,VYY,VZZ     
       DO WHILE (N.NE.0)
          IF(N.GT.PROJECT%NUMNP) THEN
             WRITE(*,*) "ERROR, NODE NUMBER .GT. NUMNP!"
             STOP
          END IF
          INCON%UI(1,N)=UXX
          INCON%UI(2,N)=UYY
          INCON%UI(3,N)=UZZ
          INCON%VI(1,N)=VXX
          INCON%VI(2,N)=VYY
          INCON%VI(3,N)=VZZ  
          NCOUNT=NCOUNT+1
          READ (IIN,*) N,UXX,UYY,UZZ,VXX,VYY,VZZ     
       END DO
       IF(NCOUNT.NE.PROJECT%NUMNP) THEN
          WRITE(*,*) "ERROR,NUMBER OF INPUTTED INITIAL DATA .LT. NUMNP!"
          STOP
       END IF
!----------------write out initial conditions
       DO N=1,PROJECT%NUMNP
          WRITE (IOUT,2030) N,INCON%UI(1,N),INCON%UI(2,N),INCON%UI(3,N),INCON%VI(1,N),INCON%VI(2,N),&
                            INCON%VI(3,N)
       END DO
!----------------data transform
       DO N=1,PROJECT%NUMNP
          DO I=1,3
            IF(NODE%ID(I,N).NE.0) THEN
              K=NODE%ID(I,N)
              INCON%U(K)=INCON%UI(I,N)
              INCON%V(K)=INCON%VI(I,N)
            END IF
          END DO
       END DO
!----------------Write new Initial conditions
       WRITE(IOUT,2040)
       DO I=1,PROJECT%NEQ
          WRITE(IOUT,2050) I,INCON%U(I),INCON%V(I)  
       END DO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 2000  FORMAT(//,' N O D A L   I N I T I A L   D A T A',/)
 2010  FORMAT(' INPUT NODAL DATA',//)
 2020  FORMAT(' NODE',14X,'UX',11X,'UY',11X,'UZ',11X,'VX',11X,'VY',11X,&
              'VZ',/,' NUMBER')
 2030  FORMAT (I5,6X,6F13.3)
 2040  FORMAT(/,'EQUATION     INITIAL_DISP           INITIAL_VEL',/&
               ' NUMBER                                        ')
 2050  FORMAT(I5,5X,F13.3,7X,F13.3)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!}}}
    END SUBROUTINE INITLD  
!
    SUBROUTINE MATERIALD
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
! *     SUBROUTINE MATERIALD
! *      READ AND STORE MATERIAL DATA                           
! *      INPUT: 
! *        MATP(7,NMAT)  :  MATERIAL PROPERTIES
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
        IMPLICIT NONE
! Internal variables        
        INTEGER :: N=0
        REAL(8) :: E=0,PR=0
! Assistant variables
        INTEGER :: NCOUNT
!
        NCOUNT=0
!
        WRITE (IOUT,2000)
        WRITE (IOUT,2010)
        WRITE (IOUT,2020)
!
        READ(IIN,*) N,E,PR
        DO WHILE(N.NE.0)
           IF(N.GT.PROJECT%NMAT) THEN
                WRITE(*,*) "ERROR, MATERIAL SET NUMBER .GT. NMAT"
                STOP
           END IF
           NCOUNT=NCOUNT+1
           MATP(1,N)=E
           MATP(2,N)=PR
           WRITE(IOUT,2030) N, E, PR
           READ(IIN,*) N,E,PR
        END DO
        IF(NCOUNT.NE.PROJECT%NMAT) THEN
            WRITE(*,*) "ERROR, NUMBER OF MATERIAL SET .LT. NMAT"
            STOP
        RETURN
!
 2000  FORMAT(//,' M A T E R I A L    D A T A',/)
 2010  FORMAT(' INPUT MATERIAL DATA',//)
 2020  FORMAT(' SET ',10X,' E ',10X,' PR ',/,&
              'NUMBER')
 2030  FORMAT(I6,10X,F13.3,F13.3)
!}}}
    END SUBROUTINE MATERIALD
            
!--------------------------------------------------
!   DEALLOCATE WAREHOUSE FOR SUBROUTINE NODED
!--------------------------------------------------
         DEALLOCATE (ID)
         DEALLOCATE (X)
         DEALLOCATE (Y)
         DEALLOCATE (Z)
!--------------------------------------------------
!   DEALLOCATE WAREHOUSE FOR SUBROUTINE LOADD
!--------------------------------------------------
         DEALLOCATE (IDIRN)
         DEALLOCATE (FLOAD)
         DEALLOCATE (NFUNC)
         DEALLOCATE (R)
         DEALLOCATE (NRF)
!--------------------------------------------------
!   DEALLOCATE WAREHOUSE FOR SUBROUTINE INITLD
!--------------------------------------------------
         DEALLOCATE (UI)
         DEALLOCATE (VI)
         DEALLOCATE (U)
         DEALLOCATE (V)
!--------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE MATERIALD
!--------------------------------------------------
         DEALLOCATE (MATP)       






       
