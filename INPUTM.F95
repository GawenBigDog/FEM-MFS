! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       .
! .                            DISCOVERY 0.1                              .
! .                                                                       .
! .    CONTACT DYNAMICS                                                   .
! .    MAIN PROGRAM                                                       .
! .                                                                       .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .
module INPUTM            
    USE INTERPRETER
    USE WAREHOUSE
    USE SOLVERM
    implicit none
    PUBLIC  :: INPUT    
    PRIVATE :: GLOBALD
    PRIVATE :: NODED
    PRIVATE :: BOUNDD
    PRIVATE :: LOADD
    PRIVATE :: INITLD
    PRIVATE :: MATERIALD
!
    PRIVATE :: BOUNDARY_IN
    PRIVATE :: LINED
!
    PRIVATE :: ELEMENT_IN
    PRIVATE :: TRUSSD
    PRIVATE :: MFS_TRUSS_ASSD
    PRIVATE :: QUADSD
    PRIVATE :: TRUSS_MFSD
    PRIVATE :: QUADS_MFSD
    PRIVATE :: MFS_2D_ASSD
!
    PRIVATE :: SURFACE_IN
    PRIVATE :: PTPTD
! 
!    PRIVATE :: SURFACE_IN
!    PRIVATE :: POINTD
    CONTAINS
!
    SUBROUTINE INPUT(IIN,IOUT,ITEM3,ITEM1,ITEM2)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *     INPUTM                                        
! *      READ AND STORE INPUT DATA                    
! *      INPUT:                                       
! *            1. GLOBAL   DATA                                 
! *            2. NODE     DATA                         
! *            3. LOAD     DATA  
! *            4. INITIAL  DATA
! *            5. MATERIAL DATA
! *            6. LINEAR ELEMENT DATA
! *            7. SURFACE  DATA
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
!
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: IIN,IOUT,ITEM3,ITEM1,ITEM2
!        
        INTEGER :: I
        INTEGER,PARAMETER :: KWNUM=10
        CHARACTER(LEN=20),DIMENSION(KWNUM) :: KWINPUT ! KEYWORDS OF DIFFERENT INPUT PHASE
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
        KWINPUT(1) = '\GLOBAL'
        KWINPUT(2) = '\NODE'
        KWINPUT(3) = '\LOAD'
        KWINPUT(4) = '\BOUND'
        KWINPUT(5) = '\INITIAL'
        KWINPUT(6) = '\MATERIAL'
        KWINPUT(7) = '\BOUNDARY'
        KWINPUT(8) = '\ELEMENT'
        KWINPUT(9) = '\SURFACE'
        KWINPUT(10) = '\END{INPUT}'     
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        I = KeyWord(KWINPUT,KWNUM)
        DO WHILE(I<=KWNUM)
           WRITE(IOUT,1000) KWINPUT(I)
           SELECT CASE(I)
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE GLOBAL DATA               *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(1)
              CALL GLOBALD(IIN,IOUT)
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE NODE   DATA               *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(2)
              CALL NODED(IIN,IOUT)
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE LOAD   DATA               *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(3)
              CALL LOADD(IIN,IOUT)
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE LOAD   DATA               *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(4)
              CALL BOUNDD(IIN,IOUT)
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE INITIAL DATA              *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(5)
              IND_INITIAL=1
              CALL INITLD(IIN,IOUT)
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE MATERIAL DATA             *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(6)
              CALL MATERIALD(IIN,IOUT)
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE BOUNDARY  DATA            *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(7)
              CALL BOUNDARY_IN(IIN,IOUT,ITEM3)
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE ELEMENT  DATA             *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(8)
              CALL ELEMENT_IN(IIN,IOUT,ITEM1)   
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE MATERIAL DATA             *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(9)
              IND_SURFACE=1
              CALL SURFACE_IN(IIN,IOUT,ITEM2) 
           CASE(10)
              RETURN          
           CASE DEFAULT
              WRITE(*,*) "ERROR, INPUT OPERATION PHASE DOES NOT EXIST!"
              STOP
           END SELECT
           I = KeyWord(KWINPUT,KWNUM)
        END DO
!
        RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 1000   FORMAT(/,A20)
!}}}
    END SUBROUTINE INPUT
!
    SUBROUTINE GLOBALD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE GLOBALD 
! *      READ AND STORE GLOBAL DATA                            
! *      INPUT:                                       
! *    	   HED    -   The master heading information for use in labeling the output
! *        NUMNP  -   Total number of nodal points
! *				0 : program stop
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT NONE 
        INTEGER,INTENT(IN)    :: IIN,IOUT
!
        INTEGER :: I,K
        INTEGER,PARAMETER :: NumKwd = 16
        CHARACTER(LEN=20),DIMENSION(NUMKWD) :: KWGLOBAL ! KEYWORDS OF DIFFERENT  PHASE 
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
        KWGLOBAL(1) = 'PROJECT:'
        KWGLOBAL(2) = 'NUMNP:'
        KWGLOBAL(3) = 'NLOAD:'
        KWGLOBAL(4) = 'NBOUND:'
        KWGLOBAL(5) = 'NBGROUP:'
        KWGLOBAL(6) = 'NMAT:'
        KWGLOBAL(7) = 'NEGROUP:'
        KWGLOBAL(8) = 'NSPAIR:'
        KWGLOBAL(9) = 'DYNAMIC:'
        KWGLOBAL(10) = 'CONTACT:'
        KWGLOBAL(11) = 'NITP_IN:'
        KWGLOBAL(12) = 'NIDIV_IN:'
        KWGLOBAL(13) = 'NITP_BC:'
        KWGLOBAL(14) = 'NIDIV_BC:'
        KWGLOBAL(15) = 'PENALTY:'
        KWGLOBAL(16) = '\END'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
        I = KeyWord(KWGLOBAL,NumKwd)
        DO WHILE(I<=NumKwd)
          SELECT CASE(I)
          CASE(1)
             CALL GetString(HEAD)
             CALL NewLine()
             WRITE(IOUT,1000) HEAD
          CASE(2)
             K = GetInt()
             NUMNP = K
          CASE(3)
             K = GetInt()
             NLOAD = K
          CASE(4)
             K = GetInt()
             NBOUND = K
          CASE(5)
             K = GetInt()
             NBGROUP = K
          CASE(6)
             K = GetInt()
             NMAT = K
          CASE(7)
             K = GetInt()
             NEGROUP = K
          CASE(8)
             K = GetInt()
             NSPAIR = K
          CASE(9)
             K = GetInt()
             DYNAMIC = K
          CASE(10)
             K = GetInt()
             CONTACT = K
          CASE(11)
             K = GetInt()
             NITP_IN=K
          CASE(12)
             K = GetInt()
             NIDIV_IN=K
          CASE(13)
             K = GetInt()
             NITP_BC=K
          CASE(14)
             K = GetInt()
             NIDIV_BC=K
          CASE(15)
             K = GetInt()
             PENALTY=K 
          CASE(16)
             EXIT
          CASE DEFAULT
             WRITE(*,*) "ERROR, UNEXPECTED GLOBAL DATA"
             STOP
          END SELECT
          I = KeyWord(KWGLOBAL,NumKwd)
        END DO 
        IF(NUMNP.EQ.0) THEN
           WRITE(*,*) "ERROR,NUMBER OF NODES .EQ. ZERO!"
           STOP
        END IF
!
        WRITE(IOUT,2000) NUMNP,NLOAD,NBOUND,NMAT,NBGROUP,NEGROUP,NSPAIR,PENALTY,DYNAMIC,CONTACT
        WRITE(IOUT,2010) NITP_IN,NIDIV_IN,NITP_BC,NIDIV_BC
!
        RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 1000   FORMAT(/,A40)
 2000   FORMAT(/,&          
          ' C O N T R O L   I N F O R M A T I O N',//,&
          '      NUMBER OF NODAL POINTS ',10(' .'),' (NUMNP)   = ',I5,//,&
          '      NUMBER OF LOAD     ',12(' .'),    ' (NLOAD)   = ',I5,//,&
          '      NUMBER OF EBC      ',12(' .'),    ' (NBOUND)  = ',I5,//,&
          '      NUMBER OF MATERIAL SET ',10(' .'),' (NMAT)    = ',I5,//,&
          '      NUMBER OF BOUNDARY GROUP',10(' .'),'(NBGROUP) = ',I5,//,&
          '      NUMBER OF ELEMENT GROUP',10(' .'),'(NEGROUP)  = ',I5,//,&
          '      NUMBER OF SURFACE PAIR ',10(' .'),'(NSPAIR)   = ',I5,//,&
          '      PENALTY FOR DIRICHLET BCS',10(' .'),'(PENALTY)= ',I5,//,&
          '                      0 --- PENALTY MODE OFF          ',/,&
          '                      1 --- PENALTY MODE ON           ',/,&
          '      DYNAMIC SIMULATION MODE',10(' .'),'(DYNAMIC)  = ',I5,//,&
          '                      0 --- DYNAMIC MODE OFF          ',/,&
          '                      1 --- DYNAMIC MODE ON           ',/,&
          '      CONTACT SIMULATION MODE',10(' .'),'(CONTACT)  = ',I5,//,&
          '                      0 --- CONTACT MODE OFF          ',/,&
          '                      1 --- CONTACT MODE ON           ',//)
 2010   FORMAT(/,&
          ' I N T E G R A T I O N   I N F O R M A T I O N (MFS)',//,&
          '      NUMBER OF GAUSS-LEGENDRE INTEGRATION POINTS (INNER) ',5(' .'),' (NITP_IN)   = ',I5,//,&
          '      NUMBER OF INTEGRATION SECTIONS (INNER)              ',5(' .'),' (NIDIV_IN)  = ',I5,//,&
          '      NUMBER OF GAUSS-LEGENDRE INTEGRATION POINTS (BC)    ',5(' .'),' (NITP_BC)   = ',I5,//,&
          '      NUMBER OF INTEGRATION SECTIONS (BC)                 ',5(' .'),' (NIDIV_BC)  = ',I5,//)  
!}}}
    END SUBROUTINE GLOBALD
!
    SUBROUTINE NODED(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *     SUBROUTINE NODED
! *     INPUT AND STORE NODE DATA                            
! *      INPUT:     
! *        ID(3,NUMNP) : Boundary condition codes (0=free,1=deleted)
! *        X(NUMNP)    : X coordinates
! *	   Y(NUMNP)    : Y coordinates
! *	   Z(NUMNP)    : Z coordinates
! *        NPOR(NUMNP) : Number of the polynomial order
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: IIN,IOUT
       INTEGER :: I, J
! Internal variables for node input data
       INTEGER :: N,KN
       INTEGER :: IDD(3)
       REAL(8) :: XX,YY,ZZ
       INTEGER :: NDOFX,NDOFY,NDOFZ
! Assistant variables
       INTEGER :: NCOUNT
       INTEGER :: KNOLD,NOLD
       REAL(8) :: XNUM, DX, DY, DZ
       INTEGER :: NUM, NUMN
       INTEGER :: K,KK
       INTEGER :: ALLocateStatus
!-------- Allocate memory
       ALLOCATE (ID(3,NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (X(NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (Y(NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (Z(NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (NDOF(3,NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! RESET NODE DATA
       DO I=1,NUMNP
          X(I)=0.0
          Y(I)=0.0
          Z(I)=0.0
          DO J=1,3
             ID(J,I)=0
             NDOF(J,I)=0
          END DO
       END DO  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       WRITE (IOUT,2000)
       WRITE (IOUT,2010)
       WRITE (IOUT,2020)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       NCOUNT=0                  !RECORD NUMBER OF INPUTTED NODES     
       KNOLD=0
       NOLD=0
!
       N = 0
       READ (IIN,*) N,(IDD(I),I=1,3),XX,YY,ZZ,NDOFX,NDOFY,NDOFZ,KN
       IF(NDOFX.GT.NPOR_MAX) NPOR_MAX=NDOFX
       IF(NDOFY.GT.NPOR_MAX) NPOR_MAX=NDOFY
       IF(NDOFZ.GT.NPOR_MAX) NPOR_MAX=NDOFZ 
!
       DO  WHILE(N.NE.0)  
          IF(N.GT.NUMNP) THEN
            WRITE(*,*) "ERROR,NODE NUMBER .GT. NUMNP!"
            STOP
          END IF
          ID(:,N)=IDD
          X(N)=XX
          Y(N)=YY
          Z(N)=ZZ
          NDOF(1,N)=NDOFX
          NDOF(2,N)=NDOFY
          NDOF(3,N)=NDOFZ
          NCOUNT=NCOUNT+1
          WRITE (IOUT,2030) N,(IDD(I),I=1,3),XX,YY,ZZ,NDOFX,NDOFY,NDOFZ,KN
            IF (KNOLD.NE.0) THEN
              NUM=(N-NOLD)/KNOLD
              NUMN=NUM - 1
              IF (NUMN.GE.1) THEN
                 XNUM=NUM
                 DX=(X(N)-X(NOLD))/XNUM
                 DY=(Y(N)-Y(NOLD))/XNUM
                 DZ=(Z(N)-Z(NOLD))/XNUM
                 K=NOLD
                 DO J=1,NUMN
                    KK=K
                    K=K + KNOLD
                    X(K) = X(KK) + DX
                    Y(K) = Y(KK) + DY
                    Z(K) = Z(KK) + DZ
                    NDOF(1,K)=NDOF(1,KK)
                    NDOF(2,K)=NDOF(2,KK)
                    NDOF(3,K)=NDOF(3,KK)
                    NCOUNT=NCOUNT+1
                    DO I=1,3
                       ID(I,K) = ID(I,KK)
                    END DO
                 END DO
              END IF
            ENDIF
!        
           NOLD=N
           KNOLD=KN
           READ (IIN,*) N,(IDD(I),I=1,3),XX,YY,ZZ,NDOFX,NDOFY,NDOFZ,KN
           IF(NDOFX.GT.NPOR_MAX) NPOR_MAX=NDOFX
           IF(NDOFY.GT.NPOR_MAX) NPOR_MAX=NDOFY
           IF(NDOFZ.GT.NPOR_MAX) NPOR_MAX=NDOFZ
!
       END DO
       IF(NCOUNT.NE.NUMNP) THEN
         WRITE(*,*) "ERROR, NUMBER OF NODES .LT. NUMNP!"
         STOP
       END IF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Write complete nodal data
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       WRITE (IOUT,2015)
       WRITE (IOUT,2020)
       DO N=1,NUMNP
         WRITE (IOUT,2030) N,(ID(I,N),I=1,3),X(N),Y(N),Z(N),(NDOF(I,N),I=1,3),KN
       END DO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Number unknowns
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       NEQ=0
       DO N=1,NUMNP
          DO I=1,3           
              IF (ID(I,N) .EQ. 0) THEN
                 NEQ=NEQ + 1
                 ID(I,N)=NEQ
                 IF(NDOF(I,N).EQ.0) THEN
                    WRITE(*,*) "ERROR, NDOF SHOULD .GE. 0!", N
                    STOP
                 ELSE
                    NEQ=NEQ+NDOF(I,N)-1
                 ENDIF
              ELSE 
                 ID(I,N)=0 
                 IF(NDOF(I,N).NE.0) THEN
                    WRITE(*,*) "ERROR, NDOF DOES NOT EQUAL TO 0!", N
                    STOP
                 ENDIF
              END IF             
          END DO
       END DO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Write equation numbers
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       WRITE (IOUT,2040) (N,(ID(I,N),I=1,3),N=1,NUMNP)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-------- Allocate memory for load vector 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ALLOCATE (R(NEQ),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (NRF(NEQ),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (RESS(3,NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (RNAT(3,NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO I=1,NEQ
          R(I)=0
          NRF(I)=1
       END DO     
!
       DO I=1,NUMNP
          DO J=1,3
             RNAT(J,I)=0.
             RESS(J,I)=0.
          ENDDO
       ENDDO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 2000  FORMAT(//,' N O D A L   P O I N T   D A T A',/)
 2010  FORMAT(' INPUT NODAL DATA',//)
 2015  FORMAT(//,' GENERATED NODAL DATA',//)
 2020  FORMAT('  NODE',10X,'BOUNDARY',25X,'NODAL POINT',35X,'NODE',/,&
       ' NUMBER     CONDITION  CODES',21X,'COORDINATES',32X,&
       'GENERATING',/,95X,'CODE',/,&
       15X,'X    Y    Z',15X,'X',12X,'Y',12X,'Z',11X,'NDOF(3,N)',8X,'KN')
 2030  FORMAT (I5,6X,3I5,6X,ES14.7,1X,ES14.7,1X,ES14.7,3X,3I5,6X,I4)
 2040  FORMAT(//,' EQUATION NUMBERS',//,'   NODE',9X,&
       'DEGREES OF FREEDOM',/,'  NUMBER',//,&
       '     N',13X,'X      Y      Z',/,(1X,I5,9X,I5,2X,I5,2X,I5))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!}}}
    END SUBROUTINE NODED
!
    SUBROUTINE LOADD(IIN,IOUT)
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
! *         NRF(NEQ)    : LOAD FUNCTION CODES                               
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: IIN,IOUT
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
!
       WRITE (IOUT,2015) NLOAD
       WRITE (IOUT,2020)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       N = 0
       I=0
       READ (IIN,*) N,IDIRNN,FLOADD,NFUNCC
       DO WHILE(N.NE.0)
          IF(N.GT.NUMNP) THEN
            WRITE(*,*) "ERROR,NODE NUMBER .GT. NUMNP IN LOAD_DATA!"
            STOP
          END IF
          I=I+1        
          IF(I.GT.NLOAD) THEN
            WRITE(*,*) "ERROR, NUMBER OF INPUTTED LOADS .GT. NLOAD!"
            STOP
          END IF
          IF(NFUNCC.EQ.0) NFUNCC=1      
          WRITE (IOUT,2030) N,IDIRNN,FLOADD,NFUNCC
!
          NumNEQ=ID(IDIRNN,N)
          RNAT(IDIRNN,N)=FLOADD
          IF (NumNEQ.NE.0) THEN
             R(NumNEQ)=FLOADD
             NRF(NumNEQ)=NFUNCC
             IF(NDOF(IDIRNN,N).GE.2) THEN
                DO J=1,NDOF(IDIRNN,N)-1
                  R(NumNEQ+J)=FLOADD
                  NRF(NumNEQ+J)=NFUNCC
                ENDDO
             ENDIF
          ENDIF
!
          READ (IIN,*) N,IDIRNN,FLOADD,NFUNCC
       END DO
       IF(I.NE.NLOAD) THEN
          WRITE(*,*) "ERROR, NUMBER OF INPUTTED LOADS .LT. NLOAD!"
          STOP
       END IF
!
       WRITE(IOUT,2040)
       DO I=1,NEQ
          WRITE(IOUT,2050) I,R(I),NRF(I)  
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
 2050  FORMAT(I5,5X,ES18.6,7X,I6)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!}}}
    END SUBROUTINE LOADD
!
    SUBROUTINE BOUNDD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
! *     SUBROUTINE BOUNDD 
! *      READ AND STORE DIRICHLET BOUNDARY DATA                            
! *      INPUT: 
! *     BOUND CASE
! *         IDIRN_BC(NBOUND): BOUNDARY DIRECTION CODES (1:X DIRECTION,2:Y DIRECTION; 3: Z DIRECTION)
! *         FBOUND(NBOUND): MAGNITUDE OF ESSENTIAL BOUNDARY
! *         NFUNC_BC(NBOUND): BOUNDARY FUNCTION CODES 
! *               1: CONSTANT
! *               2: LINEAR
! *               3: TRIGONOMETRIC                            
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: IIN,IOUT
       INTEGER :: I,J
!--------------------------------------------------
!      Internal variables
!--------------------------------------------------
       INTEGER :: N
       INTEGER :: IDIRNN_BC,NFUNCC_BC
       REAL(8) :: FBOUNDD
!--------------------------------------------------
!      Assistant variables
!--------------------------------------------------
       INTEGER :: NumNEQ=0    
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       WRITE (IOUT,2015) NBOUND
       WRITE (IOUT,2020)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       N = 0
       I=0
       READ (IIN,*) N,IDIRNN_BC,FBOUNDD,NFUNCC_BC
       DO WHILE(N.NE.0)
         IF(N.GT.NUMNP) THEN
           WRITE(*,*) "ERROR,NODE NUMBER .GT. NUMNP IN LOAD_DATA!"
           STOP
         END IF
         I=I+1        
         IF(I.GT.NBOUND) THEN
            WRITE(*,*) "ERROR, NUMBER OF INPUTTED LOADS .GT. NLOAD!"
            STOP
         END IF
         IF(NFUNCC_BC.EQ.0) NFUNCC_BC=1
!       
         WRITE (IOUT,2030) N,IDIRNN_BC,FBOUNDD,NFUNCC_BC
!
         NumNEQ=ID(IDIRNN_BC,N)         
         RESS(IDIRNN_BC,N)=FBOUNDD
         IF (NumNEQ.NE.0) THEN
            R(NumNEQ)=FBOUNDD
            NRF(NumNEQ)=NFUNCC_BC
            IF(NDOF(IDIRNN_BC,N).GE.2) THEN
               DO J=1,NDOF(IDIRNN_BC,N)-1
                  R(NumNEQ+J)=FBOUNDD
                  NRF(NumNEQ+J)=NFUNCC_BC
               ENDDO
            ENDIF
         END IF
!
         READ (IIN,*) N,IDIRNN_BC,FBOUNDD,NFUNCC_BC
       END DO
       IF(I.NE.NBOUND) THEN
          WRITE(*,*) "ERROR, NUMBER OF INPUTTED LOADS .LT. NLOAD!"
          STOP
       END IF
!
       WRITE(IOUT,2040)
       DO I=1,NEQ
          WRITE(IOUT,2050) I,R(I),NRF(I)  
       END DO
!
       RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 2015  FORMAT(/,'NUMBER OF BOUNDARY CONDITIONS =',I5,/)
 2020  FORMAT(' NODE        IDRIN_BC     BOUNDARY        NFUNC_BC',/,&
              'NUMBER        NUMBER      MAGNITUDE        NUMBER')
 2030  FORMAT(I5,9X,I5,3X,F13.3,5X,I6)
 2040  FORMAT(/,'EQUATION        LOAD           FUNCTION',/&
              ' NUMBER       MAGNITUDE          CODES ')
 2050  FORMAT(I5,5X,F13.3,7X,I6)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!}}}
    END SUBROUTINE BOUNDD
!
    SUBROUTINE INITLD(IIN,IOUT)
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
       INTEGER,INTENT(IN) :: IIN,IOUT
!
       INTEGER :: I,K,J,KK
       REAL(8) :: UXX,UYY,UZZ,VXX,VYY,VZZ
       REAL(8) :: DUX,DUY,DUZ,DVX,DVY,DVZ
       INTEGER :: KN,KNOLD,N,NOLD
       INTEGER :: NUMN,NUM
!--------------------------------------------------
!       Assistant variables
!--------------------------------------------------
       INTEGER :: NCOUNT
       INTEGER :: ALLocateStatus
!-------- Allocate memory
       ALLOCATE (UI(3,NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (VI(3,NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (U(NEQ),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (V(NEQ),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (ACC(NEQ),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"       
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO I=1,NUMNP
          DO J=1,3
             UI(J,I)=0
             VI(J,I)=0
          END DO
       END DO
!
       DO I=1,NEQ
          U(I)=0
          V(I)=0
          ACC(I)=0
       END DO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       WRITE (IOUT,2000)
       WRITE (IOUT,2010)
       WRITE (IOUT,2020)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
       N = 0
       NCOUNT=0                  !RECORD NUMBER OF INPUTTED NODES     
       KNOLD=0
       NOLD=0

       READ (IIN,*) N,UXX,UYY,UZZ,VXX,VYY,VZZ,KN  
       DO WHILE (N.NE.0)
          IF(N.GT.NUMNP) THEN
             WRITE(*,*) "ERROR, NODE NUMBER .GT. NUMNP!"
             STOP
          END IF
          UI(1,N)=UXX
          UI(2,N)=UYY
          UI(3,N)=UZZ
          VI(1,N)=VXX
          VI(2,N)=VYY
          VI(3,N)=VZZ  
          NCOUNT=NCOUNT+1
!----------------write out initial conditions
          WRITE (IOUT,2030) N,UI(1,N),UI(2,N),UI(3,N),VI(1,N),VI(2,N),VI(3,N)
           IF (KNOLD.NE.0) THEN
              NUM=(N-NOLD)/KNOLD
              NUMN=NUM - 1
              IF (NUMN.GE.1) THEN
                 DUX=(UI(1,N)-UI(1,NOLD))/NUM
                 DUY=(UI(2,N)-UI(2,NOLD))/NUM
                 DUZ=(UI(3,N)-UI(3,NOLD))/NUM
                 DVX=(VI(1,N)-VI(1,NOLD))/NUM
                 DVY=(VI(2,N)-VI(2,NOLD))/NUM
                 DVZ=(VI(3,N)-VI(3,NOLD))/NUM
                 K=NOLD
                 DO J=1,NUMN
                    KK=K
                    K=K + KNOLD
                    UI(1,K) = UI(1,KK) + DUX
                    UI(2,K) = UI(2,KK) + DUY
                    UI(3,K) = UI(3,KK) + DUZ                
                    VI(1,K) = VI(1,KK) + DVX
                    VI(2,K) = VI(2,KK) + DVY
                    VI(3,K) = VI(3,KK) + DVZ
                    NCOUNT=NCOUNT+1
                 END DO
              END IF
            ENDIF
!        
          NOLD=N
          KNOLD=KN
          READ (IIN,*) N,UXX,UYY,UZZ,VXX,VYY,VZZ,KN   
       END DO
       IF(NCOUNT.NE.NUMNP) THEN
          WRITE(*,*) "ERROR,NUMBER OF INPUTTED INITIAL DATA .LT. NUMNP!"
          STOP
       END IF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Write complete initial data
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       WRITE (IOUT,2015)
       WRITE (IOUT,2020)
       DO N=1,NUMNP
         WRITE (IOUT,2030) N,UI(1,N),UI(2,N),UI(3,N),VI(1,N),VI(2,N),VI(3,N)
       END DO       
!----------------data transform
       DO N=1,NUMNP
          DO I=1,3
            IF(ID(I,N).NE.0) THEN
              K=ID(I,N)
              U(K)=UI(I,N)
              V(K)=VI(I,N)
            END IF
          END DO
       END DO
!----------------Write new Initial conditions
       WRITE(IOUT,2040)
       DO I=1,NEQ
          WRITE(IOUT,2050) I,U(I),V(I)  
       END DO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 2000  FORMAT(//,' N O D A L   I N I T I A L   D A T A',/)
 2010  FORMAT(' INPUT NODAL DATA',//)
 2015  FORMAT(//,' GENERATED INITIAL DATA',//)
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
    SUBROUTINE MATERIALD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
! *     SUBROUTINE MATERIALD
! *      READ AND STORE MATERIAL DATA                           
! *      INPUT: 
! *        MATP(7,NMAT)  :  MATERIAL PROPERTIES
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: IIN,IOUT
! Internal variables        
        INTEGER :: N=0
        REAL(8) :: E=0,PR=0,DEN=0,PAR1=0
! Assistant variables
        INTEGER :: I,J
        INTEGER :: NCOUNT
        INTEGER :: AllocateStatus
!-------- Allocate memory
        ALLOCATE (MATP(7,NMAT),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!
        DO I=1,NMAT
           DO J=1,7
              MATP(J,I)=0
           END DO
        END DO    
! SET DEFAULT VALUE
        NCOUNT=0
!
        WRITE (IOUT,2000)
        WRITE (IOUT,2010)
        WRITE (IOUT,2020)
!
        READ(IIN,*) N,E,PR,DEN,PAR1
        DO WHILE(N.NE.0)
           IF(N.GT.NMAT) THEN
                WRITE(*,*) "ERROR, MATERIAL SET NUMBER .GT. NMAT"
                STOP
           END IF
           NCOUNT=NCOUNT+1
           MATP(1,N)=E
           MATP(2,N)=PR
           MATP(3,N)=DEN
           MATP(4,N)=PAR1
           WRITE(IOUT,2030) N, E, PR,DEN,PAR1
           READ(IIN,*) N,E,PR,DEN,PAR1
        END DO
        IF(NCOUNT.NE.NMAT) THEN
            WRITE(*,*) "ERROR, NUMBER OF MATERIAL SET .LT. NMAT"
            STOP
        END IF
        RETURN
!
 2000  FORMAT(//,' M A T E R I A L    D A T A',/)
 2010  FORMAT(' INPUT MATERIAL DATA',//)
 2020  FORMAT(' SET ',22X,' E ',14X,' PR ',14X,'DEN',14X,'PAR1',/,&
              'NUMBER')
 2030  FORMAT(I6,10X,ES13.3,5X,F13.3,5X,F13.3,5X,F13.3)
!}}}
    END SUBROUTINE MATERIALD
!
    SUBROUTINE BOUNDARY_IN(IIN,IOUT,ITEM3)
!{{{
!--------------------------------------------------
!   SELECT BOUNDARY TYPE
!   NPAR()
!   LINE
!       NPAR(1) = 1
!       NPAR(2) = NUMBC
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT,ITEM3
!   ASSISTANT VARIABLES
    INTEGER :: I,J,L,M
    INTEGER,PARAMETER :: KwBC=1
    CHARACTER(20),DIMENSION(KwBC) :: KwBoundary 
!
!---  ASSISTANT VARIABLES    
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KwBoundary(1)='\LINE'
!
    WRITE(IOUT,2000)
!
    IF(NBGROUP<1e-7) THEN
       WRITE(*,*) "ERROR, NBGROUP EQUALS TO ZEROS"
       STOP
    ENDIF
!
    DO J=1,NBGROUP
       I = KeyWord(KwBoundary,KwBC)
       WRITE(IOUT,1000) KwBoundary(I)
       SELECT CASE(I)
       CASE(1)
          NPAR=0
          NPAR(1)=1
          CALL LINED(IIN,IOUT)
          WRITE(ITEM3) (NPAR(L),L=1,7)
          WRITE(ITEM3) ((BCED(L,M),L=1,2),M=1,NPAR(2)),((BCIND(L,M),L=1,4),M=1,NPAR(2)),((BCPAR(M,L),M=1,2),L=1,NPAR(2))
          ! Deallocate variables
          DEALLOCATE (BCED)
          DEALLOCATE (BCIND)
          DEALLOCATE (BCPAR) 
       CASE DEFAULT
          WRITE(*,*) "ERROR, INPUT ERROR IN ELEMENTD"
          STOP
       END SELECT 
    END DO
!
    RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 1000 FORMAT(/,A20)
 2000 FORMAT(//,' B O U N D A R Y   I N F O R M A T I O N',///)
!}}}
    END SUBROUTINE BOUNDARY_IN
!
    SUBROUTINE LINED(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE LINED
! *      READ AND STORE LINE DATA
! *      LINE
! *         NPAR(1) = 1
! *         NPAR(2) = NUMBC
! *      INPUT:
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: NumBc
    INTEGER :: N,II,JJ,KN,INDI,INDJ,FUNCI,FUNCJ
    INTEGER :: OINDI,OINDJ,OFUNCI,OFUNCJ
    REAL(8) :: FLOADI,FLOADJ
    REAL(8) :: OFLOADI,OFLOADJ
    INTEGER :: KNOLD,NOLD,NUM,NUMN
    INTEGER :: IOLD,JOLD,NII,NJJ
    INTEGER :: NCOUNT
    INTEGER :: DII
!
    INTEGER :: AllocateStatus
!
    CHARACTER(20),DIMENSION(2) :: KWLINE
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWLINE(1) = 'NUMBC:'
    KWLINE(2) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
!
    I = KeyWord(KWLINE,2)
    DO WHILE(I<=2)
       SELECT CASE(I)
       CASE(1)
           J = GetInt()
           NPAR(2)=J
       CASE(2)
           EXIT
       CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN LINED"
           STOP
       END SELECT
       I = KeyWord(KWLINE,2)
    END DO
!
    WRITE(IOUT,1000) NPAR(2)
!   
    IF(NPAR(2).EQ.0) THEN
            WRITE(*,*) "ERROR, NUMBER OF LINE BOUNDARY .EQ.0"
            STOP
    END IF
    NumBc=NPAR(2)
!
!-------- Allocate memory
!
        ALLOCATE (BCED(2,NumBc),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
        ALLOCATE (BCIND(4,NumBc),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
        ALLOCATE (BCPAR(2,NumBc),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"        
! ----------------------------------------------------------------
        DO I=1,NumBc
           DO J=1,4
              BCIND(J,I)=0
           END DO
           BCED(1,I)=0
           BCED(2,I)=0
           BCPAR(1,I)=0.
           BCPAR(2,I)=0.
        END DO
!
        WRITE(IOUT,2000)
        NCOUNT=0
        READ(IIN,*) N, II, JJ, INDI, INDJ, FLOADI, FLOADJ, FUNCI, FUNCJ, KN
        DO WHILE(N.NE.0)
           IF(N.GT.NumBc) THEN
              WRITE(*,*) "ERROR N .GT. NumEle"
              STOP
           END IF
           NCOUNT=NCOUNT+1
!
           BCED(1,N)=II 
           BCED(2,N)=JJ
!          
           BCPAR(1,N)=FLOADI
           BCPAR(2,N)=FLOADJ
!
           BCIND(1,N)=INDI
           BCIND(2,N)=INDJ
           BCIND(3,N)=FUNCI
           BCIND(4,N)=FUNCJ
!
           IF (KNOLD.NE.0) THEN
              NUM=(N-NOLD)/KNOLD
              NUMN=NUM - 1
              IF (NUMN.GE.1) THEN
                  K=NOLD
                  NII=IOLD
                  NJJ=JOLD
                  DII=(II-IOLD)/NUM
                 DO J=1,NUMN
                    K=K + KNOLD
!
                    NII=NII+DII
                    NJJ=NJJ+DII
!
                    BCED(1,K)=NII 
                    BCED(2,K)=NJJ
!
                    BCPAR(1,K)=OFLOADI
                    BCPAR(2,K)=OFLOADJ
!
                    BCIND(1,K)=OINDI
                    BCIND(2,K)=OINDJ
                    BCIND(3,K)=OFUNCI
                    BCIND(4,K)=OFUNCJ
!
                   WRITE(IOUT,2050) K,NII,NJJ,OINDI,OINDJ,OFLOADI,OFLOADJ,OFUNCI,OFUNCJ                 
                   NCOUNT=NCOUNT+1
                 END DO
              END IF
            ENDIF
            WRITE(IOUT,2050) N, II, JJ, INDI, INDJ, FLOADI, FLOADJ, FUNCI, FUNCJ
!        
            NOLD=N
            KNOLD=KN
            IOLD=II
            JOLD=JJ
!            
            OINDI=INDI
            OINDJ=INDJ
            OFLOADI=FLOADI
            OFLOADJ=FLOADJ
            OFUNCI=FUNCI
            OFUNCJ=FUNCJ
!
           READ(IIN,*)  N, II, JJ, INDI, INDJ, FLOADI, FLOADJ, FUNCI, FUNCJ, KN
        END DO
!
        IF(NCOUNT.NE.NumBc) THEN
              WRITE(*,*) "ERROR, BOUNDARY NUMBER .LT. NUMBC",NCOUNT
              STOP
        END IF
!
    RETURN
!--------------------------------------------------
 1000 Format(/,&
              '    THE BOUNDARY TYPE: LINE ',/,&
              '    THE NUMBER OF LINE BOUNDARY      =',I5,/)
 2000 FORMAT(/,&
              ' BOUNDARY   ELEMENT  ELEMENT  IND(2,N)           FLOAD(2,N)           FUNC(2,N)',/,&
              ' NUMBER-N      I        J       ',/)
 2050 FORMAT(I5,6X,I5,4X,I5,2X,2I5,2X,2ES13.3,2X,2I5)
!}}}
    END SUBROUTINE LINED
!
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
!   QUADS
!       NPAR(1) = 2
!       NPAR(2) = ELEMENT TYPE
!                 0: AXISYMMETRIC
!                 1: PLANE STRAIN
!                 2: PLANE STRESS
!       NPAR(3) = MTYPE
!       NPAR(4) = NUMELE
!       NPAR(5) = NINTG
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT,ITEM1
!
    INTEGER :: I,J,L,M
    INTEGER,PARAMETER :: KwEle=4
    CHARACTER(20),DIMENSION(KwEle) :: KwElement
!---  ASSISTANT VARIABLES
    INTEGER :: ALLocateStatus
!---- ALLOCATE VARIABLES
    ALLOCATE (MHT(NEQ),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    DO I=1,NEQ
       MHT(I)=0
    END DO       
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KwElement(1)='\TRUSS'
    KwElement(2)='\QUADS' 
    KwElement(3)='\TRUSS_MFS'
    KwElement(4)='\QUADS_MFS'
!
    WRITE(IOUT,2000)
!
    IF(NEGROUP<1e-7) THEN
       WRITE(*,*) "ERROR, NEGROUP EQUALS TO ZEROS"
       STOP
    ENDIF
!
    DO J=1,NEGROUP
       I = KeyWord(KwElement,KwEle)
       WRITE(IOUT,1000) KwElement(I)
       SELECT CASE(I)
       CASE(1)
          NPAR=0
          NPAR(1)=1
          CALL TRUSSD(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,7)
          WRITE(ITEM1) ((XYZ(L,M),L=1,6),M=1,NPAR(3)),((LM(L,M),L=1,6),M=1,NPAR(3)),(ELEPAR(1,L),L=1,NPAR(3))
          ! Deallocate variables
          DEALLOCATE (XYZ)
          DEALLOCATE (LM)
          DEALLOCATE (ELEPAR)
       CASE(2)
          NPAR=0
          NPAR(1)=2
          CALL QUADSD(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,7)
          WRITE(ITEM1) ((XYZ(L,M),L=1,8),M=1,NPAR(4)),((LM(L,M),L=1,8),M=1,NPAR(4)),(ELEPAR(1,L),L=1,NPAR(4))
          ! Deallocate variables
          DEALLOCATE (XYZ)
          DEALLOCATE (LM)
          DEALLOCATE (ELEPAR)
        CASE(3)
          NPAR=0
          NPAR(1)=3
          CALL TRUSS_MFSD(IIN,IOUT)         
          WRITE(ITEM1) (NPAR(L),L=1,7)
          WRITE(ITEM1) ((NDMFS(M,L),M=1,3),L=1,NPAR(3)),((ELEPAR(M,L),M=1,2),L=1,NPAR(3)),((LM(M,L),M=1,NPAR(4)),L=1,NPAR(3))
          ! Deallocate variables
          DEALLOCATE(NDMFS)
          DEALLOCATE(ELEPAR)
          DEALLOCATE(LM)
       CASE(4)
          NPAR=0
          NPAR(1)=4
          CALL QUADS_MFSD(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,7)
          WRITE(ITEM1) ((NDMFS(M,L),M=1,8),L=1,NPAR(4)),((ELEPAR(M,L),M=1,6),L=1,NPAR(4)),((LM(M,L),M=1,NPAR(5)),L=1,NPAR(4))
          DEALLOCATE(NDMFS)
          DEALLOCATE(ELEPAR)
          DEALLOCATE(LM)
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
    INTEGER :: I,J,K
    INTEGER :: NumEle
    INTEGER :: N,II,JJ,KN
    INTEGER :: KNOLD,NOLD,NUM,NUMN
    INTEGER :: IOLD,JOLD,NII,NJJ
    REAL(8) :: AREA,AOLD
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
       I = KeyWord(KWTRUSS,3)
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
        ALLOCATE (ELEPAR(1,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"        
! ----------------------------------------------------------------
        DO I=1,NumEle
           DO J=1,6
              XYZ(J,I)=0
              LM(J,I)=0
           END DO
           ELEPAR(1,I)=0.
        END DO
! ----------------------------------------------------------------
        WRITE(IOUT,2000)
        NCOUNT=0
        READ(IIN,*) N, II, JJ, AREA, KN
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
           ELEPAR(1,N)=AREA
           DO I=1,3
             LM(I,N)=ID(I,II)
             LM(I+3,N)=ID(I,JJ)
           END DO
           CALL COLHT(6,LM(:,N))
!
           IF (KNOLD.NE.0) THEN
              NUM=(N-NOLD)/KNOLD
              NUMN=NUM - 1
              IF (NUMN.GE.1) THEN
                  K=NOLD
                  NII=IOLD
                  NJJ=JOLD
                 DO J=1,NUMN
                    K=K + KNOLD
!
                    NII=NII+KNOLD
                    NJJ=NJJ+KNOLD
!
                    XYZ(1,K)=X(NII)
                    XYZ(2,K)=Y(NII)
                    XYZ(3,K)=Z(NII)
!
                    XYZ(4,K)=X(NJJ)
                    XYZ(5,K)=Y(NJJ)
                    XYZ(6,K)=Z(NJJ)
!
                    ELEPAR(1,K)=AOLD
!
                   DO I=1,3
                     LM(I,K)=ID(I,NII)
                     LM(I+3,K)=ID(I,NJJ)
                   END DO
                   CALL COLHT(6,LM(:,K))
                   WRITE(IOUT,2050) K,NII,NJJ,AOLD                   
                   NCOUNT=NCOUNT+1
                 END DO
              END IF
            ENDIF
            WRITE(IOUT,2050) N,II,JJ,AREA
!        
            NOLD=N
            KNOLD=KN
            IOLD=II
            JOLD=JJ
            AOLD=AREA
!
           READ(IIN,*) N,II,JJ,AREA,KN
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
              ' ELEMENT     NODE     NODE     AREA   ',/,&
              ' NUMBER-N      I        J       ',/)
 2050 FORMAT(I5,6X,I5,4X,I5,2X,ES13.3)
!}}}
    END SUBROUTINE TRUSSD
!
    SUBROUTINE QUADSD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE QUADSD
! *      READ AND STORE TRUSS DATA
! *  QUADS
! *      NPAR(1) = 2
! *      NPAR(2) = ELEMENT TYPE
! *                0: AXISYMMETRIC
! *                1: PLANE STRAIN
! *                2: PLANE STRESS
! *      NPAR(3) = MTYPE
! *      NPAR(4) = NUMELE
! *      NPAR(5) = NINTG
! *      INPUT:
! *        MTYPE    -  MATERIAL SET NUMBER
! *    	   NUMELE   -  THE NUMBER OF TRUSS ELEMENTS
! *        NINTG    -  NUMBER OF NUMERICAL INTEGRATION POINTS
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J
    INTEGER :: NumEle
    REAL(8) :: THIC
    INTEGER :: N,II,JJ,LL,MM
    INTEGER :: NCOUNT
    INTEGER :: AllocateStatus
!
    CHARACTER(20),DIMENSION(5) :: KWQUADS
!- - - - - - - - - - - - - - -  
!-     VARIABLE DEFINITION   -
!- - - - - - - - - - - - - - -
    KWQUADS(1) = 'ITYPE:'
    KWQUADS(2) = 'MTYPE:'
    KWQUADS(3) = 'NUMELE:'
    KWQUADS(4) = 'NINTG:'
    KWQUADS(5) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=0
    NPAR(3)=1
    NPAR(4)=0
    NPAR(5)=2
!  
    I = KeyWord(KWQUADS,5)
    DO WHILE(I<=5)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(3)=J
         CASE(3)
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           J = GetInt()
           NPAR(5)=J
         CASE(5)
           EXIT
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN TRUSSD"
           STOP
       END SELECT
       I = KeyWord(KWQUADS,5)
    END DO
!
    WRITE(IOUT,1000) NPAR(2),NPAR(3),NPAR(4),NPAR(5)
!
    IF(NPAR(4).EQ.0) THEN
            WRITE(*,*) "ERROR, NUMBER OF TRUSS ELEMENT .EQ.0"
            STOP
    END IF
    NumEle=NPAR(4)
!
!-------- Allocate memory
!
        ALLOCATE (XYZ(8,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
        ALLOCATE (LM(8,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
        ALLOCATE (ELEPAR(1,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"        
! ----------------------------------------------------------------
        DO I=1,NumEle
           DO J=1,8
              XYZ(J,I)=0
              LM(J,I)=0
           END DO
           ELEPAR(1,I)=0.
        END DO
! ----------------------------------------------------------------
        WRITE(IOUT,2000)
        NCOUNT=0
        READ(IIN,*) N, II, JJ, LL, MM, THIC
        DO WHILE(N.NE.0)
           IF(N.GT.NumEle) THEN
              WRITE(*,*) "ERROR N .GT. NumEle"
              STOP
           END IF
           NCOUNT=NCOUNT+1
!
           XYZ(1,N)=X(II)
           XYZ(2,N)=Y(II)
!
           XYZ(3,N)=X(JJ)
           XYZ(4,N)=Y(JJ)
!
           XYZ(5,N)=X(LL)
           XYZ(6,N)=Y(LL)
!
           XYZ(7,N)=X(MM)
           XYZ(8,N)=Y(MM)
!
           ELEPAR(1,N)=THIC
           DO I=1,2
             LM(I  ,N)=ID(I,II)
             LM(I+2,N)=ID(I,JJ)
             LM(I+4,N)=ID(I,LL)
             LM(I+6,N)=ID(I,MM)
           END DO
           CALL COLHT(8,LM(:,N))
           WRITE(IOUT,2050) N,II,JJ,LL,MM,THIC
!
           READ(IIN,*) N,II,JJ,LL,MM,THIC
        END DO
!
        IF(NCOUNT.NE.NumEle) THEN
              WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumEle"
              STOP
        END IF
    RETURN
!--------------------------------------------------
 1000 Format(/,&
              '    THE ELELMENT TYPE: QUADS          =',I5,/,&
              '               0: AXISYMMETRIC            ',/,&
              '               1: PLANE STRAIN            ',/,&
              '               2: PLANE STRESS            ',/,&                  
              '    THE MATERIAL SET NUMBER           =',I5,/,&
              '    THE NUMBER OF QUADS ELEMENTS      =',I5,/,&
              '    GAUSS NUMERICAL INTEGRATION ORDER =',I5,/)
 2000 FORMAT(/,&
              ' ELEMENT     NODE     NODE     NODE     NODE     THIC   ',/,&
              ' NUMBER-N      I        J        L        M       ',/)
 2050 FORMAT(I5,6X,I5,4X,I5,4X,I5,4X,I5,2X,ES10.3)
!}}}
    END SUBROUTINE QUADSD
!
    SUBROUTINE TRUSS_MFSD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE TRUSS_MFSD
! *      READ AND STORE TRUSS DATA
! *      TRUSS_MFS
! *         NPAR(1) = 3
! *         NPAR(2) = MTYPE
! *         NPAR(3) = NUMELE
! *         NPAR(4) = ELEASS
! *      INPUT:
! *        MTYPE    -  MATERIAL SET NUMBER
! *    	   NUMELE   -  THE NUMBER OF TRUSS ELEMENTS  
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: NumEle
    INTEGER :: N,M,II,KN,IND
    INTEGER :: KNOLD,NOLD,NUM,NUMN,IND_OLD
    REAL(8) :: DIRN,DIRN_OLD
    INTEGER :: IOLD,NII
    REAL(8) :: RAD,ROLD
    INTEGER :: NCOUNT
    INTEGER :: NEleAss
    INTEGER :: INDSP,INDSP_OLD
    INTEGER :: AllocateStatus
!
    CHARACTER(20),DIMENSION(4) :: KWTRUSS_MFS
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWTRUSS_MFS(1) = 'MTYPE:'
    KWTRUSS_MFS(2) = 'NUMELE:'
    KWTRUSS_MFS(3) = 'ELEASS:'
    KWTRUSS_MFS(4) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=0
    NPAR(4)=NUMNP
!   
    I = KeyWord(KWTRUSS_MFS,4)
    DO WHILE(I<=5)
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
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           EXIT
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN TRUSS_MFSD"
           STOP
       END SELECT
       I = KeyWord(KWTRUSS_MFS,4)
    END DO
!
    WRITE(IOUT,1000) NPAR(2),NPAR(3),NPAR(4)
!
    IF(NPAR(3).EQ.0) THEN
            WRITE(*,*) "ERROR, NUMBER OF TRUSS ELEMENT .EQ.0"
            STOP
    END IF
!
    NumEle  = NPAR(3)
    NEleAss = NPAR(4)
!
!-------- Allocate memory
!
        ALLOCATE (ELEPAR(2,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"        
        ALLOCATE (NDMFS(3,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
        ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
! ----------------------------------------------------------------
        DO I=1,NumEle
           DO J=1,3
              NDMFS(J,I)=0
           ENDDO
           ELEPAR(1,I)=0.
           ELEPAR(2,I)=0.
!
           DO J=1,NEleAss
              LM(J,I)=0
           ENDDO
        END DO
! ----------------------------------------------------------------
        WRITE(IOUT,2000)
        NCOUNT=0
        READ(IIN,*) N, II, IND,RAD,DIRN,INDSP,KN
        DO WHILE(N.NE.0)
!           IF(II.GT.NumEle) THEN
!              WRITE(*,*) "ERROR N .GT. NumEle"
!              STOP
!           END IF
           NCOUNT=NCOUNT+1
!          
           ELEPAR(2,N)=RAD
           ELEPAR(1,N)=DIRN
           NDMFS(1,N)=II
           NDMFS(2,N)=IND
           NDMFS(3,N)=INDSP
!
           IF (KNOLD.NE.0) THEN
              NUM=(N-NOLD)/KNOLD
              NUMN=NUM - 1
              IF (NUMN.GE.1) THEN
                  K=NOLD
                  NII=IOLD
                 DO J=1,NUMN
                    K=K + KNOLD
!
                    NII=NII+KNOLD
                    NDMFS(1,K)=NII
                    NDMFS(2,K)=IND_OLD
                    NDMFS(3,K)=INDSP_OLD                 
!
                    ELEPAR(2,K)=ROLD
                    ELEPAR(1,K)=DIRN_OLD
!
                   WRITE(IOUT,2050) K,NII,IND_OLD,ROLD,DIRN_OLD,INDSP_OLD              
                   NCOUNT=NCOUNT+1
                 END DO
              END IF
            ENDIF
            WRITE(IOUT,2050) N,II,IND,RAD,DIRN,INDSP
!        
            NOLD=N
            KNOLD=KN
            IOLD=II
            IND_OLD=IND
            ROLD=RAD
            DIRN_OLD=DIRN
            INDSP_OLD=INDSP
!
            READ(IIN,*) N, II,IND,RAD,DIRN,INDSP,KN
        END DO
!
        IF(NCOUNT.NE.NumEle) THEN
              WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumEle"
              STOP
        END IF
        CALL MFS_TRUSS_ASSD(3,2,NumEle,NEleAss,NDMFS(1,1),ELEPAR(1,1),LM)
!       CONSTRUCT MHT
        DO N=1,NumEle
           DO I=1,NEleAss
             IF(LM(I,N).NE.0) THEN
              M=LM(I,N)
              CALL COLHT_MFS(NDMFS(1,N),NDMFS(1,M))
             ENDIF
           ENDDO
        ENDDO
!
    RETURN
!--------------------------------------------------
 1000 Format(/,&
              '    THE ELELMENT TYPE: TRUSS',/,&
              '    THE MATERIAL SET NUMBER                   =',I5,/,&
              '    THE NUMBER OF TRUSS ELEMENTS              =',I5,/,&
              '    THE MAXIMUM ELEMENT ASSOCIATION NUMBER    =',I5)
 2000 FORMAT(/,&
              ' ELEMENT     NODE        IND       RAD        DIRN     INDSP',/,&
              ' NUMBER-N   NUMBER               ',/)
 2050 FORMAT(I5,6X,I5,6X,I5,2X,ES13.3,2X,ES10.3,2X,I5)
!}}}
    END SUBROUTINE TRUSS_MFSD
!
    SUBROUTINE QUADS_MFSD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE QUADS_MFSD
! *      READ AND STORE TRUSS DATA
! *      TRUSS_MFS
! *         NPAR(1) = 4
! *      NPAR(2) = ELEMENT TYPE
! *                1: PLANE STRAIN
! *                2: PLANE STRESS
! *      NPAR(3) = MTYPE
! *      NPAR(4) = NUMELE
! *      NPAR(5) = ELEASS
! *      INPUT:
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J
    INTEGER :: NumEle,NEleAss
    INTEGER :: NCOUNT 
    INTEGER :: N,II,KN
    REAL(8) :: R_IN1,R_IN2,R_OUT1,R_OUT2
    INTEGER :: NOLD,IOLD,KNOLD
    INTEGER :: INDSP,INDSP_OLD
    INTEGER :: NUM,NUMN,K,NII
    INTEGER :: IND1,IND2,IND1O,IND2O
    REAL(8) :: R_IN_OLD1,R_IN_OLD2,R_OUT_OLD1,R_OUT_OLD2
    INTEGER :: ALLocateStatus 
    INTEGER :: M
    INTEGER :: DII
!
    CHARACTER(20),DIMENSION(5) :: KWQUADS_MFS
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWQUADS_MFS(1) = 'ITYPE:'
    KWQUADS_MFS(2) = 'MTYPE:'
    KWQUADS_MFS(3) = 'NUMELE:'
    KWQUADS_MFS(4) = 'ELEASS:'
    KWQUADS_MFS(5) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=1
    NPAR(4)=1
    NPAR(5)=NUMNP
!
    I = KeyWord(KWQUADS_MFS,6)
    DO WHILE(I<=5)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           IF (J .GT. 3) THEN
               WRITE(*,*) "ERROR ITYPE .GT. 3"
               STOP
           ENDIF
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(3)=J
         CASE(3)
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           J = GetInt()
           NPAR(5)=J
         CASE(5)
           EXIT
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN QUADS_MFSD"
           STOP
       END SELECT
       I = KeyWord(KWQUADS_MFS,5)
    END DO
!
    WRITE(IOUT,1000) NPAR(2),NPAR(3),NPAR(4),NITP_IN,NPAR(5)
!
    IF(NPAR(4).EQ.0) THEN
            WRITE(*,*) "ERROR, NUMBER OF QUADS ELEMENT .EQ.0"
            STOP
    END IF
!
    NumEle  = NPAR(4)
    NEleAss = NPAR(5)
!
!-------- Allocate memory
!
        ALLOCATE (NDMFS(8,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
        ALLOCATE (ELEPAR(6,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"        
        ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!
        DO I=1,NumEle
           DO J=1,8
              NDMFS(J,I)=0
           ENDDO
!
           DO J=1,6
              ELEPAR(J,I)=0.
           ENDDO
!
           DO J=1,NEleAss
              LM(J,I)=0
           ENDDO
        END DO
! ----------------------------------------------------------------
        WRITE(IOUT,2000)
        NCOUNT=0
        READ(IIN,*) N,II,IND1,IND2,R_IN1,R_IN2,R_OUT1,R_OUT2,INDSP,KN
        DO WHILE(N.NE.0)
           NCOUNT=NCOUNT+1
           IF(N.GT.NumEle) THEN
              WRITE(*,*) "ERROR, ELEMENT NUMBER .GT. NumEle"
              STOP
           END IF
!          
           ELEPAR(1,N)=R_IN1
           ELEPAR(2,N)=R_IN2
           ELEPAR(3,N)=R_OUT1
           ELEPAR(4,N)=R_OUT2
           NDMFS(1,N)=II
           NDMFS(2,N)=IND1
           NDMFS(3,N)=IND2
           NDMFS(6,N)=INDSP
!
           IF (KNOLD.NE.0) THEN
              NUM=(N-NOLD)/KNOLD
              NUMN=NUM - 1
              IF (NUMN.GE.1) THEN
                  K=NOLD
                  NII=IOLD
                  DII=(II-IOLD)/NUM
                 DO J=1,NUMN
                    K=K + KNOLD
!
                    NII=NII+DII
                    NDMFS(1,K)=NII
                    NDMFS(2,K)=IND1O
                    NDMFS(3,K)=IND2O 
                    NDMFS(6,K)=INDSP_OLD                   
!
                    ELEPAR(1,K)=R_IN_OLD1
                    ELEPAR(2,K)=R_IN_OLD2
                    ELEPAR(3,K)=R_OUT_OLD1
                    ELEPAR(4,K)=R_OUT_OLD2
!
                   WRITE(IOUT,2050) K,NII,IND1O,IND2O,R_IN_OLD1,R_IN_OLD2,R_OUT_OLD1,R_OUT_OLD2,INDSP_OLD             
                   NCOUNT=NCOUNT+1
                 END DO
              END IF
           ENDIF
           WRITE(IOUT,2050) N,II,IND1,IND2,R_IN1,R_IN2,R_OUT1,R_OUT2,INDSP
!        
           NOLD=N
           KNOLD=KN
           IOLD=II
           R_IN_OLD1=R_IN1
           R_IN_OLD2=R_IN2
           R_OUT_OLD1=R_OUT1
           R_OUT_OLD2=R_OUT2
           IND1O=IND1
           IND2O=IND2
           INDSP_OLD=INDSP
!
           READ(IIN,*) N,II,IND1,IND2,R_IN1,R_IN2,R_OUT1,R_OUT2,INDSP,KN
        END DO
!
        IF(NCOUNT.NE.NumEle) THEN
              WRITE(*,*) "ERROR, ELEMENT NUMBER .LT. NumEle", NCOUNT
              STOP
        END IF
!
        CALL MFS_2D_ASSD(8,6,NumEle,NEleAss,NDMFS(1,1),ELEPAR(1,1),LM(1,1))
!       CONSTRUCT MHT
        DO N=1,NumEle
           DO I=1,NEleAss
             M=LM(I,N)
             IF((M.NE.0).AND.(M.GE.N)) THEN
              CALL COLHT_MFS(NDMFS(1,N),NDMFS(1,M))
             ENDIF
           ENDDO
        ENDDO
!
        RETURN
!
 1000 Format(/,&
              '    THE ELELMENT TYPE: QUADS          =',I5,/,&
              '               1: PLANE STRAIN            ',/,&
              '               2: PLANE STRESS            ',/,&                  
              '    THE MATERIAL SET NUMBER           =',I5,/,&
              '    THE NUMBER OF QUADS ELEMENTS      =',I5,/,&
              '    THE NUMBER OF INTEGRATION POINTS  =',I5,/,&
              '    THE MAXIMUM ELEMENT ASSOCIATION NUMBER    =',I5)
 2000 FORMAT(/,&
              ' ELEMENT     NODE      IND(2,N)                            RAD(4,N)                       INDSP ',/,&
              ' NUMBER-N   NUMBER                 ',/)
 2050 FORMAT(I5,6X,I5,6X,2I5,2X,4ES13.3,2X,I5)
!}}}
    END SUBROUTINE QUADS_MFSD
!
    SUBROUTINE MFS_TRUSS_ASSD(NEle1,NEle2,NumEle,NEleAss,MFS_NODE,MFS_PAR,MFS_ASS)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE FOR CALCULATING NODE ASSOCIATION LM      
! *      NumEle   :  Number of elements
! *      EleAss   :  Maximum element association number
! *      MFS_NODE :  
! *      MFS_PAR  :
! *      MFS_ASS  : 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
!
    INTEGER,INTENT(IN) :: NumEle, NEleAss, NEle1,NEle2
    INTEGER,DIMENSION(NEle1,NumEle)   :: MFS_NODE
    REAL(8),DIMENSION(NEle2,NumEle)   :: MFS_PAR
    INTEGER,DIMENSION(NEleAss,NumEle) :: MFS_ASS
!
!   ASSISTANT VARIABLES
!
    INTEGER :: I,J
    INTEGER,DIMENSION(NumEle) :: K
    INTEGER :: N,M
    INTEGER :: KK
    REAL(8) :: RAD1,RAD2
    REAL(8) :: DIST  
    LOGICAL :: EXPRN  
!
    KK=NEleAss+1
!
    DO I=1,NumEle
       K(I)=1
    ENDDO
    RAD1=0
    RAD2=0
    DO I=1,NumEle
       N=MFS_NODE(1,I)
       RAD1=MFS_PAR(NEle2,I) 
       DO J=I,NumEle              
          M=MFS_NODE(1,J)
          RAD2=MFS_PAR(NEle2,J)
          DIST=(RAD1+RAD2)
!
          DIST=ABS(X(N)-X(M))-DIST
!
          EXPRN=DIST>=(-1E-8) 
          IF(EXPRN) CYCLE 
!          
          MFS_ASS(K(I),I)=J
          K(I)=K(I)+1
          IF(I.NE.J) THEN
             MFS_ASS(K(J),J)=I
             K(J)=K(J)+1
          ENDIF
          IF(K(I).GT.KK) THEN
             WRITE(*,*) "ERROR, NEleAss too small!"
             STOP
          ENDIF
       ENDDO
    ENDDO
!
    DO I=1,NumEle
       IF(K(I).LE.2) THEN
          WRITE(*,*) "ERROR, ELEMENT DOES NOT OVERLAP!"
          STOP
       ENDIF
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE MFS_TRUSS_ASSD
!
    SUBROUTINE MFS_2D_ASSD(NEle1,NEle2,NumEle,NEleAss,MFS_NODE,MFS_PAR,MFS_ASS)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE FOR CALCULATING NODE ASSOCIATION LM      
! *      NumEle   :  Number of elements
! *      EleAss   :  Maximum element association number
! *      MFS_NODE :  
! *      MFS_PAR  :
! *      MFS_ASS  : 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
!
    INTEGER,INTENT(IN) :: NumEle, NEleAss, NEle1,NEle2
    INTEGER,DIMENSION(NEle1,NumEle)   :: MFS_NODE
    REAL(8),DIMENSION(NEle2,NumEle)   :: MFS_PAR
    INTEGER,DIMENSION(NEleAss,NumEle) :: MFS_ASS
!
!   ASSISTANT VARIABLES
!
    INTEGER :: I,J
    INTEGER,DIMENSION(NumEle) :: K
    INTEGER :: N,M
    INTEGER :: KK
    REAL(8) :: RAD_OUT1,RAD_OUT2,RAD2_OUT1,RAD2_OUT2
    REAL(8) :: DIST1,DIST2
    LOGICAL :: EXPRN1,EXPRN2
!
    KK=NEleAss+1
!
    DO I=1,NumEle
       K(I)=1
    ENDDO
!
    DO I=1,NumEle
       RAD_OUT1=MFS_PAR(3,I)
       RAD_OUT2=MFS_PAR(4,I)
       N=MFS_NODE(1,I)
       DO J=I,NumEle              
          M=MFS_NODE(1,J)
!
          RAD2_OUT1=MFS_PAR(3,J)
          RAD2_OUT2=MFS_PAR(4,J)
!
          DIST1=RAD_OUT1+RAD2_OUT1
          DIST1=ABS(X(N)-X(M))-DIST1
          EXPRN1=DIST1.GT.(-1E-8)
          IF(EXPRN1) CYCLE
          DIST2=RAD_OUT2+RAD2_OUT2
          DIST2=ABS(Y(N)-Y(M))-DIST2
          EXPRN2=DIST2.GT.(-1E-8)
          IF(EXPRN2) CYCLE
!          
          MFS_ASS(K(I),I)=J
          K(I)=K(I)+1
!
          IF(I.NE.J) THEN
             MFS_ASS(K(J),J)=I
             K(J)=K(J)+1
          ENDIF
          IF(K(I).GT.KK) THEN
             WRITE(*,*) "ELleAss should .GE.",K(I)
             WRITE(*,*) "N AND M ARE:"
             WRITE(*,*) N,M
             WRITE(*,*) "X(N),Y(M),RADI1,RADI2=", X(N),Y(N),RAD_OUT1,RAD_OUT2
             WRITE(*,*) "X(N),Y(M),RADJ1,RADJ2=", X(M),Y(M),RAD2_OUT1,RAD2_OUT2
             write(*,*) "LM is:"
             WRITE(*,*) MFS_NODE(1,MFS_ASS(:,I))
             WRITE(*,*) "ERROR, NEleAss too small!"
             STOP
          ENDIF
       ENDDO
    ENDDO
!
    DO I=1,NumEle
       IF(K(I).LE.4) THEN
          WRITE(*,*) "ERROR, ELEMENT DOES NOT OVERLAP!"
          STOP
       ENDIF
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE MFS_2D_ASSD 
!**************************************************
!
!   SUBROUTINES FOR SURFACE INPUT
!
!**************************************************
    SUBROUTINE SURFACE_IN(IIN,IOUT,ITEM2)
!{{{
!--------------------------------------------------
!   ALLOCATE SURFACE MEMORY AND SELECT SURFACE TYPE
!   NPAR()
!   POINT-POINT
!       NPAR(1) = 1
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT,ITEM2
!
    INTEGER :: I,J,L,M
!
    INTEGER,PARAMETER :: KwSurf=1
    CHARACTER(20),DIMENSION(KwSurf) :: KwSurface
!
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KwSurface(1)='\Pt-Pt'
!
    WRITE(IOUT,2000)
!
    IF (NSPAIR<1E-7) THEN
            WRITE(*,*) " ERROR, NSPAIR EQUALS TO ZEROS"
            STOP
    ENDIF
!
    NCQ=0
    DO J=1,NSPAIR
       I = KeyWord(KwSurface,KwSurf)
       WRITE(IOUT,1000) KwSurface(I)
       SELECT CASE(I)
       CASE(1)
          NSFPAR=0
          NSFPAR(1)=1
          CALL PTPTD(IIN,IOUT)
          WRITE(ITEM2) (NSFPAR(L),L=1,7)
          WRITE(ITEM2) ((XYZC(L,M),L=1,3),M=1,NSFPAR(2)),((NSFC(L,M),L=1,1),M=1,NSFPAR(2)),((XYZT(L,M),L=1,3),M=1,NSFPAR(3)),&
          ((NSFT(L,M),L=1,1),M=1,NSFPAR(3))
          NCQ=NCQ+NSFPAR(2)
          ! Deallocate variables
          DEALLOCATE (NSFC) 
          DEALLOCATE (NSFT)
          DEALLOCATE (XYZC)
          DEALLOCATE (XYZT)      
       CASE DEFAULT
          WRITE(*,*) "ERROR, INPUT ERROR IN ELEMENTD"
          STOP
       END SELECT 
    END DO
    RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 1000   FORMAT(/,A20)
 2000   FORMAT(//,' S U R F A C E   I N F O R M A T I O N',///)

!}}}
    END SUBROUTINE SURFACE_IN
!
    SUBROUTINE PTPTD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE PTPTD
! *      READ AND STORE POINT-POINT SURFACE PAIR
! *      PTPTD
! *         NSFPAR(1) =   1
! *         NSFPAR(2) --- POINT NUMBER IN CONTRACTOR
! *         NSFPAR(3) --- POINT NUMBER IN TARGET 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!
    INTEGER :: I,J,K
!
    INTEGER :: N,NODE
    INTEGER :: NCOUNT
    INTEGER :: AllocateStatus
!
    CHARACTER(20),DIMENSION(3) :: KWPTPT
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWPTPT(1) = 'Contactor:'
    KWPTPT(2) = 'Target:'
    KWPTPT(3) = '\END'
!
    DO I=1,3
       J = KeyWord(KWPTPT,3)
       IF(J.LE.2) THEN
         K = GetInt()
         WRITE(IOUT,2000) KWPTPT(J),K
       ENDIF
         
       SELECT CASE(J)
       CASE(1)
          NSFPAR(2)=K
!
!-------- Allocate memory
!
          ALLOCATE (NSFC(1,NSFPAR(2)),STAT = ALLocateStatus)
          IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
          ALLOCATE (XYZC(3,NSFPAR(2)),STAT = ALLocateStatus)
          IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!
          WRITE(IOUT,2050)
          N=0
          NODE=0
          NCOUNT=0
          READ(IIN,*) N,NODE
          DO WHILE(N.NE.0) 
             IF(NODE.GT.NUMNP) THEN
                    WRITE(*,*) "ERROR, NODE .GT. NUMNP IN SURFACE INPUT DATA"
                    STOP
             ENDIF
             NCOUNT=NCOUNT+1
             NSFC(1,N)=NODE
             XYZC(1,N)=X(NODE)
             XYZC(2,N)=Y(NODE)
             XYZC(3,N)=Z(NODE)
!
             WRITE(IOUT,3000) N,NODE
             READ(IIN,*) N,NODE

          END DO
          IF(NCOUNT.NE.NSFPAR(2)) THEN
                  WRITE(*,*) "ERROR, NCOUNT.NE.NSFPAR(2) IN SURFACE INPUT DATA"
                  STOP
          END IF
       CASE(2)
          NSFPAR(3)=K
!
!-------- Allocate memory
!
          ALLOCATE (NSFT(1,NSFPAR(3)),STAT = ALLocateStatus)
          IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
          ALLOCATE (XYZT(3,NSFPAR(3)),STAT = ALLocateStatus)
          IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!
          WRITE(IOUT,2050)
          N=0
          NODE=0
          NCOUNT=0
          READ(IIN,*) N,NODE
          DO WHILE(N.NE.0) 
             IF(NODE.GT.NUMNP) THEN
                    WRITE(*,*) "ERROR, NODE .GT. NUMNP IN SURFACE INPUT DATA"
                    STOP
             ENDIF
             NCOUNT=NCOUNT+1
             NSFT(1,N)=NODE
             XYZT(1,N)=X(NODE)
             XYZT(2,N)=Y(NODE)
             XYZT(3,N)=Z(NODE)
!
             WRITE(IOUT,3000) N,NODE
             READ(IIN,*) N,NODE
          END DO
          IF(NCOUNT.NE.NSFPAR(3)) THEN
                  WRITE(*,*) "ERROR, NCOUNT.NE.NSFPAR(3) IN SURFACE INPUT DATA"
                  STOP
          END IF
       CASE(3)
          EXIT
       CASE DEFAULT
          WRITE(*,*) "ERROR, INPUT ERROR IN ELEMENTD"
          STOP
       END SELECT 
!   
    END DO
    RETURN    
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 1000   FORMAT(/,A20)
 2000   FORMAT(/,A20,2X,I5)
 2050   FORMAT(/,&
              ' SURFACE     NODE ',/,&
              ' NUMBER-N      I  ',/)
 3000   FORMAT(I5,6X,I5)
!}}}
    END SUBROUTINE PTPTD
!
end module INPUTM
