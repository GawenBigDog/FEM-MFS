program DISCOVERY
    !--------------------------------------------------
    !   FEM PROGRAM (MAIN)
    !--------------------------------------------------
    USE INTERPRETER
    USE FILEM
    USE INPUTM
    USE GENERATEM
    USE SOLVERM
    USE OUTPUTM
!    USE ODEM
    IMPLICIT NONE
    !
    INTEGER :: IIN
    INTEGER :: IOUT,ADBC,ADELE,ADSURF,ADK,ADM
    INTEGER :: ADDISP,ADSTRE,ADDISPD
    !
    INTEGER :: NUM_STEP
    !
    REAL(8) :: TT(5)
    INTEGER,PARAMETER :: NumKW=6
    CHARACTER(20),DIMENSION(NumKW) :: KeyCtrl
    !   ASSISTANT VARIABLES
    INTEGER :: I,J
    !--------------------------------------------------
    !   VARIABLES FOR DYNAMIC SOLVER
    !--------------------------------------------------
    LOGICAL :: OUTPUT_ODE
    LOGICAL :: HALF_ODE
    !
    !------------------------
    !   Variables definition
    !------------------------
    KeyCtrl(1) = '\START'
    KeyCtrl(2) = '\INPUT'
    KeyCtrl(3) = '\GENERATE'
    KeyCtrl(4) = '\STATIC_SOLVER'
    KeyCtrl(5) = '\DYNAMIC_SOLVER'
    KeyCtrl(6) = '\END'
    !   OPEN FILES
    CALL OPENFILES()  
    !   READ CONTROL FILE ADDRESS
    IIN=InputAD()
    CALL ReadFile(IIN)  !LineNumber&WordNumber are resetted to 0
    !
    I=KeyWord(KeyCtrl,NumKW)
    DO WHILE(I<=NumKW) 
        SELECT CASE(I)
         CASE(1)
                IOUT=OutputAD()
                ADBC=OutputBC()
                ADELE=OutputEle()
                ADSURF=OutputSurf()
                ADK=OutputK()
                ADM=OutputMa()
                !
                ADDISP =OutputDisp()
                ADDISPD=OutputDispD()
                ADSTRE =OutputStress()
         CASE(2)
                CALL CPU_TIME(TT(1))
                WRITE(*,*) "INPUT PROCESSING"
                WRITE(*,*) "*********************"
                CALL INPUT(IIN,IOUT,ADBC,ADELE,ADSURF)
                CALL ALLOCATEM()
                CALL COL_ALLOCATE()
         CASE(3)
                CALL CPU_TIME(TT(2))
                WRITE(*,*) "GENERATE PROCESSING"
                WRITE(*,*) "*********************"
                IF (DYNAMIC.EQ.0) THEN
                    CALL GENERATE(ADK,ADM,ADELE,ADBC)
                ELSE
                    IF (IND_INITIAL.EQ.0) THEN
                        WRITE(*,*) "ERROR, INITIAL CONDITIONS REQUIRED!"
                        STOP
                    ENDIF
                    CALL GENERATE(ADK,ADM,ADELE,ADBC)
                ENDIF
         CASE(4)
         ! ---- STATIC_SOLVER
         !{{
                CALL CPU_TIME(TT(3))
                WRITE(*,*) "STATIC SOLVER"  
                WRITE(*,*) "*********************"  
                CALL COLSOL_PRE()
                CALL COLSOL(K_COLSOL(1),DISP(1),MAXA_COLSOL(1),NEQ,NWK_COLSOL,NEQ+1,1)
                CALL COLSOL(K_COLSOL(1),DISP(1),MAXA_COLSOL(1),NEQ,NWK_COLSOL,NEQ+1,2)
                CALL CPU_TIME(TT(4))
                WRITE(*,*) "OUTPUT PROCESSING"
                WRITE(*,*) "*********************"
                CALL OUTPUT(ADDISP,ADSTRE,ADELE,ADBC)
         !}}}
         CASE(5)
           WRITE(*,*) "ERROR, ODEM IS NOT INCLUDED YET!"
           STOP
!         ! ---- DYNAMIC_SOLVER
!         !{{{
!                IF (DYNAMIC.NE.1) THEN
!                    WRITE(*,*) "ERROR, DYNAMIC IS .NE. 1"
!                    STOP
!                ENDIF
!                CALL ODE_IN(IOUT,R(1),RODE(1),ACC(1))
!                OUTPUT_ODE=.FALSE.
!                HALF_ODE=.FALSE.
!                TIME=TPARA(1)
!                DTIME=TPARA(3)
!                NUM_STEP=0
!                WRITE(*,2000) NUM_STEP,DTIME,TIME
!                IF (CONTACT.EQ.0) THEN
!                  WRITE(ADDISPD,1000)
!                  CALL OUTPUTD(IND=2,IDISP=ADDISPD)
!                  DO WHILE(.NOT.OUTPUT_ODE)
!                    IF((TPARA(2)-TIME).GT.2.0*DTIME) THEN
!                      CONTINUE
!                    ELSEIF (((TPARA(2)-TIME).GT.DTIME).AND.(.NOT.HALF_ODE)) THEN
!                      DTIME=0.5*DTIME
!                      HALF_ODE=.TRUE.
!                    ELSE
!                      DTIME=TPARA(2)-TIME
!                      OUTPUT_ODE=.TRUE.
!                    END IF
!                    TIME=TIME+DTIME
!                    !
!                    CALL ODE_SOLVER(KMATRIX(1),MMATRIX(1),R(1),NRF(1),U(1),V(1),ACC(1),TIME,DTIME)
!                    !
!                    CALL OUTPUTD(IND=2,IDISP=ADDISPD)
!                    NUM_STEP=NUM_STEP+1
!                    WRITE(*,2000) NUM_STEP,DTIME,TIME
!                  END DO
!                ELSE
!                  CONTINUE
!                ENDIF
!         !}}}
         CASE(6)
              CALL CLOSEFILES
              CALL RELEASE
              CALL CPU_TIME(TT(5))
              WRITE(*,*)  "---------------------------------------------------"
              WRITE(*,'(A20,ES14.7)')  "INPUT TIME =" ,TT(2)-TT(1)
              WRITE(*,*)  "---------------------------------------------------"
              WRITE(*,'(A20,ES14.7)')  "GENERATION TIME =" ,TT(3)-TT(2)
              WRITE(*,*)  "---------------------------------------------------"
              WRITE(*,'(A20,ES14.7)')  "SOLUTION TIME =" ,TT(4)-TT(3)
              WRITE(*,*)  "---------------------------------------------------"
              WRITE(*,'(A20,ES14.7)')  "OUTPUT TIME =" ,TT(5)-TT(4)
              WRITE(*,*)  "---------------------------------------------------"
              STOP
         CASE DEFAULT
            WRITE(*,*) "ERROR, WRONG INPUT IN CONTROL FILE"
            STOP
         END SELECT
         I=KeyWord(KeyCtrl,NumKW)
    END DO
!
1000 FORMAT(6X,"TIME",6X,"DISPLACEMENT",5X,"VELOCITY",3X,"ACCELERATION",&
    /,4X,"X_DATA",12X,"Y_DATA",6X,"Y1_DATA",6X,"Y2_DATA")
2000 FORMAT("STEP NUMBER =",I5,5X,"TIME STEP SIZE =",ES10.3,5X,"SIMULATION TIME =",F10.3)
end program DISCOVERY
