module OUTPUTM
    USE WAREHOUSE
    USE GENERATEM
    USE OPERATIONM
    IMPLICIT NONE
    PUBLIC  :: OUTPUT
    PUBLIC  :: OUTPUTD
    PRIVATE :: OUTDISPD1
    PRIVATE :: OUTDISPD2
    PRIVATE :: OUTDISP_STRESS
    PRIVATE :: OUTDISP_FEM
    PRIVATE :: OUTDISP_2D_MFS
    PRIVATE :: OUTDISP_TRUSS_MFS
    PRIVATE :: TRUSS_OUT
    CONTAINS
!   
    SUBROUTINE OUTPUT(IDISP,ISTRE,ADELE,ADBC)
!{{{
!----------------------------------------
!   OUTPUT DISPLACEMENT AND STRESS
!----------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IDISP,ISTRE,ADELE,ADBC
!   OUTPUT DISPLACEMENT AND STRESS
    CALL OUTDISP_STRESS(IDISP,ISTRE,ADELE,ADBC)
!
    RETURN
!}}}
    END SUBROUTINE OUTPUT
!
    SUBROUTINE OUTPUTD(IND,IDISP,NODE,IDIRN)
!{{{
!----------------------------------------
!   OUTPUT DATA FOR PLOTTING
!----------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: IND,IDISP
    INTEGER, INTENT(IN),OPTIONAL :: NODE,IDIRN
!   OUTPUT DISPLACEMENT FOR PLOTTING
    SELECT CASE(IND)
    CASE(1)
      CALL OUTDISPD1(IDISP,NODE,IDIRN)
    CASE(2)
      CALL OUTDISPD2(IDISP)
    CASE DEFAULT
      WRITE(*,*) "ERROR, WRONG INPUT IN OUTPUTD"
      STOP
    END SELECT
    RETURN
!}}}
    END SUBROUTINE OUTPUTD
!
    SUBROUTINE OUTDISPD1(IDISP,NODE,IDIRN)
!{{{
!----------------------------------------
!   OUTPUT DISPLACEMENT FOR PLOTTING
!----------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: IDISP
    INTEGER, INTENT(IN) :: NODE,IDIRN
!
    INTEGER :: M
    REAL(8) :: DD,VV,AA
!
    M=ID(IDIRN,NODE)
    IF(M.EQ.0) THEN
       DD = 0
       VV = 0
       AA = 0
    ELSE
       DD = U(M)
       VV = V(M)
       AA = ACC(M)
    ENDIF
!
    WRITE(IDISP,2000) TIME,DD,VV,AA
    RETURN
!
 2000 FORMAT(F10.5,5X,3F13.6)
!}}}
    END SUBROUTINE OUTDISPD1
!
    SUBROUTINE OUTDISPD2(IDISP)
!{{{
!----------------------------------------
!   OUTPUT DISPLACEMENT FOR PLOTTING
!----------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: IDISP
!
    INTEGER :: I
    REAL(8) :: DD,VV,AA
!
    WRITE(IDISP,1000) TIME
    DO I=1,NEQ
       DD=U(I)
       VV=V(I)
       AA=ACC(I)
       WRITE(IDISP,2000) I,DD,VV,AA
    END DO
!
    RETURN
!
 1000 FORMAT(//,"TIME=",F10.6)
 2000 FORMAT(2X,'NEQ',1X,I5,5X,ES13.6,2X,ES13.6,2X,ES13.6)
!}}}
    END SUBROUTINE OUTDISPD2
!
    SUBROUTINE OUTDISP_STRESS(IDISP,ISTRE,ADELE,ADBC)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To print displacements                                          .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IDISP,ISTRE,ADELE,ADBC 
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
!   TRUSS
    INTEGER :: EleType,NumEle,INDSP
    INTEGER :: BcType,NumBc
    INTEGER :: L1,L2
!   QUADS
!   MFS
    INTEGER :: NEleAss
!   
    INTEGER :: ALLocateStatus
!
    REWIND(ADELE)
    REWIND(ADBC)
    DO I=1,NEGROUP
       READ(ADELE) (NPAR(J),J=1,7)
       EleType=NPAR(1)
       SELECT CASE(EleType)
       CASE(1)
!--------------------------------------------------
!      TRUSS
!--------------------------------------------------
!{{{
           NumEle = NPAR(3)
!--------------------------------------------------
!   ALLOCATE VARIABLES
!--------------------------------------------------
           ALLOCATE (XYZ(6,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (LM(6,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (ELEPAR(1,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!--------------------------------------------------
           DO J=1,NumEle
             DO K=1,6
                XYZ(K,J)=0
                LM(K,J)=0
             END DO
             ELEPAR(1,J)=0.
           END DO
!--------------------------------------------------
           READ(ADELE) ((XYZ(K,J),K=1,6),J=1,NumEle),((LM(K,J),K=1,6),J=1,NumEle),(ELEPAR(1,J),J=1,NPAR(3))
!--------------------------------------------------
!}}}
           CALL OUTDISP_FEM(IDISP)
           CALL TRUSS_OUT(ISTRE,ADELE,I)
!--------------------------------------------------
           DEALLOCATE(XYZ)
           DEALLOCATE(LM)
           DEALLOCATE(ELEPAR)
       CASE(2)
!--------------------------------------------------
!      QUADS
!--------------------------------------------------
!{{{
!      READ
!--------------------------------------------------
           NumEle = NPAR(4)
!--------------------------------------------------
!      ALLOCATE VARIABLES
!--------------------------------------------------
           ALLOCATE (XYZ(8,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (LM(8,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (ELEPAR(1,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"             
!--------------------------------------------------
           DO J=1,NumEle
              DO K=1,8
                 XYZ(K,J)=0
                 LM(K,J)=0
              END DO
              ELEPAR(1,J)=0.
           END DO  
!--------------------------------------------------
           READ(ADELE) ((XYZ(K,J),K=1,8),J=1,NumEle),((LM(K,J),K=1,8),J=1,NumEle),(ELEPAR(1,K),K=1,NumEle)
!}}}
           CALL OUTDISP_FEM(IDISP)
!--------------------------------------------------
           DEALLOCATE(XYZ)
           DEALLOCATE(LM)
           DEALLOCATE(ELEPAR)
       CASE(3)
!--------------------------------------------------
!      MFS_BAR
!--------------------------------------------------
!{{{
           NEleAss= NPAR(4)
           NumEle = NPAR(3)  
!--------------------------------------------------
!   ALLOCATE VARIABLES
!--------------------------------------------------
           ALLOCATE (NDMFS(3,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (ELEPAR(2,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!--------------------------------------------------
           DO K=1,NumEle
              NDMFS(1,K)=0
              NDMFS(2,K)=0
              NDMFS(3,K)=0
              ELEPAR(1,K)=0.
              ELEPAR(2,K)=0.
              DO J=1,NEleAss
                 LM(J,K)=0
              ENDDO
           END DO
!
!--------------------------------------------------
           READ(ADELE) ((NDMFS(K,J),K=1,3),J=1,NumEle),((ELEPAR(J,K),J=1,2),K=1,NumEle),((LM(J,K),J=1,NEleAss),K=1,NumEle)
!}}}
           CALL OUTDISP_TRUSS_MFS(IDISP,NEleAss,NumEle)
           DEALLOCATE(NDMFS)
           DEALLOCATE(ELEPAR)
           DEALLOCATE(LM)
       CASE(4)
!--------------------------------------------------
!      MFS_QUADS
!--------------------------------------------------
!{{{
           NEleAss= NPAR(5)
           NumEle = NPAR(4)  
!--------------------------------------------------
!   ALLOCATE VARIABLES
!--------------------------------------------------
           ALLOCATE (NDMFS(8,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (ELEPAR(6,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!--------------------------------------------------
           DO K=1,NumEle
              DO J=1,8
                 NDMFS(J,K)=0
              ENDDO
!
              DO J=1,6
                 ELEPAR(J,K)=0.
              ENDDO
!
              DO J=1,NEleAss
                 LM(J,K)=0
              ENDDO
           END DO
!--------------------------------------------------
           READ(ADELE) ((NDMFS(K,J),K=1,8),J=1,NumEle),((ELEPAR(J,K),J=1,6),K=1,NumEle),((LM(J,K),J=1,NEleAss),K=1,NumEle)
!}}}
           CALL OUTDISP_2D_MFS(IDISP,NEleAss,NumEle)
           DEALLOCATE(NDMFS)
           DEALLOCATE(ELEPAR)
           DEALLOCATE(LM)
       CASE DEFAULT
            WRITE(*,*) "ERROR, WRONG INPUT IN OUTDISP SUBROUTINE!"
            STOP
       END SELECT
    ENDDO
!}}}
    END SUBROUTINE OUTDISP_STRESS
!
    SUBROUTINE OUTDISP_FEM(IDISP)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To print displacements                                          .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: IDISP
!
      REAL(8) :: D(3)
      INTEGER :: IC,II,I,KK,IL
      REAL(8) :: H(NEQ)
      LOGICAL :: EXPN1,EXPN2
!
!   Calculate strain energy 
!
      ESTRAIN=0
      H=MATRIXB_VECTOR(KMATRIX,MAXA_COLSOL,MHT,DISP,NWK_COLSOL,NEQ)
      DO I=1,NEQ
         ESTRAIN=ESTRAIN+DISP(I)*H(I)
      ENDDO
      ESTRAIN=ESTRAIN/2.0
!
!     Print displacements
!
      WRITE (IDISP,1000) TIME
      WRITE (IDISP,1050) ESTRAIN
      WRITE (IDISP,2000)
      IC=4
!
      DO II=1,NUMNP
         IC=IC + 1
         IF (IC.GE.56) THEN
            WRITE (IDISP,1000) TIME
            WRITE (IDISP,2000)
            IC=4
         END IF
!
         DO I=1,3
            D(I)=0.
         END DO
!
         DO I=1,3
            KK=ID(I,II)
            EXPN2=KK.NE.0
            IL=I
            EXPN1=NDOF(I,II).EQ.0
            IF (EXPN1.AND.EXPN2) THEN
                D(IL)=DISP(KK)
            ENDIF
         END DO
!
         WRITE (IDISP,2010) II,D
      END DO
!
      RETURN
!
 1000 FORMAT (/,' TIME =',5X,ES14.7) 
 1050 FORMAT (' TOTAL STRAIN ENERGY =',5X,ES18.6)
 2000 FORMAT (' D I S P L A C E M E N T S',//,'  NODE ',10X,&
              'X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT')
 2010 FORMAT (1X,I3,8X,3ES18.6)
!}}}
    END SUBROUTINE OUTDISP_FEM
!
    SUBROUTINE OUTDISP_2D_MFS(IDISP,NEleAss,NumEle)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To print displacements                                          .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IDISP
    INTEGER :: NEleAss,NumEle
!
    REAL(8) :: D(3)
    REAL(8) :: H(NEQ)
    REAL(8) :: HI(NPOR_MAX)
    INTEGER :: IC,II,I,KK,IL
    INTEGER :: I1
    INTEGER :: J,JJ,K,KMAX
    INTEGER :: N,M,NUM
    INTEGER :: INDD,NPOR
    REAL(8) :: RAD1,RAD2,SX,SY,S
    LOGICAL :: EXPRN1
    INTEGER :: INC
    INTEGER :: INDSP
!

!
!   Calculate strain energy 
!
      ESTRAIN=0
      H=MATRIXB_VECTOR(KMATRIX,MAXA_COLSOL,MHT,DISP,NWK_COLSOL,NEQ)
      DO I=1,NEQ
         ESTRAIN=ESTRAIN+DISP(I)*H(I)
      ENDDO
      ESTRAIN=ESTRAIN/2.0
!
!     Print displacements
!
      WRITE (IDISP,1000) TIME
      WRITE (IDISP,1050) ESTRAIN
      WRITE (IDISP,2000)
      IC=4
!
      DO II=1,NumEle
         IC=IC + 1
         IF (IC.GE.56) THEN
            WRITE (IDISP,1000) TIME
            WRITE (IDISP,2000)
            IC=4
         END IF
!
         DO I=1,3
            D(I)=0.
         END DO
!
         DO I=1,2
            I1=I+1
            N=NDMFS(1,II)
            INDSP=NDMFS(6,II)
            KK=NDMFS(I1,II)
            SELECT CASE(KK)
            CASE(1)
               D(I)=RESS(I,N)
            CASE(2)
               NUM=ID(I,N)
               D(I)=DISP(NUM)
            CASE(3)
               NUM=ID(I,N)
               D(I)=DISP(NUM)
            CASE(0)
               DO J=1,NEleAss
                  JJ=LM(J,II)
                  IF(JJ.EQ.0) CYCLE
                  M=NDMFS(1,JJ)
                  INDSP=NDMFS(6,JJ)
!
                  RAD1=ELEPAR(3,JJ)
                  RAD1=1.0/RAD1
                  IF(INDSP.EQ.3) THEN
                     RAD2=ELEPAR(4,JJ)
                     RAD2=1.0/RAD2
!
                     SX=1.0-ABS(X(N)-X(M))*RAD1
                     EXPRN1=SX.LT.1E-8
                     IF(EXPRN1) CYCLE
                     SY=1.0-ABS(Y(N)-Y(M))*RAD2
                     EXPRN1=SY.LT.1E-8
                     IF(EXPRN1) CYCLE
                  ELSE
                     S=(X(N)-X(M))*(X(N)-X(M))+(Y(N)-Y(M))*(Y(N)-Y(M))
                     S=SQRT(S)*RAD1
                     EXPRN1=S.GT.1.0
                     IF(EXPRN1) CYCLE
                  ENDIF
!
                  INDD=NDMFS(I+1,JJ)
                  IF(INDD.NE.1) THEN 
                     INC=0
                  ELSE
                     INC=1
                  ENDIF
                  KMAX=NDOF(I,M)
                  HI=QUADS_MFS_HIHJ(NEleAss,JJ,KMAX+INC,X(N),Y(N))                      
                  IF(INDD.EQ.1) D(I)=D(I)+RESS(I,M)*HI(1)
                  DO K=1,KMAX                             
                     NUM=ID(I,M)+K-1
                     NPOR=K+INC
                     D(I)=D(I)+DISP(NUM)*HI(K) 
                  ENDDO
               ENDDO
            CASE DEFAULT 
               WRITE(*,*) "ERROR, IN OUTDISP_MFS!"
               STOP
            END SELECT
         ENDDO 
!
         WRITE (IDISP,2010) N,D
      END DO    
!
    RETURN
!
 1000 FORMAT (/,' TIME =',5X,F13.3)
 1050 FORMAT (' TOTAL STRAIN ENERGY =',5X,ES18.6)
 2000 FORMAT (' D I S P L A C E M E N T S',//,'  NODE ',10X,&
              'X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT')
 2010 FORMAT (1X,I3,8X,3ES18.6)

!}}}
    END SUBROUTINE OUTDISP_2D_MFS
!
    SUBROUTINE OUTDISP_TRUSS_MFS(IDISP,NEleAss,NumEle)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To print displacements                                          .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IDISP
    INTEGER :: NEleAss,NumEle
!
    REAL(8) :: D(3)
    REAL(8) :: H(NEQ)
    REAL(8) :: HI(NPOR_MAX)
    INTEGER :: IC,II,I,KK,IL
    INTEGER :: I1
    INTEGER :: J,JJ,K,KMAX
    INTEGER :: N,M,NUM
    INTEGER :: INDD,NPOR
    REAL(8) :: S,Radius
    INTEGER :: INC
    INTEGER :: INDSP
!
!   Calculate strain energy 
!
      ESTRAIN=0
      H=MATRIXB_VECTOR(KMATRIX,MAXA_COLSOL,MHT,DISP,NWK_COLSOL,NEQ)
      DO I=1,NEQ
         ESTRAIN=ESTRAIN+DISP(I)*H(I)
      ENDDO
      ESTRAIN=ESTRAIN/2.0
!
!     Print displacements
!
      WRITE (IDISP,1200) 
      WRITE (IDISP,1000) TIME
      WRITE (IDISP,1050) ESTRAIN
      WRITE (IDISP,2000)
      IC=4
!
      DO II=1,NumEle
         IC=IC + 1
         IF (IC.GE.56) THEN
            WRITE (IDISP,1000) TIME
            WRITE (IDISP,2000)
            IC=4
         END IF
!
         DO I=1,3
            D(I)=0.
         END DO
!
            N=NDMFS(1,II)
            KK=NDMFS(2,II)
            SELECT CASE(KK)
            CASE(1)
               D(1)=RESS(1,N)
            CASE(2)
               NUM=ID(1,N)
               D(1)=DISP(NUM)
            CASE(3)
               NUM=ID(1,N)
               D(1)=DISP(NUM)
            CASE(0)
               DO J=1,NEleAss
                  JJ=LM(J,II)
                  IF(JJ.EQ.0) CYCLE
                     M=NDMFS(1,JJ)
                     Radius=ELEPAR(2,JJ)
                     S=ABS(X(N)-X(M))
                     S=S/Radius
                     IF(S.LE.1) THEN
                        INDD=NDMFS(2,JJ)
                        IF(INDD.NE.1) THEN 
                           INC=0
                        ELSE
                           INC=1
                        ENDIF
                        KMAX=NDOF(1,M)
                        INDSP=NDMFS(3,JJ)
                        HI=TRUSS_MFS_HIHJ(NEleAss,JJ,KMAX+INC,X(N),INDSP)       
                        IF(INDD.EQ.1) D(1)=D(1)+RESS(1,M)*HI(1)
                        DO K=1,KMAX                             
                             NUM=ID(1,M)+K-1
                             NPOR=K+INC
                             D(1)=D(1)+DISP(NUM)*HI(K) 
                        ENDDO
                     ENDIF
               ENDDO
            CASE DEFAULT 
               WRITE(*,*) "ERROR, IN OUTDISP_MFS!"
               STOP
            END SELECT
!
         WRITE (IDISP,2010) N,D
      END DO    
!
    RETURN
!
 1200 FORMAT (/,'TRUSS(METHOD OF FINITE SPHERES)')
 1000 FORMAT (/,' TIME =',5X,F13.3)
 1050 FORMAT (' TOTAL STRAIN ENERGY =',5X,ES18.6)
 2000 FORMAT (' D I S P L A C E M E N T S',//,'  NODE ',10X,&
              'X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT')
 2010 FORMAT (1X,I3,8X,3ES18.6)

!}}}
    END SUBROUTINE OUTDISP_TRUSS_MFS
!
    SUBROUTINE TRUSS_OUT(ISTRE,ADELE,NG)
!{{{
!--------------------------------------------------
!   CALCULATE AND OUTPUT STRESS OF TRUSS ELEMENT
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ISTRE,ADELE,NG
!
    INTEGER :: NumEle,Mtype
!   ASSISTANT VARIABLES
    REAL(8) :: E,DEN,AREA
    INTEGER :: ALLocateStatus
!   VARIABLES FOR STRESS CALCULATION
    INTEGER :: IPRINT
    INTEGER :: N,L,I,J
    REAL(8) :: XL2,D(3),ST(6),STR,P

!   READ 
    NumEle = NPAR(3)
    Mtype  = NPAR(2)
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    DEN  = MATP(3,Mtype)
!
!   CALCULATE AND OUTPUT STRESS
!   
      IPRINT=0
      DO N=1,NumEle
         AREA=ELEPAR(1,N)
         IPRINT=IPRINT + 1
         IF (IPRINT.GT.50) IPRINT=1
         IF (IPRINT.EQ.1) THEN
            WRITE (ISTRE,1000) TIME
            WRITE (ISTRE,2060) NG
         END IF
         XL2=0.
         DO L=1,3
            D(L) = XYZ(L,N) - XYZ(L+3,N)
            XL2=XL2 + D(L)*D(L)
         END DO
         DO L=1,3
            ST(L)=(D(L)/XL2)*E
            ST(L+3)=-ST(L)
         END DO
         STR=0.0         
         DO L=1,3
            I=LM(L,N)
            IF (I.GT.0) STR=STR + ST(L)*DISP(I)
            J=LM(L+3,N)
            IF (J.GT.0) STR=STR + ST(L+3)*DISP(J)
         END DO
         P=STR*AREA
         WRITE (ISTRE,2070) N,P,STR
      END DO 
      RETURN
!
 1000 FORMAT (/,' TIME =',5X,F13.3)
 2060 FORMAT (' S T R E S S  C A L C U L A T I O N S  F O R  ',&
              'E L E M E N T  G R O U P',I4,//,&
              '  ELEMENT',13X,'FORCE',12X,'STRESS',/,&
              '  NUMBER',/)
 2070 FORMAT (1X,I5,11X,E13.6,4X,E13.6)
!}}}
    END SUBROUTINE TRUSS_OUT
end module OUTPUTM
