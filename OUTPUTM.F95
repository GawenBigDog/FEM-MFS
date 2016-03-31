module OUTPUTM
    USE WAREHOUSE
    USE GENERATEM
    IMPLICIT NONE
    PUBLIC  :: OUTPUT
    PUBLIC  :: OUTPUTD
    PRIVATE :: OUTDISPD1
    PRIVATE :: OUTDISPD2
    PRIVATE :: OUTDISP_STRESS
    PRIVATE :: OUTDISP_FEM
    PRIVATE :: OUTDISP_MFS
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
    INTEGER :: EleType,NumEle
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
           ALLOCATE (NDMFS(2,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (ELEPAR(2,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!--------------------------------------------------
           DO K=1,NumEle
              NDMFS(1,K)=0
              NDMFS(2,K)=0
              ELEPAR(1,K)=0.
              ELEPAR(2,K)=0.
              DO J=1,NEleAss
                 LM(J,K)=0
              ENDDO
           END DO
!
!--------------------------------------------------
           READ(ADELE) ((NDMFS(K,J),K=1,2),J=1,NumEle),((ELEPAR(J,K),J=1,2),K=1,NumEle),((LM(J,K),J=1,NEleAss),K=1,NumEle)
!}}}
           CALL OUTDISP_MFS(IDISP,NEleAss,NumEle,1)
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
           ALLOCATE (NDMFS(5,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (ELEPAR(3,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!--------------------------------------------------
           DO K=1,NumEle
              DO J=1,5
                 NDMFS(J,K)=0
              ENDDO
!
              DO J=1,3
                 ELEPAR(J,K)=0.
              ENDDO
!
              DO J=1,NEleAss
                 LM(J,K)=0
              ENDDO
           END DO
!
!--------------------------------------------------
           READ(ADELE) ((NDMFS(K,J),K=1,5),J=1,NumEle),((ELEPAR(J,K),J=1,3),K=1,NumEle),((LM(J,K),J=1,NEleAss),K=1,NumEle)
!}}}
!{{{
           READ(ADBC) (NPAR(J),J=1,7)
           BcType=NPAR(1)
           SELECT CASE(BcType)
           CASE(1)
             NumBc  = NPAR(2)
!   ALLOCATE VARIABLES
             ALLOCATE (BCED(2,NumBc),STAT = ALLocateStatus)
             IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
             ALLOCATE (BCIND(4,NumBc),STAT = ALLocateStatus)
             IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
             ALLOCATE (BCPAR(2,NumBc),STAT = ALLocateStatus)
             IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!
             DO K=1,NumBc
                DO J=1,4
                  BCIND(J,K)=0
                ENDDO
                DO J=1,2
                  BCED(J,K)=0
                  BCPAR(J,K)=0.
                ENDDO
             ENDDO
!
             READ(ADBC) ((BCED(J,K),J=1,2),K=1,NumBc),((BCIND(J,K),J=1,4),K=1,NumBc),((BCPAR(J,K),J=1,2),K=1,NumBc)
!
             DO J=1,NumBc
               L1=BCED(1,J)
               L2=BCED(2,J)
               NDMFS(5,L1)=L2
               NDMFS(4,L2)=L1
             ENDDO
           CASE DEFAULT
              WRITE(*,*), "ERROR, WRONG INPUT IN GENERATE BOUNDARY PROCEDURE"
              STOP
           END SELECT
!}}}
           CALL OUTDISP_MFS(IDISP,NEleAss,NumEle,2)
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
!
!     Print displacements
!
      WRITE (IDISP,1000) TIME
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
            IL=I
            IF (KK.NE.0)  THEN
                D(IL)=DISP(KK)
            ENDIF
         END DO
!
         WRITE (IDISP,2010) II,D
      END DO
!
      RETURN
!
 1000 FORMAT (/,' TIME =',5X,F13.3)
 2000 FORMAT (' D I S P L A C E M E N T S',//,'  NODE ',10X,&
              'X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT')
 2010 FORMAT (1X,I3,8X,3ES18.6)
!}}}
    END SUBROUTINE OUTDISP_FEM
!
    SUBROUTINE OUTDISP_MFS(IDISP,NEleAss,NumEle,IND)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To print displacements                                          .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IDISP,IND
    INTEGER :: NEleAss,NumEle
!
    REAL(8) :: D(3)
    INTEGER :: IC,II,I,KK,IL
    INTEGER :: J,K
    INTEGER :: N,M,NUM
    INTEGER :: INDD,NPOR
    INTEGER :: JJ
    REAL(8) :: S,Radius
    INTEGER :: DIMEN
!
!     Print displacements
!
      WRITE (IDISP,1000) TIME
      WRITE (IDISP,2000)
      IC=4
!
      IF(IND.EQ.1) THEN
         DIMEN=1
      ELSE
         DIMEN=2
      ENDIF

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
         DO I=1,DIMEN
            N=NDMFS(1,II)
            KK=ID(I,N)
            IL=I
            IF (KK.NE.0)  THEN
               DO J=1,NEleAss
                  JJ=LM(J,II)
                  IF (JJ.NE.0) THEN
                      M=NDMFS(1,JJ)
                      IF(IND.EQ.1) THEN
                         Radius=ELEPAR(2,JJ)
                      ELSE
                         Radius=ELEPAR(3,JJ)
                      ENDIF
                      S=(X(N)-X(M))*(X(N)-X(M))+(Y(N)-Y(M))*(Y(N)-Y(M))
                      S=SQRT(S)/Radius
                      IF(S.LE.1) THEN
                         IF(NDOF(I,M).NE.0) THEN
                           DO K=1,NDOF(I,M)                             
                              INDD=NDMFS(I+1,JJ)
                              IF((INDD.EQ.1).OR.(INDD.EQ.3)) THEN
                                 NPOR=K+1
                              ELSE
                                 NPOR=K
                              ENDIF
                              NUM=ID(I,M)+K-1 
                              IF(IND.EQ.1) THEN                        
                                 D(IL)=D(IL)+DISP(NUM)*TRUSS_MFS_HIHJ(NEleAss,JJ,NPOR,X(N),1)
                              ELSE
                                 D(IL)=D(IL)+DISP(NUM)*QUADS_MFS_HIHJ(NEleAss,JJ,NPOR,X(N),Y(N),1) 
                              ENDIF
                           ENDDO
                         ENDIF
                         INDD=NDMFS(I+1,JJ)
                         IF ((INDD.EQ.1).OR.(INDD.EQ.3)) THEN
                            IF(IND.EQ.1) THEN                        
                               D(IL)=D(IL)+RESS(I,M)*TRUSS_MFS_HIHJ(NEleAss,JJ,1,X(N),1)
                            ELSE
                               D(IL)=D(IL)+RESS(I,M)*QUADS_MFS_HIHJ(NEleAss,JJ,1,X(N),Y(N),1) 
                            ENDIF
                         ENDIF
                      ENDIF
                  ENDIF
               ENDDO
            ELSE
               D(IL)=RESS(I,N)                  
            ENDIF
         END DO
         write(*,*) 
!
         WRITE (IDISP,2010) N,D
      END DO    
!
    RETURN
!
 1000 FORMAT (/,' TIME =',5X,F13.3)
 2000 FORMAT (' D I S P L A C E M E N T S',//,'  NODE ',10X,&
              'X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT')
 2010 FORMAT (1X,I3,8X,3ES18.6)

!}}}
    END SUBROUTINE OUTDISP_MFS
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
