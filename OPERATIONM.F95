        MODULE OPERATIONM
!--------------------------------------------------
!       USER DEFINED OPERATIONS
!--------------------------------------------------
        USE WAREHOUSE
        IMPLICIT NONE
! SUBROUTINES
        PUBLIC SHEPARD_2D
        PUBLIC SHEPARD_TRUSS
! FUNCTIONS
        PUBLIC MATRIXB_VECTOR
        PUBLIC QUADS_MFS_HIHJ
        PUBLIC QUADS_MFS_HIHJ_V2
! * * * * * * * * * * * * * * * * * * 
        PUBLIC TRUSS_MFS_HIHJ
        PUBLIC TRUSS_MFS_HIHJ_V2
        PUBLIC areClockwise
!
        CONTAINS
!
        FUNCTION MATRIXB_VECTOR(A,NMAXA,NMHT,V,NWKK,NDOFF) RESULT(H)
!{{{
          IMPLICIT NONE
          INTEGER :: NWKK,NDOFF
          INTEGER :: NMAXA(NDOFF+1),NMHT(NDOFF)
          REAL(8) :: A(NWKK),V(NDOFF)
          REAL(8) :: H(NDOFF)
!         ASSISTANT VARIABLES
          INTEGER :: I,J
          INTEGER :: JJ
!
          DO I=1,NDOFF
             H(I)=0.0D0
          ENDDO
!
          DO I=1,NDOFF
             JJ=NMAXA(I)
             H(I)=H(I)+A(JJ)*V(I)
             DO J=1,NMHT(I)
                JJ=JJ+1
                H(I)=H(I)+A(JJ)*V(I-J)
                H(I-J)=H(I-J)+A(JJ)*V(I)
             ENDDO
          ENDDO
!}}}
        END FUNCTION MATRIXB_VECTOR
!
        FUNCTION areClockwise(X1,X2,Y1,Y2) RESULT(REST)
!{{{
        IMPLICIT NONE
        REAL(8) :: X1,X2,Y1,Y2
        LOGICAL :: REST
!
        REST= -X1*Y2+Y1*X2>=-1e-8  
!}}}
        END FUNCTION areClockwise
!
      SUBROUTINE SHEPARD_2D(NEleAss,II,JJ,XX,YY,W,DWDX,DWDY,IND_BLEND)
!{{{
      IMPLICIT NONE
      INTEGER :: NEleAss, II,JJ
      REAL(8) :: XX,YY
      REAL(8) :: W(5),DWDX(5),DWDY(5)
      INTEGER :: IND_BLEND
!     ASSISTANT VARIABLES
      INTEGER :: I,J
      INTEGER :: M
      INTEGER :: INDSP
      REAL(8) :: RAD1,RAD2
      REAL(8) :: S,SS(4),SX,SY,S1
      REAL(8) :: CONSTX,CONSTY
      REAL(8) :: CONSTXD,CONSTYD
      REAL(8) :: WW,DWDXX,DWDYY
      LOGICAL :: EXPRN1
      REAL(8) :: W3,DWDX3,DWDY3
      INTEGER :: INDSPI,INDSPJ
!     
      DO I=1,5
         W(I)=0.0D0
         DWDX(I)=0.0D0
         DWDY(I)=0.0D0
      ENDDO
!
      W3=0.0D0
      DWDX3=0.0D0
      DWDY3=0.0D0
      DO I=1,NEleAss
         J=LM(I,II)
         IF (J.NE.0) THEN
             M=NDMFS(1,J)
             INDSP=NDMFS(6,J)
             RAD1=ELEPAR(3,J)
             RAD1=1/RAD1
!
             IF(INDSP.GE.3) THEN
                RAD2=ELEPAR(4,J)
                RAD2=1.0/RAD2
                CONSTX=(XX-X(M))
!
                SX=1.0-ABS(CONSTX)*RAD1
                EXPRN1=SX.LT.1E-8
                IF(EXPRN1) CYCLE
                CONSTY=(YY-Y(M))
                SY=1.0-ABS(CONSTY)*RAD2
                EXPRN1=SY.LT.1E-8
                IF(EXPRN1) CYCLE
             ELSE
                CONSTX=(XX-X(M))
                CONSTY=(YY-Y(M))
!
                S1=SQRT(CONSTX*CONSTX+CONSTY*CONSTY)
                S=S1*RAD1
                EXPRN1=S.GT.1.0
                IF(EXPRN1) CYCLE
                CONSTXD=CONSTX/S1*RAD1
                CONSTYD=CONSTY/S1*RAD1
             ENDIF
!
                SELECT CASE(INDSP)
                CASE(1)
                  SS(1)=S*S
                  SS(2)=SS(1)*S
                  SS(3)=SS(2)*S
                  SS(4)=-12*S+24*SS(1)-12*SS(2)
!
                  WW=1-6*SS(1)+8*SS(2)-3*SS(3)
                  DWDXX=SS(4)*CONSTXD
                  DWDYY=SS(4)*CONSTYD
                CASE(2)
                  SS(1)=S
                  SS(2)=S*S
                  SS(3)=SS(2)*S
                  IF(S.LE.0.5) THEN
                     SS(4)=-8*SS(1)+12*SS(2)
!
                     WW=2/3.0-4*SS(2)+4*SS(3)
                     DWDXX=SS(4)*CONSTXD
                     DWDYY=SS(4)*CONSTYD
                  ELSE
                     SS(4)=-4+8*SS(1)-4*SS(2)
!
                     WW=4/3.0-4*SS(1)+4*SS(2)-4/3.0*SS(3)
                     DWDXX=SS(4)*CONSTXD
                     DWDYY=SS(4)*CONSTYD
                  ENDIF
                CASE(3)
                     WW=SX*SY
                     DWDXX=-SY*RAD1*SIGN(1.0D0,CONSTX)
                     DWDYY=-SX*RAD2*SIGN(1.0D0,CONSTY)
                CASE(4)
                     IND_BLEND=1
!
                     WW=SX*SY
                     DWDXX=-SY*RAD1*SIGN(1.0D0,CONSTX)
                     DWDYY=-SX*RAD2*SIGN(1.0D0,CONSTY)
!
                     W3=W3+WW
                     DWDX3=DWDX3+DWDXX
                     DWDY3=DWDY3+DWDYY
                CASE DEFAULT
                   WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
                   STOP
                END SELECT                                 
!
             IF(II.EQ.J) THEN
                W(1)=WW
                DWDX(1)=DWDXX
                DWDY(1)=DWDYY
             ENDIF
!
             IF(JJ.EQ.J) THEN
                W(2)=WW
                DWDX(2)=DWDXX
                DWDY(2)=DWDYY
             ENDIF
!
             W(3)=W(3)+WW
             DWDX(3)=DWDX(3)+DWDXX
             DWDY(3)=DWDY(3)+DWDYY
         ENDIF
      ENDDO
!
      IF(IND_BLEND.EQ.0) RETURN
!
      INDSPI=NDMFS(6,II)
      INDSPJ=NDMFS(6,JJ)
!
      IF(INDSPI.EQ.4) THEN
         W(4)=W(1)*W3
         DWDX(4)=DWDX(1)*W3+W(1)*DWDX3
         DWDY(4)=DWDY(1)*W3+W(1)*DWDY3
      ENDIF
      IF(INDSPI.EQ.3) THEN
         M=NDMFS(1,II)
         RAD1=ELEPAR(3,II)
         RAD1=1.0/RAD1  
         RAD2=ELEPAR(4,II)
         RAD2=1.0/RAD2
         CONSTX=(XX-X(M))
         CONSTY=(YY-Y(M))
         SX=1.0-ABS(CONSTX)*RAD1
         SY=1.0-ABS(CONSTY)*RAD2
!
         W(4)=W(1)+W3*SX*SY
         DWDX(4)=DWDX(1)+DWDX3*SX*SY-W3*SY*RAD1*SIGN(1.0D0,CONSTX)
         DWDY(4)=DWDY(1)+DWDY3*SX*SY-W3*SX*RAD2*SIGN(1.0D0,CONSTY)
      ENDIF
!
      SX=1.0-ABS(CONSTX)*RAD1
      SY=1.0-ABS(CONSTY)*RAD2
      IF(INDSPJ.EQ.4) THEN
         W(5)=W(2)*W3
         DWDX(5)=DWDX(2)*W3+W(2)*DWDX3
         DWDY(5)=DWDY(2)*W3+W(2)*DWDY3
      ENDIF
      IF(INDSPJ.EQ.3) THEN
         M=NDMFS(1,JJ)
         RAD1=ELEPAR(3,JJ)
         RAD1=1/RAD1  
         RAD2=ELEPAR(4,JJ)
         RAD2=1.0/RAD2
         CONSTX=(XX-X(M))
         CONSTY=(YY-Y(M))
         SX=1.0-ABS(CONSTX)*RAD1
         SY=1.0-ABS(CONSTY)*RAD2
!   
         W(5)=W(2)+W3*SX*SY
         DWDX(5)=DWDX(2)+DWDX3*SX*SY-W3*SY*RAD1*SIGN(1.0D0,CONSTX)
         DWDY(5)=DWDY(2)+DWDY3*SX*SY-W3*SX*RAD2*SIGN(1.0D0,CONSTY)
      ENDIF
!
      RETURN
!}}}
      END SUBROUTINE SHEPARD_2D
!
      SUBROUTINE SHEPARD_TRUSS(NEleAss,II,XX,W,DWDX,IND)
!{{{
!---------------------------------------------------------------
!     W(1):  WEIGHTING FUNCTION CORRESPONDING TO NODE II
!     DWDX(1): ...
!     W(2):  SUM OF WEIGHTING FUNCTIONS
!     DWDX(2): ...
!     W(3):  SUM OF WEIGHTING FUNCTIONS FOR BLENDING ELEMENT
!     DWDX(3): ...
!---------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: NEleAss, II
      REAL(8) :: XX
      REAL(8) :: W(3),DWDX(3)
      INTEGER :: IND
!     ASSISTANT VARIABLES
      INTEGER :: I,J
      INTEGER :: M
      INTEGER :: INDSP
      REAL(8) :: Radius
      REAL(8) :: S,SS(4)
      REAL(8) :: CONSTM1,CONSTMX
      REAL(8) :: W2,DWDX2,W3,DWDX3
!
      W(2)=0.0D0
      DWDX(2)=0.0D0
      W3=0.0D0
      DWDX3=0.0D0
      DO I=1,NEleAss
         J=LM(I,II)
         IF (J.NE.0) THEN
             M=NDMFS(1,J)
             INDSP=NDMFS(3,J)
             Radius=ELEPAR(2,J)
!
             CONSTM1=1/Radius
             CONSTMX=(XX-X(M))
!
             S=ABS(CONSTMX)*CONSTM1
!
             IF(S.LE.1) THEN
                SELECT CASE(INDSP)
                CASE(1)
                  SS(1)=S*S
                  SS(2)=SS(1)*S
                  SS(3)=SS(2)*S
                  SS(4)=-12*S+24*SS(1)-12*SS(2)
!
                  W2=1-6*SS(1)+8*SS(2)-3*SS(3)
                  DWDX2=SS(4)*CONSTM1*SIGN(1.0D0,CONSTMX)
                CASE(2)
                  SS(1)=S
                  SS(2)=S*S
                  SS(3)=SS(2)*S
                  IF(S.LE.0.5) THEN
                     SS(4)=-8*SS(1)+12*SS(2)
                     W2=2/3.0-4*SS(2)+4*SS(3)
                     DWDX2=SS(4)*CONSTM1*SIGN(1.0D0,CONSTMX)
                  ELSE
                     SS(4)=-4+8*SS(1)-4*SS(2)
                     W2=4/3.0-4*SS(1)+4*SS(2)-4/3.0*SS(3)
                     DWDX2=SS(4)*CONSTM1*SIGN(1.0D0,CONSTMX) 
                  ENDIF
                CASE(3)
                  W2=1-S
                  DWDX2=-CONSTM1*SIGN(1.0D0,CONSTMX)
                CASE(4)
                  W2=1-S
                  DWDX2=-CONSTM1*SIGN(1.0D0,CONSTMX)
!
                  W3=W3+W2
                  DWDX3=DWDX3+DWDX2
                CASE DEFAULT
                   WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
                   STOP
                END SELECT
!
                IF(J.EQ.II) THEN
                   W(1)=W2
                   DWDX(1)=DWDX2
                ENDIF
!
                W(2)=W(2)+W2 
                DWDX(2)=DWDX(2)+DWDX2
             ENDIF
         ENDIF
      ENDDO
!
      INDSP=NDMFS(3,II)
      Radius=ELEPAR(2,II)
!
      IF(INDSP.EQ.4) THEN
         W(3)=W(1)*W3
         DWDX(3)=DWDX(1)*W3+W(1)*DWDX3
      ELSE
         M=NDMFS(1,II)
         CONSTM1=1/Radius
         CONSTMX=(XX-X(M))
         S=ABS(CONSTMX)*CONSTM1
!
         W(3)=W(1)+W3*(1-S)
         DWDX(3)=DWDX(1)+DWDX3*(1-S)+W3*(-CONSTM1*SIGN(1.0D0,CONSTMX))
      ENDIF
!
      RETURN
!}}}
      END SUBROUTINE SHEPARD_TRUSS
!
      SUBROUTINE QUADS_MFS_HIHJ_V2(NEleAss,II,JJ,XX,YY,NPORMAX1,NPORMAX2,HI,HJ) 
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE VALUE OF HI,HJ,DHI,DHJ                                .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: NEleAss,II,JJ,NPORMAX1,NPORMAX2
        REAL(8),INTENT(IN) :: XX,YY 
        REAL(8)            :: HI(3,NPORMAX1),HJ(3,NPORMAX2)
!   ASSISTANT VARIABLES
        INTEGER :: I,J
        INTEGER :: N,M
        REAL(8) :: W(5),DWDX(5),DWDY(5)
        REAL(8) :: RAD_IX,RAD_IY,RAD_JX,RAD_JY
        REAL(8) :: CONSTW2,CONSTW2_X,CONSTW2_Y
        REAL(8) :: CONSTX_I,CONSTX_J,CONSTY_I,CONSTY_J
        INTEGER :: IND_BLEND
        INTEGER :: INDSPI,INDSPJ
        INTEGER :: INDN,INDM
        REAL(8) :: NX,NY
        REAL(8) :: X_BAR,Y_BAR
!
        IND_BLEND=0
!
        RAD_IX=ELEPAR(3,II)
        RAD_IX=1.0/RAD_IX
!
        RAD_IY=ELEPAR(4,II)
        RAD_IY=1.0/RAD_IY
!
        RAD_JX=ELEPAR(3,JJ)
        RAD_JX=1.0/RAD_JX
!
        RAD_JY=ELEPAR(4,JJ)
        RAD_JY=1.0/RAD_JY
!
        N=NDMFS(1,II)
        M=NDMFS(1,JJ)
!
        CALL SHEPARD_2D(NEleAss,II,JJ,XX,YY,W,DWDX,DWDY,IND_BLEND)
!
        CONSTW2=1/W(3)
        CONSTW2_X=DWDX(3)*CONSTW2*CONSTW2
        CONSTW2_Y=DWDY(3)*CONSTW2*CONSTW2
!
        CONSTX_I=(XX-X(N))*RAD_IX
        CONSTY_I=(YY-Y(N))*RAD_IY
        CONSTX_J=(XX-X(M))*RAD_JX
        CONSTY_J=(YY-Y(M))*RAD_JY
!
        INDN=NDMFS(7,II)
        INDM=NDMFS(7,JJ)
        IF(INDN.EQ.1) THEN
          NX=ELEPAR(5,II)
          NY=ELEPAR(6,II)
          X_BAR=(NX*CONSTX_I+NY*CONSTY_I)
          Y_BAR=(-NY*CONSTX_I+NX*CONSTY_I)
          SELECT CASE(NPORMAX1)
!{{{
          CASE(1)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
          CASE(2)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
           HI(1,2)=HI(1,1)*X_BAR
           HI(2,2)=HI(2,1)*X_BAR+HI(1,1)*RAD_IX*NX
           HI(3,2)=HI(3,1)*X_BAR+HI(1,1)*RAD_IY*NY
          CASE(3)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
           HI(1,2)=HI(1,1)*X_BAR
           HI(2,2)=HI(2,1)*X_BAR+HI(1,1)*RAD_IX*NX
           HI(3,2)=HI(3,1)*X_BAR+HI(1,1)*RAD_IY*NY
!   HIXY,DHIXYDX,DHIXYDY
           HI(1,3)=HI(1,2)*Y_BAR
           HI(2,3)=HI(2,2)*Y_BAR+HI(1,2)*(-NY)*RAD_IX
           HI(3,3)=HI(3,2)*Y_BAR+HI(1,2)*NX*RAD_IY
         CASE(4)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
           HI(1,2)=HI(1,1)*X_BAR
           HI(2,2)=HI(2,1)*X_BAR+HI(1,1)*RAD_IX*NX
           HI(3,2)=HI(3,1)*X_BAR+HI(1,1)*RAD_IY*NY
!   HIXY,DHIXYDX,DHIXYDY
           HI(1,3)=HI(1,2)*Y_BAR
           HI(2,3)=HI(2,2)*Y_BAR+HI(1,2)*(-NY)*RAD_IX
           HI(3,3)=HI(3,2)*Y_BAR+HI(1,2)*NX*RAD_IY
!   HIXX,DHIXXDX,DHIXXDY
           HI(1,4)=HI(1,2)*X_BAR
           HI(2,4)=HI(2,2)*X_BAR+HI(1,2)*RAD_IX*NX
           HI(3,4)=HI(3,2)*X_BAR+HI(1,2)*RAD_IY*NY
          CASE DEFAULT
            WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
            STOP
!}}}
          END SELECT
        ELSE
          SELECT CASE(NPORMAX1)
!{{{
          CASE(1)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
          CASE(3)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
              HI(1,2)=HI(1,1)*CONSTX_I
              HI(2,2)=HI(2,1)*CONSTX_I+HI(1,1)*RAD_IX
              HI(3,2)=HI(3,1)*CONSTX_I
!   HIY,DHIYDX,DHIYDY
              HI(1,3)=HI(1,1)*CONSTY_I
              HI(2,3)=HI(2,1)*CONSTY_I
              HI(3,3)=HI(3,1)*CONSTY_I+HI(1,1)*RAD_IY
          CASE(4)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
              HI(1,2)=HI(1,1)*CONSTX_I
              HI(2,2)=HI(2,1)*CONSTX_I+HI(1,1)*RAD_IX
              HI(3,2)=HI(3,1)*CONSTX_I
!   HIY,DHIYDX,DHIYDY
              HI(1,3)=HI(1,1)*CONSTY_I
              HI(2,3)=HI(2,1)*CONSTY_I
              HI(3,3)=HI(3,1)*CONSTY_I+HI(1,1)*RAD_IY
!   HIXY,DHIXYDX,DHIXYDY
                 HI(1,4)=HI(1,2)*CONSTY_I
                 HI(2,4)=HI(2,2)*CONSTY_I
                 HI(3,4)=HI(3,3)*CONSTX_I
          CASE(6)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
              HI(1,2)=HI(1,1)*CONSTX_I
              HI(2,2)=HI(2,1)*CONSTX_I+HI(1,1)*RAD_IX
              HI(3,2)=HI(3,1)*CONSTX_I
!   HIY,DHIYDX,DHIYDY
              HI(1,3)=HI(1,1)*CONSTY_I
              HI(2,3)=HI(2,1)*CONSTY_I
              HI(3,3)=HI(3,1)*CONSTY_I+HI(1,1)*RAD_IY
!   HIXY,DHIXYDX,DHIXYDY
                 HI(1,4)=HI(1,2)*CONSTY_I
                 HI(2,4)=HI(2,2)*CONSTY_I
                 HI(3,4)=HI(3,3)*CONSTX_I
!   HIXX,DHIXXDX,DHIXXDY
                    HI(1,5)=HI(1,2)*CONSTX_I
                    HI(2,5)=HI(2,2)*CONSTX_I+HI(1,2)*RAD_IX
                    HI(3,5)=HI(3,2)*CONSTX_I
!   HIYY,DHIYYDX,DHIYYDY
                    HI(1,6)=HI(1,3)*CONSTY_I
                    HI(2,6)=HI(2,3)*CONSTY_I
                    HI(3,6)=HI(3,3)*CONSTY_I+HI(1,3)*RAD_IY
          CASE DEFAULT
            WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
            STOP
!}}}
          END SELECT
        ENDIF
!
        IF(INDM.EQ.1) THEN
          NX=ELEPAR(5,JJ)
          NY=ELEPAR(6,JJ)
          X_BAR=(NX*CONSTX_J+NY*CONSTY_J)
          Y_BAR=(-NY*CONSTX_J+NX*CONSTY_J)
          SELECT CASE(NPORMAX2)
!{{{
          CASE(1)
           HJ(1,1)=W(2)*CONSTW2
           HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
           HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
          CASE(2)
           HJ(1,1)=W(2)*CONSTW2
           HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
           HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
           HJ(1,2)=HJ(1,1)*X_BAR
           HJ(2,2)=HJ(2,1)*X_BAR+HJ(1,1)*RAD_JX*NX
           HJ(3,2)=HJ(3,1)*X_BAR+HJ(1,1)*RAD_JY*NY
          CASE(3)
           HJ(1,1)=W(2)*CONSTW2
           HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
           HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
           HJ(1,2)=HJ(1,1)*X_BAR
           HJ(2,2)=HJ(2,1)*X_BAR+HJ(1,1)*RAD_JX*NX
           HJ(3,2)=HJ(3,1)*X_BAR+HJ(1,1)*RAD_JY*NY
!   HIXY,DHIXYDX,DHIXYDY
           HJ(1,3)=HJ(1,2)*Y_BAR
           HJ(2,3)=HJ(2,2)*Y_BAR+HJ(1,2)*(-NY)*RAD_JX
           HJ(3,3)=HJ(3,2)*Y_BAR+HJ(1,2)*NX*RAD_JY
         CASE(4)
           HJ(1,1)=W(2)*CONSTW2
           HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
           HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
           HJ(1,2)=HJ(1,1)*X_BAR
           HJ(2,2)=HJ(2,1)*X_BAR+HJ(1,1)*RAD_JX*NX
           HJ(3,2)=HJ(3,1)*X_BAR+HJ(1,1)*RAD_JY*NY
!   HIXY,DHIXYDX,DHIXYDY
           HJ(1,3)=HJ(1,2)*Y_BAR
           HJ(2,3)=HJ(2,2)*Y_BAR+HJ(1,2)*(-NY)*RAD_JX
           HJ(3,3)=HJ(3,2)*Y_BAR+HJ(1,2)*NX*RAD_JY
!   HIXX,DHIXXDX,DHIXXDY
           HJ(1,4)=HJ(1,2)*X_BAR
           HJ(2,4)=HJ(2,2)*X_BAR+HJ(1,2)*RAD_JX*NX
           HJ(3,4)=HJ(3,2)*X_BAR+HJ(1,2)*RAD_JY*NY
          CASE DEFAULT
            WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
            STOP
!}}}
          END SELECT
        ELSE
         SELECT CASE(NPORMAX2)
!{{{
        CASE(1)
           HJ(1,1)=W(2)*CONSTW2
           HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
           HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
        CASE(3)
           HJ(1,1)=W(2)*CONSTW2
           HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
           HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
              HJ(1,2)=HJ(1,1)*CONSTX_J
              HJ(2,2)=HJ(2,1)*CONSTX_J+HJ(1,1)*RAD_JX
              HJ(3,2)=HJ(3,1)*CONSTX_J
!   HIY,DHIYDX,DHIYDY
              HJ(1,3)=HJ(1,1)*CONSTY_J
              HJ(2,3)=HJ(2,1)*CONSTY_J
              HJ(3,3)=HJ(3,1)*CONSTY_J+HJ(1,1)*RAD_JY
        CASE(4)
           HJ(1,1)=W(2)*CONSTW2
           HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
           HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
              HJ(1,2)=HJ(1,1)*CONSTX_J
              HJ(2,2)=HJ(2,1)*CONSTX_J+HJ(1,1)*RAD_JX
              HJ(3,2)=HJ(3,1)*CONSTX_J
!   HIY,DHIYDX,DHIYDY
              HJ(1,3)=HJ(1,1)*CONSTY_J
              HJ(2,3)=HJ(2,1)*CONSTY_J
              HJ(3,3)=HJ(3,1)*CONSTY_J+HJ(1,1)*RAD_JY
!   HIXY,DHIXYDX,DHIXYDY
                 HJ(1,4)=HJ(1,2)*CONSTY_J
                 HJ(2,4)=HJ(2,2)*CONSTY_J
                 HJ(3,4)=HJ(3,3)*CONSTX_J
        CASE(6)
           HJ(1,1)=W(2)*CONSTW2
           HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
           HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
              HJ(1,2)=HJ(1,1)*CONSTX_J
              HJ(2,2)=HJ(2,1)*CONSTX_J+HJ(1,1)*RAD_JX
              HJ(3,2)=HJ(3,1)*CONSTX_J
!   HIY,DHIYDX,DHIYDY
              HJ(1,3)=HJ(1,1)*CONSTY_J
              HJ(2,3)=HJ(2,1)*CONSTY_J
              HJ(3,3)=HJ(3,1)*CONSTY_J+HJ(1,1)*RAD_JY
!   HIXY,DHIXYDX,DHIXYDY
                 HJ(1,4)=HJ(1,2)*CONSTY_J
                 HJ(2,4)=HJ(2,2)*CONSTY_J
                 HJ(3,4)=HJ(3,3)*CONSTX_J
!   HIXX,DHIXXDX,DHIXXDY
                    HJ(1,5)=HJ(1,2)*CONSTX_J
                    HJ(2,5)=HJ(2,2)*CONSTX_J+HJ(1,2)*RAD_JX
                    HJ(3,5)=HJ(3,2)*CONSTX_J
!   HIYY,DHIYYDX,DHIYYDY
                    HJ(1,6)=HJ(1,3)*CONSTY_J
                    HJ(2,6)=HJ(2,3)*CONSTY_J
                    HJ(3,6)=HJ(3,3)*CONSTY_J+HJ(1,3)*RAD_JY
        CASE DEFAULT
            WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
            STOP
!}}}
         END SELECT
        ENDIF
!
        IF(IND_BLEND.EQ.0) RETURN
!
           INDSPI=NDMFS(6,II)
           INDSPJ=NDMFS(6,JJ)
           IF(INDSPI.GE.3) THEN
              HI(1,1)=W(4)*CONSTW2
              HI(2,1)=DWDX(4)*CONSTW2-W(4)*CONSTW2_X
              HI(3,1)=DWDY(4)*CONSTW2-W(4)*CONSTW2_Y
           ENDIF
!
           IF(INDSPJ.GE.3) THEN
              HJ(1,1)=W(5)*CONSTW2
              HJ(2,1)=DWDX(5)*CONSTW2-W(5)*CONSTW2_X
              HJ(3,1)=DWDY(5)*CONSTW2-W(5)*CONSTW2_Y
           ENDIF
!
        RETURN
!}}}
      END SUBROUTINE QUADS_MFS_HIHJ_V2
!
      FUNCTION QUADS_MFS_HIHJ(NEleAss,II,NPORMAX,XX,YY) RESULT(H)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE VALUE OF HI,HJ,DHI,DHJ                                .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: NEleAss,II,NPORMAX
        REAL(8),INTENT(IN) :: XX,YY 
        REAL(8) :: H(NPORMAX)
!   ASSISTANT VARIABLES
        INTEGER :: I,J
        INTEGER :: N,M
        INTEGER :: INDSP
        REAL(8) :: W(3)
        REAL(8) :: Radius1,Radius2,RAD1,RAD2
        REAL(8) :: S,SS(3),SX,SY
        REAL(8) :: CONSTX,CONSTY
        REAL(8) :: CONST1
        REAL(8) :: WW,W3
        LOGICAL :: EXPRN1
        INTEGER :: IND_BLEND
        INTEGER :: INDN
        REAL(8) :: NX,NY,X_BAR,Y_BAR
!
        IND_BLEND=0
!
        Radius1=ELEPAR(3,II)
        Radius2=ELEPAR(4,II)
!
        N=NDMFS(1,II)
        CONSTX=(XX-X(N))/Radius1
        CONSTY=(YY-Y(N))/Radius2
!
        W(1)=0.0D0
        W(2)=0.0D0
        W3=0.0D0
        DO I=1,NEleAss
           J=LM(I,II)
           IF (J.NE.0) THEN
             M=NDMFS(1,J)
             INDSP=NDMFS(6,J)
             RAD1=ELEPAR(3,J)
             RAD1=1.0/RAD1
             IF(INDSP.GE.3) THEN
               RAD2=ELEPAR(4,J)
               RAD2=1.0/RAD2
!
               SX=1.0-ABS(XX-X(M))*RAD1
               EXPRN1=SX.LT.1E-8
               IF(EXPRN1) CYCLE
               SY=1.0-ABS(YY-Y(M))*RAD2
               EXPRN1=SY.LT.1E-8
               IF(EXPRN1) CYCLE
             ELSE
               S=(XX-X(M))*(XX-X(M))+(YY-Y(M))*(YY-Y(M))
               S=SQRT(S)*RAD1
               EXPRN1=S.GT.1.0
               IF(EXPRN1) CYCLE
             ENDIF
!
               SELECT CASE(INDSP)
               CASE(1)
                  SS(1)=S*S
                  SS(2)=SS(1)*S
                  SS(3)=SS(2)*S
!
                  WW=1-6*SS(1)+8*SS(2)-3*SS(3)
               CASE(2)
                  SS(1)=S
                  SS(2)=S*S
                  SS(3)=SS(2)*S
                  IF(S.LE.0.5) THEN
                     WW=2/3.0-4*SS(2)+4*SS(3)
                  ELSE
                     WW=4/3.0-4*SS(1)+4*SS(2)-4/3.0*SS(3)
                  ENDIF
               CASE(3)
                  WW=SX*SY
               CASE(4)
                  IND_BLEND=1
                  WW=SX*SY                  
!
                  W3=W3+WW                  
               CASE DEFAULT
                  WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
                  WRITE(*,*) INDSP
                  STOP
               END SELECT 
!
             IF(J.EQ.II) W(1)=WW
             W(2)=W(2)+WW                
           ENDIF
        ENDDO
!
        CONST1=1/W(2)
        INDN=NDMFS(7,II)
        IF(INDN.EQ.1) THEN
          NX=ELEPAR(5,II)
          NY=ELEPAR(6,II)
          X_BAR=(NX*CONSTX+NY*CONSTY)
          Y_BAR=(-NY*CONSTX+NX*CONSTY)
          SELECT CASE(NPORMAX) 
          CASE(1)
           H(1)=W(1)*CONST1
          CASE(2)
           H(1)=W(1)*CONST1
           H(2)=H(1)*X_BAR
          CASE(3)
           H(1)=W(1)*CONST1
           H(2)=H(1)*X_BAR
           H(3)=H(2)*Y_BAR
          CASE(4)
           H(1)=W(1)*CONST1
           H(2)=H(1)*X_BAR
           H(3)=H(2)*Y_BAR
           H(4)=H(2)*X_BAR
          CASE DEFAULT
           WRITE(*,*) "ERROR, NPOR SHOULD LT 4!"
           STOP
          END SELECT
        ELSE
!
          SELECT CASE(NPORMAX) 
          CASE(1)
           H(1)=W(1)*CONST1
          CASE(3)
           H(1)=W(1)*CONST1
           H(2)=H(1)*CONSTX
           H(3)=H(1)*CONSTY
          CASE(4)
           H(1)=W(1)*CONST1
           H(2)=H(1)*CONSTX
           H(3)=H(1)*CONSTY
           H(4)=H(2)*CONSTY
          CASE(6)
           H(1)=W(1)*CONST1
           H(2)=H(1)*CONSTX
           H(3)=H(1)*CONSTY
           H(4)=H(2)*CONSTY
           H(5)=H(2)*CONSTX
           H(6)=H(3)*CONSTY
          CASE DEFAULT
           WRITE(*,*) "ERROR, NPOR SHOULD LT 4!"
           STOP
          END SELECT
        ENDIF
!
        IF(IND_BLEND.EQ.0) RETURN
!
        INDSP=NDMFS(6,II)
        IF(INDSP.EQ.4) W(3)=W(1)*W3
        IF(INDSP.EQ.3) THEN
           SX=1.0-ABS(CONSTX)
           SY=1.0-ABS(CONSTY)
!
           W(3)=W(1)+W3*SX*SY
        ENDIF
!
        H(1)=W(3)*CONST1   
!
        RETURN
!}}}
      END FUNCTION QUADS_MFS_HIHJ
!    
      FUNCTION TRUSS_MFS_HIHJ(NEleAss,II,NPORMAX,XX,IND) RESULT(H)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE VALUE OF HI,HJ,DHI,DHJ                                .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NEleAss
      INTEGER, INTENT(IN) :: II,NPORMAX
      INTEGER, INTENT(IN) :: IND
      REAL(8), INTENT(IN) :: XX
      REAL(8)             :: H(NPOR_MAX)
!   ASSISTANT VARIABLES
      INTEGER :: I,J
      INTEGER :: N,M
      INTEGER :: INDSP
      REAL(8) :: W(3)
      REAL(8) :: Radius1,Radius
      REAL(8) :: S,SS(3)
      REAL(8) :: CONSTX,CONST2
      REAL(8) :: W2
      REAL(8) :: W3
!
      Radius1=ELEPAR(2,II)
      N=NDMFS(1,II)
      CONSTX=(XX-X(N))/Radius1
!
        W(2)=0.0D0
        W3=0.0D0
        DO I=1,NEleAss
           J=LM(I,II)
           IF (J.NE.0) THEN
             M=NDMFS(1,J)
             INDSP=NDMFS(3,J)
             Radius=ELEPAR(2,J)
             S=ABS(XX-X(M))/Radius
             IF(S.LE.1) THEN
                SELECT CASE(INDSP)
                CASE(1)
                  SS(1)=S*S
                  SS(2)=SS(1)*S
                  SS(3)=SS(2)*S
!
                  W2=1-6*SS(1)+8*SS(2)-3*SS(3)
                CASE(2)
                  SS(1)=S
                  SS(2)=S*S
                  SS(3)=SS(2)*S
                  IF(S.LE.0.5) THEN
                     W2=2/3.0-4*SS(2)+4*SS(3)
                  ELSE
                     W2=4/3.0-4*SS(1)+4*SS(2)-4/3.0*SS(3)
                  ENDIF
                CASE(3)
                  W2=1-S
                CASE(4)
                  W2=1-S
!
                  W3=W3+W2
                CASE DEFAULT
                   WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
                   STOP
                END SELECT 
                IF(J.EQ.II) W(1)=W2
                W(2)=W(2)+W2
             ENDIF
           ENDIF
        ENDDO
!
        INDSP=NDMFS(3,II)
        S=ABS(CONSTX)
        IF(INDSP.EQ.4) THEN
           W(3)=W(1)*W3
        ELSE
           W(3)=W(1)+W3*(1-S)
        ENDIF
!
        CONST2=1/W(2)
        SELECT CASE(NPORMAX)
        CASE(1)
           H(1)=W(3)*CONST2
        CASE(2)
           H(1)=W(1)*CONST2
           H(2)=H(1)*CONSTX
           H(1)=W(3)*CONST2
        CASE(3)
           H(1)=W(1)*CONST2
           H(2)=H(1)*CONSTX
           H(3)=H(2)*CONSTX
           H(1)=W(3)*CONST2
        CASE(4)
           H(1)=W(1)*CONST2
           H(2)=H(1)*CONSTX
           H(3)=H(2)*CONSTX
           H(4)=H(3)*CONSTX
           H(1)=W(3)*CONST2
        CASE DEFAULT
             WRITE(*,*) "NPORMAX.GT.4 IN TRUSS_MFS_HIHJ"
             STOP
        END SELECT     
!
      RETURN
!}}}
      END FUNCTION TRUSS_MFS_HIHJ
!
      FUNCTION TRUSS_MFS_HIHJ_V2(NEleAss,II,NPORMAX,XX,IND) RESULT(H)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE VALUE OF HI,HJ,DHI,DHJ                                .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: NEleAss,II,NPORMAX
        REAL(8),INTENT(IN) :: XX
        REAL(8)            :: H(NPORMAX)
        INTEGER :: IND
!   ASSISTANT VARIABLES
        INTEGER :: I,J
        INTEGER :: N,M
        REAL(8) :: HH1,HH2,HH3
        REAL(8) :: YY
        REAL(8) :: W(3),DWDX(3)
        REAL(8) :: Radius1
        REAL(8) :: CONST1,CONST2,CONSTX
!
        Radius1=ELEPAR(2,II)
        CONST1=1/Radius1
        N=NDMFS(1,II)
!
        CALL SHEPARD_TRUSS(NEleAss,II,XX,W,DWDX,IND)
        CONST2=1/W(2)
!
        CONSTX=(XX-X(N))*CONST1
!
        SELECT CASE(NPORMAX)
        CASE(1)
           H(1)=DWDX(3)*CONST2-W(3)*DWDX(2)*CONST2*CONST2
        CASE(2)
           HH1=W(1)*CONST2
           H(1)=DWDX(1)*CONST2-W(1)*DWDX(2)*CONST2*CONST2           
           H(2)=H(1)*CONSTX+HH1*CONST1
           H(1)=DWDX(3)*CONST2-W(3)*DWDX(2)*CONST2*CONST2
        CASE(3)
           HH1=W(1)*CONST2
           H(1)=DWDX(1)*CONST2-W(1)*DWDX(2)*CONST2*CONST2           
           HH2=HH1*CONSTX
           H(2)=H(1)*CONSTX+HH1*CONST1
           H(3)=H(2)*CONSTX+HH2*CONST1
           H(1)=DWDX(3)*CONST2-W(3)*DWDX(2)*CONST2*CONST2
        CASE(4)
           HH1=W(1)*CONST2
           H(1)=DWDX(1)*CONST2-W(1)*DWDX(2)*CONST2*CONST2           
           HH2=HH1*CONSTX
           H(2)=H(1)*CONSTX+HH1*CONST1
           HH3=HH2*CONSTX
           H(3)=H(2)*CONSTX+HH2*CONST1
           H(4)=H(3)*CONSTX+HH3*CONST1
           H(1)=DWDX(3)*CONST2-W(3)*DWDX(2)*CONST2*CONST2
        CASE DEFAULT
             WRITE(*,*) "NPORMAX.GT.4 IN QUADS_MFS_HIHJ2"
             STOP
        END SELECT
!
        RETURN
!}}}
      END FUNCTION TRUSS_MFS_HIHJ_V2
!
      END MODULE OPERATIONM
