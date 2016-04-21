        MODULE OPERATIONM
!--------------------------------------------------
!       USER DEFINED OPERATIONS
!--------------------------------------------------
        USE WAREHOUSE
        IMPLICIT NONE
! SUBROUTINES
        PUBLIC SHEPARD
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
      SUBROUTINE SHEPARD(NEleAss,II,JJ,XX,YY,W,DWDX,DWDY,IND)
!{{{
      IMPLICIT NONE
      INTEGER :: NEleAss, II,JJ
      REAL(8) :: XX,YY
      REAL(8) :: W(3),DWDX(3),DWDY(3)
      INTEGER :: IND
!     ASSISTANT VARIABLES
      INTEGER :: I,J
      INTEGER :: M
      REAL(8) :: Radius
      REAL(8) :: S,SS(4)
      REAL(8) :: CONSTM1,CONSTM2,CONSTMX,CONSTMY,CONSTMX3,CONSTMY3
      REAL(8) :: WW,DWDXX,DWDYY
!
      W(3)=0.0D0
      DWDX(3)=0.0D0
      DWDY(3)=0.0D0
      DO I=1,NEleAss
         J=LM(I,II)
         IF (J.NE.0) THEN
             M=NDMFS(1,J)
             Radius=ELEPAR(3,J)
!
             CONSTM1=1/Radius
             CONSTM2=CONSTM1*CONSTM1
             CONSTMX=(XX-X(M))
             CONSTMY=(YY-Y(M))
!
             S=CONSTMX*CONSTMX+CONSTMY*CONSTMY
             S=SQRT(S)*CONSTM1
             CONSTMX3=CONSTMX/S
             CONSTMY3=CONSTMY/S
!
             IF(S.LE.1) THEN
                SELECT CASE(IND)
                CASE(1)
                  SS(1)=S*S
                  SS(2)=SS(1)*S
                  SS(3)=SS(2)*S
                  SS(4)=-12*S+24*SS(1)-12*SS(2)
!
                  WW=1-6*SS(1)+8*SS(2)-3*SS(3)
                  DWDXX=SS(4)*CONSTMX3*CONSTM2
                  DWDYY=SS(4)*CONSTMY3*CONSTM2
!
                  IF(II.EQ.J) THEN
                     W(1)=WW
                     DWDX(1)=DWDXX
                     DWDY(1)=DWDYY
                  ENDIF
                  IF(JJ.EQ.J) THEN
                     W(2)=WW
                     DWDX(2)=DWDXX
                     DWDY(2)=DWDYY
                  ENDIF
!
                  W(3)=W(3)+WW
                  DWDX(3)=DWDX(3)+DWDXX
                  DWDY(3)=DWDY(3)+DWDYY
                CASE(2)
                  SS(1)=S
                  SS(2)=S*S
                  SS(3)=SS(2)*S
                  IF(S.LE.0.5) THEN
                     SS(4)=-8*SS(1)+12*SS(2)
!
                     WW=2/3.0-4*SS(2)+4*SS(3)
                     DWDXX=SS(4)*CONSTMX3*CONSTM2
                     DWDYY=SS(4)*CONSTMY3*CONSTM2
!
                     IF(II.EQ.J) THEN
                        W(1)=WW
                        DWDX(1)=DWDXX
                        DWDY(1)=DWDYY
                     ENDIF
                     IF(JJ.EQ.J) THEN
                        W(2)=WW
                        DWDX(2)=DWDXX
                        DWDY(2)=DWDYY
                     ENDIF
!
                     W(3)=W(3)+WW
                     DWDX(3)=DWDX(3)+DWDXX
                     DWDY(3)=DWDY(3)+DWDYY
                  ELSE
                     SS(4)=-4+8*SS(1)-4*SS(2)
!
                     WW=4/3.0-4*SS(1)+4*SS(2)-4/3.0*SS(3)
                     DWDXX=SS(4)*CONSTMX3*CONSTM2
                     DWDYY=SS(4)*CONSTMY3*CONSTM2
!
                     IF(II.EQ.J) THEN
                        W(1)=WW
                        DWDX(1)=DWDXX
                        DWDY(1)=DWDYY
                     ENDIF
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
                CASE DEFAULT
                   WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
                   STOP
                END SELECT 
             ENDIF
         ENDIF
      ENDDO
!
      RETURN
!}}}
      END SUBROUTINE SHEPARD
!
      SUBROUTINE SHEPARD1(NEleAss,II,XX,W,DWDX,IND)
!{{{
      IMPLICIT NONE
      INTEGER :: NEleAss, II
      REAL(8) :: XX
      REAL(8) :: W(2),DWDX(2)
      INTEGER :: IND
!     ASSISTANT VARIABLES
      INTEGER :: I,J
      INTEGER :: N,M
      REAL(8) :: Radius1,Radius
      REAL(8) :: S1,S,SS1(4),SS(4)
      REAL(8) :: CONST1,CONSTX
      REAL(8) :: CONSTM1,CONSTMX
!     
      Radius1=ELEPAR(1,II)
      N=NDMFS(1,II)
      CONST1=1/Radius1
      CONSTX=(XX-X(N))
      S1=ABS(CONSTX)*CONST1  
!
      SELECT CASE(IND)
      CASE(1)
         SS1(1)=S1*S1
         SS1(2)=SS1(1)*S1
         SS1(3)=SS1(2)*S1
         SS1(4)=-12*S1+24*SS1(1)-12*SS1(2)
!   CALCULATE W1,W2,DW1,DW2
         W(1)  = 1-6*SS1(1)+8*SS1(2)-3*SS1(3)
         DWDX(1)=SS1(4)*CONST1*SIGN(1.0D0,CONSTX)
      CASE(2)
         SS1(1)=S1
         SS1(2)=S1*S1
         SS1(3)=SS1(2)*S1
         IF(S1.LE.0.5) THEN
            SS1(4)=-8*SS1(1)+12*SS1(2)
            W(1)=2/3.0-4*SS1(2)+4*SS1(3)
            DWDX(1)=SS1(4)*CONST1*SIGN(1.0D0,CONSTX)
         ELSE
            SS1(4)=-4+8*SS1(1)-4*SS1(2)
            W(1)=4/3.0-4*SS1(1)+4*SS1(2)-4/3.0*SS1(3)
            DWDX(1)=SS1(4)*CONST1*SIGN(1.0D0,CONSTX)
         ENDIF
      CASE DEFAULT
          WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
          STOP
      END SELECT
!
      W(2)=0.0D0
      DWDX(2)=0.0D0
      DO I=1,NEleAss
         J=LM(I,II)
         IF (J.NE.0) THEN
             M=NDMFS(1,J)
             Radius=ELEPAR(1,J)
!
             CONSTM1=1/Radius
             CONSTMX=(XX-X(M))
!
             S=ABS(CONSTMX)*CONSTM1
!
             IF(S.LE.1) THEN
                SELECT CASE(IND)
                CASE(1)
                  SS(1)=S*S
                  SS(2)=SS(1)*S
                  SS(3)=SS(2)*S
                  SS(4)=-12*S+24*SS(1)-12*SS(2)
!
                  W(2)=W(2)+1-6*SS(1)+8*SS(2)-3*SS(3)
                  DWDX(2)=DWDX(2)+SS(4)*CONSTM1*SIGN(1.0D0,CONSTMX)
                CASE(2)
                  SS(1)=S
                  SS(2)=S*S
                  SS(3)=SS(2)*S
                  IF(S.LE.0.5) THEN
                     SS(4)=-8*SS(1)+12*SS(2)
                     W(2)=W(2)+2/3.0-4*SS(2)+4*SS(3)
                     DWDX(2)=DWDX(2)+SS(4)*CONSTM1*SIGN(1.0D0,CONSTMX)
                  ELSE
                     SS(4)=-4+8*SS(1)-4*SS(2)
                     W(2)=W(2)+4/3.0-4*SS(1)+4*SS(2)-4/3.0*SS(3)
                     DWDX(2)=DWDX(2)+SS(4)*CONSTM1*SIGN(1.0D0,CONSTMX) 
                  ENDIF
                CASE DEFAULT
                   WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
                   STOP
                END SELECT 
             ENDIF
         ENDIF
      ENDDO
!
      RETURN
!}}}
      END SUBROUTINE SHEPARD1
!
      SUBROUTINE QUADS_MFS_HIHJ_V2(NEleAss,II,JJ,XX,YY,NPORMAX1,NPORMAX2,IND,HI,HJ) 
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE VALUE OF HI,HJ,DHI,DHJ                                .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: NEleAss,II,JJ,NPORMAX1,NPORMAX2
        REAL(8),INTENT(IN) :: XX,YY 
        REAL(8)            :: HI(NPORMAX1,3),HJ(NPORMAX2,3)
        INTEGER :: IND
!   ASSISTANT VARIABLES
        INTEGER :: I,J
        INTEGER :: N,M
        REAL(8) :: W(3),DWDX(3),DWDY(3)
        REAL(8) :: RAD_I,RAD_J
        REAL(8) :: CONSTR_I,CONSTR_J,CONSTX_I,CONSTX_J,CONSTY_I,CONSTY_J
        REAL(8) :: CONSTW2,CONSTW21,CONSTW2_X,CONSTW2_Y
!
        RAD_I=ELEPAR(3,II)
        CONSTR_I=1/RAD_I
        N=NDMFS(1,II)
        RAD_J=ELEPAR(3,JJ)
        CONSTR_J=1/RAD_J
        M=NDMFS(1,JJ)
!
        CALL SHEPARD(NEleAss,II,JJ,XX,YY,W,DWDX,DWDY,IND)
        CONSTW2=1/W(3)
        CONSTW21=W(3)*CONSTW2*CONSTW2
        CONSTW2_X=DWDX(3)*CONSTW2*CONSTW2
        CONSTW2_Y=DWDY(3)*CONSTW2*CONSTW2
!
        CONSTX_I=XX-X(N)
        CONSTY_I=YY-Y(N)
        CONSTX_J=XX-X(M)
        CONSTY_J=YY-Y(M)
!   HI,DHIDX,DHIDY
        IF(NPORMAX1.GE.1) THEN
           HI(1,1)=W(1)*CONSTW2
           HI(1,2)=DWDX(1)*CONSTW21-W(1)*CONSTW2_X
           HI(1,3)=DWDY(1)*CONSTW21-W(1)*CONSTW2_Y
           IF(NPORMAX1.GE.3) THEN
!   HIX,DHIXDX,DHIXDY
              HI(2,1)=HI(1,1)*CONSTX_I*CONSTR_I
              HI(2,2)=HI(1,2)*CONSTX_I*CONSTR_I+HI(1,1)*CONSTR_I
              HI(2,3)=HI(1,3)*CONSTX_I*CONSTR_I
!   HIY,DHIYDX,DHIYDY
              HI(3,1)=HI(1,1)*CONSTY_I*CONSTR_I
              HI(3,2)=HI(1,2)*CONSTY_I*CONSTR_I
              HI(3,3)=HI(1,3)*CONSTY_I*CONSTR_I+HI(1,1)*CONSTR_I
              IF(NPORMAX1.GE.4) THEN
!   HIXY,DHIXYDX,DHIXYDY
                 HI(4,1)=HI(2,1)*CONSTY_I*CONSTR_I
                 HI(4,2)=HI(2,2)*CONSTY_I*CONSTR_I
                 HI(4,3)=HI(3,3)*CONSTX_I*CONSTR_I
                 IF(NPORMAX1.GE.6) THEN
!   HIXX,DHIXXDX,DHIXXDY
                    HI(5,1)=HI(2,1)*CONSTX_I*CONSTR_I
                    HI(5,2)=(HI(2,2)*CONSTX_I+HI(2,1))*CONSTR_I
                    HI(5,3)=HI(2,3)*CONSTX_I*CONSTR_I
!   HIYY,DHIYYDX,DHIYYDY
                    HI(6,1)=HI(3,1)*CONSTY_I*CONSTR_I
                    HI(6,2)=HI(3,2)*CONSTY_I*CONSTR_I
                    HI(6,3)=(HI(3,3)*CONSTY_I+HI(3,1))*CONSTR_I
                    IF(NPORMAX1.GT.6) THEN
                       WRITE(*,*) "NPORMAX.GT.6 IN QUADS_MFS_HIHJ2"
                       STOP
                    ENDIF
                  ENDIF
               ENDIF
            ENDIF
        ENDIF
!   HJ,DHJDX,DHJDY
        IF(NPORMAX2.GE.1) THEN
           HJ(1,1)=W(2)*CONSTW2
           HJ(1,2)=DWDX(2)*CONSTW21-W(2)*CONSTW2_X
           HJ(1,3)=DWDY(2)*CONSTW21-W(2)*CONSTW2_Y
           IF(NPORMAX2.GE.3) THEN
!   HIX,DHIXDX,DHIXDY
              HJ(2,1)=HJ(1,1)*CONSTX_J*CONSTR_J
              HJ(2,2)=HJ(1,2)*CONSTX_J*CONSTR_J+HJ(1,1)*CONSTR_J
              HJ(2,3)=HJ(1,3)*CONSTX_J*CONSTR_J
!   HIY,DHIYDX,DHIYDY
              HJ(3,1)=HJ(1,1)*CONSTY_J*CONSTR_J
              HJ(3,2)=HJ(1,2)*CONSTY_J*CONSTR_J
              HJ(3,3)=HJ(1,3)*CONSTY_J*CONSTR_J+HJ(1,1)*CONSTR_J
              IF(NPORMAX2.GE.4) THEN
!   HIXY,DHIXYDX,DHIXYDY
                 HJ(4,1)=HJ(2,1)*CONSTY_J*CONSTR_J
                 HJ(4,2)=HJ(2,2)*CONSTY_J*CONSTR_J
                 HJ(4,3)=HJ(3,3)*CONSTX_J*CONSTR_J
                 IF(NPORMAX2.GE.6) THEN
!   HIXX,DHIXXDX,DHIXXDY
                    HJ(5,1)=HJ(2,1)*CONSTX_J*CONSTR_J
                    HJ(5,2)=(HJ(2,2)*CONSTX_J+HJ(2,1))*CONSTR_J
                    HJ(5,3)=HJ(2,3)*CONSTX_J*CONSTR_J
!   HIYY,DHIYYDX,DHIYYDY
                    HJ(6,1)=HJ(3,1)*CONSTY_J*CONSTR_J
                    HJ(6,2)=HJ(3,2)*CONSTY_J*CONSTR_J
                    HJ(6,3)=(HJ(3,3)*CONSTY_J+HJ(3,1))*CONSTR_J
                    IF(NPORMAX2.GT.6) THEN
                       WRITE(*,*) "NPORMAX.GT.6 IN QUADS_MFS_HIHJ2"
                       STOP
                    ENDIF
                  ENDIF
               ENDIF
            ENDIF
        ENDIF
!
        RETURN
!}}}
      END SUBROUTINE QUADS_MFS_HIHJ_V2
!
      FUNCTION QUADS_MFS_HIHJ(NEleAss,II,NPORMAX,XX,YY,IND) RESULT(H)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE VALUE OF HI,HJ,DHI,DHJ                                .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: NEleAss,II,NPORMAX
        REAL(8),INTENT(IN) :: XX,YY 
        INTEGER :: IND
        REAL(8) :: H(NPORMAX)
!   ASSISTANT VARIABLES
        INTEGER :: I,J
        INTEGER :: N,M
        REAL(8) :: W(2)
        REAL(8) :: Radius1,Radius
        REAL(8) :: S1,S,SS1(3),SS(3)
        REAL(8) :: CONST1,CONSTX,CONSTY
        REAL(8) :: CONSTM1
        REAL(8) :: CONST2
!
!      
        Radius1=ELEPAR(3,II)
        N=NDMFS(1,II)
        CONSTX=(XX-X(N))
        CONSTY=(YY-Y(N))
        CONST1=1/Radius1
        S1=CONSTX*CONSTX+CONSTY*CONSTY
        S1=SQRT(S1)*CONST1
!
        SELECT CASE(IND)
        CASE(1)
           SS1(1)=S1*S1
           SS1(2)=SS1(1)*S1
           SS1(3)=SS1(2)*S1
!   CALCULATE W1,W2,DW1,DW2
           W(1)  = 1-6*SS1(1)+8*SS1(2)-3*SS1(3)
        CASE(2)
           SS1(1)=S1
           SS1(2)=S1*S1
           SS1(3)=SS1(2)*S1
           IF(S1.LE.0.5) THEN
              W(1)=2/3.0-4*SS1(2)+4*SS1(3)
           ELSE
              W(1)=4/3.0-4*SS1(1)+4*SS1(2)-4/3.0*SS1(3)
           ENDIF
        CASE DEFAULT
            WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
            STOP
        END SELECT
!
        W(2)=0.0D0
        DO I=1,NEleAss
           J=LM(I,II)
           IF (J.NE.0) THEN
             M=NDMFS(1,J)
             Radius=ELEPAR(3,J)
             CONSTM1=1/Radius
             S=(XX-X(M))*(XX-X(M))+(YY-Y(M))*(YY-Y(M))
             S=SQRT(S)*CONSTM1
             IF(S.LE.1) THEN
                SELECT CASE(IND)
                CASE(1)
                  SS(1)=S*S
                  SS(2)=SS(1)*S
                  SS(3)=SS(2)*S
!
                  W(2)=W(2)+1-6*SS(1)+8*SS(2)-3*SS(3)
                CASE(2)
                  SS(1)=S
                  SS(2)=S*S
                  SS(3)=SS(2)*S
                  IF(S.LE.0.5) THEN
                     W(2)=W(2)+2/3.0-4*SS(2)+4*SS(3)
                  ELSE
                     W(2)=W(2)+4/3.0-4*SS(1)+4*SS(2)-4/3.0*SS(3)
                  ENDIF
                CASE DEFAULT
                   WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
                   STOP
                END SELECT 
             ENDIF
           ENDIF
        ENDDO
!
        CONST2=1/W(2)
        IF(NPORMAX.GE.1) THEN        
           H(1)=W(1)*CONST2
           IF(NPORMAX.GE.3) THEN
              H(2)=H(1)*CONSTX*CONST1
              H(3)=H(1)*CONSTY*CONST1
              IF(NPORMAX.GE.4) THEN
                 H(4)=H(2)*CONSTY*CONST1
                 IF(NPORMAX.GE.6) THEN
                    H(5)=H(2)*CONSTX*CONST1
                    H(6)=H(3)*CONSTY*CONST1
                    IF(NPORMAX.GT.6) THEN
                       WRITE(*,*) "ERROR, NPOR SHOULD LT 4!"
                       STOP
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
        ENDIF      
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
      REAL(8) :: W(2)
      REAL(8) :: Radius1,Radius
      REAL(8) :: S1,SS1(3),S,SS(3)
      REAL(8) :: CONSTX,CONSTN1,CONSTM1,CONST2
!
      Radius1=ELEPAR(1,II)
      N=NDMFS(1,II)
      CONSTX=XX-X(N)
      CONSTN1=1/Radius1
      S1=ABS(CONSTX)*CONSTN1
!   CALCULATE W1,W2,DW1,DW2
        SELECT CASE(IND)
        CASE(1)
           SS1(1)=S1*S1
           SS1(2)=SS1(1)*S1
           SS1(3)=SS1(2)*S1
!   CALCULATE W1,W2,DW1,DW2
           W(1)  = 1-6*SS1(1)+8*SS1(2)-3*SS1(3)
        CASE(2)
           SS1(1)=S1
           SS1(2)=S1*S1
           SS1(3)=SS1(2)*S1
           IF(S1.LE.0.5) THEN
              W(1)=2/3.0-4*SS1(2)+4*SS1(3)
           ELSE
              W(1)=4/3.0-4*SS1(1)+4*SS1(2)-4/3.0*SS1(3)
           ENDIF
        CASE DEFAULT
            WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
            STOP
        END SELECT
!
        W(2)=0.0D0
        DO I=1,NEleAss
           J=LM(I,II)
           IF (J.NE.0) THEN
             M=NDMFS(1,J)
             Radius=ELEPAR(1,J)
             CONSTM1=1/Radius
             S=(XX-X(M))
             S=ABS(S)*CONSTM1
             IF(S.LE.1) THEN
                SELECT CASE(IND)
                CASE(1)
                  SS(1)=S*S
                  SS(2)=SS(1)*S
                  SS(3)=SS(2)*S
!
                  W(2)=W(2)+1-6*SS(1)+8*SS(2)-3*SS(3)
                CASE(2)
                  SS(1)=S
                  SS(2)=S*S
                  SS(3)=SS(2)*S
                  IF(S.LE.0.5) THEN
                     W(2)=W(2)+2/3.0-4*SS(2)+4*SS(3)
                  ELSE
                     W(2)=W(2)+4/3.0-4*SS(1)+4*SS(2)-4/3.0*SS(3)
                  ENDIF
                CASE DEFAULT
                   WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
                   STOP
                END SELECT 
             ENDIF
           ENDIF
        ENDDO
!
        CONST2=1/W(2)
        IF(NPORMAX.GE.1) THEN
           H(1)=W(1)*CONST2
           IF(NPORMAX.GE.2) THEN
              H(2)=H(1)*CONSTX*CONSTN1
              IF(NPORMAX.GE.3) THEN
                 H(3)=H(2)*CONSTX*CONSTN1
                 IF(NPORMAX.GE.4) THEN
                    H(4)=H(3)*CONSTX*CONSTN1
                    IF(NPORMAX.GT.4) THEN
                       WRITE(*,*) "NPORMAX.GT.4 IN TRUSS_MFS_HIHJ"
                       STOP
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
        ENDIF      
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
        REAL(8) :: W(2),DWDX(2),DWDY(2)
        REAL(8) :: Radius1
        REAL(8) :: CONST1,CONST2,CONSTX
!
        Radius1=ELEPAR(1,II)
        CONST1=1/Radius1
        N=NDMFS(1,II)
!
        CALL SHEPARD1(NEleAss,II,XX,W,DWDX,IND)
        CONST2=1/W(2)
!
        CONSTX=XX-X(N)
!   DHIDX
        IF(NPORMAX.GE.1) THEN
           HH1=W(1)*CONST2
           H(1)=(DWDX(1)*W(2)-W(1)*DWDX(2))*CONST2*CONST2
           IF(NPORMAX.GE.2) THEN
!   DHIXDX
              HH2=HH1*CONSTX*CONST1
              H(2)=H(1)*CONSTX*CONST1+HH1*CONST1
              IF(NPORMAX.GE.3) THEN
!   DHIXXDX
                 HH3=HH2*CONSTX*CONST1
                 H(3)=H(2)*CONSTX*CONST1+HH2*CONST1
                 IF(NPORMAX.GE.4) THEN
!   DHIXXXDX
                    H(4)=H(3)*CONSTX*CONST1+HH3*CONST1
                    IF(NPORMAX.GT.4) THEN
                       WRITE(*,*) "NPORMAX.GT.4 IN QUADS_MFS_HIHJ2"
                       STOP
                    ENDIF
                  ENDIF
               ENDIF
            ENDIF
        ENDIF
!
        RETURN
!}}}
      END FUNCTION TRUSS_MFS_HIHJ_V2
!
      END MODULE OPERATIONM
