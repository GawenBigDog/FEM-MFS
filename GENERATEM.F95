module GENERATEM
    USE WAREHOUSE
    USE SOLVERM
    USE GAUSSM
    USE OPERATIONM
    implicit none
    PUBLIC  GENERATE
!   GENERATE LOAD VECTOR(MFS)
    PRIVATE  MFS_LINE_G
    PRIVATE  MFS_LINE_STIFF_WEAK
    PRIVATE  MFS_LINE_STIFF_PENALTY
    PRIVATE  MFS_LINE_ADDBAN
!   GENERATE STIFF MATRIX(FEM)
    PRIVATE TRUSS_G
    PRIVATE QUADS_G
!   GENERATE STIFF MATRIX(MFS)
    PRIVATE TRUSS_MFS_G
    PRIVATE TRUSS_MFS_STIFF
    PRIVATE TRUSS_MFS_BODY
!
    PRIVATE QUADS_MFS_G
    PRIVATE QUADS_MFS_STIFF_V2
    PRIVATE QUADS_MFS_GAUSSG_V2
    PRIVATE QUADS_MFS_ADDBAN_V2
! GENERATE CONSTITUTIVE MATRIX
    PRIVATE G_KCS 
!    
    PRIVATE ADDBAN
    PRIVATE ADDBAN_ZLB
    PRIVATE OUTPUTKMR
!
    CONTAINS
!
    SUBROUTINE GENERATE(ADK,ADM,ADELE,ADBC)
!{{{
!--------------------------------------------------
!   GENERATE 
!     SELECT GENERATOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADK,ADM,ADELE,ADBC   
!
    INTEGER :: I,J,BcType,EleType
!   ASSISTANT VARIABLES
    INTEGER :: ALLocateStatus
!-------------------------------------------------------------
    REWIND(ADELE)
    REWIND(ADBC)
!
    DO I=1,NEGROUP
       READ(ADELE) (NPAR(J),J=1,7)
       EleType=NPAR(1)
       SELECT CASE(EleType)
       CASE(1)
            CALL TRUSS_G(ADELE)
       CASE(2)
            CALL QUADS_G(ADELE)    
       CASE(3)
            CALL TRUSS_MFS_G(ADELE)
       CASE(4)
            CALL QUADS_MFS_G(ADELE,ADBC)
       CASE DEFAULT
            WRITE(*,*), "ERROR, WRONG INPUT IN GENERATE ELEMENT PROCEDURE"
            STOP
       END SELECT       
    END DO
!
!     CALL OUTPUTKMR(ADK,ADM)
!
    RETURN
!}}}
    END SUBROUTINE GENERATE
!
    SUBROUTINE TRUSS_G(ADELE)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      MASS      MATRIX
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE
    INTEGER :: NumEle,Mtype   
!   VARIABLES FOR K & M
    INTEGER :: N,L,KL,K
    REAL(8) :: XL2,XL,D(3),XX,YY
    REAL(8) :: ST(6)
    REAL(8) :: S(21)
    REAL(8) :: SMASS(21)
!   ASSISTANT VARIABLE FOR K*U
    REAL(8) :: FORCE
!   ASSISTANT VARIABLES
    INTEGER :: I,J
    REAL(8) :: E,DEN,AREA
    INTEGER :: ALLocateStatus
!   READ 
    NumEle = NPAR(3)
    Mtype  = NPAR(2)
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    DEN  = MATP(3,Mtype)
!   ALLOCATE VARIABLES
    ALLOCATE (XYZ(6,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (LM(6,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(1,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"  
!--------------------------------------------------
        DO I=1,NumEle
           DO J=1,6
              XYZ(J,I)=0
              LM(J,I)=0
           END DO
           ELEPAR(1,I)=0.
        END DO
!--------------------------------------------------
    READ(ADELE) ((XYZ(I,J),I=1,6),J=1,NumEle),((LM(I,J),I=1,6),J=1,NumEle),(ELEPAR(1,I),I=1,NumEle)
!   GENERATE
!   STIFFNESS MATRIX
    DO N=1,NumEle
       AREA=ELEPAR(1,N)
       XL2=0.
       DO L=1,3
          D(L)=XYZ(L,N) - XYZ(L+3,N)
          XL2=XL2 + D(L)*D(L)
       END DO
       XL=SQRT(XL2)
       XX=E*AREA*XL
       DO L=1,3
          ST(L)=D(L)/XL2
          ST(L+3)=-ST(L)
       END DO
!
       KL=0
       DO L=1,6
          YY=ST(L)*XX
          DO K=L,6
             KL=KL + 1
             S(KL)=ST(K)*YY                
          END DO
       END DO  
!  
       CALL ADDBAN(S(1),21,LM(1,N),6,1)
!   MASS MATRIX: CONSISTENT MASS 
       DO L=1,21
          SMASS(L)=0
       END DO
       SMASS(1)  = AREA*DEN*XL/3
       SMASS(4)  = AREA*DEN*XL/6
       SMASS(7)  = AREA*DEN*XL/3
       SMASS(10) = AREA*DEN*XL/6
       SMASS(12) = AREA*DEN*XL/3
       SMASS(15) = AREA*DEN*XL/6
       SMASS(16) = AREA*DEN*XL/3
       SMASS(19) = AREA*DEN*XL/3
       SMASS(21) = AREA*DEN*XL/3
       CALL ADDBAN(SMASS(1),21,LM(1,N),6,2)
!   CALCULATE K*U
       IF(DYNAMIC.EQ.1) THEN
         FORCE=0        
         DO L=1,3
            I=LM(L,N)
            IF ((I.GT.0).AND.(U(I).NE.0)) FORCE = FORCE + ST(L)*U(I)*E*AREA
            J=LM(L+3,N)
            IF ((J.GT.0).AND.(U(J).NE.0)) FORCE = FORCE + ST(L+3)*U(J)*E*AREA
         END DO 
         IF(FORCE.NE.0) THEN     
           DO L=1,3
              I=LM(L,N)
              IF(I.GT.0) RODE(I)=RODE(I)+FORCE*D(L)/XL
              J=LM(L+3,N)
              IF(J.GT.0) RODE(J)=RODE(J)-FORCE*D(L)/XL
           END DO
         END IF
       END IF
    END DO
!
    DEALLOCATE(XYZ)
    DEALLOCATE(LM)
    DEALLOCATE(ELEPAR)
    RETURN
!}}}
    END SUBROUTINE TRUSS_G
!
    SUBROUTINE QUADS_G(ADELE)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      MASS      MATRIX
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE
    INTEGER :: ITYPE,Mtype,NumEle,NINTG  
!   VARIABLES FOR K & M
    INTEGER :: N
    REAL(8) :: S(8,8)=0
    REAL(8) :: SMASS(8,8)=0
!   ASSISTANT VARIABLE FOR K*U
    REAL(8) :: FORCE
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,KK
    REAL(8) :: E,PR,DEN,THIC
    INTEGER :: ALLocateStatus
!
    REAL(8), DIMENSION(2,4) :: XX
!   READ 
    ITYPE  = NPAR(2)
    Mtype  = NPAR(3)
    NumEle = NPAR(4)
    NINTG  = NPAR(5)
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    PR   = MATP(2,Mtype)
    DEN  = MATP(3,Mtype)
!   ALLOCATE VARIABLES
    ALLOCATE (XYZ(8,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (LM(8,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(1,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"  
!--------------------------------------------------
        DO I=1,NumEle
           DO J=1,8
              XYZ(J,I)=0
              LM(J,I)=0
           END DO
           ELEPAR(1,I)=0.
        END DO
!--------------------------------------------------
    READ(ADELE) ((XYZ(I,J),I=1,8),J=1,NumEle),((LM(I,J),I=1,8),J=1,NumEle),(ELEPAR(1,I),I=1,NumEle)
!   GENERATE
!   STIFFNESS MATRIX
    DO N=1,NumEle
       THIC=ELEPAR(1,N)
       DO I=1,8
          DO J=1,8
             S(I,J)=0.
             SMASS(I,J)=0.
          ENDDO
       ENDDO
       DO I=1,4
          XX(1,I)=XYZ(2*I-1,N)
          XX(2,I)=XYZ(2*I,N)
       ENDDO
       CALL QUADS(NumEle,ITYPE,NINTG,THIC,E,PR,XX,S,SMASS,DEN)     
       CALL ADDBAN_ZLB(8,S,LM(:,N),1)  
!   MASS MATRIX: LUMPED MASS  
       DO I=1,8
          DO J=1,8
             IF (J.NE.I) THEN
               SMASS(I,I)=SMASS(I,I)+SMASS(I,J)
               SMASS(I,J)=0
             ENDIF
          ENDDO
       ENDDO
       CALL ADDBAN_ZLB(8,SMASS,LM(:,N),2)
!   CALCULATE K*U
       IF(DYNAMIC.EQ.1) THEN
          DO I=1,8
             FORCE=0
             KK=LM(I,N)
             IF(KK.GT.0) THEN
               DO J=1,8
                 K=LM(J,N)
                 IF(K.GT.0) THEN
                    FORCE=FORCE+S(I,J)*U(K)
                 ENDIF
               ENDDO
               RODE(KK)=RODE(KK)+FORCE
             ENDIF
          ENDDO 
       END IF
    END DO
!
    DEALLOCATE(XYZ)
    DEALLOCATE(LM)
    DEALLOCATE(ELEPAR)
    RETURN
!}}}
    END SUBROUTINE QUADS_G
!
    SUBROUTINE QUADS(NEL,ITYPE,NINT,THIC,YM,PR,XX,S,SMASS,DEN)              
!{{{                                                                      
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
! .                                                                   . 
! .   P R O G R A M                                                   . 
! .        TO CALCULATE ISOPARAMETRIC QUADRILATERAL ELEMENT STIFFNESS . 
! .        MATRIX FOR AXISYMMETRIC, PLANE STRESS, AND PLANE STRAIN    . 
! .        CONDITIONS                                                 . 
! .                                                                   . 
! .  - - INPUT VARIABLES - -                                          . 
! .        NEL       = NUMBER OF ELEMENT                              . 
! .        ITYPE     = ELEMENT TYPE                                   . 
! .                        EQ.0 = AXISYMMETRIC                        . 
! .                        EQ.1 = PLANE STRAIN                        . 
! .                        EQ.2 = PLANE STRESS                        . 
! .        NINT      = GAUSS NUMERICAL INTEGRATION ORDER              . 
! .        THIC      = THICKNESS OF ELEMENT                           . 
! .        YM        = YOUNG'S MODULUS                                . 
! .        PR        = POISSON'S RATIO                                . 
! .        XX(2,4)   = ELEMENT NODE COORDINATES                       . 
! .        S(8,8)    = STORAGE FOR STIFFNESS MATRIX                   . 
! .                                                                   . 
! .  - - OUTPUT - -                                                   . 
! .        S(8,8)    = CALCULATED STIFFNESS MATRIX                    . 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      IMPLICIT NONE                               
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
! .   THIS PROGRAM IS USED IN SINGLE PRECISION ARITHMETIC ON CRAY     . 
! .   EQUIPMENT AND DOUBLE PRECISION ARITHMETIC ON IBM MACHINES,      . 
! .   ENGINEERING WORKSTATIONS AND PCS. DEACTIVATE ABOVE LINE FOR     . 
! .   SINGLE PRECISION ARITHMETIC.                                    . 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      INTEGER :: NEL,ITYPE,NINT
      REAL(8) :: DEN
      REAL(8) :: YM,PR,XBAR,THIC
      REAL(8) :: F,G,H
      REAL(8) :: A,RI,SI,DET,WT,STIFF
      INTEGER :: I,J,K,IST,LX,LY,L,M
      REAL(8) :: D(4,4),B(4,8),HH(8,8),XX(2,4),S(8,8),SMASS(8,8),XG(4,4),WGT(4,4),DB(4)     
!                                                                       
!     MATRIX XG STORES GAUSS - LEGENDRE SAMPLING POINTS                 
!                                                                       
      DATA XG/   0.D0,   0.D0,   0.D0,   0.D0,   -.5773502691896D0,&    
      .5773502691896D0,   0.D0,   0.D0,   -.7745966692415D0,   0.D0,&   
      .7745966692415D0,   0.D0,   -.8611363115941D0,&                   
      -.3399810435849D0,   .3399810435849D0,   .8611363115941D0 /       
!                                                                       
!     MATRIX WGT STORES GAUSS - LEGENDRE WEIGHTING FACTORS              
!                                                                       
      DATA WGT /  2.D0,   0.D0,   0.D0,   0.D0,   1.D0,   1.D0,&        
      0.D0,   0.D0,   .5555555555556D0,   .8888888888889D0,&            
      .5555555555556D0,   0.D0,   .3478548451375D0,   .6521451548625D0,&
      .6521451548625D0,   .3478548451375D0 /                            
!                                                                       
!     O B T A I N  S T R E S S - S T R A I N  L A W                     
!                                                                       
      F=YM/(1.+PR)                                                      
      G=F*PR/(1.-2.*PR)                                                 
      H=F + G                                                           
!                                                                       
!     PLANE STRAIN ANALYSIS                                             
!                                                                       
      D(1,1)=H                                                          
      D(1,2)=G                                                          
      D(1,3)=0.                                                         
      D(2,1)=G                                                          
      D(2,2)=H                                                          
      D(2,3)=0.                                                         
      D(3,1)=0.                                                         
      D(3,2)=0.                                                         
      D(3,3)=F/2.                                                       
      IF (ITYPE.EQ.1) THEN                                              
      THIC=1.                                                           
      GO TO 20                                                          
      ENDIF                                                             
!                                                                       
!     AXISYMMETRIC ANALYSIS                                             
!                                                                       
      D(1,4)=G                                                          
      D(2,4)=G                                                          
      D(3,4)=0.                                                         
      D(4,1)=G                                                          
      D(4,2)=G                                                          
      D(4,3)=0.                                                         
      D(4,4)=H                                                          
      IF (ITYPE.EQ.0) GO TO 20                                          
!                                                                       
!     FOR PLANE STRESS ANALYSIS CONDENSE STRESS-STRAIN MATRIX           
!                                                                       
      DO 10 I=1,3                                                       
      A=D(I,4)/D(4,4)                                                   
      DO 10 J=I,3                                                       
      D(I,J)=D(I,J) - D(4,J)*A                                          
   10 D(J,I)=D(I,J)                                                     
!                                                                       
!     C A L C U L A T E  E L E M E N T  S T I F F N E S S               
!                                                                       
   20 DO 30 I=1,8                                                       
      DO 30 J=1,8                                                       
      S(I,J)=0.
   30 SMASS(I,J)=0. 
!   
      IST=3                                                             
      IF (ITYPE.EQ.0) IST=4                                             
      DO 80 LX=1,NINT                                                   
      RI=XG(LX,NINT)                                                    
      DO 80 LY=1,NINT                                                   
      SI=XG(LY,NINT)                                                    
!                                                                       
!     EVALUATE DERIVATIVE OPERATOR B AND THE JACOBIAN DETERMINANT DET   
!                                                                       
      CALL STDM (XX,B,HH,DET,RI,SI,XBAR,NEL,ITYPE)                         
!                                                                       
!     ADD CONTRIBUTION TO ELEMENT STIFFNESS                             
!                                                                       
      IF (ITYPE.GT.0) XBAR=THIC                                         
      WT=WGT(LX,NINT)*WGT(LY,NINT)*XBAR*DET                             
      DO 70 J=1,8                                                       
      DO 40 K=1,IST                                                     
      DB(K)=0.0                                                         
      DO 40 L=1,IST                                                     
   40 DB(K)=DB(K) + D(K,L)*B(L,J)                                       
      DO 60 I=J,8                                                       
      STIFF=0.0                                                         
      DO 50 L=1,IST                                                     
   50 STIFF=STIFF + B(L,I)*DB(L)                                        
   60 S(I,J)=S(I,J) + STIFF*WT                                          
   70 CONTINUE
       DO L=1,8
         DO M=1,8
            SMASS(L,M)=SMASS(L,M)+DEN*HH(L,M)*WT
         ENDDO
       ENDDO
   80 CONTINUE                                                          
!
      DO J=1,8                                                       
         DO I=J,8                                                       
           S(J,I)=S(I,J)                                                     
         ENDDO
      ENDDO
!   LUMPED MASS MATRIX
      RETURN                                                            
!                                                                       
      END                                                               
      SUBROUTINE STDM (XX,B,HH,DET,R,S,XBAR,NEL,ITYPE)                     
!                                                                       
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
! .                                                                   . 
! .   P R O G R A M                                                   . 
! .     TO EVALUATE THE STRAIN-DISPLACEMENT TRANSFORMATION MATRIX B   . 
! .     AT POINT (R,S) FOR A QUADRILATERAL ELEMENT                    . 
! .                                                                   . 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      IMPLICIT NONE
      INTEGER :: NEL,ITYPE
      INTEGER :: I,J,K2,K
      REAL(8) :: R,S,RP,SP,RM,SM
      REAL(8) :: DET,XBAR,DUM     
      REAL(8) :: XX(2,4),B(4,8),HH(8,8),H(4),P(2,4),XJ(2,2),XJI(2,2)             
!                                                                       
      RP = 1.0 + R                                                      
      SP = 1.0 + S                                                      
      RM = 1.0 - R                                                      
      SM = 1.0 - S                                                      
!                                                                       
!     INTERPOLATION FUNCTIONS                                           
!                                                                       
      H(1) = 0.25* RP* SP                                               
      H(2) = 0.25* RM* SP                                               
      H(3) = 0.25* RM* SM                                               
      H(4) = 0.25* RP* SM                                               
!
      HH=0.
      DO I=1,8
         If(mod(I,2) == 0) then
           J=I/2
           HH(I,2)=H(J)*H(1)
           HH(I,4)=H(J)*H(2)
           HH(I,6)=H(J)*H(3)
           HH(I,8)=H(J)*H(4)
         ELSE
           J=(I+1)/2
           HH(I,1)=H(J)*H(1)
           HH(I,3)=H(J)*H(2)
           HH(I,5)=H(J)*H(3)
           HH(I,7)=H(J)*H(4)
         ENDIF
      ENDDO
! 
!     NATURAL COORDINATE DERIVATIVES OF THE INTERPOLATION FUNCTIONS     
!                                                                       
!        1. WITH RESPECT TO R                                           
!                                                                       
      P(1,1) = 0.25* SP                                                 
      P(1,2) = - P(1,1)                                                 
      P(1,3) = - 0.25* SM                                               
      P(1,4) = - P(1,3)                                                 
!                                                                       
!        2. WITH RESPECT TO S                                           
!                                                                       
      P(2,1) = 0.25* RP                                                 
      P(2,2) = 0.25* RM                                                 
      P(2,3) = - P(2,2)                                                 
      P(2,4) = - P(2,1)                                                 
!                                                                       
!     EVALUATE THE JACOBIAN MATRIX AT POINT (R,S)                       
!                                                                       
   10 DO 30 I=1,2                                                       
      DO 30 J=1,2                                                       
      DUM = 0.0                                                         
      DO 20 K=1,4                                                       
   20 DUM=DUM + P(I,K)*XX(J,K)                                          
   30 XJ(I,J)=DUM                                                       
!                                                                       
!     COMPUTE THE DETERMINANT OF THE JACOBIAN MATRIX AT POINT (R,S)     
!                                                                       
      DET = XJ(1,1)* XJ(2,2) - XJ(2,1)* XJ(1,2)                         
      IF (DET.GT.0.00000001) GO TO 40                                   
      WRITE (*,2000) NEL                                                
      GO TO 800                                                         
!                                                                       
!     COMPUTE INVERSE OF THE JACOBIAN MATRIX                            
!                                                                       
   40 DUM=1./DET                                                        
      XJI(1,1) = XJ(2,2)* DUM                                           
      XJI(1,2) =-XJ(1,2)* DUM                                           
      XJI(2,1) =-XJ(2,1)* DUM                                           
      XJI(2,2) = XJ(1,1)* DUM                                           
!                                                                       
!     EVALUATE GLOBAL DERIVATIVE OPERATOR B                             
!                                                                       
      K2=0                                                              
      DO 60 K=1,4                                                       
      K2=K2 + 2                                                         
      B(1,K2-1) = 0.                                                    
      B(1,K2  ) = 0.                                                    
      B(2,K2-1) = 0.                                                    
      B(2,K2  ) = 0.                                                    
      DO 50 I=1,2                                                       
      B(1,K2-1) = B(1,K2-1) + XJI(1,I) * P(I,K)                         
   50 B(2,K2  ) = B(2,K2  ) + XJI(2,I) * P(I,K)                         
      B(3,K2  ) = B(1,K2-1)                                             
   60 B(3,K2-1) = B(2,K2  )                                             
!                                                                       
!     IN CASE OF PLANE STRAIN OR PLANE STRESS ANALYSIS DO NOT INCLUDE   
!     THE NORMAL STRAIN COMPONENT                                       
!                                                                       
      IF (ITYPE.GT.0) GO TO 900                                         
!                                                                       
!     COMPUTE THE RADIUS AT POINT (R,S)                                 
!                                                                       
      XBAR=0.0                                                          
      DO 70 K=1,4                                                       
   70 XBAR=XBAR + H(K)*XX(1,K)                                          
!                                                                       
!     EVALUATE THE HOOP STRAIN-DISPLACEMENT RELATION                    
!                                                                       
      IF (XBAR.GT.0.00000001) GO TO 90                                  
!                                                                       
!     FOR THE CASE OF ZERO RADIUS EQUATE RADIAL TO HOOP STRAIN          
!                                                                       
      DO 80 K=1,8                                                       
   80 B(4,K)=B(1,K)                                                     
      GO TO 900                                                         
!                                                                       
!     NON-ZERO RADIUS                                                   
!                                                                       
   90 DUM=1./XBAR                                                       
      K2=0                                                              
      DO 100 K=1,4                                                      
      K2=K2 + 2                                                         
      B(4,K2  ) = 0.                                                    
  100 B(4,K2-1) = H(K)*DUM                                              
      GO TO 900                                                         
!                                                                       
  800 STOP                                                              
  900 RETURN                                                            
!                                                                       
 2000 FORMAT (//,' *** ERROR *** ',&                                     
         ' ZERO OR NEGATIVE JACOBIAN DETERMINANT FOR ELEMENT (',I8,')')
!}}}                                                                     
    END 
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *   
!  
!     SUBROUTINES FOR TRUSS_MFS
!
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!{{{
!
    SUBROUTINE TRUSS_MFS_G(ADELE)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      MASS      MATRIX
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE      
!   ASSISTANT VARIABLES
    INTEGER :: N,M
    INTEGER :: I,J,K,L,II,JJ
    INTEGER :: NumEle,Mtype,NEleAss
    REAL(8) :: E,DEN,AREA
    INTEGER :: ALLocateStatus
!   ASSISTANT VARIABLES
    INTEGER :: NUM,NUM1
    INTEGER :: INDSP
    INTEGER :: KI
    REAL(8) :: Radius
    REAL(8) :: S        
    REAL(8) :: HI(NPOR_MAX),DHI(NPOR_MAX)
    REAL(8) :: XX,YY
    INTEGER :: NPOR
!   READ 
    NEleAss= NPAR(4)
    NumEle = NPAR(3)
    Mtype  = NPAR(2)
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    DEN  = MATP(3,Mtype)
    AREA = MATP(4,Mtype)
!   ALLOCATE VARIABLES
    ALLOCATE (NDMFS(3,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(2,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
    ALLOCATE (XABSC_WEIG_IN(NITP,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!--------------------------------------------------
        DO I=1,NumEle
           DO J=1,3
              NDMFS(J,I)=0
           ENDDO
           ELEPAR(1,I)=0.
           ELEPAR(2,I)=0.
           DO J=1,NEleAss
              LM(J,I)=0
           ENDDO
        END DO
!--------------------------------------------------
    READ(ADELE) ((NDMFS(I,J),I=1,3),J=1,NumEle),((ELEPAR(J,I),J=1,2),I=1,NumEle),((LM(J,I),J=1,NEleAss),I=1,NumEle)
!   GENERATE
    CALL GAULEG(NITP,XABSC_WEIG_IN(1,1),XABSC_WEIG_IN(1,2))
!
    ALLOCATE (STIFF_IJ(1,NPOR_MAX*NPOR_MAX),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!
!   STIFFNESS MATRIX
    DO I=1,NumEle
!
!   GENERATE LOAD VECTOR
!
!   Neumann boundary
!
       INDSP=NDMFS(3,I)
       IF(NDMFS(2,I).EQ.2) THEN
          N=NDMFS(1,I)
          KI=ID(1,N)
          XX=X(N)
          HI=TRUSS_MFS_HIHJ(NEleAss,I,NDOF(1,N),XX,INDSP)
          DO JJ=1,NDOF(1,N)
             NUM=KI+JJ-1
             NPOR=JJ
             R(NUM)=RNAT(1,N)*HI(NPOR)*AREA
          ENDDO
       ENDIF
!       
!   Dirichlet boundary
!
       IF(NDMFS(2,I).EQ.1) THEN
          N=NDMFS(1,I)
          KI=ID(1,N)
          IF(KI.NE.0) THEN
             XX=X(N)
             DHI=TRUSS_MFS_HIHJ_V2(NEleAss,I,NDOF(1,N)+1,XX,INDSP)
             DO JJ=1,NDOF(1,N)
                NPOR=JJ+1
                NUM=KI+JJ-1
                R(NUM) = RESS(1,N)*DHI(NPOR)*AREA*E
             ENDDO
          ENDIF
       ENDIF
!
!   BODY FORCE
!
       CALL TRUSS_MFS_BODY(I,NEleAss,AREA,INDSP)           
!
!   CALCULATE STIFFNESS MATRIX
!
       DO J=1,NEleAss
          K=LM(J,I)
          IF (K.GE.I) THEN        
             CALL TRUSS_MFS_STIFF(I,K,E,DEN,AREA,NEleAss)
          ENDIF
       ENDDO 
!
    ENDDO              
!
    DEALLOCATE(NDMFS)
    DEALLOCATE(ELEPAR)
    DEALLOCATE(LM)
    DEALLOCATE(XABSC_WEIG_IN)
    DEALLOCATE(STIFF_IJ)
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_MFS_G
!
    SUBROUTINE TRUSS_MFS_STIFF(II,JJ,E,DEN,AREA,NEleAss)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE STIFFNESS MATRIX OF TRUSS_MFS                         .
! .          A(NWK)       :  GLOBAL STIFFNESS MATRIX                  .
! .          AMAXA(NEQ+1) :  ADDRESS MATRIX                           .
! .          N,M          :  NODE NUMBER                              .
! .          AREA1,RAD1,AREA2,RAD2                                    .
! .          E            :  YOUNG'S MODULUS                          .
! .          DEN          :  DENSITY                                  .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II,JJ
    REAL(8)            :: E,DEN,AREA
    INTEGER            :: NEleAss
!   ASSISTANT VARIABLES
    INTEGER            :: N,M
    INTEGER            :: KI,KJ,KK2
!
    INTEGER            :: L,I,J,K,LL
    INTEGER            :: INDK,INDL
    REAL(8)            :: XM,XL,XTMP,WEIGHT
    INTEGER            :: IND1,IND2
    LOGICAL            :: IND1_X,IND2_X,EXPN1,EXPN2
    INTEGER            :: NPOR1,NPOR2
    REAL(8)            :: HIHJ(NPOR_MAX*NPOR_MAX)
    REAL(8)            :: HI(NPOR_MAX),HJ(NPOR_MAX),DHI(NPOR_MAX),DHJ(NPOR_MAX)
    REAL(8)            :: S2,CONST2,CONST1,CONST3,DIRNI,DIRNJ
    REAL(8)            :: STIFF(NPOR_MAX*NPOR_MAX)
    INTEGER            :: NUM
    INTEGER            :: INDSPI,INDSPJ
!
    N       =  NDMFS(1,II)
    M       =  NDMFS(1,JJ)
    INDSPI  =  NDMFS(3,II)
    INDSPJ  =  NDMFS(3,JJ)
!
    IND1 = NDMFS(2,II)
    IND2 = NDMFS(2,JJ)
    IND1_X = (IND1.EQ.1)
    IND2_X = (IND2.EQ.1)
!
    EXPN1=(IND1.NE.0).AND.(IND2.NE.0)
!
    IF (IND1_X) THEN
        NPOR1=NDOF(1,N)+1
    ELSE
        NPOR1=NDOF(1,N)
    ENDIF
!
    IF (IND2_X) THEN
        NPOR2=NDOF(1,M)+1
    ELSE
        NPOR2=NDOF(1,M)
    ENDIF
!
    XL = ELEPAR(2,II)
    XM = X(N)
    DIRNI=ELEPAR(1,II)
!
    L=NPOR1*NPOR2 
    DO I=1,L
       HIHJ(I)=0.0D0
    ENDDO 
!
    CONST2=1.0/ELEPAR(2,JJ)
    DO I=1,NITP
       XTMP=XM+XL*XABSC_WEIG_IN(I,1)
       S2=ABS(XTMP-X(M))*CONST2
       IF(S2.GT.1) CYCLE
       WEIGHT=XABSC_WEIG_IN(I,2)
       IF(EXPN1) THEN
          CONST3=(XTMP-XM)*DIRNI
          EXPN2=CONST3.LT.0.0
          IF(EXPN2) CYCLE
!
          DHI=TRUSS_MFS_HIHJ_V2(NEleAss,II,NPOR1,XTMP,INDSPI)
          DHJ=TRUSS_MFS_HIHJ_V2(NEleAss,JJ,NPOR2,XTMP,INDSPJ)
!
          DO J=1,NPOR1
            DO K=1,NPOR2
               LL=(J-1)*NPOR2+K
               HIHJ(LL)=HIHJ(LL)+DHI(J)*DHJ(K)*WEIGHT
            ENDDO
          ENDDO             
       ELSE
!
          DHI=TRUSS_MFS_HIHJ_V2(NEleAss,II,NPOR1,XTMP,INDSPI)
          DHJ=TRUSS_MFS_HIHJ_V2(NEleAss,JJ,NPOR2,XTMP,INDSPJ)
!
          DO J=1,NPOR1
            DO K=1,NPOR2
               LL=(J-1)*NPOR2+K
               HIHJ(LL)=HIHJ(LL)+DHI(J)*DHJ(K)*WEIGHT
            ENDDO
          ENDDO    
       ENDIF
    ENDDO
!
    CONST1=E*AREA*XL
    DO I=1,L
       STIFF_IJ(1,I)=HIHJ(I)*CONST1
       STIFF(I)=HIHJ(I)*CONST1
    ENDDO                     
!
    IF(IND1.EQ.1) THEN
       IF(IND2.EQ.1) THEN
          XTMP = X(N)
          HI = TRUSS_MFS_HIHJ(NEleAss,II,NPOR1,XTMP,INDSPI)
          HJ = TRUSS_MFS_HIHJ(NEleAss,JJ,NPOR2,XTMP,INDSPJ)
          DHI = TRUSS_MFS_HIHJ_V2(NEleAss,II,NPOR1,XTMP,INDSPI)
          DHJ = TRUSS_MFS_HIHJ_V2(NEleAss,JJ,NPOR2,XTMP,INDSPJ)
          DO I=1,NPOR1
             DO J=1,NPOR2
                K=(I-1)*NPOR2+J
                STIFF_IJ(1,K)=STIFF_IJ(1,K)+DHI(I)*HJ(J)+HI(I)*DHJ(J)
             ENDDO
          ENDDO
       ENDIF
    ENDIF
!
!   UPDATE LOAD VECTORS
!
    IF(IND1_X) THEN
        IF(ID(1,M).NE.0) THEN
           DO I=1,NDOF(1,M)
              IF(IND2_X) THEN
                 K=I+1
              ELSE
                 K=I
              ENDIF
              LL=K
              NUM=ID(1,M)+I-1
              R(NUM)=R(NUM)-STIFF_IJ(1,LL)*RESS(1,N)
           ENDDO           
        ENDIF
    ENDIF
!
    IF(IND2_X.AND.(N.NE.M)) THEN
       DO I=1,NDOF(1,N)
          LL=(I-1)*NPOR2+1
          NUM=ID(1,N)+I-1
          R(NUM)=R(NUM)-STIFF_IJ(1,LL)*RESS(1,M)
       ENDDO
    ENDIF 
!
!   ASSEMBLY
!
    K=1
    INDK=1
    IF(IND1_X) THEN
       K=2
       INDK=2
    ENDIF
    DO WHILE(K.LE.NPOR1)
       L=1
       INDL=1
       IF(IND2_X) THEN
          L=2
          INDL=2
       ENDIF
       DO WHILE(L.LE.NPOR2)
          KI=ID(1,N)+K-INDK
          KJ=ID(1,M)+L-INDL
          KK2=(K-1)*NPOR2+L
          IF((N.NE.M).OR.(KJ.GE.KI)) THEN
             IF(KJ.LE.KI) THEN
                LL=MAXA_COLSOL(KI)+ABS(KJ-KI)
             ELSE
                LL=MAXA_COLSOL(KJ)+ABS(KJ-KI)
             ENDIF
             K_COLSOL(LL)=K_COLSOL(LL)+STIFF_IJ(1,KK2)
             KMATRIX(LL)=KMATRIX(LL)+STIFF(KK2)
          ENDIF
         L=L+1
       ENDDO
       K=K+1
    ENDDO          
! 
    RETURN
!}}}
    END SUBROUTINE TRUSS_MFS_STIFF
!
    SUBROUTINE TRUSS_MFS_BODY(II,NEleAss,AREA,INDSP)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE  BODY FORCE CONTRIBUTION TO LOAD VECTOR               .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER :: II
    INTEGER :: NEleAss
    REAL(8) :: AREA
    INTEGER :: INDSP
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: N,NUM,NPOR,NPORMAX
    REAL(8) :: HI(NPOR_MAX),STIFF(NPOR_MAX)
    REAL(8) :: XM,XL,XTMP,WEIGHT
    INTEGER :: IND1
    INTEGER :: INC
    LOGICAL :: EXPN2,EXPN3
    REAL(8) :: CONST1,CONST2,DIRN
!
       N   = NDMFS(1,II)
       IND1= NDMFS(2,II)
       IF(IND1.EQ.1) THEN
          NPORMAX=NDOF(1,N)+1
          INC=1
       ELSE
          NPORMAX=NDOF(1,N)
          INC=0
       ENDIF
! 
       DO I=1,NPOR_MAX
          STIFF(I)=0.0D0
       ENDDO
!
       XM = X(N)
       XL = ELEPAR(2,II)
       NUM=ID(1,N)
       DIRN=ELEPAR(1,II)
       IF(NUM.NE.0) THEN
          DO J=1,NITP
             XTMP = XM+XL*XABSC_WEIG_IN(J,1)
             WEIGHT = XABSC_WEIG_IN(J,2)
!
             EXPN2=(IND1.NE.0)
             IF(EXPN2) THEN
                CONST2=(XTMP-XM)*DIRN
                EXPN3=CONST2.LT.0.0
                IF(EXPN3) CYCLE
                   HI=TRUSS_MFS_HIHJ(NEleAss,II,NPORMAX,XTMP,INDSP)
                   DO K=1,NDOF(1,N)
                      NPOR=K+INC
                      STIFF(K)=STIFF(K)+HI(NPOR)*WEIGHT
                   ENDDO
             ELSE
                   HI=TRUSS_MFS_HIHJ(NEleAss,II,NPORMAX,XTMP,INDSP)
                   DO K=1,NDOF(1,N)
                      NPOR=K+INC
                      STIFF(K)=STIFF(K)+HI(NPOR)*WEIGHT
                   ENDDO
             ENDIF
!               
          ENDDO
       ENDIF
!
       CONST1=AREA*XL
       DO I=1,NDOF(1,N)
          NUM   = ID(1,N)+I-1
          R(NUM)= R(NUM) + STIFF(I)*CONST1
       ENDDO                                 
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_MFS_BODY
!
!}}}
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *   
    SUBROUTINE QUADS_MFS_G(ADELE,ADBC)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      MASS      MATRIX
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE,ADBC
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L
    INTEGER :: NEleAss,NumEle,Mtype,ITYPE
    REAL(8) :: E,PR,DEN,THIC
    INTEGER :: BcType
!   ASSISTANT VARIABLES
    INTEGER :: KI
    INTEGER :: ALLocateStatus
    REAL(8) :: TT1,TT2
!   READ 
    NEleAss= NPAR(5)
    NumEle = NPAR(4)
    Mtype  = NPAR(3)
    ITYPE  = NPAR(2)
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    PR   = MATP(2,Mtype)
    DEN  = MATP(3,Mtype)
    THIC = MATP(4,Mtype)
!   ALLOCATE VARIABLES
    ALLOCATE (NDMFS(6,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(4,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!
    ALLOCATE (XABSC_WEIG_BC(NITP_BC,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
    ALLOCATE (XABSC_WEIG_IN(NITP_IN,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!--------------------------------------------------
!        DO I=1,NumEle
!           DO J=1,5
!              NDMFS(J,I)=0
!           ENDDO
!!
!           DO J=1,3
!              ELEPAR(J,I)=0.
!           ENDDO
!!
!           DO J=1,NEleAss
!              LM(J,I)=0
!           ENDDO
!        END DO
!!
!        DO I=1,NITP_BC
!           XABSC_WEIG_BC(I,1)=0
!           XABSC_WEIG_BC(I,2)=0
!        ENDDO
!!
!        DO I=1,NITP_IN
!           XABSC_WEIG_IN(I,1)=0
!           XABSC_WEIG_IN(I,2)=0
!        ENDDO
!--------------------------------------------------
    READ(ADELE) ((NDMFS(I,J),I=1,6),J=1,NumEle),((ELEPAR(J,I),J=1,4),I=1,NumEle),((LM(J,I),J=1,NEleAss),I=1,NumEle)
!   GENERATE CONSTITUTIVE MATRIX
    CALL G_KCS(2,ITYPE,E,PR)
!   GENERATE LOAD VECTOR  
    CALL GAULEG(NITP_BC,XABSC_WEIG_BC(1,1),XABSC_WEIG_BC(1,2))
!   CALL CPU_TIME(TT1)
    CALL GAULEG(NITP_IN,XABSC_WEIG_IN(1,1),XABSC_WEIG_IN(1,2)) 
!   CALL CPU_TIME(TT2)
!   WRITE(*,*) TT2-TT1
!
    READ(ADBC) (NPAR(J),J=1,7)
!
    ALLOCATE (STIFF_IJ(NPOR_MAX*NPOR_MAX,4),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!
    BcType=NPAR(1)
    SELECT CASE(BcType)
    CASE(1)
        CALL MFS_LINE_G(ADBC,NEleAss,THIC)
    CASE DEFAULT
        WRITE(*,*), "ERROR, WRONG INPUT IN GENERATE BOUNDARY PROCEDURE"
        STOP
    END SELECT
!   GENERATE STIFFNESS MATRIX 
    WRITE(*,*) "Generating Stiff matrix OF ELEMENT"    
    DO I=1,NumEle
!   CALCULATE STIFFNESS MATRIX
       DO J=1,NEleAss
          K=LM(J,I)
          IF (K.GE.I) THEN
!   CALCULATE ENTRIES OF STIFFNESS MATRIX  
             CALL QUADS_MFS_STIFF_V2(I,K,NEleAss,THIC)
          ENDIF
       ENDDO              
    ENDDO
!
    DEALLOCATE(NDMFS)
    DEALLOCATE(ELEPAR)
    DEALLOCATE(LM)
    DEALLOCATE(XABSC_WEIG_BC)
    DEALLOCATE(XABSC_WEIG_IN)
    DEALLOCATE(STIFF_IJ)
    DEALLOCATE(KCS)
!    
    RETURN
!}}}
    END SUBROUTINE QUADS_MFS_G   
!
    SUBROUTINE MFS_LINE_G(ADBC,NEleAss,THIC)
!{{{
!--------------------------------------------------
!   GENERATE LOAD VECTOR
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADBC,NEleAss
    REAL(8),INTENT(INOUT) :: THIC
!   ASSISTANT VARIABLES
    INTEGER :: INDSP
    INTEGER :: I,J,K,JJ,KK,JJJ
    INTEGER :: NumBc
    INTEGER :: ALLocateStatus
!  
    INTEGER :: L1,L2,N1,N2,L,N
    REAL(8) :: XM,XL,YM,YL,XX,YY
    REAL(8) :: FLOAD(2),FUNC(2)
    INTEGER :: INDBC,IND
    LOGICAL :: INDBCX,INDBCY
    INTEGER :: NUM
    INTEGER :: NPOR,NPORMAX 
    REAL(8) :: DIV
    REAL(8) :: S(NPOR_MAX),HI(NPOR_MAX),LENGTH
    REAL(8) :: S1,S2,S3,DHIDX,DHIDY
    REAL(8) :: NX,NY
!
!   READ 
    NumBc  = NPAR(2)
!   ALLOCATE VARIABLES
    ALLOCATE (BCED(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (BCIND(4,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (BCPAR(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!
!    DO I=1,NumBc
!       DO J=1,4
!          BCIND(J,I)=0
!       ENDDO
!       DO J=1,2
!          BCED(J,I)=0
!          BCPAR(J,I)=0.
!       ENDDO
!    ENDDO
!
    READ(ADBC) ((BCED(J,I),J=1,2),I=1,NumBc),((BCIND(J,I),J=1,4),I=1,NumBc),((BCPAR(J,I),J=1,2),I=1,NumBc)
!
    WRITE(*,*) "Generating Load vector OF BOUNDARY "
!
    DO I=1,NumBc
!      CALCULATE LOAD VECTOR
       L1=BCED(1,I)
       L2=BCED(2,I)
       N1=NDMFS(1,L1)
       N2=NDMFS(1,L2)
!
       NDMFS(5,L1)=L2
       NDMFS(4,L2)=L1
!
       LENGTH=(X(N1)-X(N2))*(X(N1)-X(N2))+(Y(N1)-Y(N2))*(Y(N1)-Y(N2))
       LENGTH=SQRT(LENGTH)
!
       XM=0.5*(X(N1)+X(N2))
       XL=0.5*(X(N1)-X(N2))
       YM=0.5*(Y(N1)+Y(N2))
       YL=0.5*(Y(N1)-Y(N2))
!
       FLOAD(1)=BCPAR(1,I)
       FLOAD(2)=BCPAR(2,I)
       FUNC(1)=BCIND(3,I)
       FUNC(2)=BCIND(4,I)
!
       NX=(Y(N2)-Y(N1))/LENGTH
       NY=-(X(N2)-X(N1))/LENGTH
!
!       DO J=1,2
!          IF(BCIND(J,I).EQ.1) THEN
!             DISP(J)=FLOAD(J)
!          ELSE
!             DISP(J)=0.
!          ENDIF
!       ENDDO
!------------------------------------------
       DO K=1,2
          L=BCED(K,I)
          N=NDMFS(1,L)
          IF(NDMFS(2,L).EQ.1) THEN
            NPORMAX=NDOF(1,N)+1
          ELSE
            NPORMAX=NDOF(1,N)
          ENDIF
!
          DO JJ=1,NPORMAX
             S(JJ)=0.0D0
          ENDDO
!
          DO JJ=1,NITP_BC
             XX = XM+XL*XABSC_WEIG_BC(JJ,1)
             YY = YM+YL*XABSC_WEIG_BC(JJ,1)
             DIV= XABSC_WEIG_BC(JJ,2)
!
             HI = QUADS_MFS_HIHJ(NEleAss,L,NPORMAX,XX,YY)
             DO KK=1,NPORMAX
                S(KK)= S(KK) + HI(KK)*DIV
             ENDDO
          ENDDO
!
          DO J=1,2
             INDBC=BCIND(J,I)
             SELECT CASE(INDBC)
             CASE(2)
                 IND=NDMFS(J+1,L)
                 DO JJ=1,NDOF(J,N)
                    NUM=ID(J,N)+JJ-1
                    IF(IND.EQ.1) THEN
                       KK=JJ+1
                    ELSE
                       KK=JJ
                    ENDIF
                    R(NUM)=R(NUM)+S(KK)*LENGTH/2.*THIC*FLOAD(J)
                 ENDDO
             CASE(1) 
                 CONTINUE
             CASE DEFAULT
                 WRITE(*,*) "ERROR, IND.GT.1 IN SUBROUTINE MFS_LINE_G!"
                 STOP
             END SELECT
          ENDDO       
       ENDDO
!
!!             DO K=1,2
!                L=BCED(K,I)
!                N=NDMFS(1,L)   
!                DO JJ=1,NDOF(J,N)
!                   NUM=ID(J,N)+JJ-1
!!
!                   NPOR=JJ+1
!                   S1=0.
!                   S2=0.
!                   S3=0.
!                   DO KK=1,NITP_BC
!                      XX = XM+XL*XABSC_WEIG_BC(KK,1)
!                      YY = YM+YL*XABSC_WEIG_BC(KK,1)
!                      DIV= XABSC_WEIG_BC(KK,2)
!!
!                      DHIDX = QUADS_MFS_HIHJ(NEleAss,L,NPOR,XX,YY,2) 
!                      DHIDY = QUADS_MFS_HIHJ(NEleAss,L,NPOR,XX,YY,3)
!!
!                      IF(J.EQ.1) THEN
!                        S1=S1+DHIDX*NX*DISP(1)*DIV
!                        S2=S2+DHIDX*NY*DISP(2)*DIV
!                        S3=S3+(DHIDY*NY*DISP(1)+DHIDY*NX*DISP(2))*DIV
!                      ELSE
!                        S2=S2+DHIDY*NX*DISP(1)*DIV
!                        S1=S1+DHIDY*NY*DISP(2)*DIV
!                        S3=S3+(DHIDX*NY*DISP(1)+DHIDX*NX*DISP(2))*DIV
!                      ENDIF
!                   ENDDO
!                   GR(NUM)=GR(NUM)-(KCS(1)*S1+KCS(2)*S2+KCS(3)*S3)*LENGTH/2.*THIC
!                ENDDO
!             ENDDO     
!******************************
! IMPOSITION OF DIRICHLET Bcs
!******************************
        INDBCX=BCIND(1,I).EQ.1
        INDBCY=BCIND(2,I).EQ.1
        IF(INDBCX.OR.INDBCY) THEN
          DO J=1,2
             DO K=J,2
                L1=BCED(J,I)
                L2=BCED(K,I)
!******************************
!  METHOD OF WEAK FORM
!******************************
!                CALL MFS_LINE_STIFF_WEAK(L1,L2,I,NEleAss,NX,NY,XM,XL,YM,YL,LENGTH,THIC)
!******************************
!  PENALTY METHOD
!******************************
                 CALL MFS_LINE_STIFF_PENALTY(L1,L2,I,NEleAss,NX,NY,XM,XL,YM,YL,LENGTH,THIC)
             ENDDO
          ENDDO
        ENDIF
!
    ENDDO
!
    DEALLOCATE(BCED)
    DEALLOCATE(BCIND)
    DEALLOCATE(BCPAR) 
!
    RETURN
!}}}
    END SUBROUTINE MFS_LINE_G
!
    SUBROUTINE MFS_LINE_STIFF_WEAK(II,JJ,KK,NEleAss,NX,NY,XM,XL,YM,YL,LENGTH,THIC)
!{{{
!--------------------------------------------------
!   GENERATE STIFFNESS MATRIX RESULT FROM BCs
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER :: II,JJ,KK,NEleAss
    REAL(8) :: NX,NY,XM,XL,YM,YL,LENGTH
    REAL(8) :: THIC
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L
    INTEGER :: LL
    INTEGER :: N1,N2
!   REAL(8) :: S(2,2)
    INTEGER :: NPOR1,NPOR2
!
    LOGICAL :: IND1_X,IND1_Y,IND2_X,IND2_Y 
    REAL(8) :: HIHJ(NPOR_MAX*NPOR_MAX,4)
    REAL(8) :: HI(3,NPOR_MAX),HJ(3,NPOR_MAX)
    REAL(8) :: XTMP,YTMP
    REAL(8) :: WEIGHT
    REAL(8) :: CONST1
    REAL(8) :: KCS1,KCS2,KCS3
!
!    INTEGER :: KI,KJ
!
    N1=NDMFS(1,II)
    N2=NDMFS(1,JJ)
!
    IND1_X  =  NDMFS(2,II).EQ.1
    IND1_Y  =  NDMFS(3,II).EQ.1
    IND2_X  =  NDMFS(2,JJ).EQ.1
    IND2_Y  =  NDMFS(3,JJ).EQ.1
!
    IF (IND1_X) THEN
        NPOR1=NDOF(1,N1)+1
    ELSE
        NPOR1=NDOF(1,N1)
    ENDIF
!
    IF (IND2_X) THEN
        NPOR2=NDOF(1,N2)+1
    ELSE
        NPOR2=NDOF(1,N2)
    ENDIF
!
!    DO I=1,NPOR_MAX*NPOR_MAX
!       STIFF_IJ(1,I)=0.0D0
!       STIFF_IJ(2,I)=0.0D0
!       STIFF_IJ(3,I)=0.0D0
!       STIFF_IJ(4,I)=0.0D0
!    ENDDO
!
    L=NPOR1*NPOR2
    DO I=1,L
       HIHJ(I,1)=0.0D0
       HIHJ(I,2)=0.0D0
       HIHJ(I,3)=0.0D0
       HIHJ(I,4)=0.0D0
    ENDDO
!
!    DO I=1,3
!       DO J=1,NPOR_MAX
!          HI(J,I)=0.0D0
!          HJ(J,I)=0.0D0
!       ENDDO
!    ENDDO
!
      DO I=1,NITP_IN
         XTMP = XM+XL*XABSC_WEIG_IN(I,1)
         YTMP = YM+YL*XABSC_WEIG_IN(I,1)
!
         CALL QUADS_MFS_HIHJ_V2(NEleAss,II,JJ,XTMP,YTMP,NPOR1,NPOR2,HI(1,1),HJ(1,1))
!
         WEIGHT=XABSC_WEIG_IN(I,2)
         DO J=1,NPOR1
            DO K=1,NPOR2
               LL=(J-1)*NPOR2+K
               HIHJ(LL,1)=HIHJ(LL,1)+HI(1,J)*HJ(2,K)*WEIGHT
               HIHJ(LL,2)=HIHJ(LL,2)+HI(1,J)*HJ(3,K)*WEIGHT
               HIHJ(LL,3)=HIHJ(LL,3)+HI(2,J)*HJ(1,K)*WEIGHT
               HIHJ(LL,4)=HIHJ(LL,4)+HI(3,J)*HJ(1,K)*WEIGHT
            ENDDO
         ENDDO
      ENDDO
!
      CONST1=LENGTH/2.0*THIC
      KCS1=KCS(1)
      KCS2=KCS(2)
      KCS3=KCS(3)
      DO I=1,L
         STIFF_IJ(I,1)=(HIHJ(I,1)+HIHJ(I,3))*NX*KCS1+(HIHJ(I,2)+HIHJ(I,4))*NY*KCS3
         STIFF_IJ(I,2)=(HIHJ(I,2)*NX+HIHJ(I,3)*NY)*KCS2+(HIHJ(I,1)*NY+HIHJ(I,4)*NX)*KCS3
         STIFF_IJ(I,3)=(HIHJ(I,1)*NY+HIHJ(I,4)*NX)*KCS2+(HIHJ(I,2)*NX+HIHJ(I,3)*NY)*KCS3
         STIFF_IJ(I,4)=(HIHJ(I,2)+HIHJ(I,4))*NY*KCS1+(HIHJ(I,1)+HIHJ(I,3))*NX*KCS3
!
         STIFF_IJ(I,1)=STIFF_IJ(I,1)*CONST1
         STIFF_IJ(I,2)=STIFF_IJ(I,2)*CONST1
         STIFF_IJ(I,3)=STIFF_IJ(I,3)*CONST1
         STIFF_IJ(I,4)=STIFF_IJ(I,4)*CONST1
      ENDDO  
!
      CALL MFS_LINE_ADDBAN(II,JJ,KK,NPOR1,NPOR2,IND1_X,IND1_Y,IND2_X,IND2_Y,STIFF_IJ(1,1))
!
    RETURN
!}}}
    END SUBROUTINE MFS_LINE_STIFF_WEAK
!
    SUBROUTINE MFS_LINE_STIFF_PENALTY(II,JJ,KK,NEleAss,NX,NY,XM,XL,YM,YL,LENGTH,THIC)
!{{{
!--------------------------------------------------
!   GENERATE STIFFNESS MATRIX RESULT FROM BCs
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER :: II,JJ,KK,NEleAss
    REAL(8) :: NX,NY,XM,XL,YM,YL,LENGTH
    REAL(8) :: THIC
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L
    INTEGER :: LL
    INTEGER :: N1,N2
!   REAL(8) :: S(2,2)
    INTEGER :: NPOR1,NPOR2
!
    LOGICAL :: IND1_X,IND1_Y,IND2_X,IND2_Y 
    REAL(8) :: HIHJ(NPOR_MAX*NPOR_MAX)
    REAL(8) :: HI(NPOR_MAX),HJ(NPOR_MAX)
    REAL(8) :: XTMP,YTMP
    REAL(8) :: WEIGHT
    REAL(8) :: CONST1
    REAL(8) :: KCS1,KCS2,KCS3
!
!    INTEGER :: KI,KJ
!
    N1=NDMFS(1,II)
    N2=NDMFS(1,JJ)
!
    IND1_X  =  NDMFS(2,II).EQ.1
    IND1_Y  =  NDMFS(3,II).EQ.1
    IND2_X  =  NDMFS(2,JJ).EQ.1
    IND2_Y  =  NDMFS(3,JJ).EQ.1
!
    IF (IND1_X) THEN
        NPOR1=NDOF(1,N1)+1
    ELSE
        NPOR1=NDOF(1,N1)
    ENDIF
!
    IF (IND2_X) THEN
        NPOR2=NDOF(1,N2)+1
    ELSE
        NPOR2=NDOF(1,N2)
    ENDIF
!
    L=NPOR1*NPOR2
    DO I=1,L
       HIHJ(I)=0.0D0
    ENDDO
!
      DO I=1,NITP_IN
         XTMP = XM+XL*XABSC_WEIG_IN(I,1)
         YTMP = YM+YL*XABSC_WEIG_IN(I,1)
         WEIGHT=XABSC_WEIG_IN(I,2)
!
         HI = QUADS_MFS_HIHJ(NEleAss,II,NPOR1,XTMP,YTMP)
         HJ = QUADS_MFS_HIHJ(NEleAss,JJ,NPOR2,XTMP,YTMP)
         DO J=1,NPOR1
            DO K=1,NPOR2
               LL=(J-1)*NPOR2+K
               HIHJ(LL)=HIHJ(LL)+HI(J)*HJ(K)*WEIGHT
            ENDDO
         ENDDO
      ENDDO
!
      CONST1=LENGTH/2.0*THIC*BETA
      DO I=1,L
         STIFF_IJ(I,1)=HIHJ(I)*CONST1
         STIFF_IJ(I,2)=0.0D0
         STIFF_IJ(I,3)=0.0D0
         STIFF_IJ(I,4)=HIHJ(I)*CONST1
      ENDDO  
!
      CALL MFS_LINE_ADDBAN(II,JJ,KK,NPOR1,NPOR2,IND1_X,IND1_Y,IND2_X,IND2_Y,STIFF_IJ(1,1))
!
    RETURN
!}}}
    END SUBROUTINE MFS_LINE_STIFF_PENALTY
!
    SUBROUTINE MFS_LINE_ADDBAN(II,JJ,KK,NPOR1,NPOR2,IND1_X,IND1_Y,IND2_X,IND2_Y,S)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II,JJ,KK,NPOR1,NPOR2
    LOGICAL            :: IND1_X,IND1_Y,IND2_X,IND2_Y
    REAL(8)            :: S(NPOR_MAX*NPOR_MAX,4)
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L
    INTEGER :: JMIN,JMAX
    INTEGER :: N,M
    LOGICAL :: INDBCX,INDBCY
    INTEGER :: INDK,INDL,KK1,KK2
    INTEGER :: KI,KJ
    INTEGER :: LL
!
    INDBCX=(BCIND(1,KK).EQ.2)
    INDBCY=(BCIND(2,KK).EQ.2)
!
    N         =  NDMFS(1,II)
    M         =  NDMFS(1,JJ)
!
    DO I=1,2
       IF(INDBCX) THEN
          JMIN=2
       ELSE
          JMIN=1
       ENDIF
       IF(INDBCY) THEN
          JMAX=1
       ELSE
          JMAX=2
       ENDIF
!
       DO J=JMIN,JMAX
          K=1
          INDK=1
          IF((IND1_X).AND.(I.EQ.1)) THEN
             K=2
             INDK=2
          ENDIF
          IF((IND1_Y).AND.(I.EQ.2)) THEN
             K=2
             INDK=2
          ENDIF          
          DO WHILE(K.LE.NPOR1)
             L=1
             INDL=1
             IF((IND2_X).AND.(J.EQ.1)) THEN
                 L=2
                 INDL=2
             ENDIF
             IF((IND2_Y).AND.(J.EQ.2)) THEN
                 L=2
                 INDL=2
             ENDIF 
             DO WHILE(L.LE.NPOR2)
                KI=ID(I,N)+K-INDK
                KJ=ID(J,M)+L-INDL
                KK1=(I-1)*2+J
                KK2=(K-1)*NPOR2+L
                IF ((N.NE.M).OR.(KJ.GE.KI)) THEN
                   IF(KJ.LE.KI) THEN
                      LL=MAXA_COLSOL(KI)+ABS(KJ-KI)
                   ELSE
                      LL=MAXA_COLSOL(KJ)+ABS(KJ-KI)
                   ENDIF
                   K_COLSOL(LL)=K_COLSOL(LL)-S(KK2,KK1)
                ENDIF
                L=L+1
             ENDDO
             K=K+1
          ENDDO
       ENDDO
    ENDDO
!
    RETURN    
!}}}
    END SUBROUTINE MFS_LINE_ADDBAN
!
    SUBROUTINE QUADS_MFS_STIFF_V2(II,JJ,NEleAss,THIC)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE STIFFNESS MATRIX OF TRUSS_MFS                         .
! .          A(NWK)       :  GLOBAL STIFFNESS MATRIX                  .
! .          AMAXA(NEQ+1) :  ADDRESS MATRIX                           .
! .          N,M          :  NODE NUMBER                              .
! .          E            :  YOUNG'S MODULUS                          .
! .          PR           :  POisson ratio                            .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II,JJ
    INTEGER            :: NEleAss
    REAL(8)            :: THIC
!   ASSISTANT VARIABLES
    INTEGER            :: I,J,K,L
    INTEGER            :: N,M
    INTEGER            :: NPOR1,NPOR2
    LOGICAL            :: IND1_X,IND1_Y,IND2_X,IND2_Y
    LOGICAL            :: EXPN1
    INTEGER            :: IND1,IND2
!
    REAL(8)            :: XM,YM,XL,YL,LENGTHN,LENGTHM,XTMP,WEIGHT
    REAL(8)            :: HIHJH(NPOR_MAX*NPOR_MAX,4),HIHJG(NPOR_MAX*NPOR_MAX,4)
    REAL(8)            :: KCS1,KCS2,KCS3,CONST1
!
    N       =  NDMFS(1,II)
    M       =  NDMFS(1,JJ)
!
    IND1    = NDMFS(2,II) 
    IND2    = NDMFS(2,JJ)
!
    IND1_X  = (IND1.EQ.1)
    IND1_Y  = (NDMFS(3,II).EQ.1)
    IND2_X  = (IND2.EQ.1)
    IND2_Y  = (NDMFS(3,JJ).EQ.1)
!
    EXPN1=((IND1.NE.0).AND.(IND2.NE.0))
!
    IF (IND1_X) THEN
        NPOR1=NDOF(1,N)+1
    ELSE
        NPOR1=NDOF(1,N)
    ENDIF
!
    IF (IND2_X) THEN
        NPOR2=NDOF(1,M)+1
    ELSE
        NPOR2=NDOF(1,M)
    ENDIF
!
    XM=X(N)
    YM=Y(N)
    XL=ELEPAR(3,II)
    YL=ELEPAR(4,II)
!
!    DO I=1,NPOR_MAX*NPOR_MAX
!       STIFF_IJ(1,I)=0.0D0
!       STIFF_IJ(2,I)=0.0D0
!       STIFF_IJ(3,I)=0.0D0
!       STIFF_IJ(4,I)=0.0D0
!    ENDDO
!
    L=NPOR1*NPOR2
    DO I=1,L
       HIHJH(I,1)=0.0D0
       HIHJH(I,2)=0.0D0
       HIHJH(I,3)=0.0D0
       HIHJH(I,4)=0.0D0
!
       HIHJG(I,1)=0.0D0
       HIHJG(I,2)=0.0D0
       HIHJG(I,3)=0.0D0
       HIHJG(I,4)=0.0D0
    ENDDO
!
    DO I=1,NITP_IN
       XTMP=XM+XL*XABSC_WEIG_IN(I,1)
       WEIGHT=XABSC_WEIG_IN(I,2)
!
       CALL QUADS_MFS_GAUSSG_V2(N,M,II,JJ,NPOR1,NPOR2,NEleAss,XTMP,YM,YL,HIHJG(1,1),EXPN1,L)
       DO J=1,L
          HIHJH(J,1)=HIHJH(J,1)+HIHJG(J,1)*WEIGHT
          HIHJH(J,2)=HIHJH(J,2)+HIHJG(J,2)*WEIGHT
          HIHJH(J,3)=HIHJH(J,3)+HIHJG(J,3)*WEIGHT
          HIHJH(J,4)=HIHJH(J,4)+HIHJG(J,4)*WEIGHT
       ENDDO
    ENDDO
!
    KCS1=KCS(1)
    KCS2=KCS(2)
    KCS3=KCS(3)
    CONST1=THIC*XL
    DO I=1,L
       STIFF_IJ(I,1)=(KCS1*HIHJH(I,1)+KCS3*HIHJH(I,4))*CONST1
       STIFF_IJ(I,2)=(KCS2*HIHJH(I,2)+KCS3*HIHJH(I,3))*CONST1
       STIFF_IJ(I,3)=(KCS2*HIHJH(I,3)+KCS3*HIHJH(I,2))*CONST1
       STIFF_IJ(I,4)=(KCS1*HIHJH(I,4)+KCS3*HIHJH(I,1))*CONST1
!      WRITE(*,*) STIFF_IJ(1,K),STIFF_IJ(2,K),STIFF_IJ(2,K),STIFF_IJ(3,K)
    ENDDO 
!
    CALL QUADS_MFS_ADDBAN_V2(II,JJ,NPOR1,NPOR2,IND1_X,IND1_Y,IND2_X,IND2_Y,STIFF_IJ(1,1))
!   
    RETURN
!!}}}
    END SUBROUTINE QUADS_MFS_STIFF_V2
!
    SUBROUTINE QUADS_MFS_GAUSSG_V2(N,M,II,JJ,NPOR1,NPOR2,NEleAss,XTMP,YM,YL,HIHJG,EXPN1,LL)
!{{{
    IMPLICIT NONE
    INTEGER :: N,M,II,JJ,NPOR1,NPOR2,NEleAss
    REAL(8) :: XTMP,YM,YL
    REAL(8) :: HIHJG(NPOR_MAX*NPOR_MAX,4)
    LOGICAL :: EXPN1
    INTEGER :: LL
!   ASSISTANT VARIABLES
    INTEGER :: INDSPI,INDSPJ
    INTEGER :: I,J,K,L
    INTEGER :: L1,L2,N1,N2
    REAL(8) :: YTMP
    REAL(8) :: HI(3,NPOR1),HJ(3,NPOR2)
    REAL(8) :: X1,X2,Y1,Y2,XX1,YY1,WEIGHT
    LOGICAL :: REST1,REST2,REST3,LOGIC1,LOGIC2
    LOGICAL :: EXPRN1
    REAL(8) :: RAD1,RAD2,SX,SY,S
!
    DO I=1,4
       DO J=1,LL
         HIHJG(J,I)=0.0D0
       ENDDO
    ENDDO
!
    INDSPI=NDMFS(6,II)
    INDSPJ=NDMFS(6,JJ)
!
    XX1=XTMP-X(N)
    S=(XTMP-X(M))*(XTMP-X(M))
    DO I=1,NITP_IN     
       YTMP=YM+YL*XABSC_WEIG_IN(I,1)
       YY1=YTMP-Y(N)
!
       RAD1=ELEPAR(3,JJ)
       RAD1=1.0/RAD1
       IF(INDSPJ.EQ.3) THEN
          RAD2=ELEPAR(4,JJ)
          RAD2=1.0/RAD2
!
          SX=1.0-ABS(XTMP-X(M))*RAD1
          EXPRN1=SX.LT.1E-8
          IF(EXPRN1) CYCLE
          SY=1.0-ABS(YTMP-Y(M))*RAD2
          EXPRN1=SY.LT.1E-8
          IF(EXPRN1) CYCLE               
       ELSE
           S=S+(YTMP-Y(M))*(YTMP-Y(M))
           S=S*RAD1*RAD1
           EXPRN1=S.GT.1.0
           IF(EXPRN1) CYCLE
       ENDIF
!
       RAD1=ELEPAR(3,II)
       RAD1=1.0/RAD1
       IF(INDSPI.EQ.3) THEN
          RAD2=ELEPAR(4,II)
          RAD2=1.0/RAD2
!
          SX=1.0-ABS(XTMP-X(N))*RAD1
          EXPRN1=SX.LT.1E-8
          IF(EXPRN1) CYCLE
          SY=1.0-ABS(YTMP-Y(N))*RAD2
          EXPRN1=SY.LT.1E-8
          IF(EXPRN1) CYCLE               
       ELSE
           S=XX1*XX1+YY1*YY1
           S=S*RAD1*RAD1
           EXPRN1=S.GT.1.0
           IF(EXPRN1) CYCLE
       ENDIF
!
       WEIGHT=XABSC_WEIG_IN(I,2)
       IF (EXPN1) THEN
          L1=NDMFS(4,II)
          L2=NDMFS(5,II)
          N1=NDMFS(1,L1)
          N2=NDMFS(1,L2)
! 
          X2=X(N2)-X(N)
          X1=X(N1)-X(N)
          Y2=Y(N2)-Y(N)
          Y1=Y(N1)-Y(N)
!
          REST1=areClockwise(X1,X2,Y1,Y2)
          REST2=areClockwise(XX1,X2,YY1,Y2)
          REST3=areClockwise(XX1,X1,YY1,Y1)
          LOGIC1=(REST1.AND.(.NOT.(REST2.AND.(.NOT.REST3))))
          LOGIC2=((.NOT.REST2).AND.(REST3.AND.(.NOT.REST2)))
          IF ((LOGIC1).OR.(LOGIC2)) CYCLE
!
          CALL QUADS_MFS_HIHJ_V2(NEleAss,II,JJ,XTMP,YTMP,NPOR1,NPOR2,HI(1,1),HJ(1,1))
!
          DO J=1,NPOR1
             DO K=1,NPOR2
                L=(J-1)*NPOR2+K
                HIHJG(L,1)=HIHJG(L,1)+HI(2,J)*HJ(2,K)*WEIGHT
                HIHJG(L,2)=HIHJG(L,2)+HI(2,J)*HJ(3,K)*WEIGHT
                HIHJG(L,3)=HIHJG(L,3)+HI(3,J)*HJ(2,K)*WEIGHT
                HIHJG(L,4)=HIHJG(L,4)+HI(3,J)*HJ(3,K)*WEIGHT
             ENDDO             
          ENDDO                
       ELSE
          CALL QUADS_MFS_HIHJ_V2(NEleAss,II,JJ,XTMP,YTMP,NPOR1,NPOR2,HI(1,1),HJ(1,1))
!
          DO J=1,NPOR1
             DO K=1,NPOR2
                L=(J-1)*NPOR2+K
                HIHJG(L,1)=HIHJG(L,1)+HI(2,J)*HJ(2,K)*WEIGHT
                HIHJG(L,2)=HIHJG(L,2)+HI(2,J)*HJ(3,K)*WEIGHT
                HIHJG(L,3)=HIHJG(L,3)+HI(3,J)*HJ(2,K)*WEIGHT
                HIHJG(L,4)=HIHJG(L,4)+HI(3,J)*HJ(3,K)*WEIGHT
             ENDDO             
          ENDDO 
       ENDIF
    ENDDO
!
    DO I=1,LL
       HIHJG(I,1)=HIHJG(I,1)*YL
       HIHJG(I,2)=HIHJG(I,2)*YL
       HIHJG(I,3)=HIHJG(I,3)*YL
       HIHJG(I,4)=HIHJG(I,4)*YL
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE QUADS_MFS_GAUSSG_V2
!
    SUBROUTINE QUADS_MFS_ADDBAN_V2(II,JJ,NPOR1,NPOR2,IND1_X,IND1_Y,IND2_X,IND2_Y,S)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II,JJ,NPOR1,NPOR2
    LOGICAL            :: IND1_X,IND1_Y,IND2_X,IND2_Y
    REAL(8)            :: S(NPOR_MAX*NPOR_MAX,4)
!   ASSISTANT VARIABLES
    INTEGER :: N,M
!
    INTEGER :: I,J,K,L
    INTEGER :: KI,KJ
    INTEGER :: LL,LL1,KK1,KK2
    INTEGER :: INDK,INDL
!
    N         =  NDMFS(1,II)
    M         =  NDMFS(1,JJ)
!
    DO I=1,2
       DO J=1,2
          K=1
          INDK=1
          IF((IND1_X).AND.(I.EQ.1)) THEN
             K=2
             INDK=2
          ENDIF
          IF((IND1_Y).AND.(I.EQ.2)) THEN
             K=2
             INDK=2
          ENDIF          
          DO WHILE(K.LE.NPOR1)
             L=1
             INDL=1
             IF((IND2_X).AND.(J.EQ.1)) THEN
                 L=2
                 INDL=2
             ENDIF
             IF((IND2_Y).AND.(J.EQ.2)) THEN
                 L=2
                 INDL=2
             ENDIF 
             DO WHILE(L.LE.NPOR2)
                KI=ID(I,N)+K-INDK
                KJ=ID(J,M)+L-INDL
                KK1=(I-1)*2+J
                KK2=(K-1)*NPOR2+L
                IF ((N.NE.M).OR.(KJ.GE.KI)) THEN
                   IF(KJ.LE.KI) THEN
!                     LL=AMAXA(KI)+ABS(KJ-KI)
                      LL1=MAXA_COLSOL(KI)+ABS(KJ-KI)
                   ELSE
!                     LL=AMAXA(KJ)+ABS(KJ-KI)
                      LL1=MAXA_COLSOL(KJ)+ABS(KJ-KI)
                   ENDIF
!                  A(LL)=A(LL)+S(KK1,KK2)
                   K_COLSOL(LL1)=K_COLSOL(LL1)+S(KK2,KK1)
                   KMATRIX(LL1)=KMATRIX(LL1)+S(KK2,KK1)
                ENDIF
                L=L+1
             ENDDO
             K=K+1
          ENDDO
       ENDDO
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE QUADS_MFS_ADDBAN_V2
!
!**************************************************   version 2
!
    SUBROUTINE ADDBAN (S,NND,LMM,ND,IND)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To assemble upper triangular element stiffness into             .
! .   compacted global stiffness                                      .
! .                                                                   .
! .         A = GLOBAL STIFFNESS                                      .
! .         S = ELEMENT STIFFNESS                                     .
! .         ND = DEGREES OF FREEDOM IN ELEMENT STIFFNESS              .
! .                                                                   .
! .                   S(1)        S(2)        S(3)        . . .       .
! .         S   =                 S(ND+1)     S(ND+2)     . . .       .
! .                                           S(2*ND)     . . .       .
! .                                                       . . .       .
! .                                                                   .
! .                                                                   .
! .                   A(1)        A(3)        A(6)        . . .       .
! .         A   =                 A(2)        A(5)        . . .       .
! .                                           A(4)        . . .       .
! .                                                       . . .       .
! .                                                                   .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: ND,NND
      REAL(8)    :: S(NND)
      INTEGER    :: LMM(ND)
      INTEGER    :: NDI, I,II, MI, KS, J, JJ, IJ, KK, KSS
      INTEGER    :: IND
!
      NDI=0
      DO I=1,ND
         II=LMM(I)
         IF (II .GT. 0) THEN
            MI=MAXA_COLSOL(II)
            KS=I
            DO J=1,ND
               JJ=LMM(J)
               IF (JJ .GT. 0) THEN
                  IJ=II - JJ
                  IF (IJ .GE. 0) THEN
                     KK=MI + IJ
                     KSS=KS
                     IF (J.GE.I) KSS=J + NDI
                     IF(IND.EQ.1) THEN
                        K_COLSOL(KK)=K_COLSOL(KK) + S(KSS)
                        KMATRIX(KK)=KMATRIX(KK)+S(KSS)
                     ELSE
                        M_COLSOL(KK)=M_COLSOL(KK) + S(KSS)
                     ENDIF
                  END IF
               END IF
               KS=KS + ND - J
            END DO
         END IF
         NDI=NDI + ND - I
      END DO
!
      RETURN
!}}}
    END SUBROUTINE ADDBAN
!
    SUBROUTINE OUTPUTKMR(ADK,ADM)
!{{{
!--------------------------------------------------
!   OUTPUT STIFFNESS AND MASS MATRIX
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADK,ADM
    INTEGER :: I,J,K
    REAL(8) :: KI(NEQ)
!--- ASSISTANT VARIABLES
    INTEGER :: II,JJ
!
    DO I=1,NEQ
       KI(I)=0.0D0
    ENDDO
!
    DO I=1,NEQ
       JJ=MHT(I)+1
       DO J=1,JJ
          II=MAXA_COLSOL(I)+J-1
          KI(J)=K_COLSOL(II)
       ENDDO
       WRITE(ADK,1000) (KI(K),K=1,JJ)
    ENDDO    
!
     DO I=1,NEQ
       WRITE(13,*) R(I)
     ENDDO
!
     RETURN
!
1000 FORMAT(*(ES16.9,1X))
!}}}
    END SUBROUTINE OUTPUTKMR
!      
    SUBROUTINE ADDBAN_ZLB(ND,S,LMM,IND)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To assemble upper triangular element stiffness into             .
! .   compacted global stiffness                                      .
! .                                                                   .
! .         A = GLOBAL STIFFNESS                                      .
! .         S = ELEMENT STIFFNESS                                     .
! .         ND = DEGREES OF FREEDOM IN ELEMENT STIFFNESS              .
! .                                                                   .
! .                   S(1)        S(2)        S(3)        . . .       .
! .         S   =                 S(ND+1)     S(ND+2)     . . .       .
! .                                           S(2*ND)     . . .       .
! .                                                       . . .       .
! .                                                                   .
! .                                                                   .
! .                   A(1)        A(3)        A(6)        . . .       .
! .         A   =                 A(2)        A(5)        . . .       .
! .                                           A(4)        . . .       .
! .                                                       . . .       .
! .                                                                   .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: ND
   REAL(8),DIMENSION(ND,ND),INTENT(IN) :: S
   INTEGER,DIMENSION(ND),INTENT(IN) :: LMM
   INTEGER :: IND
!
   INTEGER :: I,J,L,M,K
!
   DO I=1,8
      L=LMM(I)
      IF(L.GT.0) THEN
        DO J=1,8
           M=LMM(J)
           IF(M.GT.0) THEN
              K=M-L
              IF(K.GE.0) THEN
                K=MAXA_COLSOL(M)+K
                IF(IND.EQ.1) THEN
                   K_COLSOL(K)=K_COLSOL(K)+S(I,J)
                   KMATRIX(K)=KMATRIX(K)+S(I,J)
                ELSE
                   M_COLSOL(K)=M_COLSOL(K)+S(I,J)
                ENDIF                
              END IF
           ENDIF
        ENDDO
      ENDIF
   ENDDO
!
   RETURN
!}}}
    END SUBROUTINE ADDBAN_ZLB
!
    SUBROUTINE G_KCS(DIMEN,ITYPE,E,PR)
!{{{
!--------------------------------------------------
!   GENERATE CONSTITUTIVE MATRIX
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: DIMEN,ITYPE
    REAL(8),INTENT(IN)  :: E,PR
!
    INTEGER :: I
    INTEGER :: ALLocateStatus
!
    SELECT CASE(DIMEN)
    CASE(2)
       ALLOCATE (KCS(3),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!
       DO I=1,3
          KCS(I)=0.
       ENDDO
!
       SELECT CASE(ITYPE)
       CASE(1)
          KCS(1)=E*(1-PR)/(1+PR)/(1-2*PR)
          KCS(2)=E*PR/(1+PR)/(1-2*PR)
          KCS(3)=E/2/(1+PR)
       CASE(2)
          KCS(1)=E/(1-PR*PR)
          KCS(2)=E*PR/(1-PR*PR)
          KCS(3)=E/2/(1+PR)
       CASE DEFAULT
          WRITE(*,*) "ERROR, ITYPE.GT.2!"
          STOP
       END SELECT
    CASE DEFAULT
         WRITE(*,*) "ERROR, DIMEN.GT.2!"
         STOP
    END SELECT
!
    RETURN
!}}}
    END SUBROUTINE G_KCS
!
end module GENERATEM
