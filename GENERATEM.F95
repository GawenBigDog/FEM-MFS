module GENERATEM
    USE WAREHOUSE
    USE SOLVERM
    USE GAUSSM
    implicit none
    PUBLIC  GENERATE
!   GENERATE LOAD VECTOR(MFS)
    PRIVATE  MFS_LINE_G
    PRIVATE  MFS_LINE_STIFF
    PRIVATE  MFS_LINE_STDM
!   GENERATE STIFF MATRIX(FEM)
    PRIVATE TRUSS_G
    PRIVATE QUADS_G
!   GENERATE STIFF MATRIX(MFS)
    PRIVATE TRUSS_MFS_G
    PRIVATE TRUSS_MFS_ITP
    PRIVATE TRUSS_MFS_STIFF
    PRIVATE TRUSS_MFS_STDM
    PRIVATE TRUSS_MFS_BODY
    PUBLIC  TRUSS_MFS_HIHJ
!
    PRIVATE QUADS_MFS_G
!   VERSION 1
    PRIVATE QUADS_MFS_STIFF
    PRIVATE QUADS_MFS_GAUSSG
    PRIVATE QUADS_MFS_STDM
    PRIVATE QUADS_MFS_ADDBAN
!   VERSION 2
    PRIVATE QUADS_MFS_STIFF_V2
    PRIVATE QUADS_MFS_GAUSSG_V2
    PRIVATE QUADS_MFS_ADDBAN_V2
! FUNCTION
    PUBLIC  QUADS_MFS_HIHJ
    PUBLIC  QUADS_MFS_HIHJ_V2
! GENERATE CONSTITUTIVE MATRIX
    PRIVATE G_KCS 
!    
    PRIVATE ADDBAN
    PRIVATE ADDBAN_ZLB
    PRIVATE OUTPUTKM

    CONTAINS
!
    SUBROUTINE GENERATE(ADK,ADM,ADELE,ADBC,GNPAR,GMATP,GKMATRIX,GMMATRIX,GMAXA,GR,GRODE,GU)
!{{{
!--------------------------------------------------
!   GENERATE 
!     SELECT GENERATOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADK,ADM,ADELE,ADBC
!
    INTEGER,DIMENSION(7),     INTENT(INOUT) :: GNPAR
    REAL(8),DIMENSION(7,NMAT),INTENT(INOUT) :: GMATP
    REAL(8),DIMENSION(NWK),   INTENT(INOUT) :: GKMATRIX,GMMATRIX
    INTEGER,DIMENSION(NEQ+1), INTENT(IN) :: GMAXA
    REAL(8),DIMENSION(NEQ),   OPTIONAL   :: GR
    REAL(8),DIMENSION(NEQ),   INTENT(INOUT),OPTIONAL :: GRODE
    REAL(8),DIMENSION(NEQ),   INTENT(IN),   OPTIONAL :: GU      
!
    INTEGER :: I,J,BcType,EleType
!   ASSISTANT VARIABLES
    INTEGER :: ALLocateStatus
!-------------------------------------------------------------
    REWIND(ADELE)
    REWIND(ADBC)
!
    DO I=1,NEGROUP
       READ(ADELE) (GNPAR(J),J=1,7)
       EleType=GNPAR(1)
       SELECT CASE(EleType)
       CASE(1)
            CALL TRUSS_G(ADELE,GNPAR(1),GMATP(1,1),GKMATRIX(1),GMMATRIX(1),GMAXA(1),GRODE(1),GU(1))
       CASE(2)
            CALL QUADS_G(ADELE,GNPAR(1),GMATP(1,1),GKMATRIX(1),GMMATRIX(1),GMAXA(1),GRODE(1),GU(1))           
       CASE(3)
            CALL TRUSS_MFS_G(ADELE,GNPAR(1),GMATP(1,1),GKMATRIX(1),GMMATRIX(1),GMAXA(1))
       CASE(4)
            CALL QUADS_MFS_G(ADELE,ADBC,GNPAR(1),GMATP(1,1),GKMATRIX(1),GMMATRIX(1),GMAXA(1),GR(1))
       CASE DEFAULT
            WRITE(*,*), "ERROR, WRONG INPUT IN GENERATE ELEMENT PROCEDURE"
            STOP
       END SELECT       
    END DO
!
    CALL OUTPUTKM(ADK,ADM,GKMATRIX(1),GMMATRIX(1),GMAXA(1))
    DO I=1,NEQ
       WRITE(13,*) GR(I)
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE GENERATE
!
    SUBROUTINE TRUSS_G(ADELE,GNPAR,GMATP,GKMATRIX,GMMATRIX,GMAXA,GRODE,GU)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      MASS      MATRIX
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE
    INTEGER :: NumEle,Mtype
!
!
    INTEGER,DIMENSION(7),     INTENT(INOUT) :: GNPAR
    REAL(8),DIMENSION(7,NMAT),INTENT(INOUT) :: GMATP
    REAL(8),DIMENSION(NWK),   INTENT(INOUT) :: GKMATRIX,GMMATRIX
    INTEGER,DIMENSION(NEQ+1), INTENT(IN) :: GMAXA
    REAL(8),DIMENSION(NEQ),   INTENT(INOUT),OPTIONAL :: GRODE
    REAL(8),DIMENSION(NEQ),   INTENT(IN),   OPTIONAL :: GU    
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
    NumEle = GNPAR(3)
    Mtype  = GNPAR(2)
!   SET MATERIAL PARAMETERS
    E    = GMATP(1,Mtype)
    DEN  = GMATP(3,Mtype)
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
       CALL ADDBAN(GKMATRIX(1),GMAXA(1),S(1),21,LM(1,N),6)
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
       CALL ADDBAN(GMMATRIX(1),GMAXA(1),SMASS(1),21,LM(1,N),6)
!   CALCULATE K*U
       IF(DYNAMIC.EQ.1) THEN
         FORCE=0        
         DO L=1,3
            I=LM(L,N)
            IF ((I.GT.0).AND.(GU(I).NE.0)) FORCE = FORCE + ST(L)*GU(I)*E*AREA
            J=LM(L+3,N)
            IF ((J.GT.0).AND.(GU(J).NE.0)) FORCE = FORCE + ST(L+3)*GU(J)*E*AREA
         END DO 
         IF(FORCE.NE.0) THEN     
           DO L=1,3
              I=LM(L,N)
              IF(I.GT.0) GRODE(I)=GRODE(I)+FORCE*D(L)/XL
              J=LM(L+3,N)
              IF(J.GT.0) GRODE(J)=GRODE(J)-FORCE*D(L)/XL
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
    SUBROUTINE QUADS_G(ADELE,GNPAR,GMATP,GKMATRIX,GMMATRIX,GMAXA,GRODE,GU)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      MASS      MATRIX
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE
    INTEGER :: ITYPE,Mtype,NumEle,NINTG
!
    INTEGER,DIMENSION(7),     INTENT(INOUT) :: GNPAR
    REAL(8),DIMENSION(7,NMAT),INTENT(INOUT) :: GMATP
    REAL(8),DIMENSION(NWK),   INTENT(INOUT) :: GKMATRIX,GMMATRIX
    INTEGER,DIMENSION(NEQ+1), INTENT(IN) :: GMAXA
    REAL(8),DIMENSION(NEQ),   INTENT(INOUT),OPTIONAL :: GRODE
    REAL(8),DIMENSION(NEQ),   INTENT(IN),   OPTIONAL :: GU    
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
    ITYPE  = GNPAR(2)
    Mtype  = GNPAR(3)
    NumEle = GNPAR(4)
    NINTG  = GNPAR(5)
!   SET MATERIAL PARAMETERS
    E    = GMATP(1,Mtype)
    PR   = GMATP(2,Mtype)
    DEN  = GMATP(3,Mtype)
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
       CALL ADDBAN_ZLB(8,S,LM(:,N),GKMATRIX(1),GMAXA(1))  
!   MASS MATRIX: LUMPED MASS  
       DO I=1,8
          DO J=1,8
             IF (J.NE.I) THEN
               SMASS(I,I)=SMASS(I,I)+SMASS(I,J)
               SMASS(I,J)=0
             ENDIF
          ENDDO
       ENDDO
       CALL ADDBAN_ZLB(8,SMASS,LM(:,N),GMMATRIX(1),GMAXA(1))
!   CALCULATE K*U
       IF(DYNAMIC.EQ.1) THEN
          DO I=1,8
             FORCE=0
             KK=LM(I,N)
             IF(KK.GT.0) THEN
               DO J=1,8
                 K=LM(J,N)
                 IF(K.GT.0) THEN
                    FORCE=FORCE+S(I,J)*GU(K)
                 ENDIF
               ENDDO
               GRODE(KK)=GRODE(KK)+FORCE
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
    SUBROUTINE TRUSS_MFS_G(ADELE,GNPAR,GMATP,GKMATRIX,GMMATRIX,GMAXA)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      MASS      MATRIX
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE
!
    INTEGER,DIMENSION(7),     INTENT(INOUT) :: GNPAR
    REAL(8),DIMENSION(7,NMAT),INTENT(INOUT) :: GMATP
    REAL(8),DIMENSION(NWK),   INTENT(INOUT) :: GKMATRIX,GMMATRIX
    INTEGER,DIMENSION(NEQ+1), INTENT(IN)    :: GMAXA        
!   ASSISTANT VARIABLES
    INTEGER :: N,M
    INTEGER :: I,J,K,L,II,JJ
    INTEGER :: NumEle,Mtype,NEleAss
    REAL(8) :: E,DEN,AREA
    INTEGER :: ALLocateStatus
!   ASSISTANT VARIABLES
    INTEGER :: NUM,NUM1
    INTEGER :: KI
    REAL(8) :: Radius
    REAL(8) :: S        
    REAL(8) :: HI,DHI
    REAL(8) :: XX
    INTEGER :: NPOR
!   READ 
    NEleAss= GNPAR(4)
    NumEle = GNPAR(3)
    Mtype  = GNPAR(2)
!   SET MATERIAL PARAMETERS
    E    = GMATP(1,Mtype)
    DEN  = GMATP(3,Mtype)
    AREA = GMATP(4,Mtype)
!   ALLOCATE VARIABLES
    ALLOCATE (NDMFS(2,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(2,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
    ALLOCATE (MFS_ITP_XYZDIV(3,NITP),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!--------------------------------------------------
        DO I=1,NumEle
           NDMFS(1,I)=0
           NDMFS(2,I)=0
           ELEPAR(1,I)=0.
           ELEPAR(2,I)=0.
           DO J=1,NEleAss
              LM(J,I)=0
           ENDDO
        END DO
!
        DO I=1,NITP
           MFS_ITP_XYZDIV(1,I)=0
           MFS_ITP_XYZDIV(2,I)=0
           MFS_ITP_XYZDIV(3,I)=0
        ENDDO
!--------------------------------------------------
    READ(ADELE) ((NDMFS(I,J),I=1,2),J=1,NumEle),((ELEPAR(J,I),J=1,2),I=1,NumEle),((LM(J,I),J=1,NEleAss),I=1,NumEle)
!   GENERATE
!   STIFFNESS MATRIX
    DO I=1,NumEle
!
!   GENERATE LOAD VECTOR
!
!   Neumann boundary
!
       IF(NDMFS(2,I).EQ.2) THEN
         N=NDMFS(1,I)
         DO II=1,3
            KI=ID(II,N)
            IF(KI.NE.0) THEN
               XX=X(N)
               DO JJ=1,NDOF(II,N)
                  NUM=KI+JJ-1
                  NPOR=JJ
                  HI=TRUSS_MFS_HIHJ(NEleAss,I,NPOR,XX,1)
                  R(NUM)=R(NUM)*HI*AREA
               ENDDO
            ENDIF
         ENDDO
       ENDIF
!       
!   Dirichlet boundary
!
       IF(NDMFS(2,I).EQ.1) THEN
          N=NDMFS(1,I)
          DO II=1,3
             KI=ID(II,N)
             IF(KI.NE.0) THEN
                XX=X(N)
                DO JJ=1,NDOF(II,N)
                   NPOR=JJ+1
                   DHI=TRUSS_MFS_HIHJ(NEleAss,I,NPOR,XX,2)
                   NUM=KI+JJ-1
                   R(NUM) = R(NUM)*DHI*AREA
                ENDDO
             ENDIF
          ENDDO
       ENDIF
!
!      CALCULATE LOAD VECTOR CONTRIBUTION FROM DIRICHLET BOUNDARY 
!
       N=NDMFS(1,I)
       DO II=1,3
          KI=ID(II,N)
          IF(KI.NE.0) THEN
             DO J=1,NEleAss
                K=LM(J,I)
                IF(K.NE.0) THEN
                   IF(NDMFS(2,K).EQ.1) THEN
                      M=NDMFS(1,K)
                      CALL TRUSS_MFS_ITP(I,K,NumEle,ELEPAR(1,1),NDMFS(1,1),MFS_ITP_XYZDIV(1,1))
                      DO L=1,NDOF(II,N)
                         IF(NDMFS(2,I).EQ.1) THEN
                            NPOR=L+1
                         ELSE
                            NPOR=L
                         ENDIF
!
                         NUM=KI+L-1
                         S=0
                         CALL TRUSS_MFS_STDM(NPOR,1,I,K,NumEle,NEleAss,NDMFS(1,1),MFS_ITP_XYZDIV(1,1),S)
                         R(NUM)=R(NUM)-S*RESS(II,M)*AREA
                      ENDDO
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
       ENDDO
!   CALCULATE STIFFNESS MATRIX
       DO J=1,NEleAss
          K=LM(J,I)
          IF ((K.NE.0).AND.(K.GE.I)) THEN
             CALL TRUSS_MFS_ITP(I,K,NumEle,ELEPAR(1,1),NDMFS(1,1),MFS_ITP_XYZDIV(1,1))          
             CALL TRUSS_MFS_STIFF(GKMATRIX,GMAXA,I,K,E,DEN,AREA,NumEle,NEleAss,NDMFS(1,1),MFS_ITP_XYZDIV(1,1))
!       
!   Body Force
!
             IF(I.EQ.K) THEN
                CALL TRUSS_MFS_BODY(R(1),I,NumEle,NEleAss,NDMFS(1,1),MFS_ITP_XYZDIV(1,1),AREA)
             ENDIF
          ENDIF
       ENDDO 
!
    ENDDO              
!
    DEALLOCATE(NDMFS)
    DEALLOCATE(ELEPAR)
    DEALLOCATE(LM)
    DEALLOCATE(MFS_ITP_XYZDIV)
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_MFS_G
!
    SUBROUTINE TRUSS_MFS_ITP(II,JJ,NumEle,AELEPAR,ANDMFS,AMFS_ITP_XYDIV)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! . CALCULATE  COORDINATES OF INTEGRATION POINTS IN THE OVERLAPPING AREA  .
! . INPUT:                                                                .
! .        II,  JJ                                                        .
! . OUTPUT:                                                               .
! .        MFS_ITPXYDIV:  COORDINATES OF INTEGRATION POINTS               .   
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II,JJ
    INTEGER            :: NumEle
    REAL(8)            :: AELEPAR(2,NumEle),AMFS_ITP_XYDIV(3,NITP)
    INTEGER            :: ANDMFS(2,NumEle)
!   ASSISTANT VARIABLES
    INTEGER            :: N,M
    INTEGER            :: I,J
    INTEGER            :: INDN
    REAL(8)            :: ANGLE,Radius1,Radius2
    REAL(8)            :: RAD,DIV,DIVX,DIVY
    REAL(8)            :: XX,YY   
!
    N         =  ANDMFS(1,II)
    M         =  ANDMFS(1,JJ)
    INDN      =  ANDMFS(2,II)
    ANGLE     =  AELEPAR(1,II)
    Radius1   =  AELEPAR(2,II)
    Radius2   =  AELEPAR(2,JJ)
!
!   DETERMINE THE OVERLAPPING AREA
!
    IF(N.EQ.M) THEN
       IF(INDN.NE.0) THEN
          RAD=ANGLE*Pi
          DIV=Radius1/NITP
          XX=X(N)
          YY=Y(N)
          DO I=1,NITP
             AMFS_ITP_XYDIV(1,I)=XX+DIV*(I-1./2)*COS(RAD)
             AMFS_ITP_XYDIV(2,I)=YY+DIV*(I-1./2)*SIN(RAD)
             AMFS_ITP_XYDIV(3,I)=DIV
          ENDDO          
       ELSE
          RAD=ANGLE*Pi
          DIV=2*Radius1/NITP
          XX=X(N)-COS(RAD)*Radius1
          YY=Y(N)-SIN(RAD)*Radius1
          DO I=1,NITP
             AMFS_ITP_XYDIV(1,I)=XX+DIV*(I-1./2)*COS(RAD)
             AMFS_ITP_XYDIV(2,I)=YY+DIV*(I-1./2)*SIN(RAD)
             AMFS_ITP_XYDIV(3,I)=DIV
          ENDDO
       ENDIF
    ELSE
       RAD=ANGLE*Pi
       DIV=(X(N)-X(M))*(X(N)-X(M))+(Y(N)-Y(M))*(Y(N)-Y(M))
       DIV=((Radius1+Radius2)-SQRT(DIV))/NITP
       DIVX=(X(M)-X(N))*COS(RAD)
       DIVY=(Y(M)-Y(N))*COS(RAD)
       IF(DIVX.GE.0) THEN
          XX=X(M)-Radius2*COS(RAD)
       ELSE
          XX=X(N)-Radius1*COS(RAD)
       ENDIF
!
       IF(DIVY.GE.0) THEN
          YY=Y(M)-Radius2*SIN(RAD)
       ELSE
          YY=Y(N)-Radius1*SIN(RAD)
       ENDIF
!
       DO I=1,NITP
          AMFS_ITP_XYDIV(1,I)=XX+DIV*(I-1./2)*COS(RAD)
          AMFS_ITP_XYDIV(2,I)=YY+DIV*(I-1./2)*SIN(RAD)
          AMFS_ITP_XYDIV(3,I)=DIV
       ENDDO       
    ENDIF    
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_MFS_ITP
!
    SUBROUTINE TRUSS_MFS_STIFF(A,AMAXA,II,JJ,E,DEN,AREA,NumEle,NEleAss,ANDMFS,AMFS_ITP_XYDIV)
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
    REAL(8)            :: A(NWK)
    INTEGER            :: AMAXA(NEQ+1)
    REAL(8)            :: E,DEN,AREA
    INTEGER            :: NumEle,NEleAss
    INTEGER            :: ANDMFS(2,NumEle)
    REAL(8)            :: AMFS_ITP_XYDIV(3,NITP)
!   ASSISTANT VARIABLES
    INTEGER            :: N,M
    INTEGER            :: KI,KJ
    INTEGER            :: I,J,K,PORJ,PORK,L
    INTEGER            :: IND
    REAL(8)            :: S
!
    N       =  ANDMFS(1,II)
    M       =  ANDMFS(1,JJ)
!
    DO I=1,3
       KI=ID(I,N)
       IF(KI.NE.0) THEN                  
         DO J=1,NDOF(I,N)
          IND=ANDMFS(2,II)
          IF(IND.EQ.1) THEN
             PORJ=J+1  
          ELSE
             PORJ=J
          ENDIF
          KI=ID(I,N)+J-1
          DO K=1,NDOF(I,M)
             IND=ANDMFS(2,JJ)
             IF(IND.EQ.1) THEN
                PORK=K+1
             ELSE
                PORK=K
             ENDIF
!   DETERMINE THE ADDRESS
             KJ=ID(I,M)+K-1
!   CALCULATE ENTRIES OF THE STIFFNESS MATRIX
             IF ((N.NE.M).OR.(KJ.GE.KI)) THEN                
                S=0
                CALL TRUSS_MFS_STDM(PORJ,PORK,II,JJ,NumEle,NEleAss,ANDMFS(1,1),AMFS_ITP_XYDIV(1,1),S)
!   ASSEMBLE STIFFNESS MATRIX
                IF(KJ.LE.KI) THEN
                  L=AMAXA(KI)+ABS(KJ-KI)
                ELSE
                  L=AMAXA(KJ)+ABS(KJ-KI)
                ENDIF
                A(L)=AREA*E*S  
             ENDIF             
!
          ENDDO
         ENDDO
       ENDIF    
    ENDDO  
! 
    RETURN
!}}}
    END SUBROUTINE TRUSS_MFS_STIFF
!
    SUBROUTINE TRUSS_MFS_STDM(JJ,KK,IE,JE,NumEle,NEleAss,ANDMFS,AMFS_ITP_XYDIV,STIFF)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE ENTRY OF THE STIFFNESS MATRIX                         .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER :: JJ,KK,IE,JE
    INTEGER :: NumEle, NEleAss
    INTEGER :: ANDMFS(2,NumEle)
    REAL(8) :: AMFS_ITP_XYDIV(3,NITP)
    REAL(8) :: STIFF
!   ASSISTANT VARIABLES
    INTEGER :: N
    INTEGER :: I
    INTEGER :: IND 
    REAL(8) :: HI,HJ
    REAL(8) :: DHIDX,DHJDX
    REAL(8) :: XX,DIV
!
    STIFF=0
!
    DO I=1,NITP
       XX    = AMFS_ITP_XYDIV(1,I)
       DIV   = AMFS_ITP_XYDIV(3,I)
       DHIDX = TRUSS_MFS_HIHJ(NEleAss,IE,JJ,XX,2)
       DHJDX = TRUSS_MFS_HIHJ(NEleAss,JE,KK,XX,2)
       STIFF = STIFF+DHIDX*DHJDX*DIV
    ENDDO
!
!   CONTRIBUTIONS OF STIFFNESS FROM BOUNDARY
!
    IND=ANDMFS(2,IE)                      
!
    IF(IND.EQ.1) THEN
       IND=ANDMFS(2,JE)
       IF(IND.EQ.1) THEN
          N  = ANDMFS(1,IE)
          XX = X(N)
          HI = TRUSS_MFS_HIHJ(NEleAss,IE,JJ,XX,1)
          HJ = TRUSS_MFS_HIHJ(NEleAss,JE,KK,XX,1)
          DHIDX = TRUSS_MFS_HIHJ(NEleAss,IE,JJ,XX,2)
          DHJDX = TRUSS_MFS_HIHJ(NEleAss,JE,KK,XX,2)
          STIFF = STIFF+DHIDX*HJ+HI*DHJDX
       ENDIF
    ENDIF                                      
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_MFS_STDM
!
    SUBROUTINE TRUSS_MFS_BODY(AR,II,NumEle,NEleAss,ANDMFS,AMFS_ITP_XYDIV,AREA)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE  BODY FORCE CONTRIBUTION TO LOAD VECTOR               .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER :: II
    INTEGER :: NumEle,NEleAss
    REAL(8) :: AR(NEQ)
    INTEGER :: ANDMFS(2,NumEle) 
    REAL(8) :: AMFS_ITP_XYDIV(3,NITP)
    REAL(8) :: AREA
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: N,NUM,NPOR
    INTEGER :: IND
    REAL(8) :: XX,DIV
    REAL(8) :: HI
!
    DO I=1,3
       N   = ANDMFS(1,II)
       NUM=ID(I,N)
       IF(NUM.NE.0) THEN
          DO J=1,NITP
             XX  = AMFS_ITP_XYDIV(1,J)
             DIV = AMFS_ITP_XYDIV(3,J)
             DO K=1,NDOF(I,N)
                IND = ANDMFS(2,II)
                IF(IND.EQ.1) THEN
                   NPOR=K+1
                ELSE
                   NPOR=K
                ENDIF
                HI    = TRUSS_MFS_HIHJ(NEleAss,II,NPOR,XX,1)
                NUM   = ID(I,N)+K-1
                R(NUM)= R(NUM) + XX*HI*DIV*AREA
             ENDDO     
          ENDDO
       ENDIF
!
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_MFS_BODY
!
    FUNCTION TRUSS_MFS_HIHJ(NEleAss,II,NPOR,XX,IND) RESULT(H)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE VALUE OF HI,HJ,DHI,DHJ                                .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NEleAss
    INTEGER, INTENT(IN) :: II,NPOR
    INTEGER, INTENT(IN) :: IND
    REAL(8), INTENT(IN) :: XX
    REAL(8)             :: H
!   ASSISTANT VARIABLES
    INTEGER :: I,J
    INTEGER :: N,M
    REAL(8) :: W1,DW1,W2,DW2
    REAL(8) :: H0,H1,H2,DH0,DH1,DH2
    REAL(8) :: Radius1,Radius
    REAL(8) :: S
!
    Radius1=ELEPAR(2,II)
    N=NDMFS(1,II)
    S=(XX-X(N))/Radius1
!   CALCULATE W1,W2,DW1,DW2
    W1  = 1-6*S*S+8*S*S*ABS(S)-3*S*S*S*S
    DW1 = (-12*S+24*S*ABS(S)-12*S*S*S)/Radius1
!
    W2=0
    DW2=0
    DO I=1,NEleAss
       J=LM(I,II)
       IF (J.NE.0) THEN
          M=NDMFS(1,J)
          Radius=ELEPAR(2,J)
          S=(XX-X(M))/Radius
          IF(ABS(S).LE.1) THEN
             W2=W2+1-6*S*S+8*S*S*ABS(S)-3*S*S*S*S
             DW2=DW2+(-12*S+24*S*ABS(S)-12*S*S*S)/Radius
          ENDIF
       ENDIF
    ENDDO
!
    SELECT CASE(IND)
!   CALCULATE H0,H1,H2
    CASE(1)
       SELECT CASE(NPOR)
       CASE(1)
          H0=W1/W2
          H=H0
       CASE(2)
          H1=W1/W2*(XX-X(N))/Radius1
          H=H1
       CASE(3)
          H2=W1/W2*(XX-X(N))*(XX-X(N))/Radius1/Radius1
          H=H2
       CASE DEFAULT
          WRITE(*,*) "ERROR, NPOR SHOULD LT 3!"
          STOP
       END SELECT
!   CALCULATE DH0,DH1,DH2
    CASE(2)
       SELECT CASE(NPOR)
       CASE(1)
          DH0=(DW1*W2-W1*DW2)/W2/W2
          H=DH0
       CASE(2)
          DH0=(DW1*W2-W1*DW2)/W2/W2
          H0 =W1/W2
          DH1=DH0*(XX-X(N))/Radius1+H0/Radius1
          H=DH1
       CASE(3)
          DH0=(DW1*W2-W1*DW2)/W2/W2
          H0 =W1/W2
          DH2=DH0*(XX-X(N))*(XX-X(N))/Radius1/Radius1+H0*2*(XX-X(N))/Radius1/Radius1
          H=DH1
       CASE DEFAULT
          WRITE(*,*) "ERROR, NPOR SHOULD LT 4!"
          STOP
       END SELECT
    CASE DEFAULT
       WRITE(*,*) "ERROR, IND SHOULD LT 3!"
       STOP
    END SELECT
!
    RETURN
!}}}
    END FUNCTION TRUSS_MFS_HIHJ
!}}}
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *   
    SUBROUTINE QUADS_MFS_G(ADELE,ADBC,GNPAR,GMATP,GKMATRIX,GMMATRIX,GMAXA,GR)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      MASS      MATRIX
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE,ADBC
!
    INTEGER,DIMENSION(7),     INTENT(INOUT) :: GNPAR
    REAL(8),DIMENSION(7,NMAT),INTENT(INOUT) :: GMATP
    REAL(8),DIMENSION(NWK),   INTENT(INOUT) :: GKMATRIX,GMMATRIX
    INTEGER,DIMENSION(NEQ+1), INTENT(IN)    :: GMAXA   
    REAL(8),DIMENSION(NEQ)                  :: GR 
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L
    INTEGER :: NEleAss,NumEle,Mtype,ITYPE
    REAL(8) :: E,PR,DEN,THIC
    INTEGER :: BcType
!   ASSISTANT VARIABLES
    INTEGER :: KI
    INTEGER :: ALLocateStatus
!
    REAL(8) :: XX=0,YY=0.5
!   READ 
    NEleAss= GNPAR(5)
    NumEle = GNPAR(4)
    Mtype  = GNPAR(3)
    ITYPE  = GNPAR(2)
!   SET MATERIAL PARAMETERS
    E    = GMATP(1,Mtype)
    PR   = GMATP(2,Mtype)
    DEN  = GMATP(3,Mtype)
    THIC = GMATP(4,Mtype)
!   ALLOCATE VARIABLES
    ALLOCATE (NDMFS(5,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(3,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!
    ALLOCATE (XABSC_WEIG_BC(NITP_BC,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
    ALLOCATE (XABSC_WEIG_IN(NITP_IN,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!--------------------------------------------------
        DO I=1,NumEle
           DO J=1,5
              NDMFS(J,I)=0
           ENDDO
!
           DO J=1,3
              ELEPAR(J,I)=0.
           ENDDO
!
           DO J=1,NEleAss
              LM(J,I)=0
           ENDDO
        END DO
!
        DO I=1,NITP_BC
           XABSC_WEIG_BC(I,1)=0
           XABSC_WEIG_BC(I,2)=0
        ENDDO
!
        DO I=1,NITP_IN
           XABSC_WEIG_IN(I,1)=0
           XABSC_WEIG_IN(I,2)=0
        ENDDO
!--------------------------------------------------
    READ(ADELE) ((NDMFS(I,J),I=1,5),J=1,NumEle),((ELEPAR(J,I),J=1,3),I=1,NumEle),((LM(J,I),J=1,NEleAss),I=1,NumEle)
!   GENERATE CONSTITUTIVE MATRIX
    CALL G_KCS(2,ITYPE,E,PR)
!   GENERATE LOAD VECTOR  
    CALL GAULEG(NITP_BC,XABSC_WEIG_BC(1,1),XABSC_WEIG_BC(1,2))
!
    READ(ADBC) (GNPAR(J),J=1,7)
!
    ALLOCATE (STIFF_IJ(4,NPOR_MAX*NPOR_MAX),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!
    BcType=GNPAR(1)
    SELECT CASE(BcType)
    CASE(1)
        CALL MFS_LINE_G(ADBC,GNPAR(1),NEleAss,THIC,GR(1))
    CASE DEFAULT
        WRITE(*,*), "ERROR, WRONG INPUT IN GENERATE BOUNDARY PROCEDURE"
        STOP
    END SELECT
!   GENERATE STIFFNESS MATRIX 
    CALL GAULEG(NITP_IN,XABSC_WEIG_IN(1,1),XABSC_WEIG_IN(1,2)) 
!
    DO I=1,NumEle
       WRITE(*,*) "Generating Stiff matrix OF ELEMENT  = ",I      
!   CALCULATE STIFFNESS MATRIX
       DO J=1,NEleAss
          K=LM(J,I)
          IF ((K.NE.0).AND.(K.GE.I)) THEN
!   CALCULATE ENTRIES OF STIFFNESS MATRIX  
!             CALL QUADS_MFS_STIFF(GKMATRIX,GMAXA,I,K,NumEle,NEleAss,ELEPAR(1,1),NDMFS(1,1),LM(1,1),THIC)
              CALL QUADS_MFS_STIFF_V2(GKMATRIX(1),GMAXA(1),I,K,NumEle,NEleAss,ELEPAR(1,1),NDMFS(1,1),THIC)
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
    SUBROUTINE MFS_LINE_G(ADBC,GNPAR,NEleAss,THIC,GR)
!{{{
!--------------------------------------------------
!   GENERATE LOAD VECTOR
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADBC,NEleAss
    REAL(8),INTENT(INOUT) :: THIC
    INTEGER,DIMENSION(7),     INTENT(INOUT) :: GNPAR
    REAL(8),DIMENSION(NEQ)                  :: GR
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,JJ,KK
    INTEGER :: NumBc
    INTEGER :: ALLocateStatus
!  
    INTEGER :: L1,L2,N1,N2,L,N
    REAL(8) :: XM,XL,YM,YL,XX,YY
    REAL(8) :: FLOAD(2),FUNC(2)
    INTEGER :: INDBC,IND
    INTEGER :: NUM
    INTEGER :: NPOR   
    REAL(8) :: DIV
    REAL(8) :: S,HI,LENGTH
    REAL(8) :: S1,S2,S3,DHIDX,DHIDY
    REAL(8) :: NX,NY,DISP(2)
!
!   READ 
    NumBc  = GNPAR(2)
!   ALLOCATE VARIABLES
    ALLOCATE (BCED(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (BCIND(4,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (BCPAR(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!
    DO I=1,NumBc
       DO J=1,4
          BCIND(J,I)=0
       ENDDO
       DO J=1,2
          BCED(J,I)=0
          BCPAR(J,I)=0.
       ENDDO
    ENDDO
!
    READ(ADBC) ((BCED(J,I),J=1,2),I=1,NumBc),((BCIND(J,I),J=1,4),I=1,NumBc),((BCPAR(J,I),J=1,2),I=1,NumBc)
!
    DO I=1,NumBc
       WRITE(*,*) "Generating Load vector OF BOUNDARY  = ",I
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
       DO J=1,2
          IF(BCIND(J,I).EQ.1) THEN
             DISP(J)=FLOAD(J)
          ELSE
             DISP(J)=0.
          ENDIF
       ENDDO
!------------------------------------------
       DO J=1,2
          INDBC=BCIND(J,I)
          SELECT CASE(INDBC)
          CASE(2)
             DO K=1,2
                L=BCED(K,I)
                N=NDMFS(1,L)
                IF(NDOF(J,N).GE.1) THEN            
                  DO JJ=1,NDOF(J,N)
                     NUM=ID(J,N)+JJ-1
                     IND=NDMFS(J+1,L)
                     IF(IND.EQ.1) THEN
                        NPOR=JJ+1
                     ELSE
                        NPOR=JJ
                     ENDIF
                     S=0.
                     DO KK=1,NITP_BC
                        XX = XM+XL*XABSC_WEIG_BC(KK,1)
                        YY = YM+YL*XABSC_WEIG_BC(KK,1)
                        DIV= XABSC_WEIG_BC(KK,2)
!
                        HI    = QUADS_MFS_HIHJ(NEleAss,L,NPOR,XX,YY,1)
                        S     = S + HI*FLOAD(J)*DIV
                     ENDDO
                     GR(NUM)=GR(NUM)+S*LENGTH/2.*THIC
                  ENDDO
                ENDIF
             ENDDO
          CASE(1)
               CONTINUE
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
          CASE DEFAULT
             WRITE(*,*) "ERROR, IND.GT.1 IN SUBROUTINE MFS_LINE_G!"
             STOP
          END SELECT
       ENDDO 
!
       DO J=1,2
          DO K=J,2
             L1=BCED(J,I)
             L2=BCED(K,I)
!             CALL MFS_LINE_STIFF(L1,L2,I,NEleAss,NX,NY,XM,XL,YM,YL,LENGTH,THIC)
             CALL MFS_LINE_STIFF_V2(L1,L2,I,NEleAss,NX,NY,XM,XL,YM,YL,LENGTH,THIC)
          ENDDO
       ENDDO
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
    SUBROUTINE MFS_LINE_STIFF(II,JJ,KK,NEleAss,NX,NY,XM,XL,YM,YL,LENGTH,THIC)
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
    REAL(8) :: S(2,2)
    INTEGER :: INDBC1,INDBC2
    INTEGER :: NPOR1,NPOR2
!
    INTEGER :: IND1_X,IND1_Y,IND2_X,IND2_Y 
    REAL(8) :: HIHJ(4,NPOR_MAX*NPOR_MAX)
    REAL(8) :: HI(NPOR_MAX,3),HJ(NPOR_MAX,3)
    REAL(8) :: XTMP,YTMP
!
     INTEGER :: KI,KJ
!
    N1=NDMFS(1,II)
    N2=NDMFS(1,JJ)
!
    INDBC1=BCIND(1,KK)
    INDBC2=BCIND(2,KK)
!
    IF(INDBC1.EQ.1) THEN
       NPOR1=NDOF(1,N1)
       NPOR2=NDOF(1,N2)
    ELSEIF(INDBC2.EQ.1) THEN
       NPOR1=NDOF(2,N1)
       NPOR2=NDOF(2,N2)
    ENDIF    
!
    DO I=1,2
       DO J=1,2
          S(J,I)=0.0
       ENDDO
    ENDDO

!
    IF((INDBC1.EQ.1).OR.(INDBC2.EQ.1)) THEN
       DO I=1,NPOR1
          DO J=1,NPOR2
             CALL MFS_LINE_STDM(I+1,J+1,II,JJ,KK,NEleAss,NX,NY,XM,XL,YM,YL,LENGTH,S,THIC) 
!     ASSEMBLY  
             IF(INDBC1.EQ.2) THEN
                S(1,1)=0.
                S(1,2)=0.
                S(2,1)=0.
             ENDIF
!
             IF(INDBC2.EQ.2) THEN
                S(1,2)=0.
                S(2,1)=0.
                S(2,2)=0.
             ENDIF
!
             DO K=1,2
                KI=ID(K,N1)+I-1
                DO L=1,2
                   KJ=ID(L,N2)+J-1
                   IF ((N1.NE.N2).OR.(KJ.GE.KI)) THEN
                      IF(KJ.LE.KI) THEN
                         LL=MAXA(KI)+ABS(KJ-KI)
                      ELSE
                         LL=MAXA(KJ)+ABS(KJ-KI)
                      ENDIF
                      KMATRIX(LL)=KMATRIX(LL)-S(K,L)
                   ENDIF                   
                ENDDO             
             ENDDO
!
          ENDDO
       ENDDO
    ENDIF
!
    RETURN
!}}}
    END SUBROUTINE MFS_LINE_STIFF
!
    SUBROUTINE MFS_LINE_STIFF_V2(II,JJ,KK,NEleAss,NX,NY,XM,XL,YM,YL,LENGTH,THIC)
!{{{
!--------------------------------------------------
!   GENERATE STIFFNESS MATRIX RESULT FROM BCs
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER :: II,JJ,KK,NEleAss
    REAL(8) :: NX,NY,XM,XL,YM,YL,LENGTH
    REAL(8) :: THIC
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: LL
    INTEGER :: N1,N2
!   REAL(8) :: S(2,2)
    INTEGER :: INDBC1,INDBC2
    INTEGER :: NPOR1,NPOR2
!
    INTEGER :: IND1_X,IND1_Y,IND2_X,IND2_Y 
    REAL(8) :: HIHJ(4,NPOR_MAX*NPOR_MAX)
    REAL(8) :: HI(NPOR_MAX,3),HJ(NPOR_MAX,3)
    REAL(8) :: XTMP,YTMP
!
!    INTEGER :: KI,KJ
!
    N1=NDMFS(1,II)
    N2=NDMFS(1,JJ)
!
    INDBC1=BCIND(1,KK)
    INDBC2=BCIND(2,KK)
!
    IND1_X  =  NDMFS(2,II)
    IND1_Y  =  NDMFS(3,II)
    IND2_X  =  NDMFS(2,JJ)
    IND2_Y  =  NDMFS(3,JJ)
!
    IF (IND1_X.EQ.1) THEN
        NPOR1=NDOF(1,N1)+1
    ELSE
        NPOR1=NDOF(1,N1)
    ENDIF
!
    IF (IND2_X.EQ.1) THEN
        NPOR2=NDOF(1,N2)+1
    ELSE
        NPOR2=NDOF(1,N2)
    ENDIF
!
    DO I=1,NPOR_MAX*NPOR_MAX
       STIFF_IJ(1,I)=0.0D0
       STIFF_IJ(2,I)=0.0D0
       STIFF_IJ(3,I)=0.0D0
       STIFF_IJ(4,I)=0.0D0
    ENDDO
!
    DO I=1,NPOR_MAX*NPOR_MAX
       HIHJ(1,I)=0.0D0
       HIHJ(2,I)=0.0D0
       HIHJ(3,I)=0.0D0
       HIHJ(4,I)=0.0D0
    ENDDO
!
    DO I=1,3
       DO J=1,NPOR_MAX
          HI(J,I)=0.0D0
          HJ(J,I)=0.0D0
       ENDDO
    ENDDO
!
    IF ((INDBC1.EQ.1).OR.(INDBC2.EQ.1)) THEN
      DO I=1,NITP_BC
         XTMP = XM+XL*XABSC_WEIG_BC(I,1)
         YTMP = YM+YL*XABSC_WEIG_BC(I,1)
!
         HI = QUADS_MFS_HIHJ_V2(NEleAss,II,XTMP,YTMP,NPOR1)
         HJ = QUADS_MFS_HIHJ_V2(NEleAss,JJ,XTMP,YTMP,NPOR2)
!
         DO J=1,NPOR1
            DO K=1,NPOR2
               LL=(J-1)*NPOR1+K
               HIHJ(1,LL)=HIHJ(1,LL)+HI(J,1)*HJ(K,2)*XABSC_WEIG_BC(I,2)
               HIHJ(2,LL)=HIHJ(2,LL)+HI(J,1)*HJ(K,3)*XABSC_WEIG_BC(I,2)
               HIHJ(3,LL)=HIHJ(3,LL)+HI(J,2)*HJ(K,1)*XABSC_WEIG_BC(I,2)
               HIHJ(4,LL)=HIHJ(4,LL)+HI(J,3)*HJ(K,1)*XABSC_WEIG_BC(I,2)
            ENDDO
         ENDDO
      ENDDO
!
      DO I=1,NPOR1
         DO J=1,NPOR2
            LL=(I-1)*NPOR1+J
            STIFF_IJ(1,LL)=(HIHJ(1,LL)+HIHJ(3,LL))*NX*KCS(1)+(HIHJ(2,LL)+HIHJ(4,LL))*NY*KCS(3)
            STIFF_IJ(2,LL)=(HIHJ(2,LL)*NX+HIHJ(3,LL)*NY)*KCS(2)+(HIHJ(1,LL)*NY+HIHJ(4,LL)*NX)*KCS(3)
            STIFF_IJ(3,LL)=(HIHJ(1,LL)*NY+HIHJ(4,LL)*NX)*KCS(2)+(HIHJ(2,LL)*NX+HIHJ(3,LL)*NY)*KCS(3)
            STIFF_IJ(4,LL)=(HIHJ(2,LL)+HIHJ(4,LL))*NY*KCS(1)+(HIHJ(1,LL)+HIHJ(3,LL))*NX*KCS(3)
!
            STIFF_IJ(1,LL)=STIFF_IJ(1,LL)*LENGTH/2.0*THIC
            STIFF_IJ(2,LL)=STIFF_IJ(2,LL)*LENGTH/2.0*THIC
            STIFF_IJ(3,LL)=STIFF_IJ(3,LL)*LENGTH/2.0*THIC
            STIFF_IJ(4,LL)=STIFF_IJ(4,LL)*LENGTH/2.0*THIC 
         ENDDO
      ENDDO  
!
      CALL MFS_LINE_ADDBAN(KMATRIX(1),MAXA(1),II,JJ,KK,NPOR1,NPOR2,IND1_X,IND1_Y,IND2_X,IND2_Y,STIFF_IJ(1,1))
    ENDIF
!
    RETURN
!}}}
    END SUBROUTINE MFS_LINE_STIFF_V2
!
    SUBROUTINE MFS_LINE_ADDBAN(A,AMAXA,II,JJ,KK,NPOR1,NPOR2,IND1_X,IND1_Y,IND2_X,IND2_Y,S)
!{{{
    IMPLICIT NONE
    REAL(8)            :: A(NWK)
    INTEGER            :: AMAXA(NEQ+1)
    INTEGER,INTENT(IN) :: II,JJ,KK,NPOR1,NPOR2
    INTEGER            :: IND1_X,IND1_Y,IND2_X,IND2_Y
    REAL(8)            :: S(4,NPOR1*NPOR2)
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L
    INTEGER :: JMIN,JMAX
    INTEGER :: N,M
    INTEGER :: INDBC1,INDBC2
    INTEGER :: INDK,INDL,KK1,KK2
    INTEGER :: KI,KJ
    INTEGER :: LL
!
    INDBC1=BCIND(1,KK)
    INDBC2=BCIND(2,KK)
!
    N         =  NDMFS(1,II)
    M         =  NDMFS(1,JJ)
!
    DO I=1,2
       IF(INDBC1.EQ.2) THEN
          JMIN=1
       ELSE
          JMIN=1
       ENDIF
       IF(INDBC2.EQ.2) THEN
          JMAX=2
       ELSE
          JMAX=2
       ENDIF
!
       DO J=JMIN,JMAX
          K=1
          INDK=1
          IF((IND1_X.EQ.1).AND.(I.EQ.1)) THEN
             K=2
             INDK=2
          ENDIF
          IF((IND1_Y.EQ.1).AND.(I.EQ.2)) THEN
             K=2
             INDK=2
          ENDIF          
          DO WHILE(K.LE.NPOR1)
             L=1
             INDL=1
             IF((IND2_X.EQ.1).AND.(J.EQ.1)) THEN
                 L=2
                 INDL=2
             ENDIF
             IF((IND2_Y.EQ.1).AND.(J.EQ.2)) THEN
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
                      LL=AMAXA(KI)+ABS(KJ-KI)
                   ELSE
                      LL=AMAXA(KJ)+ABS(KJ-KI)
                   ENDIF
                   A(LL)=A(LL)-S(KK1,KK2)
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
    SUBROUTINE MFS_LINE_STDM(NPORM,NPORN,II,JJ,KK,NEleAss,NX,NY,XM,XL,YM,YL,LENGTH,S,THIC)
!{{{
!--------------------------------------------------
!   GENERATE STIFFNESS MATRIX RESULT FROM BCs
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER :: NPORM,NPORN,II,JJ,KK,NEleAss
    REAL(8) :: NX,NY,XM,XL,YM,YL,LENGTH
    REAL(8) :: S(2,2),THIC
!   ASSISTANT VARIABLES
    INTEGER :: I,J
    REAL(8) :: XX,YY,DIV
    REAL(8) :: HI,HJ,HIX,HIY,HJX,HJY
    REAL(8) :: HIHJ(4)
!
    DO I=1,2
       DO J=1,2
          S(J,I)=0.
       ENDDO
    ENDDO
!
    DO I=1,4
       HIHJ(I)=0.
    ENDDO 
!    
    DO I=1,NITP_BC
       XX = XM+XL*XABSC_WEIG_BC(I,1)
       YY = YM+YL*XABSC_WEIG_BC(I,1)
       DIV= XABSC_WEIG_BC(I,2)
!
       HI = QUADS_MFS_HIHJ(NEleAss,II,NPORM,XX,YY,1)
       HIX= QUADS_MFS_HIHJ(NEleAss,II,NPORM,XX,YY,2)
       HIY= QUADS_MFS_HIHJ(NEleAss,II,NPORM,XX,YY,3)
!
       HJ = QUADS_MFS_HIHJ(NEleAss,JJ,NPORN,XX,YY,1)
       HJX= QUADS_MFS_HIHJ(NEleAss,JJ,NPORN,XX,YY,2)
       HJY= QUADS_MFS_HIHJ(NEleAss,JJ,NPORN,XX,YY,3)
!
       HIHJ(1)=HIHJ(1)+HI*HJX*DIV
       HIHJ(2)=HIHJ(2)+HI*HJY*DIV
       HIHJ(3)=HIHJ(3)+HIX*HJ*DIV
       HIHJ(4)=HIHJ(4)+HIY*HJ*DIV
    ENDDO
!
    S(1,1)=(HIHJ(1)+HIHJ(3))*NX*KCS(1)+(HIHJ(2)+HIHJ(4))*NY*KCS(3)
    S(1,2)=(HIHJ(2)*NX+HIHJ(3)*NY)*KCS(2)+(HIHJ(1)*NY+HIHJ(4)*NX)*KCS(3)
    S(2,1)=(HIHJ(1)*NY+HIHJ(4)*NX)*KCS(2)+(HIHJ(2)*NX+HIHJ(3)*NY)*KCS(3)
    S(2,2)=(HIHJ(2)+HIHJ(4))*NY*KCS(1)+(HIHJ(1)+HIHJ(3))*NX*KCS(3)
!
    S(1,1)=S(1,1)*LENGTH/2.0*THIC
    S(1,2)=S(1,2)*LENGTH/2.0*THIC
    S(2,1)=S(2,1)*LENGTH/2.0*THIC
    S(2,2)=S(2,2)*LENGTH/2.0*THIC      
!
    RETURN
!}}}
    END SUBROUTINE MFS_LINE_STDM
!
!**************************************************   version 1
!
    SUBROUTINE QUADS_MFS_STIFF(A,AMAXA,II,JJ,NumEle,NEleAss,AELEPAR,ANDMFS,ALM,THIC)
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
    REAL(8)            :: A(NWK)
    INTEGER            :: AMAXA(NEQ+1)
    INTEGER            :: NumEle,NEleAss
    INTEGER            :: ANDMFS(5,NumEle),ALM(NEleAss,NumEle)
    REAL(8)            :: AELEPAR(3,NumEle)
    REAL(8)            :: THIC
!   ASSISTANT VARIABLES
    INTEGER            :: I,J,K
    INTEGER            :: N,M
    INTEGER            :: IND1,IND2,NPOR1,NPOR2
    INTEGER            :: IND1_X,IND1_Y,IND2_X,IND2_Y
    REAL(8)            :: S(2,2)
!
    N       =  ANDMFS(1,II)
    M       =  ANDMFS(1,JJ)
!
    IND1_X  =  NDMFS(2,II)
    IND1_Y  =  NDMFS(3,II)
    IND2_X  =  NDMFS(2,JJ)
    IND2_Y  =  NDMFS(3,JJ)
!
    DO I=1,2
       DO J=1,2
          S(J,I)=0.
       ENDDO
    ENDDO
!
    IF (IND1_X.EQ.1) THEN
        NPOR1=NDOF(1,N)+1
    ELSE
        NPOR1=NDOF(1,N)
    ENDIF
!
    IF (IND2_X.EQ.1) THEN
        NPOR2=NDOF(1,M)+1
    ELSE
        NPOR2=NDOF(1,M)
    ENDIF
!
    DO I=1,NPOR1   
       DO J=1,NPOR2          
          CALL QUADS_MFS_STDM(I,J,II,JJ,NumEle,NEleAss,AELEPAR,ANDMFS,ALM,XABSC_WEIG_IN,S,THIC)     
!   ASSEMBLY
          CALL QUADS_MFS_ADDBAN(A,AMAXA,II,JJ,I,J,S,IND1_X,IND1_Y,IND2_X,IND2_Y)
       ENDDO
    ENDDO     
!
    RETURN
!!}}}
    END SUBROUTINE QUADS_MFS_STIFF
!
    SUBROUTINE QUADS_MFS_STDM (NPORM,NPORN,II,JJ,NumEle,NEleAss,AELEPAR,ANDMFS,ALM,AXABSC_WEIG_IN,S,THIC)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE ENTRY OF THE STIFFNESS MATRIX                         .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NPORM,NPORN,II,JJ,NumEle,NEleAss
    REAL(8),INTENT(IN) :: AELEPAR(3,NumEle)
    INTEGER,INTENT(IN) :: ANDMFS(5,NumEle),ALM(NEleAss,NumEle)
    REAL(8),INTENT(IN) :: AXABSC_WEIG_IN(NITP_IN,2)
    REAL(8),INTENT(INOUT) :: S(2,2) 
    REAL(8),INTENT(IN) :: THIC
!   ASSISTANT VARIABLES
    INTEGER :: I,J
    INTEGER :: N
    REAL(8) :: LENGTH,XM,XL,YM,YL
    REAL(8) :: XTMP
    REAL(8) :: HIHJH(4),HIHJG(4)
!
    DO I=1,2
       DO J=1,2
          S(J,I)=0.
       ENDDO
    ENDDO
!
    DO I=1,4
       HIHJH(I)=0.0D0
    ENDDO
!
    N=NDMFS(1,II)
    LENGTH=AELEPAR(3,II)
!
    XM=X(N)
    YM=Y(N)
    XL=LENGTH
    YL=LENGTH
!
    DO I=1,NITP_IN
       XTMP=XM+XL*AXABSC_WEIG_IN(I,1)
       CALL QUADS_MFS_GAUSSG(II,JJ,NPORM,NPORN,NEleAss,XTMP,YM,YL,HIHJG(1),AXABSC_WEIG_IN(1,1))
       HIHJH(1)=HIHJH(1)+HIHJG(1)*AXABSC_WEIG_IN(I,2)
       HIHJH(2)=HIHJH(2)+HIHJG(2)*AXABSC_WEIG_IN(I,2)
       HIHJH(3)=HIHJH(3)+HIHJG(3)*AXABSC_WEIG_IN(I,2)
       HIHJH(4)=HIHJH(4)+HIHJG(4)*AXABSC_WEIG_IN(I,2)
    ENDDO
!
    HIHJH(1)=HIHJH(1)*XL
    HIHJH(2)=HIHJH(2)*XL
    HIHJH(3)=HIHJH(3)*XL
    HIHJH(4)=HIHJH(4)*XL
!
    S(1,1)=(KCS(1)*HIHJH(1)+KCS(3)*HIHJH(4))*THIC
    S(1,2)=(KCS(2)*HIHJH(2)+KCS(3)*HIHJH(3))*THIC
    S(2,1)=(KCS(2)*HIHJH(3)+KCS(3)*HIHJH(2))*THIC
    S(2,2)=(KCS(1)*HIHJH(4)+KCS(3)*HIHJH(1))*THIC
!
    RETURN
!}}}
    END SUBROUTINE QUADS_MFS_STDM
!
    SUBROUTINE QUADS_MFS_GAUSSG(II,JJ,NPORM,NPORN,NEleAss,XTMP,YM,YL,HIHJG,AXABSC_WEIG_IN)
!{{{
    IMPLICIT NONE
    INTEGER :: II,JJ,NPORM,NPORN,NEleAss
    REAL(8) :: XTMP,YM,YL
    REAL(8) :: HIHJG(4)
    REAL(8),INTENT(IN) :: AXABSC_WEIG_IN(NITP_IN,2)
!   ASSISTANT VARIABLES
    INTEGER :: I
    INTEGER :: N,M,IND1,IND2
    INTEGER :: L1,L2,N1,N2
    REAL(8) :: LENGTHM
    REAL(8) :: YTMP
    REAL(8) :: HIX,HIY,HJX,HJY
    REAL(8) :: S1,S2,X1,X2,Y1,Y2,XX1,YY1
    LOGICAL :: REST1,REST2,REST3,LOGIC1,LOGIC2
!
    N=NDMFS(1,II)
    M=NDMFS(1,JJ)
    LENGTHM=ELEPAR(3,JJ)
!
    IND1=NDMFS(2,II)
    IND2=NDMFS(2,JJ)
!
    DO I=1,4
       HIHJG(I)=0.0D0
    ENDDO
!
    DO I=1,NITP_IN       
       YTMP=YM+YL*AXABSC_WEIG_IN(I,1)
       XX1=XTMP-X(N)
       YY1=YTMP-Y(N)
       S1=XX1*XX1+YY1*YY1
       S1=SQRT(S1)/YL
       S2=(XTMP-X(M))*(XTMP-X(M))+(YTMP-Y(M))*(YTMP-Y(M))
       S2=SQRT(S2)/LENGTHM     
       IF((S1.GT.1).OR.(S2.GT.1)) THEN
          CONTINUE
       ELSEIF ((IND1.NE.0).AND.(IND2.NE.0)) THEN
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
          IF ((LOGIC1).OR.(LOGIC2)) THEN
              CONTINUE
          ELSE
             HIX =QUADS_MFS_HIHJ(NEleAss,II,NPORM,XTMP,YTMP,2)
             HIY =QUADS_MFS_HIHJ(NEleAss,II,NPORM,XTMP,YTMP,3)
             HJX =QUADS_MFS_HIHJ(NEleAss,JJ,NPORN,XTMP,YTMP,2)
             HJY =QUADS_MFS_HIHJ(NEleAss,JJ,NPORN,XTMP,YTMP,3)
!
             HIHJG(1)=HIHJG(1)+HIX*HJX*AXABSC_WEIG_IN(I,2)
             HIHJG(2)=HIHJG(2)+HIX*HJY*AXABSC_WEIG_IN(I,2)
             HIHJG(3)=HIHJG(3)+HIY*HJX*AXABSC_WEIG_IN(I,2)
             HIHJG(4)=HIHJG(4)+HIY*HJY*AXABSC_WEIG_IN(I,2)
          ENDIF
       ELSE
          HIX =QUADS_MFS_HIHJ(NEleAss,II,NPORM,XTMP,YTMP,2)
          HIY =QUADS_MFS_HIHJ(NEleAss,II,NPORM,XTMP,YTMP,3)
          HJX =QUADS_MFS_HIHJ(NEleAss,JJ,NPORN,XTMP,YTMP,2)
          HJY =QUADS_MFS_HIHJ(NEleAss,JJ,NPORN,XTMP,YTMP,3)
!
          HIHJG(1)=HIHJG(1)+HIX*HJX*AXABSC_WEIG_IN(I,2)
          HIHJG(2)=HIHJG(2)+HIX*HJY*AXABSC_WEIG_IN(I,2)
          HIHJG(3)=HIHJG(3)+HIY*HJX*AXABSC_WEIG_IN(I,2)
          HIHJG(4)=HIHJG(4)+HIY*HJY*AXABSC_WEIG_IN(I,2)
       ENDIF
    ENDDO
    HIHJG(1)=HIHJG(1)*YL
    HIHJG(2)=HIHJG(2)*YL
    HIHJG(3)=HIHJG(3)*YL
    HIHJG(4)=HIHJG(4)*YL
!
    RETURN
!}}}
    END SUBROUTINE QUADS_MFS_GAUSSG
!
    SUBROUTINE QUADS_MFS_ADDBAN(A,AMAXA,II,JJ,NPORM,NPORN,S,IND1_X,IND1_Y,IND2_X,IND2_Y)
!{{{
    IMPLICIT NONE
    REAL(8)            :: A(NWK)
    INTEGER            :: AMAXA(NEQ+1)
    INTEGER,INTENT(IN) :: II,JJ,NPORM,NPORN
    INTEGER            :: IND1_X,IND1_Y,IND2_X,IND2_Y
    REAL(8)            :: S(2,2)
!   ASSISTANT VARIABLES
    INTEGER :: N,M
    INTEGER :: NPORMM(2),NPORNN(2)
!
    INTEGER :: I,J
    INTEGER :: KI,KJ
    INTEGER :: L
!
    N         =  NDMFS(1,II)
    M         =  NDMFS(1,JJ)
    NPORMM(1) =  NPORM
    NPORMM(2) =  NPORM
    NPORNN(1) =  NPORN
    NPORNN(2) =  NPORN
!
       IF(IND1_X.EQ.1) NPORMM(1)=NPORM-1
       IF(IND1_Y.EQ.1) NPORMM(2)=NPORM-1
       IF(IND2_X.EQ.1) NPORNN(1)=NPORN-1
       IF(IND2_Y.EQ.1) NPORNN(2)=NPORN-1
!
       IF(NPORMM(1).EQ.0) THEN
          S(1,1)=0.
          S(1,2)=0.
       ENDIF
       IF(NPORMM(2).EQ.0) THEN
          S(2,1)=0.
          S(2,2)=0.
       ENDIF
       IF(NPORNN(1).EQ.0) THEN
          S(1,1)=0.
          S(2,1)=0.
       ENDIF
       IF(NPORNN(2).EQ.0) THEN
          S(1,2)=0.
          S(2,2)=0.
       ENDIF
!
    DO I=1,2
       KI=ID(I,N)+NPORMM(I)-1
       IF(KI.GT.0) THEN
          DO J=1,2
             KJ=ID(J,M)+NPORNN(J)-1
             IF(KJ.GT.0) THEN
               IF ((N.NE.M).OR.(KJ.GE.KI)) THEN
                  IF(KJ.LE.KI) THEN
                     L=AMAXA(KI)+ABS(KJ-KI)
                  ELSE
                     L=AMAXA(KJ)+ABS(KJ-KI)
                  ENDIF
                  A(L)=A(L)+S(I,J)
               ENDIF
             ENDIF
          ENDDO
       ENDIF          
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE QUADS_MFS_ADDBAN
!
!**************************************************   version 2
!
    SUBROUTINE QUADS_MFS_STIFF_V2(A,AMAXA,II,JJ,NumEle,NEleAss,AELEPAR,ANDMFS,THIC)
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
    REAL(8)            :: A(NWK)
    INTEGER            :: AMAXA(NEQ+1)
    INTEGER            :: NumEle,NEleAss
    INTEGER            :: ANDMFS(5,NumEle)
    REAL(8)            :: AELEPAR(3,NumEle)
    REAL(8)            :: THIC
!   ASSISTANT VARIABLES
    INTEGER            :: I,J,K,L
    INTEGER            :: N,M
    INTEGER            :: NPOR1,NPOR2
    INTEGER            :: IND1_X,IND1_Y,IND2_X,IND2_Y 
!
    REAL(8)            :: XM,YM,XL,YL,LENGTH,XTMP
    REAL(8)            :: HIHJH(4,NPOR_MAX*NPOR_MAX),HIHJG(4,NPOR_MAX*NPOR_MAX)
!
    N       =  ANDMFS(1,II)
    M       =  ANDMFS(1,JJ)
!
    IND1_X  =  NDMFS(2,II)
    IND1_Y  =  NDMFS(3,II)
    IND2_X  =  NDMFS(2,JJ)
    IND2_Y  =  NDMFS(3,JJ)
!
    IF (IND1_X.EQ.1) THEN
        NPOR1=NDOF(1,N)+1
    ELSE
        NPOR1=NDOF(1,N)
    ENDIF
!
    IF (IND2_X.EQ.1) THEN
        NPOR2=NDOF(1,M)+1
    ELSE
        NPOR2=NDOF(1,M)
    ENDIF
!
    LENGTH=AELEPAR(3,II)
    XM=X(N)
    YM=Y(N)
    XL=LENGTH
    YL=LENGTH
!
    DO I=1,NPOR_MAX*NPOR_MAX
       STIFF_IJ(1,I)=0.0D0
       STIFF_IJ(2,I)=0.0D0
       STIFF_IJ(3,I)=0.0D0
       STIFF_IJ(4,I)=0.0D0
    ENDDO
!
    DO I=1,NPOR_MAX*NPOR_MAX
       HIHJH(1,I)=0.0D0
       HIHJH(2,I)=0.0D0
       HIHJH(3,I)=0.0D0
       HIHJH(4,I)=0.0D0
    ENDDO
!
    DO I=1,NITP_IN
       XTMP=XM+XL*XABSC_WEIG_IN(I,1)
       CALL QUADS_MFS_GAUSSG_V2(II,JJ,NPOR1,NPOR2,NEleAss,XTMP,YM,YL,HIHJG(1,1))
       DO J=1,NPOR1
          DO K=1,NPOR2
             L=(J-1)*NPOR2+K
             HIHJH(1,L)=HIHJH(1,L)+HIHJG(1,L)*XABSC_WEIG_IN(I,2)
             HIHJH(2,L)=HIHJH(2,L)+HIHJG(2,L)*XABSC_WEIG_IN(I,2)
             HIHJH(3,L)=HIHJH(3,L)+HIHJG(3,L)*XABSC_WEIG_IN(I,2)
             HIHJH(4,L)=HIHJH(4,L)+HIHJG(4,L)*XABSC_WEIG_IN(I,2)
          ENDDO
       ENDDO
    ENDDO
!
    DO I=1,NPOR1
       DO J=1,NPOR2
          K=(I-1)*NPOR2+J
          STIFF_IJ(1,K)=(KCS(1)*HIHJH(1,K)+KCS(3)*HIHJH(4,K))*THIC*XL
          STIFF_IJ(2,K)=(KCS(2)*HIHJH(2,K)+KCS(3)*HIHJH(3,K))*THIC*XL
          STIFF_IJ(3,K)=(KCS(2)*HIHJH(3,K)+KCS(3)*HIHJH(2,K))*THIC*XL
          STIFF_IJ(4,K)=(KCS(1)*HIHJH(4,K)+KCS(3)*HIHJH(1,K))*THIC*XL
!          WRITE(*,*) STIFF_IJ(1,K),STIFF_IJ(2,K),STIFF_IJ(2,K),STIFF_IJ(3,K)
       ENDDO
    ENDDO 
!
    CALL QUADS_MFS_ADDBAN_V2(A(1),AMAXA(1),II,JJ,NPOR1,NPOR2,IND1_X,IND1_Y,IND2_X,IND2_Y,STIFF_IJ(1,1))
!   
    RETURN
!!}}}
    END SUBROUTINE QUADS_MFS_STIFF_V2
!
    SUBROUTINE QUADS_MFS_GAUSSG_V2(II,JJ,NPOR1,NPOR2,NEleAss,XTMP,YM,YL,HIHJG)
!{{{
    IMPLICIT NONE
    INTEGER :: II,JJ,NPOR1,NPOR2,NEleAss
    REAL(8) :: XTMP,YM,YL
    REAL(8) :: HIHJG(4,NPOR1*NPOR2)
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L
    INTEGER :: N,M,IND1,IND2
    INTEGER :: L1,L2,N1,N2
    REAL(8) :: LENGTHM
    REAL(8) :: YTMP
    REAL(8) :: HI(NPOR1,3),HJ(NPOR2,3)
    REAL(8) :: S1,S2,X1,X2,Y1,Y2,XX1,YY1
    LOGICAL :: REST1,REST2,REST3,LOGIC1,LOGIC2
!
    N=NDMFS(1,II)
    M=NDMFS(1,JJ)
    LENGTHM=ELEPAR(3,JJ)
!
    IND1=NDMFS(2,II)
    IND2=NDMFS(2,JJ)
!
    DO I=1,NPOR1*NPOR2
       DO J=1,4
         HIHJG(J,I)=0.0D0
         HIHJG(J,I)=0.0D0
       ENDDO
    ENDDO
!
    DO I=1,NITP_IN       
       YTMP=YM+YL*XABSC_WEIG_IN(I,1)
       XX1=XTMP-X(N)
       YY1=YTMP-Y(N)
       S1=XX1*XX1+YY1*YY1
       S1=SQRT(S1)/YL
       S2=(XTMP-X(M))*(XTMP-X(M))+(YTMP-Y(M))*(YTMP-Y(M))
       S2=SQRT(S2)/LENGTHM     
       IF((S1.GT.1).OR.(S2.GT.1)) THEN
          CONTINUE
       ELSEIF ((IND1.NE.0).AND.(IND2.NE.0)) THEN
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
          IF ((LOGIC1).OR.(LOGIC2)) THEN
              CONTINUE
          ELSE
             HI =QUADS_MFS_HIHJ_V2(NEleAss,II,XTMP,YTMP,NPOR1)
             HJ =QUADS_MFS_HIHJ_V2(NEleAss,JJ,XTMP,YTMP,NPOR2)
!
             DO J=1,NPOR1
                DO K=1,NPOR2
                   L=(J-1)*NPOR2+K
                   HIHJG(1,L)=HIHJG(1,L)+HI(J,2)*HJ(K,2)*XABSC_WEIG_IN(I,2)
                   HIHJG(2,L)=HIHJG(2,L)+HI(J,2)*HJ(K,3)*XABSC_WEIG_IN(I,2)
                   HIHJG(3,L)=HIHJG(3,L)+HI(J,3)*HJ(K,2)*XABSC_WEIG_IN(I,2)
                   HIHJG(4,L)=HIHJG(4,L)+HI(J,3)*HJ(K,3)*XABSC_WEIG_IN(I,2)
                ENDDO             
             ENDDO                
          ENDIF
       ELSE
          HI =QUADS_MFS_HIHJ_V2(NEleAss,II,XTMP,YTMP,NPOR1)
          HJ =QUADS_MFS_HIHJ_V2(NEleAss,JJ,XTMP,YTMP,NPOR2)
!
          DO J=1,NPOR1
             DO K=1,NPOR2
                L=(J-1)*NPOR2+K
                HIHJG(1,L)=HIHJG(1,L)+HI(J,2)*HJ(K,2)*XABSC_WEIG_IN(I,2)
                HIHJG(2,L)=HIHJG(2,L)+HI(J,2)*HJ(K,3)*XABSC_WEIG_IN(I,2)
                HIHJG(3,L)=HIHJG(3,L)+HI(J,3)*HJ(K,2)*XABSC_WEIG_IN(I,2)
                HIHJG(4,L)=HIHJG(4,L)+HI(J,3)*HJ(K,3)*XABSC_WEIG_IN(I,2)
             ENDDO             
          ENDDO 
       ENDIF
    ENDDO
!
    DO I=1,NPOR1*NPOR2
       HIHJG(1,I)=HIHJG(1,I)*YL
       HIHJG(2,I)=HIHJG(2,I)*YL
       HIHJG(3,I)=HIHJG(3,I)*YL
       HIHJG(4,I)=HIHJG(4,I)*YL
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE QUADS_MFS_GAUSSG_V2
!
    SUBROUTINE QUADS_MFS_ADDBAN_V2(A,AMAXA,II,JJ,NPOR1,NPOR2,IND1_X,IND1_Y,IND2_X,IND2_Y,S)
!{{{
    IMPLICIT NONE
    REAL(8)            :: A(NWK)
    INTEGER            :: AMAXA(NEQ+1)
    INTEGER,INTENT(IN) :: II,JJ,NPOR1,NPOR2
    INTEGER            :: IND1_X,IND1_Y,IND2_X,IND2_Y
    REAL(8)            :: S(4,NPOR1*NPOR2)
!   ASSISTANT VARIABLES
    INTEGER :: N,M
!
    INTEGER :: I,J,K,L
    INTEGER :: KI,KJ
    INTEGER :: LL,KK1,KK2
    INTEGER :: INDK,INDL
!
    N         =  NDMFS(1,II)
    M         =  NDMFS(1,JJ)
!
    DO I=1,2
       DO J=1,2
          K=1
          INDK=1
          IF((IND1_X.EQ.1).AND.(I.EQ.1)) THEN
             K=2
             INDK=2
          ENDIF
          IF((IND1_Y.EQ.1).AND.(I.EQ.2)) THEN
             K=2
             INDK=2
          ENDIF          
          DO WHILE(K.LE.NPOR1)
             L=1
             INDL=1
             IF((IND2_X.EQ.1).AND.(J.EQ.1)) THEN
                 L=2
                 INDL=2
             ENDIF
             IF((IND2_Y.EQ.1).AND.(J.EQ.2)) THEN
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
                      LL=AMAXA(KI)+ABS(KJ-KI)
                   ELSE
                      LL=AMAXA(KJ)+ABS(KJ-KI)
                   ENDIF
                   A(LL)=A(LL)+S(KK1,KK2)
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
    SUBROUTINE ADDBAN (A,AMAXA,S,NND,LMM,ND)
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
      REAL(8)    :: A(NWK), S(NND)
      INTEGER    :: LMM(ND),AMAXA(NEQ+1)
      INTEGER    :: NDI, I,II, MI, KS, J, JJ, IJ, KK, KSS
!
      NDI=0
      DO I=1,ND
         II=LMM(I)
         IF (II .GT. 0) THEN
            MI=AMAXA(II)
            KS=I
            DO J=1,ND
               JJ=LMM(J)
               IF (JJ .GT. 0) THEN
                  IJ=II - JJ
                  IF (IJ .GE. 0) THEN
                     KK=MI + IJ
                     KSS=KS
                     IF (J.GE.I) KSS=J + NDI
                     A(KK)=A(KK) + S(KSS)
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
    SUBROUTINE OUTPUTKM(ADK,ADM,OKMATRIX,OMMATRIX,OMAXA)
!{{{
!--------------------------------------------------
!   OUTPUT STIFFNESS AND MASS MATRIX
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADK,ADM
    INTEGER :: I,J,K
    REAL(8) :: KI(NEQ),MI(NEQ)
!
    REAL(8),DIMENSION(NWK),   INTENT(IN) :: OKMATRIX,OMMATRIX
    INTEGER,DIMENSION(NEQ+1), INTENT(IN) :: OMAXA
!--- ASSISTANT VARIABLES
    INTEGER :: II
!
    DO I=1,NEQ
       DO J=1,NEQ
          IF(I.LE.J) THEN
             II=OMAXA(J)+J-I
          ELSE
             II=OMAXA(I)+I-J
          END IF
             KI(J)=OKMATRIX(II)
             MI(J)=OMMATRIX(II)
        END DO
        WRITE(ADK,1000) (KI(K),K=1,NEQ)
        WRITE(ADM,1000) (MI(K),K=1,NEQ)
     END DO
     RETURN
!
1000 FORMAT(*(ES16.9,1X))
!}}}
    END SUBROUTINE OUTPUTKM
!      
    SUBROUTINE ADDBAN_ZLB(ND,S,LMM,A,AMAXA)
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
   REAL(8),DIMENSION(NWK),INTENT(INOUT) :: A
   INTEGER,DIMENSION(NEQ+1),INTENT(IN) :: AMAXA
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
                K=AMAXA(M)+K
                A(K)=A(K)+S(I,J)
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
!**************************************************  FUNCTIONS
!
    FUNCTION QUADS_MFS_HIHJ(NEleAss,II,NPOR,XX,YY,IND) RESULT(H)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE VALUE OF HI,HJ,DHI,DHJ                                .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NEleAss,II,NPOR,IND
    REAL(8),INTENT(IN) :: XX,YY 
    REAL(8)            :: H
!   ASSISTANT VARIABLES
    INTEGER :: I,J
    INTEGER :: N,M
    REAL(8) :: W1,W2,DW1DX,DW1DY,DW2DX,DW2DY
    REAL(8) :: H0,HX,HY,HXY,HXX,HYY
    REAL(8) :: DH0DX,DHXDX,DHYDX,DHXYDX,DHXXDX,DHYYDX
    REAL(8) :: DH0DY,DHXDY,DHYDY,DHXYDY,DHXXDY,DHYYDY
    REAL(8) :: Radius1,Radius
    REAL(8) :: S1,S
!
    Radius1=ELEPAR(3,II)
    N=NDMFS(1,II)
    S1=(XX-X(N))*(XX-X(N))+(YY-Y(N))*(YY-Y(N))
    S1=SQRT(S1)/Radius1
!   CALCULATE W1,W2,DW1,DW2
    W1  = 1-6*S1*S1+8*S1*S1*S1-3*S1*S1*S1*S1
!
    W2=0.0D0
    DO I=1,NEleAss
       J=LM(I,II)
       IF (J.NE.0) THEN
          M=NDMFS(1,J)
          Radius=ELEPAR(3,J)
          S=(XX-X(M))*(XX-X(M))+(YY-Y(M))*(YY-Y(M))
          S=SQRT(S)/Radius
          IF(S.LE.1) THEN
             W2=W2+1-6*S*S+8*S*S*S-3*S*S*S*S
          ENDIF
       ENDIF
    ENDDO
!
    SELECT CASE(IND)
    CASE(1)
       SELECT CASE(NPOR)
       CASE(1)
           H0=W1/W2
           H=H0
       CASE(2)
           HX=W1/W2*(XX-X(N))/Radius1
           H=HX
       CASE(3)
           HY=W1/W2*(YY-Y(N))/Radius1
           H=HY
       CASE(4)
           HXY=W1/W2*(XX-X(N))*(YY-Y(N))/Radius1/Radius1
           H=HXY
       CASE(5)
           HXX=W1/W2*(XX-X(N))*(XX-X(N))/Radius1/Radius1
           H=HXX
       CASE(6)
           HYY=W1/W2*(YY-Y(N))*(YY-Y(N))/Radius1/Radius1
           H=HYY
       CASE DEFAULT
          WRITE(*,*) "ERROR, NPOR SHOULD LT 4!"
          STOP
       END SELECT        
    CASE(2)
       DW1DX=(-12*S1+24*S1*S1-12*S1*S1*S1)*(XX-X(N))/S1/Radius1/Radius1
       DW2DX=0.0D0
       DO I=1,NEleAss
          J=LM(I,II)
          IF (J.NE.0) THEN
             M=NDMFS(1,J)
             Radius=ELEPAR(3,J)
             S=(XX-X(M))*(XX-X(M))+(YY-Y(M))*(YY-Y(M))
             S=SQRT(S)/Radius
             IF(S.LE.1) THEN
                DW2DX=DW2DX+(-12*S+24*S*S-12*S*S*S)*(XX-X(M))/S/Radius/Radius
             ENDIF
          ENDIF
       ENDDO
!
       SELECT CASE(NPOR)
       CASE(1)
           DH0DX=(DW1DX*W2-W1*DW2DX)/W2/W2
           H=DH0DX
       CASE(2)
           DH0DX=(DW1DX*W2-W1*DW2DX)/W2/W2
           H0=W1/W2
           DHXDX=DH0DX*(XX-X(N))/Radius1+H0/Radius1
           H=DHXDX
       CASE(3)
           DH0DX=(DW1DX*W2-W1*DW2DX)/W2/W2
           DHYDX=DH0DX*(YY-Y(N))/Radius1
           H=DHYDX
       CASE(4)               
           DH0DX = (DW1DX*W2-W1*DW2DX)/W2/W2
           H0    = W1/W2           
           DHXYDX= DH0DX*(XX-X(N))*(YY-Y(N))/Radius1/Radius1+H0*(YY-Y(N))/Radius1/Radius1
           H     = DHXYDX
       CASE(5)
           DH0DX = (DW1DX*W2-W1*DW2DX)/W2/W2
           H0    = W1/W2           
           DHXXDX= DH0DX*(XX-X(N))*(XX-X(N))/Radius1/Radius1+2*H0*(XX-X(N))/Radius1/Radius1
           H     = DHXXDX
       CASE(6)
           DH0DX = (DW1DX*W2-W1*DW2DX)/W2/W2         
           DHYYDX= DH0DX*(YY-Y(N))*(YY-Y(N))/Radius1/Radius1
           H     = DHYYDX
       CASE DEFAULT
          WRITE(*,*) "ERROR, NPOR SHOULD LT 4!"
          STOP
       END SELECT
    CASE(3)
       DW1DY=(-12*S1+24*S1*S1-12*S1*S1*S1)*(YY-Y(N))/S1/Radius1/Radius1
       DW2DY=0
       DO I=1,NEleAss
          J=LM(I,II)
          IF (J.NE.0) THEN
             M=NDMFS(1,J)
             Radius=ELEPAR(3,J)
             S=(XX-X(M))*(XX-X(M))+(YY-Y(M))*(YY-Y(M))
             S=SQRT(S)/Radius
             IF(S.LE.1) THEN
                DW2DY=DW2DY+(-12*S+24*S*S-12*S*S*S)*(YY-Y(M))/S/Radius/Radius                
             ENDIF
          ENDIF
       ENDDO
!
       SELECT CASE(NPOR)
       CASE(1)
           DH0DY=(DW1DY*W2-W1*DW2DY)/W2/W2
           H=DH0DY
       CASE(2)
           DH0DY=(DW1DY*W2-W1*DW2DY)/W2/W2
           DHXDY=DH0DY*(XX-X(N))/Radius1
           H=DHXDY
       CASE(3)
           DH0DY=(DW1DY*W2-W1*DW2DY)/W2/W2
           H0=W1/W2
           DHYDY=DH0DY*(YY-Y(N))/Radius1+H0/Radius1
           H=DHYDY
       CASE(4)
           DH0DY=(DW1DY*W2-W1*DW2DY)/W2/W2
           H0    = W1/W2           
           DHXYDY= DH0DY*(XX-X(N))*(YY-Y(N))/Radius1/Radius1+H0*(XX-X(N))/Radius1/Radius1
           H     = DHXYDY
       CASE(5)  
           DH0DY =(DW1DY*W2-W1*DW2DY)/W2/W2
           DHXXDY= DH0DY*(XX-X(N))*(XX-X(N))/Radius1/Radius1  
           H     = DHXXDY    
       CASE(6)
           DH0DY =(DW1DY*W2-W1*DW2DY)/W2/W2
           H0    = W1/W2 
           DHYYDY= DH0DY*(YY-Y(N))*(YY-Y(N))/Radius1/Radius1+2*H0*(YY-Y(N))/Radius1/Radius1             
           H     = DHYYDY               
       CASE DEFAULT
          WRITE(*,*) "ERROR, NPOR SHOULD LT 4!"
          STOP
       END SELECT
    CASE DEFAULT
       WRITE(*,*) "ERROR, IND SHOULD LT 4!"
       STOP
    END SELECT
!
    RETURN
!}}}
    END FUNCTION QUADS_MFS_HIHJ
!
    FUNCTION QUADS_MFS_HIHJ_V2(NEleAss,II,XX,YY,NPORMAX) RESULT(H)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE VALUE OF HI,HJ,DHI,DHJ                                .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NEleAss,II,NPORMAX
    REAL(8),INTENT(IN) :: XX,YY 
    REAL(8)            :: H(NPORMAX,3)
!   ASSISTANT VARIABLES
    INTEGER :: I,J
    INTEGER :: N,M
    REAL(8) :: W1,W2,DW1DX,DW1DY,DW2DX,DW2DY
    REAL(8) :: H0
    REAL(8) :: Radius1,Radius
    REAL(8) :: S1,S,SS1(4),SS(4)
!
    Radius1=ELEPAR(3,II)
    N=NDMFS(1,II)
    S1=(XX-X(N))*(XX-X(N))+(YY-Y(N))*(YY-Y(N))
    S1=SQRT(S1)/Radius1
    SS1(1)=S1*S1
    SS1(2)=S1*S1*S1
    SS1(3)=S1*S1*S1*S1
    SS1(4)=-12*S1+24*SS1(1)-12*SS1(2)
!   CALCULATE W1,W2,DW1,DW2
    W1  = 1-6*SS1(1)+8*SS1(2)-3*SS1(3)
    DW1DX=SS1(4)*(XX-X(N))/S1/Radius1/Radius1
    DW1DY=SS1(4)*(YY-Y(N))/S1/Radius1/Radius1
!
    W2=0.0D0
    DW2DX=0.0D0
    DW2DY=0.0D0
    DO I=1,NEleAss
       J=LM(I,II)
       IF (J.NE.0) THEN
          M=NDMFS(1,J)
          Radius=ELEPAR(3,J)
          S=(XX-X(M))*(XX-X(M))+(YY-Y(M))*(YY-Y(M))
          S=SQRT(S)/Radius
!
          SS(1)=S*S
          SS(2)=S*S*S
          SS(3)=S*S*S*S
          SS(4)=-12*S+24*SS(1)-12*SS(2)
!
          IF(S.LE.1) THEN
             W2=W2+1-6*SS(1)+8*SS(2)-3*SS(3)
             DW2DX=DW2DX+SS(4)*(XX-X(M))/S/Radius/Radius
             DW2DY=DW2DY+SS(4)*(YY-Y(M))/S/Radius/Radius 
          ENDIF
       ENDIF
    ENDDO
!   HI,DHIDX,DHIDY
    IF(NPORMAX.GE.1) THEN
      H(1,1)=W1/W2
      H(1,2)=(DW1DX*W2-W1*DW2DX)/W2/W2
      H(1,3)=(DW1DY*W2-W1*DW2DY)/W2/W2
      IF(NPORMAX.GE.3) THEN
!   HIX,DHIXDX,DHIXDY
        H(2,1)=H(1,1)*(XX-X(N))/Radius1
        H(2,2)=H(1,2)*(XX-X(N))/Radius1+H(1,1)/Radius1
        H(2,3)=H(1,3)*(XX-X(N))/Radius1
!   HIY,DHIYDX,DHIYDY
        H(3,1)=H(1,1)*(YY-Y(N))/Radius1
        H(3,2)=H(1,2)*(YY-Y(N))/Radius1
        H(3,3)=H(1,3)*(YY-Y(N))/Radius1+H(1,1)/Radius1
        IF(NPORMAX.GE.4) THEN
!   HIXY,DHIXYDX,DHIXYDY
           H(4,1)=H(2,1)*(YY-Y(N))/Radius1
           H(4,2)=(H(3,2)*(XX-X(N))+H(3,1))/Radius1
           H(4,3)=(H(2,3)*(YY-Y(N))+H(2,1))/Radius1
           IF(NPORMAX.GE.6) THEN
!   HIXX,DHIXXDX,DHIXXDY
              H(5,1)=H(2,1)*(XX-X(N))/Radius1
              H(5,2)=(H(2,2)*(XX-X(N))+H(2,1))/Radius1
              H(5,3)=H(2,3)*(XX-X(N))/Radius1
!   HIYY,DHIYYDX,DHIYYDY
              H(6,1)=H(3,1)*(YY-Y(N))/Radius1
              H(6,2)=H(3,2)*(YY-Y(N))/Radius1
              H(6,3)=(H(3,3)*(YY-Y(N))+H(3,1))/Radius1
              IF(NPORMAX.GT.6) THEN
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
    END FUNCTION QUADS_MFS_HIHJ_V2
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

end module GENERATEM
