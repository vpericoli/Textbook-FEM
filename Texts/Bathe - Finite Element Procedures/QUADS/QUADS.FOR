      SUBROUTINE QUADS (NEL,ITYPE,NINT,THIC,YM,PR,XX,S,IOUT)            
C                                                                       
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C .                                                                   . 
C .   P R O G R A M                                                   . 
C .        TO CALCULATE ISOPARAMETRIC QUADRILATERAL ELEMENT STIFFNESS . 
C .        MATRIX FOR AXISYMMETRIC, PLANE STRESS, AND PLANE STRAIN    . 
C .        CONDITIONS                                                 . 
C .                                                                   . 
C .  - - INPUT VARIABLES - -                                          . 
C .        NEL       = NUMBER OF ELEMENT                              . 
C .        ITYPE     = ELEMENT TYPE                                   . 
C .                        EQ.0 = AXISYMMETRIC                        . 
C .                        EQ.1 = PLANE STRAIN                        . 
C .                        EQ.2 = PLANE STRESS                        . 
C .        NINT      = GAUSS NUMERICAL INTEGRATION ORDER              . 
C .        THIC      = THICKNESS OF ELEMENT                           . 
C .        YM        = YOUNG'S MODULUS                                . 
C .        PR        = POISSON'S RATIO                                . 
C .        XX(2,4)   = ELEMENT NODE COORDINATES                       . 
C .        S(8,8)    = STORAGE FOR STIFFNESS MATRIX                   . 
C .        IOUT      = UNIT NUMBER USED FOR OUTPUT                    . 
C .                                                                   . 
C .  - - OUTPUT - -                                                   . 
C .        S(8,8)    = CALCULATED STIFFNESS MATRIX                    . 
C .                                                                   . 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C .   THIS PROGRAM IS USED IN SINGLE PRECISION ARITHMETIC ON CRAY     . 
C .   EQUIPMENT AND DOUBLE PRECISION ARITHMETIC ON IBM MACHINES,      . 
C .   ENGINEERING WORKSTATIONS AND PCS. DEACTIVATE ABOVE LINE FOR     . 
C .   SINGLE PRECISION ARITHMETIC.                                    . 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      DIMENSION D(4,4),B(4,8),XX(2,4),S(8,8),XG(4,4),WGT(4,4),DB(4)     
C                                                                       
C     MATRIX XG STORES GAUSS - LEGENDRE SAMPLING POINTS                 
C                                                                       
      DATA XG/   0.D0,   0.D0,   0.D0,   0.D0,   -.5773502691896D0,     
     1 .5773502691896D0,   0.D0,   0.D0,   -.7745966692415D0,   0.D0,   
     2 .7745966692415D0,   0.D0,   -.8611363115941D0,                   
     3 -.3399810435849D0,   .3399810435849D0,   .8611363115941D0 /      
C                                                                       
C     MATRIX WGT STORES GAUSS - LEGENDRE WEIGHTING FACTORS              
C                                                                       
      DATA WGT /  2.D0,   0.D0,   0.D0,   0.D0,   1.D0,   1.D0,         
     1 0.D0,   0.D0,   .5555555555556D0,   .8888888888889D0,            
     2 .5555555555556D0,   0.D0,   .3478548451375D0,   .6521451548625D0,
     3 .6521451548625D0,   .3478548451375D0 /                           
C                                                                       
C     O B T A I N  S T R E S S - S T R A I N  L A W                     
C                                                                       
      F=YM/(1.+PR)                                                      
      G=F*PR/(1.-2.*PR)                                                 
      H=F + G                                                           
C                                                                       
C     PLANE STRAIN ANALYSIS                                             
C                                                                       
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
C                                                                       
C     AXISYMMETRIC ANALYSIS                                             
C                                                                       
      D(1,4)=G                                                          
      D(2,4)=G                                                          
      D(3,4)=0.                                                         
      D(4,1)=G                                                          
      D(4,2)=G                                                          
      D(4,3)=0.                                                         
      D(4,4)=H                                                          
      IF (ITYPE.EQ.0) GO TO 20                                          
C                                                                       
C     FOR PLANE STRESS ANALYSIS CONDENSE STRESS-STRAIN MATRIX           
C                                                                       
      DO 10 I=1,3                                                       
      A=D(I,4)/D(4,4)                                                   
      DO 10 J=I,3                                                       
      D(I,J)=D(I,J) - D(4,J)*A                                          
   10 D(J,I)=D(I,J)                                                     
C                                                                       
C     C A L C U L A T E  E L E M E N T  S T I F F N E S S               
C                                                                       
   20 DO 30 I=1,8                                                       
      DO 30 J=1,8                                                       
   30 S(I,J)=0.                                                         
      IST=3                                                             
      IF (ITYPE.EQ.0) IST=4                                             
      DO 80 LX=1,NINT                                                   
      RI=XG(LX,NINT)                                                    
      DO 80 LY=1,NINT                                                   
      SI=XG(LY,NINT)                                                    
C                                                                       
C     EVALUATE DERIVATIVE OPERATOR B AND THE JACOBIAN DETERMINANT DET   
C                                                                       
      CALL STDM (XX,B,DET,RI,SI,XBAR,NEL,ITYPE,IOUT)                    
C                                                                       
C     ADD CONTRIBUTION TO ELEMENT STIFFNESS                             
C                                                                       
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
   80 CONTINUE                                                          
C                                                                       
      DO 90 J=1,8                                                       
      DO 90 I=J,8                                                       
   90 S(J,I)=S(I,J)                                                     
C                                                                       
      RETURN                                                            
C                                                                       
      END                                                               
      SUBROUTINE STDM (XX,B,DET,R,S,XBAR,NEL,ITYPE,IOUT)                
C                                                                       
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C .                                                                   . 
C .   P R O G R A M                                                   . 
C .     TO EVALUATE THE STRAIN-DISPLACEMENT TRANSFORMATION MATRIX B   . 
C .     AT POINT (R,S) FOR A QUADRILATERAL ELEMENT                    . 
C .                                                                   . 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION XX(2,4),B(4,8),H(4),P(2,4),XJ(2,2),XJI(2,2)             
C                                                                       
      RP = 1.0 + R                                                      
      SP = 1.0 + S                                                      
      RM = 1.0 - R                                                      
      SM = 1.0 - S                                                      
C                                                                       
C     INTERPOLATION FUNCTIONS                                           
C                                                                       
      H(1) = 0.25* RP* SP                                               
      H(2) = 0.25* RM* SP                                               
      H(3) = 0.25* RM* SM                                               
      H(4) = 0.25* RP* SM                                               
C                                                                       
C     NATURAL COORDINATE DERIVATIVES OF THE INTERPOLATION FUNCTIONS     
C                                                                       
C        1. WITH RESPECT TO R                                           
C                                                                       
      P(1,1) = 0.25* SP                                                 
      P(1,2) = - P(1,1)                                                 
      P(1,3) = - 0.25* SM                                               
      P(1,4) = - P(1,3)                                                 
C                                                                       
C        2. WITH RESPECT TO S                                           
C                                                                       
      P(2,1) = 0.25* RP                                                 
      P(2,2) = 0.25* RM                                                 
      P(2,3) = - P(2,2)                                                 
      P(2,4) = - P(2,1)                                                 
C                                                                       
C     EVALUATE THE JACOBIAN MATRIX AT POINT (R,S)                       
C                                                                       
   10 DO 30 I=1,2                                                       
      DO 30 J=1,2                                                       
      DUM = 0.0                                                         
      DO 20 K=1,4                                                       
   20 DUM=DUM + P(I,K)*XX(J,K)                                          
   30 XJ(I,J)=DUM                                                       
C                                                                       
C     COMPUTE THE DETERMINANT OF THE JACOBIAN MATRIX AT POINT (R,S)     
C                                                                       
      DET = XJ(1,1)* XJ(2,2) - XJ(2,1)* XJ(1,2)                         
      IF (DET.GT.0.00000001) GO TO 40                                   
      WRITE (IOUT,2000) NEL                                             
      GO TO 800                                                         
C                                                                       
C     COMPUTE INVERSE OF THE JACOBIAN MATRIX                            
C                                                                       
   40 DUM=1./DET                                                        
      XJI(1,1) = XJ(2,2)* DUM                                           
      XJI(1,2) =-XJ(1,2)* DUM                                           
      XJI(2,1) =-XJ(2,1)* DUM                                           
      XJI(2,2) = XJ(1,1)* DUM                                           
C                                                                       
C     EVALUATE GLOBAL DERIVATIVE OPERATOR B                             
C                                                                       
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
C                                                                       
C     IN CASE OF PLANE STRAIN OR PLANE STRESS ANALYSIS DO NOT INCLUDE   
C     THE NORMAL STRAIN COMPONENT                                       
C                                                                       
      IF (ITYPE.GT.0) GO TO 900                                         
C                                                                       
C     COMPUTE THE RADIUS AT POINT (R,S)                                 
C                                                                       
      XBAR=0.0                                                          
      DO 70 K=1,4                                                       
   70 XBAR=XBAR + H(K)*XX(1,K)                                          
C                                                                       
C     EVALUATE THE HOOP STRAIN-DISPLACEMENT RELATION                    
C                                                                       
      IF (XBAR.GT.0.00000001) GO TO 90                                  
C                                                                       
C     FOR THE CASE OF ZERO RADIUS EQUATE RADIAL TO HOOP STRAIN          
C                                                                       
      DO 80 K=1,8                                                       
   80 B(4,K)=B(1,K)                                                     
      GO TO 900                                                         
C                                                                       
C     NON-ZERO RADIUS                                                   
C                                                                       
   90 DUM=1./XBAR                                                       
      K2=0                                                              
      DO 100 K=1,4                                                      
      K2=K2 + 2                                                         
      B(4,K2  ) = 0.                                                    
  100 B(4,K2-1) = H(K)*DUM                                              
      GO TO 900                                                         
C                                                                       
  800 STOP                                                              
  900 RETURN                                                            
C                                                                       
 2000 FORMAT (//,' *** ERROR *** ',                                     
     1    ' ZERO OR NEGATIVE JACOBIAN DETERMINANT FOR ELEMENT (',I8,')')
C                                                                       
      END                                                               
