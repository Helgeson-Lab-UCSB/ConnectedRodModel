      SUBROUTINE TEFIT(Y,X,BB,RES,NARRAY,IBB,PARMFILE)
CTEFIT      NLLSQ - NON LINEAR LEAST SQUARES FITTING PROGRAM
C     MODIFIED TO GIVE THE OPTION OF A TAYLOR EXPANSION FIT ALONE
C     NONLINEAR LEAST SQUARES FITTING ALGORITHM BY D W MARQUARD
C     ORIGINAL PROGRAM REWRITTEN BY W A BURNETTE  BTL JULY  1967
C     MODIFIED TO GIVE THE OPTION OF A TAYLOR EXPANSION FIT ALONE
C     NARRAY CONTAINS PROGRAM PARAMETERS  ARRAY CONTAINS STATISTICAL
C     CONSTANTS   SET ARRAY EQUAL TO 0.FOR STANDARD SET OF CONSTANTS
C     MAXIMUM NUMBER OF PARAMETERS IS 20 THIS MAY BE CHANGED BY ALTERING
C     DIMENSION STATEMENTS AND MATRIX STORING STATEMENTS


C...	NOTES ON VARIABLES:
C...	   NARRAY(1) = # OF DATA POINTS =N
C...	   NARRAY(2) = # OF INDEPENDENT VARIABLES =M
C...	   NARRAY(3) = TOTAL # OF PARAMETERS (FIXED & VARIED)=K
C...	   NARRAY(4) = # OF PARAMETERS FIXED =IP
C...	   NARRAY(5) = INTERMEDIATE PRINTOUT SPECIFICATION= INTP
C...	   NARRAY(6) = FINAL PRINTOUT SPECIFICATION	=IFP
C...	   NARRAY(7) = LOGICAL UNIT NUMBER FOR OUTPUT =IK
C...	   NARRAY(8) = # OF ITERATIONS = KITER
C...     BB= PARAMETER VALUES
C...     IBB= FIXED PARAMETER NUMBERS. FOR EXAMPLE IF THERE ARE 10 PARAMETERS AND 3 ARE FIXED THE INPUT MAY LOOK LIKE "3,7,9" FOR FIXED PARAMETERS 3 AND 7 AND 9



      IMPLICIT NONE
      COMMON/BLK1/B(20),P(20),N,M,K
      COMMON/BLK2/A(40,20),SA(20),K2,IK
      COMMON/BLK3/BS(20),DB(20),G(20),K3
      COMMON/BLK4/AL,DELTA,E,FF,GAMCR,T,TAU,ZETA,PHI,SE,PHICR
      COMMON/BLK5/IB(20),IP
      COMMON/BLK6/WEIGHT(2000)
!      COMMON /BLK7/ ITAYL
      COMMON /BLK8/LJ
      COMMON /R_AND_N_BLK/RSQAVG,NCYLAVG 
      REAL*8 Y(2000),X(2000),RES(2000)
      CHARACTER*60 PARMFILE
      INTEGER IBB(20),NARRAY(8)
      REAL*8 BB(20)
      REAL*8 CONST(8),SCONST(8)
      REAL*8  B,P
      REAL*8  RE(NARRAY(1)),F(NARRAY(1))
      REAL*8  PHI,PHIOLD
      REAL*8  A,SA,BS,DB,G,AL,DELTA,E,FF,GAMCR,T,TAU,ZETA,SE,PHICR
      REAL*8  WEIGHT
      REAL*8  DD,DG,GG,CHISQ
      INTEGER K2,IK,K3,IB,IP,LJ
      REAL*8  GAMMA,XL,CGAM,WS
      INTEGER N,M,K,L,MS
      INTEGER INTP,IFP,KITER,J,I6,I7                        !itayl
      INTEGER NLJ,II,I,FNKW,FKW
      REAL*8  STE,HJTD,OPU,OPL,SPL,SPU
      REAL*8  RSQAVG,NCYLAVG
      DATA    SCONST/1.0D-02,5.0D-3,5.00D-5,4.0D0,45.0D0,2.0D0,1.0D-3,1.00D-31/
      !SCONST VALUES ARE THE CONSTANTS USED IN THE CALCULATIONS
      !AL (THE LAMBDA IN THE lEVENBERG-mARQUARDT FITTIN),DELTA (THE PERCENTAGE OF STEP SIZES IN THE PARTIAL DERIVATIVES)
      !,E(=precision of fit),FF (USED IN SUPPORT PLANE CLACULATIONS),GAMCR (USED TO DETERMINE WHEN TO QUIT),T,TAU(PRECISION USED IN THE EPSILON TEST),ZETA(=machine epsilon)
   
      N=NARRAY(1)
      M=NARRAY(2)
      K=NARRAY(3)
      K2=K
      K3=K
      IP=NARRAY(4)
      INTP=NARRAY(5)
      IFP=NARRAY(6)
      IK=NARRAY(7)
      KITER=NARRAY(8)
      IF(KITER.EQ.(-1))NARRAY(8)=0
      
      F=0.0D0
      
C     SET THE CONSTANTS FOR THE L-M FITTING ROUTINE WHICH ARE READ IN IN THE DATA STATEMENT ABOVE

      CONST=SCONST
      AL=SCONST(1)
      DELTA=SCONST(2)
      E=SCONST(3)
      FF=SCONST(4)
      GAMCR=SCONST(5)
      T=SCONST(6)
      TAU=SCONST(7)
      ZETA=SCONST(8)
      IF(KITER.LE.0)KITER=30
      IF(IK.LE.0)IK=5

 26   IF(IFP.EQ.(-1))GO TO 100
      
      WRITE(IK,2090)
 2090 FORMAT(1H1)

      WRITE(IK,2001)N,K,M,DELTA,E,FF,GAMCR,T,TAU,ZETA,AL
2001  FORMAT(21H0NO OF DATA POINTS IS,I4,22H   NO OF PARAMETERS IS,I3,3X
     1,30HNO OF INDEPENDENT VARIABLES IS,I3/7H0DELTA=,E15.8,5H   E=,E15.
     28,6H   FF=,E15.8,9H   GAMCR=,E15.8/5H   T=,E15.8,7H   TAU=,E15.8,
     38H   ZETA=,E15.8,6H   AL=,E15.8 )


 100  BS=BB		!SET BS AND B TO THE INITIAL GUESS PARAMETERS
      B=BB
      IB=IBB	! SET IB TO LIST OF FIXED PARAMETER NUMBERS
      LJ=0
      
      !*******IF ALL PARAMETERS ARE FIXED DO A CALCULATION OF CHI SQUARE********************************************
      IF (IP.EQ.K) THEN                           
          CALL MODEL(F,Y,X,RE,RES,4)
          CHISQ= SUM(RE*RE*WEIGHT(1:N)) / DBLE(N-K+IP+1)
          WRITE(IK,*)"CHISQ=",CHISQ
          RETURN
      ENDIF    
      !*************************************************************************************************************
      
      CALL SUMSQ(PHI,Y,X,RES)
 130  IF(LJ.GE.KITER)GOTO 404	 !end program if iteration number is GREATER THAN total
C                               number of iterations allowed
      LJ=LJ+1
      WRITE(6,*)'PARMFILE=',PARMFILE
C     BEGIN LJTH ITERATION*************************

      CALL NEWA(Y,X,RES)

 1311 IF(AL.LT..1E-07) GO TO 131
      AL=AL/10.
 131  CALL SCALE
      PHIOLD=PHI

C     STORE MATRIX
13111     DO 132 I=1,K
          II=I+20
              DO 132 J=1,K
              A(II,J)=A(I,J)
132       CONTINUE
      CALL  SOLVE

135       B=BS+DB

C     COMPUTE GAMMA
150   DD=0.
      DG=0.
      GG=0.
          DO 152 J=1,K
          IF(SA(J).EQ.0.)GOTO 152
          GG=GG+G(J)*G(J)/(SA(J)*SA(J))
          DD=DD+DB(J)*DB(J)*SA(J)*SA(J)
  152     DG=DG+DB(J)*G(J)
      XL=SQRT(DD)
      IF(DD*GG.GT.0.)GOTO 160   
      
 155  GAMMA=0.
      GOTO 170
 160  CGAM=DG/SQRT(DD*GG)

	IF(CGAM .GT. 1.)CGAM=.999999
      WS=SQRT(1.-CGAM*CGAM)
      GAMMA=57.2957795*ATAN2(WS,CGAM)

 170  CALL SUMSQ(PHI,Y,X,RES)

      IF(PHI.LE.PHIOLD)GOTO 175  
      IF(GAMMA-GAMCR)300,300,180   
  175 BS=B
      IF (GAMMA.LT.90.) GO TO 190
      
C     GAMMA LAMBDA TEST
 178  IF(AL-1.)190,403,403
  180 AL=AL*10.
      CALL SOLVE
      GOTO 135
      
C     EPSILON TEST
  190 CALL  ETEST(L)
      GOTO(401,200),L
      
C     BEGIN INTERMEDIATE OUTPUT ROUTINE
 200  INTP=1                    !INTP IS FIXED TO 1 COULD BE CHANGED IF YOU WANT MORE INTERMEDIATE OUTPUT
      !*************WRITE TO OUTPUT FILE********************************
      WRITE(IK,2000)
      CHISQ= PHI / (N-K+IP+1)
      WRITE(IK,2002)LJ,CHISQ,AL,(B(J),J=1,K)
      WRITE(IK,*)  'RSQAVG,NCYLAVG=',RSQAVG,NCYLAVG
	write(IK,*)	 'GAMMA INTER=',GAMMA
	write(IK,*)	 'LENGTH OF DB INTER=',XL
	write(IK,*)	 'DB CORRECTION VECTOR INTER=',(DB(J),J=1,K)
    
      !*********WRITE TO SCREEN**************************************
      WRITE(6,2000)
      WRITE(6,2002)LJ,CHISQ,AL,(B(J),J=1,K)
      WRITE(6,*)  'RSQAVG,NCYLAVG=',RSQAVG,NCYLAVG
	write(6,*)	 'GAMMA INTER=',GAMMA
	write(6,*)	 'LENGTH OF DB INTER=',XL
	write(6,*)	 'DB CORRECTION VECTOR INTER=',(DB(J),J=1,K)
      

      IF(INTP.EQ.1)GOTO 130
      CALL NEWA(Y,X,RES)
C     STORE MATRIX
          DO 205 I=1,K
          II=I+20
              DO 205 J=1,K
 205      A(II,J)=A(I,J)
          
      CALL GJR(MS)

      GO TO (207,130),MS
 207  IF(INTP.EQ.2)GOTO 210
      WRITE(IK,2004)
      CALL PRINT1

 210  CALL SCALE
      WRITE(IK,2006)
      CALL PRINT2

C     GET MATRIX FROM STORAGE
          DO 220 I=1,K
          II=I+20
              DO 220 J=1,K
 220      A(I,J)=A(II,J)
          
      IF(LJ.GE.KITER)GO TO 404
      LJ=LJ+1
      GO TO 1311
 300  DB=DB/2.
      B=BS+DB
      
C     GAMMA EPSILON TEST
      CALL ETEST(L)

      GO TO (402,321),L
  321 CALL  SUMSQ(PHI,Y,X,RES)

      IF(PHIOLD.LT.PHI) GO TO 300

 330  BS=B
      GO TO 200
      
C     BEGIN FINAL  PRINTOUT ROUTINE
 401  CONTINUE
      WRITE(IK,2090)
      WRITE(IK,2010)
      GOTO 405
  402 IF(IFP.EQ.(-1))GO TO 4025
      WRITE(IK,2090)
      WRITE(IK,2011)
 4025 IF(IFP.NE.(-1).AND.PHIOLD.LT.PHI)WRITE(IK,2092)
      IF(PHIOLD.GE.PHI)GO TO 4029
      PHI=PHIOLD

 4027 B=BS
4029  CONTINUE    

      GOTO 405
 403  CONTINUE
 
      WRITE(IK,2012)
      GOTO 405
  404 WRITE(IK,2013)
 2013 FORMAT(10H1FORCE OFF)

      NARRAY(8)=-1
       
 405  BS=B
      BB=B
      CHISQ= PHI / (N-K+IP+1)
      WRITE(IK,2002)LJ,CHISQ,AL,(B(J),J=1,K)
	write(ik,*)	 'GAMMA1=',GAMMA
	write(ik,*)	 'LENGTH OF DB1=',XL
	write(ik,*)	 'DB CORRECTION VECTOR=1',(DB(J),J=1,K)
      
      CALL NEWA(Y,X,RES)
      IF(IFP.LE.1)GOTO430
          DO 410 I=1,K
          II=I+20
              DO 410 J=1,K
  410         A(II,J)=A(I,J)
      WRITE(IK,2022)
      CALL PRINT1
      CALL SCALE
      WRITE(IK,2023)
      CALL PRINT2
C     GET MATRIX FROM STORAGE
          DO 420 I=1,K
          II=I+20
              DO 420 J=1,K
 420          A(I,J)=A(II,J)
 430  CALL GJR(MS)

      GO TO (440,435),MS
 435  WRITE(IK,2060)
      GO TO 455
 440  IF(IFP.EQ.0)GOTO450
      WRITE(IK,2024)
      CALL PRINT1
 
 450  CALL SCALE
      WRITE(IK,2025)
      CALL PRINT2

 455  CONTINUE

C--------------------------
      WRITE(IK,20301)
3005  CONTINUE  
      
      CALL MODEL(F,Y,X,RE,RES,4)
 
      DO 460 I=1,N
C      RESIDUAL ARRAY OPTION SATISFIED HERE
460   WRITE(IK,20311)I,X(I),WEIGHT(I),Y(I),F(I),RE(I)
C--------------------------
C     ONE PARAMETER SUPPORT PLANE COMPUTATIONS
      FNKW=N-K+IP
      IF(FNKW.LE.0.)GOTO 589
      FKW=K-IP
      SE=SQRT(PHI/FNKW)
      WRITE(IK,2040)
          DO 470 I=1,K
C     CHECK FOR OMITTED PARAMETERS
          IF(IP.EQ.0)GOTO464
              DO 462 J=1,IP
              IF (I.EQ.IB(J)) GO TO 469
 462          CONTINUE
  464     STE=SA(I)*SE
          HJTD=SQRT(FF*FKW)*STE
          OPL=BS(I)-STE*T
          OPU=BS(I)+STE*T
          SPL=BS(I)-HJTD
          SPU=BS(I)+HJTD
          WRITE(IK,2041)I,STE,OPL,OPU,SPL,SPU
          GOTO 470
 469      WRITE(IK,2042)I
 470      CONTINUE
      IF (IFP.EQ.1) GO TO 602
C     NONLINEAR CONFIDENCE REGION CALCULATIONS
      WS=FKW/FNKW
      PHICR=PHI*(1.0D0+WS*FF)
      WRITE(IK,2049)PHICR
      CALL  CONFRG(Y,X,RES)
      IF(IFP.GE.0)WRITE(IK,2090)
      RETURN
 589  WRITE(IK,2060)
!  590 IF(IFP.EQ.0) GO TO 602
 599  IF(IFP.GE.0)WRITE(IK,2090)
      RETURN
C     RETURNING PARAMETERS WITH NO OUTPUT
! 600  DO601J=1,K
!601  BB(J)=B(J)
C      RESIDUAL ARRAY OPTION WITH NO OUTPUT
 602  J=4

 2000 FORMAT(30H0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX/1H0)
 2002 FORMAT(1X/19H0NO OF ITERATIONS =I3/9H0CHI**2 =,E15.8,4X,8HLAMBDA =,E15.8/11H0PARAMETERS/(1X,7E17.8/))
 2003 FORMAT(8H0GAMMA =,E15.8,4X,14HLENGTH OF DB =,E15.8/21H0DB CORRECTION VECTOR/(1X,7E17.8/))
 2004 FORMAT(1X/12H0PTP INVERSE)
 2006 FORMAT(1X/25H0CORRELATION COEFFICIENTS)
 2010 FORMAT(28H0CONVERGENCE BY EPSILON TEST)
 2011 FORMAT(34H0CONVERGENCE BY GAMMA EPSILON TEST)
 2012 FORMAT(33H0CONVERGENCE BY GAMMA LAMBDA TEST)
 2022 FORMAT(1X/11H0PTP MATRIX)
 2023 FORMAT(1X/29H0PTP CORRELATION COEFFICIENTS)
 2024 FORMAT(1X/12H0PTP INVERSE)
 2025 FORMAT(1X/35H0PARAMETER CORRELATION COEFFICIENTS)
! 2030 FORMAT(1X/1H0,6X,11H   OBSERVED,11X,9HPREDICTED,10X,8HRESIDUAL/1X)
20301 FORMAT(1X/'   I        Q(I)        WEIGHT(I) INTENSITY MEASURED(I) INTENSITY CALC(I)   RESIDUAL(I)'/1X)
 2031 FORMAT(1X,I3,3X,2(E15.8,4X),E15.8)
20311 FORMAT(1X,I3,3X,F12.6,2X,F12.6,2X,3(E15.8,4X),E15.8)
 2040 FORMAT(1X/1H0,12X,4H STD,18X,15HONE - PARAMETER,22X,13HSUPPORT PLANE/2X,4HPARA,7X,5HERROR,13X,5HLOWER,13X,5HUPPER,13X,
     1 5HLOWER,13X,5HUPPER)
 2041 FORMAT(2X,I4,5E18.8)
 2042 FORMAT(2X,I4,5X,23HPARAMETER HELD CONSTANT)
 2049 FORMAT(1X/29H0 NONLINEAR CONFIDENCE LIMITS/15H0PHI CRITICAL =E15.8/6H0 PARA,6X,8H LOWER B,8X,10H LOWER PHI,10X,8H UPPER B,
     1 8X,10H UPPER PHI)
 2060 FORMAT(57H0OUTPUT IS ABRIEVIATED DUE TO MATHEMATICAL CONSIDERATIONS)
 2092 FORMAT(50H0CORRECTION VECTOR FOR LAST ITERATION WAS NOT USED)
      RETURN      
      END



      SUBROUTINE SUMSQ(PHI,Y,X,RES)
C      COMPUTES SUM OF SQUARES
      IMPLICIT NONE
      COMMON/BLK1/B(20),P(20),N,M,K
      COMMON/BLK6/WEIGHT(2000)
      REAL*8      Y(2000),X(2000),RES(2000)
      REAL*8      PHI,B,P,WEIGHT
      REAL*8      RE(N),F(N)
      INTEGER     I,N,M,K
      
      PHI=0.
      F=0.
      CALL MODEL(F,Y,X,RE,RES,1)
      PHI=SUM(RE*RE)
      RETURN
      END







      SUBROUTINE CONFRG(Y,X,RES)
CCONFRG               CONFRG - NON LINEAR CONFIDENCE REGION CALCULATIONS
      IMPLICIT NONE
      COMMON/BLK1/B(20),P(20),N,M,K
      COMMON/BLK2/A(40,20),SA(20),K2,IK
      COMMON/BLK3/BS(20),DB(20),G(20),K3
      COMMON/BLK4/AL,DELTA,E,FF,GAMCR,T,TAU,ZETA,PHI,SE,PHICR
      COMMON/BLK5/IB(20),IP
      REAL*8  Y(2000),X(2000),RES(2000)
      REAL*8  B,P
      REAL*8  BS,DB,G
      INTEGER I,IP,IB,IK,J,K,N,M,K2,K3
      REAL*8  DDS,D,DJ,SE,SA,PH,PHICR,PPH,Q,XK1,PHI,XK2,XK3,BC,BL
      REAL*8  PL,BU,PU,A,AL,DELTA,E,FF,GAMCR,T,TAU,ZETA
      LOGICAL NOLO
           
          DO 580 J=1,K
          NOLO=.FALSE.
C         CHECK FOR OMITTED PARAMETERS
          IF(IP.EQ.0) GO TO 509
              DO 504 I=1,IP
              IF(J.EQ.IB(I))GOTO 506
 504          CONTINUE
          GOTO 509
 506      WRITE(IK,2042)J
          GO TO 580
 509      DDS=-1.
 510      D=DDS
          DJ=SE*SA(J)
          B(J)=BS(J)+D*DJ
          CALL SUMSQ(PH,Y,X,RES)
          IF(PH.LT.PHICR)GOTO530
 520      D=D/2.
          IF(ABS(D).LE..001)GOTO570
          B(J)=BS(J)+D*DJ
          CALL SUMSQ(PPH,Y,X,RES)
          IF(PPH-PHICR)540,540,520
 530      D=D+DDS
          IF(ABS(D).GE.5.0) GO TO 570
          B(J)=BS(J)+D*DJ
          CALL SUMSQ(PPH,Y,X,RES)
          IF(PPH.LT.PHICR)GOTO530
 540      Q=1.-D
          XK1=PHI/D+PH/Q-PPH/(D*Q )
          XK2=-PHI*(1.+D)/D-PH*D/Q+PPH/(D*Q)
          XK3=PHI-PHICR
          BC=(-XK2+SQRT(XK2**2 -4.*XK1*XK3))/(2.*XK1)
          IF(DDS.GT.0.) GO TO 550
          B(J)=BS(J)-BC*DJ
          BL=B(J)
          CALL SUMSQ(PL,Y,X,RES)
 548      DDS=1.
          GOTO 510
  550     B(J)=BS(J)+BC*DJ
          BU=B(J)
          CALL SUMSQ(PU,Y,X,RES)
          GOTO 576
  570     IF(DDS.GT.0.) GO TO 571
          NOLO=.TRUE.
          GOTO 548
 571      IF(NOLO)GOTO 575
C         OMITTING UPPER LIMITS
          WRITE(IK,2055)J,BL,PL
          GOTO 580
C         OMITTING BOTH
 575      WRITE(IK,2056)J
          GOTO 580
 576      IF(NOLO)GOTO578
          WRITE(IK,2052)J,BL,PL,BU,PU
          GOTO 580
C         OMITTING LOWER LIMITS
 578      WRITE(IK,2053)J,BU,PU
 580      B(J)=BS(J)
 2042 FORMAT(2X,I3,5X,23HPARAMETER HELD CONSTANT)
 2052 FORMAT(2X,I3,4E18.8)
 2055 FORMAT(2X,I3,2E18.8,11H  NOT FOUND)
 2053 FORMAT(2X,I3,11H  NOT FOUND,25X,2E18.8)
 2056 FORMAT(2X,I3,18X,11H  NOT FOUND)
      RETURN
      END







      SUBROUTINE NEWA(Y,X,RES)
CNEWA         NEWA - CALCULATES THE ALPHA MATRIX FROM THE L-M ALGORITM, A, AND GRADIENT VECTOR, G.
      IMPLICIT NONE
      COMMON/BLK1/B(20),P(20),N,M,K   
      COMMON/BLK2/A(40,20),SA(20),K2,IK
      COMMON/BLK3/BS(20),DB(20),G(20),K3
      COMMON/BLK4/AL,DELTA,E,FF,GAMCR,T,TAU,ZETA,PHI,SE,PHICR
      COMMON/BLK5/IB(20),IP
!      COMMON /BLK7/ ITAYL
      REAL*8      Y(2000),X(2000),RES(2000)
      REAL*8      B,P
      REAL*8      F(N),FDEL(N)
      REAL*8      RE(N),RD(N)
      REAL*8      A,SA,BS,DB,G,AL,DELTA,E,FF,GAMCR,T,TAU,ZETA,PHI,SE,PHICR
      REAL*8      PHIP,AB,DELTAB
      REAL*8      P1(K,N),G1(K,N)
      INTEGER     KS,IK,K3,IB,IP                              !,ITAYL
      INTEGER     J,K,I,II,N,JJ,M,K2,NN

C...		N IS THE NUMBER OF DATA POINTS
C...		M IS THE NUMBER OF INDEPENDENT VARIABLES
C...      IP IS THE NUMBER OF FIXED PARAMETERS
C...		K IS THE TOTAL NUMBER OF PARAMETERS
C...		RE IS THE RESIDUAL DIVIDED BY SIGMA
C...		B ARE THE INITIAL VALUES FOR THE PARAMETERS
C...		P ARE PARTIAL DERIVATIVES
C...		Y AND X ARE THE DATA
C...		RES IS THE RESULT
C...		A IS MIXED PARTIAL DERIVATIVE MATRIX 
C...		IB ARE THE ARRAY ELEMENT NUMBERS OF THE FIXED PARAMETERS     
      
    !  WRITE(6,99997)
99997	FORMAT(/' *** SUBROUTINE NEWA ***')
          DO 1 J=1,K
          G(J)=0.0D0
          P(J)=0.0D0
              DO 3 NN=1,N
              P1(J,NN)=0.0D0
              G1(J,NN)=0.0D0
3             CONTINUE              
              DO 1 I=1,K
 1            A(J,I)=0.0D0
          PHIP=0.0D0

C             LOOK FOR PARTIALS
              J=2 
              F=0.
              RE=0.
              CALL MODEL(F,Y,X,RE,RES,J)          !J=3 AFTER THIS CALL

              RD=RE
              PHIP=SUM(RE*RE)
     	            DO 30 JJ=1,K
C                 CHECK FOR OMITTED PARAMETERS
                  IF(IP.GT.0)GOTO 25      !ip is always gt 0 unless you LET all parameters VARY- i DON'T UNDERSTAND THIS ONE
   10             GOTO(20,31,20),J
C                 COMPUTE PARTIALS IF NECESSARY
 20               AB=B(JJ)
                      IF (B(JJ).NE.0.) GOTO 118	! CORRECTION FOR ZERO
                      DELTAB=DELTA	                ! PARAMETER
                      GOTO 117
                    
118                DELTAB=DELTA*AB            
117                B(JJ)=AB+DELTAB		!CALCULATE NEW VALUE FOR PARAMETER BY
									!ADDING DELTA * LAST PARAMETER VALUE 
									!UNLESS THE PARAMETER IS ZERO INWHICH CASE THE 
									!DELTA IS SIMPLY DELTA

                  J=1
                  FDEL=0.
                  RE=0.
                  CALL MODEL(FDEL,Y,X,RE,RES,J)               !J RETURNS AS 1
                  RE=RD              
                  P1(JJ,1:N)=(FDEL-F)/(DELTAB)			!CRUDE FORM OF THE DERIVATIVE DF/DA  
                  B(JJ)=AB
                  GOTO 31
   25                 DO 26 I=1,IP
                      IF(JJ.EQ.IB(I)) GO TO 29	!CHECK TO SEE IF PARAMETER JJ IS FIXED											
 26                   CONTINUE
                  GOTO 10
  29              P1(JJ,1:N)=0.
C                 USING PARTIALS AT ITH DATA POINT
  31              G1(JJ,1:N)=RE*P1(JJ,1:N)	!GRADIENT IN CHI: I THINK ITS DCHI/DF
										!TIMES DF/DA(I)
!                 ********************NOTE P IS ONLY USED IN THIS ROUTINE SO i'LL LET THE P'S STAY AT 0.0D0
                  G(JJ)=SUM(G1(JJ,1:N))
   30             CONTINUE
              

              DO 40 I=1,K
                  DO 40 J=I,K
   40             A(I,J)=SUM(P1(I,1:N)*P1(J,1:N))	    !MIXED PARTIAL DERIVATIVE MATRIX
									 	            !THE I,JTH ELEMENT IS PARTAL OF F
										            !W/R THE ITH PARAMETER TIMES THE 
										            !PARTIAL OF F W/R ITH PARAMETER
										            ! SUMMED OVER THE EVALUATION AT ALL 
										            ! VALUES OF X 

          DO 55 I=1,K
              DO 55 J=I,K
   55         A(J,I)=A(I,J)
C         A(I,I)=1.0 FOR OMITTED PARAMETER I
      IF(IP.EQ.0)RETURN
          DO 60 I=1,IP
              DO 60 J=1,K
 60       IF(J.EQ.IB(I))A(J,J)=1.
         RETURN
      END







      SUBROUTINE SCALE
CSCALE      SUBROUTINE SCALE
C     SCALES ACCORDING TO DIAGONAL ELEMENTS
      IMPLICIT NONE
      COMMON/BLK2/A(40,20),SA(20),K2,IK
      REAL*8 A,SA,WS
      INTEGER K2,IK,K,I,J

      K=K2
          DO 20 I=1,K
          IF(A(I,I).GT.0.) GO TO 15
          SA(I)=0.
          GOTO 20
 15       SA(I)=SQRT(A(I,I))
 20       CONTINUE


          DO 50 I=1,K
              DO 40 J=1,I
              WS=SA(I)*SA(J)
              IF(WS.GT.0.)GOTO 30
              A(I,J)=0.
              GOTO 40
 30           A(I,J)=A(I,J)/WS
 40           A(J,I)=A(I,J)
 50       A(I,I)=1.0
      RETURN
      END



      SUBROUTINE SOLVE
CSOLVE      SOLVES (PTP)(DB)=(G)  WHERE PTP IS STORED IN A(I+20,J)
C     SOLVES  A SET OF LINEAR EQUATIONS IN DB DETERMINED BY MATRIX
C     A AND VECTOR G. USES SUBROUTINE GJR TO INVERT MATRIX
      IMPLICIT NONE
      COMMON/BLK2/A(40,20),SA(20),K2,IK
      COMMON/BLK3/BS(20), DB(20), G(20),K3
      COMMON/BLK4/AL,DELTA,E,FF,GAMCR,T,TAU,ZETA,PHI,SE,PHICR
      REAL*8 A,SA,BS,DB,G,AL,DELTA,E,FF,GAMCR,T,TAU,ZETA,PHI,SE,PHICR
      INTEGER K2,IK,K3
      INTEGER K,L,I,II,J,MS
      
      K=K2
      L=1
C     GET MATRIX FROM STORAGE
 1        DO 10 I=1,K
          II=I+20
              DO 9 J=1,K
 9            A(I,J)=A(II,J)
 10       A(I,I)=1.+AL
 20       CALL GJR(MS)

 25       DO 40 I=1,K
          DB(I)=0.
          IF(SA(I).LE.0.)GOTO40
              DO 30 J=1,K
              IF(SA(J).LE.0.)GOTO30
              DB(I)=A(I,J)*G(J)/SA(J) +DB(I)
 30           CONTINUE
          DB(I)=DB(I)/SA(I)
 40       CONTINUE
      RETURN
 100  AL=AL*10.
      L=L+1
      IF(L.GE.6)STOP
      GOTO1
      END



      SUBROUTINE ETEST(ML)
CETEST      SUBROUTINE ETEST
      IMPLICIT NONE
      COMMON/BLK1/B(20),P(20),N,M,K
      COMMON/BLK3/BS(20),DB(20),G(20),K3
      COMMON/BLK4/AL,DELTA,E,FF,GAMCR,T,TAU,ZETA,PHI,SE,PHICR
      REAL*8      B,P,BS,DB,G,AL,DELTA,E,FF,GAMCR,T,TAU,ZETA,PHI,SE,PHICR,EPS,W
      INTEGER     N,M,K,K3,I,ML
      
      EPS=E
      ML=1
          DO 20 I=1,K
          W=ABS(DB(I))/(TAU+ABS(B(I)))
          IF (W.GE.EPS) GO TO 30
   20     CONTINUE
      GO TO 40
   30 ML=2
   40 RETURN
      END


      SUBROUTINE GJR(MSING)
C     GJR - INVERTS A MATRIX IN A(I,J), I=1,20,J=1,20
C     GAUSS-JORDAN-RUTISHAUSER MATRIX INVERSION WITH DOUBLE PIVOTING
      IMPLICIT NONE
      COMMON/BLK2/A(40,20),SA(20),K2,IK
      COMMON/BLK4/AL,DELTA,E,FF,GAMCR,T,TAU,ZETA,PHI,SE,PHICR
      REAL*8      B(20),C(20)
      INTEGER     P(20),Q(20)
      REAL*8      A,SA,AL,DELTA,E,FF,GAMCR,T,TAU,ZETA,PHI,SE,PHICR
      INTEGER     K2,IK
      REAL*8      EPS,PIVOT,Z
      INTEGER     N,K,MSING,I,J,L,M
      
      EPS=ZETA
      N=K2
      MSING=1
          DO 10 K=1,N
C         DETERMINATION OF PIVOT ELEMENT
          PIVOT=0.
              DO 20 I=K,N
                  DO 20 J=K,N
                  IF(ABS(A(I,J))-ABS(PIVOT))20,20,30
 30               PIVOT=A(I,J)
                  P(K)=I
                  Q(K)=J
 20               CONTINUE
              IF(ABS(PIVOT)-EPS)40,40,50
C     EXCHANGE OF PIVOTAL ROW WITH KTH ROW
   50         IF (P(K).EQ.K) GO TO 80
                  DO 70 J=1,N
                  L=P(K)
                  Z=A(L,J)
                  A(L,J)=A(K,J)
 70               A(K,J)=Z
C     EXCHANGE OF COLUMN
 80           IF(Q(K).EQ.K)GOTO 90
                  DO 100 I=1,N
                  L=Q(K)
                  Z=A(I,L)
                  A(I,L)=A(I,K)
  100             A(I,K)=Z
 90               CONTINUE
C             JORDAN STEP
                  DO 110 J=1,N
                  IF(J.EQ.K)GOTO120
                  B(J)=-A(K,J)/PIVOT
                  C(J)=A(J,K)
                  GOTO 140
 120              B(J)=1./PIVOT
                  C(J)=1.
 140              A(K,J)=0.
 110              A(J,K)=0.
                  
                  DO 10 I=1,N
                      DO 10 J=1,N
   10                 A(I,J)=A(I,J)+C(I)*B(J)
C     REORDERING THE MATRIX
          DO 155 M=1,N
          K=N-M+1
          IF (P(K).EQ.K) GOTO 170
              DO 180 I=1,N
              L=P(K)
              Z=A(I,L)
              A(I,L)=A(I,K)
 180          A(I,K)=Z
 170      IF(Q(K).EQ.K)GOTO 155
              DO 150 J=1,N
              L=Q(K)
              Z=A(L,J)
              A(L,J)=A(K,J)
 150          A(K,J)=Z
 155      CONTINUE
      RETURN
 40   WRITE(IK,45)P(K),Q(K),PIVOT
 45   FORMAT(20H0SINGULAR MATRIX  I=,I3,4H  J=,I3,8H  PIVOT=E16.8/)
      MSING=2
      RETURN
      END


      SUBROUTINE PRINT1
CPRINT1     SUBROUTINE PRINT1
C     PRINTS A: K BY K MATRIX
      IMPLICIT NONE
      COMMON/BLK2/A(40,20),SA(20),K2,IK
      REAL*8      A,SA
      INTEGER     K2,IK,K,L,JJ,LL,I,J
      
      K=K2
      L=1
 5    JJ=7*L
      LL=JJ-6
      IF(K.LT.LL)GOTO 30
      IF(K.LT.JJ)GOTO 20
      WRITE(IK,105)LL,JJ
          DO 15 I=1,K
 15       WRITE(IK,106)(A(I,J),J=LL,JJ)
      L=L+1
      GOTO 5
 20   WRITE(IK,105)LL,K
          DO 25 I=1,K
 25       WRITE(IK,106)(A(I,J),J=LL,K)
 30   RETURN
 105  FORMAT(1X/8H0COLUMNS,I4,9H  THROUGH,I4)
 106  FORMAT(1X,7E17.8)
      END



CPRINT2     SUBROUTINE PRINT2
      SUBROUTINE PRINT2
C     PRINTS A K BY K CORRELATION COEFFICIENT MATRIX
      COMMON/BLK2/A(40,20),SA(20),K2,IK
      REAL*8      A,SA
      INTEGER     K2,IK
      
      L=1
      K=K2
    5 JJ=13*L
      LL=JJ-12
      IF(K.LT.LL)GOTO30
      IF(K.LT.JJ)GOTO20
      WRITE(IK,105)LL,JJ
          DO 15 I=1,K
 15       WRITE(IK,107)(A(I,J),J=LL,JJ)
      L=L+1
      GOTO 5
 20   WRITE(IK,105)LL,K
          DO 25 I=1,K
 25       WRITE(IK,107)(A(I,J),J=LL,K)
 30   RETURN
 105  FORMAT(1X/8H0COLUMNS,I4,9H  THROUGH,I4)
  107 FORMAT(1X,13F9.4)
      END
