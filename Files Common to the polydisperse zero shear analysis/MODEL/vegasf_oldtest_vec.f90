
!***************************************THIS USES THE ORIGINAL lAPAGE VERSION OF VEGAS*******************************

      SUBROUTINE VEGASF_OLD(FUNCNAME,AVGI,ACC,NCALL,ITMX,NDMX,NCYL,NDIM,XU,XL,COSPSI,SINPSI,Q_VEC,NQ)
!
!   SUBROUTINE PERFORMS NDIM-DIMENSIONAL MONTE CARLO INTEG'N
!      - BY G.P. LEPAGE    SEPT 1976/(REV)AUG 1979
!      - ALGORITHM DESCRIBED IN J COMP PHYS 27,192(1978)
!       THIS IS A VECTORIZED VERSION OF THE ORIGINAL TO SIMULTANEOUSLY DO MULTIPLE Q VALUES; SPEEDS UP THE CALCULATION A LOT      
!
        IMPLICIT NONE
        EXTERNAL      FUNCNAME
        INTEGER       NCYL       
        INTEGER       ITMX,NDEV,NCALL,NDIM,NPRN,NDMX,MDS,IT,NQ,INIT
        INTEGER       NDO,J,ND,NG,NPG,K,NDM,I
        INTEGER       KG(ndim),IA(ndim)   
        REAL*8        D(ndmx,ndim),DI(ndmx,ndim),XIN(ndmx),R(ndmx),DX(ndim),DT(ndim),X(ndim)
        REAL*8        RAND(ndim)
        REAL          RAND1(ndim)
        REAL*8        ACC,CALLS,DXG,DV2G,XND,XJAC
        REAL*8        Q(nq),FB(nq),F2B(nq),FVEC(nq),F(nq),F2(nq),TI(nq),TSI(nq),TI2(nq),WGT1(nq),SI(nq),SWGT(nq),SCHI(nq)
        REAL*8        AVGI(nq),CHI2A(nq),SD(nq),FXN(nq)
        REAL*8        RC,XN,DR,XO,WGT
        REAL*8        XI(ndmx,ndim)
        REAL*8        XI2(ndmx,ndim)            ! DUPLICATE FOR DEBUGGING
        REAL*8        ALPH,ONE,PI
        REAL*8        Q_VEC(nq)
        REAL*8        XU(ndim)
        REAL*8        XL(ndim)
        REAL*8        COSPSI,SINPSI
        REAL*8        SD_RATIO(NQ)
        INTEGER       L,LL,MAXLO(NQ)
        REAL*8        GRIDVEC(NDMX)
        REAL*8        XN1(NDIM),XO1(NDIM),RC1(NDIM)
       
        INIT=0
1       NPRN=0                               !NPRN CONTROLS INTERMEDIATE OUTPUT FROM VEGAS, 1 MEANS FULL OUTPUT FROM EACH ITERATION
        ONE=1.0D0     
        ALPH=1.0D0                            !THIS IS THE RATE AT WHICH THE GRID IS MODIFIED-- i RESET FROM ALPHA = 1.5 TO =1.2 AND SEEMS MORE STABLE
        SD=0.0d0
        MDS=0                                !THIS ALLOWS PROGRAM TO CHOOSE BETWEEN IMPORTANCE AND STRATIFIED SAMPLING- LOOKS LIKE KEEPING IT AT ZERO KEEPS IT FROM BLOWING UP
        
        
!       INIT=0 IS THE NORMAL ENTRY POINT         
            IF (INIT.LE.0) THEN
            NDO=1
            XI=1.0D0                            !arbitrary choice of initialization to 1
            ENDIF
            
!       INIT =1 INITIALIZES CUMMULATIVE VARIABLES, BUT NOT GRID
            IF (INIT.LE.1) THEN
            IT=0   
            SI=0.d0
            SWGT=0.0d0
            SCHI=0.0d0
            ENDIF
                
 !       INIT=2 - INHERIT PREVIOUS GRID AND ITS ANSWERS
            IF (INIT.LE.2) THEN      
            ND=NDMX                               !NDMAX IS MAXIMUM NUMBER OF GRID POINTS ALONG EACH AXIS
            NG=1
            IF (MDS.EQ.0) GO TO 2                  !MDS =0 IMPORTANCE SAMPLING ONLY, MDS NE 0 USES IMPORTANCE SAMPLING
            NG=(NCALL/2.d0+0.25)**(1.d0/NDIM)          !NCALLS IS THE APPROXIMATE NUMBER OF EVALUATIONS OF THE INTEGRAND PER ITERATION; NG the average INTEGER NUMBER OF CALCULATED INTEGRANDS
                                                    ! IN EACH DIMENSION A TWO POINT INTEGRATION IS USED SO NCALLS/2 IS USED FOR NG
            MDS=-1
            IF ((2*NG-NDMX).LT.0) GO TO 2          !NG IS THE CRITERION FOR SWITCHING FROM IMPORTANCE TO STRATIFIED SAMPLING; 
                                                    !IF WE SKIP TO LINE 2 WE USE IMPORTANCE SAMPLING
            MDS=1
            NPG=NG/NDMX+1                         !AVERAGE NUMBER OF CALLS PER BOX ALONG A SINGLE DIMENSION PLUS 1 ( I.E. NUMBER PER GRID)
            ND=NG/NPG                             !MAXIMUM NUMBER OF SLICES ALONG AN AXIS IN ANY ONE DIMENSION
            NG=NPG*ND                             !NEW NUMBER OF GRID POINTS ALONG ANY AXIS
            ENDIF

 2      K=NG**NDIM                            !INTEGER NUMBER OF INTEGRAND EVALUATIONS ALONG EACH AXIS
        NPG=NCALL/K                           !
        IF (NPG.LT.2) NPG=2                    !         
        CALLS=DBLE(NPG*K)                           !NEW NUMBER OF ACTUAL CALLS
        DXG=DBLE(1.0d0/NG)                            ! 1/(INTEGER NUMBER OF BOXES PER DIRECTION)
        DV2G=(CALLS*DXG**DBLE(NDIM))**2/DBLE(NPG)/DBLE(NPG)/(DBLE(NPG)-ONE)           !TOTAL NUMBER OF CALLS PER GRID VOLUME/
        XND=DBLE(ND)                                !NUMBER OF GRIDS IN ANY DIRECTION
        NDM=ND-1                              ! NUMBER OF GRID INCREMENTS IN ANY DIRECTION MINUS 1-- FOR DO LOOPS
        DXG=DXG*XND                           ! % OF NDMAX BOXES USED- WILL BE USED TO SCALE HOW MANY ELEMENTS OF THE XVECTOR
        XJAC=ONE/CALLS
            DO 3 J=1,NDIM                       !XJAC IS THE VOLUME OF THE TOTAL HYPERCUBE IN "REAL" SPACE DIVIDED BY THE NUMBER CALLS
            DX(J)=XU(J)-XL(J)
3           XJAC=XJAC*DX(J)
            
            IF (ND.EQ.NDO) GO TO 8

        GRIDVEC=dble([1:NDMX])/NDMX
        XI=reshape(GRIDVEC,[NDMX,NDIM],GRIDVEC)       !CONSTRUCT AN NDMX x NDIM MATRIX WITH VALUES FROM 0 TO ONE IN STEPS OF 1/NDMX
        XI(NDMX,1:NDIM)=1.0D0
        NDO=ND

!8       IF (NPRN.GE.0) THEN
!          WRITE(20,200) NDIM,CALLS,IT,ITMX,ACC,NPRN,ALPH,MDS,ND,(XL(J),XU(J),J=1,NDIM)
!        ENDIF 

 !     ENTRY VEGAS3(NDIM,FXN,AVGI,SD,CHI2A)
!         - MAIN INTEGRATION LOOP

8        CALL RANDOM_SEED  
        DT=0.0D0
        
9       IT=IT+1  

        TI=0.0d0                                   !INITIALIZATION
        TSI=0.0d0                                    !INITIALIZATION
        KG=1
        D=0.0D0
!        DI=0.0D0
!
 11     FB= 0.0d0                                   !
        F2B=0.0d0                                   !
        
            DO 12 K=1,NPG                             
            CALL RANDOM_NUMBER(RAND)
            WGT=XJAC                                  !WEIGHT FUNCTION 

            XN1=0.0D0
            RC1=0.0D0
            XO1=0.0D0
            XN1=(DBLE(KG)-RAND)*DXG+ONE
            IA=INT(XN1)            
                DO 15 J=1,NDIM
                IF (IA(J).GT.1)THEN
                    XO1(J)=XI(IA(J),J) - XI(IA(J)-1,J)
                    RC1(J)=XI(IA(J)-1,J) + (XN1(J)-IA(J))*XO1(J)
                ELSE
                    XO1(J)=XI(IA(J),J)
                    RC1(J)=(XN1(J)-IA(J)) *XO1(J)
                ENDIF
15              CONTINUE
            X=XL+RC1*DX
            WGT=WGT*PRODUCT(XO1)*(XND**NDIM)            
            F=WGT

            CALL FUNCNAME(X,FXN,NCYL,COSPSI,SINPSI,Q_VEC,NQ,K,IT)     

            F=F*FXN                                     !THE FUNCTION TIMES 1/P IN LEPAGE PAPER
            F2=F*F
            FB=FB+F
            F2B=F2B+F2
            !**********************************
                DO 16 J=1,NDIM
                !**************ONLY USE THE middle value of the q's to redefine the grid for all q          
!                DI(IA(J),J)=DI(IA(J),J)+F(NQ/2)                           !THIS IS THE INTEGRAL OF THE FUNCTION FOR THE 
 16         IF (MDS.GE.0) D(IA(J),J)=D(IA(J),J)+F2(NQ/2)               !D IS SET AS A VECTOR OF SIZE ND WITH EACH ELEMENT BEING F*WEIGHT SQUARE
12          CONTINUE
        
        F2B=DSQRT(F2B*NPG)
        F2B=(F2B-FB)*(F2B+FB)
        TI=TI+FB
        TSI=TSI+F2B
        IF (MDS.GE.0) GO TO 18
            DO 17 J=1,NDIM
 17         D(IA(J),J)=D(IA(J),J)+F2B(NQ/2)
 18     K=NDIM
 19     KG(K)=INT(MOD(KG(K),NG))+1
        IF (KG(K).NE.1) GO TO 11
        K=K-1
        IF (K.GT.0) GO TO 19
        
!       COMPUTE FINAL RESULTS FOR THIS ITERATION
        TSI=TSI*DV2G
        TI2=TI*TI
        WGT1=ONE/TSI
        SI=SI+TI*WGT1
        SWGT=SWGT+WGT1
        SCHI=SCHI+TI2*WGT1
        AVGI=SI/SWGT
        CHI2A=(SCHI-SI*AVGI)/(DBLE(IT)-.9999d0)
        SD=DSQRT(ONE/SWGT)
        TSI=DSQRT(TSI)

!   REFINE GRID ***************************************************     
21      DO 23 J=1,NDIM
        XO=D(1,J)
        XN=D(2,J)
        D(1,J)=(XO+XN)/2.0d0
        DT(J)=D(1,J)
            DO 22 I=2,NDM
            D(I,J)=XO+XN
            XO=XN
            XN=D(I+1,J)
            D(I,J)=(D(I,J)+XN)/3.d0
 22         DT(J)=DT(J)+D(I,J)
        D(ND,J)=(XN+XO)/2.d0
 23     DT(J)=DT(J)+D(ND,J)
!
            DO 28 J=1,NDIM
            RC=0.0d0
                DO 24 I=1,ND
                R(I)=0.0d0
                IF (D(I,J).LE.0.d0) GO TO 24
                XO=DT(J)/D(I,J)
                R(I)=((XO-ONE)/XO/DLOG(XO))**ALPH
 24             RC=RC+R(I)
            RC=RC/XND
            K=0
            XN=0.d0
            DR=XN
            I=K
 25         K=K+1
            DR=DR+R(K)
            XO=XN
            XN=XI(K,J)
 26         IF (RC.GT.DR) GO TO 25
            I=I+1
            DR=DR-RC
            XIN(I)=XN-(XN-XO)*DR/R(K)
            IF (I.LT.NDM) GO TO 26
                DO 27 I=1,NDM
 27             XI(I,J)=XIN(I)
            XI(ND,J)=ONE
28          CONTINUE
      
         !**************************************************************
!        WRITE(6,*)'IT,ITMX,NCYL,NCALL=',IT,ITMX,NCYL,NCALL

!        PAUSE
        SD_RATIO=SD/AVGI
        MAXLO=FINDLOC(SD_RATIO,DBLE(MAXVAL(SD_RATIO)),KIND=4)                     !CALCULATE THE POSITION OF THE MAXIUM ERROR IN AVGI IN ALL FOR THE Q'S AND QUIT ONCE THAT MAX IS LOWER THAN THE STANDARD DEVIATION

        IF (IT.LT.ITMX.AND.ACC*AVGI(MAXLO(1)).LT.SD(MAXLO(1))) GO TO 9
            IF(IT.EQ.ITMX) THEN
!            WRITE(6,*)'VEGAS FAILED TO CONVERGE- INCREASING NCALL AND STARTING OVER'
!            WRITE (6,*)' OLD NCALL=',NCALL
            NCALL=NCALL+50000
!            WRITE(6,*)'NEW NCALL=',NCALL
            GOTO 1
            ENDIF
               

200     FORMAT(/35H INPUT PARAMETERS FOR VEGAS:  NDIM=,I3,8H  NCALL=,F8.0/28X,5H  IT=,I5,7H  ITMX=,I5/28X, &
            6H  ACC=,G10.3/28X,7H  NPRN=,I3,7H  ALPH=,F5.2/28X,6H  MDS=,I3,6H   ND=,I4/28X,10H  (XL,XU)=,(T40,2H( ,G13.6,3H , ,G13.6,2H )))

201     FORMAT(///21H INTEGRATION BY VEGAS//14H ITERATION NO.,I3, 14H:   INTEGRAL =,3G15.8/21X,10HSTD DEV  =,&
            3G15.8/ 34H ACCUMULATED RESULTS:   INTEGRAL =,3G15.8/ 24X,10HSTD DEV  =,3G15.8/24X,17HCHI**2 PER IT'N =,3G15.8)
 
202          FORMAT(/15H DATA FOR AXIS ,I2/25H    X       DELTA I       ,24H   X       DELTA I       ,18H   X       DELTA I /   &
            (1H ,F7.6,1X,G11.4,5X,F7.6,1X,G11.4,5X,F7.6,1X,G11.4))
             
             

      RETURN
    END

