!THIS IS THE MAIN PROGRAM TO CALCULATE THE SCATTERED INTENSITY AS A FUNC OF Q 
!FOR N CONNECTED CYCLINDERS UNDER SHEAR- WRITTEN BY GREG SMITH JANUARY 2021
	COMMON /BLK1/B(20),P(20),RE,N,M,K

					!
	PARAMETER (MAXDAT=50000) !MAXIMUM NUMBER OF DATA POINTS
	PARAMETER (MAXPAR=20)	! MAXIMUM NUMBER OF PARAMETERS
					!
					!
      INTEGER NARRAY,STEP,II,NS,NCYL,LCNT,I,SHEARCNT         	
      DOUBLE PRECISION YOUT(MAXDAT) ,QX(MAXDAT),QY(MAXDAT),Q(MAXDAT),COSPSI(MAXDAT),SINPSI(MAXDAT)
      DOUBLE PRECISION QXNEW(MAXDAT),QYNEW(MAXDAT)
      DOUBLE PRECISION QI,CS,SS

	DOUBLE PRECISION QZ,DQZ,DFLAG,QV
      DOUBLE PRECISION QX1,QY1
      DOUBLE PRECISION CSYQYL
      DOUBLE PRECISION EXPTPHI(16),St_List(16)
      CHARACTER(110) FILENAME(16) !changed from CHARACTER(80) to CHARACTER(40) to accomodate for long file name
      CHARACTER(10) TIME1,TIME2,TIME3
      INTEGER(2)  TVALS1(8),TVALS2(8)
      
      
      CALL DATE_AND_TIME(TIME1,TIME2,TIME3,TVALS1)
      
      FILENAME(1)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.27_4cyl_L431_Pq_3det.dat'
      FILENAME(2)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.275_4cyl_L431_Pq_3det.dat'
      FILENAME(3)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.28_4cyl_L431_Pq_3det.dat'
      FILENAME(4)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.282_4cyl_L431_Pq_3det.dat'
      FILENAME(5)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.284_4cyl_L431_Pq_3det.dat'
      FILENAME(6)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.286_4cyl_L431_Pq_3det.dat'
      FILENAME(7)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.288_4cyl_L431_Pq_3det.dat'
      FILENAME(8)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.29_4cyl_L431_Pq_3det.dat'
      FILENAME(9)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.295_4cyl_L431_Pq_3det.dat'
      FILENAME(10)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.298_4cyl_L431_Pq_3det.dat'
      FILENAME(11)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.3_4cyl_L431_Pq_3det.dat'
      FILENAME(12)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.302_4cyl_L431_Pq_3det.dat'
      FILENAME(13)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.304_4cyl_L431_Pq_3det.dat'
      FILENAME(14)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.306_4cyl_L431_Pq_3det.dat'
      FILENAME(15)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.308_4cyl_L431_Pq_3det.dat'
      FILENAME(16)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.31_4cyl_L431_Pq_3det.dat'
      !FILENAME(1)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.736_4cyl_L431_Pq_3det.dat'
      !FILENAME(2)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.7365_4cyl_L431_Pq_3det.dat'
      !FILENAME(3)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.737_4cyl_L431_Pq_3det.dat'
      !FILENAME(4)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.7375_4cyl_L431_Pq_3det.dat'
      !FILENAME(5)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.738_4cyl_L431_Pq_3det.dat'
      !FILENAME(6)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.7385_4cyl_L431_Pq_3det.dat'
      !FILENAME(7)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.739_4cyl_L431_Pq_3det.dat'
      !FILENAME(8)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.7395_4cyl_L431_Pq_3det.dat'
      !FILENAME(9)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.74_4cyl_L431_Pq_3det.dat'
      !FILENAME(10)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.7405_4cyl_L431_Pq_3det.dat'
      !FILENAME(11)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.741_4cyl_L431_Pq_3det.dat'
      !FILENAME(12)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.7415_4cyl_L431_Pq_3det.dat'
      !FILENAME(13)='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\PHIO_0_St0.742_4cyl_L431_Pq_3det.dat'

      
      !EXPTPHI=(/0.785,0.676,0.4605,0.199,0.133,0.118,0.04326,0.015445,0.001218/)
      !EXPTPHI=(/0.785398,0.75,0.7,0.65,0.6,0.5,0.4,0.3,0.2,0.1/)
      !EXPTPHI=(/0.1989,0.1180,0.04326,0.001873,0./)
      !St_List=(/0.29,0.39,0.62,0.74,0.765/)
      !EXPTPHI=(/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)!13 files
      !St_List=(/0.736,0.7365,0.737,0.7375,0.738,0.7385,0.739,0.7395,0.74,0.7405,0.741,0.7415,0.742/)!13 files
      EXPTPHI=(/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./) !16 files
      St_List=(/0.27,0.275,0.28,0.282,0.284,0.286,0.288,0.29,0.295,0.298,0.3,0.302,0.304,0.306,0.308,0.31/)!16 files

      DO 1000 SHEARCNT=1,16 
        B(1)=1. !THIS IS AN ARBITRARY AMPLITUDE correction
	  B(2)=EXPTPHI(SHEARCNT)      !PHIZERO in rad from experiment
	  B(3)=22.837				        !THIS IS THE RADIUS OF THE CYLINDERS !20.3
	  B(4)=431.35			            !THIS IS THE CYLINDER LENGTH (THE KUHN LENGTH) !250
	  B(5)=St_List(SHEARCNT)		!THE EXTENSIONAL STRETCH
	  B(6)=0.14998			        !AN ADdITIVE BACKGROUND TERM; ONLY IMPORTANT AT HIGH Q !0.1511
        B(7)= 8.577                    !VOLUME FRACTION OF MICELLES TIMES THE SCATTERING LENGTH DENISTY SQUARE /10**19 CM-4
	WRITE(6,*)B
	WRITE(6,*)
	WRITE(6,*)

******************************************************************************************************
	NCYL=4				!NCYL IS THE NUMBER OF CYLINDERS
******************************************************************************************************      

      !OPEN(UNIT=4,FILE='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\Midq3m.dat',FORM='FORMATTED',ERR=999,STATUS='OLD')
      OPEN(UNIT=4,FILE='D:\Research\SANSpaper\Code\041421_PTU_NORM_2D\IO\LowQMidQHighQ.dat',FORM='FORMATTED',ERR=999,STATUS='OLD')

c***********THIS LOOP READS IN THE DATA FILE WITH FORMAT QX - QY
C*********** THE MAXIMUM NUMBER Of DATA POINTS THAT CAN BE READ IN IS 20000 (ARBITRARY NUMBER)
          DO 100 LCNT=1,MAXDAT
          READ(4,970,END=960)QX1,QY1
          QX(LCNT)=QX1
          QY(LCNT)=QY1
 !         WRITE (6,*)'QX,QY',QX(LCNT),QY(LCNT)
  970     FORMAT(2E16.7)   
  100     CONTINUE
      
  960 NARRAY=LCNT-1
      CLOSE (UNIT=4)  
      
      DO 121 I=1,NARRAY
      Q(I)=SQRT(QX(I)*QX(I)+QY(I)*QY(I))
      COSPSI(I)=QX(I)/Q(I)
      SINPSI(I)=QY(I)/Q(I)
 121  CONTINUE

      
 	    KK=1
	    STEP=1 !THIS USES EVERY STEP'TH POINT IN Q FOR THE CALCULATION

    	    DO 800 I=1,NARRAY,STEP
		QXNEW(KK)=QX(I)
          QYNEW(KK)=QY(I)
          QI=Q(I)
          CS=COSPSI(I)
          SS=SINPSI(I)
************************************************************CALL THE CALCULATION ROUTINES**************
		CALL CSY_MULTI(QI,CS,SS,CSYQYL,NCYL)
*******************************************************************************************************

          YOUT(KK)=CSYQYL
		KK=KK+1
  800     CONTINUE
      
   
      OPEN(UNIT=4,FILE=TRIM(FILENAME(SHEARCNT)),ERR=999)

		DO 200 II=1,KK-1
		WRITE(4,980)QXNEW(II),QYNEW(II),YOUT(II)
  980	    FORMAT(6E15.7)   
  200	    CONTINUE

      CLOSE(UNIT=4)
1000  CONTINUE
      GOTO 998
  999 WRITE(6,*)'ERROR'
  998 CALL DATE_AND_TIME(TIME1,TIME2,TIME3,TVALS2)
      WRITE(6,*)'START DATE=',TVALS1(3),'HOUR=',TVALS1(5),'MINUTE=',TVALS1(6),'SEC=',TVALS1(7)
      WRITE(6,*)'END DATE=',TVALS2(3),'HOUR=',TVALS2(5),'MINUTE=',TVALS2(6),'SEC=',TVALS2(7)
      END
