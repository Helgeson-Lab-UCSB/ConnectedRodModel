C*******************************************************************************
C   THIS IS A SUBROUTINE TO SET UP THE PARAMETERS FOR THE CALCULATION AND INTEGRATION OVER THE CYCLINDER ANGLES 
C   GREG SMITH 1/10/2021
C*******************************************************************************

	SUBROUTINE	CSY_MULTI(QZ,COSPSI,SINPSI,CSYQYL,NCYL,NQ,NCALL)
	IMPLICIT NONE					!
	REAL*8 CSYQYL(NQ),B,P,RE,CQ(NQ)
	integer ndim, ncall, maxiter, NGRID,i
      INTEGER MAXPAR,L,LL,N,M,K
	INTEGER NCYL,NQ
      INTEGER NCOMP
      PARAMETER (NCOMP=1)
      REAL*8 upper((NCYL*2)+1)                        !ADD ONE FOR THE DISTRIBUTION IN THE RADIUS
      REAL*8 lower((NCYL*2)+1)
      REAL*8 UPPER_TR(NCYL*2)
      REAL*8 LOWER_TR(NCYL*2)
	COMMON  /BLK1/B(20),P(20),RE,N,M,K
	COMMON /STUFF/RCYL,LCYL,SIGMA_R
	COMMON /PIBLOCK/ PI
      COMMON  /LBLS/LABEL(20)     
	REAL*8 COSPSI(NQ),SINPSI(NQ)
	REAL*8 RCYL,LCYL,Z,SIGMA_R
	REAL*8 PQ(NQ)
      REAL*8 QZ(nq)
	REAL*8 ACC, RESUL(NQ),PI,RESUL_TR(NQ)
	!EXTERNAL SCALEDINTEGRAND
      EXTERNAL Ncylffac_SECT_FIT

	CHARACTER *10	MODID
	PARAMETER		(MAXPAR=20)				
      CHARACTER *25	LABEL 

      DATA        LABEL/'AMPL','CYLINDER RADIUS','CYLINDER LENGTH',
     1	        'BACKGROUND','VOL FRACTION x SLD','S(Q=0)','AVG CYLINDER NUMBER(LC)', 'NUMBER CYL MAX',
     2            'SIGMA RADIUS','NCALL_START','THIN ROD LENGTH',9*' '/

						!
C...
	PI=4.0D0*DATAN(1.0D0)      
      WRITE(6,*)'NCYL=',NCYL
C
C*****************************************************

C			 FUNCTION fivecyl_FIT_PROB_SECT_funct(Q,PARAM,NPARAM)

! PARAM IS A VECTOR WITH [RCYL,LCYL,LHYDRO,AMPLITUDE]
! THIS HAS AN ARBITRAY MULIPLIER OF 1000.  tHE RESULT IS THAT THE mc INTEGRATION IS ALMOST 1/3 THE TIME BUT GIVES THE SAME ANSWER!!!!!!!!!!!!!


! THIS IS THE SANS INTENSITY CALCULATION FOR NCYL JOINTED CYLINDERS
! EACH CYLINDER CAN ROTATE FREELY ABOUT THE AXES AND Q IS FIXED
! I(Q)=INTEGRAL(FFxFF*)
! GAMMA IS THE ANGLE BETWEEN THE Q VECTOR AND THE AXIS OF EACH CYLINDER AND IS DETERMINED BY THETA(I) AND PHI(I)
! S(I) IS THE VECTOR FROM THE ORIGIN TO THE CENTER OF THE CYLINDERS 
!! THETA(I) IS THE ROTATION OF A CYLINDER ABOUT THE Y AXIS AND PHI(I) IS THE ROTATION ABOUT THE Z-AXIS
!THE FLOW DIRECTION IS ALONG THE X-AXIS (SEE Hayter and Penfold)
! Q LIES IN THE X-Y PLANE AND IS DEFINED BY MAGNITUDE (Q) AND ROTATION WITH RESPECT TO X-AXIS (PSI)  - THIS VERSION IS FOR THE FLOW-GRADIENT CASE
	

		RCYL=B(2)                   
		LCYL=B(3)                  
          SIGMA_R=B(9)
		Z=DBLE(0.)


c	These vectors have the lower and upper angular bounds for the integrals.
c     E.G. For 5 cylinders there are 10 integrals- ONE FOR EACH THETA AND PHI
c     and 5 integrals of the psi angles for the cylinders and 5 for the phi angles of the cylinders  

          DO 60 L=2,(2*NCYL)+1
              LOWER(L)=Z
   60     CONTINUE
          LOWER(1)=B(2)-(4.0D0*B(9))                   !INTEGRATE FROM B(3)-4SIGMA TO B(3) + 4SIGMA
          IF (LOWER(1).LT.0.0D0)LOWER(1)=0.0D0
          DO 61 LL=2,NCYL+1
              UPPER(LL)=PI
              UPPER(LL+NCYL)=2.0D0*PI
   61     CONTINUE
          UPPER(1)=B(2)+(4.0D0*B(9)) 
          LOWER_TR=LOWER(2:((2*NCYL)+1))
          UPPER_TR=UPPER(2:((2*NCYL)+1))
          
          
          
          
	acc=0.01D0              !ACC IS THE ACCURACY OF THE INTEGRAL CALCULATION
                              !WAS .001----FEB 23 2002 INCREASED TO SEE IF I CAN GET MORE CYLINDERS- GREATLY AFFECTS COMPUTATION TIME
	resul=0.0D0
	ndim=((2*NCYL)+1)
      maxiter=5	            !SUGGESTED IN lEPAGE PAPER	
	NGRID=50                !SUGGESTED IN lAPAGE PAPER

      
!******************************************************CALL THE MONTE CARLO INTEGRATION ROUTINES VEGAS********

      CALL VEGASF_OLD(Ncylffac_SECT_FIT,RESUL,ACC,NCALL,MAXITER,NGRID,NCYL,NDIM,UPPER,LOWER,COSPSI,SINPSI,QZ,NQ)      


!********************************************************************************************************      

		PQ=RESUL/NCYL/NCYL/1000.      !using normailization factor for the sum of the form factors square for ncyl cylinders which is 1/(ncyl**2) SO Pq GOES TO 1 AT Q=0
                                      
     
       CSYQYL=PQ			        


	RETURN

	END












