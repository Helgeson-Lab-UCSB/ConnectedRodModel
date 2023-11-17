C*******************************************************************************
C   THIS IS A SUBROUTINE TO SET UP THE PARAMETERS FOR THE CALCULATION AND INTEGRATION OVER THE CYCLINDER ANGLES 
C   GREG SMITH 1/10/2021
C*******************************************************************************

	SUBROUTINE	CSY_MULTI(QZ,COSPSI,SINPSI,CSYQYL,NCYL)
						!
	DOUBLE PRECISION CSYQYL
	integer ndim, ncall, maxiter, neval
	INTEGER NCYL
      INTEGER NCOMP
      PARAMETER (NCOMP=1)
      DOUBLE PRECISION upper(NCYL*2)
      DOUBLE PRECISION lower(NCYL*2)
	COMMON		/BLK1/B(20),P(20),RE,N,M,K
	COMMON /STUFF/QPASS,PHIO,RCYL,LCYL,STRETCH
	COMMON /PIBLOCK/ PI
	DOUBLE PRECISION SRATEIN,COSPSI,SINPSI
	DOUBLE PRECISION QPASS,PHIO,RCYL,LCYL,STRETCH,Z
	DOUBLE PRECISION PQ
      DOUBLE PRECISION QZ
	double precision ACC,RESUL(NCOMP),PI
	EXTERNAL ScaledIntegrand

	CHARACTER *10	MODID
	PARAMETER		(MAXPAR=20)				


						!
C...
	PI=4.0D0*DATAN(1.0D0)         

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
	

		PHIO=DBLE(B(2)) !phi0 from experiment
		RCYL=DBLE(B(3))                   
		LCYL=DBLE(B(4))                  
		STRETCH=DBLE(B(5))
		Z=DBLE(0.)


c	These vectors have the lower and upper angular bounds for the integrals.
c     E.G. For 5 cylinders there are 10 integrals- ONE FOR EACH THETA AND PHI
c     and 5 integrals of the psi angles for the cylinders and 5 for the phi angles of the cylinders  

          DO 60 L=1,(2*NCYL)
              LOWER(L)=Z
   60     CONTINUE
          DO 61 LL=1,NCYL
              UPPER(LL)=PI
              UPPER(LL+NCYL)=2.*PI
   61     CONTINUE
!      WRITE(6,*)'LOWER,UPPER', LOWER,UPPER
 !     STOP

		  QPASS=QZ       ! THIS IS ONLY A SINGLE VALUE FOR QZ

	acc=0.05D0      
	resul=0.0D0
	ndim=(2*NCYL)
	ncall=10000				!originally 10000
	maxiter=100000		
	neval=0

!******************************************************CALL THE MONTE CARLO INTEGRATION ROUTINES VEGAS********

	CALL VEGASf(ScaledIntegrand,resul,acc,NDIM,NCALL,MAXITER,NEVAL,NCYL,UPPER,LOWER,COSPSI,SINPSI)
      
!********************************************************************************************************      
      PQ=RESUL(1)/NCYL/NCYL/1000.     !using normailization factor for the sum of the form factors square for ncyl cylinders which is 1/(ncyl**2) SO Pq GOES TO 1 AT Q=0
          !PQ=B(7)*PQ                      !SCALED TO PHI*(DELTA BETA)**2 WHICH IS PARAMETER B(7)IN 10-19 CM-1
          !PQ=PI*RCYL*RCYL*LCYL*NCYL*PQ*1.E-5   !SCALED TO VOLUME OF A "MICELLE" WHICH IS THE VOLUME OF A CYLINDER X NUMBER OF CYLINDER ; THE LAST FACTOR CONVERTS THE DIMESNIONS IN ANGSTROMS TO CM-1

       !CSYQYL=B(1)*PQ +B(6)			                    !B(1) is an arbitrary scaling factor and b(6) is an additive background
      CSYQYL=PQ
	RETURN

	END












