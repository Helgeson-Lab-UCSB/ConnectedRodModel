       	SUBROUTINE	NORMALCALCS(R01,SIGMA1,FUNCNUM1,FUNC) 
    
        IMPLICIT NONE
	    double precision accuracy
	    DOUBLE PRECISION PI
        DOUBLE PRECISION RESUL1(1),SIGMA,R0,FUNC,FUNCNUM
        DOUBLE PRECISION SIGMA1,R01,FUNCNUM1
        integer maxdim,NCOMP,MAXITER,NCALL,NEVAL,J,NDIM
        parameter (maxdim = 1)
        double precision upper(maxdim)
        common /ubound/ upper
        double precision lower(maxdim)
        common /lbound/ lower
        COMMON /NORMS/ SIGMA, R0,FUNCNUM

	    external ScaledIntegrand2

        PI=4.D0*DATAN(1.D0)
	    maxiter=5000
	    ncall=100000
        SIGMA=SIGMA1
        R0=R01
        FUNCNUM=FUNCNUM1
        NCOMP=1. 
        ndim=1
	    accuracy=.00003
	    LOWER(1)=R0-4.0D0*SIGMA
        IF(LOWER(1).LT.0.0D0)LOWER(1)=0.0D0
        UPPER(1)=R0+4.0D0*SIGMA
!**********************THIS VERSION STILL CALLS THE NON-VECTORIZED VERSION OF VEGAS SINCE IT IS A SIMPLER INTEGRATION 
        CALL vegasF2(ScaledIntegrand2,resul1,accuracy, ndim,ncall, maxiter, neval,NCOMP)  
!********************** 
        FUNC=RESUL1(1)
    
      END
