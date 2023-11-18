      SUBROUTINE mcint(NDIM,X,NCOMP,RESUL)
		IMPLICIT NONE
		double precision X(ndim)
		double precision RESUL(1),R0,SIGMA,A,FUNCNUM
        INTEGER NCOMP,NDIM
		COMMON /NORMS/ SIGMA,R0,FUNCNUM      

		A=X(1)**FUNCNUM
	  RESUL(1:NCOMP)=A*A*EXP(-(X(1)-R0)*(X(1)-R0)/2.0D0/SIGMA/SIGMA)
	  RETURN
	  END

