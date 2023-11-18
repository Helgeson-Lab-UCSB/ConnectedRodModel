		subroutine CFF2 (GAM,Q,R,L,fq,NQ)
! THIS ROUTINE CALCULATES THE FORM FACTOR OF A SINGLE CYLINDER
        USE IFPORT		
        IMPLICIT NONE
!		Use  PORTLIB
		COMMON /PIBLOCK/ PI
        INTEGER NQ,I
		REAL*8 TERM1(NQ),TERM2(NQ),FQ(NQ)
		REAL*8 GAM

		REAL*8 Q(NQ),R,L
		REAL*8 PI
        
		TERM1=Q*R*SIN(GAM)
		TERM2=(Q*L/2.)*COS(GAM)
		   IF (GAM.EQ.0) THEN 
		   FQ=SIN(Q*L/2.)/(Q*L/2.)
           ELSEIF (GAM.EQ.PI/2.) THEN 
				DO 10 I=1,NQ 
				FQ(I)=2.0D0*DBESJ1(Q(I)*R)/(Q(I)*R)
10				CONTINUE
           ELSE
				DO 20 I=1,NQ   
				FQ(I)=2.0D0*(DBESJ1(TERM1(I))/TERM1(I))*(SIN(TERM2(I))/TERM2(I))
20				CONTINUE
		   ENDIF
		END