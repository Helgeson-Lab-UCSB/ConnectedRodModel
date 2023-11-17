		subroutine CFF2 (GAM,Q,R,L,fq)
! THIS IS THE FORM FACTOR FOR A CYLINDER
!		IMPLICIT NONE
!		Use  PORTLIB
		COMMON /PIBLOCK/ PI
		double precision TERM1,TERM2,FQ
		double precision GAM

		double precision Q,R,L
		double precision PI
!		REAL BESJ1
!       WRITE(6,*)"GAM,Q,R,L",GAM,Q,R,L


!		write(6,*)'besj1( 
		TERM1=Q*R*SIN(GAM)
		TERM2=(Q*L/2.)*COS(GAM)
 !       write(6,*)'q,r,gam,sin(gam),term1,term2=',q,r,gam,sin(gam),term1,term2




		   IF (GAM.EQ.0) THEN 
		   FQ=SIN(Q*L/2.)/(Q*L/2.)
!		   write(6,*)'GAM EQ 0:gam,fq=',gam,fq
		   ELSEIF (GAM.EQ.PI/2.) THEN 
		   FQ=2.*BESSEL_J1(Q*R)/(Q*R) !changed from DBESJ1 to BESSEL_J1
!		   write(6,*)'GAM EQ PI/2: gam,fq=',gam,fq
		   ELSE
		   FQ=2.*(BESSEL_J1(TERM1)/TERM1)*(SIN(TERM2)/TERM2) !changed from DBESJ1 to BESSEL_J1
!		   write(6,*)'gam,fq,term1,BESSEL_J1(term1)=',gam,fq,term1,BESSEL_J1(term1) !changed from DBESJ1 to BESSEL_J1
		   ENDIF
!      WRITE(6,*)"FQ=",FQ

!RETURN,FQ
!		CFF2=FQ
		END