		subroutine p_th_ph1a_Ncyl(theta,phi,prob,NCYL)
! This is the Hayter Penfold expression for shear
! 
!not elog below is the natural log
!P_THETA_PHI=FLTARR(N_ELEMENTS(THETA),N_ELEMENTS(PHI))
!		IMPLICIT NONE
		COMMON /PIBLOCK/ PI
        COMMON /STUFF/QPASS,PHIO,RCYL,LCYL,STRETCH
 		integer NCYL
!        parameter (NCYL= 7)			!NUMBER OF CYLINDERS


		DOUBLE PRECISION theta(NCYL),phi(NCYL)
		double precision A(NCYL),B(NCYL)
		DOUBLE PRECISION phizero,R
		DOUBLE PRECISION QPASS,PHIO,RCYL,LCYL,STRETCH

		double precision prob(NCYL)
		double precision PI
!		write(6,*)'theta,phi=',theta,phi
        phizero=PHIO
		r=rcyl
		
!		WRITE(6,*)'theta=',theta
!		WRITE(6,*)'phi=',phi
!		WRITE(6,*)'shear_rate,L_HYD,R=',shear_rate,L_HYD,R

!!!! use phizero directly in the calculation of OPDF, don't use shear rate
		!a=(1.-cos(2.*phizero))*(1.+sin(theta)*sin(theta)*cos(2.*phizero))**1.5
		!b=4.*pi*(1.-sin(theta)*sin(theta)*cos(2.*phizero)*cos(2.*(phi-phizero)))**2.
!		WRITE(6,*)'A,B=',A,B
        IF (stretch.NE.0) THEN
		prob=exp((3.*stretch-stretch**3.)/(1.-stretch**2.)*cos(phi-phizero)*sin(theta))
        prob=prob/(4.*pi*(1-stretch**2.)/(3.*stretch-stretch**3.)*sinh((3.*stretch-stretch**3.)/(1.-stretch**2.)))
        ELSE
            prob=1./(4.*pi) !prob expression has 0 on the denominator if stretch=0
        ENDIF
        ! already normalized
        !prob=1./4.*pi !for equilibrium calculation only
		end