		subroutine Ncylffac_SECT_FIT (X,IQ,NCYL,COSPSI,SINPSI,Q,NQ)
!**************************
!**********************************************************************************************
    
!    THIS IS FOR THE GRADIENT-VORTICITY DIRECTION
    
!**********************************************************************************************    
!  THIS IS THE ROUTINE THAT HAS MOST OF THE CALULATIONS FOR THE FUNCTIONS NEEDED FOR 
!  N CYLINDER CALCULATIONS


! THIS HAS AN ARBITRAY MULIPLIER OF 1000.  THE RESULT IS THAT THE MONTE CARLO INTEGRATION 
! IS ALMOST 1/3 THE TIME BUT GIVES THE SAME ANSWER!!!!!!!!!!!!!

! This is a N cylinder form factor: the probablity distribution for each cylinder is NOT in here but could be added- it can be edited to include 
! both orientations of flow for a couette cell as noted

		IMPLICIT NONE
		COMMON /STUFF/RCYL,LCYL,SIGMA_R
		integer NCYL,J,II,KK,JJ,NQ
		REAL*8 GAMMAI(NCYL)													!,GAMMAI_MINUS(NCYL)		add this for couette geometry where you go through the cell twice
		REAL*8 X((2*NCYL)+1)
		REAL*8 THETA(NCYL), PHI(NCYL)
		REAL*8 COSPSI, SINPSI, PI
		REAL*8 COSTHETA(NCYL), SINTHETA(NCYL),COSPHI(NCYL),SINPHI(NCYL)
		REAL*8 SX(NCYL),SY(NCYL)											!,SZ_MINUS(NCYL),SX_MINUS(NCYL)		add this for couette geometry where you go through the cell twice
		REAL*8 RCYL,LCYL,FQ(2,NQ),RCYL1,SIGMA_R
        REAL*8 Q(NQ)
        REAL*8 QDOTS(NCYL,NQ)												!,QDOTS_MINUS(NCYL)			add this for couette geometry where you go through the cell twice
		REAL*8 F1(NCYL,NQ)													!,F1_MINUS(NCYL)							add this for couette geometry where you go through the cell twice
		REAL*8 PROB0(NCYL)
		REAL*8 TEM1(NQ), TEM2(NQ),IQ(NQ),IQ1(NQ)							!, TEM1_MINUS,TEM2_MINUS				add this for couette geometry where you go through the cell twice

        PI=4.0D0*DATAN(1.0D0)  
        
        RCYL1=X(1)       
		THETA=X(2:NCYL+1)
		PHI=X(NCYL+2:((NCYL*2)+1))
		COSTHETA=COS(THETA)
		SINTHETA=SIN(THETA)
		COSPHI=COS(PHI)
		SINPHI=SIN(PHI)
		GAMMAI=ACOS(COSPSI*COSPHI*SINTHETA+SINPSI*SINTHETA*SINPHI)   
																			!GAMMAI_MINUS=ACOS(-COSPSI*COSPHI*SINTHETA+SINPSI*COSTHETA)   add this for couette geometry where 
																			!you go through the cell twice
		SX(1)=0.0D0
		SY(1)=0.0D0
        																	!SX_MINUS(1)=DBLE(0.)    add this for couette geometry where you go through the cell twice
																			!SZ_MINUS(1)=DBLE(0.)    add this for couette geometry where you go through the cell twice
		QDOTS=0.0D0        
																			!QDOTS_MINUS(1)=DBLE(0.)        add this for couette geometry where you go through the cell twice
			DO 20 J=2,NCYL
			SX(J)=SX(J-1)+(COSPHI(J)*SINTHETA(J)+COSPHI(J-1)*SINTHETA(J-1))	
            SY(J)=SY(J-1)+(SINPHI(J)*SINTHETA(J)+SINPHI(J-1)*SINTHETA(J-1))
 			QDOTS(J,1:NQ)= Q(1:NQ)*LCYL/2.*(COSPSI*SX(J)+SINPSI*SY(J))            
                                                                            !QDOTS_MINUS(J)=QPASS*LCYL/2.*(-COSPSI*SX(J)+SINPSI*SZ(J))	add this for couette geometry where you go through the cell twice
20          CONTINUE

		F1=0.0D0
		DO 10 II=1,NCYL
		call CFF2(GAMMAI(II),Q,RCYL1,LCYL,FQ(1,1:NQ),NQ)       
																			!call CFF2(GAMMAI_MINUS(II),QPASS,RCYL,LCYL,FQ(2))  add this for couette geometry where you go through the cell twice       
        F1(II,1:NQ)= FQ(1,1:NQ)      
																			!F1_MINUS(II)=FQ(2) add this for couette geometry where you go through the cell twice
10      CONTINUE

        PROB0=1.0D0/4.0D0/PI        !-------------PROBABILITY SET TO A UNIFORM VALUES FOR ZERO SHEAR
             
		TEM1=F1(1,1:NQ)
																			!		TEM1_MINUS=F1_MINUS(1)  add this for couette geometry where you go through the cell twice
		TEM2=0.0D0
																			!		TEM2_MINUS=DBLE(0.)     add this for couette geometry where you go through the cell twice
		DO 40 KK=2,NCYL
		TEM1(1:NQ)=TEM1(1:NQ)+F1(KK,1:NQ)*COS(QDOTS(KK,1:NQ))
		TEM2(1:NQ)=TEM2(1:NQ)+F1(KK,1:NQ)*SIN(QDOTS(KK,1:NQ))
																			!TEM1_MINUS=TEM1_MINUS+F1_MINUS(KK)*COS(QDOTS_MINUS(KK))		add this for couette geometry where you go through the cell twice
																			!TEM2_MINUS=TEM2_MINUS+F1_MINUS(KK)*SIN(QDOTS_MINUS(KK))		add this for couette geometry where you go through the cell twice
40		CONTINUE

		IQ=(TEM1*TEM1+TEM2*TEM2)											!*2.  !+(TEM1_MINUS*TEM1_MINUS+TEM2_MINUS*TEM2_MINUS)  add this for couette geometry where you go through the cell twice									
		IQ1=1.0D0
			DO 50 JJ=1,NCYL
			IQ1(1:NQ)=IQ1(1:NQ)*SINTHETA(JJ)*PROB0(JJ)
50			CONTINUE
		IQ=IQ*IQ1
		IQ=IQ*1000.
		IQ=IQ*RCYL1*RCYL1*DEXP(-(RCYL1-RCYL)*(RCYL1-RCYL)/2.0D0/SIGMA_R/SIGMA_R)
		end