        SUBROUTINE ScaledIntegrand(ndim, x, result,NCYL,UPPER,LOWER,COSPSI,SINPSI)

!**************** ROUTINE TO RESCALE THE X-VALUES FROM VEGAS FROM A HYPERCUBE THAT HAS UNIT DIMENSIONS TO REAL ANGLES IN THE INTEGRAND RANGES    
        integer ndim
        double precision x(ndim), result(1)
        integer NCYL
        double precision upper((NCYL*2))
        double precision lower((NCYL*2))
        integer dim, comp
        double precision range, jacobian, scaledx(NCYL*2),COSPSI,SINPSI

        jacobian = 1.0D0
        do dim = 1, ndim
          range = upper(dim) - lower(dim)
          jacobian = jacobian*range
          scaledx(dim) = lower(dim) + x(dim)*range
        enddo

		CALL Ncylffac_SECT_FIT(SCALEDX,RESULT(1),NCYL,COSPSI,SINPSI)


        result(1) = result(1)*jacobian
        end
