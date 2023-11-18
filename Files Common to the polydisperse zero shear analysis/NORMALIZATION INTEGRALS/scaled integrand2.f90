        SUBROUTINE ScaledIntegrand2(ndim, x, ncomp, resul)
        implicit none
        integer ndim, ncomp
        double precision x(ndim), resul(1)

        integer maxdim
        parameter (maxdim = 1)

        double precision upper(maxdim)
        common /ubound/ upper

        double precision lower(maxdim)
        common /lbound/ lower

        integer dim, comp
        double precision range, jacobian, scaledx(maxdim)
        jacobian = 1
        do dim = 1, ndim
          range = upper(dim) - lower(dim)
          jacobian = jacobian*range
          scaledx(dim) = lower(dim) + x(dim)*range
        enddo
		CALL MCINT(ndim, scaledx, ncomp, resul)

        do comp = 1, ncomp
          resul(comp) = resul(comp)*jacobian
        enddo
        end
