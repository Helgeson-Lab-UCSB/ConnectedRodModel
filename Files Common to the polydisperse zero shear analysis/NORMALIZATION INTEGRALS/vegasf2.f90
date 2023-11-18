! vegas.F
! VEGAS Monte Carlo integration of a vector either with ordinary random
! numbers or with quasi-random numbers
! based on code by M. Martinez, J. Illana, J. Bossert, and A. Vicini
! this file is part of FormCalc
! last modified 14 Sep 01 th

! NDIM is the maximum number of dimensions the integrand may have.
! Note: if your Fortran compiler allows, replace NDIM by ndim
! (i.e. #define NDIM ndim) to dimension the arrays dynamically.
!$define NDIM 3

! NCOMP is the number of components of the integrand vector.
!$define NCOMP 2

! If DEBUG is defined, the error is printed out after each iteration
! so that one can tune the parameters for a particular integral.
!#define DEBUG

! RNG determines which random-number generator is used:
! - 1 uses Fortran's ran(...) function,
! - 2 uses a Faure quasi-random sequence,
! - 3 uses a Sobol quasi-random sequence (default).
!      #define RNG 3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! vegas integrates a vector in the ndim-dimensional unit hypercube
!! using the Monte-Carlo method. The NCOMP-dimensional integrand is
!! invoked via the subroutine func(ndim, x, NCOMP, RESUL). After
!! sampling ncall points, the grid is refined. This is called one
!! iteration. The iterations loop terminates if the relative error is
!! below accuracy, or after maxiter iterations.

!! Caution: the refinement of the grid is done only with respect to
!! the highest component of the integrand vector.

	  subroutine vegasf2(func,resul,accuracy, ndim,ncall, maxiter, neval,NCOMP)
	  implicit none
	  integer ncomp
!      parameter (ncomp =1)
	  double precision resul(2000)
	  double precision accuracy
	  integer ndim, ncall, maxiter, neval
	  external func

	  double precision fun(NCOMP), sfun(NCOMP), sfun2(NCOMP)
	  double precision sint(NCOMP), sweight(NCOMP)
	  double precision fun2, sint2, weight
	  double precision r, dr, xo, xn, err
	  integer iter, call, dim, grid, g, c

	  integer ngrid
	  parameter (ngrid = 100)

	  double precision xi(ngrid, NDIM), d(ngrid, NDIM)
	  double precision x(NDIM), imp(ngrid), tmp(ngrid - 1)
	  integer pos(NDIM)
!	  ndim=3
!      write(6,*)'accu,ndim,ncall,maxi,nev=',accuracy,ndim,ncall,maxiter,neval    

	  !call inirandom(maxiter*ncall, ndim)
      CALL RANDOM_SEED
	  iter = 0
	  neval = 0
	  do c = 1, NCOMP
	    sint(c) = 0
	    sweight(c) = 0
	  enddo
	  sint2 = 0

! define the initial distribution of intervals
!$OMP PARALLEL
!$OMP DO        
        do grid = 1, ngrid
	    r = dble(grid)/ngrid
	    do dim = 1, ndim
	      xi(grid, dim) = r
	    enddo
        enddo
!$OMP END DO
!$OMP END PARALLEL
        
        
! iterations loop
1	  continue
	  iter = iter + 1

! initialize iteration variables
  	  do c = 1, NCOMP
	    sfun(c) = 0
	    sfun2(c) = 0
        enddo

        
        do dim = 1, ndim
	    do grid = 1, ngrid
	      d(grid, dim) = 0
	    enddo
        enddo

        
        
	  do call = 1, ncall
	    weight = 1D0/ncall

! compute the point position
	!    call getrandom(x)
        CALL RANDOM_NUMBER(X)
          

          do dim = 1, ndim
	      r = x(dim)*ngrid + 1
	      grid = int(r)
	      xo = 0
	      if(grid .gt. 1) xo = xi(grid - 1, dim)
	      xn = xi(grid, dim) - xo
	      x(dim) = xo + (r - grid)*xn
	      pos(dim) = grid
	      weight = weight*xn*ngrid
          enddo

          
! compute the function value
!		write(6,*)'ndim,x',ndim,x
!		write(6,*)'function value=', fun									!********
	    call func(ndim, x, NCOMP, fun)

          
          do c = 1, NCOMP
	      fun2 = fun(c)*weight
	      sfun(c) = sfun(c) + fun2
	      fun2 = fun2**2
	      sfun2(c) = sfun2(c) + fun2
          enddo
         
          

          do dim = 1, ndim
	      d(pos(dim), dim) = d(pos(dim), dim) + fun2
          enddo
         
          
        enddo
       

        
        neval = neval + ncall
!      write(6,*)'im here 1'																!*****
! compute the integral and error values
!	write(6,*)'ncomp=',ncomp
	  do c = 1, NCOMP
	    fun2 = sfun(c)**2
	    r = sfun2(c)*ncall - fun2
!		write(6,*)'r=',r								!******
!		write(6,*)'sfun=',sfun								!******
!		write(6,*)'sfun2=',sfun2								!******
!		write(6,*)'fun2=',fun2								!******



	    if(r .ne. 0) then
	      weight = fun2/abs(r)*(ncall - 1)
	      sweight(c) = sweight(c) + weight
	      sint(c) = sint(c) + sfun(c)*weight
	    endif
	    if(sweight(c) .eq. 0) then
	      resul(c) = 0
	    else
	      RESUL(c) = sint(c)/sweight(c)
	    endif
	  enddo
	  sint2 = sint2 + fun2
      err = sqrt(sint2/(sweight(NCOMP)*iter))/abs(RESUL(NCOMP))
!      write(6,*)'im here 3'														!***
!ifdef DEBUG
!	  print *, "iteration ", iter, "  error ", err
!	  print *,"resul=", resul
!endif
	  if(RESUL(NCOMP) .eq. 0 .or. err .lt. accuracy) return
	  if(iter .gt. maxiter) then
	    print *,"Warning: VEGAS failed to reach the desired accuracy."
	    print *, "Remaining relative error: ", err
	    return
        endif

! redefine the grid (importance sampling)
! - smooth the f^2 value stored for each interval


        
!$OMP PARALLEL
        do dim = 1, ndim
	    xo = d(1, dim)
	    xn = d(2, dim)
	    d(1, dim) = .5D0*(xo + xn)
	    x(dim) = d(1, dim)

!$OMP DO           
	    do grid = 2, ngrid - 1
	      r = xo + xn
	      xo = xn
	      xn = d(grid + 1, dim)
	      d(grid, dim) = (r + xn)/3D0
	      x(dim) = x(dim) + d(grid, dim)
	    enddo
!$OMP END DO

	    d(ngrid, dim) = .5D0*(xo + xn)
	    x(dim) = x(dim) + d(ngrid, dim)
        enddo
!$OMP END PARALLEL

        
   
! - compute the importance function of each interval
	  do dim = 1, ndim
	    r = 0

           
          do grid = 1, ngrid
	      imp(grid) = 0
	      if(d(grid, dim) .gt. 0) then
	        xo = x(dim)/d(grid, dim)
	        imp(grid) = ((xo - 1)/xo/log(xo))**1.5D0
	      endif
	      r = r + imp(grid)
	    enddo
 
          
          
	    r = r/ngrid
!      write(6,*)'im here 2'															     	   !***
! - redefine the size of each interval
	    dr = 0
	    xn = 0
	    g = 0
	    do grid = 1, ngrid - 1
	      do while(dr .lt. r)
	        g = g + 1
	        dr = dr + imp(g)
	        xo = xn
	        xn = xi(g, dim)
	      enddo
	      dr = dr - r
	      tmp(grid) = xn - (xn - xo)*dr/imp(g)
	    enddo
	    do grid = 1, ngrid - 1
	      xi(grid, dim) = tmp(grid)
	    enddo
	    xi(ngrid, dim) = 1
	  enddo
 
        
	  goto 1

      end


