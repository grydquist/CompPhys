PROGRAM hw2
USE interp
implicit none

TYPE(interpType) :: interpolator
REAL(KIND=8), ALLOCATABLE :: xin(:), yin(:)
REAL(KIND=8) ::  diff
INTEGER i,m

! Set number of points to use
m = 21

! Allocate x's/y's to be given to interpolator
ALLOCATE(xin(m),yin(m))

! Evaluating those x's and y's
xin(1) = -1D0
yin(1) = 1D0/(1D0+25D0*xin(1)**2)
diff = 2D0/(m-1D0)
DO i = 2,m
    xin(i) = xin(i-1) + diff
    yin(i) = 1D0/(1D0+25D0*xin(i)**2)
ENDDO

! Initializa the interpolation structure with 300 points
interpolator = interpType(m,xin,yin,300)

! Use a polynomial interpolation and write to file
CALL interpolator%poly
CALL interpolator%write('p21.txt')

! Use a cubic spline interpolation and write to file
CALL interpolator%cspline
CALL interpolator%write('c21.txt')


END PROGRAM hw2