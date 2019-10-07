MODULE interp
    USE solve

!   Set up interp structure. This will interpolate a given set of points
!   on the interval of those points onto the number of points specified
    TYPE interpType
!       Number of points in (degree of polynomial-1) and out
        INTEGER :: m,pts
!       x and y points in
        REAL(KIND=8), ALLOCATABLE :: xin(:), yin(:)
!       x and y points (interp'd) out
        REAL(KIND=8), ALLOCATABLE :: xout(:), yout(:)
    
!       Procedures contained in structure
        CONTAINS
!           Polynomial interpolation procedure
            PROCEDURE :: poly => polyinterp
!           cubic spline interpolation procedure
            PROCEDURE :: cspline => csplineinterp
!           procedure that writes txt file with xout and yout
            PROCEDURE :: write => writeinterp
    END TYPE interpType

!   Initialization procedure
    INTERFACE interpType
        PROCEDURE :: newint
    END INTERFACE

    CONTAINS

!=======================================================================
!   Initialize structure
    FUNCTION newint(m,xin,yin,pts) RESULT(interp)
        TYPE(interpType) :: interp
        INTEGER,INTENT(IN) :: m,pts
        REAL(KIND=8), INTENT(IN) :: xin(m), yin(m)
        REAL(KIND=8) :: diff

!       Allocating x and y out vectors with given arguments
        ALLOCATE(interp%xout(pts),interp%yout(pts))

!       Some temporary variables to make things easier
        interp%xin = xin
        interp%yin = yin
        interp%m = m
        interp%pts = pts

!       Calculating x values of points that will be given out
        diff = (xin(m) - xin(1))/(pts-1D0)
        interp%xout(1) = xin(1)
        DO i = 2,pts
            interp%xout(i) = interp%xout(i-1) + diff
        ENDDO
    END FUNCTION newint
    
!=======================================================================
!   Polynomial interpolation
    SUBROUTINE polyinterp(interp)
        CLASS(interpType),INTENT(INOUT) :: interp
        REAL(KIND=8) :: P(interp%m,interp%m),xtmp,xin(interp%m),&
        & yin(interp%m)
        INTEGER :: i,j,k

!       Temporary variables
        xin = interp%xin
        yin = interp%yin
        m = interp%m

!       Polynomial of degree zero
        P(1,:) = yin

!       Loop over all points I want to interpolate with
        DO k = 1,interp%pts
!           Temporary point (to be interpolated with)
            xtmp = interp%xout(k)
!           Loop over the columns in the tableau
            DO j = 2,m
!               Loop over the rows in the jth column of the tableau
                DO i = 1,m-j+1
!                   Calculate polynomial of this column/row of tableau
                    P(j,i) = ((xtmp - xin(i    ))*P(j-1,i+1) &
                    &      -  (xtmp - xin(i+j-1))*P(j-1,i  )) &
                    &      /  (xin(i+j-1)-xin(i))
                ENDDO
            ENDDO
!           y value is the far right element of the tableau
            interp%yout(k) = P(m,1)
        ENDDO
    END SUBROUTINE polyinterp

!=======================================================================
!   Cubic spline interpolation
    SUBROUTINE csplineinterp(interp)
        CLASS(interpType),INTENT(INOUT) :: interp
        REAL(KIND=8) :: y2(interp%m),xin(interp%m), &
        & Am(interp%m,interp%m),bm(interp%m), yin(interp%m), A, B, C, D&
        & ,xm, xp, xo, ym, yp, y2m,y2p
        INTEGER :: i,j,m

!       Temporary variables
        m = interp%m
        xin = interp%xin
        yin = interp%yin

!       First we need to solve for y2. We have the matrix problem:
!       Am y2 = bm (since A and b are used later)
        Am = 0D0
        bm = 0D0
        Am(1,1) = 1D0
        Am(m,m) = 1D0
        DO i = 2,m-1
!           Calculating the Am and bm terms in the tridiagonal matrix
            Am(i,i-1) = (xin(i)   - xin(i-1))/6D0
            Am(i,i)   = (xin(i+1) - xin(i-1))/3D0
            Am(i,i+1)   = (xin(i+1) - xin(i)  )/6D0
            bm(i) = (yin(i+1) - yin(i))/(xin(i+1)-xin(i)) &
            &    - (yin(i)-yin(i-1))/(xin(i)-xin(i-1))
        ENDDO

!       Solve tridiagonal Am y2 = bm
        CALL tridiag(Am,y2,bm,m)

!       xm = lower point in current interval
!       xp = upper point in current interval
!       Same for y2
        xm = xin(1)
        xp = xin(2)
        ym = yin(1)
        yp = yin(2)
        y2m = y2(1)
        y2p = y2(2)
        j = 2

!       Loop over points to interpolate
        DO i = 1, interp%pts
            xo = interp%xout(i)
!           Have we moved to the next interval? Move xp, xm, etc if so
!           Not a perfect way to do this, but works for the problem
            IF (xo.gt.xp) THEN
                xm = xin(j)
                xp = xin(j+1)
                ym = yin(j)
                yp = yin(j+1)
                y2m = y2(j)
                y2p = y2(j+1)
                j=j+1
            ENDIF

!           Calculate constants in equation
            A = (xp - xo)/(xp-xm)
            B = 1 - A
            C = 1D0/6D0*(A**3D0 - A)*(xp-xm)**2D0
            D = 1D0/6D0*(B**3D0 - B)*(xp-xm)**2D0

!           Use constants to get interpolated y
            interp%yout(i) = A*ym + B*yp + C*y2m + D*y2p
        ENDDO

    END SUBROUTINE

!=======================================================================
!   Write out the data
    SUBROUTINE writeinterp(interp,name)
        CLASS(interpType), INTENT(IN) :: interp
        CHARACTER(len = 7), INTENT(IN) :: name
        INTEGER :: i

!       Write into file name
        OPEN(10, file = name)
        DO i = 1,interp%pts
            write(10,'((f10.3))') (/interp%xout(i),interp%yout(i)/)
        ENDDO
        CLOSE(10)

    END SUBROUTINE writeinterp

END MODULE