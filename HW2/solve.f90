MODULE solve

CONTAINS

! Solves a tridiagonal matrix according to numerical recipes
SUBROUTINE tridiag(A,x,b,n)
    INTEGER, INTENT(IN) :: n
    REAL(KIND = 8), INTENT(IN) :: A(n,n), b(n)
    REAL(KIND = 8), INTENT(OUT) :: x(n)
    REAL(KIND = 8) :: gam(n), bet
    INTEGER :: i

    bet = 0D0
    gam = 0D0
    x = 0D0
    bet = A(1,1)
    x(1) = b(1)/bet
    DO i =2,n
        gam(i) = A(i-1,i)/bet
        bet = A(i,i) - A(i,i-1)*gam(i)
        x(i) = (b(i)-A(i,i-1)*x(i-1))/bet
    ENDDO

    i = n-1
    DO WHILE(i.ge.1)
        x(i) = x(i) - gam(i+1)*x(i+1)
        i = i-1
    ENDDO

END SUBROUTINE tridiag

END MODULE solve