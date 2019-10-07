PROGRAM hw1

IMPLICIT NONE
REAL (KIND = 8), ALLOCATABLE :: A(:,:), b(:), x(:), y(:),bo(:)
REAL (KIND = 8) ::  xe(2), x1(2), x2(2) &
& , r1(2), r2(2), sol1(2), sol2(2), &
&  Ao(2,2),Ai(1600,1600), Ain(1600,1600), bi(1600), ci(1600), Ai2(2,2)
REAL(KIND = 8) ::  x3(1600), Aio(1600,1600), bio(1600), be(1600),&
& dx(1600)
INTEGER :: n,i,j

! Setting numbers/allocating/reading in
n = 2
ALLOCATE(A(n,n), b(n), y(n), x(n))
A(1,1) = 0.7073725D0
A(1,2) = 0.5204556D0
A(2,1) = 0.8158208D0
A(2,2) = 0.6002474D0
Ao = A
open(12,file = 'A.txt')
read(12,*) Ai
close(12)
Ai = TRANSPOSE(Ai)
Aio = Ai
open(12,file = 'b.txt')
read(12,*) bi
close(12)
bio = bi
open(12,file = 'c.txt')
read(12,*) ci
close(12)
b(1) = 0.1869169D0
b(2) = 0.2155734D0
bo = b
xe(1) = 1D0
xe(2) = -1D0
x1(1) = 0.9999999D0
x1(2) =-1.0000001D0
x2(1) = 0.4073666D0
x2(2) =-0.1945277D0

! Calculating solutions with the given vectors
sol1 = matmul(A,x1)
sol2 = matmul(A,x2)

r1 = b - sol1
r2 = b - sol2

print *, '=========Problem 1========='
print *, 'Residual for x1: ',r1
print *, 'Residual for x2: ',r2

! Here comes the LU decomposition of A
CALL LUdcmp(A,n,b)
CALL LUsolve(A,n,b,x)
print *, 'LUdcmp solution: ', x
print *, 'LUdcmp error:    ', x - xe
print *, 'LUdcmp residual: ', bo - matmul(Ao,x)

! Getting the condition number by multiplying norm by norm of inverse
A = Ao
call matrixinv(A,Ai2,n)
print *, 'Condition number:',norm2(A)*norm2(Ai2) 

print *, '=========Problem 2========='
! Get inverse and mltiply by b
CALL matrixinv(Ai,Ain,1600)
Ai = Aio
print *, 'Condition number:',norm2(Ai)*norm2(Ain) 
x3 = matmul(Ain,bio)
print *, 'GJ error:        ',sum(abs(x3-ci))

! LU dcmp and solve
CALL LUdcmp(Ai,1600,bi)
CALL LUsolve(Ai,1600,bi,x3)

print *, 'LUdcmp error:    ', sum(abs(x3-ci))

! Iterative improvement
be = matmul(Aio,x3) - bio
CALL LUsolve(Ai,1600,be,dx)
x3 = x3 + dx
print *, 'Error after iter:',sum(abs(x3-ci))

END PROGRAM

!=======================================================================
! Returns pivoted LU decomp and pivoted solutoin vector

SUBROUTINE LUdcmp(A,n,b)
INTEGER n,i,j,k, irow
REAL(KIND = 8) A(n,n), sm, big, b(n), tmp(n)
LOGICAL piv

!   pivot first?
!   Column Loop
    DO j = 1,n
        piv = .false.
        DO i =j,n
!           Get top left element
            big = A(j,j)
!           Are other elements in this column bigger than the diagonal?
            IF (abs(A(i,j)).gt.abs(big)) THEN
                irow = i
                big = A(i,j)
                piv = .true.
            ENDIF
        ENDDO

!       Did we find an element bigger than the diagonal? Swap if so
        IF (piv) THEN
            tmp = A(j,:)
            A(j,:) = A(irow,:)
            A(irow,:) = tmp
            tmp(1) = b(j)
            b(j) = b(irow)
            b(irow) = tmp(1)
        ENDIF
    ENDDO

!   Loop over all columns
    DO j = 1, n
!       Loop rows
        DO i = 1, n
            sm = 0D0
    !       Calculate beta
            IF (i.le.j) THEN
    !           First row shenaningans (in first, beta = A)
                IF (i.eq.1)then
                     A(i,j) = A(i,j)
                    CYCLE
                endif
                DO k = 1,i-1
    !               Dealing with alpha = 1 on diag
                    IF (i.eq.k) THEN
                        sm = sm + A(k,j)
                    ELSE
                        sm = sm + A(i,k)*A(k,j)
                    ENDIF
                ENDDO
                A(i,j) = A(i,j) - sm
    !       Calculate alpha
            ELSE
    !           In first column alpa = 1/Beta_jj
                IF(j.eq.1) THEN
                    A(i,j) = A(i,j)/A(j,j)
                    CYCLE
                ENDIF
                DO k=1,j-1
                    sm = sm + A(i,k)*A(k,j)
                ENDDO
                A(i,j) = 1/A(j,j)*(A(i,j)-sm)
            END IF
        ENDDO 
    ENDDO
END SUBROUTINE

!=======================================================================
! Solves LU dcmp of A

SUBROUTINE LUsolve(A,n,b,x)
INTEGER :: n,i,j
REAL(KIND = 8) ::A(n,n), x(n), b(n), y(n), sm

y(1) = b(1)
DO i =2,n
    sm = 0d0
    DO j = 1,i-1
!       Deal with alph diagonals
        IF (i.eq.j) THEN
            sm = sm + y(j)
        ELSE
            sm = sm + A(i,j)*y(j)
        ENDIF
    ENDDO
    y(i) = b(i) - sm
ENDDO

x(n) = y(n)/(A(n,n))

i=n-1
DO WHILE (i.ge.1)
    sm = 0D0
    DO j = i+1,n
        sm = sm + A(i,j)*x(j)
    ENDDO
    x(i) = 1/A(i,i)*(y(i)-sm)
    i = i-1
ENDDO

END SUBROUTINE LUsolve

!=======================================================================
! Returns inverse matrix using Gauss Jordan

SUBROUTINE Matrixinv(A,Ain,n)
INTEGER :: i,j,k,l,n,irow
REAL(KIND = 8):: big,A(n,n),Ain(n,n),tmp

! Identity matrix
DO i = 1,n
    DO j = 1,n
        Ain(i,j) = 0.0D0
    ENDDO
    Ain(i,i) = 1.0D0
ENDDO

! Loop over columns
DO i = 1,n
    big = A(i,i)
    DO j = i,n
        IF (A(j,i).gt.big) THEN
            big = A(j,i)
            irow = j
        ENDIF
    ENDDO

!   Swap rows
    IF (big.gt.A(i,i)) THEN
        DO k = 1,n
            tmp = A(i,k)
            A(i,k) = A(irow,k)
            A(irow,k) = tmp
            tmp = Ain(i,k) 
            Ain(i,k) = Ain(irow,k)
            Ain(irow,k) = tmp
        ENDDO
    ENDIF

!   Divide all entries in row i from A(i,j) by the value A(i,i)
    tmp = A(i,i)
    DO j = 1,n
        A(i,j) = A(i,j)/tmp
        Ain(i,j) = Ain(i,j)/tmp
    ENDDO

!   Zero out
    DO j = i+1,n
        tmp = A(j,i)
        DO k = 1,n
            A(j,k) = A(j,k) - tmp*A(i,k)
            Ain(j,k) = Ain(j,k) - tmp*Ain(i,k)
        ENDDO
    ENDDO
ENDDO

! Subtract rows
DO i = 1,n-1
    DO j = i+1,n
        tmp = A(i,j)
        DO l = 1,n
            A(i,l) = A(i,l)-tmp*A(j,l)
            Ain(i,l) = Ain(i,l)-tmp*Ain(j,l)
        ENDDO
    ENDDO
ENDDO
END SUBROUTINE Matrixinv

