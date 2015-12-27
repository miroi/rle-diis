
module equations

integer, parameter :: k  = 10, m = 5, verb=3, max_iter=100
real*8, allocatable :: a(:),a_s(:), resid_vect(:),a_m(:)
real*8, allocatable :: T(:,:),t_s(:)
real*8, allocatable :: B(:,:),B_s(:,:)
integer, allocatable :: ipiv(:)
real*8 :: resid
real*8, external :: dnrm2

contains

subroutine allocate_vars
 allocate(a(k),a_s(k),t_s(k),a_m(k))
 allocate(B(k,k),B_s(k,k))
 allocate(T(k,m))
 allocate(ipiv(k))
end subroutine

subroutine init_t
  integer :: i,j
  do j=1,m
  do i=1,k
    T(i,j)=0
  enddo
  enddo
end subroutine

subroutine update_a_B
  integer :: i
  do i=1,k
    !a(i)=-3.0d0*dfloat(i)+3.0d0 + 0.5d0*T(i,m)
     !a(i)=-log( dfloat(i) + 3.0d0 ) + 0.005d0*T(i,m)
     a(i)=-log( dfloat(i) + 3.0d0 ) 
    a_s(i) = a(i) ! store a(i)
    do j=1,k
      !B(i,j)=dfloat(i)-0.08d0*dfloat(j)+0.05d0*T(i,m)-0.02d0*T(j,m)
      !B(i,j)=log( dfloat(i)+dfloat(j+2)) / log(dfloat(i+1)*dfloat(j+2))+0.05d0*T(i,m)-0.02d0*T(j,m)

      !B(i,j)=log( dfloat(i)+dfloat(j+2)) / log(dfloat(i+1)*dfloat(j+2)) + log(0.005d0*T(i,m)*T(i,m)+0.001d0) + log(0.002d0*T(j,m)*T(j,m)+0.02d0)

     !OK this !
     ! B(i,j)=log( dfloat(i)+dfloat(j+2)) / log(dfloat(i+1)*dfloat(j+2)) + log(0.0005d0*T(i,m)*T(i,m)+0.001d0) & 
     !  - log(0.0002d0*T(j,m)*T(j,m)+0.002d0)

      B(i,j)=log( dfloat(i)+dfloat(j+2)) / log(dfloat(i+1)*dfloat(j+2)) +  & 
             log(0.0005d0*T(i,m)*T(i,m)+0.001d0)/log(0.0002d0*T(j,m)*T(j,m)+0.002d0) 

     ! B(i,j)=log( dfloat(i)+dfloat(j+2)) / log(dfloat(i+1)*dfloat(j+2)) + 0.001*log (0.00005d0*T(i,m)*T(i,m) + 0.0001)  < --- OK !!!

     ! B(i,j)=log( dfloat(i)+dfloat(j+2)) / log(dfloat(i+1)*dfloat(j+2))  ! this is fine ...

      B_s(i,j)=B(i,j) ! store B(i,j)
    enddo
  enddo
end subroutine

subroutine solve_a_Bt_lineq
!
! solves system of linear equations,  a(t) + B(t).t = 0
!
! vector a(t) and matrix B(t) are allocated and filled, get the t vector - solutions
!
  integer :: info,i,j

  if (verb >= 3) then
    ! print *,'B:',b
    ! print *,'a:',a
    ! print *,'-a:',-a
  endif

  !  B.t = -a
  a_m=-a; 
  !print *,"a:",a
  !print *,"-a:",-a
  call dgesv( k, 1, B , k, ipiv, a_m, k , info )
  if (info.gt.0) stop "error in dgesv routine !"

  if (verb >= 3) then
     !print *,'roots of a+Bt=0, t:',a_m
  endif

! fill roots : -a to T(k,m), shift columns to the left
 do j=2,m
 do i=1,k
   T(i,j-1)=T(i,j)
 enddo
 enddo

! 1st column: a_m to T(:,m)
 do i=1,k
   T(i,m) = a_m(i); t_s(i)=T(i,m)
 enddo
 !print *,'t_s',t_s

end subroutine

subroutine get_resid
! calculates the norm of the residuum vector resid =  |a+B.t|

  !print *,"get_resid:  get B_s.t_s + a_s into a_s"

 ! get vector B_s.t_s + a_s into a_s
  call dgemv('n', k, k, 1.0D0, B_s, k, t_s, 1, +1.0D0, a_s, 1)
  !print *,'zero matrix, B.t+a',a_s

   ! the a_s should be zero vector
   ! get euklid.norm |a+B.t| - must be zero
   !resid=dnrm2(k,a_s,1)/(dfloat(k))
   resid=dnrm2(k,a_s,1)

   !print *,dnrm2(k,a_s,1)

end subroutine

end module equations


Program Test_RLE_DIIS
use equations
logical :: do_iter=.true.

integer :: iter = 0

print *,"Hello !"
call allocate_vars
call init_t
call update_a_B

do while ( do_iter) 
call solve_a_Bt_lineq
call update_a_B

call get_resid
print *,'iter=',iter,' resid=', resid

iter = iter + 1

do_iter = resid > 0.00002 .and. iter < max_iter
enddo








End Program
