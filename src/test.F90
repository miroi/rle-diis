
module equations

integer, parameter :: k  = 4, m = 3, verb=3, max_iter=30
real*8, allocatable :: a(:),a_s(:), resid_vect(:),a_m(:),sigma(:),bta(:)
real*8, allocatable :: T(:,:),t_s(:),tta(:),ttbta(:),ttbta_m(:),t_mp1(:)
real*8, allocatable :: B(:,:),B_s(:,:),BT(:,:),BTBT(:,:),TTBTBT(:,:),TTB(:,:)
integer, allocatable :: ipiv(:)
real*8 :: resid
real*8, external :: dnrm2

contains

subroutine allocate_vars
 allocate(a(k),a_s(k),t_s(k),a_m(k),sigma(k),bta(k),ttbta(m),ttbta_m(m))
 allocate(B(k,k),B_s(k,k),TTB(m,k),t_mp1(m))
 allocate(T(k,m),BT(k,m),BTBT(k,m),TTBTBT(m,m))
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

  if (verb >= 4) then
     print *,'roots of a+Bt=0, t:',a_m
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

subroutine solve_diis
  integer :: info
! Solve  T^T.B^T.a + T^T.B^T.B.T.sigma = 0 for sigma vector
! the B matrix and the a vector are given

! get vector B^T.a = bta
  call dgemv('t', k, k, 1.0D0, B_s, k, a_s, 1, +0.0D0, bta, 1)
  print *,'bta=B^T.a:',bta

! get T^T.bta = ttbta
  call dgemv('t', k, m, 1.0D0, T, k, bta, 1, +0.0D0, ttbta, 1)
  print *,'T^T.bta = ttbta:',ttbta

! get B.T into BT
  call dgemm('n', 'n', k, k, m, 1.0D0, B_s, k, T, k, 0.0D0, BT, k)
  print *,'B(k,k).T(k,m)=BT(k,m)=',BT

! get B^T . BT into BTBT
  call dgemm('t', 'n', k, k, m, 1.0D0, B_s, k, BT, k, 0.0D0, BTBT, k)
  print *,'B^T(k,k) . BT(k,m) = BTBT(k,m):',BTBT

! get T^T.BTBT into TTBTBT
  call dgemm('t', 'n', k, k, m, 1.0D0, T, k, BTBT, k, 0.0D0, TTBTBT, m)
  print *,'T^T(m,k).BTBT(k,m) = TTBTBT(m,m)',TTBTBT

! we left with ttbta + TTBTBT.sigma = 0 to get sigma
  ttbta_m = -ttbta
  call dgesv( m, 1, TTBTBT , m, ipiv, ttbta_m, m, info )
  if (info.gt.0) print *,"error in dgesv routine !"
  sigma = ttbta_m
  print *,'sigma vector:',sigma

end subroutine


subroutine solve_rle
! solve  T^T.a + T^T.B.T.sigma = 0 
! to get sigma (m) vector

! get T^T.a into ttbta
  call dgemv('t', k, m, 1.0D0, T, k, a_s, 1, +0.0D0, ttbta, 1)
  print *,'T^T(m,k).a(k) = ttbta(m):',ttbta

! get T^T.B into TTB
  call dgemm('t', 'n', m, k, k, 1.0D0, T, k, B_s, k, 0.0D0, TTB, m)
  print *,'T^T(m,k) .B(k,k) =  TTB(m,k) :',TTB

! get TTB.T into TTBTBT
  call dgemm('n', 'n', m, m, k, 1.0D0, TTB, k, T, k, 0.0D0, TTBTBT, m)
  print *,'TTB(m,k).T(k,m) = TTBTBT(m,m):',TTBTBT

! we left with ttbta(m) + TTBTBT(m,m).sigma(m) = 0 to get sigma
  ttbta_m = -ttbta
  call dgesv( m, 1, TTBTBT , m, ipiv, ttbta_m, m, info )
  if (info.gt.0) print *,"error in dgesv routine !"
  sigma = ttbta_m
  print *,'sigma vector:',sigma


end subroutine

subroutine update_sigma_t
! using calculated sigma vector and the T matrix get new estimate of sigma.T, T[m+1]
! get  T(k,m).sigma(m) =  t_mp1(m)
  call dgemv('n', k, m, 1.0D0, T, k, sigma, 1, +0.0D0, t_mp1, 1)
  print *,'new t[m+1] obtained as T(k,m).sigma(m) = t_mp1(m):',t_mp1

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
iter = iter + 1

call solve_a_Bt_lineq

if (iter > k) then 
  !call solve_diis ! ... some error in routine ...
  call solve_rle
  call update_sigma_t
endif

call update_a_B
call get_resid
print *,'iter=',iter,' resid=', resid

do_iter = resid > 0.00002 .and. iter < max_iter

enddo








End Program
