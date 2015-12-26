
module equations_params

integer, parameter :: k  = 10, m = 5
real*8, allocatable :: a(:)
real*8, allocatable :: t(:)
real*8, allocatable :: B(:,:)

contains

subroutine allocate_a_B
 allocate(a(k))
 allocate(B(k,k))
 allocate(t(k))
end subroutine

subroutine init_t
  integer :: i
  do i=1,k
    t(i)=0
  enddo
end subroutine

subroutine update_a_B
  integer :: i
  do i=1,k
    a(i)=dfloat(i)+2.0d0 + 0.05d0*t(i)
    do j=1,k
      B(i,j)=dfloat(i)+0.5d0*dfloat(j)+0.05d0*t(i)-0.02d0*t(j)
    enddo
  enddo
end subroutine

subroutine solve_a_Bt_lineq
  integer :: info
  !call dgesv( n, n, a, n, ipiv, b, n , info )
  !if (info.gt.0) stop "error in dgesv routine !"
end subroutine

end module equations_params


Program Test_RLE_DIIS
use equations_params
integer :: iter = 0
print *,"Hello !"
call allocate_a_B
call init_t

call update_a_B
call solve_a_Bt_lineq

End Program
