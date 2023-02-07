!------- OMP scheme:
!====================================================
!nthreads =16 A double do loop is organized in 16
! blocks of equal size with addresses [i1:i2][j1:j2]
! defined, for each thread by the brp(rank)%i1,i2,j1,j2
! struct. Each thread reads on (private) pt(:) copy of data
!         and writes on the (shared) ac(:) output data
!   cpu time : a scalar execution takes 40-44 cpu sec
!   cpu time : OMP execution using 16 threads takes
!    3 seconds for the 4 threads executing the first do loop (see below)
!    9 seconds for the 12 threads executing the second do loop
!    total 12 sec
!    NB each of the sixteen double loops has equal size [i1:i2]x[j1:j2]
!==================================
!    Credo che molto dipenda dalle clausole private/shared naturalmente
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program OpenMP_test1
use OMP_LIB

subroutine omp_test_grav(pt,bpr,ac)

type(vec4),intent(in) :: pt(:)
type(pair),intent(in) :: bpr(:)
type(vec4),intent(inout) :: ac(:)
integer :: rank,p,q,i1,i2,j1,j2
real :: phi0(0:3),grv0(0:3),xc(1:3)

!$OMP parallel default(private) shared(ac,bpr)
   rank=1+OMP_get_thread_num()
   i1=bpr(rank)%ib(1)
   i2=bpr(rank)%ib(2)
   j1=bpr(rank)%jb(1)
   j2=bpr(rank)%jb(2)
   if(i1==j1)then
    do q=i1,i2-1
     grv0=0.0
     do p=q+1,i2
      xc(1:3)=pt(q)%pos(1:3)-pt(p)%pos(1:3)
      call two_body(xc,phi0)
      phi0=pt(p)%pos(0)*phi0
      grv0=grv0+phi0
      ac(p)%pos(1:3)=ac(p)%pos(1:3)+phi0(1:3)
     end do
    end do
    !This loop are correctely executed in parallel
   else
    do q=i1,i2
     grv0=0.0
     do p=j1,j2
      xc(1:3)=pt(q)%pos(1:3)-pt(p)%pos(1:3)
      call two_body(xc,phi0)
      phi0=pt(p)%pos(0)*phi0
      grv0=grv0+phi0
      ac(p)%pos(1:3)=ac(p)%pos(1:3)+phi0(1:3)
     end do
    end do
    !write(*,*)rank,i1,i2,j1,j2
   endif
!$OMP end parallel

end subroutine omp_test_grav
