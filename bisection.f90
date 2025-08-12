#include "../parallel.h"
!*************************************************************************
 PROGRAM bisection
!*************************************************************************
   use io
   implicit none
   double precision :: r,th,z
   integer :: n
   type (coll) :: c11,c21,c31
   type (coll) :: c12,c22,c32
   _loop_km_vars


   call mpi_precompute()
   call par_precompute()
   call mes_precompute()
   call var_precompute()
   call tra_precompute()
   call tim_precompute()
   call vel_precompute()
   call  io_precompute()
   
   if(mpi_sze/=1) stop 'set _Np 1'

   tim_t  = 0d0
   tim_dt = 0.001d0
 
   print*, 'Enter input statefile #1:'
   read(*,*) io_statefile   
    
   
   call io_load_state()
   call var_coll_copy(vel_ur, c11)
   call var_coll_copy(vel_ut, c21)
   call var_coll_copy(vel_uz, c31)
   
   print*, 'Enter input statefile #2:'
   read(*,*) io_statefile     
   
   call io_load_state()
   call var_coll_copy(vel_ur, c12)
   call var_coll_copy(vel_ut, c22)
   call var_coll_copy(vel_uz, c32)  
   
   
   
   call var_coll_init(vel_ur)
   call var_coll_init(vel_ut)
   call var_coll_init(vel_uz)
 
   
! **** SET IC IN REAL SPACE ****
!   do n = 1, i_N
!      do m = 0, i_Th-1
!         do k = 0, i_Z-1
!            r  = mes_D%r(n,1)
!            th = (2d0*d_PI/dble(i_Th)) * m
!            z  = (2d0*d_PI/(d_alpha*i_Z)) * k
!            vel_r%Re(k,m,n) =  (1-r**2)/10
!            vel_t%Re(k,m,n) =  0 
!            vel_z%Re(k,m,n) =  0 
!         end do   
!      end do
!   end do
!   call tra_phys2coll(vel_r,vel_ur, vel_t,vel_ut, vel_z,vel_uz)

! **** SET IC IN COLLOCATED SPACE ****
! k and m are set by the _loop_km macro,
! see parallel.h

   _loop_km_begin
      do n = 1, i_N
         r = mes_D%r(n,1)
         vel_ur%Re(n,nh) = (c11%Re(n,nh) + c12%Re(n,nh))/2
         vel_ur%Im(n,nh) = (c11%Im(n,nh) + c12%Im(n,nh))/2 
         vel_ut%Re(n,nh) = (c21%Re(n,nh) + c22%Re(n,nh))/2 
         vel_ut%Im(n,nh) = (c21%Im(n,nh) + c22%Im(n,nh))/2
         vel_uz%Re(n,nh) = (c31%Re(n,nh) + c32%Re(n,nh))/2 
         vel_uz%Im(n,nh) = (c31%Im(n,nh) + c32%Im(n,nh))/2 
      end do
   _loop_km_end

! **** SAVE TO state1000.cdf.dat ****
   call io_save_state()


!*************************************************************************
 END PROGRAM bisection
!*************************************************************************

