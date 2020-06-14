!************************************************************************                            
!*Cavity Flow in a Porous Medium of Binary Mixture with the Soret Effect*
!************************************************************************
!Parameters:
!	  Ra - Rayleigh number
!	  Le - Lewise number
!     m	 - porosity
!  soret - Soret parameter
!   ampl - vibration amplitude
!  omega - vibration frequincy

program sq_convection
!use GVLRG_int


!---------------------variable declaration(beg)---------------------------
implicit none
integer, parameter :: N=50    ! must be odd

real(kind(1d0)) Ra, Ra1, Ra2, Le, m, soret, ampl, omega, ampl_psi           ! manag. parameters
real(kind(1d0)) eps_poisson, eps_threshold     ! accuracy
real(kind(1d0)) psi_new(0:N,0:N), psi(0:N,0:N), T(0:N,0:N), eta(0:N,0:N)            ! arrays

real(kind(1d0)) time_for_count, h, dt, time, relativ, relax, period    
real(kind(1d0)) max_psi, max_psi_new, max_psi_old , increment_new, increment, max_eta_new, freq , time_in, time_out ! inc
integer i, j, k, inc_counter, time_counter, N_half
logical first_time

integer  Max_counter_poisson

real(kind(1d0)) PI/3.14159265/
real(kind(1d0)) a_sq, b_sq, c_sq, x1, x2, x3, y1, y2, y3, R_sub


character(len=13) filename_psi, filename_T, filename_eta, filename_C
character(len=2) filename_counter
integer period_counter, period_Num_steps/10/
real(kind(1d0)) time_step_period 
!---------------------variable declaration(end)---------------------------

!--------------------Reading data from file(beg)--------------------------
open(1,file='input1.dat')
read(1,*) Ra 
read(1,*) Le 
read(1,*) m	 
read(1,*) soret
read(1,*) ampl
read(1,*) omega
read(1,*) time_for_count
read(1,*) eps_poisson
read(1,*) Max_counter_poisson
read(1,*) eps_threshold
read(1,*) relativ
read(1,*) relax
!---------------------Reading data from file(end)--------------------------

open(5,file='out_psimax_t.dat')
open(9,file='eps_Ra_sub.dat')

h=1./N		 	! grid step
N_half=N/2

!R_sub=180. 
!do while (soret<0.1)
period=0.5

!----------------------Initial and boundary conditions(beg)----------------
do i=1, N-1
 do j=1, N-1 
    psi(i,j)=10.*sin(PI*i*h)*sin(PI*j*h)*cos(PI*sqrt(2.)*(i*h-0.5))       
    T(i,j)=10.*sin(PI*i*h)*sin(PI*j*h)*sin(PI*sqrt(2.)*(i*h-0.5))/sqrt(Ra)
    eta(i,j)= 0. !-0.01*cos(PI*i*h)                                   
 enddo
enddo
eta(20,20)=0.01 ! 25

do i=0, N
     psi(i,0)=0.
     psi(i,N)=0.
     psi(0,i)=0.
     psi(N,i)=0.
     T(i,0)=0.
     T(i,N)=0.
     T(0,i)=0.
     T(N,i)=0.
     eta(i,0)=eta(i,1)
     eta(i,N)=eta(i,N-1)
     eta(0,i)=eta(1,i)
     eta(N,i)=eta(N-1,i)
enddo
!----------------------Initial and boundary conditions(end)----------------

!Ra=R_sub+25.
time_in=0.
time_out=period
call fields_comp_SIP (N, Ra, Le, m, soret, ampl, omega, time_in, time_out, ampl_psi, psi, T, eta, '.')
x1=ampl_psi
y1=Ra
write(*,'(x, F6.2, x, F8.5 )') Ra, ampl_psi
!write(9,'(x, F6.2, x, F8.5 )') Ra, ampl_psi

!Ra=R_sub+15.
!time_in=0.
!time_out=period
!call fields_comp_SIP (Ra, Le, m, soret, ampl, omega, psi, T, eta, time_in, time_out, eps_poisson, Max_counter_poisson, relativ, relax, ampl_psi)
!x2=ampl_psi
!y2=Ra
!write(*,'(x, F6.2, x, F8.5 )') Ra, ampl_psi
!
!Ra=R_sub+2.
!time_in=0.
!time_out=period
!call fields_comp_SIP (Ra, Le, m, soret, ampl, omega, psi, T, eta, time_in, time_out, eps_poisson, Max_counter_poisson, relativ, relax, ampl_psi)
!x3=ampl_psi
!y3=Ra
!write(*,'(x, F6.2, x, F8.5 )') Ra, ampl_psi
!
!!write(9,'(x, F6.2, x, F8.5 )') Ra, ampl_psi
!!write(*,'(x, F6.2, x, F8.5 )') Ra, ampl_psi
!
!a_sq = (-x3 * y1 + x1 * y3 + x2 * y1 - x2 * y3 + x3 * y2 - x1 * y2) / (-x2 * x3 ** 2 + x2 * x1 ** 2 + x1 * x3 ** 2 - x3 * x1 ** 2 + x3 * x2 ** 2 - x1 * x2 ** 2)
!b_sq = -(x1 ** 2 * y3 - x1 ** 2 * y2 - x2 ** 2 * y3 + y2 * x3 ** 2 + y1 * x2 ** 2 - y1 * x3 ** 2) / (-x2 * x3 ** 2 + x2 * x1 ** 2 + x1 * x3 ** 2 - x3 * x1 ** 2 + x3 * x2 ** 2 - x1 * x2 ** 2)
!c_sq = (x1 ** 2 * x2 * y3 - x1 ** 2 * x3 * y2 - x2 ** 2 * x1 * y3 + y2 * x1 * x3 ** 2 - y1 * x2 * x3 ** 2 + x2 ** 2 * x3 * y1) / (-x2* x3 ** 2 + x2 * x1 ** 2 + x1 * x3 ** 2 - x3 * x1 ** 2 + x3 * x2 ** 2 - x1 * x2 ** 2)
!R_sub=(4.*a_sq*c_sq - b_sq*b_sq)/(4.*a_sq)
!
!write(9,'(x, F6.3, x, F8.3, x, F8.3 )') soret,  R_sub, -b_sq/(2.*a_sq)
!write(*,'(x, F6.3, x, F8.3, x, F8.3 )') soret,  R_sub, -b_sq/(2.*a_sq)
!
!soret=soret+0.01
!!Ra=Ra-1.
!enddo

!period_counter=1
!time_in=period
!time_step_period=period/period_Num_steps
!do while(period_counter<=period_Num_steps+1)
!
!!-------------------------Write data to file(beg)-------------------------- 
! write(filename_counter,'(I2)') period_counter
! filename_psi='out_psi'//filename_counter//'.dat'
! filename_T='out_T'//filename_counter//'.dat'
! filename_eta='out_eta'//filename_counter//'.dat'
! filename_C='out_C'//filename_counter//'.dat'
! 
! open(period_counter+100,file=filename_psi)
! open(period_counter+200,file=filename_T)
! open(period_counter+300,file=filename_eta)
! open(period_counter+400,file=filename_C)
! 
! do i=0, N
!  do j=0, N
!   write(period_counter+100,'(x,f10.5,x,f10.5,x,f10.5)') i*h,  j*h, psi(i,j) 
!   write(period_counter+200,'(x,f10.5,x,f10.5,x,f10.5)') i*h,  j*h, T(i,j) 
!   write(period_counter+300,'(x,f10.5,x,f10.5,x,f10.5)') i*h,  j*h, eta(i,j) 
!   write(period_counter+400,'(x,f10.5,x,f10.5,x,f10.5)') i*h,  j*h, eta(i,j)+soret*T(i,j) 
!  enddo
! enddo
!!-------------------------Write data to file(end)--------------------------
!
!time_out=time_out+time_step_period
!call fields_comp_SIP (Ra, Le, m, soret, ampl, omega, psi, T, eta, time_in, time_out, eps_poisson, Max_counter_poisson, relativ, relax)
!write(*,*) time_in, time_out 
!time_in=time_out
!period_counter=period_counter+1
! 
!enddo

!-------------------------Write data to file(beg)--------------------------  
 open(2,file='out_psi.dat')
 open(3,file='out_T.dat')
 open(7,file='out_eta.dat')
 open(8,file='out_C.dat')
 do i=0, N
  do j=0, N
   write(2,'(x,f10.5,x,f10.5,x,f10.5)') i*h,  j*h, psi(i,j) 
   write(3,'(x,f10.5,x,f10.5,x,f10.5)') i*h,  j*h, T(i,j) 
   write(7,'(x,f10.5,x,f10.5,x,f10.5)') i*h,  j*h, eta(i,j) 
   write(8,'(x,f10.5,x,f10.5,x,f10.5)') i*h,  j*h, eta(i,j)+soret*T(i,j) 
  enddo
 enddo
!-------------------------Write data to file(end)--------------------------

end program sq_convection
