subroutine init_fields(N, psi, T, eta)

! =====================================
! Initial and boundary conditions
! =====================================
implicit none
integer, intent(in) :: N    ! must be odd
real*8, intent(out) :: psi(0:N,0:N)
real*8, intent(out) :: T(0:N,0:N)
real*8, intent(out) :: eta(0:N,0:N)
integer i, j, k 
real*8 PI/3.14159265/, h

h=1./N
do i=1, N-1
  do j=1, N-1 
     psi(i,j)=0. !10.*sin(PI*i*h)*sin(PI*j*h)*cos(PI*sqrt(2.)*(i*h-0.5))       
     T(i,j)=0.! 10.*sin(PI*i*h)*sin(PI*j*h)*sin(PI*sqrt(2.)*(i*h-0.5))/sqrt(78.)
     eta(i,j)= 0. !-0.01*cos(PI*i*h)                                   
  enddo
enddo

psi(20,20)=0.01 
T(20,20)=0.01
eta(20,20)=0.01 

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
end

subroutine fields_comp_sip(N, Ra, Le, m, soret, ampl, omega, time_in, time_out, ampl_psi, psi, T, eta, case_path)
! ======================================================================
! Cavity Flow in a Porous Medium of Binary Mixture with the Soret Effect
! SIP: strongly implicite procedure with Stony approximation
! ======================================================================
! Parameters:
!	  Ra - Rayleigh number
!	  Le - Lewise number
!     m	 - porosity
!  soret - Soret parameter
!   ampl - vibration amplitude
!  omega - vibration frequincy
! ======================================================================

!---------------------variable declaration(beg)---------------------------
implicit none
integer, intent(in) :: N
integer i, j, k
real*8, intent(in) :: Ra, Le, m, soret, ampl, omega           ! manag. parameters
real*8, intent(in) :: time_in, time_out
real*8, intent(out) :: ampl_psi
real*8, intent(inout), dimension(0:N, 0:N) :: psi
real*8, intent(inout), dimension(0:N, 0:N) :: T 
real*8, intent(inout), dimension(0:N, 0:N) :: eta 
real*8 psi_new(0:N,0:N)
real*8 h, dt, time, max_psi, max_psi_old, max_eta, max_eta_old
real*8 alpha /0.92/
integer icounter_poisson, out_iter_counter, inner_iter_counter 
logical stac_flg
real*8 Un(0:N,0:N), Ue(0:N,0:N), Lp(0:N,0:N), Ls(0:N,0:N), Lw(0:N,0:N), res(0:N,0:N)
real*8 Ap(0:N,0:N), Ae(0:N,0:N), Aw(0:N,0:N), As(0:N,0:N), An(0:N,0:N), Q(0:N,0:N), resN
real*8 eps_SIP/0.000001/, eps_stac/0.000001/, relativ/1.0/
character (len = *) case_path
!---------------------variable declaration(end)---------------------------

h=1./N    ! grid step
dt=relativ*h*h  ! time step
time=time_in
stac_flg=.true.

open(5, file = case_path//'/'//'out_psimax_t.dat')

do i=0, N
  Un(0,i)=0.; Un(i,0)=0.; Un(N,i)=0.; Un(i,N)=0.; 
  Ue(0,i)=0.; Ue(i,0)=0.; Ue(N,i)=0.; Ue(i,N)=0.;
  Lp(0,i)=0.; Lp(i,0)=0.; Lp(N,i)=0.; Lp(i,N)=0.;
  Lw(0,i)=0.; Lw(i,0)=0.; Lw(N,i)=0.; Lw(i,N)=0.;
  Ls(0,i)=0.; Ls(i,0)=0.; Ls(N,i)=0.; Ls(i,N)=0.;
  res(0,i)=0.; res(i,0)=0.; res(N,i)=0.; res(i,N)=0.; 
enddo

max_psi_old=0.
max_eta_old=0.
stac_flg=.true.

!--------------------------Main cycle on time(beg)-------------------------
do while (time < time_out)
!do while ((time < time_out) .and. stac_flg)

!out_iter_counter=0.
!do while (out_iter_counter<1) 
!out_iter_counter=out_iter_counter+1
!write (*,'(I3, x, F6.3, x, F6.3, x, F6.3, x, I3)') out_iter_counter, psi(20,20), T(20,20), eta(20,20), inner_iter_counter

!------------------------Stream function computing(beg)--------------------
 do j=1, N-1
  do i=1, N-1
   Ap(i,j)=-4.
   As(i,j)=1.
   Ae(i,j)=1.
   An(i,j)=1.
   Aw(i,j)=1.
   Q(i,j)=-h*(Ra+ampl*sin(omega*time))*((1.+soret)*(T(i+1,j)-T(i-1,j))+(eta(i+1,j)-eta(i-1,j)))/2.
  enddo
 enddo

  
 do j=1, N-1
  do i=1, N-1
   Lw(i,j)=Aw(i,j)/(1.+alpha*Un(i-1,j))
   Ls(i,j)=As(i,j)/(1.+alpha*Ue(i,j-1))
   Lp(i,j)=Ap(i,j) + alpha*(Lw(i,j)*Un(i-1,j) + Ls(i,j)*Ue(i,j-1)) - Lw(i,j)*Ue(i-1,j) - Ls(i,j)*Un(i,j-1)
   Un(i,j)=(An(i,j) - alpha*Lw(i,j)*Un(i-1,j))/Lp(i,j)
   Ue(i,j)=(Ae(i,j) - alpha*Ls(i,j)*Ue(i,j-1))/Lp(i,j)
  enddo
 enddo

resN = 1.
inner_iter_counter=0.

do while (resN>eps_SIP)
inner_iter_counter=inner_iter_counter+1
 
 resN=0.
 do j=1, N-1
  do i=1, N-1
   res(i,j)=Q(i,j)-Ae(i,j)*psi(i+1,j)-Aw(i,j)*psi(i-1,j)-An(i,j)*psi(i,j+1)-As(i,j)*psi(i,j-1)-Ap(i,j)*psi(i,j)

   resN=resN+ABS(res(i,j))
   enddo
 enddo
 

 do j=1, N-1
  do i=1, N-1
   res(i,j)=(res(i,j) - Lw(i,j)*res(i-1,j) - Ls(i,j)*res(i,j-1))/Lp(i,j)
  enddo
 enddo


 do j=N-1, 1, -1 
  do i=N-1, 1, -1 
   res(i,j)=res(i,j) - Un(i,j)*res(i,j+1) - Ue(i,j)*res(i+1,j)
   psi(i,j)=psi(i,j)+res(i,j)
  enddo
 enddo
enddo
!------------------------Stream function computing(end)--------------------

!------------------------Temperature computing(beg)------------------------
 do j=1, N-1
  do i=1, N-1
   Ap(i,j)=-(4.*relativ+1.)
   As(i,j)=relativ*(1. - (psi(i+1,j) - psi(i-1,j))/4.)
   Ae(i,j)=relativ*(1. - (psi(i,j+1) - psi(i,j-1))/4.)
   An(i,j)=relativ*(1. + (psi(i+1,j) - psi(i-1,j))/4.)
   Aw(i,j)=relativ*(1. + (psi(i,j+1) - psi(i,j-1))/4.)
   Q(i,j)=-T(i,j)+ relativ*h*(psi(i+1,j)-psi(i-1,j))/2.
  enddo
 enddo

  
 do j=1, N-1
  do i=1, N-1
   Lw(i,j)=Aw(i,j)/(1.+alpha*Un(i-1,j))
   Ls(i,j)=As(i,j)/(1.+alpha*Ue(i,j-1))
   Lp(i,j)=Ap(i,j) + alpha*(Lw(i,j)*Un(i-1,j) + Ls(i,j)*Ue(i,j-1)) - Lw(i,j)*Ue(i-1,j) - Ls(i,j)*Un(i,j-1)
   Un(i,j)=(An(i,j) - alpha*Lw(i,j)*Un(i-1,j))/Lp(i,j)
   Ue(i,j)=(Ae(i,j) - alpha*Ls(i,j)*Ue(i,j-1))/Lp(i,j)
  enddo
 enddo

resN = 1.
inner_iter_counter=0.

do while (resN>eps_SIP)
inner_iter_counter=inner_iter_counter+1
 
 resN=0.
 do j=1, N-1
  do i=1, N-1
   res(i,j)=Q(i,j)-Ae(i,j)*T(i+1,j)-Aw(i,j)*T(i-1,j)-An(i,j)*T(i,j+1)-As(i,j)*T(i,j-1)-Ap(i,j)*T(i,j)

   resN=resN+ABS(res(i,j))
   enddo
 enddo
 

 do j=1, N-1
  do i=1, N-1
   res(i,j)=(res(i,j) - Lw(i,j)*res(i-1,j) - Ls(i,j)*res(i,j-1))/Lp(i,j)
  enddo
 enddo


 do j=N-1, 1, -1 
  do i=N-1, 1, -1 
   res(i,j)=res(i,j) - Un(i,j)*res(i,j+1) - Ue(i,j)*res(i+1,j)
   T(i,j)=T(i,j)+res(i,j)
  enddo
 enddo
enddo
!------------------------Temperature computing(end)------------------------


!------------------------Chemical potential computing(beg)-----------------
 do j=1, N-1
  do i=1, N-1
   Ap(i,j)=-(4.*Le*relativ/m+1.)
   As(i,j)=relativ/m*(Le - (psi(i+1,j) - psi(i-1,j))/4.)
   Ae(i,j)=relativ/m*(Le - (psi(i,j+1) - psi(i,j-1))/4.)
   An(i,j)=relativ/m*(Le + (psi(i+1,j) - psi(i-1,j))/4.)
   Aw(i,j)=relativ/m*(Le + (psi(i,j+1) - psi(i,j-1))/4.)
   Q(i,j)=-eta(i,j)+soret*relativ*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)-4.*T(i,j))+ soret*(1./m-1.)*relativ*h*(psi(i+1,j)-psi(i-1,j))/2.
   Q(i,j)= Q(i,j)+soret*(1./m-1.)*relativ*( (psi(i,j+1)-psi(i,j-1))*(T(i+1,j)-T(i-1,j)) - (psi(i+1,j)-psi(i-1,j))*(T(i,j+1)-T(i,j-1)) )/4.
  enddo
 enddo

  
 do j=1, N-1
  do i=1, N-1
   Lw(i,j)=Aw(i,j)/(1.+alpha*Un(i-1,j))
   Ls(i,j)=As(i,j)/(1.+alpha*Ue(i,j-1))
   Lp(i,j)=Ap(i,j) + alpha*(Lw(i,j)*Un(i-1,j) + Ls(i,j)*Ue(i,j-1)) - Lw(i,j)*Ue(i-1,j) - Ls(i,j)*Un(i,j-1)
   Un(i,j)=(An(i,j) - alpha*Lw(i,j)*Un(i-1,j))/Lp(i,j)
   Ue(i,j)=(Ae(i,j) - alpha*Ls(i,j)*Ue(i,j-1))/Lp(i,j)
  enddo
 enddo

resN = 1.
inner_iter_counter=0.

do while (resN>eps_SIP)
inner_iter_counter=inner_iter_counter+1
 
 resN=0. 
 do j=1, N-1
  do i=1, N-1
   res(i,j)=Q(i,j)-Ae(i,j)*eta(i+1,j)-Aw(i,j)*eta(i-1,j)-An(i,j)*eta(i,j+1)-As(i,j)*eta(i,j-1)-Ap(i,j)*eta(i,j)

   resN=resN+ABS(res(i,j))
   enddo
 enddo
 

 do j=1, N-1
  do i=1, N-1
   res(i,j)=(res(i,j) - Lw(i,j)*res(i-1,j) - Ls(i,j)*res(i,j-1))/Lp(i,j)
  enddo
 enddo


 do j=N-1, 1, -1 
  do i=N-1, 1, -1 
   res(i,j)=res(i,j) - Un(i,j)*res(i,j+1) - Ue(i,j)*res(i+1,j)
   eta(i,j)=eta(i,j)+res(i,j)
  enddo
 enddo
 
 do i=0, N
  eta(i,0)=eta(i,1)
  eta(i,N)=eta(i,N-1)
  eta(0,i)=eta(1,i)
  eta(N,i)=eta(N-1,i)
 enddo

enddo
!------------------------Chemical potential computing(end)-----------------

!enddo

 time=time+dt

!---------------Finding maximum value of stream function(beg)--------------
max_psi=abs(psi(1,1))
max_eta=abs(eta(1,1))

do i=1, N-1
 do j=1, N-1
  if (abs(psi(i,j))>max_psi) max_psi=abs(psi(i,j))
  if (abs(eta(i,j))>max_eta) max_eta=abs(eta(i,j))
 enddo
enddo
!---------------Finding maximum value of stream function(end)--------------

if ((abs(max_psi_old-max_psi)<eps_stac) .and. (abs(max_eta_old-max_eta)<eps_stac)) stac_flg=.false.
max_psi_old=max_psi
max_eta_old=max_eta

!write(*,'(x,A,f10.6,x,A,f10.6,x,A,i5,x,A,f10.6,x,A,f10.6)') 'time=', time, 'psi_max=', max_psi , 'iter=', inner_iter_counter, 'eta_max', max_eta
!write(5,'(x,f10.6,x,f10.6,x,f10.6,x,f10.6)') time, max_psi, eta(15,15) !max_eta_new      

enddo
ampl_psi=max_psi

!-------------------------Write data to file(beg)--------------------------  
open(2, file = case_path//'/'//'out_psi.dat')
open(3, file = case_path//'/'//'out_T.dat')
open(7, file = case_path//'/'//'out_eta.dat')
open(8, file = case_path//'/'//'out_C.dat')
do i = 0, N
 do j = 0, N
  write(2,'(x,f10.5,x,f10.5,x,f10.5)') i*h,  j*h, psi(i,j) 
  write(3,'(x,f10.5,x,f10.5,x,f10.5)') i*h,  j*h, T(i,j) 
  write(7,'(x,f10.5,x,f10.5,x,f10.5)') i*h,  j*h, eta(i,j) 
  write(8,'(x,f10.5,x,f10.5,x,f10.5)') i*h,  j*h, eta(i,j)+soret*T(i,j) 
 enddo
enddo
!-------------------------Write data to file(end)--------------------------

end