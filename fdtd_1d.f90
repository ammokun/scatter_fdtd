program fdtd_1d

    use module_param
    implicit none

    real(8) :: t,dtime,nemax
    real(8) :: a_e,alp_e,bet_e,ome_e,gam_e,temg_in,temv_in
    real(8), dimension(0:imax) :: ez_ib,ez_i,ez_sb,ez_s,ez,ez_if
    real(8), dimension(0:imax) :: hy_s,hy_i,hy,vele_z,jdotE
    real(8), dimension(0:imax) :: e_rmi,e_rms,e_eff,xx,xx_in

    real(8), dimension(0:imax) :: nea,neb,ne,ne_m
    real(8), dimension(0:imax) :: nia,nib
    real(8), dimension(0:imax) :: nna,nnb
    real(8), dimension(0:imax) :: col,mu_e,mu_a,mu_i,mu_n,mu_in

    real(8), dimension(0:imax) ::x


    do i=0,imax
        x(i)=dx * i
        nea(i)=1.0d21 &
        *dexp(-(dble(i-imax*0.9d0)*dx)**2/(2.0d0*(1.0d-3*fb/f)**2)) 
        col(i)=0.0d0
    end do

    do i=0,imax

        ez_i(i)=0.0d0
        ez_s(i)=0.0d0
        ez(i)=0.0d0
        hy_i(i)=0.0d0
        hy_s(i)=0.0d0
        hy(i)=0.0d0

     end do
  !!!---------------------!!!
  !!!-----FDTD Scheme-----!!!
  !!!---------------------!!!
        do n=1,nmax
           !do m=1,divt !divt loop start
              
  !!---Electric Field---!!
  
  !-Incidnt-!
  
              ez_i(0)=amp*dsin(2.0d0*pi*f*t)
              
              do i=1,imax-1
                 ez_ib(i)=ez_i(i)
                 ez_i(i)=ez_i(i)+dt_e/eps0/dx*(hy_i(i)-hy_i(i-1))
              end do
                          
              ez_i(imax)=ez_ib(imax-1)+(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_i(imax-1)-ez_ib(imax))            
              ez_ib(imax)=ez_i(imax)
  
  !-[END] Incidnt-!
  
  !-Add Scatter-!
  
              do i=1,imax-1
                 
                 a_e=col(i)*dt_e/2.0d0
                 gam_e=1.0d0+a_e
                 alp_e=(1.0d0-a_e)/(1.0d0+a_e)               
                 ome_e=dsqrt(nea(i)*ee**2/me/eps0)
                 bet_e=ome_e**2*dt_e**2/4.0d0/gam_e               
                 
                 ez_sb(i)=ez_s(i)
                 
                 ez_s(i)=(1.0d0-bet_e)/(1.0d0+bet_e)*ez_s(i) &
                     +ee*nea(i)*dt_e/2.0d0/eps0 &
                     *(1.0d0+alp_e)/(1.0d0+bet_e)*vele_z(i) &
                     -bet_e/(1.0d0+bet_e)*(ez_i(i)+ez_ib(i)) &
                     +dt_e/(1.0d0+bet_e)/eps0/dx &
                     *(hy_s(i)-hy_s(i-1))
                 
                 vele_z(i)=alp_e*vele_z(i) &
                     -ee*dt_e/me/gam_e/2.0d0 &
                     *(ez_s(i)+ez_sb(i)+ez_i(i)+ez_ib(i))
                 
                 ez(i)=ez_i(i)+ez_s(i)                 
                 
                 jdotE=ee*vele_z(i)*ez(i) !joule heating by EM wave
                 
              end do
  
              ez_s(0)=ez_sb(1) &
                 +(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_s(1)-ez_sb(0))
              ez(0)=ez_s(0)+ez_i(0)
            
              !if(mirror==1)then
               !  ez_s(imax)=-ez_i(imax)
              !else 
                 ez_s(imax)=ez_sb(imax-1) &
                     +(cc*dt_e-dx)/(cc*dt_e+dx) &
                     *(ez_s(imax-1)-ez_sb(imax)) 
                 ez(imax)=ez_s(imax)+ez_i(imax) !!!!!          
              !end if
  
              ez_sb(0)=ez_s(0)
              ez_sb(imax)=ez_s(imax)
  
  !-[END] Add Scatter-!               
  
              t=t+0.5d0*dt_e
  
  !!---[END] Electric Field---!!
  
  !!---Magnetic Field---!!
  
              do i=0,imax-1
                 hy_i(i)=hy_i(i)+dt_e/mu0/dx*(ez_i(i+1)-ez_i(i))
              end do
              
              do i=0,imax-1
                 hy_s(i)=hy_s(i)+dt_e/mu0/dx*(ez_s(i+1)-ez_s(i))
                 hy(i)=hy_i(i)+hy_s(i)
              end do
  
  !!---[END] Magnetic Field---!!
  
              t=t+0.5d0*dt_e
              
              do i=0,imax
                 e_rmi(i)=e_rmi(i)+ez(i)**2
              end do
  
           !end do!divt loop end
           do i=0,imax
            e_rms(i)=dsqrt(e_rmi(i)/dble(divt))
            e_rmi(i)=0.0d0
         end do

         if(mod(n,sout)==0) then
            write(*,*) "n:",n
            write(filename,'(a,i5.5,a)') "./output/",n,".dat"
            open(10, file = filename, form = 'formatted')
            do i = 0, imax
               write(10,100) i,ez(i),hy(i),ez_i(i),hy_i(i),ez_s(i),hy_s(i)               
            end do
            close(10)
         end if
        end do !end main loop  
  !!!---------------------------!!!
  !!!-----[END] FDTD Scheme-----!!!
  !!!---------------------------!!!

        100     format(i15.4,1h,100(e15.7,1h,))

end program