program main

    use module_param
    implicit none

    real(8) :: t,dtime,nemax
    real(8) :: a_e,alp_e,bet_e,ome_e,gam_e,temg_in,temv_in
    real(8), dimension(0:jmax,0:imax) :: ez_ib,ez_i,ez_sb,ez_s,ez,ez_if
    real(8), dimension(0:jmax,0:imax) :: hy_s,hy_i,hy,hx_s,hx_i,hx,vele_z,jdotE
    real(8), dimension(0:jmax,0:imax) :: e_rmi,e_rms,e_eff,xx,xx_in
    real(8), dimension(0:jmax,0:imax) :: nea,neb,ne,ne_m
    real(8), dimension(0:jmax,0:imax) :: nia,nib
    real(8), dimension(0:jmax,0:imax) :: nna,nnb
    real(8), dimension(0:jmax,0:imax) :: col,mu_e,mu_a,mu_i,mu_n,mu_in

    real(8), dimension(0:imax) ::x
    real(8), dimension(0:jmax) ::y



    do i=0,imax
        x(i)=dx * i
    end do
    do j=0,jmax
        y(j)=dx * j
    end do

    t=0.0d0

    do i=0,imax
        do j=0,jmax

            nea(j,i)=0.0d0
            col(j,i)=0.0d0

            ez_i(j,i)=0.0d0
            ez_s(j,i)=0.0d0
            ez(j,i)=0.0d0
            hx_i(j,i)=0.0d0
            hx_s(j,i)=0.0d0
            hx(j,i)=0.0d0
            hy_i(j,i)=0.0d0
            hy_s(j,i)=0.0d0
            hy(j,i)=0.0d0

        end do
    end do
        do n=1,10 !main loop start
    !!!---------------------!!!
    !!!-----FDTD Scheme-----!!!
    !!!---------------------!!!
           do m=1,divt !divt loop start
              
  !!---Electric Field---!!
  
  !-Incidnt-!
            do j=1,jmax-1
              ez_i(j,0)=amp*dsin(2.0d0*pi*f*t)
              hy_i(j,0)=-ez_i(j,0)/377d0
              hx_i(j,0)=ez_i(j,0)/377d0
            end do

            do i=0,imax
                do j=0,jmax
                    ez_ib(j,i)=ez_i(j,i)
                end do
            end do

            do i=1,imax-1
              do j=1,jmax-1
                 ez_i(j,i)=ez_i(j,i)+dt_e/eps0/dx*(hy_i(j,i)-hy_i(j,i-1))-dt_e/eps0/dx*(hx_i(j,i)-hx_i(j-1,i))  !ez time integrate
              end do                          
            end do

            do j=1,jmax-1
                ez_i(j,imax)=ez_ib(j,imax-1)+(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_i(j,imax-1)-ez_ib(j,imax)) !mur absorption for incidnet field           
                ez_ib(j,imax)=ez_i(j,imax)
            end do
!
            do i=0,imax
                ez_i(0,i)=ez_ib(1,i)+(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_i(1,i)-ez_ib(0,i)) !mur absorption for incidnet field           
                ez_ib(0,i)=ez_i(0,i)
                
                ez_i(jmax,i)=ez_ib(jmax-1,i)+(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_i(jmax-1,i)-ez_ib(jmax,i)) !mur absorption for incidnet field           
                ez_ib(jmax,i)=ez_i(jmax,i)
            end do
  !-[END] Incidnt-!
  
  !-Add Scatter-!
  
              do i=1,imax-1
                do j=1,jmax-1
                 
                 a_e=col(j,i)*dt_e/2.0d0
                 gam_e=1.0d0+a_e
                 alp_e=(1.0d0-a_e)/(1.0d0+a_e)               
                 ome_e=dsqrt(nea(j,i)*ee**2/me/eps0)
                 bet_e=ome_e**2*dt_e**2/4.0d0/gam_e               
                 
                 ez_sb(j,i)=ez_s(j,i)
                 
                 ez_s(j,i)=(1.0d0-bet_e)/(1.0d0+bet_e)*ez_s(j,i) &
                     +ee*nea(j,i)*dt_e/2.0d0/eps0 &
                     *(1.0d0+alp_e)/(1.0d0+bet_e)*vele_z(j,i) &
                     -bet_e/(1.0d0+bet_e)*(ez_i(j,i)+ez_ib(j,i)) &
                     +dt_e/(1.0d0+bet_e)/eps0/dx &
                     *(hy_s(j,i)-hy_s(j,i-1)) &
                     -dt_e/(1.0d0+bet_e)/eps0/dx &
                     *(hy_s(j,i)-hy_s(j-1,i))
                 
                 vele_z(j,i)=alp_e*vele_z(j,i) &
                     -ee*dt_e/me/gam_e/2.0d0 &
                     *(ez_s(j,i)+ez_sb(j,i)+ez_i(j,i)+ez_ib(j,i))
                 
                 ez(j,i)=ez_i(j,i)+ez_s(j,i)                 
                 
                 jdotE(j,i)=ee*vele_z(j,i)*ez(j,i) !joule heating by EM wave
                end do
              end do
  
              !Mur absorption boundary for scattered field
              do j=1,jmax-1
                ez_s(j,0)=ez_sb(j,1) +(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_s(j,1)-ez_sb(j,0))
                ez(j,0)=ez_s(j,0)+ez_i(j,0)
              end do

              do i=0,imax
                ez_s(0,i)=ez_sb(1,i) +(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_s(1,i)-ez_sb(0,i))
                ez(0,i)=ez_s(0,i)+ez_i(0,i)

                ez_s(jmax,i)=ez_sb(jmax-1,i) +(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_s(jmax-1,i)-ez_sb(jmax,i))
                ez(jmax,i)=ez_s(jmax,i)+ez_i(jmax,i)

              end do
            
            if(mirror==1)then
                do j=1,jmax-1
                     ez_s(j,imax)=-ez_i(j,imax)
                end do
            else 
                do j=1,jmax-1
                    ez_s(j,imax)=ez_sb(j,imax-1)+(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_s(j,imax-1)-ez_sb(j,imax)) 
                    ez(j,imax)=ez_s(j,imax)+ez_i(j,imax) !!!!!
                end do
                                
            end if
            
            do j=0,jmax
              ez_sb(j,0)=ez_s(j,0)
              ez_sb(j,imax)=ez_s(j,imax)
            end do
  !-[END] Add Scatter-!               
  
              t=t+0.5d0*dt_e
  
  !!---[END] Electric Field---!!
  
  !!---Magnetic Field---!!
  
            do i=0,imax-1
                do j=0,jmax
                    hy_i(j,i)=hy_i(j,i)+dt_e/mu0/dx*(ez_i(j,i+1)-ez_i(j,i))
                end do
            end do

            do i=0,imax
                do j=0,jmax-1
                    hx_i(j,i)=hx_i(j,i)-dt_e/mu0/dx*(ez_i(j+1,i)-ez_i(j,i))
                end do
            end do
              
              do i=0,imax-1
                do j=0,jmax
                 hy_s(j,i)=hy_s(j,i)+dt_e/mu0/dx*(ez_s(j,i+1)-ez_s(j,i))
                 hy(j,i)=hy_i(j,i)+hy_s(j,i)
                end do
              end do
            do i=0,imax
                do j=0,jmax-1
                    hx_s(j,i)=hx_s(j,i)+dt_e/mu0/dx*(ez_s(j+1,i)-ez_s(j,i))
                    hx(j,i)=hx_i(j,i)+hx_s(j,i)
                end do
            end do
  
  !!---[END] Magnetic Field---!!
  
              t=t+0.5d0*dt_e
              
              do i=0,imax
                do j=0,jmax
                    e_rmi(j,i)=e_rmi(j,i)+ez(j,i)**2
                end do
              end do
  
           end do!divt loop end
            do i=0,imax
                do j=0,jmax
                e_rms(j,i)=dsqrt(e_rmi(j,i)/dble(divt))
                e_rmi(j,i)=0.0d0
                end do
            end do
        end do !end main loop  
  !!!---------------------------!!!
  !!!-----[END] FDTD Scheme-----!!!
  !!!---------------------------!!!


           open(10, file = 'output.dat', form = 'formatted')
           do i = 0, imax
                do j=0,jmax
                    write(10,*) i,j,ez_i(j,i),hx_i(j,i),hy_i(j,i)
                end do
                write(10,*) " "
           end do
           close(10)

end program