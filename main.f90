program main

    use module_param
    implicit none

    real(8) :: t,dtime,nemax
    real(8) :: a_e,alp_e,bet_e,ome_e,gam_e,temg_in,temv_in
    real(8), dimension(0:jmax,0:imax) :: ez_ib,ez_i,ez_sb,ez_s,ez,ez_if
    real(8), dimension(0:jmax,0:imax) :: hy_s,hy_i,hy,hx_s,hx_i,hx,vele_z,jdotE,hy_sb,hy_ib,hx_sb,hx_ib
    real(8), dimension(0:jmax,0:imax) :: e_rmi,e_rms,e_eff,xx,xx_in
    real(8), dimension(0:jmax,0:imax) :: nea,neb,ne,ne_m
    real(8), dimension(0:jmax,0:imax) :: ne_fdtd,col_fdtd
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
            ne_fdtd(j,i)=0.0d0
            col_fdtd(j,i)=0.0d0

            ez_i(j,i)=0.0d0
            ez_s(j,i)=0.0d0
            ez(j,i)=0.0d0
            hx_i(j,i)=0.0d0
            hx_s(j,i)=0.0d0
            hx(j,i)=0.0d0
            hy_i(j,i)=0.0d0
            hy_s(j,i)=0.0d0
            hy(j,i)=0.0d0
            e_rms(j,i)=0.0d0
            e_rmi(j,i)=0.0d0
            vele_z(j,i)=0.0d0
        end do
    end do

   do i=200,220
        do j=140,160
            nea(j,i)=1.0d021 &
            *dexp(-(dble(i-200)*dx)**2/(2.0d0*(1.0d-3*fb/f)**2)) &
            *dexp(-(dble(j-149)*dx)**2/(2.0d0*(1.0d-3*fb/f)**2)) !!!!!
            col(j,i)=1.0d12 &
            *dexp(-(dble(i-200)*dx)**2/(2.0d0*(1.0d-3*fb/f)**2)) &
            *dexp(-(dble(j-149)*dx)**2/(2.0d0*(1.0d-3*fb/f)**2)) !!!!!
        end do
    end do          

    do i=1,imax-1
      do j=1,jmax-1
          ne_fdtd(j,i)=nea(j,i)
          col_fdtd(j,i)=col(j,i)
      end do
  end do                
        do n=1,nmax !main loop start
    !!!---------------------!!!
    !!!-----FDTD Scheme-----!!!
    !!!---------------------!!!
           !do m=1,divt !loop start
              
  !!---Electric Field---!!
  
  !-Incidnt-!
            do j=0,jmax
              ez_i(j,0)=amp*dsin(2.0d0*pi*f*t)
              hy_i(j,0)=-ez_i(j,0)/377d0
              hx_i(j,0)=ez_i(j,0)/377d0
            end do

!            do j=1,jmax-1
                !ez_i(150,0)=amp*dsin(2.0d0*pi*f*t)
!                hy_i(150,0)=-ez_i(j,0)/377d0
!                hx_i(150,0)=ez_i(j,0)/377d0
!            end do

            do i=0,imax
                do j=0,jmax
                    ez_ib(j,i)=ez_i(j,i)
                end do
            end do

            do i=1,imax
              do j=1,jmax
                 ez_i(j,i)=ez_i(j,i)+dt_e/eps0/dx*(hy_i(j,i)-hy_i(j,i-1))-dt_e/eps0/dx*(hx_i(j,i)-hx_i(j-1,i))  !ez time integrate
              end do                          
            end do 

            do j=0,jmax
                ez_i(j,imax)=ez_ib(j,imax-1)+(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_i(j,imax-1)-ez_ib(j,imax)) !mur absorption for incidnet field           
                ez_ib(j,imax)=ez_i(j,imax)
                ez_ib(j,imax-1)=ez_i(j,imax-1)
                
                ez_i(j,0)=ez_ib(j,1)+(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_i(j,1)-ez_ib(j,0)) !mur absorption for incidnet field           
                ez_ib(j,0)=ez_i(j,0)
                ez_ib(j,1)=ez_i(j,1)
            end do

            !do i=0,imax
            !    ez_i(0,i)=ez_ib(1,i)+(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_i(1,i)-ez_ib(0,i)) !mur absorption for incidnet field           
            !    ez_ib(0,i)=ez_i(0,i)
            !    ez_ib(1,i)=ez_i(1,i)
            !    
            !    ez_i(jmax,i)=ez_ib(jmax-1,i)+(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_i(jmax-1,i)-ez_ib(jmax,i)) !mur absorption for incidnet field           
            !    ez_ib(jmax,i)=ez_i(jmax,i)
            !    ez_ib(jmax-1,i)=ez_i(jmax-1,i)
            !    !jmaxあたりで減衰する
            !end do
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
                 
                 ez_s(j,i)=(1.0d0-bet_e)/(1.0d0+bet_e)*ez_sb(j,i) &
                     +ee*nea(j,i)*dt_e/2.0d0/eps0 &
                     *(1.0d0+alp_e)/(1.0d0+bet_e)*vele_z(j,i) &
                     -bet_e/(1.0d0+bet_e)*(ez_i(j,i)+ez_ib(j,i)) &
                     +dt_e/(1.0d0+bet_e)/eps0/dx &
                     *(hy_s(j,i)-hy_s(j,i-1)) &
                     -dt_e/(1.0d0+bet_e)/eps0/dx &
                     *(hx_s(j,i)-hx_s(j-1,i))
                ! ez_s(j,i)=(1.0d0-bet_e)/(1.0d0+bet_e)*ez_s(j,i) &
                !     +dt_e/(1.0d0+bet_e)/eps0/dx &
                !     *(hy_s(j,i)-hy_s(j,i-1)) &
                !     -dt_e/(1.0d0+bet_e)/eps0/dx &
                !     *(hx_s(j,i)-hx_s(j-1,i))
 
                 vele_z(j,i)=alp_e*vele_z(j,i) &
                     -ee*dt_e/me/gam_e/2.0d0 &
                     *(ez_s(j,i)+ez_sb(j,i)+ez_i(j,i)+ez_ib(j,i))
                end do
              end do
              
              !do i=200,250
              !  do j=125,175
              !    ez_s(j,i)=-ez_i(j,i)
              !  end do
              !end do

              do i=1,imax-1
                do j=1,jmax-1
                 ez(j,i)=ez_i(j,i)+ez_s(j,i)                 
                 
                 jdotE(j,i)=ee*vele_z(j,i)*ez(j,i) !joule heating by EM wave
                end do
              end do


  
              !Mur absorption boundary for scattered field
              do j=1,jmax-1
                ez_s(j,0)=ez_sb(j,1) +(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_s(j,1)-ez_sb(j,0))
                ez_sb(j,0)=ez_s(j,0)
                ez_sb(j,1)=ez_s(j,1)
                
                ez(j,0)=ez_s(j,0)+ez_i(j,0)
                ez_s(j,imax)=ez_sb(j,imax-1)+(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_s(j,imax-1)-ez_sb(j,imax)) 
                ez_sb(j,imax)=ez_s(j,imax)
                ez_sb(j,imax-1)=ez_s(j,imax-1)
                
                ez(j,imax)=ez_s(j,imax)+ez_i(j,imax) !!!!!
              end do

              do i=1,imax-1
                ez_s(0,i)=ez_sb(1,i) +(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_s(1,i)-ez_sb(0,i))
                ez_sb(0,i)=ez_s(0,i)
                ez_sb(1,i)=ez_s(1,i)

                ez(0,i)=ez_s(0,i)+ez_i(0,i)
                
                ez_s(jmax,i)=ez_sb(jmax-1,i) +(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_s(jmax-1,i)-ez_sb(jmax,i))
                ez_sb(jmax,i)=ez_s(jmax,i)
                ez_sb(jmax-1,i)=ez_s(jmax-1,i)


                ez(jmax,i)=ez_s(jmax,i)+ez_i(jmax,i)
              end do
!            
!            !if(mirror==1)then
!            !    do j=0,jmax
!            !         ez_s(j,imax)=-ez_i(j,imax)
!            !    end do
!            !else 
!                do j=0,jmax
!                    ez_s(j,imax)=ez_sb(j,imax-1)+(cc*dt_e-dx)/(cc*dt_e+dx)*(ez_s(j,imax-1)-ez_sb(j,imax)) 
!                    ez(j,imax)=ez_s(j,imax)+ez_i(j,imax) !!!!!
!                end do
!                                
!            !end if
!            
!            do j=0,jmax
!              ez_sb(j,0)=ez_s(j,0)
!              ez_sb(j,imax)=ez_s(j,imax)
!            end do
  !-[END] Add Scatter-!               
  
              t=t+0.5d0*dt_e
  
  !!---[END] Electric Field---!!
  
  !!---Magnetic Field---!!
  
              do i=0,imax-1
                do j=0,jmax
                    hy_i(j,i)=hy_i(j,i)+dt_e/mu0/dx*(ez_i(j,i+1)-ez_i(j,i))
                    hy_s(j,i)=hy_s(j,i)+dt_e/mu0/dx*(ez_s(j,i+1)-ez_s(j,i))
                    hy(j,i)=hy_i(j,i)+hy_s(j,i)
                end do
              end do
            do i=1,imax-1
                hy_i(0,i)=hy_ib(1,i)+(cc*dt_e-dx)/(cc*dt_e+dx)*(hy_i(1,i)-hy_ib(0,i)) !mur absorption for incidnet field           
                hy_ib(0,i)=hy_i(0,i)
                hy_ib(1,i)=hy_i(1,i)
                
                hy_i(jmax,i)=hy_ib(jmax-1,i)+(cc*dt_e-dx)/(cc*dt_e+dx)*(hy_i(jmax-1,i)-hy_ib(jmax,i)) !mur absorption for incidnet field           
                hy_ib(jmax,i)=hy_i(jmax,i)
                hy_ib(jmax-1,i)=hy_i(jmax-1,i)
                
                hy_s(0,i)=hy_sb(1,i)+(cc*dt_e-dx)/(cc*dt_e+dx)*(hy_s(1,i)-hy_sb(0,i)) !mur absorption for incidnet field           
              hy_sb(0,i)=hy_s(0,i)
              hy_sb(1,i)=hy_s(1,i)
                
                hy_s(jmax,i)=hy_sb(jmax-1,i)+(cc*dt_e-dx)/(cc*dt_e+dx)*(hy_s(jmax-1,i)-hy_sb(jmax,i)) !mur absorption for incidnet field           
                hy_sb(jmax,i)=hy_s(jmax,i)
                hy_sb(jmax-1,i)=hy_s(jmax-1,i)

                hy(0,i)=hy_i(0,i)+hy_s(0,i)
                hy(jmax,i)=hy_i(jmax,i)+hy_s(jmax,i)
            end do

            do j=1,jmax-1
              hy_i(j,0)=hy_ib(j,1)+(cc*dt_e-dx)/(cc*dt_e+dx)*(hy_i(j,1)-hy_ib(j,0)) !mur absorption for incidnet field           
              hy_ib(j,0)=hy_i(j,0)
              hy_ib(j,1)=hy_i(j,1)
              
              hy_i(j,imax)=hy_ib(j,imax-1)+(cc*dt_e-dx)/(cc*dt_e+dx)*(hy_i(j,imax-1)-hy_ib(j,imax)) !mur absorption for incidnet field           
              hy_ib(j,imax)=hy_i(j,imax)
              hy_ib(j,imax-1)=hy_i(j,imax-1)

              hy_s(j,0)=hy_sb(j,1)+(cc*dt_e-dx)/(cc*dt_e+dx)*(hy_s(j,1)-hy_sb(j,0)) !mur absorption for incidnet field           
              hy_sb(j,0)=hy_s(j,0)
              hy_sb(j,1)=hy_s(j,1)
              
              hy_s(j,imax)=hy_sb(j,imax-1)+(cc*dt_e-dx)/(cc*dt_e+dx)*(hy_s(j,imax-1)-hy_sb(jmax,i)) !mur absorption for incidnet field           
              hy_sb(j,imax)=hy_s(j,imax)
              hy_sb(j,imax-1)=hy_s(j,imax-1)

              hy(j,0)=hy_i(j,0)+hy_s(j,0)
              hy(j,imax)=hy_i(j,imax)+hy_s(j,imax)
          end do

            do i=0,imax
                do j=0,jmax-1
                    hx_i(j,i)=hx_i(j,i)-dt_e/mu0/dx*(ez_i(j+1,i)-ez_i(j,i))
                    hx_s(j,i)=hx_s(j,i)-dt_e/mu0/dx*(ez_s(j+1,i)-ez_s(j,i))
                    hx(j,i)=hx_i(j,i)+hx_s(j,i)
                end do
            end do
            do i=0,imax
              hx_i(0,i)=hx_ib(1,i)+(cc*dt_e-dx)/(cc*dt_e+dx)*(hx_i(1,i)-hx_ib(0,i)) !mur absorption for incidnet field           
              hx_ib(0,i)=hx_i(0,i)
              hx_ib(1,i)=hx_i(1,i)
              
              hx_i(jmax,i)=hx_ib(jmax-1,i)+(cc*dt_e-dx)/(cc*dt_e+dx)*(hx_i(jmax-1,i)-hx_ib(jmax,i)) !mur absorption for incidnet field           
              hx_ib(jmax,i)=hx_i(jmax,i)
              hx_ib(jmax-1,i)=hx_i(jmax-1,i)

              hx_s(0,i)=hx_sb(1,i)+(cc*dt_e-dx)/(cc*dt_e+dx)*(hx_s(1,i)-hx_sb(0,i)) !mur absorption for incidnet field           
              hx_sb(0,i)=hx_s(0,i)
              hx_sb(1,i)=hx_s(1,i)
              
              hx_s(jmax,i)=hx_sb(jmax-1,i)+(cc*dt_e-dx)/(cc*dt_e+dx)*(hx_s(jmax-1,i)-hx_sb(jmax,i)) !mur absorption for incidnet field           
              hx_sb(jmax,i)=hx_s(jmax,i)
              hx_sb(jmax-1,i)=hx_s(jmax-1,i)

              hx(0,i)=hx_i(0,i)+hx_s(0,i)
              hx(jmax,i)=hx_i(jmax,i)+hx_s(jmax,i)
            end do

            do j=0,jmax
              hx_i(j,0)=hx_ib(j,1)+(cc*dt_e-dx)/(cc*dt_e+dx)*(hx_i(j,1)-hx_ib(j,0)) !mur absorption for incidnet field           
              hx_ib(j,0)=hx_i(j,0)
              hx_ib(j,1)=hx_i(j,1)
              
              hx_i(j,imax)=hx_ib(j,imax-1)+(cc*dt_e-dx)/(cc*dt_e+dx)*(hx_i(j,imax-1)-hx_ib(j,imax)) !mur absorption for incidnet field           
              hx_ib(j,imax)=hx_i(j,imax)
              hx_ib(j,imax-1)=hx_i(j,imax-1)

              hx_s(j,0)=hx_sb(j,1)+(cc*dt_e-dx)/(cc*dt_e+dx)*(hx_s(j,1)-hx_sb(j,0)) !mur absorption for incidnet field           
              hx_sb(j,0)=hx_s(j,0)
              hx_sb(j,1)=hx_s(j,1)
              
              hx_s(j,imax)=hx_sb(j,imax-1)+(cc*dt_e-dx)/(cc*dt_e+dx)*(hx_s(j,imax-1)-hx_sb(j,imax)) !mur absorption for incidnet field           
              hx_sb(j,imax)=hx_s(j,imax)
              hx_sb(j,imax-1)=hx_s(j,imax-1)

              hx(j,0)=hx_i(j,0)+hx_s(j,0)
              hx(j,imax)=hx_i(j,imax)+hx_s(j,imax)              
            end do


            
  
  !!---[END] Magnetic Field---!!
  
              t=t+0.5d0*dt_e
              
              do i=0,imax
                do j=0,jmax
                    e_rmi(j,i)=e_rmi(j,i)+ez(j,i)**2
                end do
              end do
              if(mod(n,sout)==0) then
                write(*,*) "n:",n
                write(filename,'(a,i5.5,a)') "./output/",n,".dat"
                open(10, file = filename, form = 'formatted')
                do i = 0, imax
                     do j=0,jmax
                         write(10,100) i,j,ez(j,i),hx(j,i),hy(j,i),ez_i(j,i),hx_i(j,i),hy_i(j,i),ez_s(j,i),hx_s(j,i),hy_s(j,i),&
                         vele_z(j,i)!e_rms(j,i)
                     end do
                     write(10,*) " "
                end do

                close(10)
                !do i=0,imax
                !  do j=0,jmax
                !    if(ez(j,i)>10.0d0) then
                !      write(*,*)"erorr!, (i,j)",i,j
                !      stop
                !    end if
                !  end do
                !end do
            
            end if
           !end do!divt loop end
            !do i=0,imax
            !    do j=0,jmax
            !    e_rms(j,i)=dsqrt(e_rmi(j,i)/dble(divt))
            !    e_rmi(j,i)=0.0d0
            !    end do
            !end do

        
        end do !end main loop  

        100     format(i15.4,i15.4,1h,100(e15.7,1h,))
  !!!---------------------------!!!
  !!!-----[END] FDTD Scheme-----!!!
  !!!---------------------------!!!



end program