! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing interface to greg thompson's microphysics code 
!
module mphys_thompson09n

  Use parameters, only : num_h_moments, num_h_bins, nspecies, nz, dt &
       , h_names, mom_units, max_char_len, nx
  Use column_variables
  Use physconst, only : p0, r_on_cp, pi

  Use module_mp_thompson09n
  Use diagnostics, only: save_dg, i_dgtime

  Implicit None

  !Logical switches 
  logical :: micro_unset=.True.
  integer:: ih, imom
  character(max_char_len) :: name, units

contains
  
  Subroutine mphys_thompson09_interfacen
    
    real :: pptrain, pptsnow, pptgraul, pptice
    real :: pptrain_2d(nx), pptsnow_2d(nx), pptgraul_2d(nx) &
         , pptice_2d(nx), pptrain_level(nz), pptrain_2d_prof(nz,nx)   
    real :: t1d(nz), p1d(nz), dz1d(nz),qv1d(nz),qc1d(nz) &
         ,qr1d(nz), qi1d(nz), ni1d(nz), qs1d(nz), qg1d(nz) &
         ,nr1d(nz), dzq(nz)
    real :: nc1d(nz), nifa1d(nz), nwfa1d(nz), w1d(nz)
    real :: rainprod1d(nz), evapprod1d(nz)
    real :: rainprod(nz), evapprod(nz)
    integer :: kts, kte, i, j, k

    
    kts=1
    kte=nz
    j=1

    qc1d(:) = 0.0
    qr1d(:) = 0.0
    nr1d(:) = 0.0
    qi1d(:) = 0.0
    ni1d(:) = 0.0
    qs1d(:) = 0.0
    qg1d(:) = 0.0

    do i=1,nx
       pptrain = 0. !taxa de precipitacao de chuva
       pptsnow = 0. !taxa e precipitacao de neve
       pptgraul = 0. !taxa de precipitacao de graupel
       pptice = 0. !taxa de precipitacao de gelo
       do k=1,nz
          t1d(k) = (theta(k,i) + (dtheta_adv(k,i)+dtheta_div(k,i))*dt )*exner(k,i) !Temp potencial em um determinado nivel vertical e ponto horizontal
          p1d(k) = p0*exner(k,i)**(1./r_on_cp) ! pressao em um determinado nivel vertical e ponto horizontal
!          dz1d(k) = dz(k) ! espessura da camada em um determinado nivel vertical e ponto horizontal 
          dzq(k) = dz(k) ! espessura da camada em um determinado nivel vertical e ponto horizontal 
          qv1d(k) = qv(k,i) + (dqv_adv(k,i)+dqv_div(k,i))*dt !umidade especifica em um determinado nivel vertical e ponto horizontal 
          
          qc1d(k) = hydrometeors(k,i,1)%moments(1,1) &       !Conteudo de agua condensada (nuvem) em um
               + (dhydrometeors_adv(k,i,1)%moments(1,1) &    !determinado nivel da vertical e ponto
               + dhydrometeors_div(k,i,1)%moments(1,1))*dt   !horizontal
          
          qr1d(k) = hydrometeors(k,i,2)%moments(1,1) &       !conteudo de agua na forma de chuva em um 
               + (dhydrometeors_adv(k,i,2)%moments(1,1) &    !determinado nivel da vertical e ponto
               + dhydrometeors_div(k,i,2)%moments(1,1))*dt   !horizontal
          
          nr1d(k) = hydrometeors(k,i,2)%moments(1,2) &       !conteudo de numero de particulas de chuva
               + (dhydrometeors_adv(k,i,2)%moments(1,2) &    !em um determinado nivel da vertical e
               + dhydrometeors_div(k,i,2)%moments(1,2))*dt   !ponto na horizontal
          
          if (.not. iiwarm) then   !Se nao for um ambiente quente
             qi1d(k) = hydrometeors(k,i,3)%moments(1,1) &     !conteudo de agua na forma de gelo em um
                  + (dhydrometeors_adv(k,i,3)%moments(1,1) &    !determinado nivel da vertical e ponto
                  + dhydrometeors_div(k,i,3)%moments(1,1))*dt    !horizontal
             
             ni1d(k) = hydrometeors(k,i,3)%moments(1,2) &      !conteudo de numero de particulas de gelo 
                  + (dhydrometeors_adv(k,i,3)%moments(1,2) &    !em um determinado nuvel da vertical e
                  + dhydrometeors_div(k,i,3)%moments(1,2))*dt   !ponto na horizontal 
         
             qs1d(k) = hydrometeors(k,i,4)%moments(1,1) &      !conteudo de agua na forma de neve em um
                  + (dhydrometeors_adv(k,i,4)%moments(1,1) &    !determinado nivel da vertical e ponto
                  + dhydrometeors_div(k,i,4)%moments(1,1))*dt    !horizontal 
             
             qg1d(k) = hydrometeors(k,i,5)%moments(1,1) &      !conteudo de agua na forma de graupel em um
                  + (dhydrometeors_adv(k,i,5)%moments(1,1) &    !determinado nivel da vertical e ponto
                  + dhydrometeors_div(k,i,5)%moments(1,1))*dt    !horizontal 
             
          end if
          
       end do
       
       ! Initialise microphysics 
       if (micro_unset)then
          call thompson_init
          micro_unset=.False.
       end if

! qv1d: vapor d'água específico.
! qc1d: água de nuvem específica.
! qi1d: gelo específico.
! qr1d: chuva específica.
! qs1d: neve específica.
! qg1d: graupel específico.
! ni1d: número de partículas de gelo.
! nr1d: número de gotas de chuva.
! nc1d: número de partículas de nuvem.
! nwfa1d: número de partículas de água líquida em nuvens.
! nifa1d: número de partículas de gelo em nuvens.
! t1d: temperatura.
! p1d: pressão.
! w1d: velocidade vertical.
! dz1d: espessura vertical da camada.
! pptrain: Acumulação de precipitação convectiva.
! pptsnow: Acumulação de precipitação de neve.
! pptgraul: Acumulação de precipitação de graupel.
! pptice: Acumulação de precipitação de gelo.
! rainprod1d: produção de chuva.
! evapprod1d: produção de evaporação.
! kts: Índice vertical de início do loop.
! kte: Índice vertical de término do loop.
! dt: Passo de tempo.
! i: Índice horizontal da grade.
! j: Índice vertical da grade.


!       call mp_thompson(qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, ni1d,     &
!       nr1d, nc1d, nwfa1d, nifa1d, t1d, p1d, w1d, dz1d,  &
!       pptrain, pptsnow, pptgraul, pptice, &
 ! KiD specific diagnostic  
 !      pptrain_level,                                        &
 ! End KiD diag                      
 !      rainprod1d, evapprod1d, &
 !      kts, kte, dt, i, j)


       call mp_thompson(qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, ni1d, &
       nr1d, nc1d, nwfa1d, nifa1d, t1d, p1d, w1d, dzq, &
       pptrain, pptsnow, pptgraul, pptice, &
 ! KiD specific diagnostic  
 !      pptrain_level,                                        &
 ! End KiD diag   
#ifdef WRF_CHEM                  
       rainprod, evapprod, &
#endif
       kts, kte, dt, i, j)


       if (nx == 1) then
          ! Save some diagnostics
          imom=1
          !rain ppt
          ih=2
          name='surface_ppt_for_'//trim(h_names(ih))
          units=trim(mom_units(imom))//' m'
          call save_dg(pptrain, name, i_dgtime,  units, dim='time')
       !ice ppt
          ih=3
          name='surface_ppt_for_'//trim(h_names(ih))
          units=trim(mom_units(imom))//' m'
          call save_dg(pptice, name, i_dgtime,  units, dim='time')
          !snow ppt
          ih=4
          name='surface_ppt_for_'//trim(h_names(ih))
          units=trim(mom_units(imom))//' m'
          call save_dg(pptsnow, name, i_dgtime,  units, dim='time')
          !graupel ppt
          ih=5
          name='surface_ppt_for_'//trim(h_names(ih))
          units=trim(mom_units(imom))//' m'
          call save_dg(pptgraul, name, i_dgtime,  units, dim='time')
          !total ppt
          name='total_surface_ppt'
          units=trim(mom_units(imom))//' m'
          call save_dg((pptice+pptrain+pptsnow+pptgraul)/nx, name, i_dgtime, &
               units, dim='time')
       else
          pptrain_2d(i) = pptrain
          pptice_2d(i) = pptice
          pptsnow_2d(i) = pptsnow
          pptgraul_2d(i) = pptgraul
          !
          ! pptrain_level is declared and calced in 
          ! the module module_mp_thompson09
          !pptrain_2d_prof(:,i) = pptrain_level(:)
       endif
       

    ! back out tendencies
       !Essa parte do código está calculando as tendências das variáveis 
       !de microfísica de nuvens e precipitação ao longo da dimensão vertical 
       do k=1,nz
          dtheta_mphys(k,i)=(t1d(k)/exner(k,i)-theta(k,i))/dt &  !Tendência da temperatura potencial para um nível vertical específico e ponto horizontal. 
               - ( dtheta_adv(k,i)+dtheta_div(k,i))    !é calculada com base na diferença entre a temperatura potencial atual (theta(k,i)) 
                                                       !e a temperatura potencial calculada a partir das variáveis t1d, exner, dtheta_adv, dtheta_div, dt
                                                       !O resultado é ajustado pela contribuição das tendências advectivas e divergentes.
          
          dqv_mphys(k,i)=(qv1d(k) - qv(k,i))/dt & !Tendência da umidade específica para um nível vertical específico e ponto horizontal. 
               - ( dqv_adv(k,i)+dqv_div(k,i))     !é calculada com base na diferença entre a umidade específica atual (qv(k,i))
                                                  !e a umidade específica calculada a partir das variáveis qv1d, dt, dqv_adv, dqv_div
                                                  !O resultado é ajustado pela contribuição das tendências advectivas e divergentes.
          
          dhydrometeors_mphys(k,i,1)%moments(1,1)= &   !Tendência do momento da fase de condensado (nuvem) para um nível vertical específico e ponto horizontal. 
               (qc1d(k)-hydrometeors(k,i,1)%moments(1,1))/dt & !A tendência é calculada com base na diferença entre o conteúdo de água condensada atual (qc1d(k)) 
               - (dhydrometeors_adv(k,i,1)%moments(1,1)  & ! e o conteúdo de água condensada da fase de condensado calculado a partir das variáveis hydrometeors,
               + dhydrometeors_div(k,i,1)%moments(1,1)) ! dt, dhydrometeors_adv, dhydrometeors_div. O resultado é ajustado pela contribuição das tendências advectivas e divergentes.

          dhydrometeors_mphys(k,i,2)%moments(1,1)= & !Tendência do momento da fase de chuva para um nível vertical específico e ponto horizontal.
               (qr1d(k)-hydrometeors(k,i,2)%moments(1,1))/dt & !Similar à dhydrometeors_mphys(k,i,1)%moments(1,1), mas para a fase de chuva.
               - (dhydrometeors_adv(k,i,2)%moments(1,1)  &
               + dhydrometeors_div(k,i,2)%moments(1,1))

          dhydrometeors_mphys(k,i,2)%moments(1,2)= & !Tendência do momento do número de partículas de chuva para um nível vertical específico e ponto horizontal.
               (nr1d(k)-hydrometeors(k,i,2)%moments(1,2))/dt & !Similar à dhydrometeors_mphys(k,i,2)%moments(1,1), mas para o número de partículas de chuva.
               - (dhydrometeors_adv(k,i,2)%moments(1,2)  &
               + dhydrometeors_div(k,i,2)%moments(1,2))

          if (.not.iiwarm)then 
            dhydrometeors_mphys(k,i,3)%moments(1,1)= &
               (qi1d(k)-hydrometeors(k,i,3)%moments(1,1))/dt & !gelo = indice 3
               - (dhydrometeors_adv(k,i,3)%moments(1,1)  &
               + dhydrometeors_div(k,i,3)%moments(1,1))

            dhydrometeors_mphys(k,i,3)%moments(1,2)= &
               (ni1d(k)-hydrometeors(k,i,3)%moments(1,2))/dt & 
               - (dhydrometeors_adv(k,i,3)%moments(1,2)  &
               + dhydrometeors_div(k,i,3)%moments(1,2))

            dhydrometeors_mphys(k,i,4)%moments(1,1)= &
               (qs1d(k)-hydrometeors(k,i,4)%moments(1,1))/dt & !neve = indice 4
               - (dhydrometeors_adv(k,i,4)%moments(1,1)  &
               + dhydrometeors_div(k,i,4)%moments(1,1))

            dhydrometeors_mphys(k,i,5)%moments(1,1)= &
               (qg1d(k)-hydrometeors(k,i,5)%moments(1,1))/dt & !graupel = indice 5
               - (dhydrometeors_adv(k,i,5)%moments(1,1)  &
               + dhydrometeors_div(k,i,5)%moments(1,1))
          end if
       end do
    end do
    
    if (nx > 1) then
       ! Save some mean scalar diagnostics
       imom=1
       !rain ppt
       ih=2
       name='surface_ppt_for_'//trim(h_names(ih))
       units=trim(mom_units(imom))//' m'
       call save_dg(pptrain_2d/nx, name, i_dgtime,  units, dim='time') !Média da taxa de precipitação de chuva ao longo de nx
       !ice ppt
       ih=3
       name='surface_ppt_for_'//trim(h_names(ih))
       units=trim(mom_units(imom))//' m'
       call save_dg(pptice_2d/nx, name, i_dgtime,  units, dim='time') !Média da taxa de precipitação de gelo ao longo de nx
       !snow ppt
       ih=4
       name='surface_ppt_for_'//trim(h_names(ih))
       units=trim(mom_units(imom))//' m'
       call save_dg(pptsnow_2d/nx, name, i_dgtime,  units, dim='time') !Média da taxa de precipitação de neve ao longo de nx
       !graupel ppt
       ih=5
       name='surface_ppt_for_'//trim(h_names(ih))
       units=trim(mom_units(imom))//' m'
       call save_dg(pptgraul_2d/nx, name, i_dgtime,  units, dim='time') !Média da taxa de precipitação de graupel ao longo de nx
       !total ppt
       name='total_surface_ppt'
       units=trim(mom_units(imom))//' m'
       call save_dg((pptice_2d+pptrain_2d+pptsnow_2d+pptgraul_2d)/nx, & !media total de precipitação ao longo de nx
            name, i_dgtime, units, dim='time') 
       
       ! save all horizontal columns
       imom=1
       !rain ppt
       ih=2
       name='surface_ppt_for_'//trim(h_names(ih))
       units=trim(mom_units(imom))//' m'
       call save_dg(pptrain_2d, name, i_dgtime,  units, dim='time')
       !ice ppt
       ih=3
       name='surface_ppt_for_'//trim(h_names(ih))
       units=trim(mom_units(imom))//' m'
       call save_dg(pptice_2d, name, i_dgtime,  units, dim='time')
       !snow ppt
       ih=4
       name='surface_ppt_for_'//trim(h_names(ih))
       units=trim(mom_units(imom))//' m'
       call save_dg(pptsnow_2d, name, i_dgtime,  units, dim='time')
       !graupel ppt
       ih=5
       name='surface_ppt_for_'//trim(h_names(ih))
       units=trim(mom_units(imom))//' m'
       call save_dg(pptgraul_2d, name, i_dgtime,  units, dim='time')
       !total ppt
       name='total_surface_ppt'
       units=trim(mom_units(imom))//' m'
       call save_dg((pptice_2d+pptrain_2d+pptsnow_2d+pptgraul_2d), &
            name, i_dgtime, units, dim='time') 
       !save precip flux at all levels and columns
       name='total_ppt_level'
       units=trim(mom_units(imom))//' m'
       call save_dg(pptrain_2d_prof, name, i_dgtime,  units, dim='z,x')
    endif
    
  end Subroutine mphys_thompson09_interfacen
     
end module mphys_thompson09n
