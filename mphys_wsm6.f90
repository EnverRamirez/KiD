! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing interface to greg thompson's microphysics code 
!
! Mudei RKIND para R8SIZE, que é como está no KID
! adicionei USE typeKind
! coloque como é no KID
! len=512 (pq é o valor do MPAS)

module mphys_wsm6

   Use parameters, only : num_h_moments, num_h_bins, nspecies, nz, dt &
        , h_names, mom_units, max_char_len, nx!, dtcld
   Use column_variables
   Use physconst, only : p0, r_on_cp, pi, g, cp, cv, R, tk0c, ep1, ep2 &
         , rhoair0, cice, rhoW, epsilon, psat, xlf, xls, xlv0
 
   Use module_mp_wsm6_2024
   Use diagnostics, only: save_dg, i_dgtime
   Use namelists, only: iiwarm, hail_opt
   
 
   Implicit None
 
   !Logical switches 
   logical :: micro_unset=.True.
   integer:: ih, imom
   character(max_char_len) :: name, units
 
   REAL(wp), PARAMETER, PRIVATE:: qmin = 1.E-15
   real(WP),parameter:: cliq    = 4190.
 contains
   
   Subroutine mphys_wsm6_interface
     
     real :: pptrain, pptsnow, pptgraul, pptice
     real :: rainncv (nx), snowncv (nx), graupelncv (nx)
     !real :: snow, graupel, sr
     real :: sr (nx)
     real :: pptrain_2d(nx), pptsnow_2d(nx), pptgraul_2d(nx) &
          , pptice_2d(nx), pptrain_level(nz), pptrain_2d_prof(nz,nx), &
          rain (nx), snow (nx), graupel (nx)
     real :: t1d(nz), p1d(nz), dz1d(nz),qv1d(nz),qc1d(nz) &
          ,qr1d(nz), qi1d(nz), ni1d(nz), qs1d(nz), qg1d(nz) &
          ,nr1d(nz), dzq(nz), delz(nz,nx), dz_p(nz)
      real :: th(nz),qc(nz,nx), &
             qr(nz,nx), qi(nz,nx), qs(nz,nx), qg(nz,nx), p(nz)
          !,nr1d(nz), dzq(nz)
      real (8) :: t_hv(nz),qc_hv(nz,nx), qv_hv(nz,nx), &
             qr_hv(nz,nx), qi_hv(nz,nx), qs_hv(nz,nx), qg_hv(nz,nx) &
            ,p_hv(nz,nx), delz_hv(nz,nx)
     real :: delt
     real :: dtcld
     real(4) :: q(nz,nx)
     real :: rvs
     real :: nc1d(nz), nifa1d(nz), nwfa1d(nz), w1d(nz)
     real :: rainprod1d(nz), evapprod1d(nz)
     real :: rainprod(nz), evapprod(nz)
     logical :: diagflag
     integer :: do_radar_ref
     integer :: errflg
     character :: errmsg
     integer :: has_reqc,has_reqi,has_reqs
     integer :: ids,ide,jds,jde,kds,kde, &
                ims,ime,jms,jme,kms,kme, &
                its,ite,jts,jte
     ! real :: re_cloud, re_ice, re_snow - ORIGINAL
     real :: re_qc_bg,re_qi_bg,re_qs_bg,re_qc_max,re_qi_max,re_qs_max
     real :: refl_10cm (nz, nx), re_cloud (nz, nx), re_ice (nz, nx), re_snow(nz,nx)
     integer :: kts, kte, i, j, k
 
     real(8) :: rho(nz), R_v, xlf0, xlf
     real(8) :: den (nz), rv
     real(8) :: cpd, rd, cpv, t0c
     real(8) :: pii(nz,nx)
     !real(8) :: pii
     logical :: dryrho
     
 
     real(8) :: den0, denr, dens, cl
     !integer :: hail_opt
 
 
     
     kts=1
     kte=nz
     j=1
 
    !qc1d(:) = 0.0
    !qr1d(:) = 0.0
    !nr1d(:) = 0.0
    !qi1d(:) = 0.0
    !ni1d(:) = 0.0
    !qs1d(:) = 0.0
    !qg1d(:) = 0.0
 
     !qc_hv(:) = 0.0
     !qr_hv(:) = 0.0
     !nr1d(:) = 0.0
     !qi_hv(:) = 0.0
     !ni1d(:) = 0.0
     !qs_hv(:) = 0.0
     !qg_hv(:) = 0.0
 
     !qc_hv(:,:) = 0.0
     !qr_hv(:,:) = 0.0
     !nr1d(:) = 0.0
     !qi_hv(:,:) = 0.0
     !ni1d(:) = 0.0
     !qs_hv(:,:) = 0.0
     !qg_hv(:,:) = 0.0

     qc(:,:) = 0.0
     qr(:,:) = 0.0
     qi(:,:) = 0.0
     qs(:,:) = 0.0
     qg(:,:) = 0.0
 
 !Era dt, mudei para dt - Marilia 25/03/2024
     !dt = dtcld ! - comentei 24/04/2024
 
     do i=1,nx
        pptrain = 0. !taxa de precipitacao de chuva
        pptsnow = 0. !taxa e precipitacao de neve
        pptgraul = 0. !taxa de precipitacao de graupel
        pptice = 0. !taxa de precipitacao de gelo
        do k=1,nz
           th(k) = (theta(k,i) + (dtheta_adv(k,i)+dtheta_div(k,i))*dt )*exner(k,i) !Temp potencial em um determinado nivel vertical e ponto horizontal
           p(k) = p0*exner(k,i)**(1./r_on_cp) ! pressao em um determinado nivel vertical e ponto horizontal
 !          dz1d(k) = dz(k) ! espessura da camada em um determinado nivel vertical e ponto horizontal 
 !          dzq(k) = dz(k) ! espessura da camada em um determinado nivel vertical e ponto horizontal 
 !          delz_hv(k) = dz(k) ! espessura da camada em um determinado nivel vertical e ponto horizontal - 25/03/2024
           delz_hv(k,i) = delz(k,i) ! espessura da camada em um determinado nivel vertical e ponto horizontal - 25/03/2024
           q(k,i) = qv(k,i) + (dqv_adv(k,i)+dqv_div(k,i))*dt !umidade especifica em um determinado nivel vertical e ponto horizontal 
           
           qc(k,i) = hydrometeors(k,i,1)%moments(1,1) &       !Conteudo de agua condensada (nuvem) em um
                + (dhydrometeors_adv(k,i,1)%moments(1,1) &    !determinado nivel da vertical e ponto
                + dhydrometeors_div(k,i,1)%moments(1,1))*dt   !horizontal
           
           qr(k,i) = hydrometeors(k,i,2)%moments(1,1) &       !conteudo de agua na forma de chuva em um 
                + (dhydrometeors_adv(k,i,2)%moments(1,1) &    !determinado nivel da vertical e ponto
                + dhydrometeors_div(k,i,2)%moments(1,1))*dt   !horizontal
           
         !temos que voltar aqui - DUVIDA!
          ! nr1d(k) = hydrometeors(k,i,2)%moments(1,2) &       !conteudo de numero de particulas de chuva
          !      + (dhydrometeors_adv(k,i,2)%moments(1,2) &    !em um determinado nivel da vertical e
          !      + dhydrometeors_div(k,i,2)%moments(1,2))*dt   !ponto na horizontal
           
           if (.not. iiwarm) then   !Se nao for um ambiente quente
              qi(k,i) = hydrometeors(k,i,3)%moments(1,1) &     !conteudo de agua na forma de gelo em um
                   + (dhydrometeors_adv(k,i,3)%moments(1,1) &    !determinado nivel da vertical e ponto
                   + dhydrometeors_div(k,i,3)%moments(1,1))*dt    !horizontal
                     
                   !temos que voltar aqui - DUVIDA!
             ! ni1d(k) = hydrometeors(k,i,3)%moments(1,2) &      !conteudo de numero de particulas de gelo 
             !      + (dhydrometeors_adv(k,i,3)%moments(1,2) &    !em um determinado nuvel da vertical e
             !      + dhydrometeors_div(k,i,3)%moments(1,2))*dt   !ponto na horizontal 
          
              qs(k,i) = hydrometeors(k,i,4)%moments(1,1) &      !conteudo de agua na forma de neve em um
                   + (dhydrometeors_adv(k,i,4)%moments(1,1) &    !determinado nivel da vertical e ponto
                   + dhydrometeors_div(k,i,4)%moments(1,1))*dt    !horizontal 
                if (.not. hail_opt) then
              qg(k,i) = hydrometeors(k,i,5)%moments(1,1) &      !conteudo de agua na forma de graupel em um
                   + (dhydrometeors_adv(k,i,5)%moments(1,1) &    !determinado nivel da vertical e ponto
                   + dhydrometeors_div(k,i,5)%moments(1,1))*dt    !horizontal 
                end if
           end if
           
        end do
 
 
        
        ! Initialise microphysics 
        if (micro_unset)then
           call mp_wsm6_init!(den0,denr,dens,cl,cpv,errmsg,errflg) !colocar as variaveis
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
 
 
 !       call wsm6(qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, ni1d, &
 !       nr1d, nc1d, nwfa1d, nifa1d, t1d, p1d, w1d, dzq, &
 !       pptrain, pptsnow, pptgraul, pptice, &
  ! KiD specific diagnostic  
  !      pptrain_level,                                        &
  ! End KiD diag   
 !#ifdef WRF_CHEM                  
 !       rainprod, evapprod, &
 !#endif
 !       kts, kte, dt, i, j)

        !Marilia 10/05/2024
        den = rho
        !pii = exner
        cpd = cp
        rv = R_v 
        R_v = 461.6
        xlf0 = xlf
        xlf     = 3.50e5
        t0c = 273.15
        rd = R
!Marilia 10/05/2024

 
       ! dryrho = .false.
 
        !if (.not. dryrho) then
        ! rho = 0.622*p/(287*th*(qv(:,:)+0.622)) !perfil vertical da densidade do ar
       !  !rho = 0.622*p/(R*th*(qv(:,:)+0.622)) !perfil vertical da densidade do ar
       ! else 
       !  rho = p/(R*th) 
       ! endif

 
        !th = t1d/exner(:,i)
        
 
        !call wsm6(t1d/exner(:,i),qv1d,qc1d,qr1d, &
        !call wsm6(th,qv1d,qc1d,qr1d, &
        !          !qi1d,qs1d,qg1d,rho,exner(:,i), &
        !          qi1d,qs1d,qg1d,rho,exner, &
        !          p1d, dzq, delt, g, cp, cv, R, rvs, tk0c, &
        !          ep1, ep2, epsilon, xls, xlv0, xlf, &
        !          rhoair0, rhoW, cliq, cice, psat, &
        !          !rain,pptrain_2d,snow,pptsnow_2d,graupel,pptgraul_2d,sr, &
        !          rain,rainncv,snow,snowncv,graupel,graupelncv,sr, &
        !          refl_10cm,diagflag,do_radar_ref,                 &
        !          has_reqc,has_reqi,has_reqs,re_qc_bg,re_qi_bg,    &
        !          re_qs_bg,re_qc_max,re_qi_max,re_qs_max,          &
        !          re_cloud,re_ice,re_snow,                         &
        !          ids,ide,jds,jde,kds,kde,                         &
        !          ims,ime,jms,jme,kms,kme,                         &
        !          its,ite,jts,jte,kts,kte,                         &
        !          errmsg,errflg                                    &
!#if(WRF_CHEM == 1)
!                  ,wetscav_on,evapprod,rainprod                    &
!#endif
!                 )
         call wsm6(th,q,qc,qr,qi,qs,qg,den,pii,p,delz,              &
         delt,g,cpd,cpv,rd,rv,t0c,ep1,ep2,qmin,           &
         xls,xlv0,xlf0,den0,denr,cliq,cice,psat,          &
         rain,rainncv,snow,snowncv,graupel,graupelncv,sr, &
         refl_10cm,diagflag,do_radar_ref,                 &
         has_reqc,has_reqi,has_reqs,re_qc_bg,re_qi_bg,    &
         re_qs_bg,re_qc_max,re_qi_max,re_qs_max,          &
         re_cloud,re_ice,re_snow,                         &
         ids,ide,jds,jde,kds,kde,                         &
         ims,ime,jms,jme,kms,kme,                         &
         its,ite,jts,jte,kts,kte,                         &
         errmsg,errflg                                    &
#if(WRF_CHEM == 1)
     ,wetscav_on,evapprod,rainprod                    &
#endif
    )
 
                    
 
 
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
        else !NAO TENHO PPTICE
           pptrain_2d(i) = pptrain
           pptice_2d(i) = pptice
           pptsnow_2d(i) = pptsnow
           pptgraul_2d(i) = pptgraul
           !
           ! pptrain_level is declared and calced in 
           ! the module module_mp_thompson09
           !pptrain_2d_prof(:,i) = pptrain_level(:)
        endif
        
        !do k=1,nz
        ! dtheta_mphys(k,i)=(th(k)/exner(k,i)-theta(k,i))/dt
        ! dqv_mphys(k,i)=(q(k,i) - qv(k,i))/dt
        ! dhydrometeors_mphys(k,i,1)%moments(1,1)= (qc_hv(k,i)-hydrometeors(k,i,1)%moments(1,1))/dt
        ! dhydrometeors_mphys(k,i,2)%moments(1,1)= (qr_hv(k,i)-hydrometeors(k,i,2)%moments(1,1))/dt
        ! !dhydrometeors_mphys(k,i,2)%moments(1,2)= nr_tend1d(k)
        ! dhydrometeors_mphys(k,i,3)%moments(1,1)= (qi_hv(k,i)-hydrometeors(k,i,3)%moments(1,1))/dt
        ! !dhydrometeors_mphys(k,i,3)%moments(1,2)= ni_tend1d(k)
        ! dhydrometeors_mphys(k,i,4)%moments(1,1)= (qs_hv(k,i)-hydrometeors(k,i,4)%moments(1,1))/dt
        ! !dhydrometeors_mphys(k,i,4)%moments(1,2)= ns_tend1d(k)
        ! dhydrometeors_mphys(k,i,5)%moments(1,1)= (qg_hv(k,i)-hydrometeors(k,i,5)%moments(1,1))/dt
        ! !dhydrometeors_mphys(k,i,5)%moments(1,2)= ng_tend1d(k)
      !end do
 
     ! back out tendencies
        !Essa parte do código está calculando as tendências das variáveis 
        !de microfísica de nuvens e precipitação ao longo da dimensão vertical 

        do k=1,nz
        !do i=1,nx
           dtheta_mphys(k,i)=(th(k)/exner(k,i)-theta(k,i))/dt &  !Tendência da temperatura potencial para um nível vertical específico e ponto horizontal. 
                - ( dtheta_adv(k,i)+dtheta_div(k,i))    !é calculada com base na diferença entre a temperatura potencial atual (theta(k,i)) 
                                                        !e a temperatura potencial calculada a partir das variáveis t1d, exner, dtheta_adv, dtheta_div, dt
                                                        !O resultado é ajustado pela contribuição das tendências advectivas e divergentes.
                                                        
           
           dqv_mphys(k,i)=(qv_hv(k,i) - q(k,i))/dt & !Tendência da umidade específica para um nível vertical específico e ponto horizontal. 
                - ( dqv_adv(k,i)+dqv_div(k,i))     !é calculada com base na diferença entre a umidade específica atual (qv(k,i))
                                                   !e a umidade específica calculada a partir das variáveis qv1d, dt, dqv_adv, dqv_div
                                                   !O resultado é ajustado pela contribuição das tendências advectivas e divergentes.
                                                   ! Era qv(k,i), mudei para qv_hv(k,i)
           dhydrometeors_mphys(k,i,1)%moments(1,1)= &   !Tendência do momento da fase de condensado (nuvem) para um nível vertical específico e ponto horizontal. 
                (qc_hv(k,i)-hydrometeors(k,i,1)%moments(1,1))/dt & !A tendência é calculada com base na diferença entre o conteúdo de água condensada atual (qc1d(k)) 
               - (dhydrometeors_adv(k,i,1)%moments(1,1)  & ! e o conteúdo de água condensada da fase de condensado calculado a partir das variáveis hydrometeors,
                + dhydrometeors_div(k,i,1)%moments(1,1)) ! dt, dhydrometeors_adv, dhydrometeors_div. O resultado é ajustado pela contribuição das tendências advectivas e divergentes.
                                                   ! Era qc(k,i), mudei para qc_hv(k,i) 
           dhydrometeors_mphys(k,i,2)%moments(1,1)= & !Tendência do momento da fase de chuva para um nível vertical específico e ponto horizontal.
                (qr_hv(k,i)-hydrometeors(k,i,2)%moments(1,1))/dt & !Similar à dhydrometeors_mphys(k,i,1)%moments(1,1), mas para a fase de chuva.
                - (dhydrometeors_adv(k,i,2)%moments(1,1)  &
                + dhydrometeors_div(k,i,2)%moments(1,1))
                                                   ! Era qr(k,i), mudei para qr_hv(k,i)
          ! dhydrometeors_mphys(k,i,2)%moments(1,2)= & !Tendência do momento do número de partículas de chuva para um nível vertical específico e ponto horizontal.
          !      (nr1d(k)-hydrometeors(k,i,2)%moments(1,2))/dt & !Similar à dhydrometeors_mphys(k,i,2)%moments(1,1), mas para o número de partículas de chuva.
          !      - (dhydrometeors_adv(k,i,2)%moments(1,2)  &
          !      + dhydrometeors_div(k,i,2)%moments(1,2))
 
           if (.not.iiwarm)then 
             dhydrometeors_mphys(k,i,3)%moments(1,1)= &
                (qi_hv(k,i)-hydrometeors(k,i,3)%moments(1,1))/dt & !gelo = indice 3
                - (dhydrometeors_adv(k,i,3)%moments(1,1)  &
                + dhydrometeors_div(k,i,3)%moments(1,1))
                                                   ! Era qi(k,i), mudei para qi_hv(k,i) 
            ! dhydrometeors_mphys(k,i,3)%moments(1,2)= &
            !    (ni1d(k)-hydrometeors(k,i,3)%moments(1,2))/dt & 
            !    - (dhydrometeors_adv(k,i,3)%moments(1,2)  &
            !    + dhydrometeors_div(k,i,3)%moments(1,2))
 
             dhydrometeors_mphys(k,i,4)%moments(1,1)= &
                (qs_hv(k,i)-hydrometeors(k,i,4)%moments(1,1))/dt & !neve = indice 4
                - (dhydrometeors_adv(k,i,4)%moments(1,1)  &
                + dhydrometeors_div(k,i,4)%moments(1,1))
                                                   ! Era qs(k,i), mudei para qs_hv(k,i)
                !if (.not.hail_opt)then !Marilia 26/04/2024
             dhydrometeors_mphys(k,i,5)%moments(1,1)= &
                (qg_hv(k,i)-hydrometeors(k,i,5)%moments(1,1))/dt & !graupel = indice 5
                - (dhydrometeors_adv(k,i,5)%moments(1,1)  &
                + dhydrometeors_div(k,i,5)%moments(1,1))
                                                   ! Era qg(k,i), mudei para qg_hv(k,i)
                !end if
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
     
   end Subroutine mphys_wsm6_interface
      
 end module mphys_wsm6 