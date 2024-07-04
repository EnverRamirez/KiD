! wsm6
! MPAS to KiD
! Marilia
! fev/2024

      module module_mp_wsm6_2024
      !use mpas_log !continua... pq eh um modulo de mensagens de erro. nao tem calculos
      !use mpas_kind_types,only: RKIND
      use typeKind

      !use mp_wsm6,only: mp_wsm6_run => coloquei dentro desse código
      !use mp_wsm6_effectrad,only: mp_wsm6_effectRad_run
      Use namelists, only: iiwarm, hail_opt
      Use switches, only: l_sediment
     !Use mp_radar
      Use module_libmassv
      Use diagnostics, only: save_dg, i_dgtime
      Use parameters, only: nx, dt

      implicit none
      private
      public:: wsm6, &
               mp_wsm6_run,   &
               mp_wsm6_init,  &
               mp_wsm6_final, &
               refl10cm_wsm6

      

      

 !MARILIA

      INTEGER, PARAMETER, PRIVATE:: R8SIZE = 8 !Marilia
      INTEGER, PARAMETER, PRIVATE:: R4SIZE = 4 ! Marilia

      real(kind=R8SIZE),public:: dt_microp

            !troquei (KIND=R8SIZE) por kind=r8size, que é o padrão do KiD

      real(KIND=R8SIZE),parameter,private:: dtcldcr = 120.    ! maximum time step for minor loops
      real(KIND=R8SIZE),parameter,private:: n0r = 8.e6        ! intercept parameter rain
!real(KIND=R8SIZE),parameter,private:: n0g = 4.e6        ! intercept parameter graupel
      real(KIND=R8SIZE),parameter,private:: avtr = 841.9      ! a constant for terminal velocity of rain
      real(KIND=R8SIZE),parameter,private:: bvtr = 0.8        ! a constant for terminal velocity of rain
      real(KIND=R8SIZE),parameter,private:: r0 = .8e-5        ! 8 microm  in contrast to 10 micro m
      real(KIND=R8SIZE),parameter,private:: peaut = .55       ! collection efficiency
      real(KIND=R8SIZE),parameter,private:: xncr = 3.e8       ! maritime cloud in contrast to 3.e8 in tc80
      real(KIND=R8SIZE),parameter,private:: xmyu = 1.718e-5   ! the dynamic viscosity kgm-1s-1
      real(KIND=R8SIZE),parameter,private:: avts = 11.72      ! a constant for terminal velocity of snow
      real(KIND=R8SIZE),parameter,private:: bvts = .41        ! a constant for terminal velocity of snow
      !real(KIND=R8SIZE),parameter,private:: avtg = 330.       ! a constant for terminal velocity of graupel
      !real(KIND=R8SIZE),parameter,private:: bvtg = 0.8        ! a constant for terminal velocity of graupel
      !real(KIND=R8SIZE),parameter,private:: deng = 500.       ! density of graupel ! set later with hail_opt
      real(KIND=R8SIZE),parameter,private:: lamdarmax = 8.e4  ! limited maximum value for slope parameter of rain
      real(KIND=R8SIZE),parameter,private:: lamdasmax = 1.e5  ! limited maximum value for slope parameter of snow
      !real(KIND=R8SIZE),parameter,private:: lamdagmax = 6.e4  ! limited maximum value for slope parameter of graupel
      real(KIND=R8SIZE),parameter,private:: dicon = 11.9      ! constant for the cloud-ice diamter
      real(KIND=R8SIZE),parameter,private:: dimax = 500.e-6   ! limited maximum value for the cloud-ice diamter
      real(KIND=R8SIZE),parameter,private:: pfrz1 = 100.      ! constant in Biggs freezing
      real(KIND=R8SIZE),parameter,private:: pfrz2 = 0.66      ! constant in Biggs freezing
      real(KIND=R8SIZE),parameter,private:: qcrmin = 1.e-9    ! minimun values for qr, qs, and qg
      real(KIND=R8SIZE),parameter,private:: eacrc = 1.0       ! Snow/cloud-water collection efficiency
      real(KIND=R8SIZE),parameter,private:: dens  =  100.0    ! Density of snow
      real(KIND=R8SIZE),parameter,private:: qs0   =  6.e-4    ! threshold amount for aggretion to occur

      real(KIND=R8SIZE),parameter,public :: n0smax =  1.e11   ! maximum n0s (t=-90C unlimited)
      real(KIND=R8SIZE),parameter,public :: n0s = 2.e6        ! temperature dependent intercept parameter snow
      real(KIND=R8SIZE),parameter,public :: alpha = .12       ! .122 exponen factor for n0s

      REAL(WP), PARAMETER :: RV     = 461.6
      real(wp), parameter :: denr   = 1000.  
      real(wp), parameter :: den0   = 1.28
      !real(wp), parameter :: dens   = 100.
      real(wp), parameter :: cl     = 4190.
      real(wp), parameter :: cpv    = 1846.4

      real(KIND=R8SIZE),save::                 &
         qc0,qck1,                                &
         bvtr1,bvtr2,bvtr3,bvtr4,g1pbr,           &
         g3pbr,g4pbr,g5pbro2,pvtr,eacrr,pacrr,    &
         bvtr6,g6pbr,                             &
         precr1,precr2,roqimax,bvts1,             &
         bvts2,bvts3,bvts4,g1pbs,g3pbs,g4pbs,     &
         n0g,avtg,bvtg,deng,lamdagmax,            & !RAS13.3 - set these in wsm6init
         g5pbso2,pvts,pacrs,precs1,precs2,pidn0r, &
         xlv1,pacrc,pi,                           &
         bvtg1,bvtg2,bvtg3,bvtg4,g1pbg,           &
         g3pbg,g4pbg,g5pbgo2,pvtg,pacrg,          &
         precg1,precg2,pidn0g,                    &
         rslopermax,rslopesmax,rslopegmax,        &
         rsloperbmax,rslopesbmax,rslopegbmax,     &
         rsloper2max,rslopes2max,rslopeg2max,     &
         rsloper3max,rslopes3max,rslopeg3max

      real(KIND=R8SIZE),public,save:: pidn0s,pidnc

      !vem do modulo mp_radar

      integer,parameter,public:: nrbins = 50
integer,parameter,public:: slen = 20
character(len=slen), public:: mixingrulestring_s, matrixstring_s, &
            inclusionstring_s, &
            hoststring_s, hostmatrixstring_s, hostinclusionstring_s, &
            mixingrulestring_g, matrixstring_g, inclusionstring_g,   &
            hoststring_g, hostmatrixstring_g, hostinclusionstring_g

complex(KIND=R8SIZE),public:: m_w_0, m_i_0

double precision,dimension(nrbins+1),public:: xxdx
double precision,dimension(nrbins),public:: xxds,xdts,xxdg,xdtg
double precision,parameter,public:: lamda_radar = 0.10           ! in meters
double precision,public:: k_w,pi5,lamda4

double precision, dimension(nrbins+1), public:: simpson
double precision, dimension(3), parameter, public:: basis = (/1.d0/3.d0, 4.d0/3.d0, 1.d0/3.d0/)

real(KIND=R8SIZE),public,dimension(4):: xcre,xcse,xcge,xcrg,xcsg,xcgg
real(KIND=R8SIZE),public:: xam_r,xbm_r,xmu_r,xobmr
real(KIND=R8SIZE),public:: xam_s,xbm_s,xmu_s,xoams,xobms,xocms
real(KIND=R8SIZE),public:: xam_g,xbm_g,xmu_g,xoamg,xobmg,xocmg
real(KIND=R8SIZE),public:: xorg2,xosg2,xogg2


!..Single melting snow/graupel particle 90% meltwater on external sfc
character(len=256):: radar_debug

double precision,parameter,public:: melt_outside_s = 0.9d0
double precision,parameter,public:: melt_outside_g = 0.9d0

!fim do que veio do modulo mp_radar
      !MARILIA
      contains
      !=================================================================================================================
      subroutine wsm6(th,q,qc,qr,qi,qs,qg,den,pii,p,delz,              &
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
!=================================================================================================================
       
!--- input arguments:
      logical,intent(in),optional:: diagflag

      integer,intent(in):: ids,ide,jds,jde,kds,kde, &
                           ims,ime,jms,jme,kms,kme, &
                           its,ite,jts,jte,kts,kte

      integer,intent(in):: has_reqc,has_reqi,has_reqs
      integer,intent(in),optional:: do_radar_ref

      !INTEGER, PARAMETER, PRIVATE:: R8SIZE = 8 !Marilia
      !INTEGER, PARAMETER, PRIVATE:: R4SIZE = 4 ! Marilia

      !troquei (kind=r8kind) por kind=r8size, que é o padrão do KiD

      real(KIND=R4SIZE),intent(in) :: delt!, rv

     ! real(KIND=R8SIZE),intent(in):: &
     !    g,rd,rv,t0c,den0,cpd,cpv,ep1,ep2,qmin,xls,xlv0,xlf0, &
     !    cliq,cice,psat,denr

      real(KIND=R8SIZE),intent(in):: &
         g,rd,rv,t0c,den0,cpd,cpv,ep1,ep2,qmin,xls,xlv0,xlf0, &
         cliq,cice,psat,denr

      real(KIND=R4SIZE),intent(in):: &  !era r8size - mudei dia 15/03/2024
         re_qc_bg,re_qi_bg,re_qs_bg,re_qc_max,re_qi_max,re_qs_max

      real(KIND=R4SIZE),intent(in),dimension(ims:ime,kms:kme,jms:jme )::&
         !den, &
         !pii, &
         p, &
         delz
      real(KIND=R8SIZE),intent(in),dimension(ims:ime,kms:kme,jms:jme )::& !Marilia
         den!, &
         !p, &
         !delz!, &
         !delt

      real(KIND=R8SIZE),intent(in),dimension(ims:ime,kms:kme,jms:jme )::&
        pii
       

      !inout arguments:
      !real(KIND=R8SIZE),intent(inout),dimension(ims:ime,jms:jme):: &
      !   rain,rainncv,sr

      real(KIND=R4SIZE),intent(inout),dimension(ims:ime,jms:jme):: &
         rain, rainncv

      real(KIND=R4SIZE),intent(inout),dimension(ims:ime,jms:jme):: &
         sr

      real(KIND=R4SIZE),intent(inout),dimension(ims:ime,jms:jme),optional:: &
         snow,snowncv

      real(KIND=R4SIZE),intent(inout),dimension(ims:ime,jms:jme),optional:: &
         graupel,graupelncv


!         real(KIND=R8SIZE),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: &
!         th, &
      real(KIND=R4SIZE),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: & !Marilia
!         th, & - Marilia
         q,  & !razao de mistura do vapor dagua
         qc, & !razao de mistura de agua de nuvem
         qi, & !razao de mistura de cristais de gelo
         qr, & !razao de mistura de chuva
         qs, & !razao de mistura de neve
         qg    !razao de mistura de graupel

         real(KIND=R4SIZE),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: &
         th !Marilia

      real(KIND=R4SIZE),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: & !Era r8size - mudei 15/03/2024
         re_cloud, &
         re_ice,   &
         re_snow

      real(KIND=R4SIZE),intent(inout),dimension(ims:ime,kms:kme,jms:jme),optional:: & !Era r8size - mudei 14/03/2024
         refl_10cm

#if(WRF_CHEM == 1) 
      logical,intent(in):: wetscav_on
      real(KIND=R8SIZE),intent(inout),dimension(ims:ime,kms:kme,jms:jme ):: &
         rainprod,evapprod
#endif

!output arguments:
      character(len=*),intent(out):: errmsg
      integer,intent(out):: errflg

!local variables and arrays:
      logical:: do_microp_re
      integer:: i,j,k

      real(KIND=R8SIZE),dimension(kts:kte):: qv1d,t1d,p1d,qr1d,qs1d,qg1d,dBZ
      real(KIND=R8SIZE),dimension(kts:kte):: den1d,qc1d,qi1d,re_qc,re_qi,re_qs
      !real(KIND=R8SIZE),dimension(kts:kte):: t_hv,den_hv,p_hv,delz_hv !adicionei 26/03/2024 - Marilia - voltei para (i,k) - 25/04/2024
      !real(KIND=R8SIZE),dimension(kts:kte):: qv_hv,qc_hv,qi_hv,qr_hv,qs_hv,qg_hv !adicionei 26/03/2024 - Marilia- voltei para (i,k) - 25/04/2024
      !real(KIND=R8SIZE),dimension(kts:kte):: re_qc_hv,re_qi_hv,re_qs_hv !adicionei 26/03/2024 - Marilia- voltei para (i,k) - 25/04/2024

      real(KIND=R8SIZE),dimension(its:ite):: rainncv_hv,rain_hv,sr_hv
      real(KIND=R8SIZE),dimension(its:ite):: snowncv_hv,snow_hv
      real(KIND=R8SIZE),dimension(its:ite):: graupelncv_hv,graupel_hv
      real(KIND=R8SIZE),dimension(its:ite,kts:kte):: t_hv,den_hv,p_hv,delz_hv !descomentei 25/04/2024
      real(KIND=R8SIZE),dimension(its:ite,kts:kte):: qv_hv,qc_hv,qi_hv,qr_hv,qs_hv,qg_hv !descomentei 25/04/2024
      real(KIND=R8SIZE),dimension(its:ite,kts:kte):: re_qc_hv,re_qi_hv,re_qs_hv !descomentei 25/04/2024

      !dt=delt
!-----------------------------------------------------------------------------------------------------------------

      do j = jts,jte

         do i = its,ite
            !input arguments:
            do k = kts,kte
               den_hv(i,k)  = den(i,k,j) ! era den_hv(i,k) - Marilia 26/03/2024 - voltei para (i,k) - 25/04/2024
               p_hv(i,k)    = p(i,k,j) ! era p_hv(i,k)  - Marilia 26/03/2024- voltei para (i,k) - 25/04/2024
               delz_hv(i,k) = delz(i,k,j) !era delz_hv(i,k) - Marilia 26/03/2024- voltei para (i,k) - 25/04/2024
            enddo

            !inout arguments:
            rain_hv(i) = rain(i,j)

            do k = kts,kte
               t_hv(i,k)  = th(i,k,j)*pii(i,k,j) ! era t_hv(i,k)  - Marilia 26/03/2024- voltei para (i,k) - 25/04/2024
               qv_hv(i,k) = q(i,k,j) ! era qv_hv(i,k) - Marilia 26/03/2024- voltei para (i,k) - 25/04/2024
               qc_hv(i,k) = qc(i,k,j) !era qc_hv(i,k) - Marilia 26/03/2024- voltei para (i,k) - 25/04/2024
               qi_hv(i,k) = qi(i,k,j) !era qi_hv(i,k) - Marilia 26/03/2024- voltei para (i,k) - 25/04/2024
               qr_hv(i,k) = qr(i,k,j) !era qr_hv(i,k) - Marilia 26/03/2024- voltei para (i,k) - 25/04/2024
               qs_hv(i,k) = qs(i,k,j) !era qs_hv(i,k) - Marilia 26/03/2024- voltei para (i,k) - 25/04/2024
               qg_hv(i,k) = qg(i,k,j) !era qg_hv(i,k) - Marilia 26/03/2024- voltei para (i,k) - 25/04/2024
            enddo
         enddo

         if(present(snow) .and. present(snowncv)) then
            do i = its,ite
               snow_hv(i) = snow(i,j)
            enddo
         endif
         if(present(graupel) .and. present(graupelncv)) then
            do i = its,ite
               graupel_hv(i) = graupel(i,j)
            enddo
         endif

!--- call to cloud microphysics scheme:
      call mp_wsm6_run(t=t_hv,q=qv_hv,qc=qc_hv,qi=qi_hv,qr=qr_hv,qs=qs_hv,qg=qg_hv, &
                 den=den_hv,p=p_hv,delz=delz_hv,delt=dt_microp,g=g,cpd=cpd,cpv=cpv,    &
                 rd=rd,rv=rv,t0c=t0c,ep1=ep1,ep2=ep2,qmin=qmin,xls=xls,xlv0=xlv0, &
                 xlf0=xlf0,den0=den0,denr=denr,cliq=cliq,cice=cice,psat=psat,     &
                 rain=rain_hv,rainncv=rainncv_hv,sr=sr_hv,snow=snow_hv,&
                 snowncv=snowncv_hv,graupel=graupel_hv,graupelncv=graupelncv_hv,  &
                 its=its,ite=ite,kts=kts,kte=kte,errmsg=errmsg,errflg=errflg      &
#if(WRF_CHEM == 1)
                 ,rainprod2d=rainprod_hv,evapprod2d=evapprodhv         &
#endif
                    )

      do i = its,ite
         !inout arguments:
         rain(i,j)    = rain_hv(i)
         rainncv(i,j) = rainncv_hv(i)
         sr(i,j)      = sr_hv(i)
      
         do k = kts,kte
            th(i,k,j) = t_hv(i,k)/pii(i,k,j) !era t_hv(i,k)/pii(i,k,j)  - Marilia 26/03/2024 - voltei para (i,k) 25/04/2024
            q(i,k,j)  = qv_hv(i,k) ! era qv_hv(i,k)  - Marilia 26/03/2024 - voltei para (i,k) 25/04/2024
            qc(i,k,j) = qc_hv(i,k) ! era qc_hv(i,k) - Marilia 26/03/2024 - voltei para (i,k) 25/04/2024
            qi(i,k,j) = qi_hv(i,k) ! era qi_hv(i,k) - Marilia 26/03/2024 - voltei para (i,k) 25/04/2024
            qr(i,k,j) = qr_hv(i,k) ! era qr_hv(i,k) - Marilia 26/03/2024 - voltei para (i,k) 25/04/2024
            qs(i,k,j) = qs_hv(i,k) ! era qs_hv(i,k) - Marilia 26/03/2024 - voltei para (i,k) 25/04/2024
            qg(i,k,j) = qg_hv(i,k) ! era qg_hv(i,k) - Marilia 26/03/2024 - voltei para (i,k) 25/04/2024
         enddo
      enddo

      if(present(snow) .and. present(snowncv)) then
         do i = its,ite
            snow(i,j)    = snow_hv(i)
            snowncv(i,j) = snowncv_hv(i)
         enddo
      endif
      if(present(graupel) .and. present(graupelncv)) then
         do i = its,ite
            graupel(i,j)    = graupel_hv(i)
            graupelncv(i,j) = graupelncv_hv(i)
         enddo
      endif

#if(WRF_CHEM == 1)
         if(wetscav_on) then
            do k = kts,kte
               do i = its, ite
                  rainprod(i,k,j) = rainprod_hv(i,k)
                  evapprod(i,k,j) = evapprod_hv(i,k)
               enddo
            enddo
         else
            do k = kts,kte
               do i = its, ite
                  rainprod(i,k,j) = 0.
                  evapprod(i,k,j) = 0.
               enddo
            enddo
         endif
#endif

!--- call to computation of effective radii for cloud water, cloud ice, and snow:
         do_microp_re = .false.
         if(has_reqc == 1 .and. has_reqi == 1 .and. has_reqs == 1) do_microp_re = .true.

         do k = kts,kte
            do i = its,ite
               t_hv(i,k) = th(i,k,j)*pii(i,k,j) !era t_hv(i,k) - voltei para (i,k) 25/04/2024
               re_qc_hv(i,k) = re_cloud(i,k,j) !era re_qc_hv(i,k) - voltei para (i,k) 25/04/2024
               re_qi_hv(i,k) = re_ice(i,k,j) !era re_qi_hv(i,k) - voltei para (i,k) 25/04/2024
               re_qs_hv(i,k) = re_snow(i,k,j) !era re_qs_hv(i,k) - voltei para (i,k) 25/04/2024
            enddo
         enddo

         call  mp_wsm6_effectRad_run(do_microp_re,t_hv,qc_hv, &
                        qi_hv,qs_hv,den_hv,qmin,t0c, &
                        re_qc_bg,re_qi_bg,re_qs_bg,re_qc_max,re_qi_max,re_qs_max,re_qc_hv, &
                        re_qi_hv,re_qs_hv,its,ite,kts,kte,errmsg,errflg)

         do k = kts,kte
            do i = its,ite
               re_cloud(i,k,j) = re_qc_hv(i,k) !era re_qc_hv(i,k)  - Marilia 26/03/2024 - voltei para (i,k) 25/04/2024
               re_ice(i,k,j)   = re_qi_hv(i,k) !era re_qi_hv(i,k) - Marilia 26/03/2024 - voltei para (i,k) 25/04/2024
               re_snow(i,k,j)  = re_qs_hv(i,k) !era re_qs_hv(i,k) - Marilia 26/03/2024 - voltei para (i,k) 25/04/2024
            enddo
         enddo

      enddo

      end subroutine wsm6

!=================================================================================================================

!=================================================================================================================
      !subroutine mp_wsm6_init(den0,denr,dens,cl,cpv,hail_opt,errmsg,errflg) ! ORIGINAL
      subroutine mp_wsm6_init !deixei assim 26/04/2024
!=================================================================================================================
   
   !input arguments:
      !integer,intent(in):: hail_opt ! RAS !ORIGINAL
      !real(KIND=R4SIZE),intent(in):: den0,denr,dens,cl,cpv
   
   !output arguments:
      !character(len=*),intent(out):: errmsg
      !integer,intent(out):: errflg
   
!-----------------------------------------------------------------------------------------------------------------
   
! RAS13.1 define graupel parameters as graupel-like or hail-like,
!         depending on namelist option
     ! if(hail_opt .eq. 1) then !Hail!
     !    n0g       = 4.e4
     !    deng      = 700.
     !    avtg      = 285.0
     !    bvtg      = 0.8
     !    lamdagmax = 2.e4
     ! else !Graupel!
     !    n0g       = 4.e6
     !    deng      = 500
     !    avtg      = 330.0
     !    bvtg      = 0.8
     !    lamdagmax = 6.e4
     ! endif

      if(.not. hail_opt) then !Graupel!
         n0g       = 4.e6
         deng      = 500
         avtg      = 330.0
         bvtg      = 0.8
         lamdagmax = 6.e4
      else !Hail!
         n0g       = 4.e4
         deng      = 700.
         avtg      = 285.0
         bvtg      = 0.8
         lamdagmax = 2.e4
      endif
      !
      pi = 4.*atan(1.)
      xlv1 = cl-cpv
!
      qc0  = 4./3.*pi*denr*r0**3*xncr/den0  ! 0.419e-3 -- .61e-3
      qck1 = .104*9.8*peaut/(xncr*denr)**(1./3.)/xmyu*den0**(4./3.) ! 7.03
      pidnc = pi*denr/6.        ! syb
!
      bvtr1 = 1.+bvtr
      bvtr2 = 2.5+.5*bvtr
      bvtr3 = 3.+bvtr
      bvtr4 = 4.+bvtr
      bvtr6 = 6.+bvtr
      g1pbr = rgmma(bvtr1)
      g3pbr = rgmma(bvtr3)
      g4pbr = rgmma(bvtr4)            ! 17.837825
      g6pbr = rgmma(bvtr6)
      g5pbro2 = rgmma(bvtr2)          ! 1.8273
      pvtr = avtr*g4pbr/6.
      eacrr = 1.0
      pacrr = pi*n0r*avtr*g3pbr*.25*eacrr
      precr1 = 2.*pi*n0r*.78
      precr2 = 2.*pi*n0r*.31*avtr**.5*g5pbro2
      roqimax = 2.08e22*dimax**8
!
      bvts1 = 1.+bvts
      bvts2 = 2.5+.5*bvts
      bvts3 = 3.+bvts
      bvts4 = 4.+bvts
      g1pbs = rgmma(bvts1)    !.8875
      g3pbs = rgmma(bvts3)
      g4pbs = rgmma(bvts4)    ! 12.0786
      g5pbso2 = rgmma(bvts2)
      pvts = avts*g4pbs/6.
      pacrs = pi*n0s*avts*g3pbs*.25
      precs1 = 4.*n0s*.65
      precs2 = 4.*n0s*.44*avts**.5*g5pbso2
      pidn0r =  pi*denr*n0r
      pidn0s =  pi*dens*n0s
!
      pacrc = pi*n0s*avts*g3pbs*.25*eacrc
!
      bvtg1 = 1.+bvtg
      bvtg2 = 2.5+.5*bvtg
      bvtg3 = 3.+bvtg
      bvtg4 = 4.+bvtg
      g1pbg = rgmma(bvtg1)
      g3pbg = rgmma(bvtg3)
      g4pbg = rgmma(bvtg4)
      pacrg = pi*n0g*avtg*g3pbg*.25
      g5pbgo2 = rgmma(bvtg2)
      pvtg = avtg*g4pbg/6.
      precg1 = 2.*pi*n0g*.78
      precg2 = 2.*pi*n0g*.31*avtg**.5*g5pbgo2
      pidn0g =  pi*deng*n0g
      !
      rslopermax = 1./lamdarmax
      rslopesmax = 1./lamdasmax
      rslopegmax = 1./lamdagmax
      rsloperbmax = rslopermax ** bvtr
      rslopesbmax = rslopesmax ** bvts
      rslopegbmax = rslopegmax ** bvtg
      rsloper2max = rslopermax * rslopermax
      rslopes2max = rslopesmax * rslopesmax
      rslopeg2max = rslopegmax * rslopegmax
      rsloper3max = rsloper2max * rslopermax
      rslopes3max = rslopes2max * rslopesmax
      rslopeg3max = rslopeg2max * rslopegmax
      
!+---+-----------------------------------------------------------------+
!.. Set these variables needed for computing radar reflectivity.  These
!.. get used within radar_init to create other variables used in the
!.. radar module.
      xam_r = PI*denr/6.
      xbm_r = 3.
      xmu_r = 0.
      xam_s = PI*dens/6.
      xbm_s = 3.
      xmu_s = 0.
      xam_g = PI*deng/6.
      xbm_g = 3.
      xmu_g = 0.
      
      call radar_init
      
      !errmsg = 'mp_wsm6_init OK'
      !errflg = 0
      
      end subroutine mp_wsm6_init
   
   !=================================================================================================================
      subroutine mp_wsm6_final(errmsg,errflg)
   !=================================================================================================================
   
   !--- output arguments:
      character(len=*),intent(out):: errmsg
      integer,intent(out):: errflg
   
   !-----------------------------------------------------------------------------------------------------------------
   
      errmsg = 'mp_wsm6_final OK'
      errflg = 0
   
      end subroutine mp_wsm6_final
   
   !=================================================================================================================
      subroutine mp_wsm6_run(t,q,qc,qi,qr,qs,qg,den,p,delz,delt,   &
                             g,cpd,cpv,rd,rv,t0c,ep1,ep2,qmin,xls, &
                             xlv0,xlf0,den0,denr,cliq,cice,psat,   &
                             rain,rainncv,sr,snow,snowncv,graupel, &
                             graupelncv,rainprod2d,evapprod2d,     &
                             its,ite,kts,kte,errmsg,errflg         &
                          )
   !=================================================================================================================!
   !  This code is a 6-class GRAUPEL phase microphyiscs scheme (WSM6) of the
   !  Single-Moment MicroPhyiscs (WSMMP). The WSMMP assumes that ice nuclei
   !  number concentration is a function of temperature, and seperate assumption
   !  is developed, in which ice crystal number concentration is a function
   !  of ice amount. A theoretical background of the ice-microphysics and related
   !  processes in the WSMMPs are described in Hong et al. (2004).
   !  All production terms in the WSM6 scheme are described in Hong and Lim (2006).
   !  All units are in m.k.s. and source/sink terms in kgkg-1s-1.
   !
   !  WSM6 cloud scheme
   !
   !  Coded by Song-You Hong and Jeong-Ock Jade Lim (Yonsei Univ.)
   !           Summer 2003
   !
   !  Implemented by Song-You Hong (Yonsei Univ.) and Jimy Dudhia (NCAR)
   !           Summer 2004
   !
   !  further modifications :
   !        semi-lagrangian sedimentation (JH,2010),hong, aug 2009
   !        ==> higher accuracy and efficient at lower resolutions
   !        reflectivity computation from greg thompson, lim, jun 2011
   !        ==> only diagnostic, but with removal of too large drops
   !        add hail option from afwa, aug 2014
   !        ==> switch graupel or hail by changing no, den, fall vel.
   !        effective radius of hydrometeors, bae from kiaps, jan 2015
   !        ==> consistency in solar insolation of rrtmg radiation
   !        bug fix in melting terms, bae from kiaps, nov 2015
   !        ==> density of air is divided, which has not been
   !
   !  Reference) Hong, Dudhia, Chen (HDC, 2004) Mon. Wea. Rev.
   !             Hong and Lim (HL, 2006) J. Korean Meteor. Soc.
   !             Dudhia, Hong and Lim (DHL, 2008) J. Meteor. Soc. Japan
   !             Lin, Farley, Orville (LFO, 1983) J. Appl. Meteor.
   !             Rutledge, Hobbs (RH83, 1983) J. Atmos. Sci.
   !             Rutledge, Hobbs (RH84, 1984) J. Atmos. Sci.
   !             Juang and Hong (JH, 2010) Mon. Wea. Rev.
   !
   
   !input arguments:
      integer,intent(in):: its,ite,kts,kte

      real(KIND=R8SIZE),intent(in),dimension(its:ite,kts:kte)::     &  !Marilia - era R8SIZE antes voltei para (i,k) 25/04/2024
                                                                        ! Mudei 27/03/2024
      !real(KIND=R8SIZE),intent(in),dimension(kts:kte)::     &  !Marilia - era R8SIZE antes
                                                                  !den, &
                                                                    p!, &
                                                                  !delz
      real(KIND=R8SIZE),intent(in),dimension(its:ite,kts:kte)::     &  !Marilia - era R8SIZE antes
                                                                  den, &
                                                                  delz
      real(KIND=R8SIZE),intent(in) :: rv
      
      real(KIND=R8SIZE),intent(in)::                                &
                                                               delt, &
                                                                  g, &
                                                                cpd, &
                                                                cpv, &
                                                                t0c, &
                                                               den0, &
                                                                 rd, &
                                                               !  rv, &
                                                                ep1, &
                                                                ep2, &
                                                               qmin, &
                                                                xls, &
                                                               xlv0, &
                                                               xlf0, &
                                                               cliq, &
                                                               cice, &
                                                               psat, &
                                                               denr
   
   !inout arguments:
      real(KIND=R8SIZE),intent(inout),dimension(its:ite,kts:kte)::  & !- mudei 27/03/2024
      ! real(KIND=R8SIZE),intent(inout),dimension(kts:kte)::  & voltei para (i,k) 25/04/2024
                                                                  t
      real(KIND=R8SIZE),intent(inout),dimension(its:ite,kts:kte)::  & ! & - mudei 27/03/2024
      !real(KIND=R8SIZE),intent(inout),dimension(kts:kte)::   & - voltei para (i,k) 25/04/2024
                                                                  q, & !razao de mistura do vapor dagua
                                                                 qc, & !razao de mistura de agua de nuvem
                                                                 qi, & !razao de mistura de cristias de gelo
                                                                 qr, & !razao de mistura de chuva
                                                                 qs, & !razao de mistura de neve
                                                                 qg    !razao de mistura de graupel
      real(KIND=R8SIZE),intent(inout),dimension(its:ite)::           &
                                                               rain, &
                                                            rainncv, &
                                                                 sr
   
      real(KIND=R8SIZE),intent(inout),dimension(its:ite),optional::  &
                                                               snow, &
                                                            snowncv
      real(KIND=R8SIZE),intent(inout),dimension(its:ite),optional::  &
                                                            graupel, &
                                                         graupelncv
   
      real(KIND=R8SIZE),intent(inout),dimension(its:ite,kts:kte),   & !- Mudei 27/03/2024 voltei para (i,k) 25/04/2024
      !real(KIND=R8SIZE),intent(inout),dimension(kts:kte),   &
       optional::                                                    &
                                                         rainprod2d, &
                                                         evapprod2d
   
   !output arguments:
      character(len=*),intent(out):: errmsg
      integer,intent(out):: errflg
   
   !local variables and arrays:
      !real(KIND=R8SIZE),dimension(kts:kte,3)::             & !adicionei 27/03/2024 !voltei para (i,k,numero) 25/04/2024
                                                            !qsat, & voltei para (i,k,numero) 25/04/2024
                                                            !rh, & voltei para (i,k,numero) 25/04/2024
                                                            !falk, &
                                                            !fall, &
                                                            !qrs_tmp, &
                                                            !work1
      !real(KIND=R8SIZE),dimension(kts:kte)::             & !adicionei 27/03/2024
                                                           ! fallc, &
                                                           ! falkc, &
                                                            !workr, &
                                                            !worka, &
                                                            !!work1c!voltei para (i,k) 25/04/2024
      real(KIND=R8SIZE),dimension(its:ite,kts:kte,3)::             &
                                                                 rh, &
                                                               qsat, &
                                                             rslope, &
                                                            rslope2, &
                                                            rslope3, &
                                                            rslopeb , &
                                                            qrs_tmp, &
                                                               falk, &
                                                               fall , &
                                                              work1 !voltei para (i,k,numero) 25/04/2024
      real(KIND=R8SIZE),dimension(its:ite,kts:kte)::                &
                                                             fallc, &
                                                              falkc, &
                                                             work1c, & !voltei para (i,k) 25/04/2024
                                                             work2c , &
                                                              workr, & !voltei para (i,k) 25/04/2024
                                                              worka !voltei para (i,k) 25/04/2024
      real(KIND=R8SIZE),dimension(its:ite,kts:kte)::                & !- Mudei 27/03/2024
      !real(KIND=R8SIZE),dimension(kts:kte)::                &
                                                            den_tmp, &
                                                           delz_tmp

      real(KIND=R8SIZE),dimension(its:ite,kts:kte)::                &   !- Mudei 27/03/2024 - voltei para (i,k) 25/04/2024                                                 
      !real(KIND=R8SIZE),dimension(kts:kte)::                &
                                                              pigen, & !taxa de producao para geracao (nucleacao) de gelo a partir de vapor
                                                              pidep, & !taxa de produção para deposição/sublimação de gelo
                                                              pcond, & !Taxa de produção para condensação/evaporação de água de nuvem
                                                              prevp, & !Taxa de produção para taxa de evaporação/condensação de chuva.
                                                              psevp, & !Taxa de produção para evaporação de neve derretendo.
                                                              pgevp, & !Taxa de produção para taxa de deposição/sublimação de graupel.
                                                              psdep, & !Taxa de produção para taxa de deposição/sublimação de neve.
                                                              pgdep, & !Taxa de produção para taxa de deposição/sublimação de graupel.
                                                              praut, & !Taxa de produção para conversão automática de água de nuvem para formar chuva.
                                                              psaut, & !Taxa de produção para conversão automática de gelo de nuvem para formar neve.
                                                              pgaut, & !Taxa de produção para conversão automática de neve para formar graupel.
                                                              piacr, & !Taxa de produção para acreção de chuva por gelo de nuvem.
                                                              pracw, & !Taxa de produção para acreção de água de nuvem por chuva.
                                                              praci, & !Taxa de produção para acreção de gelo de nuvem por chuva.
                                                              pracs, & !Taxa de produção para acreção de neve por chuva.
                                                              psacw, & !Taxa de produção para acreção de água de nuvem por neve.
                                                              psaci, & !Taxa de produção para acreção de gelo de nuvem por neve.
                                                              psacr, & !Taxa de produção para acreção de chuva por neve.
                                                              pgacw, & !Taxa de produção para acreção de água de nuvem por graupel.
                                                              pgaci, & !Taxa de produção para acreção de gelo de nuvem por graupel.
                                                              pgacr, & !Taxa de produção para acreção de chuva por graupel.
                                                              pgacs, & !Taxa de produção para acreção de neve por graupel.
                                                              paacw, & !accretion of cloud water by averaged snow/graupel
                                                              psmlt, & !Taxa de produção para derretimento de neve para formar água de nuvem.
                                                              pgmlt, & !Taxa de produção induzida pela taxa de derretimento aprimorada de graupel.
                                                              pseml, & !Taxa de produção induzida pelo derretimento aprimorado de neve.
                                                              pgeml    !Taxa de produção induzida pela taxa de derretimento aprimorada de graupel.
      real(KIND=R8SIZE),dimension(its:ite,kts:kte)::                & !- Mudei 27/03/2024
                                                            denfac

      real(KIND=R8SIZE),dimension(its:ite,kts:kte)::                & !- Mudei 27/03/2024 
      !real(KIND=R8SIZE),dimension(kts:kte)::                & voltei para (i,k) 25/04/2024
                                                               qsum, &
                                                                 xl, &
                                                                cpm, &
                                                              work2, &
                                                             !denfac, & - Mudei 27/03/2024
                                                                xni, &
                                                            denqrs1, &
                                                            denqrs2, &
                                                            denqrs3, &
                                                             denqci, &
                                                             n0sfac
      real(KIND=R8SIZE),dimension(its:ite)::                        &
                                                            delqrs1, &
                                                            delqrs2, &
                                                            delqrs3, &
                                                              delqi
      real(KIND=R8SIZE),dimension(its:ite)::                        &
                                                          tstepsnow, &
                                                         tstepgraup
      integer,dimension(its:ite)::                                     &
                                                              mstep, &
                                                              numdt
      logical,dimension(its:ite)::                             flgcld
      real(KIND=R8SIZE)::                                           &
               cpmcal, xlcal, diffus,                                &
               viscos, xka, venfac, conden, diffac,                  &
               x, y, z, a, b, c, d, e,                               &
               qdt, holdrr, holdrs, holdrg, supcol, supcolt, pvt,    &
               coeres, supsat, dtcld, xmi, eacrs, satdt,             &
               qimax, diameter, xni0, roqi0,                         &
               fallsum, fallsum_qsi, fallsum_qg,                     &
               vt2i,vt2r,vt2s,vt2g,acrfac,egs,egi,                   &
               xlwork2, factor, source, value,                       &
               xlf, pfrzdtc, pfrzdtr, supice, alpha2, delta2, delta3
      real(KIND=R8SIZE):: vt2ave
      real(KIND=R8SIZE):: holdc, holdci
      integer:: i, j, k, mstepmax,                                     &
              iprt, latd, lond, loop, loops, ifsat, n, idim, kdim
   
   !Temporaries used for inlining fpvs function
      real(KIND=R8SIZE):: dldti, xb, xai, tr, xbi, xa, hvap, cvap, hsub, dldt, ttp
   
   ! variables for optimization
      real(KIND=R8SIZE),dimension(its:ite):: dvec1,tvec1
      real(KIND=R8SIZE):: temp
   
   !-----------------------------------------------------------------------------------------------------------------
   
   ! compute internal functions
   !
      cpmcal(x) = cpd*(1.-max(x,qmin))+max(x,qmin)*cpv
      xlcal(x) = xlv0-xlv1*(x-t0c)
   !----------------------------------------------------------------
   ! diffus: diffusion coefficient of the water vapor
   ! viscos: kinematic viscosity(m2s-1)
   ! Optimizatin : A**B => exp(log(A)*(B))
   !
      diffus(x,y) = 8.794e-5 * exp(log(x)*(1.81)) / y   ! 8.794e-5*x**1.81/y
      viscos(x,y) = 1.496e-6 * (x*sqrt(x)) /(x+120.)/y  ! 1.496e-6*x**1.5/(x+120.)/y
      xka(x,y) = 1.414e3*viscos(x,y)*y
      diffac(a,b,c,d,e) = d*a*a/(xka(c,d)*rv*c*c)+1./(e*diffus(c,b))
      venfac(a,b,c) = exp(log((viscos(b,c)/diffus(b,a)))*((.3333333))) &
                  /sqrt(viscos(b,c))*sqrt(sqrt(den0/c))
      conden(a,b,c,d,e) = (max(b,qmin)-c)/(1.+d*d/(rv*e)*c/(a*a))
   !
   !
      idim = ite-its+1
      kdim = kte-kts+1
   !
   !----------------------------------------------------------------
   ! paddint 0 for negative values generated by dynamics
   !
      do k = kts, kte
         do i = its, ite
            qc(i,k) = max(qc(i,k),0.0) !era qc(i,k) = max(qc(i,k),0.0) - Marilia 27/03/2024
            qr(i,k) = max(qr(i,k),0.0) !era qr(i,k) = max(qr(i,k),0.0) - Marilia 27/03/2024
            qi(i,k) = max(qi(i,k),0.0) !era qi(i,k) = max(qi(i,k),0.0) - Marilia 27/03/2024
            qs(i,k) = max(qs(i,k),0.0) !era qs(i,k) = max(qs(i,k),0.0) - Marilia 27/03/2024
            qg(i,k) = max(qg(i,k),0.0) !era qg(i,k) = max(qg(i,k),0.0) - Marilia 27/03/2024
         enddo
      enddo
   !
   !----------------------------------------------------------------
   ! latent heat for phase changes and heat capacity. neglect the
   ! changes during microphysical process calculation emanuel(1994)
   !
      do k = kts, kte
         do i = its, ite !- comentei 27/03/2024 - descomentei 25/04/2024
         cpm(i,k) = cpmcal(q(i,k)) ! era cpm(i,k) = cpmcal(q(i,k)) - voltei para (i,k) 25/04/2024
         xl(i,k) = xlcal(t(i,k)) ! era xl(i,k) = xlcal(t(i,k)) - voltei para (i,k) 25/04/2024
         enddo !- comentei 27/03/2024 - descomentei 25/04/2024
      enddo
      do k = kts, kte
         do i = its, ite  !- comentei 27/03/2024
         delz_tmp(i,k) = delz(i,k) !era delz_tmp(i,k) = delz(i, k)
         den_tmp(i,k) = den(i,k) !era den_tmp(i,k) = den(i,k)
         enddo !- comentei 27/03/2024
      enddo
   !
   !----------------------------------------------------------------
   ! initialize the surface rain, snow, graupel
   !
      do i = its, ite
         rainncv(i) = 0.
         if(present(snowncv) .and. present(snow)) snowncv(i) = 0.
         if(present(graupelncv) .and. present(graupel)) graupelncv(i) = 0.
         sr(i) = 0.
      ! new local array to catch step snow and graupel
         tstepsnow(i) = 0.
         tstepgraup(i) = 0.
      enddo
   !
   !----------------------------------------------------------------
   ! compute the minor time steps. - Calcula o numero de loops com base
   ! na divisão entre 'delt' e 'dtcldcr', garantindo que haja pelo menos
   ! um loop, mesmo que o resultado da divisão seja menor que 1.
   !
      loops = max(nint(delt/dtcldcr),1)
      dtcld = delt/loops
      if(delt.le.dtcldcr) dtcld = delt
      !
      do loop = 1,loops
   !
   !----------------------------------------------------------------
   ! initialize the large scale variables
   !
         do i = its, ite
         mstep(i) = 1
         flgcld(i) = .true.
         enddo
      !
      !  do k = kts, kte
      !    do i = its, ite
      !      denfac(i,k) = sqrt(den0/den(i,k))
      !    enddo
      !  enddo
         do k = kts, kte
         do i = its,ite
            dvec1(i) = den(i,k)
         enddo
         call vrec(tvec1,dvec1,ite-its+1)
         do i = its, ite
            tvec1(i) = tvec1(i)*den0
         enddo
         call vsqrt(dvec1,tvec1,ite-its+1)
         do i = its,ite
            denfac(i,k) = dvec1(i) !era denfac(i,k) = dvec1(i)
         enddo
         enddo
   !
   ! Inline expansion for fpvs
   !  qsat(i,k,1) = fpvs(t(i,k),0,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
   !  qsat(i,k,2) = fpvs(t(i,k),1,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
      hsub = xls
      hvap = xlv0
      cvap = cpv
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
      do k = kts, kte
        do i = its, ite !- mudei todos os (i,k,numero) -> ficou (k,numero) - voltei para (i,k,numero) 25/04/2024
          tr=ttp/t(i,k)
          qsat(i,k,1)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          qsat(i,k,1) = min(qsat(i,k,1),0.99*p(i,k))
          qsat(i,k,1) = ep2 * qsat(i,k,1) / (p(i,k) - qsat(i,k,1))
          qsat(i,k,1) = max(qsat(i,k,1),qmin)
          rh(i,k,1) = max(q(i,k) / qsat(i,k,1),qmin)
          tr=ttp/t(i,k)
          if(t(i,k).lt.ttp) then
            qsat(i,k,2)=psat*exp(log(tr)*(xai))*exp(xbi*(1.-tr))
          else
            qsat(i,k,2)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          endif
          qsat(i,k,2) = min(qsat(i,k,2),0.99*p(i,k))
          qsat(i,k,2) = ep2 * qsat(i,k,2) / (p(i,k) - qsat(i,k,2))
          qsat(i,k,2) = max(qsat(i,k,2),qmin)
          rh(i,k,2) = max(q(i,k) / qsat(i,k,2),qmin)
        enddo
      enddo
   !
   !----------------------------------------------------------------
   ! initialize the variables for microphysical physics
   !
   !
      do k = kts, kte
        do i = its, ite
          prevp(i,k) = 0.
          psdep(i,k) = 0.
          pgdep(i,k) = 0.
          praut(i,k) = 0.
          psaut(i,k) = 0.
          pgaut(i,k) = 0.
          pracw(i,k) = 0.
          praci(i,k) = 0.
          piacr(i,k) = 0.
          psaci(i,k) = 0.
          psacw(i,k) = 0.
          pracs(i,k) = 0.
          psacr(i,k) = 0.
          pgacw(i,k) = 0.
          paacw(i,k) = 0.
          pgaci(i,k) = 0.
          pgacr(i,k) = 0.
          pgacs(i,k) = 0.
          pigen(i,k) = 0.
          pidep(i,k) = 0.
          pcond(i,k) = 0.
          psmlt(i,k) = 0.
          pgmlt(i,k) = 0.
          pseml(i,k) = 0.
          pgeml(i,k) = 0.
          psevp(i,k) = 0.
          pgevp(i,k) = 0.
          falk(i,k,1) = 0.
          falk(i,k,2) = 0.
          falk(i,k,3) = 0.
          fall(i,k,1) = 0.
          fall(i,k,2) = 0.
          fall(i,k,3) = 0.
          fallc(i,k) = 0.
          falkc(i,k) = 0.
          xni(i,k) = 1.e3
        enddo
      enddo
   !-------------------------------------------------------------
   ! Ni: ice crystal number concentraiton   [HDC 5c]
   !-------------------------------------------------------------
      do k = kts, kte
        do i = its, ite !- Marilia 27/03/2024 - descomentei 24/04/2024
          temp = (den(i,k)*max(qi(i,k),qmin)) !mudei: era temp = (den(i,k)*max(qi(i,k),qmin))
          temp = sqrt(sqrt(temp*temp*temp))  
          xni(i,k) = min(max(5.38e7*temp,1.e3),1.e6) !mudei: era xni(i,k) = min(max(5.38e7*temp,1.e3),1.e6)
        enddo !- Marilia 27/03/2024
      enddo
      
   ! a sedimentação dos hidrometeoros em queda é calculada primeiro, 
   ! juntamente com o cálculo da fusão da neve e granizo dentro do 
   ! loop de sedimentação.
   ! Só vai fazer o loop de sedimentação se l_sediment for true

      if (l_sediment) then !- 24/04/2024
   !
   !----------------------------------------------------------------
   ! compute the fallout term:
   ! first, vertical terminal velosity for minor loops
   !----------------------------------------------------------------
      do k = kts, kte
        do i = its, ite
          qrs_tmp(i,k,1) = qr(i,k) !mudei: era qrs_tmp(i,k,1) = qr(i,k)
          qrs_tmp(i,k,2) = qs(i,k) !mudei: era qrs_tmp(i,k,2) = qs(i,k)
          qrs_tmp(i,k,3) = qg(i,k) !mudei: era qrs_tmp(i,k,3) = qg(i,k)
        enddo
      enddo
      call slope_wsm6(qrs_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2,rslope3, &
                      work1,its,ite,kts,kte)
   !
      do k = kte, kts, -1
        do i = its, ite !- Marilia 27/03/2024 - tirei todos os (i,k) - ficou só (k)- descomentei 24/04/2024
          workr(i,k) = work1(i,k,1)
          qsum(i,k) = max( (qs(i,k)+qg(i,k)), 1.E-15)
          if( qsum(i,k) .gt. 1.e-15 ) then
            worka(i,k) = (work1(i,k,2)*qs(i,k) + work1(i,k,3)*qg(i,k)) &
                       / qsum(i,k)
          else
            worka(i,k) = 0.
          endif
          denqrs1(i,k) = den(i,k)*qr(i,k)
          denqrs2(i,k) = den(i,k)*qs(i,k)
          denqrs3(i,k) = den(i,k)*qg(i,k)
          if(qr(i,k).le.0.0) workr(i,k) = 0.0
        enddo
      enddo
      call nislfv_rain_plm(idim,kdim,den_tmp,denfac,t,delz_tmp,workr,denqrs1,  &
                           delqrs1,dtcld,1,1)
      call nislfv_rain_plm6(idim,kdim,den_tmp,denfac,t,delz_tmp,worka,         &
                            denqrs2,denqrs3,delqrs2,delqrs3,dtcld,1,1)
      do k = kts, kte
        do i = its, ite !- mudei (i,k) para (k) - Marilia 27/03/2024 - descomentei 24/04/2024 e voltei ao que era
          qr(i,k) = max(denqrs1(i,k)/den(i,k),0.)
          qs(i,k) = max(denqrs2(i,k)/den(i,k),0.)
          qg(i,k) = max(denqrs3(i,k)/den(i,k),0.)
          fall(i,k,1) = denqrs1(i,k)*workr(i,k)/delz(i,k)
          fall(i,k,2) = denqrs2(i,k)*worka(i,k)/delz(i,k)
          fall(i,k,3) = denqrs3(i,k)*worka(i,k)/delz(i,k)
        enddo
      enddo
      do i = its, ite
        !fall(i,1,1) = delqrs1(i)/delz(i,1)/dtcld
        !fall(i,1,2) = delqrs2(i)/delz(i,1)/dtcld
        !fall(i,1,3) = delqrs3(i)/delz(i,1)/dtcld - Era assim antes. Deixei como está abaixo - Marilia 15/04/2024
        fall(i,1,1) = delqrs1(i)/delz(i,1)/dtcld !comentado para testes 23/04/2024
        fall(i,1,2) = delqrs2(i)/delz(i,1)/dtcld
        fall(i,1,3) = delqrs3(i)/delz(i,1)/dtcld
      enddo
      do k = kts, kte
        do i = its, ite !mudei (i,k,numero) para (k,numero) e (i,k) para (k) - Marilia 27/03/2024 
          qrs_tmp(i,k,1) = qr(i,k)
          qrs_tmp(i,k,2) = qs(i,k)
          qrs_tmp(i,k,3) = qg(i,k)
        enddo
      enddo
   end if !Marilia - termina o loop l_sediment - comentei 24/04/2024
      call slope_wsm6(qrs_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2,rslope3, &
                      work1,its,ite,kts,kte)
   !end if !Marilia - termina o loop l_sediment
      do k = kte, kts, -1
        do i = its, ite !- mudei (i,k) para (k) - Marilia 27/03/2024 - voltei para o que era antes 25/04/2024
          supcol = t0c-t(i,k) !
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(t(i,k).gt.t0c) then
   !---------------------------------------------------------------
   ! psmlt: melting of snow [HL A33] [RH83 A25]
   !       (T>T0: S->R)
   !---------------------------------------------------------------
            xlf = xlf0
            work2(i,k) = venfac(p(i,k),t(i,k),den(i,k)) !- mudei (i,k) para (k) - Marilia 27/03/2024
            if(qs(i,k).gt.0.) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psmlt(i,k) = xka(t(i,k),den(i,k))/xlf*(t0c-t(i,k))*pi/2.         &
                         *n0sfac(i,k)*(precs1*rslope2(i,k,2)                   &
                         +precs2*work2(i,k)*coeres)/den(i,k)
              psmlt(i,k) = min(max(psmlt(i,k)*dtcld/mstep(i),                  &
                         -qs(i,k)/mstep(i)),0.)
              qs(i,k) = qs(i,k) + psmlt(i,k)
              qr(i,k) = qr(i,k) - psmlt(i,k)
              t(i,k) = t(i,k) + xlf/cpm(i,k)*psmlt(i,k)
            endif
   !---------------------------------------------------------------
   ! pgmlt: melting of graupel [HL A23]  [LFO 47]
   !       (T>T0: G->R)
   !---------------------------------------------------------------
            if(qg(i,k).gt.0.) then !- mudei (i,k) para (k) - Marilia 27/03/2024
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgmlt(i,k) = xka(t(i,k),den(i,k))/xlf                            &
                         *(t0c-t(i,k))*(precg1*rslope2(i,k,3)                  &
                         +precg2*work2(i,k)*coeres)/den(i,k)
              pgmlt(i,k) = min(max(pgmlt(i,k)*dtcld/mstep(i),                  &
                        -qg(i,k)/mstep(i)),0.)
              qg(i,k) = qg(i,k) + pgmlt(i,k)
              qr(i,k) = qr(i,k) - pgmlt(i,k)
              t(i,k) = t(i,k) + xlf/cpm(i,k)*pgmlt(i,k)
            endif
          endif
        enddo
      enddo
      !endif  !Marilia - termina o loop
   !---------------------------------------------------------------
   ! Vice [ms-1] : fallout of ice crystal [HDC 5a] !queda de cristais de gelo
   !---------------------------------------------------------------
      do k = kte, kts, -1
        do i = its, ite !- mudei (i,k) para (k) - Marilia 27/03/2024
          if(qi(i,k).le.0.) then
            work1c(i,k) = 0.
          else
            xmi = den(i,k)*qi(i,k)/xni(i,k)
            diameter  = max(min(dicon * sqrt(xmi),dimax), 1.e-25)
            work1c(i,k) = 1.49e4*exp(log(diameter)*(1.31))
          endif
        enddo
      enddo
   !
   !  forward semi-laglangian scheme (JH), PCM (piecewise constant),  (linear)
   !
      do k = kte, kts, -1
        do i = its, ite !- mudei (i,k) para (k) - Marilia 27/03/2024, com excecao do den - voltei para (i,k)
          denqci(i,k) = den(i,k)*qi(i,k)
        enddo
      enddo
      call nislfv_rain_plm(idim,kdim,den_tmp,denfac,t,delz_tmp,work1c,denqci,  &
                           delqi,dtcld,1,0)
      do k = kts, kte
        do i = its, ite !- mudei (i,k) para (k) - Marilia 27/03/2024, com excecao do den - voltei para (i,k)
          qi(i,k) = max(denqci(i,k)/den(i,k),0.)
        enddo
      enddo
      do i = its, ite
      !do k = kts, kte
        !fallc(i,1) = delqi(i)/delz(i,1)/dtcld - mudei para como está na linha abaixo - Marili 15/04/2024
        !fallc(k) = delqi(i)/delz(i)/dtcld
        fallc(i,1) = delqi(i)/delz(i,1)/dtcld
      enddo
   !
   !----------------------------------------------------------------
   !      rain (unit is mm/sec;kgm-2s-1: /1000*delt ===> m)==> mm for wrf
   !
   
      do i = its, ite !- mudei (i,kts,numero) para (kts,numero) e (i,kts) para (kts) - Marilia 27/03/2024
        fallsum = fall(i,kts,1)+fall(i,kts,2)+fall(i,kts,3)+fallc(i,kts) !era fallc(i,kts), ficou fallc(kts) - Marilia 15/04/2024
        fallsum_qsi = fall(i,kts,2)+fallc(i,kts) !voltei ao que era inicialmente fall, fallc, delz - 24/04/2024
        fallsum_qg = fall(i,kts,3)
        if(fallsum.gt.0.) then
          rainncv(i) = fallsum*delz(i,kts)/denr*dtcld*1000. + rainncv(i)
          rain(i) = fallsum*delz(i,kts)/denr*dtcld*1000. + rain(i)
        endif
        if(fallsum_qsi.gt.0.) then
          tstepsnow(i) = fallsum_qsi*delz(i,kts)/denr*dtcld*1000.              &
                       + tstepsnow(i)
          if(present(snowncv) .and. present(snow)) then
            snowncv(i) = fallsum_qsi*delz(i,kts)/denr*dtcld*1000.              &
                       + snowncv(i)
            snow(i) = fallsum_qsi*delz(i,kts)/denr*dtcld*1000. + snow(i)
          endif
        endif
        if(fallsum_qg.gt.0.) then
          tstepgraup(i)  = fallsum_qg*delz(i,kts)/denr*dtcld*1000.             &
                         + tstepgraup(i)
          if(present (graupelncv) .and. present (graupel)) then
            graupelncv(i) = fallsum_qg*delz(i,kts)/denr*dtcld*1000.            &
                          + graupelncv(i)
            graupel(i) = fallsum_qg*delz(i,kts)/denr*dtcld*1000. + graupel(i)
          endif
        endif
        if(present (snowncv)) then
           if(fallsum.gt.0.)sr(i)=(snowncv(i) + graupelncv(i))/(rainncv(i)+1.e-12)
        else
           if(fallsum.gt.0.)sr(i)=(tstepsnow(i) + tstepgraup(i))/(rainncv(i)+1.e-12)
        endif
      enddo
   !
   !---------------------------------------------------------------
   ! pimlt: instantaneous melting of cloud ice [HL A47] [RH83 A28]
   !       (T>T0: I->C)
   !---------------------------------------------------------------
      do k = kts, kte 
        do i = its, ite !- mudei (i,k) para (k) - Marilia 27/03/2024, com excecao do denc- voltei para (i,k)
          supcol = t0c-t(i,k)
          xlf = xls-xl(i,k)
          if(supcol.lt.0.) xlf = xlf0
          if(supcol.lt.0.and.qi(i,k).gt.0.) then
            qc(i,k) = qc(i,k) + qi(i,k)
            t(i,k) = t(i,k) - xlf/cpm(i,k)*qi(i,k)
            qi(i,k) = 0.
          endif
   !---------------------------------------------------------------
   ! pihmf: homogeneous freezing of cloud water below -40c [HL A45]
   !        (T<-40C: C->I)
   !---------------------------------------------------------------
          if(supcol.gt.40..and.qc(i,k).gt.0.) then !- mudei (i,k) para (k) - Marilia 27/03/2024, com excecao do den - voltei para (i,k)
            qi(i,k) = qi(i,k) + qc(i,k)
            t(i,k) = t(i,k) + xlf/cpm(i,k)*qc(i,k)
            qc(i,k) = 0.
          endif
   !---------------------------------------------------------------
   ! pihtf: heterogeneous freezing of cloud water [HL A44]
   !        (T0>T>-40C: C->I)
   !---------------------------------------------------------------
          if(supcol.gt.0..and.qc(i,k).gt.qmin) then !- mudei (i,k) para (k) - Marilia 27/03/2024, com excecao do den - voltei para (i,k)
   !        pfrzdtc = min(pfrz1*(exp(pfrz2*supcol)-1.)                         &
   !                * den(i,k)/denr/xncr*qc(i,k)**2*dtcld,qc(i,k))
            supcolt=min(supcol,50.)
            pfrzdtc = min(pfrz1*(exp(pfrz2*supcolt)-1.)                        &
                    * den(i,k)/denr/xncr*qc(i,k)*qc(i,k)*dtcld,qc(i,k))
            qi(i,k) = qi(i,k) + pfrzdtc
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtc
            qc(i,k) = qc(i,k)-pfrzdtc
          endif
   !---------------------------------------------------------------
   ! pgfrz: freezing of rain water [HL A20] [LFO 45]
   !        (T<T0, R->G)
   !---------------------------------------------------------------
          if(supcol.gt.0..and.qr(i,k).gt.0.) then !- mudei (i,k) para (k) - Marilia 27/03/2024, com excecao do den - voltei para (i,k)
   !        pfrzdtr = min(20.*pi**2*pfrz1*n0r*denr/den(i,k)                    &
   !                * (exp(pfrz2*supcol)-1.)*rslope3(i,k,1)**2                 &
   !                * rslope(i,k,1)*dtcld,qr(i,k))
            temp = rslope3(i,k,1)
            temp = temp*temp*rslope(i,k,1)
            supcolt=min(supcol,50.)
            pfrzdtr = min(20.*(pi*pi)*pfrz1*n0r*denr/den(i,k)                  &
                    *(exp(pfrz2*supcolt)-1.)*temp*dtcld,                       &
                      qr(i,k))
            qg(i,k) = qg(i,k) + pfrzdtr
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtr
            qr(i,k) = qr(i,k)-pfrzdtr
          endif
        enddo
      enddo
   !
   !
   !----------------------------------------------------------------
   ! update the slope parameters for microphysics computation
   !
      do k = kts, kte 
        !do i = its, ite !- mudei (i,k, numero) para (k, numero) e (i,k) para (k) - - voltei para (i,k)
          qrs_tmp(i,k,1) = qr(i,k)
          qrs_tmp(i,k,2) = qs(i,k)
          qrs_tmp(i,k,3) = qg(i,k)
        !enddo
      enddo
      call slope_wsm6(qrs_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2,rslope3, &
                      work1,its,ite,kts,kte)
   !------------------------------------------------------------------
   ! work1: the thermodynamic term in the denominator associated with
   !        heat conduction and vapor diffusion (ry88, y93, h85)
   ! work1: termo termodinâmico no denominador associado à condução de 
   !        calor e difusão de vapor (ry88, y93, h85)
   !        
   ! work2: parameter associated with the ventilation effects(y93)
   ! work2: parâmetro associado aos efeitos de ventilação (y93)
   !
      do k = kts, kte
        do i = its, ite !- mudei (i,k, numero) para (k, numero) e (i,k) para (k) - Marilia 27/03/2024, com excecao do den
          work1(i,k,1) = diffac(xl(i,k),p(i,k),t(i,k),den(i,k),qsat(i,k,1))
          work1(i,k,2) = diffac(xls,p(i,k),t(i,k),den(i,k),qsat(i,k,2))
          work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
        enddo
      enddo
   !
   !===============================================================
   !
   ! warm rain processes
   !
   ! - follows the processes in RH83 and LFO except for autoconcersion
   !
   !===============================================================
   !
      !if (iiwarm) then !marilia
      do k = kts, kte
        do i = its, ite !- mudei (i,k, numero) para (k, numero) e (i,k) para (k) - Marilia 27/03/2024, com excecao do den - voltei para (i,k) 25/04/2024
          supsat = max(q(i,k),qmin)-qsat(i,k,1)
          satdt = supsat/dtcld
   !---------------------------------------------------------------
   ! praut: auto conversion rate from cloud to rain [HDC 16]
   !        (C->R)
   !---------------------------------------------------------------
          if(qc(i,k).gt.qc0) then !- voltei para (i,k) 25/04/2024
            praut(i,k) = qck1*qc(i,k)**(7./3.)
            praut(i,k) = min(praut(i,k),qc(i,k)/dtcld)
          endif
   !---------------------------------------------------------------
   ! pracw: accretion of cloud water by rain [HL A40] [LFO 51]
   !        (C->R)
   !---------------------------------------------------------------
          if(qr(i,k).gt.qcrmin.and.qc(i,k).gt.qmin) then !- voltei para (i,k) 25/04/2024
            pracw(i,k) = min(pacrr*rslope3(i,k,1)*rslopeb(i,k,1)               &
                       * qc(i,k)*denfac(i,k),qc(i,k)/dtcld)
          endif
   !---------------------------------------------------------------
   ! prevp: evaporation/condensation rate of rain [HDC 14]
   !        (V->R or R->V)
   !---------------------------------------------------------------
          if(qr(i,k).gt.0.) then !- voltei para (i,k) 25/04/2024
            coeres = rslope2(i,k,1)*sqrt(rslope(i,k,1)*rslopeb(i,k,1))
            prevp(i,k) = (rh(i,k,1)-1.)*(precr1*rslope2(i,k,1)                 &
                       + precr2*work2(i,k)*coeres)/work1(i,k,1)
            if(prevp(i,k).lt.0.) then
              prevp(i,k) = max(prevp(i,k),-qr(i,k)/dtcld)
              prevp(i,k) = max(prevp(i,k),satdt/2)
            else
              prevp(i,k) = min(prevp(i,k),satdt/2)
            endif
          endif
        enddo
      enddo

      !else  !marilia
   !   
      if(.not. iiwarm) then
   !===============================================================
   !
   ! cold rain processes
   !
   ! - follows the revised ice microphysics processes in HDC
   ! - the processes same as in RH83 and RH84  and LFO behave
   !   following ice crystal hapits defined in HDC, inclduing
   !   intercept parameter for snow (n0s), ice crystal number
   !   concentration (ni), ice nuclei number concentration
   !   (n0i), ice diameter (d)
   !
   !===============================================================
   !
      do k = kts, kte
        do i = its, ite !- comentei e modifiquei todos os (i,k) para (k) 15/04/2024 !- voltei para (i,k) 25/04/2024
          supcol = t0c-t(i,k) 
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          supsat = max(q(i,k),qmin)-qsat(i,k,2) !era (q(i,k))
          satdt = supsat/dtcld
          ifsat = 0
   !-------------------------------------------------------------
   ! Ni: ice crystal number concentraiton   [HDC 5c]
   !-------------------------------------------------------------
   !      xni(i,k) = min(max(5.38e7*(den(i,k)                                  &
   !               * max(qi(i,k),qmin))**0.75,1.e3),1.e6) !- voltei para (i,k) 25/04/2024
          temp = (den(i,k)*max(qi(i,k),qmin))
          temp = sqrt(sqrt(temp*temp*temp))
          xni(i,k) = min(max(5.38e7*temp,1.e3),1.e6)
          eacrs = exp(0.07*(-supcol))
   !
          xmi = den(i,k)*qi(i,k)/xni(i,k) !- voltei para (i,k) 25/04/2024
          diameter  = min(dicon * sqrt(xmi),dimax)
          vt2i = 1.49e4*diameter**1.31
          vt2r=pvtr*rslopeb(i,k,1)*denfac(i,k)
          vt2s=pvts*rslopeb(i,k,2)*denfac(i,k)
          vt2g=pvtg*rslopeb(i,k,3)*denfac(i,k)
          qsum(i,k) = max( (qs(i,k)+qg(i,k)), 1.E-15)
          if(qsum(i,k) .gt. 1.e-15) then
            vt2ave=(vt2s*qs(i,k)+vt2g*qg(i,k))/(qsum(i,k))
          else
            vt2ave=0.
          endif
          if(supcol.gt.0.and.qi(i,k).gt.qmin) then
            if(qr(i,k).gt.qcrmin) then
   !-------------------------------------------------------------
   ! praci: accretion of cloud ice by rain [HL A15] [LFO 25]
   !          (T<T0: I->R)
   !-------------------------------------------------------------
              acrfac = 2.*rslope3(i,k,1)+2.*diameter*rslope2(i,k,1)            & !- voltei para (i,k) 25/04/2024
                     + diameter**2*rslope(i,k,1)
              praci(i,k) = pi*qi(i,k)*n0r*abs(vt2r-vt2i)*acrfac/4.
   ! reduce collection efficiency (suggested by B. Wilt)
              praci(i,k) = praci(i,k)*min(max(0.0,qr(i,k)/qi(i,k)),1.)**2
              praci(i,k) = min(praci(i,k),qi(i,k)/dtcld)
   !-------------------------------------------------------------
   ! piacr: accretion of rain by cloud ice [HL A19] [LFO 26]
   !        (T<T0: R->S or R->G)
   !-------------------------------------------------------------
              piacr(i,k) = pi**2*avtr*n0r*denr*xni(i,k)*denfac(i,k)            & !- voltei para (i,k) 25/04/2024
                         * g6pbr*rslope3(i,k,1)*rslope3(i,k,1)                 &
                         * rslopeb(i,k,1)/24./den(i,k)
   ! reduce collection efficiency (suggested by B. Wilt)
              piacr(i,k) = piacr(i,k)*min(max(0.0,qi(i,k)/qr(i,k)),1.)**2
              piacr(i,k) = min(piacr(i,k),qr(i,k)/dtcld)
            endif
   !-------------------------------------------------------------
   ! psaci: accretion of cloud ice by snow [HDC 10]
   !        (T<T0: I->S)
   !-------------------------------------------------------------
            if(qs(i,k).gt.qcrmin) then
              acrfac = 2.*rslope3(i,k,2)+2.*diameter*rslope2(i,k,2)            &
                     + diameter**2*rslope(i,k,2)
              psaci(i,k) = pi*qi(i,k)*eacrs*n0s*n0sfac(i,k)                    &
                         * abs(vt2ave-vt2i)*acrfac/4.
              psaci(i,k) = min(psaci(i,k),qi(i,k)/dtcld)
            endif
   !-------------------------------------------------------------
   ! pgaci: accretion of cloud ice by graupel [HL A17] [LFO 41]
   !        (T<T0: I->G)
   !-------------------------------------------------------------
            if(qg(i,k).gt.qcrmin) then !- voltei para (i,k) 25/04/2024
              egi = exp(0.07*(-supcol))
              acrfac = 2.*rslope3(i,k,3)+2.*diameter*rslope2(i,k,3)            &
                     + diameter**2*rslope(i,k,3)
              pgaci(i,k) = pi*egi*qi(i,k)*n0g*abs(vt2ave-vt2i)*acrfac/4.
              pgaci(i,k) = min(pgaci(i,k),qi(i,k)/dtcld)
            endif
          endif
   !-------------------------------------------------------------
   ! psacw: accretion of cloud water by snow  [HL A7] [LFO 24]
   !        (T<T0: C->S, and T>=T0: C->R)
   !-------------------------------------------------------------
          if(qs(i,k).gt.qcrmin.and.qc(i,k).gt.qmin) then !- voltei para (i,k) 25/04/2024
            psacw(i,k) = min(pacrc*n0sfac(i,k)*rslope3(i,k,2)*rslopeb(i,k,2)   &
   ! reduce collection efficiency (suggested by B. Wilt)
                       * min(max(0.0,qs(i,k)/qc(i,k)),1.)**2                   &
                       * qc(i,k)*denfac(i,k),qc(i,k)/dtcld)
          endif
   !-------------------------------------------------------------
   ! pgacw: accretion of cloud water by graupel [HL A6] [LFO 40]
   !        (T<T0: C->G, and T>=T0: C->R)
   !-------------------------------------------------------------
          if(qg(i,k).gt.qcrmin.and.qc(i,k).gt.qmin) then !- voltei para (i,k) 25/04/2024
            pgacw(i,k) = min(pacrg*rslope3(i,k,3)*rslopeb(i,k,3)               &
   ! reduce collection efficiency (suggested by B. Wilt)
                       * min(max(0.0,qg(i,k)/qc(i,k)),1.)**2                   &
                       * qc(i,k)*denfac(i,k),qc(i,k)/dtcld)
          endif
   !-------------------------------------------------------------
   ! paacw: accretion of cloud water by averaged snow/graupel
   !        (T<T0: C->G or S, and T>=T0: C->R)
   !-------------------------------------------------------------
          if(qsum(i,k) .gt. 1.e-15) then !- voltei para (i,k) 25/04/2024
            paacw(i,k) = (qs(i,k)*psacw(i,k)+qg(i,k)*pgacw(i,k))               &
                       /(qsum(i,k))
          endif
   !-------------------------------------------------------------
   ! pracs: accretion of snow by rain [HL A11] [LFO 27]
   !         (T<T0: S->G)
   !-------------------------------------------------------------
          if(qs(i,k).gt.qcrmin.and.qr(i,k).gt.qcrmin) then !- voltei para (i,k) 25/04/2024
            if(supcol.gt.0) then
              acrfac = 5.*rslope3(i,k,2)*rslope3(i,k,2)*rslope(i,k,1)          &
                     + 2.*rslope3(i,k,2)*rslope2(i,k,2)*rslope2(i,k,1)         &
                     + .5*rslope2(i,k,2)*rslope2(i,k,2)*rslope3(i,k,1)
              pracs(i,k) = pi**2*n0r*n0s*n0sfac(i,k)*abs(vt2r-vt2ave)          &
                         * (dens/den(i,k))*acrfac
   ! reduce collection efficiency (suggested by B. Wilt)
              pracs(i,k) = pracs(i,k)*min(max(0.0,qr(i,k)/qs(i,k)),1.)**2
              pracs(i,k) = min(pracs(i,k),qs(i,k)/dtcld)
            endif
   !-------------------------------------------------------------
   ! psacr: accretion of rain by snow [HL A10] [LFO 28]
   !         (T<T0:R->S or R->G) (T>=T0: enhance melting of snow)
   !-------------------------------------------------------------
            acrfac = 5.*rslope3(i,k,1)*rslope3(i,k,1)*rslope(i,k,2)            & !- voltei para (i,k) 25/04/2024
                   + 2.*rslope3(i,k,1)*rslope2(i,k,1)*rslope2(i,k,2)           &
                   +.5*rslope2(i,k,1)*rslope2(i,k,1)*rslope3(i,k,2)
            psacr(i,k) = pi**2*n0r*n0s*n0sfac(i,k)*abs(vt2ave-vt2r)            &
                       * (denr/den(i,k))*acrfac
   ! reduce collection efficiency (suggested by B. Wilt)
            psacr(i,k) = psacr(i,k)*min(max(0.0,qs(i,k)/qr(i,k)),1.)**2
            psacr(i,k) = min(psacr(i,k),qr(i,k)/dtcld)
          endif
   !-------------------------------------------------------------
   ! pgacr: accretion of rain by graupel [HL A12] [LFO 42]
   !         (T<T0: R->G) (T>=T0: enhance melting of graupel)
   !-------------------------------------------------------------
          if(qg(i,k).gt.qcrmin.and.qr(i,k).gt.qcrmin) then !- voltei para (i,k) 25/04/2024
            acrfac = 5.*rslope3(i,k,1)*rslope3(i,k,1)*rslope(i,k,3)            &
                   + 2.*rslope3(i,k,1)*rslope2(i,k,1)*rslope2(i,k,3)           &
                   + .5*rslope2(i,k,1)*rslope2(i,k,1)*rslope3(i,k,3)
            pgacr(i,k) = pi**2*n0r*n0g*abs(vt2ave-vt2r)*(denr/den(i,k))        &
                       * acrfac
   ! reduce collection efficiency (suggested by B. Wilt)
            pgacr(i,k) = pgacr(i,k)*min(max(0.0,qg(i,k)/qr(i,k)),1.)**2
            pgacr(i,k) = min(pgacr(i,k),qr(i,k)/dtcld)
          endif
   !
   !-------------------------------------------------------------
   ! pgacs: accretion of snow by graupel [HL A13] [LFO 29]
   !        (S->G): This process is eliminated in V3.0 with the
   !        new combined snow/graupel fall speeds
   !-------------------------------------------------------------
          if(qg(i,k).gt.qcrmin.and.qs(i,k).gt.qcrmin) then !- voltei para (i,k) 25/04/2024
            pgacs(i,k) = 0.
          endif
          if(supcol.le.0) then
            xlf = xlf0
   !-------------------------------------------------------------
   ! pseml: enhanced melting of snow by accretion of water [HL A34]
   !        (T>=T0: S->R)
   !-------------------------------------------------------------
            if(qs(i,k).gt.0.)                                                   &!- voltei para (i,k) 25/04/2024
              pseml(i,k) = min(max(cliq*supcol*(paacw(i,k)+psacr(i,k))         &
                         / xlf,-qs(i,k)/dtcld),0.)
   !-------------------------------------------------------------
   ! pgeml: enhanced melting of graupel by accretion of water [HL A24] [RH84 A21-A22]
   !        (T>=T0: G->R)
   !-------------------------------------------------------------
            if(qg(i,k).gt.0.)                                                  & !- voltei para (i,k) 25/04/2024
              pgeml(i,k) = min(max(cliq*supcol*(paacw(i,k)+pgacr(i,k))         &
                         / xlf,-qg(i,k)/dtcld),0.)
          endif
          if(supcol.gt.0) then
   !-------------------------------------------------------------
   ! pidep: deposition/Sublimation rate of ice [HDC 9]
   !       (T<T0: V->I or I->V)
   !-------------------------------------------------------------
            if(qi(i,k).gt.0.and.ifsat.ne.1) then !- voltei para (i,k) 25/04/2024
               pidep(i,k) = 4.*diameter*xni(i,k)*(rh(i,k,2)-1.)/work1(i,k,2)
               supice = satdt-prevp(i,k)
               if(pidep(i,k).lt.0.) then
                pidep(i,k) = max(max(pidep(i,k),satdt/2),supice)
                pidep(i,k) = max(pidep(i,k),-qi(i,k)/dtcld)
               else
                pidep(i,k) = min(min(pidep(i,k),satdt/2),supice)
               endif
               if(abs(prevp(i,k)+pidep(i,k)).ge.abs(satdt)) ifsat = 1
               endif
   !-------------------------------------------------------------
   ! psdep: deposition/sublimation rate of snow [HDC 14]
   !        (T<T0: V->S or S->V)
   !-------------------------------------------------------------
            if(qs(i,k).gt.0..and.ifsat.ne.1) then !- voltei para (i,k) 25/04/2024
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psdep(i,k) = (rh(i,k,2)-1.)*n0sfac(i,k)*(precs1*rslope2(i,k,2)   &
                         + precs2*work2(i,k)*coeres)/work1(i,k,2)
              supice = satdt-prevp(i,k)-pidep(i,k)
              if(psdep(i,k).lt.0.) then
                psdep(i,k) = max(psdep(i,k),-qs(i,k)/dtcld)
                psdep(i,k) = max(max(psdep(i,k),satdt/2),supice)
              else
                psdep(i,k) = min(min(psdep(i,k),satdt/2),supice)
              endif
              if(abs(prevp(i,k)+pidep(i,k)+psdep(i,k)).ge.abs(satdt)) &
                ifsat = 1
            endif
   !-------------------------------------------------------------
   ! pgdep: deposition/sublimation rate of graupel [HL A21] [LFO 46]
   !        (T<T0: V->G or G->V)
   !-------------------------------------------------------------
            if(qg(i,k).gt.0..and.ifsat.ne.1) then !- voltei para (i,k) 25/04/2024
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgdep(i,k) = (rh(i,k,2)-1.)*(precg1*rslope2(i,k,3)       &
                         + precg2*work2(i,k)*coeres)/work1(i,k,2)
              supice = satdt-prevp(i,k)-pidep(i,k)-psdep(i,k)
              if(pgdep(i,k).lt.0.) then
                pgdep(i,k) = max(pgdep(i,k),-qg(i,k)/dtcld)
                pgdep(i,k) = max(max(pgdep(i,k),satdt/2),supice)
              else
                pgdep(i,k) = min(min(pgdep(i,k),satdt/2),supice)
              endif
              if(abs(prevp(i,k)+pidep(i,k)+psdep(i,k)+pgdep(i,k)).ge.  &
                abs(satdt)) ifsat = 1
            endif
   !-------------------------------------------------------------
   ! pigen: generation(nucleation) of ice from vapor [HL 50] [HDC 7-8]
   !       (T<T0: V->I)
   !-------------------------------------------------------------
            if(supsat.gt.0.and.ifsat.ne.1) then !- voltei para (i,k) 25/04/2024
              supice = satdt-prevp(i,k)-pidep(i,k)-psdep(i,k)-pgdep(i,k)
              xni0 = 1.e3*exp(0.1*supcol)
              roqi0 = 4.92e-11*xni0**1.33
              pigen(i,k) = max(0.,(roqi0/den(i,k)-max(qi(i,k),0.))/dtcld)
              pigen(i,k) = min(min(pigen(i,k),satdt),supice)
            endif
   !
   !-------------------------------------------------------------
   ! psaut: conversion(aggregation) of ice to snow [HDC 12]
   !        (T<T0: I->S)
   !-------------------------------------------------------------
            if(qi(i,k).gt.0.) then !- voltei para (i,k) 25/04/2024
              qimax = roqimax/den(i,k)
              psaut(i,k) = max(0.,(qi(i,k)-qimax)/dtcld)
            endif
   !
   !-------------------------------------------------------------
   ! pgaut: conversion(aggregation) of snow to graupel [HL A4] [LFO 37]
   !        (T<T0: S->G)
   !-------------------------------------------------------------
            if(qs(i,k).gt.0.) then !- voltei para (i,k) 25/04/2024
              alpha2 = 1.e-3*exp(0.09*(-supcol))
              pgaut(i,k) = min(max(0.,alpha2*(qs(i,k)-qs0)),qs(i,k)/dtcld)
            endif
          endif
   !
   !-------------------------------------------------------------
   ! psevp: evaporation of melting snow [HL A35] [RH83 A27]
   !       (T>=T0: S->V)
   !-------------------------------------------------------------
          if(supcol.lt.0.) then !- voltei para (i,k) 25/04/2024
            if(qs(i,k).gt.0..and.rh(i,k,1).lt.1.) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psevp(i,k) = (rh(i,k,1)-1.)*n0sfac(i,k)*(precs1                  &
                         * rslope2(i,k,2)+precs2*work2(i,k)                    &
                         * coeres)/work1(i,k,1)
              psevp(i,k) = min(max(psevp(i,k),-qs(i,k)/dtcld),0.)
            endif
   !-------------------------------------------------------------
   ! pgevp: evaporation of melting graupel [HL A25] [RH84 A19]
   !       (T>=T0: G->V)
   !-------------------------------------------------------------
            !- mudei: (i,k) para (k) e (i,k,numero) para (k,numero) - Marilia 27/03/2024 !- voltei para (i,k) 25/04/2024
            if(qg(i,k).gt.0..and.rh(i,k,1).lt.1.) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgevp(i,k) = (rh(i,k,1)-1.)*(precg1*rslope2(i,k,3)               &
                         + precg2*work2(i,k)*coeres)/work1(i,k,1)
              pgevp(i,k) = min(max(pgevp(i,k),-qg(i,k)/dtcld),0.)
            endif
          endif
        enddo
      enddo
      endif  !marilia
   !
   !
   !----------------------------------------------------------------
   !     check mass conservation of generation terms and feedback to the
   !     large scale
   !
      !- mudei: (i,k) para (k) e (i,k,numero) para (k,numero) - Marilia 27/03/2024 !- voltei para (i,k) 25/04/2024
         do k = kts, kte
         do i = its, ite
      !
            delta2=0.
            delta3=0.
            if(qr(i,k).lt.1.e-4.and.qs(i,k).lt.1.e-4) delta2=1.
            if(qr(i,k).lt.1.e-4) delta3=1.
            if(t(i,k).le.t0c) then
      !
      ! cloud water
      !
               !- mudei: (i,k) para (k) e (i,k,numero) para (k,numero) - Marilia 27/03/2024 !- voltei para (i,k) 25/04/2024
               value = max(qmin,qc(i,k))
               source = (praut(i,k)+pracw(i,k)+paacw(i,k)+paacw(i,k))*dtcld
               if (source.gt.value) then
               factor = value/source
               praut(i,k) = praut(i,k)*factor
               pracw(i,k) = pracw(i,k)*factor
               paacw(i,k) = paacw(i,k)*factor
               endif
      !
      ! cloud ice
      !
               !- mudei: (i,k) para (k) e (i,k,numero) para (k,numero) - Marilia 27/03/2024 !- voltei para (i,k) 25/04/2024
               value = max(qmin,qi(i,k))
               source = (psaut(i,k)-pigen(i,k)-pidep(i,k)+praci(i,k)+psaci(i,k)   &
                     + pgaci(i,k))*dtcld
               if (source.gt.value) then
               factor = value/source
               psaut(i,k) = psaut(i,k)*factor
               pigen(i,k) = pigen(i,k)*factor
               pidep(i,k) = pidep(i,k)*factor
               praci(i,k) = praci(i,k)*factor
               psaci(i,k) = psaci(i,k)*factor
               pgaci(i,k) = pgaci(i,k)*factor
               endif
      !
      ! rain
      !
               !- mudei: (i,k) para (k) e (i,k,numero) para (k,numero) - Marilia 27/03/2024 !- voltei para (i,k) 25/04/2024
               value = max(qmin,qr(i,k))
               source = (-praut(i,k)-prevp(i,k)-pracw(i,k)+piacr(i,k)+psacr(i,k)  &
                     + pgacr(i,k))*dtcld
               if (source.gt.value) then
               factor = value/source
               praut(i,k) = praut(i,k)*factor
               prevp(i,k) = prevp(i,k)*factor
               pracw(i,k) = pracw(i,k)*factor
               piacr(i,k) = piacr(i,k)*factor
               psacr(i,k) = psacr(i,k)*factor
               pgacr(i,k) = pgacr(i,k)*factor
               endif
      !
      ! snow
      !
               !- mudei: (i,k) para (k) e (i,k,numero) para (k,numero) - Marilia 27/03/2024 !- voltei para (i,k) 25/04/2024
               value = max(qmin,qs(i,k))
               source = -(psdep(i,k)+psaut(i,k)-pgaut(i,k)+paacw(i,k)+piacr(i,k)  &
                     * delta3+praci(i,k)*delta3-pracs(i,k)*(1.-delta2)           &
                     + psacr(i,k)*delta2+psaci(i,k)-pgacs(i,k) )*dtcld
               if (source.gt.value) then
               factor = value/source
               psdep(i,k) = psdep(i,k)*factor
               psaut(i,k) = psaut(i,k)*factor
               pgaut(i,k) = pgaut(i,k)*factor
               paacw(i,k) = paacw(i,k)*factor
               piacr(i,k) = piacr(i,k)*factor
               praci(i,k) = praci(i,k)*factor
               psaci(i,k) = psaci(i,k)*factor
               pracs(i,k) = pracs(i,k)*factor
               psacr(i,k) = psacr(i,k)*factor
               pgacs(i,k) = pgacs(i,k)*factor
               endif
      !
      ! graupel
      !
               !- mudei: (i,k) para (k) e (i,k,numero) para (k,numero) - Marilia 27/03/2024 !- voltei para (i,k) 25/04/2024
               value = max(qmin,qg(i,k))
               source = -(pgdep(i,k)+pgaut(i,k)                                   &
                     + piacr(i,k)*(1.-delta3)+praci(i,k)*(1.-delta3)             &
                     + psacr(i,k)*(1.-delta2)+pracs(i,k)*(1.-delta2)             &
                     + pgaci(i,k)+paacw(i,k)+pgacr(i,k)+pgacs(i,k))*dtcld
               if (source.gt.value) then
               factor = value/source
               pgdep(i,k) = pgdep(i,k)*factor
               pgaut(i,k) = pgaut(i,k)*factor
               piacr(i,k) = piacr(i,k)*factor
               praci(i,k) = praci(i,k)*factor
               psacr(i,k) = psacr(i,k)*factor
               pracs(i,k) = pracs(i,k)*factor
               paacw(i,k) = paacw(i,k)*factor
               pgaci(i,k) = pgaci(i,k)*factor
               pgacr(i,k) = pgacr(i,k)*factor
               pgacs(i,k) = pgacs(i,k)*factor
               endif
      !
                  !- mudei: (i,k) para (k) e (i,k,numero) para (k,numero) - Marilia 27/03/2024
               work2(i,k)=-(prevp(i,k)+psdep(i,k)+pgdep(i,k)+pigen(i,k)+pidep(i,k))
      ! update
               q(i,k) = q(i,k)+work2(i,k)*dtcld
               qc(i,k) = max(qc(i,k)-(praut(i,k)+pracw(i,k)                       &
                     + paacw(i,k)+paacw(i,k))*dtcld,0.)
               qr(i,k) = max(qr(i,k)+(praut(i,k)+pracw(i,k)                       &
                     + prevp(i,k)-piacr(i,k)-pgacr(i,k)                         &
                     - psacr(i,k))*dtcld,0.)
               qi(i,k) = max(qi(i,k)-(psaut(i,k)+praci(i,k)                       &
                     + psaci(i,k)+pgaci(i,k)-pigen(i,k)-pidep(i,k))             &
                     * dtcld,0.)
               qs(i,k) = max(qs(i,k)+(psdep(i,k)+psaut(i,k)+paacw(i,k)            &
                     - pgaut(i,k)+piacr(i,k)*delta3                             &
                     + praci(i,k)*delta3+psaci(i,k)-pgacs(i,k)                  &
                     - pracs(i,k)*(1.-delta2)+psacr(i,k)*delta2)                &
                     * dtcld,0.)
               qg(i,k) = max(qg(i,k)+(pgdep(i,k)+pgaut(i,k)                       &
                     + piacr(i,k)*(1.-delta3)                                   &
                     + praci(i,k)*(1.-delta3)+psacr(i,k)*(1.-delta2)            &
                     + pracs(i,k)*(1.-delta2)+pgaci(i,k)+paacw(i,k)             &
                     + pgacr(i,k)+pgacs(i,k))*dtcld,0.)
               xlf = xls-xl(i,k)
               xlwork2 = -xls*(psdep(i,k)+pgdep(i,k)+pidep(i,k)+pigen(i,k))       &
                        -xl(i,k)*prevp(i,k)-xlf*(piacr(i,k)+paacw(i,k)           &
                        +paacw(i,k)+pgacr(i,k)+psacr(i,k))
               t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
            else
      !
      ! cloud water
      !
                  !- mudei: (i,k) para (k) e (i,k,numero) para (k,numero) - Marilia 27/03/2024 !- voltei para (i,k) 25/04/2024
               value = max(qmin,qc(i,k))
               source=(praut(i,k)+pracw(i,k)+paacw(i,k)+paacw(i,k))*dtcld
               if (source.gt.value) then
               factor = value/source
               praut(i,k) = praut(i,k)*factor
               pracw(i,k) = pracw(i,k)*factor
               paacw(i,k) = paacw(i,k)*factor
               endif
      !
      ! rain
      !
                  !- mudei: (i,k) para (k) e (i,k,numero) para (k,numero) - Marilia 27/03/2024 !- voltei para (i,k) 25/04/2024
               value = max(qmin,qr(i,k))
               source = (-paacw(i,k)-praut(i,k)+pseml(i,k)+pgeml(i,k)-pracw(i,k)  &
                        -paacw(i,k)-prevp(i,k))*dtcld
               if (source.gt.value) then
               factor = value/source
               praut(i,k) = praut(i,k)*factor
               prevp(i,k) = prevp(i,k)*factor
               pracw(i,k) = pracw(i,k)*factor
               paacw(i,k) = paacw(i,k)*factor
               pseml(i,k) = pseml(i,k)*factor
               pgeml(i,k) = pgeml(i,k)*factor
               endif
      !
      ! snow
      !
                  !- mudei: (i,k) para (k) e (i,k,numero) para (k,numero) - Marilia 27/03/2024!- voltei para (i,k) 25/04/2024
               value = max(qcrmin,qs(i,k))
               source=(pgacs(i,k)-pseml(i,k)-psevp(i,k))*dtcld
               if (source.gt.value) then
               factor = value/source
               pgacs(i,k) = pgacs(i,k)*factor
               psevp(i,k) = psevp(i,k)*factor
               pseml(i,k) = pseml(i,k)*factor
               endif
      !
      ! graupel
      !
               !- mudei: (i,k) para (k) e (i,k,numero) para (k,numero) - Marilia 27/03/2024 !- voltei para (i,k) 25/04/2024
               value = max(qcrmin,qg(i,k))
               source=-(pgacs(i,k)+pgevp(i,k)+pgeml(i,k))*dtcld
               if (source.gt.value) then
               factor = value/source
               pgacs(i,k) = pgacs(i,k)*factor
               pgevp(i,k) = pgevp(i,k)*factor
               pgeml(i,k) = pgeml(i,k)*factor
               endif
      !
               work2(i,k)=-(prevp(i,k)+psevp(i,k)+pgevp(i,k)) !- mudei: (i,k) para (k) e (i,k,numero) para (k,numero) - Marilia 27/03/2024
      ! update !- voltei para (i,k) 25/04/2024
               q(i,k) = q(i,k)+work2(i,k)*dtcld
               qc(i,k) = max(qc(i,k)-(praut(i,k)+pracw(i,k)                       &
                     + paacw(i,k)+paacw(i,k))*dtcld,0.)
               qr(i,k) = max(qr(i,k)+(praut(i,k)+pracw(i,k)                       &
                     + prevp(i,k)+paacw(i,k)+paacw(i,k)-pseml(i,k)              &
                     - pgeml(i,k))*dtcld,0.)
               qs(i,k) = max(qs(i,k)+(psevp(i,k)-pgacs(i,k)                       &
                     + pseml(i,k))*dtcld,0.)
               qg(i,k) = max(qg(i,k)+(pgacs(i,k)+pgevp(i,k)                       &
                     + pgeml(i,k))*dtcld,0.)
               xlf = xls-xl(i,k)
               xlwork2 = -xl(i,k)*(prevp(i,k)+psevp(i,k)+pgevp(i,k))              &
                        -xlf*(pseml(i,k)+pgeml(i,k))
               t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
            endif

            if (nx == 1) then

               if (.not. iiwarm) then

                  call save_dg(k, prevp(i,k), 'prevp', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, praut(i,k), 'praut', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, pracw(i,k), 'pracw', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, pcond(i,k), 'pcond', i_dgtime, units='/kg&
                        &/s',dim='z') 
                  call save_dg(k, pgmlt(i,k), 'pgmlt', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, pseml(i,k), 'pseml', i_dgtime, units='/kg&
                        &/s',dim='z')  

               end if

                  call save_dg(k, psdep(i,k), 'psdep', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, pgdep(i,k), 'pgdep', i_dgtime, units='/kg&
                        &/s',dim='z') 
                  call save_dg(k, psaut(i,k), 'psaut', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, pgaut(i,k), 'pgaut', i_dgtime, units='/kg&
                        &/s',dim='z')                     
                  call save_dg(k, praci(i,k), 'praci', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, piacr(i,k), 'piacr', i_dgtime, units='/kg&
                        &/s',dim='z') 
                  call save_dg(k, psaci(i,k), 'psaci', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, psacw(i,k), 'psacw', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, pracs(i,k), 'pracs', i_dgtime, units='/kg&
                        &/s',dim='z')    
                  call save_dg(k, psacr(i,k), 'psacr', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, pgacw(i,k), 'pgacw', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, paacw(i,k), 'paacw', i_dgtime, units='/kg&
                        &/s',dim='z') 
                  call save_dg(k, pgaci(i,k), 'pgaci', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, pgacr(i,k), 'pgacr', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, pgacs(i,k), 'pgacs', i_dgtime, units='/kg&
                        &/s',dim='z')                          
                  call save_dg(k, pigen(i,k), 'pigen', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, pidep(i,k), 'pidep', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, psmlt(i,k), 'psmlt', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, pgeml(i,k), 'pgeml', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, psevp(i,k), 'psevp', i_dgtime, units='/kg&
                        &/s',dim='z')
                  call save_dg(k, pgevp(i,k), 'pgevp', i_dgtime, units='/kg&
                        &/s',dim='z')    
            end if

         enddo
         enddo
      !
      ! Inline expansion for fpvs
      !  qsat(i,k,1) = fpvs(t(i,k),0,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
      !  qsat(i,k,2) = fpvs(t(i,k),1,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
         hsub = xls
         hvap = xlv0
         cvap = cpv
         ttp=t0c+0.01
         dldt=cvap-cliq
         xa=-dldt/rv
         xb=xa+hvap/(rv*ttp)
         dldti=cvap-cice
         xai=-dldti/rv
         xbi=xai+hsub/(rv*ttp)
         do k = kts, kte
         do i = its, ite !- mudei: (i,k) para (k) e (i,k,numero) para (k,numero) - Marilia 27/03/2024 !- voltei para (i,k) 25/04/2024
            tr=ttp/t(i,k)
            qsat(i,k,1)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
            qsat(i,k,1) = min(qsat(i,k,1),0.99*p(i,k))
            qsat(i,k,1) = ep2 * qsat(i,k,1) / (p(i,k) - qsat(i,k,1))
            qsat(i,k,1) = max(qsat(i,k,1),qmin)
            tr=ttp/t(i,k)
            if(t(i,k).lt.ttp) then
               qsat(i,k,2)=psat*exp(log(tr)*(xai))*exp(xbi*(1.-tr))
            else
               qsat(i,k,2)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
            endif
            qsat(i,k,2) = min(qsat(i,k,2),0.99*p(i,k))
            qsat(i,k,2) = ep2 * qsat(i,k,2) / (p(i,k) - qsat(i,k,2))
            qsat(i,k,2) = max(qsat(i,k,2),qmin)
         enddo
         enddo
      !
      !----------------------------------------------------------------
      ! pcond: condensational/evaporational rate of cloud water [HL A46] [RH83 A6]
      ! if there exists additional water vapor condensated/if
      ! evaporation of cloud water is not enough to remove subsaturation
      !
         do k = kts, kte
          do i = its, ite !- mudei: (i,k) para (k) e (i,k,numero) para (k,numero) - Marilia 27/03/2024 
            work1(i,k,1) = conden(t(i,k),q(i,k),qsat(i,k,1),xl(i,k),cpm(i,k))
            work2(i,k) = qc(i,k)+work1(i,k,1)
            pcond(i,k) = min(max(work1(i,k,1)/dtcld,0.),max(q(i,k),0.)/dtcld)
            if(qc(i,k).gt.0..and.work1(i,k,1).lt.0.)                          &
               pcond(i,k) = max(work1(i,k,1),-qc(i,k))/dtcld
               q(i,k) = q(i,k)-pcond(i,k)*dtcld
               qc(i,k) = max(qc(i,k)+pcond(i,k)*dtcld,0.)
               t(i,k) = t(i,k)+pcond(i,k)*xl(i,k)/cpm(i,k)*dtcld
          enddo
         enddo
      !
      !
      !----------------------------------------------------------------
      ! padding for small values
      !
         do k = kts, kte
            do i = its, ite !- mudei: (i,k) para (k) - Marilia 27/03/2024 !- voltei para (i,k) 25/04/2024
               if(qc(i,k).le.qmin) qc(i,k) = 0.0
               if(qi(i,k).le.qmin) qi(i,k) = 0.0
            enddo
         enddo
      enddo                  ! big loops
   
      if(present(rainprod2d) .and. present(evapprod2d)) then
         do k = kts, kte
          do i = its,ite !- mudei: (i,k) para (k) !- Marilia 27/03/2024 !- voltei para (i,k) 25/04/2024
            rainprod2d(i,k) = praut(i,k)+pracw(i,k)+praci(i,k)+psaci(i,k)+pgaci(i,k) &
                           + psacw(i,k)+pgacw(i,k)+paacw(i,k)+psaut(i,k)
            evapprod2d(i,k) = -(prevp(i,k)+psevp(i,k)+pgevp(i,k)  &
                           + psdep(i,k)+pgdep(i,k))
          enddo
         enddo
      endif


   !
   !----------------------------------------------------------------
   ! CCPP checks:
   !
   
      errmsg = 'mp_wsm6_run OK'
      errflg = 0
   
      end subroutine mp_wsm6_run
   
   !=================================================================================================================
      real(KIND=R8SIZE) function rgmma(x)
   !=================================================================================================================
   !rgmma function:  use infinite product form
   
      real(KIND=R8SIZE),intent(in):: x
      
      integer:: i
      real(KIND=R8SIZE),parameter:: euler=0.577215664901532
      real(KIND=R8SIZE):: y
      
      !-----------------------------------------------------------------------------------------------------------------
      
      if(x.eq.1.)then
         rgmma=0.
      else
         rgmma=x*exp(euler*x)
         do i = 1,10000
            y = float(i)
            rgmma=rgmma*(1.000+x/y)*exp(-x/y)
         enddo
         rgmma=1./rgmma
      endif
      
      end function rgmma
   
   !=================================================================================================================
      real(KIND=R8SIZE) function fpvs(t,ice,rd,rv,cvap,cliq,cice,hvap,hsub,psat,t0c)
   !=================================================================================================================
   
      integer,intent(in):: ice
      real(KIND=R8SIZE),intent(in):: cice,cliq,cvap,hsub,hvap,psat,rd,rv,t0c
      real(KIND=R8SIZE),intent(in):: t
      
      real(KIND=R8SIZE):: tr,ttp,dldt,dldti,xa,xb,xai,xbi
      
      !-----------------------------------------------------------------------------------------------------------------
      
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
      tr=ttp/t
      if(t.lt.ttp.and.ice.eq.1) then
         fpvs=psat*(tr**xai)*exp(xbi*(1.-tr))
      else
         fpvs=psat*(tr**xa)*exp(xb*(1.-tr))
      endif
      
      end function fpvs
   
   !=================================================================================================================
      subroutine slope_wsm6(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,vt,its,ite,kts,kte)
   !=================================================================================================================
   
   !--- input arguments:
      integer:: its,ite,kts,kte
      
      real(KIND=R8SIZE),intent(in),dimension(its:ite,kts:kte):: den,denfac,t
      real(KIND=R8SIZE),intent(in),dimension(its:ite,kts:kte,3):: qrs
      
      !--- inout arguments:
      real(KIND=R8SIZE),intent(inout),dimension(its:ite,kts:kte,3):: &
         rslope,rslopeb,rslope2,rslope3,vt
      
      !--- local variables and arrays:
      integer:: i,k
      
      real(KIND=R8SIZE),parameter:: t0c = 273.15
      real(KIND=R8SIZE):: lamdar,lamdas,lamdag,x,y,z,supcol
      real(KIND=R8SIZE),dimension(its:ite,kts:kte):: n0sfac
      
      !-----------------------------------------------------------------------------------------------------------------
      
      !size distributions: (x=mixing ratio, y=air density):
      !valid for mixing ratio > 1.e-9 kg/kg.
      lamdar(x,y)=   sqrt(sqrt(pidn0r/(x*y)))      ! (pidn0r/(x*y))**.25
      lamdas(x,y,z)= sqrt(sqrt(pidn0s*z/(x*y)))    ! (pidn0s*z/(x*y))**.25
      lamdag(x,y)=   sqrt(sqrt(pidn0g/(x*y)))      ! (pidn0g/(x*y))**.25
      
      do k = kts, kte
         do i = its, ite
         supcol = t0c-t(i,k)
      !---------------------------------------------------------------
      ! n0s: Intercept parameter for snow [m-4] [HDC 6]
      !---------------------------------------------------------------
         n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
         if(qrs(i,k,1).le.qcrmin)then
            rslope(i,k,1) = rslopermax
            rslopeb(i,k,1) = rsloperbmax
            rslope2(i,k,1) = rsloper2max
            rslope3(i,k,1) = rsloper3max
         else
            rslope(i,k,1) = 1./lamdar(qrs(i,k,1),den(i,k))
            rslopeb(i,k,1) = rslope(i,k,1)**bvtr
            rslope2(i,k,1) = rslope(i,k,1)*rslope(i,k,1)
            rslope3(i,k,1) = rslope2(i,k,1)*rslope(i,k,1)
         endif
         if(qrs(i,k,2).le.qcrmin)then
            rslope(i,k,2) = rslopesmax
            rslopeb(i,k,2) = rslopesbmax
            rslope2(i,k,2) = rslopes2max
            rslope3(i,k,2) = rslopes3max
         else
            rslope(i,k,2) = 1./lamdas(qrs(i,k,2),den(i,k),n0sfac(i,k))
            rslopeb(i,k,2) = rslope(i,k,2)**bvts
            rslope2(i,k,2) = rslope(i,k,2)*rslope(i,k,2)
            rslope3(i,k,2) = rslope2(i,k,2)*rslope(i,k,2)
         endif
         if(qrs(i,k,3).le.qcrmin)then
            rslope(i,k,3) = rslopegmax
            rslopeb(i,k,3) = rslopegbmax
            rslope2(i,k,3) = rslopeg2max
            rslope3(i,k,3) = rslopeg3max
         else
            rslope(i,k,3) = 1./lamdag(qrs(i,k,3),den(i,k))
            rslopeb(i,k,3) = rslope(i,k,3)**bvtg
            rslope2(i,k,3) = rslope(i,k,3)*rslope(i,k,3)
            rslope3(i,k,3) = rslope2(i,k,3)*rslope(i,k,3)
         endif
         vt(i,k,1) = pvtr*rslopeb(i,k,1)*denfac(i,k)
         vt(i,k,2) = pvts*rslopeb(i,k,2)*denfac(i,k)
         vt(i,k,3) = pvtg*rslopeb(i,k,3)*denfac(i,k)
         if(qrs(i,k,1).le.0.0) vt(i,k,1) = 0.0
         if(qrs(i,k,2).le.0.0) vt(i,k,2) = 0.0
         if(qrs(i,k,3).le.0.0) vt(i,k,3) = 0.0
         enddo
      enddo
      
      end subroutine slope_wsm6
   
   !=================================================================================================================
      subroutine slope_rain(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,vt,its,ite,kts,kte)
   !=================================================================================================================
   
   !--- input arguments:
      integer:: its,ite,kts,kte
      
      real(KIND=R8SIZE),intent(in),dimension(its:ite,kts:kte):: qrs,den,denfac,t
      
      !--- inout arguments:
      real(KIND=R8SIZE),intent(inout),dimension(its:ite,kts:kte):: &
         rslope,rslopeb,rslope2,rslope3,vt
      
      !--- local variables and arrays:
      integer:: i,k
      
      real(KIND=R8SIZE),parameter:: t0c = 273.15
      real(KIND=R8SIZE):: lamdar,x,y
      real(KIND=R8SIZE),dimension(its:ite,kts:kte):: n0sfac
      
      !-----------------------------------------------------------------------------------------------------------------
      
      !size distributions: (x=mixing ratio, y=air density):
      !valid for mixing ratio > 1.e-9 kg/kg.
      lamdar(x,y)=   sqrt(sqrt(pidn0r/(x*y)))      ! (pidn0r/(x*y))**.25
      
      do k = kts, kte
         do i = its, ite
         if(qrs(i,k).le.qcrmin)then
            rslope(i,k) = rslopermax
            rslopeb(i,k) = rsloperbmax
            rslope2(i,k) = rsloper2max
            rslope3(i,k) = rsloper3max
         else
            rslope(i,k) = 1./lamdar(qrs(i,k),den(i,k))
            rslopeb(i,k) = rslope(i,k)**bvtr
            rslope2(i,k) = rslope(i,k)*rslope(i,k)
            rslope3(i,k) = rslope2(i,k)*rslope(i,k)
         endif
         vt(i,k) = pvtr*rslopeb(i,k)*denfac(i,k)
         if(qrs(i,k).le.0.0) vt(i,k) = 0.0
         enddo
      enddo
      
      end subroutine slope_rain
   
   !=================================================================================================================
      subroutine slope_snow(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,vt,its,ite,kts,kte)
   !=================================================================================================================
   
   !--- input arguments:
      integer:: its,ite,kts,kte
   
      real(KIND=R8SIZE),intent(in),dimension(its:ite,kts:kte):: qrs,den,denfac,t
   
   !--- inout arguments:
      real(KIND=R8SIZE),intent(inout),dimension(its:ite,kts:kte):: &
       rslope,rslopeb,rslope2,rslope3,vt
   
   !--- local variables and arrays:
      integer:: i,k
   
      real(KIND=R8SIZE),parameter:: t0c = 273.15
      real(KIND=R8SIZE):: lamdas,x,y,z,supcol
      real(KIND=R8SIZE),dimension(its:ite,kts:kte):: n0sfac
   
   !-----------------------------------------------------------------------------------------------------------------
   
   !size distributions: (x=mixing ratio, y=air density):
   !valid for mixing ratio > 1.e-9 kg/kg.
      lamdas(x,y,z)= sqrt(sqrt(pidn0s*z/(x*y)))    ! (pidn0s*z/(x*y))**.25
   !
      do k = kts, kte
         do i = its, ite
         supcol = t0c-t(i,k)
      !---------------------------------------------------------------
      ! n0s: Intercept parameter for snow [m-4] [HDC 6]
      !---------------------------------------------------------------
         n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
         if(qrs(i,k).le.qcrmin)then
            rslope(i,k) = rslopesmax
            rslopeb(i,k) = rslopesbmax
            rslope2(i,k) = rslopes2max
            rslope3(i,k) = rslopes3max
         else
            rslope(i,k) = 1./lamdas(qrs(i,k),den(i,k),n0sfac(i,k))
            rslopeb(i,k) = rslope(i,k)**bvts
            rslope2(i,k) = rslope(i,k)*rslope(i,k)
            rslope3(i,k) = rslope2(i,k)*rslope(i,k)
         endif
         vt(i,k) = pvts*rslopeb(i,k)*denfac(i,k)
         if(qrs(i,k).le.0.0) vt(i,k) = 0.0
         enddo
      enddo
      
      end subroutine slope_snow
   
   !=================================================================================================================
      subroutine slope_graup(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,vt,its,ite,kts,kte)
   !=================================================================================================================
   
   !--- input arguments:
      integer:: its,ite,kts,kte
      
      real(KIND=R8SIZE),intent(in),dimension(its:ite,kts:kte):: qrs,den,denfac,t
      
      !--- inout arguments:
      real(KIND=R8SIZE),intent(inout),dimension(its:ite,kts:kte):: &
         rslope,rslopeb,rslope2,rslope3,vt
      
      !--- local variables and arrays:
      integer:: i,k
      
      real(KIND=R8SIZE),parameter:: t0c = 273.15
      real(KIND=R8SIZE):: lamdag,x,y
      real(KIND=R8SIZE),dimension(its:ite,kts:kte):: n0sfac
   
   !-----------------------------------------------------------------------------------------------------------------
   
   !size distributions: (x=mixing ratio, y=air density):
   !valid for mixing ratio > 1.e-9 kg/kg.
      lamdag(x,y)=   sqrt(sqrt(pidn0g/(x*y)))      ! (pidn0g/(x*y))**.25
      
      do k = kts, kte
         do i = its, ite
      !---------------------------------------------------------------
      ! n0s: Intercept parameter for snow [m-4] [HDC 6]
      !---------------------------------------------------------------
         if(qrs(i,k).le.qcrmin)then
            rslope(i,k) = rslopegmax
            rslopeb(i,k) = rslopegbmax
            rslope2(i,k) = rslopeg2max
            rslope3(i,k) = rslopeg3max
         else
            rslope(i,k) = 1./lamdag(qrs(i,k),den(i,k))
            rslopeb(i,k) = rslope(i,k)**bvtg
            rslope2(i,k) = rslope(i,k)*rslope(i,k)
            rslope3(i,k) = rslope2(i,k)*rslope(i,k)
         endif
         vt(i,k) = pvtg*rslopeb(i,k)*denfac(i,k)
         if(qrs(i,k).le.0.0) vt(i,k) = 0.0
         enddo
      enddo
      
      end subroutine slope_graup
   
   !=================================================================================================================
      subroutine nislfv_rain_plm(im,km,denl,denfacl,tkl,dzl,wwl,rql,precip,dt,id,iter)
   !=================================================================================================================
   !
   ! for non-iteration semi-Lagrangain forward advection for cloud
   ! with mass conservation and positive definite advection
   ! 2nd order interpolation with monotonic piecewise linear method
   ! this routine is under assumption of decfl < 1 for semi_Lagrangian
   !
   ! dzl    depth of model layer in meter
   ! wwl    terminal velocity at model layer m/s
   ! rql    cloud density*mixing ration
   ! precip precipitation
   ! dt     time step
   ! id     kind of precip: 0 test case; 1 raindrop
   ! iter   how many time to guess mean terminal velocity: 0 pure forward.
   !        0 : use departure wind for advection
   !        1 : use mean wind for advection
   !        > 1 : use mean wind after iter-1 iterations
   !
   ! author: hann-ming henry juang <henry.juang@noaa.gov>
   !         implemented by song-you hong
   !
   
   !--- input arguments:
      integer,intent(in):: im,km,id,iter
      
      real(KIND=R8SIZE),intent(in):: dt
      real(KIND=R8SIZE),intent(in),dimension(im,km):: dzl,denl,denfacl,tkl
      
      !--- inout arguments:
      real(KIND=R8SIZE),intent(inout),dimension(im):: precip
      real(KIND=R8SIZE),intent(inout),dimension(im,km):: rql,wwl
      
      !---- local variables and arrays:
      integer:: i,k,n,m,kk,kb,kt
      real(KIND=R8SIZE):: tl,tl2,qql,dql,qqd
      real(KIND=R8SIZE):: th,th2,qqh,dqh
      real(KIND=R8SIZE):: zsum,qsum,dim,dip,c1,con1,fa1,fa2
      real(KIND=R8SIZE):: allold,allnew,zz,dzamin,cflmax,decfl
      real(KIND=R8SIZE),dimension(km):: dz,ww,qq,wd,wa,was
      real(KIND=R8SIZE),dimension(km):: den,denfac,tk
      real(KIND=R8SIZE),dimension(km):: qn,qr,tmp,tmp1,tmp2,tmp3
      real(KIND=R8SIZE),dimension(km+1):: wi,zi,za
      real(KIND=R8SIZE),dimension(km+1):: dza,qa,qmi,qpi
      
      !-----------------------------------------------------------------------------------------------------------------
      
      precip(:) = 0.0
      
      i_loop: do i=1,im
         dz(:) = dzl(i,:)
         qq(:) = rql(i,:)
         ww(:) = wwl(i,:)
         den(:) = denl(i,:)
         denfac(:) = denfacl(i,:)
         tk(:) = tkl(i,:)
      ! skip for no precipitation for all layers
         allold = 0.0
         do k=1,km
            allold = allold + qq(k)
         enddo
         if(allold.le.0.0) then
            cycle i_loop
         endif
      !
      ! compute interface values
         zi(1)=0.0
         do k=1,km
            zi(k+1) = zi(k)+dz(k)
         enddo
      !
      ! save departure wind
         wd(:) = ww(:)
         n=1
      100  continue
      ! plm is 2nd order, we can use 2nd order wi or 3rd order wi
      ! 2nd order interpolation to get wi
         wi(1) = ww(1)
         wi(km+1) = ww(km)
         do k=2,km
            wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
         enddo
      ! 3rd order interpolation to get wi
         fa1 = 9./16.
         fa2 = 1./16.
         wi(1) = ww(1)
         wi(2) = 0.5*(ww(2)+ww(1))
         do k=3,km-1
            wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
         enddo
         wi(km) = 0.5*(ww(km)+ww(km-1))
         wi(km+1) = ww(km)
      !
      ! terminate of top of raingroup
         do k=2,km
            if( ww(k).eq.0.0 ) wi(k)=ww(k-1)
         enddo
      !
      ! diffusivity of wi
         con1 = 0.05
         do k=km,1,-1
            decfl = (wi(k+1)-wi(k))*dt/dz(k)
            if( decfl .gt. con1 ) then
               wi(k) = wi(k+1) - con1*dz(k)/dt
            endif
         enddo
      ! compute arrival point
         do k=1,km+1
            za(k) = zi(k) - wi(k)*dt
         enddo
      !
         do k=1,km
            dza(k) = za(k+1)-za(k)
         enddo
         dza(km+1) = zi(km+1) - za(km+1)
      !
      ! computer deformation at arrival point
         do k=1,km
            qa(k) = qq(k)*dz(k)/dza(k)
            qr(k) = qa(k)/den(k)
         enddo
         qa(km+1) = 0.0
      !     call maxmin(km,1,qa,' arrival points ')
      !
      ! compute arrival terminal velocity, and estimate mean terminal velocity
      ! then back to use mean terminal velocity
         if( n.le.iter ) then
            call slope_rain(qr,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa, &
                            1,1,1,km)
            if( n.ge.2 ) wa(1:km)=0.5*(wa(1:km)+was(1:km))
            do k=1,km
      !#ifdef DEBUG
      !        print*,' slope_wsm3 ',qr(k)*1000.,den(k),denfac(k),tk(k),tmp(k),tmp1(k),tmp2(k),ww(k),wa(k)
      !#endif
      ! mean wind is average of departure and new arrival winds
               ww(k) = 0.5* ( wd(k)+wa(k) )
            enddo
            was(:) = wa(:)
            n=n+1
            go to 100
         endif
      !
      ! estimate values at arrival cell interface with monotone
         do k=2,km
            dip=(qa(k+1)-qa(k))/(dza(k+1)+dza(k))
            dim=(qa(k)-qa(k-1))/(dza(k-1)+dza(k))
            if( dip*dim.le.0.0 ) then
               qmi(k)=qa(k)
               qpi(k)=qa(k)
            else
               qpi(k)=qa(k)+0.5*(dip+dim)*dza(k)
               qmi(k)=2.0*qa(k)-qpi(k)
               if( qpi(k).lt.0.0 .or. qmi(k).lt.0.0 ) then
                  qpi(k) = qa(k)
                  qmi(k) = qa(k)
               endif
            endif
         enddo
         qpi(1)=qa(1)
         qmi(1)=qa(1)
         qmi(km+1)=qa(km+1)
         qpi(km+1)=qa(km+1)
      !
      ! interpolation to regular point
         qn = 0.0
         kb=1
         kt=1
         intp : do k=1,km
            kb=max(kb-1,1)
            kt=max(kt-1,1)
      ! find kb and kt
            if( zi(k).ge.za(km+1) ) then
               exit intp
            else
               find_kb : do kk=kb,km
                  if( zi(k).le.za(kk+1) ) then
                     kb = kk
                     exit find_kb
                  else
                     cycle find_kb
                  endif
               enddo find_kb
               find_kt : do kk=kt,km
                  if( zi(k+1).le.za(kk) ) then
                     kt = kk
                     exit find_kt
                  else
                     cycle find_kt
                  endif
               enddo find_kt
               kt = kt - 1
      ! compute q with piecewise constant method
               if( kt.eq.kb ) then
                  tl=(zi(k)-za(kb))/dza(kb)
                  th=(zi(k+1)-za(kb))/dza(kb)
                  tl2=tl*tl
                  th2=th*th
                  qqd=0.5*(qpi(kb)-qmi(kb))
                  qqh=qqd*th2+qmi(kb)*th
                  qql=qqd*tl2+qmi(kb)*tl
                  qn(k) = (qqh-qql)/(th-tl)
               else if( kt.gt.kb ) then
                  tl=(zi(k)-za(kb))/dza(kb)
                  tl2=tl*tl
                  qqd=0.5*(qpi(kb)-qmi(kb))
                  qql=qqd*tl2+qmi(kb)*tl
                  dql = qa(kb)-qql
                  zsum  = (1.-tl)*dza(kb)
                  qsum  = dql*dza(kb)
                  if( kt-kb.gt.1 ) then
                     do m=kb+1,kt-1
                        zsum = zsum + dza(m)
                        qsum = qsum + qa(m) * dza(m)
                     enddo
                  endif
                  th=(zi(k+1)-za(kt))/dza(kt)
                  th2=th*th
                  qqd=0.5*(qpi(kt)-qmi(kt))
                  dqh=qqd*th2+qmi(kt)*th
                  zsum  = zsum + th*dza(kt)
                  qsum  = qsum + dqh*dza(kt)
                  qn(k) = qsum/zsum
               endif
               cycle intp
            endif
      !
         enddo intp
      !
      ! rain out
         sum_precip: do k=1,km
            if( za(k).lt.0.0 .and. za(k+1).lt.0.0 ) then
               precip(i) = precip(i) + qa(k)*dza(k)
               cycle sum_precip
            else if ( za(k).lt.0.0 .and. za(k+1).ge.0.0 ) then
               precip(i) = precip(i) + qa(k)*(0.0-za(k))
               exit sum_precip
            endif
            exit sum_precip
         enddo sum_precip
      !
      ! replace the new values
         rql(i,:) = qn(:)
      enddo i_loop
      
      end subroutine nislfv_rain_plm
   
   !=================================================================================================================
      subroutine nislfv_rain_plm6(im,km,denl,denfacl,tkl,dzl,wwl,rql,rql2,precip1,precip2,dt,id,iter)
   !=================================================================================================================
   !
   ! for non-iteration semi-Lagrangain forward advection for cloud
   ! with mass conservation and positive definite advection
   ! 2nd order interpolation with monotonic piecewise linear method
   ! this routine is under assumption of decfl < 1 for semi_Lagrangian
   !
   ! dzl    depth of model layer in meter
   ! wwl    terminal velocity at model layer m/s
   ! rql    cloud density*mixing ration
   ! precip precipitation
   ! dt     time step
   ! id     kind of precip: 0 test case; 1 raindrop
   ! iter   how many time to guess mean terminal velocity: 0 pure forward.
   !        0 : use departure wind for advection
   !        1 : use mean wind for advection
   !        > 1 : use mean wind after iter-1 iterations
   !
   ! author: hann-ming henry juang <henry.juang@noaa.gov>
   !         implemented by song-you hong
   !
   
   !--- input arguments:
      integer,intent(in):: im,km,id,iter
      
      real(KIND=R8SIZE),intent(in):: dt
      real(KIND=R8SIZE),intent(in),dimension(im,km):: dzl,denl,denfacl,tkl
      
      !--- inout arguments:
      real(KIND=R8SIZE),intent(inout),dimension(im):: precip1,precip2
      real(KIND=R8SIZE),intent(inout),dimension(im,km):: rql,rql2,wwl
      
      !---- local variables and arrays:
      integer:: i,ist,k,n,m,kk,kb,kt
      real(KIND=R8SIZE):: tl,tl2,qql,dql,qqd
      real(KIND=R8SIZE):: th,th2,qqh,dqh
      real(KIND=R8SIZE):: zsum,qsum,dim,dip,c1,con1,fa1,fa2
      real(KIND=R8SIZE):: allold,allnew,zz,dzamin,cflmax,decfl
      real(KIND=R8SIZE),dimension(km):: dz,ww,qq,qq2,wd,wa,wa2,was
      real(KIND=R8SIZE),dimension(km):: den,denfac,tk
      real(KIND=R8SIZE),dimension(km):: qn,qr,qr2,tmp,tmp1,tmp2,tmp3
      real(KIND=R8SIZE),dimension(km+1):: wi,zi,za
      real(KIND=R8SIZE),dimension(km+1):: dza,qa,qa2,qmi,qpi
      real(KIND=R8SIZE),dimension(im):: precip
      
      !-----------------------------------------------------------------------------------------------------------------
      
      precip(:) = 0.0
      precip1(:) = 0.0
      precip2(:) = 0.0
      
      i_loop: do i=1,im
         dz(:) = dzl(i,:)
         qq(:) = rql(i,:)
         qq2(:) = rql2(i,:)
         ww(:) = wwl(i,:)
         den(:) = denl(i,:)
         denfac(:) = denfacl(i,:)
         tk(:) = tkl(i,:)
      ! skip for no precipitation for all layers
         allold = 0.0
         do k=1,km
            allold = allold + qq(k) + qq2(k)
         enddo
         if(allold.le.0.0) then
            cycle i_loop
         endif
      !
      ! compute interface values
         zi(1)=0.0
         do k=1,km
            zi(k+1) = zi(k)+dz(k)
         enddo
      !
      ! save departure wind
         wd(:) = ww(:)
         n=1
      100  continue
      ! plm is 2nd order, we can use 2nd order wi or 3rd order wi
      ! 2nd order interpolation to get wi
         wi(1) = ww(1)
         wi(km+1) = ww(km)
         do k=2,km
            wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
         enddo
      ! 3rd order interpolation to get wi
         fa1 = 9./16.
         fa2 = 1./16.
         wi(1) = ww(1)
         wi(2) = 0.5*(ww(2)+ww(1))
         do k=3,km-1
            wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
         enddo
         wi(km) = 0.5*(ww(km)+ww(km-1))
         wi(km+1) = ww(km)
      !
      ! terminate of top of raingroup
         do k=2,km
            if( ww(k).eq.0.0 ) wi(k)=ww(k-1)
         enddo
      !
      ! diffusivity of wi
         con1 = 0.05
         do k=km,1,-1
            decfl = (wi(k+1)-wi(k))*dt/dz(k)
            if( decfl .gt. con1 ) then
               wi(k) = wi(k+1) - con1*dz(k)/dt
            endif
         enddo
      ! compute arrival point
         do k=1,km+1
            za(k) = zi(k) - wi(k)*dt
         enddo
      !
         do k=1,km
            dza(k) = za(k+1)-za(k)
         enddo
         dza(km+1) = zi(km+1) - za(km+1)
      !
      ! computer deformation at arrival point
         do k=1,km
            qa(k) = qq(k)*dz(k)/dza(k)
            qa2(k) = qq2(k)*dz(k)/dza(k)
            qr(k) = qa(k)/den(k)
            qr2(k) = qa2(k)/den(k)
         enddo
         qa(km+1) = 0.0
         qa2(km+1) = 0.0
      !   call maxmin(km,1,qa,' arrival points ')
      !
      ! compute arrival terminal velocity, and estimate mean terminal velocity
      ! then back to use mean terminal velocity
         if( n.le.iter ) then
            call slope_snow(qr,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa, &
                           1,1,1,km)
            call slope_graup(qr2,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa2,1,1,1,km)
            do k = 1, km
               tmp(k) = max((qr(k)+qr2(k)), 1.E-15)
               if( tmp(k) .gt. 1.e-15 ) then
                  wa(k) = (wa(k)*qr(k) + wa2(k)*qr2(k))/tmp(k)
               else
                  wa(k) = 0.
               endif
            enddo
            if( n.ge.2 ) wa(1:km)=0.5*(wa(1:km)+was(1:km))
            do k=1,km
      !#ifdef DEBUG
      !         print*,' slope_wsm3 ',qr(k)*1000.,den(k),denfac(k),tk(k),tmp(k),tmp1(k),tmp2(k), &
      !            ww(k),wa(k)
      !#endif
      ! mean wind is average of departure and new arrival winds
               ww(k) = 0.5* ( wd(k)+wa(k) )
            enddo
            was(:) = wa(:)
            n=n+1
            go to 100
         endif
      
         ist_loop : do ist = 1, 2
            if (ist.eq.2) then
               qa(:) = qa2(:)
            endif
      !
            precip(i) = 0.
      !
      ! estimate values at arrival cell interface with monotone
            do k=2,km
               dip=(qa(k+1)-qa(k))/(dza(k+1)+dza(k))
               dim=(qa(k)-qa(k-1))/(dza(k-1)+dza(k))
               if( dip*dim.le.0.0 ) then
                  qmi(k)=qa(k)
                  qpi(k)=qa(k)
               else
                  qpi(k)=qa(k)+0.5*(dip+dim)*dza(k)
                  qmi(k)=2.0*qa(k)-qpi(k)
                  if( qpi(k).lt.0.0 .or. qmi(k).lt.0.0 ) then
                     qpi(k) = qa(k)
                     qmi(k) = qa(k)
                  endif
               endif
            enddo
            qpi(1)=qa(1)
            qmi(1)=qa(1)
            qmi(km+1)=qa(km+1)
            qpi(km+1)=qa(km+1)
      !
      ! interpolation to regular point
            qn = 0.0
            kb=1
            kt=1
            intp : do k=1,km
               kb=max(kb-1,1)
               kt=max(kt-1,1)
      ! find kb and kt
               if( zi(k).ge.za(km+1) ) then
                  exit intp
               else
                  find_kb : do kk=kb,km
                     if( zi(k).le.za(kk+1) ) then
                        kb = kk
                        exit find_kb
                     else
                        cycle find_kb
                     endif
                  enddo find_kb
                  find_kt : do kk=kt,km
                     if( zi(k+1).le.za(kk) ) then
                        kt = kk
                        exit find_kt
                     else
                        cycle find_kt
                     endif
                  enddo find_kt
                  kt = kt - 1
      ! compute q with piecewise constant method
                  if( kt.eq.kb ) then
                     tl=(zi(k)-za(kb))/dza(kb)
                     th=(zi(k+1)-za(kb))/dza(kb)
                     tl2=tl*tl
                     th2=th*th
                     qqd=0.5*(qpi(kb)-qmi(kb))
                     qqh=qqd*th2+qmi(kb)*th
                     qql=qqd*tl2+qmi(kb)*tl
                     qn(k) = (qqh-qql)/(th-tl)
                  else if( kt.gt.kb ) then
                     tl=(zi(k)-za(kb))/dza(kb)
                     tl2=tl*tl
                     qqd=0.5*(qpi(kb)-qmi(kb))
                     qql=qqd*tl2+qmi(kb)*tl
                     dql = qa(kb)-qql
                     zsum  = (1.-tl)*dza(kb)
                     qsum  = dql*dza(kb)
                     if( kt-kb.gt.1 ) then
                        do m=kb+1,kt-1
                           zsum = zsum + dza(m)
                           qsum = qsum + qa(m) * dza(m)
                        enddo
                     endif
                     th=(zi(k+1)-za(kt))/dza(kt)
                     th2=th*th
                     qqd=0.5*(qpi(kt)-qmi(kt))
                     dqh=qqd*th2+qmi(kt)*th
                     zsum  = zsum + th*dza(kt)
                     qsum  = qsum + dqh*dza(kt)
                     qn(k) = qsum/zsum
                  endif
                  cycle intp
               endif
      !
            enddo intp
      !
      ! rain out
            sum_precip: do k=1,km
               if( za(k).lt.0.0 .and. za(k+1).lt.0.0 ) then
                  precip(i) = precip(i) + qa(k)*dza(k)
                  cycle sum_precip
               else if ( za(k).lt.0.0 .and. za(k+1).ge.0.0 ) then
                  precip(i) = precip(i) + qa(k)*(0.0-za(k))
                  exit sum_precip
               endif
               exit sum_precip
            enddo sum_precip
      !
      ! replace the new values
            if(ist.eq.1) then
               rql(i,:) = qn(:)
               precip1(i) = precip(i)
            else
               rql2(i,:) = qn(:)
               precip2(i) = precip(i)
            endif
         enddo ist_loop
      
      enddo i_loop
      
      end subroutine nislfv_rain_plm6
   
   !=================================================================================================================
      subroutine refl10cm_wsm6(qv1d,qr1d,qs1d,qg1d,t1d,p1d,dBZ,kts,kte)
      implicit none
   !=================================================================================================================
   
   !..Sub arguments
      integer,intent(in):: kts,kte
      real(KIND=R8SIZE),intent(in),dimension(kts:kte):: qv1d,qr1d,qs1d,qg1d,t1d,p1d
      real(KIND=R8SIZE),intent(inout),dimension(kts:kte):: dBz
   
   !..Local variables
      logical:: melti
      logical,dimension(kts:kte):: l_qr,l_qs,l_qg
   
      INTEGER:: i,k,k_0,kbot,n
   
      real(KIND=R8SIZE),parameter:: R=287.
      real(KIND=R8SIZE):: temp_c
      real(KIND=R8SIZE),dimension(kts:kte):: temp,pres,qv,rho
      real(KIND=R8SIZE),dimension(kts:kte):: rr,rs,rg
      real(KIND=R8SIZE),dimension(kts:kte):: ze_rain,ze_snow,ze_graupel
   
      double precision:: fmelt_s,fmelt_g
      double precision:: cback,x,eta,f_d
      double precision,dimension(kts:kte):: ilamr,ilams,ilamg
      double precision,dimension(kts:kte):: n0_r, n0_s, n0_g
      double precision:: lamr,lams,lamg
   
   !-----------------------------------------------------------------------------------------------------------------
   
      do k = kts, kte
         dBZ(k) = -35.0
      enddo
   
   !+---+-----------------------------------------------------------------+
   !..Put column of data into local arrays.
   !+---+-----------------------------------------------------------------+
      do k = kts, kte
         temp(k) = t1d(k)
         temp_c = min(-0.001, temp(k)-273.15)
         qv(k) = max(1.e-10, qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
      
         if (qr1d(k) .gt. 1.e-9) then
            rr(k) = qr1d(k)*rho(k)
            n0_r(k) = n0r
            lamr = (xam_r*xcrg(3)*n0_r(k)/rr(k))**(1./xcre(1))
            ilamr(k) = 1./lamr
            l_qr(k) = .true.
         else
            rr(k) = 1.e-12
            l_qr(k) = .false.
         endif
      
         if (qs1d(k) .gt. 1.e-9) then
            rs(k) = qs1d(k)*rho(k)
            n0_s(k) = min(n0smax, n0s*exp(-alpha*temp_c))
            lams = (xam_s*xcsg(3)*n0_s(k)/rs(k))**(1./xcse(1))
            ilams(k) = 1./lams
            l_qs(k) = .true.
         else
            rs(k) = 1.e-12
            l_qs(k) = .false.
         endif
      
         if (qg1d(k) .gt. 1.e-9) then
            rg(k) = qg1d(k)*rho(k)
            n0_g(k) = n0g
            lamg = (xam_g*xcgg(3)*n0_g(k)/rg(k))**(1./xcge(1))
            ilamg(k) = 1./lamg
            l_qg(k) = .true.
         else
            rg(k) = 1.e-12
            l_qg(k) = .false.
         endif
      enddo
   
   !+---+-----------------------------------------------------------------+
   !..Locate K-level of start of melting (k_0 is level above).
   !+---+-----------------------------------------------------------------+
            melti = .false.
            k_0 = kts
            do k = kte-1, kts, -1
               if ( (temp(k).gt.273.15) .and. L_qr(k)                         &
                                       .and. (L_qs(k+1).or.L_qg(k+1)) ) then
                  k_0 = MAX(k+1, k_0)
                  melti=.true.
                  goto 195
               endif
            enddo
      195  continue
      
   !+---+-----------------------------------------------------------------+
   !..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
   !.. and non-water-coated snow and graupel when below freezing are
   !.. simple. Integrations of m(D)*m(D)*N(D)*dD.
   !+---+-----------------------------------------------------------------+
   
      do k = kts, kte
         ze_rain(k) = 1.e-22
         ze_snow(k) = 1.e-22
         ze_graupel(k) = 1.e-22
         if (l_qr(k)) ze_rain(k) = n0_r(k)*xcrg(4)*ilamr(k)**xcre(4)
         if (l_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/pi)*(6.0/pi)     &
                                 * (xam_s/900.0)*(xam_s/900.0)          &
                                 * n0_s(k)*xcsg(4)*ilams(k)**xcse(4)
         if (l_qg(k)) ze_graupel(k) = (0.176/0.93) * (6.0/pi)*(6.0/pi)  &
                                 * (xam_g/900.0)*(xam_g/900.0)       &
                                 * n0_g(k)*xcgg(4)*ilamg(k)**xcge(4)
      enddo
   
   
   !+---+-----------------------------------------------------------------+
   !..Special case of melting ice (snow/graupel) particles.  Assume the
   !.. ice is surrounded by the liquid water.  Fraction of meltwater is
   !.. extremely simple based on amount found above the melting level.
   !.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
   !.. routines).
   !+---+-----------------------------------------------------------------+
   
      if (melti .and. k_0.ge.kts+1) then
         do k = k_0-1, kts, -1
      
      !..Reflectivity contributed by melting snow
            if (L_qs(k) .and. L_qs(k_0) ) then
               fmelt_s = MAX(0.005d0, MIN(1.0d0-rs(k)/rs(k_0), 0.99d0))
               eta = 0.d0
               lams = 1./ilams(k)
               do n = 1, nrbins
                  x = xam_s * xxDs(n)**xbm_s
                  call rayleigh_soak_wetgraupel (x,dble(xocms),dble(xobms), &
                        fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                        cback, mixingrulestring_s, matrixstring_s,          &
                        inclusionstring_s, hoststring_s,                    &
                        hostmatrixstring_s, hostinclusionstring_s)
                  f_d = n0_s(k)*xxds(n)**xmu_s * dexp(-lams*xxds(n))
                  eta = eta + f_d * cback * simpson(n) * xdts(n)
               enddo
               ze_snow(k) = sngl(lamda4 / (pi5 * k_w) * eta)
            endif
      
      
      !..Reflectivity contributed by melting graupel
      
            if (l_qg(k) .and. l_qg(k_0) ) then
               fmelt_g = max(0.005d0, min(1.0d0-rg(k)/rg(k_0), 0.99d0))
               eta = 0.d0
               lamg = 1./ilamg(k)
               do n = 1, nrbins
                  x = xam_g * xxdg(n)**xbm_g
                  call rayleigh_soak_wetgraupel (x,dble(xocmg),dble(xobmg), &
                        fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
                        cback, mixingrulestring_g, matrixstring_g,          &
                        inclusionstring_g, hoststring_g,                    &
                        hostmatrixstring_g, hostinclusionstring_g)
                  f_d = n0_g(k)*xxdg(n)**xmu_g * dexp(-lamg*xxdg(n))
                  eta = eta + f_d * cback * simpson(n) * xdtg(n)
               enddo
               ze_graupel(k) = sngl(lamda4 / (pi5 * k_w) * eta)
            endif
      
         enddo
      endif
      
      do k = kte, kts, -1
         dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
      enddo
      
      
      end subroutine refl10cm_wsm6

!================================================================================================================

      !=================================================================================================================
      subroutine mp_wsm6_effectRad_run(do_microp_re,t,qc,qi,qs,rho,qmin,t0c,re_qc_bg,re_qi_bg,re_qs_bg, &
         re_qc_max,re_qi_max,re_qs_max,re_qc,re_qi,re_qs,its,ite,kts,kte, &
         errmsg,errflg)
!=================================================================================================================
!  Compute radiation effective radii of cloud water, ice, and snow for
!  single-moment microphysics.
!  These are entirely consistent with microphysics assumptions, not
!  constant or otherwise ad hoc as is internal to most radiation
!  schemes.
!  Coded and implemented by Soo ya Bae, KIAPS, January 2015.
!-----------------------------------------------------------------------------------------------------------------


!..Sub arguments
logical,intent(in):: do_microp_re
integer,intent(in):: its,ite,kts,kte
real(KIND=R8SIZE),intent(in):: qmin
real(KIND=R8SIZE),intent(in):: t0c
real(KIND=R4SIZE),intent(in):: re_qc_bg,re_qi_bg,re_qs_bg !era r8size - mudei 15/03/2024
real(KIND=R4SIZE),intent(in):: re_qc_max,re_qi_max,re_qs_max !era r8size - mudei 15/03/2024
real(KIND=R8SIZE),dimension(its:ite,kts:kte),intent(in)::  t
real(KIND=R8SIZE),dimension(its:ite,kts:kte),intent(in)::  qc
real(KIND=R8SIZE),dimension(its:ite,kts:kte),intent(in)::  qi
real(KIND=R8SIZE),dimension(its:ite,kts:kte),intent(in)::  qs
real(KIND=R8SIZE),dimension(its:ite,kts:kte),intent(in)::  rho
real(KIND=R8SIZE),dimension(its:ite,kts:kte),intent(inout):: re_qc
real(KIND=R8SIZE),dimension(its:ite,kts:kte),intent(inout):: re_qi
real(KIND=R8SIZE),dimension(its:ite,kts:kte),intent(inout):: re_qs

!...Output arguments:
character(len=*),intent(out):: errmsg
integer,intent(out):: errflg

!..Local variables
integer:: i,k
integer:: inu_c
real(KIND=R8SIZE),dimension(its:ite,kts:kte):: ni
real(KIND=R8SIZE),dimension(its:ite,kts:kte):: rqc
real(KIND=R8SIZE),dimension(its:ite,kts:kte):: rqi
real(KIND=R8SIZE),dimension(its:ite,kts:kte):: rni
real(KIND=R8SIZE),dimension(its:ite,kts:kte):: rqs
real(KIND=R8SIZE):: temp
real(KIND=R8SIZE):: lamdac
real(KIND=R8SIZE):: supcol,n0sfac,lamdas
real(KIND=R8SIZE):: diai      ! diameter of ice in m
logical:: has_qc, has_qi, has_qs
!..Minimum microphys values
real(KIND=R8SIZE),parameter:: R1 = 1.E-12
real(KIND=R8SIZE),parameter:: R2 = 1.E-6
!..Mass power law relations:  mass = am*D**bm
real(KIND=R8SIZE),parameter:: bm_r = 3.0
real(KIND=R8SIZE),parameter:: obmr = 1.0/bm_r
real(KIND=R8SIZE),parameter:: nc0  = 3.E8

!-----------------------------------------------------------------------------------------------------------------

if(.not. do_microp_re) return

!--- initialization of effective radii of cloud water, cloud ice, and snow to background values:
do k = kts,kte
do i = its,ite
re_qc(i,k) = re_qc_bg
re_qi(i,k) = re_qi_bg
re_qs(i,k) = re_qs_bg
enddo
enddo

!--- computation of effective radii:
has_qc = .false.
has_qi = .false.
has_qs = .false.

do k = kts,kte
do i = its,ite
! for cloud
rqc(i,k) = max(R1,qc(i,k)*rho(i,k))
if (rqc(i,k).gt.R1) has_qc = .true.
! for ice
rqi(i,k) = max(R1,qi(i,k)*rho(i,k))
temp = (rho(i,k)*max(qi(i,k),qmin))
temp = sqrt(sqrt(temp*temp*temp))
ni(i,k) = min(max(5.38e7*temp,1.e3),1.e6)
rni(i,k)= max(R2,ni(i,k)*rho(i,k))
if (rqi(i,k).gt.R1 .and. rni(i,k).gt.R2) has_qi = .true.
! for snow
rqs(i,k) = max(R1,qs(i,k)*rho(i,k))
if (rqs(i,k).gt.R1) has_qs = .true.
enddo
enddo

if (has_qc) then
do k = kts,kte
do i = its,ite
if (rqc(i,k).le.R1) CYCLE
lamdac = (pidnc*nc0/rqc(i,k))**obmr
re_qc(i,k) =  max(2.51E-6,min(1.5*(1.0/lamdac),re_qc_max))
enddo
enddo
endif

if (has_qi) then
do k = kts,kte
do i = its,ite
if (rqi(i,k).le.R1 .or. rni(i,k).le.R2) CYCLE
diai = 11.9*sqrt(rqi(i,k)/ni(i,k))
re_qi(i,k) = max(10.01E-6,min(0.75*0.163*diai,re_qi_max))
enddo
enddo
endif

if (has_qs) then
do i = its,ite
do k = kts,kte
if (rqs(i,k).le.R1) CYCLE
supcol = t0c-t(i,k)
n0sfac = max(min(exp(alpha*supcol),n0smax/n0s),1.)
lamdas = sqrt(sqrt(pidn0s*n0sfac/rqs(i,k)))
re_qs(i,k) = max(25.E-6,min(0.5*(1./lamdas),re_qs_max))
enddo
enddo
endif

!--- limit effective radii of cloud water, cloud ice, and snow to maximum values:
do k = kts,kte
do i = its,ite
re_qc(i,k) = max(re_qc_bg,min(re_qc(i,k),re_qc_max))
re_qi(i,k) = max(re_qi_bg,min(re_qi(i,k),re_qi_max))
re_qs(i,k) = max(re_qs_bg,min(re_qs(i,k),re_qs_max))
enddo
enddo

errmsg = 'mp_wsm6_effectRad_run OK'
errflg = 0

end subroutine mp_wsm6_effectRad_run

!=================================================================================================================
! Vem tudo do modulo mp_Radar
!=================================================================================================================


!contains


!=================================================================================================================
subroutine radar_init
!implicit none
!=================================================================================================================

integer:: n

!-----------------------------------------------------------------------------------------------------------------

pi5 = 3.14159*3.14159*3.14159*3.14159*3.14159
lamda4 = lamda_radar*lamda_radar*lamda_radar*lamda_radar
m_w_0 = m_complex_water_ray (lamda_radar, 0.0d0)
m_i_0 = m_complex_ice_maetzler (lamda_radar, 0.0d0)
k_w = (abs( (m_w_0*m_w_0 - 1.0) /(m_w_0*m_w_0 + 2.0) ))**2

do n = 1, nrbins+1
   simpson(n) = 0.0d0
enddo
do n = 1, nrbins-1, 2
   simpson(n) = simpson(n) + basis(1)
   simpson(n+1) = simpson(n+1) + basis(2)
   simpson(n+2) = simpson(n+2) + basis(3)
enddo

do n = 1, slen
   mixingrulestring_s(n:n) = char(0)
   matrixstring_s(n:n) = char(0)
   inclusionstring_s(n:n) = char(0)
   hoststring_s(n:n) = char(0)
   hostmatrixstring_s(n:n) = char(0)
   hostinclusionstring_s(n:n) = char(0)
   mixingrulestring_g(n:n) = char(0)
   matrixstring_g(n:n) = char(0)
   inclusionstring_g(n:n) = char(0)
   hoststring_g(n:n) = char(0)
   hostmatrixstring_g(n:n) = char(0)
   hostinclusionstring_g(n:n) = char(0)
enddo

mixingrulestring_s = 'maxwellgarnett'
hoststring_s = 'air'
matrixstring_s = 'water'
inclusionstring_s = 'spheroidal'
hostmatrixstring_s = 'icewater'
hostinclusionstring_s = 'spheroidal'

mixingrulestring_g = 'maxwellgarnett'
hoststring_g = 'air'
matrixstring_g = 'water'
inclusionstring_g = 'spheroidal'
hostmatrixstring_g = 'icewater'
hostinclusionstring_g = 'spheroidal'

!..Create bins of snow (from 100 microns up to 2 cm).
xxdx(1) = 100.d-6
xxdx(nrbins+1) = 0.02d0
do n = 2, nrbins
   xxdx(n) = dexp(real(n-1,KIND=R8SIZE)/real(nrbins,KIND=R8SIZE) &
            * dlog(xxdx(nrbins+1)/xxdx(1)) +dlog(xxdx(1)))
enddo
do n = 1, nrbins
   xxds(n) = dsqrt(xxdx(n)*xxdx(n+1))
   xdts(n) = xxdx(n+1) - xxdx(n)
enddo

!..create bins of graupel (from 100 microns up to 5 cm).
xxdx(1) = 100.d-6
xxdx(nrbins+1) = 0.05d0
do n = 2, nrbins
   xxdx(n) = dexp(real(n-1,KIND=R8SIZE)/real(nrbins,KIND=R8SIZE) &
            * dlog(xxdx(nrbins+1)/xxdx(1)) +dlog(xxdx(1)))
enddo
do n = 1, nrbins
   xxdg(n) = dsqrt(xxdx(n)*xxdx(n+1))
   xdtg(n) = xxdx(n+1) - xxdx(n)
enddo


!.. The calling program must set the m(D) relations and gamma shape
!.. parameter mu for rain, snow, and graupel.  Easily add other types
!.. based on the template here.  For majority of schemes with simpler
!.. exponential number distribution, mu=0.

xcre(1) = 1. + xbm_r
xcre(2) = 1. + xmu_r
xcre(3) = 4. + xmu_r
xcre(4) = 7. + xmu_r
do n = 1, 4
   xcrg(n) = wgamma(xcre(n))
enddo
xorg2 = 1./xcrg(2)

xcse(1) = 1. + xbm_s
xcse(2) = 1. + xmu_s
xcse(3) = 4. + xmu_s
xcse(4) = 7. + xmu_s
do n = 1, 4
   xcsg(n) = wgamma(xcse(n))
enddo
xosg2 = 1./xcsg(2)

xcge(1) = 1. + xbm_g
xcge(2) = 1. + xmu_g
xcge(3) = 4. + xmu_g
xcge(4) = 7. + xmu_g
do n = 1, 4
   xcgg(n) = wgamma(xcge(n))
enddo
xogg2 = 1./xcgg(2)

xobmr = 1./xbm_r
xoams = 1./xam_s
xobms = 1./xbm_s
xocms = xoams**xobms
xoamg = 1./xam_g
xobmg = 1./xbm_g
xocmg = xoamg**xobmg

end subroutine radar_init

!=================================================================================================================
subroutine rayleigh_soak_wetgraupel(x_g,a_geo,b_geo,fmelt, &
                                    meltratio_outside, &
                                    m_w,m_i,lambda,c_back, &
                                    mixingrule,matrix, &
                                    inclusion, &
                                    host,hostmatrix,hostinclusion)
!implicit none
!=================================================================================================================

!--- input arguments:
character(len=*), intent(in):: mixingrule, matrix, inclusion, &
                              host, hostmatrix, hostinclusion

complex(KIND=R8SIZE),intent(in):: m_w, m_i

double precision, intent(in):: x_g, a_geo, b_geo, fmelt, lambda, meltratio_outside

!--- output arguments:
double precision,intent(out):: c_back

!--- local variables:
integer:: error

complex(KIND=R8SIZE):: m_core, m_air

double precision, parameter:: pix=3.1415926535897932384626434d0
double precision:: d_large, d_g, rhog, x_w, xw_a, fm, fmgrenz, &
                        volg, vg, volair, volice, volwater,    &
                        meltratio_outside_grenz, mra

!-----------------------------------------------------------------------------------------------------------------

!refractive index of air:
m_air = (1.0d0,0.0d0)

!Limiting the degree of melting --- for safety: 
fm = dmax1(dmin1(fmelt, 1.0d0), 0.0d0)
!Limiting the ratio of (melting on outside)/(melting on inside):
mra = dmax1(dmin1(meltratio_outside, 1.0d0), 0.0d0)

!The relative portion of meltwater melting at outside should increase
!from the given input value (between 0 and 1)
!to 1 as the degree of melting approaches 1,
!so that the melting particle "converges" to a water drop.
!Simplest assumption is linear:
mra = mra + (1.0d0-mra)*fm

x_w = x_g * fm

d_g = a_geo * x_g**b_geo

if(D_g .ge. 1d-12) then

   vg = PIx/6. * D_g**3
   rhog = DMAX1(DMIN1(x_g / vg, 900.0d0), 10.0d0)
   vg = x_g / rhog
      
   meltratio_outside_grenz = 1.0d0 - rhog / 1000.

   if (mra .le. meltratio_outside_grenz) then
      !..In this case, it cannot happen that, during melting, all the
      !.. air inclusions within the ice particle get filled with
      !.. meltwater. This only happens at the end of all melting.
      volg = vg * (1.0d0 - mra * fm)

   else
      !..In this case, at some melting degree fm, all the air
      !.. inclusions get filled with meltwater.
      fmgrenz=(900.0-rhog)/(mra*900.0-rhog+900.0*rhog/1000.)

      if (fm .le. fmgrenz) then
         !.. not all air pockets are filled:
         volg = (1.0 - mra * fm) * vg
      else
         !..all air pockets are filled with meltwater, now the
         !.. entire ice sceleton melts homogeneously:
         volg = (x_g - x_w) / 900.0 + x_w / 1000.
      endif

   endif

   d_large  = (6.0 / pix * volg) ** (1./3.)
   volice = (x_g - x_w) / (volg * 900.0)
   volwater = x_w / (1000. * volg)
   volair = 1.0 - volice - volwater
      
   !..complex index of refraction for the ice-air-water mixture
   !.. of the particle:
   m_core = get_m_mix_nested (m_air, m_i, m_w, volair, volice,      &
                     volwater, mixingrule, host, matrix, inclusion, &
                     hostmatrix, hostinclusion, error)
   if (error .ne. 0) then
      c_back = 0.0d0
      return
   endif

   !..rayleigh-backscattering coefficient of melting particle: 
   c_back = (abs((m_core**2-1.0d0)/(m_core**2+2.0d0)))**2           &
         * pi5 * d_large**6 / lamda4

else
   c_back = 0.0d0
endif

end subroutine rayleigh_soak_wetgraupel

!=================================================================================================================
real(KIND=R8SIZE) function wgamma(y)
implicit none
!=================================================================================================================

!--- input arguments:
real(KIND=R8SIZE),intent(in):: y

!-----------------------------------------------------------------------------------------------------------------

wgamma = exp(gammln(y))

end function wgamma

!=================================================================================================================
real(KIND=R8SIZE) function gammln(xx)
implicit none
!(C) Copr. 1986-92 Numerical Recipes Software 2.02
!=================================================================================================================

!--- inout arguments: 
real(KIND=R8SIZE),intent(in):: xx

!--- local variables:
integer:: j

double precision,parameter:: stp = 2.5066282746310005d0
double precision,dimension(6), parameter:: &
   cof = (/76.18009172947146d0, -86.50532032941677d0, &
            24.01409824083091d0, -1.231739572450155d0, &
            .1208650973866179d-2, -.5395239384953d-5/)
double precision:: ser,tmp,x,y

!-----------------------------------------------------------------------------------------------------------------

!--- returns the value ln(gamma(xx)) for xx > 0.
x = xx
y = x
tmp = x+5.5d0
tmp = (x+0.5d0)*log(tmp)-tmp
ser = 1.000000000190015d0
do j = 1,6
   y=y+1.d0
   ser=ser+cof(j)/y
enddo

gammln=tmp+log(stp*ser/x)

end function gammln
      
!=================================================================================================================
complex(KIND=R8SIZE) function get_m_mix_nested (m_a, m_i, m_w, volair,      &
               volice, volwater, mixingrule, host, matrix,        &
               inclusion, hostmatrix, hostinclusion, cumulerror)
implicit none
!=================================================================================================================

!--- input arguments:
character(len=*),intent(in):: mixingrule, host, matrix,          &
                  inclusion, hostmatrix, hostinclusion

complex(KIND=R8SIZE),intent(in):: m_a, m_i, m_w

double precision,intent(in):: volice, volair, volwater

!--- output arguments:
integer,intent(out):: cumulerror

!--- local variables:
integer:: error

complex(KIND=R8SIZE):: mtmp

double precision:: vol1, vol2

!-----------------------------------------------------------------------------------------------------------------

!..Folded: ( (m1 + m2) + m3), where m1,m2,m3 could each be air, ice, or water
cumulerror = 0
get_m_mix_nested = cmplx(1.0d0,0.0d0)

if (host .eq. 'air') then
   if (matrix .eq. 'air') then
      write(radar_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
      !call physics_message(radar_debug)
      cumulerror = cumulerror + 1
   else
      vol1 = volice / MAX(volice+volwater,1d-10)
      vol2 = 1.0d0 - vol1
      mtmp = get_m_mix (m_a, m_i, m_w, 0.0d0, vol1, vol2,              &
                        mixingrule, matrix, inclusion, error)
      cumulerror = cumulerror + error
         
      if (hostmatrix .eq. 'air') then
         get_m_mix_nested = get_m_mix (m_a, mtmp, 2.0*m_a,              &
                           volair, (1.0d0-volair), 0.0d0, mixingrule,  &
                           hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
      elseif (hostmatrix .eq. 'icewater') then
         get_m_mix_nested = get_m_mix (m_a, mtmp, 2.0*m_a,              &
                           volair, (1.0d0-volair), 0.0d0, mixingrule,  &
                           'ice', hostinclusion, error)
         cumulerror = cumulerror + error
      else
         write(radar_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ', hostmatrix
         !call physics_message(radar_debug)
         cumulerror = cumulerror + 1
      endif
   endif

elseif (host .eq. 'ice') then

   if (matrix .eq. 'ice') then
      write(radar_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
      !call physics_message(radar_debug)
      cumulerror = cumulerror + 1
   else
      vol1 = volair / MAX(volair+volwater,1d-10)
      vol2 = 1.0d0 - vol1
      mtmp = get_m_mix (m_a, m_i, m_w, vol1, 0.0d0, vol2,               &
                        mixingrule, matrix, inclusion, error)
      cumulerror = cumulerror + error

      if (hostmatrix .eq. 'ice') then
         get_m_mix_nested = get_m_mix (mtmp, m_i, 2.0*m_a,              &
                           (1.0d0-volice), volice, 0.0d0, mixingrule,  &
                           hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
      elseif (hostmatrix .eq. 'airwater') then
         get_m_mix_nested = get_m_mix (mtmp, m_i, 2.0*m_a,              &
                           (1.0d0-volice), volice, 0.0d0, mixingrule,  &
                           'air', hostinclusion, error)
         cumulerror = cumulerror + error          
      else
         write(radar_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ', hostmatrix
         !call physics_message(radar_debug)
         cumulerror = cumulerror + 1
      endif
   endif

elseif (host .eq. 'water') then

   if (matrix .eq. 'water') then
      write(radar_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
      !call physics_message(radar_debug)
      cumulerror = cumulerror + 1
   else
      vol1 = volair / MAX(volice+volair,1d-10)
      vol2 = 1.0d0 - vol1
      mtmp = get_m_mix (m_a, m_i, m_w, vol1, vol2, 0.0d0,               &
                        mixingrule, matrix, inclusion, error)
      cumulerror = cumulerror + error

      if (hostmatrix .eq. 'water') then
         get_m_mix_nested = get_m_mix (2*m_a, mtmp, m_w,                &
                        0.0d0, (1.0d0-volwater), volwater, mixingrule,  &
                        hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
      elseif (hostmatrix .eq. 'airice') then
         get_m_mix_nested = get_m_mix (2*m_a, mtmp, m_w,                &
                        0.0d0, (1.0d0-volwater), volwater, mixingrule,  &
                        'ice', hostinclusion, error)
         cumulerror = cumulerror + error          
      else
         write(radar_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ', hostmatrix
         !call physics_message(radar_debug)
         cumulerror = cumulerror + 1
      endif
   endif

elseif (host .eq. 'none') then

   get_m_mix_nested = get_m_mix (m_a, m_i, m_w,                         &
                     volair, volice, volwater, mixingrule,             &
                     matrix, inclusion, error)
   cumulerror = cumulerror + error
      
else
   write(radar_debug,*) 'GET_M_MIX_NESTED: unknown matrix: ', host
   !call physics_message(radar_debug)
   cumulerror = cumulerror + 1
endif

if (cumulerror .ne. 0) then
   write(radar_debug,*) 'get_m_mix_nested: error encountered'
   !call physics_message(radar_debug)
   get_m_mix_nested = cmplx(1.0d0,0.0d0)    
endif

end function get_m_mix_nested

!=================================================================================================================
complex(KIND=R8SIZE) function get_m_mix (m_a, m_i, m_w, volair, volice, &
                     volwater, mixingrule, matrix, inclusion, &
                     error)
implicit none
!=================================================================================================================

!--- input arguments:
character(len=*),intent(in):: mixingrule, matrix, inclusion

complex(KIND=R8SIZE), intent(in):: m_a, m_i, m_w

double precision, intent(in):: volice, volair, volwater

!--- output arguments:
integer,intent(out):: error

!-----------------------------------------------------------------------------------------------------------------
error = 0
get_m_mix = cmplx(1.0d0,0.0d0)

if (mixingrule .eq. 'maxwellgarnett') then
   if (matrix .eq. 'ice') then
      get_m_mix = m_complex_maxwellgarnett(volice, volair, volwater,  &
                  m_i, m_a, m_w, inclusion, error)
   elseif (matrix .eq. 'water') then
      get_m_mix = m_complex_maxwellgarnett(volwater, volair, volice,  &
                  m_w, m_a, m_i, inclusion, error)
   elseif (matrix .eq. 'air') then
      get_m_mix = m_complex_maxwellgarnett(volair, volwater, volice,  &
                  m_a, m_w, m_i, inclusion, error)
   else
      write(radar_debug,*) 'GET_M_MIX: unknown matrix: ', matrix
      !call physics_message(radar_debug)
      error = 1
   endif

else
   write(radar_debug,*) 'GET_M_MIX: unknown mixingrule: ', mixingrule
   !call physics_message(radar_debug)
   error = 2
endif

if (error .ne. 0) then
   write(radar_debug,*) 'GET_M_MIX: error encountered'
   !call physics_message(radar_debug)
endif

end function get_m_mix

!=================================================================================================================
complex(KIND=R8SIZE) function m_complex_maxwellgarnett(vol1, vol2, vol3, &
                                 m1, m2, m3, inclusion, error)
implicit none
!=================================================================================================================

!--- input arguments:
character(len=*),intent(in):: inclusion

complex(KIND=R8SIZE),intent(in):: m1,m2,m3

double precision,intent(in):: vol1,vol2,vol3


!--- output arguments:
integer,intent(out):: error

!--- local variables:
complex(KIND=R8SIZE) :: beta2, beta3, m1t, m2t, m3t

!-----------------------------------------------------------------------------------------------------------------

error = 0

if (dabs(vol1+vol2+vol3-1.0d0) .gt. 1d-6) then
   write(radar_debug,*) 'M_COMPLEX_MAXWELLGARNETT: sum of the ', &
                  'partial volume fractions is not 1...ERROR'
   !call physics_message(radar_debug)
   m_complex_maxwellgarnett = CMPLX(-999.99d0,-999.99d0)
   error = 1
   return
endif

m1t = m1**2
m2t = m2**2
m3t = m3**2

if (inclusion .eq. 'spherical') then
   beta2 = 3.0d0*m1t/(m2t+2.0d0*m1t)
   beta3 = 3.0d0*m1t/(m3t+2.0d0*m1t)
elseif (inclusion .eq. 'spheroidal') then
   beta2 = 2.0d0*m1t/(m2t-m1t) * (m2t/(m2t-m1t)*LOG(m2t/m1t)-1.0d0)
   beta3 = 2.0d0*m1t/(m3t-m1t) * (m3t/(m3t-m1t)*LOG(m3t/m1t)-1.0d0)
else
   write(radar_debug,*) 'M_COMPLEX_MAXWELLGARNETT: ', 'unknown inclusion: ', inclusion
   !call physics_message(radar_debug)
   m_complex_maxwellgarnett=cmplx(-999.99d0,-999.99d0,KIND=R8SIZE)
   error = 1
   return
endif

m_complex_maxwellgarnett = sqrt(((1.0d0-vol2-vol3)*m1t + vol2*beta2*m2t + vol3*beta3*m3t) / &
                           (1.0d0-vol2-vol3+vol2*beta2+vol3*beta3))

end function m_complex_maxwellgarnett

!=================================================================================================================
complex(KIND=R8SIZE) function m_complex_water_ray(lambda,t)
implicit none
!=================================================================================================================

!complex refractive Index of Water as function of Temperature T
![deg C] and radar wavelength lambda [m]; valid for
!lambda in [0.001,1.0] m; T in [-10.0,30.0] deg C
!after Ray (1972)

!--- input arguments:
double precision,intent(in):: t,lambda

!--- local variables:
double precision,parameter:: pix=3.1415926535897932384626434d0
double precision:: epsinf,epss,epsr,epsi
double precision:: alpha,lambdas,sigma,nenner
complex(KIND=R8SIZE),parameter:: i = (0d0,1d0)

!-----------------------------------------------------------------------------------------------------------------

epsinf  = 5.27137d0 + 0.02164740d0 * T - 0.00131198d0 * T*T
epss    = 78.54d+0 * (1.0 - 4.579d-3 * (T - 25.0)                 &
         + 1.190d-5 * (T - 25.0)*(T - 25.0)                        &
         - 2.800d-8 * (T - 25.0)*(T - 25.0)*(T - 25.0))
alpha   = -16.8129d0/(T+273.16) + 0.0609265d0
lambdas = 0.00033836d0 * exp(2513.98d0/(T+273.16)) * 1e-2

nenner = 1.d0+2.d0*(lambdas/lambda)**(1d0-alpha)*sin(alpha*PIx*0.5) &
      + (lambdas/lambda)**(2d0-2d0*alpha)
epsr = epsinf + ((epss-epsinf) * ((lambdas/lambda)**(1d0-alpha)   &
      * sin(alpha*PIx*0.5)+1d0)) / nenner
epsi = ((epss-epsinf) * ((lambdas/lambda)**(1d0-alpha)            &
      * cos(alpha*PIx*0.5)+0d0)) / nenner                          &
      + lambda*1.25664/1.88496
      
m_complex_water_ray = sqrt(cmplx(epsr,-epsi))
      
end function m_complex_water_ray

!=================================================================================================================
complex(KIND=R8SIZE) function m_complex_ice_maetzler(lambda,t)
implicit none
!=================================================================================================================
      
!complex refractive index of ice as function of Temperature T
![deg C] and radar wavelength lambda [m]; valid for
!lambda in [0.0001,30] m; T in [-250.0,0.0] C
!Original comment from the Matlab-routine of Prof. Maetzler:
!Function for calculating the relative permittivity of pure ice in
!the microwave region, according to C. Maetzler, "Microwave
!properties of ice and snow", in B. Schmitt et al. (eds.) Solar
!System Ices, Astrophys. and Space Sci. Library, Vol. 227, Kluwer
!Academic Publishers, Dordrecht, pp. 241-257 (1998). Input:
!TK = temperature (K), range 20 to 273.15
!f = frequency in GHz, range 0.01 to 3000
         
!--- input arguments:
double precision,intent(in):: t,lambda

!--- local variables:
double precision:: f,c,tk,b1,b2,b,deltabeta,betam,beta,theta,alfa

!-----------------------------------------------------------------------------------------------------------------

c = 2.99d8
tk = t + 273.16
f = c / lambda * 1d-9

b1 = 0.0207
b2 = 1.16d-11
b = 335.0d0
deltabeta = exp(-10.02 + 0.0364*(tk-273.16))
betam = (b1/tk) * ( exp(b/tk) / ((exp(b/tk)-1)**2) ) + b2*f*f
beta = betam + deltabeta
theta = 300. / tk - 1.
alfa = (0.00504d0 + 0.0062d0*theta) * exp(-22.1d0*theta)
m_complex_ice_maetzler = 3.1884 + 9.1e-4*(tk-273.16)
m_complex_ice_maetzler = m_complex_ice_maetzler                   &
                        + cmplx(0.0d0, (alfa/f + beta*f)) 
m_complex_ice_maetzler = sqrt(conjg(m_complex_ice_maetzler))
      
end function m_complex_ice_maetzler
      

!=================================================================================================================
      end module module_mp_wsm6_2024
!=================================================================================================================
