module scm_gls_DGNG

   implicit none

contains

!----------------------------------------------------------------------------------------
subroutine gls_stp( Hz      ,        &            ! Cell thickness           [m]
                    u       ,        &            ! zonal velocity [m/s]
                    v       ,        &            ! meridional velocity [m/s]
                    varT    ,        &            ! Temperature variance (mean of T'^2)  [Celsius^2]
                    bvf     ,        &            ! Brunt-Vaisala frequency
                    alpha   ,        &            ! Thermal expansion coefficient in linear EOS
                    trb     ,        &            ! GLS variables TKE + length scale
                    lmix    ,        &            ! mixing length scale [m]
                    eps     ,        &            ! TKE dissipation [m2/s3]
                    Akv     ,        &            ! Turbulent viscosity  [m2/s]
                    Akt     ,        &            ! Turbulent diffusivity [m2/s]
                    gamma_GLS,       &            ! "Possibly CounterGradient" (PCG) term for tracers turbulent flux [m2/s]
                    c_mu    ,        &            ! Stability function for u,v
                    c_mu_prime ,     &            ! Stability function for T,S (DownGradient part)
                    c_mu_prime_star, &            ! Stability function for T,S ("Possibly CounterGradient" part, i.e. gamma)
                    alpha_n,         &            ! Relative to the stratification (Coeff in stability functions)
                    alpha_m,         &            ! Relative to the shear of mean currents (Coeff in stability functions)
                    alpha_v,         &            ! Relative to the temperature variance (Coeff in stability functions)
                    r_D     ,        &
                    sustr   ,        &            ! zonal wind stress [m2/s2 = (N/m2) / (kg/m3) ]
                    svstr   ,        &            ! meridional wind stress  [m2/s2 = (N/m2) / (kg/m3)]
                    gls_scheme ,     &
                    sfunc_opt  ,     &
                    dt      ,        &            ! Time-step
                    Zob     ,        &            ! bottom roughness length
                    Neu_bot ,        &            ! Nature of bottom boundary condition
                    kt      ,        &            ! time step of the external 1D simulation performed in Python
                    nstp    ,        &            ! time n
                    nnew    ,        &            ! time n+1
                    N       ,        &            ! Number of vertical grid points        )
                    ntra    ,        &
                    ngls    ,        &
                    ntime          )
       !------------------------------------------------------------------------
       integer,                              intent(in   ) :: N,ntra,ntime,ngls
       integer,                              intent(in   ) :: nstp
       integer,                              intent(in   ) :: nnew
       integer,                              intent(in   ) :: kt
       integer,                              intent(in   ) :: gls_scheme
       integer,                              intent(in   ) :: sfunc_opt
       real(8),dimension( 1:N, ntime      ), intent(in   ) :: u
       real(8),dimension( 1:N, ntime      ), intent(in   ) :: v
       real(8),dimension( 0:N, ntime      ), intent(in   ) :: varT
       real(8),dimension( 0:N             ), intent(in   ) :: bvf
       real(8),dimension( 0:N             ), intent(inout) :: Akv
       real(8),dimension( 0:N, ntra       ), intent(inout) :: Akt
       real(8),dimension( 0:N, ntra       ), intent(inout) :: gamma_GLS
       real(8),dimension( 0:N, ntime, ngls), intent(inout) :: trb
       real(8),dimension( 0:N             ), intent(inout) :: lmix
       real(8),dimension( 0:N             ), intent(inout) :: eps
       real(8),dimension( 0:N             ), intent(inout) :: c_mu
       real(8),dimension( 0:N             ), intent(inout) :: c_mu_prime
       real(8),dimension( 0:N             ), intent(inout) :: c_mu_prime_star
       real(8),dimension( 0:N             ), intent(inout) :: alpha_n
       real(8),dimension( 0:N             ), intent(inout) :: alpha_m
       real(8),dimension( 0:N             ), intent(inout) :: alpha_v
       real(8),                              intent(in   ) :: alpha
! Grid variables
       real(8),dimension( 1:N      ),        intent(in   ) ::  Hz
       real(8),                              intent(in   ) ::  sustr
       real(8),                              intent(in   ) ::  svstr
       real(8),                              intent(in   ) ::  dt
       real(8),                              intent(in   ) ::  r_D,Zob
       logical,                              intent(in   ) ::  Neu_bot
       !------------------------------------------------------------------------
! local variables
       integer   ::  k,ig,ig1,ig2,igls,itke,tind
       real(8)   :: diss  (1:N-1)
       real(8)   :: shear2(1:N-1)
       real(8)   :: FC    (1:N  )
       real(8)   :: DC    (1:N-1)
       real(8)   :: CF(1:N-1)
       real(8)   :: RH(1:N-1)
       real(8)   :: rp,    rm,    rn                !<-- n,m and p exponents
       real(8)   :: beta1, beta2, beta3m, beta3p    !<-- beta terms for the psi equation
       real(8)   :: OneOverSig(2)                   !<-- inverse of Schmidt number for tke and psi
       real(8)   :: e1,e2,e3
       real(8)   :: c1   ,c2    ,c3    ,c4    ,c5
       real(8)   :: c1T   ,c2T   ,c3T   ,c4T   ,cT
       real(8)   :: sf_d0 ,sf_d1 ,sf_d2 ,sf_d3 ,sf_d4 , sf_d5
       real(8)   :: sf_n0 ,sf_n1 ,sf_n2,sf_n3
       real(8)   :: sf_n0T,sf_n1T,sf_n2T
       real(8)   :: sf_n0T_star,sf_n1T_star,sf_n2T_star
       real(8)   :: lim_am0,lim_am1,lim_am2,lim_am3,lim_am4,lim_am5,lim_am6
       real(8)   :: z0_s,ustar_sfc_sq,ustar_bot_sq,L_lim,z0_b,trb_min(2)
       real(8)   :: cff,cff1,cff2,cff3m,cff3p,lgthsc,flux_top,flux_bot,trb_sfc,trb_bot
       real(8)   :: invk       , invG       , Bprod , Sprod    , epsilon
       real(8)   :: denom, gls_min
       real(8)   :: alpha_n_min, alpha_m_max, cm0   , cm0inv2  , gls, du, dv
       real(8),parameter :: vonKar = 0.4
       real(8),parameter :: nuwm   =   1e-04
       real(8),parameter :: nuws   = 0.1e-04
       real(8),parameter :: eps_bvf   =   1e-14
       real(8),parameter :: eps_min   = 1.0e-12
       real(8),parameter :: tke_min   = 1.0e-06
       real(8),parameter :: galp   =  0.53
       real(8),parameter ::  chk   =  1400./9.81
       real(8),parameter ::  g   =  9.81
       !--------------------------------------------------
       igls = 2
       itke = 1
       !--------------------------------------------------
       SELECT CASE( gls_scheme )
       CASE(1)     ! k-omega
          rp    = -1.0 ; rm    = 0.5  ; rn     = -1.0
          beta1 = 0.555; beta2 = 0.833; beta3m = -0.6; beta3p = 1.0
          OneOverSig = (/ 0.5, 0.5 /)
       CASE(2)     ! k-epsilon
          rp    = 3.0 ; rm    = 1.5 ; rn     = -1.0
          beta1 = 1.44; beta2 = 1.92; beta3m = -1.83; beta3p = -1.83
          OneOverSig = (/ 1.0, 0.8333 /)
          ! NEW
          ! 1) beta3m: -0.4 => -1.83, value -1.83 calculated from my side according to Umlauf 2003b (Title: "Extending the K-omega turbulence...")
          ! 2) beta3p = 1.0 => -1.83, better to take beta3p = beta3m, cf Umlauf 2005
          ! 3) OneOverSig = (/ 1.0, 0.7692 /)  =>  (/ 1.0, 0.8333 /)
          !    =>  0.7692 = 1/1.3    (sigma_epsilon = 1.3 is the value reported in Warner2005)
          !    =>  0.8333 = 1/1.2    (sigma_epsilon = 1.2 is the value found with (14) of Umlauf2003 with our value of c_mu0)
       CASE(3) ! gen-model
          rp    = 0.0; rm    = 1.0 ; rn     = -0.67
          beta1 = 1.0; beta2 = 1.22; beta3m =  0.05; beta3p = 1.0
          OneOverSig = (/ 1.25, 0.9345 /)
       CASE DEFAULT
          print*,'Error in the definition of the closure scheme'
          stop
       END SELECT
       e1 =  3.0 + 1.*rp / rn     ! k-eps e1 = 0
       e2 =  1.5 + 1.*rm / rn     ! k-eps e2 = 0
       e3 = -1.0 / rn             ! k-eps e3 = 1
       !--------------------------------------------------
       Call stab_func(sfunc_opt,c1,c2,c3,c4,c5,c1T,c2T,c3T,c4T,cT)
       !--------------------------------------------------
       ! From my Mathematica derivation
       ! Case 9 equations

       sf_n0 = (4.-4.*c2+3.*c4) / (6.*c1)
       sf_n1 = ( (1.-c3) * (2.*c1T*(c4-c5) - c1*(2.-2.*c2T-c4T)) ) / (3.*c1**2*c1T**2)
       sf_n2 = - c4T*(4.-4.*c2+3.*c4)*(2.-2.*c2T-c4T) / (24*c1*c1T**2)
       sf_n3 = ( (1.-c3)*(1.-c3T)*(2.*c1T*(4.-4.*c2+3.*c5) + 3.*c1*(2.-2.*c2T-c4T)) ) / (6.*c1**2*c1T**2)

       sf_n0T = 2. / (3.*c1T)                   ! Downgradient part, prop to dThetadz
       sf_n1T = (2.*(1.-c3)) / (3.*c1*c1T**2)   ! Downgradient part, prop to dThetadz
       sf_n2T = ( c1*c4T*(4.-4.*c2+3.*c4) + 8*c5*c1T*(1.-c2+c5) - 2.*c4*c1T*(2.-2.*c2+3.*c5) ) / (12.*c1**2*c1T**2)

       sf_n0T_star = (1.-c3T) / (c1T)                     ! "Possibly countergradient" part, prop to Theta^2 * k*g*alpha/eps
       sf_n1T_star = ((1.-c3)*(1.-c3T)) / (c1*c1T**2)     ! "Possibly countergradient" part, prop to Theta^2 * k*g*alpha/eps
       sf_n2T_star = ( (1.-c3T) * (3.*c5**2 + 6.*c5*(1.-c2) + 2.*(1.-c2)**2) ) / (3.*c1**2*c1T)

       sf_d0 = 1.
       sf_d1 = (7.*(1.-c3)) / (3.*c1*c1T)
       sf_d2 = (3.*c5**2 + 6.*c5*(1.-c2) + 2.*(1.-c2)**2) / (3.*c1**2)  -  (c4T*(2.-2.*c2T-c4T)) / (4.*c1T**2)
       sf_d3 = ( (1.-c3) * (3.*c1*c4T*(1.-c2+c5) + c5*c1T*(2.-2.*c2+c5) - c1*(1.-c2T)*(2.-2.*c2+3.*c5)) )/(3.*c1**3*c1T**2)
       sf_d4 = (4.*(1.-c3)**2) / (3.*c1**2*c1T**2)
       sf_d5 = ( -c4T*(2.-2.*c2T-c4T)*(3.*c5**2 + 6.*c5*(1.-c2) + 2.*(1.-c2)**2) ) / (12.*c1**2*c1T**2)


       lim_am0 = sf_d0*sf_n0
       lim_am1 = sf_d0*sf_n1 + sf_d1*sf_n0
       lim_am2 = sf_d1*sf_n1 + sf_d4*sf_n0
       lim_am3 = sf_d4*sf_n1
       lim_am4 = sf_d2*sf_n0
       lim_am5 = sf_d2*sf_n1+sf_d3*sf_n0
       lim_am6 = sf_d3*sf_n1
       !--------------------------------------------------
       ! Initialization of various constants
       cm0     =  (  (2 * (sf_d5 - sf_n2))  &
                     /  ( -1 * (sf_d2 - sf_n0) - sqrt( (sf_d2 - sf_n0)**2 - 4*sf_d0 * (sf_d5 - sf_n2) ) )  )**0.25  ! Compute cmu0


       cm0inv2 = 1./cm0**2                                               ! inverse of cmu0 squared
       ! minimum value of alpha_n (corresponding to the free convective case at equilibrium, G = epsilon)
       alpha_n_min = 0.5*( - ( sf_d1 + sf_n0T )  + sqrt(  ( sf_d1 + sf_n0T )**2     &
                - 4. * sf_d0 *( sf_d4 + sf_n1T ) ) ) / ( sf_d4 + sf_n1T )
       cff     = (cm0**3 )*(tke_min**1.5) / eps_min                      ! Compute gls_min consistently
       gls_min = (cm0**rp)*(tke_min**rm ) * ( cff**rn )                  !  with eps_min/tke_min

       trb_min(itke) = tke_min
       trb_min(igls) = gls_min

       !--------------------------------------------------
       ! Compute the vertical shear
       tind = nstp
       DO k=1,N-1
          cff = 2. / ( Hz(k ) + Hz(k+1 ) )
          du  = cff*( u(k+1, tind)-u(k, tind) )
          dv  = cff*( v(k+1, tind)-v(k, tind) )
          shear2(k) = du*du + dv*dv
       ENDDO

       !--------------------------------------------------
       ! Compute ustar squared at the surface and at the bottom
       ustar_sfc_sq = sqrt( sustr**2+svstr**2 )
       ustar_bot_sq = r_D * sqrt( u(1,tind)**2 + v(1,tind)**2  )
       ! Compute the dissipation rate
       DO k=1,N-1
          cff       = (cm0**e1) * ( trb( k,nstp,itke )**e2 )  &
                                * ( trb( k,nstp,igls )**e3 )
          diss(k)   = MAX( cff , eps_min )
       ENDDO

       !--------------------------------------------------
       DO ig = 1,ngls     ! ig = 2 for gls and = 1 for tke
       !--------------------------------------------------
          ! Off-diagonal terms for the tridiagonal problem
          cff=-0.5*dt
          DO k=2,N-1
             FC(k) = cff*OneOverSig(ig)*( Akv(k)+Akv(k-1) ) / Hz(k)
          ENDDO

          IF(Neu_bot) THEN
             FC(1) = 0.
          ELSE
             FC(1) = cff*OneOverSig(ig)*( Akv(1)+Akv(0) ) / Hz(1)
          END IF

          FC(N)=0.    ! The Neumann condition is added afterwards via a term in RH(N-1)
          ! Production/Dissipation terms and diagonal term
          DO k=1,N-1
             ig1   = (igls-ig); ig2 = (ig-itke)  ! tke: (ig1 = 1,ig2 = 0) ; gls: (ig1 = 0,ig2 = 1)
             invk  =     1. / trb( k,nstp,itke ) ! 1/tke
             gls   =          trb( k,nstp,igls )
             invG  =  ig1*invk+ig2*(1./gls)          ! invG = 1/tke for tke ; invG = 1/psi for gls     ! NEW : added the missing division by tke : ig1+ig2*(1./gls) => ig1*invk+ig2*(1./gls)
             cff1  =  ig1+ig2*beta1   * invk*gls     ! Coeff for P (shear production)
             cff2  = (ig1+ig2*beta2 ) * invk         ! Coeff for the dissipation
             cff3m =  ig1+ig2*beta3m  * invk*gls     ! Coeff for G (buoyancy). Ex for k-eps : beta3m = -0.4
             cff3p =  ig1+ig2*beta3p  * invk*gls     ! Coeff for G (buoyancy).                beta3p = 1. for all the schemes
             ! Shear and buoyancy production
             Sprod =  cff1*Akv(k) * shear2(k)
             Bprod = cff3p * (-Akt(k,1)*bvf(k)  &
                              +trb( k,nstp,itke )/diss(k) * (g*alpha)**2 * c_mu_prime_star(k) * varT(k,nstp))
             ! for the tke equation: cff3p = 1, Bprod can be >0 or <0
             ! for the gls equation: cff3p = beta3p * psi/tke, Bprod can be >0 or <0
             ! Patankar trick to ensure non-negative solutions
             cff   =       0.5*(Hz(k)+Hz(k+1))
             IF( (Bprod + Sprod) .gt. 0.) THEN
                RH(k) = cff*( trb(k,nstp,ig) + dt*(Bprod+Sprod) )     ! NEW : Correction of the mistake : trb(k,nnew,ig) => trb(k,nstp,ig)
                DC(k) = cff*(1.+dt*cff2*diss(k))-FC(k)-FC(k+1)
             ELSE
                RH(k) = cff*( trb(k,nstp,ig) + dt*       Sprod  )     ! Sprod > 0, Bprod < 0; using Patankar trick for treating the Bprod contribution
                DC(k) = cff*(1.+dt*(cff2*diss(k)                    &
                                  -invG*Bprod)) - FC(k) - FC(k+1)     ! NEW : Correction of the mistake : trb(k,nnew,ig) => trb(k,nstp,ig)
             ENDIF
          ENDDO

          ! Boundary conditions
          IF( ig == itke ) THEN
             ! surface
             trb_sfc      = MAX( tke_min, cm0inv2*ustar_sfc_sq )
             flux_top     = 0.
             ! bottom
             trb_bot      = MAX( tke_min, cm0inv2*ustar_bot_sq )
             flux_bot     = 0.
             ! finalize
             IF(Neu_bot) THEN
                RH(1   ) = RH(  1) + dt*flux_bot
             ELSE
                RH(1   ) = RH(  1) - FC(1)*trb_bot
             ENDIF
             RH(N-1 ) = RH(N-1) + dt*flux_top
             trb(N,nnew,ig ) = trb_sfc
             trb(0,nnew,ig ) = trb_bot
          ELSE
             ! surface
             z0_s = MAX( 1.e-2 , chk*ustar_sfc_sq )   !<-- Charnock
!                  cff     = 30.*tanh( 0.6 / (28.*sqrt( ustar_sfc_sq(i,j) )) )
!                  z0_s    = MAX( 1.e-2 ,
!     &                       1.3*( 782.353/g )*ustar_sfc_sq(i,j)*(cff**1.5) )
              cff = 0.5*( trb(N-1,nnew,itke )+trb( N  ,nnew,itke ) )
              lgthsc      = vonKar*(0.5*Hz(N)+z0_s)
              trb_sfc     = MAX(gls_min,(cm0**rp)*(lgthsc**rn)*(cff**rm))
              flux_top    = -rn*cm0**(rp+1.)*vonKar*OneOverSig(igls)  &
                                       *(cff**(rm+0.5))*(lgthsc**rn)
              ! bottom
              z0_b        = MAX( Zob , 1.E-04 )
              cff         = 0.5*( trb(1,nnew,itke ) + trb(0,nnew,itke ) )
              lgthsc      = vonKar*(0.5*Hz(1)+z0_b)
              trb_bot     = MAX(gls_min,(cm0**rp)*(lgthsc**rn)*(cff**rm))
              flux_bot    =-rn*cm0**(rp+1.)                &
                                *vonKar*OneOverSig(igls)   &
                                *(cff**(rm+0.5))*(lgthsc**rn)

              IF( ustar_bot_sq == 0. ) THEN
                 flux_bot = 0.
                 trb_bot  = gls_min
              ENDIF
              ! finalize
              IF(Neu_bot) THEN
                 RH(  1  ) = RH(  1) + dt*flux_bot
              ELSE
                 RH(  1  ) = RH(  1) - FC(1)*trb_bot
              END IF
              RH( N-1 ) = RH(N-1) + dt*flux_top
              trb( N,nnew,ig ) = trb_sfc
              trb( 0,nnew,ig ) = trb_bot
            ENDIF

            ! tridiagonal resolution
            cff       =  1./DC(N-1)
            CF(N-1) = cff*FC(N-1)
            RH(N-1) = cff*RH(N-1)

            DO k=N-2,1,-1
               cff   =   1./(DC(k)-CF(k+1)*FC(k+1))
               CF(k) = cff*FC(k)
               RH(k) = cff*( RH(k)-FC(k+1)*RH(k+1))
            ENDDO

            trb(1,nnew,ig ) = MAX( RH(1), trb_min(ig) )

            DO k=2,N-1
               RH(k) = RH(k)-CF(k)*RH(k-1)
               trb(k,nnew,ig ) = MAX( RH(k), trb_min(ig) )
            ENDDO
         !--------------------------------------------------
         ENDDO     ! ig loop
         !--------------------------------------------------

         DO k=1,N-1
            !
            ! Galperin limitation : l <= l_lim
            L_lim = galp * sqrt( 2.* trb(k,nnew,itke)) /        &
                                  ( sqrt(max(eps_bvf, bvf(k)))  )
            !
            ! Limitation on psi (use MAX because rn is negative)
            cff = (cm0**rp) * (L_lim**rn) * (trb(k,nnew,itke)**rm)
            trb( k,nnew,igls ) = MAX( trb( k,nnew,igls ),cff )
            !
            ! Dissipation rate
            epsilon = (cm0**e1) * ( trb(k,nnew,itke )**e2 )   &
                                * ( trb(k,nnew,igls )**e3 )
            epsilon = MAX(epsilon,eps_min)
            eps(k) = epsilon
            !
            ! Compute alpha_n and alpha_m
            cff     = ( trb(k,nnew,itke)/epsilon )**2
            alpha_m(k)     = cff*  shear2(k)    ! Caution, mix of previous time step (shear2) and new time step (cff)
            alpha_n(k)     = cff*     bvf(k)    ! Same
            alpha_v(k)     = trb(k,nnew,itke) / epsilon**2 * (g * alpha)**2 * varT(k,nstp)    ! Same
            !
            ! Limitation of alpha_n and alpha_m
            ! alpha_n_min    = -10000.
            alpha_n(k)     = MIN(  MAX( 0.73*alpha_n_min , alpha_n(k) ) , 1.0e10 )



            alpha_m_max = ( lim_am0 + lim_am1 * alpha_n(k)                         &
                           +  lim_am2 * alpha_n(k)**2 + lim_am3 * alpha_n(k)**3) /   &
                      ( lim_am4 + lim_am5 * alpha_n(k) + lim_am6 * alpha_n(k)**2)


            ! alpha_m_max = 10000.
            alpha_m(k) = MIN(alpha_m(k) , alpha_m_max)
            !
            ! Compute stability functions
            denom = sf_d0  +  sf_d1*alpha_n(k) +  sf_d2*alpha_m(k)   &
                 + sf_d3*alpha_n(k)*alpha_m(k) + sf_d4*alpha_n(k)**2 + sf_d5*alpha_m(k)**2
            c_mu(k)      = MAX( (sf_n0  +  sf_n1*alpha_n(k) +  sf_n2*alpha_m(k) + sf_n3*alpha_v(k))/denom   , 0.)
            c_mu_prime(k) = MAX( (sf_n0T + sf_n1T*alpha_n(k) + sf_n2T*alpha_m(k))/denom   , 0.)
            c_mu_prime_star(k) = MAX( (sf_n0T_star + sf_n1T_star*alpha_n(k) + sf_n2T_star*alpha_m(k))/denom   , 0.)
            !

            ! Finalize the computation of Akv and Akt
            cff1 = trb( k,nnew,itke )**2 / epsilon
            cff2 = trb( k,nnew,itke ) / epsilon * (g*alpha) * varT(k,nstp)
            Akv(k  )= MAX( cff1*c_mu(k)     ,nuwm )
            Akt(k,1)= MAX( cff1*c_mu_prime(k),nuws )
            gamma_GLS(k,1)= cff2*c_mu_prime_star(k)

            Akt(k,2:ntra)= Akt(k,1)
            gamma_GLS(k,2:ntra)= 0. ! for the moment, we ignore the salinity effect
            !sh Akt(k,2)= Akt(k,1)
            lmix( k ) =  cm0**3 * cff / sqrt( trb( k,nnew,itke ) )
         ENDDO

         Akv(0)   = Akv(1)    ! Needed for the temperature variance equation
         Akt(0,1) = Akt(1,1)  ! Needed for applying the bottom boundary condition of the temperature
         Akt(0,2:ntra) = Akt(0,1)

         if (kt == 0) then
            ! print*,'sf_d0            = ',sf_d0
            ! print*,'sf_d1            = ',sf_d1
            ! print*,'sf_d2            = ',sf_d2
            ! print*,'sf_d3            = ',sf_d3
            ! print*,'sf_d4            = ',sf_d4
            ! print*,'sf_d5            = ',sf_d5
            ! print*,'sf_n0            = ',sf_n0
            ! print*,'sf_n1            = ',sf_n1
            ! print*,'sf_n2            = ',sf_n2
            ! print*,'sf_n3            = ',sf_n3
            ! print*,'sf_n0T           = ',sf_n0T
            ! print*,'sf_n1T           = ',sf_n1T
            ! print*,'sf_n2T           = ',sf_n2T
            ! print*,'sf_n0T_star     = ',sf_n0T_star
            ! print*,'sf_n1T_star     = ',sf_n1T_star
            ! print*,'sf_n2T_star     = ',sf_n2T_star
            print*,'cm0              = ',cm0
            print*,'alpha_n_min      = ',alpha_n_min
            print*,'alpha_m_max      = ',alpha_m_max
            print*,'keps_DGNG'
         ENDIF

         return

end subroutine gls_stp


subroutine stab_func(sfunc_opt,c1,c2,c3,c4,c5,c1T,c2T,c3T,c4T,cT)
implicit none
integer,intent(in)    :: sfunc_opt
real(8),intent(out)   :: c1,c2,c3,c4,c5,c1T,c2T,c3T,c4T,cT
SELECT CASE(sfunc_opt)

CASE(1)
c1=0.;c2=0.;c3=0.;c4=0.;c5=0.
c1T=0.;c2T=0.;c3T=0.;c4T=0.;cT=0.
CASE DEFAULT ! Canuto A
c1=2.5;c2=0.984;c3=0.5;c4=0.512;c5=0.416
c1T=5.95;c2T=0.6;c3T=0.333333;c4T=0.4;cT=1.44
END SELECT
end subroutine stab_func

end module scm_gls_DGNG
