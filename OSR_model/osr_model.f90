!  Idea :
!    - endurance surface beta(sigma, alpha) = 0 (Drucker-Prager)
!    - cinematic translation alpha (back-stress déviatorique)
!    - scalar damage variable D, fracture at D = 1
!    - incremental criterion beta >= 0 AND dbeta > 0  => evolution of alpha and D
!
!  mains equations :
!    (1)  beta     = (sigma_eff + A*I1 - sigma_oe) / sigma_oe
!    (2)  sigma_eff = sqrt(3/2 (s-alpha):(s-alpha))
!    (6)  d alpha = dbeta * C * (s - alpha)
!    (7)  dD/dbeta = K * exp(L * beta)
!    (33) dbeta = [ 3/2 (s-alpha):ds / sigma_eff + A*tr(dsigma) ]
!                 / ( sigma_oe + C*sigma_eff )

module osr_model_mod
   implicit none
   private

   integer, parameter, public :: dp = kind(1.0d0)

   ! material parameters (5 parameters in the OSR model)

   type, public :: osr_params
      real(dp) :: A         ! hydrostatic pressure sensibility
      real(dp) :: sigma_oe  ! endurance limite for null average [MPa]
      real(dp) :: C         ! back-stress evolution rate
      real(dp) :: K         ! prefactor damage rate
      real(dp) :: L         ! exponential factor 
   end type osr_params

   ! public subroutines 
   public :: osr_beta          ! evaluate beta and sigma_eff 
   public :: osr_dbeta_loading ! time derivative of beta during a loading step 
   public :: osr_rhs           ! right-hand side of the ODE system 
   public :: osr_step_rk4      ! one fixed-step RK4 update 
   public :: osr_integrate_rk4 ! full integration with rupture detection 
   public :: osr_step_dopri 
   public :: osr_integrate_dopri 
   public :: osr_active_state

   ! Dormand-Prince coefficents 
   real(dp), parameter :: C2 = 1.0_dp/5.0_dp
   real(dp), parameter :: C3 = 3.0_dp/10.0_dp
   real(dp), parameter :: C4 = 4.0_dp/5.0_dp
   real(dp), parameter :: C5 = 8.0_dp/9.0_dp
   real(dp), parameter :: A21 = 1.0_dp/5.0_dp
   real(dp), parameter :: A31 = 3.0_dp/40.0_dp
   real(dp), parameter :: A32 = 9.0_dp/40.0_dp
   real(dp), parameter :: A41 = 44.0_dp/45.0_dp
   real(dp), parameter :: A42 = -56.0_dp/15.0_dp
   real(dp), parameter :: A43 = 32.0_dp/9.0_dp
   real(dp), parameter :: A51 = 19372.0_dp/6561.0_dp
   real(dp), parameter :: A52 = -25360.0_dp/2187.0_dp
   real(dp), parameter :: A53 = 64448.0_dp/6561.0_dp
   real(dp), parameter :: A54 = -212.0_dp/729.0_dp
   real(dp), parameter :: A61 = 9017.0_dp/3168.0_dp
   real(dp), parameter :: A62 = -355.0_dp/33.0_dp
   real(dp), parameter :: A63 = 46732.0_dp/5247.0_dp
   real(dp), parameter :: A64 = 49.0_dp/176.0_dp
   real(dp), parameter :: A65 = -5103.0_dp/18656.0_dp

   ! weight for the order five solution 
   real(dp), parameter :: B1 = 35.0_dp/384.0_dp
   real(dp), parameter :: B3 = 500.0_dp/1113.0_dp
   real(dp), parameter :: B4 = 125.0_dp/192.0_dp
   real(dp), parameter :: B5 = -2187.0_dp/6784.0_dp
   real(dp), parameter :: B6 = 11.0_dp/84.0_dp

   ! error coefficient e_i = b_i (order 5) - b_hat_i (order 4)
   real(dp), parameter :: E1 = 71.0_dp/57600.0_dp
   real(dp), parameter :: E3 = -71.0_dp/16695.0_dp
   real(dp), parameter :: E4 = 71.0_dp/1920.0_dp
   real(dp), parameter :: E5 = -17253.0_dp/339200.0_dp
   real(dp), parameter :: E6 = 22.0_dp/525.0_dp
   real(dp), parameter :: E7 = -1.0_dp/40.0_dp

   ! largest argument passed to exp() in the damage law
   ! exp overflows above ~709, and with L = 14.4 that is reached for beta ~ 49.
   ! Clamping keeps dD/dt finite (and the failure detection meaningful) instead of
   ! letting an Inf propagate into alpha and D through the RK stages.
   real(dp), parameter :: EXP_ARG_MAX = 300.0_dp

contains

   ! split a stress tensor into its trace I1 and deviator s
   pure subroutine deviatoric(sigma, s, I1)
      real(dp), intent(in)  :: sigma(3,3)
      real(dp), intent(out) :: s(3,3)
      real(dp), intent(out) :: I1
      integer :: i

      I1 = sigma(1,1) + sigma(2,2) + sigma(3,3)
      s = sigma
      do i = 1, 3
         s(i,i) = s(i,i) - I1 / 3.0_dp
      end do
   end subroutine deviatoric

   ! Double contraction A : B for 3x3 tensors 
   pure function ddot(A, B) result(res)
      real(dp), intent(in) :: A(3,3), B(3,3)
      real(dp) :: res
      res = sum(A * B)
   end function ddot

   ! evaluate the endurance function beta(sigma, alpha)
   ! Return also sigma_eff et d = s - alpha (useful elsewhere).
   pure subroutine osr_beta(p, sigma, alpha, beta, sigma_eff, d)
      type(osr_params), intent(in)  :: p
      real(dp),         intent(in)  :: sigma(3,3), alpha(3,3)
      real(dp),         intent(out) :: beta, sigma_eff
      real(dp),         intent(out) :: d(3,3)
      real(dp)                      :: s(3,3), I1

      call deviatoric(sigma, s, I1)
      d = s - alpha
      sigma_eff = sqrt(1.5_dp * ddot(d, d))
      beta = (sigma_eff + p%A * I1 - p%sigma_oe) / p%sigma_oe
   end subroutine osr_beta

   ! beta temporal derived during the loading (eq. 33)
   ! this expresion already inject alpha law evolution,
   ! his sign is the one of dbeta during loading phase.
   pure subroutine osr_dbeta_loading(p, sigma, alpha, dsigma, dbeta, beta, sigma_eff, d)
      type(osr_params), intent(in)  :: p
      real(dp),         intent(in)  :: sigma(3,3), alpha(3,3), dsigma(3,3)
      real(dp),         intent(out) :: dbeta, beta, sigma_eff
      real(dp),         intent(out) :: d(3,3)
      real(dp)                      :: ds(3,3), dI1, num, den

      call osr_beta(p, sigma, alpha, beta, sigma_eff, d)

      ! Deviator and trace of dsigma
      dI1 = dsigma(1,1) + dsigma(2,2) + dsigma(3,3)
      ds  = dsigma
      ds(1,1) = ds(1,1) - dI1 / 3.0_dp
      ds(2,2) = ds(2,2) - dI1 / 3.0_dp
      ds(3,3) = ds(3,3) - dI1 / 3.0_dp

      ! guard against a null effective stress(sigma_eff appears in the denominator of the numerator). This only happens at the very start (alpha=0)
      if (sigma_eff < 1.0e-12_dp) then
         dbeta = 0.0_dp
         return
      end if

      num   = 1.5_dp * ddot(d, ds) / sigma_eff + p%A * dI1
      den   = p%sigma_oe + p%C * sigma_eff
      dbeta = num / den
   end subroutine osr_dbeta_loading
  
   ! rhs : (alpha, D) -> (dalpha/dt, dD/dt)
   pure subroutine osr_rhs(p, sigma, alpha, D_dummy, dsigma, dalpha_dt, dD_dt)
      type(osr_params), intent(in)  :: p
      real(dp),         intent(in)  :: sigma(3,3), alpha(3,3), dsigma(3,3)
      real(dp),         intent(in)  :: D_dummy   ! D non utilisé ici (g indep. de D)
      real(dp),         intent(out) :: dalpha_dt(3,3), dD_dt
      real(dp) :: dbeta, beta, sigma_eff, d(3,3)

      ! Note : D_dummy kept to prepare a potential extension 
      if (D_dummy >= 0.0_dp) continue   ! suppress warning unused

      call osr_dbeta_loading(p, sigma, alpha, dsigma, dbeta, beta, sigma_eff, d)

      ! loading criterion (eq. 5)
      if (beta < 0.0_dp .or. dbeta <= 0.0_dp) then
         dalpha_dt = 0.0_dp
         dD_dt     = 0.0_dp
         return 
      end if

      dalpha_dt = dbeta * p%C * d
      dD_dt     = dbeta * p%K * exp(min(p%L * beta, EXP_ARG_MAX))
   end subroutine osr_rhs
  
   ! a classic RK4 step 
   ! callbacks sigma_fn(t) and dsigma_fn(t) resend 3x3 tensors (at time t)
   ! this callback design permit the same integrator handle very different loading histories without changing any line 
   subroutine osr_step_rk4(p, sigma_fn, dsigma_fn, t, dt, alpha, D)
      type(osr_params), intent(in)    :: p
      real(dp),         intent(in)    :: t, dt
      real(dp),         intent(inout) :: alpha(3,3), D ! because enter with value at time t and finish with value at t+dt

      ! any function matching this signature can be passed in 
      interface
         function sigma_fn(tloc) result(s)
            import :: dp
            real(dp), intent(in) :: tloc
            real(dp) :: s(3,3)
         end function sigma_fn
         function dsigma_fn(tloc) result(s)
            import :: dp
            real(dp), intent(in) :: tloc
            real(dp) :: s(3,3)
         end function dsigma_fn 
      end interface 

      real(dp) :: k1a(3,3), k2a(3,3), k3a(3,3), k4a(3,3)
      real(dp) :: k1D, k2D, k3D, k4D
      real(dp) :: th, t_end

      th    = t + 0.5_dp * dt ! midpoint 
      t_end = t + dt ! end of the step 

      ! RK4 
      call osr_rhs(p, sigma_fn(t), alpha, D, dsigma_fn(t), k1a, k1D)
      call osr_rhs(p, sigma_fn(th), alpha + 0.5_dp*dt*k1a, D + 0.5_dp*dt*k1D, dsigma_fn(th), k2a, k2D)
      call osr_rhs(p, sigma_fn(th), alpha + 0.5_dp*dt*k2a, D + 0.5_dp*dt*k2D, dsigma_fn(th), k3a, k3D)
      call osr_rhs(p, sigma_fn(t_end), alpha + dt*k3a, D + dt*k3D, dsigma_fn(t_end), k4a, k4D)

      alpha = alpha + (dt / 6.0_dp) * (k1a + 2.0_dp*k2a + 2.0_dp*k3a + k4a)
      D     = D     + (dt / 6.0_dp) * (k1D + 2.0_dp*k2D + 2.0_dp*k3D + k4D)
   end subroutine osr_step_rk4

   ! fixed step RK4 integration with fracture detection (D = 1)
   ! integrates from t0 to tf with n_steps uniform steps
   ! linear interpolation to estimate t_failure
   subroutine osr_integrate_rk4(p, sigma_fn, dsigma_fn, t0, tf, n_steps, alpha, D, t_failure, failed)
      type(osr_params), intent(in)    :: p
      real(dp),         intent(in)    :: t0, tf
      integer,          intent(in)    :: n_steps
      real(dp),         intent(inout) :: alpha(3,3), D
      real(dp),         intent(out)   :: t_failure
      logical,          intent(out)   :: failed

      interface
         function sigma_fn(tloc) result(s)
            import :: dp
            real(dp), intent(in) :: tloc
            real(dp) :: s(3,3)
         end function sigma_fn
         function dsigma_fn(tloc) result(s)
            import :: dp
            real(dp), intent(in) :: tloc
            real(dp) :: s(3,3)
         end function dsigma_fn
      end interface

      real(dp) :: t, dt, D_prev
      integer  :: i 

      dt        = (tf - t0) / real(n_steps, dp)
      t         = t0
      failed    = .false.
      t_failure = -1.0_dp

      do i = 1, n_steps
         D_prev = D
         call osr_step_rk4(p, sigma_fn, dsigma_fn, t, dt, alpha, D)
         t = t + dt

         if (D >= 1.0_dp) then
            if (D > D_prev) then
               ! linear interp. to locate the crossing 
               t_failure = t - dt + dt * (1.0_dp - D_prev) / (D - D_prev)
            else
               t_failure = t
            end if
            D= 1.0_dp
            failed = .true.
            return
         end if
      end do
   end subroutine osr_integrate_rk4

   ! implementation of DOPRI scheme (useful for complex loading history)
   subroutine osr_step_dopri(p, sigma_fn, dsigma_fn, t, h, alpha, D, alpha_new, D_new, err_norm, atol, rtol)
      type(osr_params), intent(in) :: p
      real(dp), intent(in)         :: t, h, atol, rtol 
      real(dp), intent(in)         :: alpha(3,3), D 
      real(dp), intent(out)        :: alpha_new(3,3), D_new, err_norm 

      interface 
         function sigma_fn(t_loc) result(s)
            import :: dp 
            real(dp), intent(in) :: t_loc 
            real(dp)             :: s(3,3)
         end function 
         function dsigma_fn(tloc) result(s)
            import :: dp 
            real(dp), intent(in) :: tloc 
            real(dp) :: s(3,3)
         end function 
      end interface 

      real(dp) :: k1a(3,3), k2a(3,3), k3a(3,3), k4a(3,3), k5a(3,3), k6a(3,3), k7a(3,3), k1D, k2D, k3D, k4D, k5D, k6D, k7D
      real(dp) :: a_stage(3,3), D_stage, err_alpha(3,3), err_D , sc_A(3,3), sc_D, sum_sq
      integer  :: i, j 

      ! stage 1 
      call osr_rhs(p, sigma_fn(t), alpha, D, dsigma_fn(t), k1a, k1D)

      ! stage 2 
      a_stage = alpha + h*(A21*k1a) ; D_stage = D + h*(A21*k1D)
      call osr_rhs(p, sigma_fn(t+C2*h), a_stage, D_stage, dsigma_fn(t+C2*h), k2a, k2D)

      ! stage 3 
      a_stage = alpha + h*(A31*k1a + A32*k2a)
      D_stage = D + h*(A31*k1D + A32*k2D)
      call osr_rhs(p, sigma_fn(t+C3*h), a_stage, D_stage, dsigma_fn(t+C3*h), k3a, k3D)

      ! stage 4 
      a_stage = alpha + h*(A41*k1a + A42*k2a + A43*k3a)
      D_stage = D + h*(A41*k1D + A42*k2D + A43*k3D)
      call osr_rhs(p, sigma_fn(t+C4*h), a_stage, D_stage, dsigma_fn(t+C4*h), k4a, k4D)

      ! stage 5 
      a_stage = alpha +h*(A51*k1a + A52*k2a + A53*k3a + A54*k4a)
      D_stage = D + h*(A51*k1D + A52*k2D + A53*k3D + A54*k4D)
      call osr_rhs(p, sigma_fn(t+C5*h), a_stage, D_stage, dsigma_fn(t+C5*h), k5a, k5D)

      ! stage 6 
      a_stage = alpha + h*(A61*k1a + A62*k2a + A63*k3a + A64*k4a +A65*k5a)
      D_stage = D + h*(A61*k1D + A62*k2D + A63*k3D + A64*k4D +A65*k5D)
      call osr_rhs(p, sigma_fn(t+h), a_stage, D_stage, dsigma_fn(t+h), k6a, k6D)

      ! order-5 solution 
      alpha_new = alpha + h*(B1*k1a + B3*k3a + B4*k4a + B5*k5a + B6*k6a)
      D_new     = D + h*(B1*k1D + B3*k3D + B4*k4D + B5*k5D + B6*k6D)

      ! stage 7 for the error estimate 
      call osr_rhs(p, sigma_fn(t+h), alpha_new, D_new, dsigma_fn(t+h), k7a, k7D)
      
      ! local error estimate : E = y5 - y_hat = h*sum(e_i*k_i)
      err_alpha = h*(E1*k1a + E3*k3a + E4*k4a + E5*k5a + E6*k6a + E7*k7a)
      err_D     = h*(E1*k1D + E3*k3D + E4*k4D + E5*k5D + E6*k6D + E7*k7D)

      ! Per-component tolerance scales 
      sum_sq=0.0_dp 

      do j = 1,3 
         do i=1,3 
            sc_a(i,j) = atol + rtol*max(abs(alpha(i,j)), abs(alpha_new(i,j)))
         end do 
      end do 
      sc_D = atol + rtol*max(abs(D), abs(D_new))

      ! root-mean-square normalised error
      do j = 1,3 
         do i = 1,3 
            sum_sq = sum_sq + (err_alpha(i,j)/sc_a(i,j))**2
         end do 
      end do 
      sum_sq = sum_sq + (err_D/sc_D)**2
      err_norm = sqrt(sum_sq/10.0_dp)

   end subroutine osr_step_dopri

   subroutine osr_integrate_dopri(p, sigma_fn, dsigma_fn, t0, tf, n_initial_steps, alpha, D, t_failure, failed, use_events, n_steps_used, &
                                  rtol_in, atol_in, n_attempts)

      type(osr_params), intent(in) :: p
      real(dp), intent(in)         :: t0, tf 
      integer, intent(in)          :: n_initial_steps
      real(dp), intent(inout)      :: alpha(3,3), D 
      real(dp), intent(out)        :: t_failure 
      logical, intent(out)         :: failed 
      logical, intent(in), optional  :: use_events 
      integer, intent(out), optional :: n_steps_used, n_attempts
      real(dp), intent(in), optional :: atol_in, rtol_in
      
      interface
         function sigma_fn(tloc) result(s)
            import :: dp
            real(dp), intent(in) :: tloc
            real(dp) :: s(3,3)
         end function sigma_fn
         function dsigma_fn(tloc) result(s)
            import :: dp
            real(dp), intent(in) :: tloc
            real(dp) :: s(3,3)
         end function dsigma_fn
      end interface

      real(dp), parameter :: SAFETY = 0.9_dp, FACMIN = 0.1_dp, FACMAX = 5.0_dp, EXPON = -0.2_dp 
      integer, parameter  :: MAX_STEPS = 10000000
      
      real(dp) :: t, h, hmin, hmax, err_norm, factor, ATOL_use, RTOL_use 
      real(dp) :: alpha_new(3,3), D_new, D_prev 
      real(dp) :: beta_start, beta_end, dbeta_start, dbeta_end
      logical  :: active_start, active_end
      integer  :: nstep, naccept
      logical  :: events_on 

      events_on = .false. 
      if (present(use_events)) events_on = use_events 

      t= t0
      h = (tf-t0)/real(max(n_initial_steps,1),dp)
      hmin = abs(h)*1.0e-8_dp 
      hmax = 10.0_dp*abs(h) 
      failed = .false. 
      t_failure = -1.0_dp 
      nstep = 0 
      naccept = 0

      ATOL_use = 1.0e-10_dp 
      RTOL_use = 1.0e-6_dp 
      if (present(atol_in)) ATOL_use = atol_in 
      if (present(rtol_in)) RTOL_use = rtol_in 

      do while (t<tf .and. nstep<MAX_STEPS) 
         nstep = nstep + 1 
         if (t+h>tf) h = tf-t 

         D_prev = D 
         call osr_step_dopri(p, sigma_fn, dsigma_fn, t, h, alpha, D, alpha_new, D_new, err_norm, ATOL_use, RTOL_use)

         if (err_norm<1.0_dp) then 
            ! event handling by step rejection 
            if (events_on) then 
               call osr_active_state(p, sigma_fn(t), alpha, dsigma_fn(t), active_start, beta_start, dbeta_start)
               call osr_active_state(p, sigma_fn(t+h), alpha_new, dsigma_fn(t+h), active_end, beta_end, dbeta_end)
               if (active_start .neqv. active_end) then
                  ! print *, "event detected, halving h from ", h, "at t = ", t -> debug 
                  if (h>2.0_dp*hmin) then 
                     h = 0.5_dp*h
                     cycle
                  end if 
                  ! if h is already extremely small, accept the step 
                  ! this avoids stagnation near the switching boundary 
               end if 
            end if 
            ! normal accepted step 
            t = t+h 
            alpha = alpha_new 
            D = D_new 
            naccept = naccept +1 

            if (D>=1.0_dp) then 
               if (D>D_prev) then 
                  t_failure = t - h + h*(1.0_dp-D_prev)/(D-D_prev) ! linear interpolation 
               else 
                  t_failure = t 
               end if 
               D = 1.0_dp 
               failed = .true. 
               if (present(n_steps_used)) n_steps_used = naccept  
               if (present(n_attempts)) n_attempts = nstep
               return 
            end if 
            if (err_norm > 0.0_dp) then 
               factor = SAFETY*err_norm**EXPON ! 0.9*err**0.2
            else 
               factor = FACMAX
            end if
            h = h*min(FACMAX,max(FACMIN,factor))
            if (h>hmax) h=hmax 
         else 
            ! step rejected : shrink the step and retry 
            ! alpha_new and D_new correspond to the old h so they must not be accepted unless the current step was already hmin  
            if (h<=1.001_dp*hmin) then 
               print *, "DOPRI accepting poor step at hmin, t = ", t, &
                        "err_norm = ", err_norm
               t = t + h 
               alpha = alpha_new 
               D = D_new 
               naccept = naccept + 1

               if (D>= 1.0_dp) then 
                  D = 1.0_dp 
                  failed = .true.
                  t_failure = t

                  if (present(n_steps_used)) n_steps_used = naccept
                  if (present(n_attempts)) n_attempts = nstep
                  return 
               end if 

               h = min(hmax, 10.0_dp*hmin)

            else 
               factor = SAFETY*err_norm**EXPON ! 0.9*err_norm^-0.2
               h = h*max(FACMIN,factor)
               if (h<hmin) h = hmin
            end if 
         end if 
      end do 

      if (present(n_steps_used)) n_steps_used = naccept
      if (present(n_attempts)) n_attempts = nstep
   
   end subroutine osr_integrate_dopri
 
   pure subroutine osr_active_state(p, sigma, alpha, dsigma, active, beta, dbeta)
      type(osr_params), intent(in) :: p 
      real(dp), intent(in)         :: sigma(3,3), alpha(3,3), dsigma(3,3)
      logical, intent(out)         :: active 
      real(dp), intent(out)        :: beta, dbeta

      real(dp) :: sigma_eff, d_dev(3,3)

      call osr_dbeta_loading(p, sigma, alpha, dsigma, dbeta, beta, sigma_eff, d_dev)

      active = (beta>=0.0_dp .and. dbeta>0.0_dp)
      
   end subroutine osr_active_state
         
end module osr_model_mod