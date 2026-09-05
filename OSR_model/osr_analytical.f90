!  semi-analytical solution of the OSR mode lfor uniaxial periodic loading 
!  implementvéquations (23), (24) et (29) of the paper 
!  provides a reference solution against which the numerical implementation of the model can be verified
!  the two solutions come from the same set of equations but coputed by independent methods 
!  strategy :
!    1. the stress oscillates between sigma2 = sigma_m + sigma_a (max)
!       and sigma4 = sigma_m - sigma_a (min).
!    2. in steady state, back-stress positions a2 and a4 at these peaks satisfy the coupled non-linear system (23)-(24) -> solve with Newton-Raphson
!    3. the damage per cycle can be find with (29) -> N = 1 / dD

module osr_analytical_mod
   use osr_model_mod, only: dp, osr_params
   implicit none
   private

   public :: solve_periodic_state
   public :: N_analytical

contains

   ! residuals of the periodic state equations (23) and (24)
   ! finds (a2, a4) in such a way that F1 = F2 = 0.
   ! valid = .false. if an argument of the log becomes <= 0.
   subroutine residuals(p, sigma2, sigma4, a2, a4, F1, F2, valid)
      type(osr_params), intent(in)  :: p
      real(dp),         intent(in)  :: sigma2, sigma4, a2, a4
      real(dp),         intent(out) :: F1, F2
      logical,          intent(out) :: valid
      real(dp) :: CA, num1, den1, num2, den2

      CA = p%C * p%A

      ! arguments of the two log terms (must all be strictly positive)
      num1 = 1.0_dp - (CA/p%sigma_oe) * (sigma2 - 1.5_dp*a2)
      den1 = 1.0_dp - (CA/(p%sigma_oe*(p%A+1.0_dp))) &
                       * (p%sigma_oe - 1.5_dp*p%A*a4)
      num2 = 1.0_dp - (CA/p%sigma_oe) * (sigma4 - 1.5_dp*a4)
      den2 = 1.0_dp - (CA/(p%sigma_oe*(p%A-1.0_dp))) &
                       * (p%sigma_oe - 1.5_dp*p%A*a2)

      valid = (num1 > 0.0_dp) .and. (den1 > 0.0_dp) .and. &
              (num2 > 0.0_dp) .and. (den2 > 0.0_dp)
      if (.not. valid) then
         F1 = 0.0_dp ;  F2 = 0.0_dp
         return
      end if

      ! equation (23)
      F1 = 1.5_dp*a2 - (p%A+1.0_dp)*sigma2 + p%sigma_oe &
           - (p%sigma_oe/CA)*(p%A+1.0_dp)*log(num1/den1)
      ! equation (24)
      F2 = -1.5_dp*a4 - (p%A-1.0_dp)*sigma4 + p%sigma_oe &
           - (p%sigma_oe/CA)*(p%A-1.0_dp)*log(num2/den2)
   end subroutine residuals

   ! solve the system (23)-(24) with Newton-Raphson
   ! the jacobian is approximated by finite differences which qvoids the algebraic burden of deriving it analytically
   ! starting point is (0,0) which is always inside the valid domain and close enough to converge for realistic parameters 
   subroutine solve_periodic_state(p, sigma2, sigma4, a2, a4, converged)
      type(osr_params), intent(in)  :: p
      real(dp),         intent(in)  :: sigma2, sigma4
      real(dp),         intent(out) :: a2, a4
      logical,          intent(out) :: converged

      integer,  parameter :: maxiter = 100
      real(dp), parameter :: tol = 1.0e-10_dp
      real(dp), parameter :: h   = 0.05_dp   ! finite difference step 
      real(dp) :: F1, F2, F1p, F2p
      real(dp) :: J11, J12, J21, J22, det, dx1, dx2, nrm, lambda
      logical  :: valid, vtmp
      integer  :: iter, k

      a2 = 0.0_dp ;  a4 = 0.0_dp  ! initial back-stress position are zero 
      converged = .false.

      do iter = 1, maxiter
         call residuals(p, sigma2, sigma4, a2, a4, F1, F2, valid)
         if (.not. valid) return

         nrm = sqrt(F1*F1 + F2*F2)
         if (nrm < tol) then ! converged
            converged = .true.
            return
         end if

         ! numerical jacobian (forward differences)
         call residuals(p, sigma2, sigma4, a2+h, a4, F1p, F2p, vtmp)
         J11 = (F1p - F1)/h ;  J21 = (F2p - F2)/h
         call residuals(p, sigma2, sigma4, a2, a4+h, F1p, F2p, vtmp)
         J12 = (F1p - F1)/h ;  J22 = (F2p - F2)/h

         det = J11*J22 - J12*J21
         if (abs(det) < 1.0e-30_dp) return ! singular jacobian 

         ! dx = -J^{-1} F (2x2 explicit inverse)
         dx1 = -( J22*F1 - J12*F2)/det
         dx2 = -(-J21*F1 + J11*F2)/det

         ! Damped line search : keep halving lambda until the trial 
         ! point (a2+lambda*dx1, a4+lambda*dx2) is inside the log validity domain 
         lambda = 1.0_dp
         do k = 1, 40
            call residuals(p, sigma2, sigma4, a2+lambda*dx1, a4+lambda*dx2, &
                           F1p, F2p, vtmp)
            if (vtmp) exit
            lambda = 0.5_dp*lambda
         end do

         a2 = a2 + lambda*dx1
         a4 = a4 + lambda*dx2
      end do
   end subroutine solve_periodic_state

   ! analytical number of cycles to failure under a uniaxial periodic loading of amplitude sigma_a and mean sigma_m  
   !   infinite_life = .true.  -> loading under the endurance limit 
   !   N < 0                   -> failure of Newton convergence 

   function N_analytical(p, sigma_m, sigma_a, infinite_life) result(N)
      type(osr_params), intent(in)  :: p
      real(dp),         intent(in)  :: sigma_m, sigma_a
      logical,          intent(out) :: infinite_life
      real(dp) :: N
      real(dp) :: sigma2, sigma4, a2, a4, b2, b4, dD
      logical  :: converged

      infinite_life = .false.
      N = huge(1.0_dp)

      ! stress peaks of the cycle 
      sigma2 = sigma_m + sigma_a
      sigma4 = sigma_m - sigma_a

      ! solve for the steady-state back-stress positions
      call solve_periodic_state(p, sigma2, sigma4, a2, a4, converged)
      if (.not. converged) then
         N = -1.0_dp
         return
      end if

      ! beta at both peaks (eq. 15 with j = +1 at peak up, -1 at peak down)
      b2 = ((p%A+1.0_dp)*sigma2 - 1.5_dp*a2 - p%sigma_oe)/p%sigma_oe
      b4 = ((p%A-1.0_dp)*sigma4 + 1.5_dp*a4 - p%sigma_oe)/p%sigma_oe

      ! per cycle damage increment (eq. 28)
      dD = (p%K/p%L)*(exp(p%L*b2) + exp(p%L*b4) - 2.0_dp)

      if (dD <= 0.0_dp) then
         infinite_life = .true.       ! no damage -> infinite life 
         N = huge(1.0_dp)
      else
         N = 1.0_dp/dD                ! eq. (29)
      end if
   end function N_analytical

end module osr_analytical_mod

! Callbacks functions sigma(t) and dsigma/dt(t) for a uniaxial sinusoidal loading
module uniaxial_loading_mod
   use osr_model_mod, only: dp
   implicit none

   real(dp), save :: u_sigma_m = 0.0_dp ! mean stress [MPa]
   real(dp), save :: u_sigma_a = 0.0_dp ! amplitude [MPa]
   real(dp), save :: u_omega   = 0.0_dp ! angular frequency [rad/time]

contains

   ! compatible with the signature of the osr_integrate_rk4 
   function uni_sigma(t) result(s)
      real(dp), intent(in) :: t
      real(dp) :: s(3,3)
      s = 0.0_dp
      s(1,1) = u_sigma_m + u_sigma_a * sin(u_omega * t)
   end function uni_sigma

   function uni_dsigma(t) result(s)
      real(dp), intent(in) :: t
      real(dp) :: s(3,3)
      s = 0.0_dp
      s(1,1) = u_sigma_a * u_omega * cos(u_omega * t)
   end function uni_dsigma

end module uniaxial_loading_mod