!  auxiliary modules for the OSR pipeline :
!    * osr_io_mod              -> reading / writting of the field files 
!    * pipeline_state_mod      -> shared state (current material point, loading parameters, load history)
!    * pipeline_callbacks_mod  -> sigma(t) and dsigma/dt(t) reconstructed by superposition of the units stress fields with the time coefficients
! pipeline_state_mod uses module variables that are set once and read-only during the parallel loop 

module osr_io_mod
   use osr_model_mod, only: dp
   implicit none
   private

   public :: stress_field_t
   public :: read_stress_field
   public :: write_damage_field
   public :: write_damage_csv
   public :: read_stress_field_binary, read_stress_field_auto 

   ! Type which stored a stress field for a loading case 
   type :: stress_field_t
      integer                :: n_points = 0
      integer, allocatable   :: elem_id(:) ! new (parent element ID of each retained node, useful for dedup)
      integer, allocatable   :: node_id(:) ! ID of the node (0 if centroid)
      real(dp), allocatable  :: coords(:,:)    ! (3, n_points)
      real(dp), allocatable  :: sigma(:,:,:)   ! (3, 3, n_points)
   end type stress_field_t

contains

   ! read a ASCII file 
   ! expected format (1 heading ligne begin by '#', puis : node_id   X   Y   Z   SXX SYY SZZ   SXY SYZ SXZ
   ! stress_scale converts the stresses of the file into the MPa expected by the
   ! material parameters (1.0 if the export is already in MPa, 1.0e-6 for an Ansys
   ! model solved in MKS, where the export is in Pa)
   subroutine read_stress_field(filename, sf, stress_scale)
      character(len=*),    intent(in)  :: filename
      type(stress_field_t), intent(out) :: sf
      real(dp),             intent(in)  :: stress_scale

      integer  :: unit, ios, i, n
      real(dp) :: r_eid, r_nid
      real(dp) :: xx, yy, zz, sxx, syy, szz, sxy, syz, sxz
      character(len=1024) :: line

      ! first step : count useful lines 
      open(newunit=unit, file=filename, status='old', action='read', &
           iostat=ios)
      if (ios /= 0) then
         print '(A,A)', "Erreur d'ouverture : ", trim(filename)
         stop 1
      end if

      n = 0
      do
         read(unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0)   cycle
         if (line(1:1) == '#')      cycle
         !debug 
         if (line(1:1) == 'e')      cycle 
         n = n + 1
      end do

      ! Allocation (now we know n)
      sf%n_points = n
      allocate(sf%coords(3, n))
      allocate(sf%sigma (3, 3, n))
      allocate(sf%elem_id(n))
      allocate(sf%node_id(n))

      ! second step : reading
      rewind(unit)
      i = 0
      do
         read(unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0)   cycle
         if (line(1:1) == '#')      cycle
         if (line(1:1) == 'e')      cycle 
         ! ignore if the line not begin with a number 

         i = i + 1
         !first try : nodes format (11 columns)
         read(line, *, iostat=ios) r_eid, r_nid, xx, yy, zz, sxx, syy, szz, sxy, syz, sxz
         if (ios/=0) then 
            ! second try  : centroïd format (10 columns)
            read(line, *, iostat=ios) r_eid, xx, yy, zz, sxx, syy, szz, sxy, syz, sxz
            r_nid = 0.0_dp 
         end if 
         sf%elem_id(i) = nint(r_eid)
         sf%node_id(i) = nint(r_nid)
         sf%coords(:, i)  = [xx, yy, zz]
         ! convert to MPa here, once, so that every consumer of stress_field_t
         ! works in the unit system the material parameters are calibrated for
         sxx = sxx*stress_scale ; syy = syy*stress_scale ; szz = szz*stress_scale
         sxy = sxy*stress_scale ; syz = syz*stress_scale ; sxz = sxz*stress_scale
         sf%sigma(1,1,i) = sxx;  sf%sigma(2,2,i) = syy;  sf%sigma(3,3,i) = szz
         sf%sigma(1,2,i) = sxy;  sf%sigma(2,1,i) = sxy
         sf%sigma(2,3,i) = syz;  sf%sigma(3,2,i) = syz
         sf%sigma(1,3,i) = sxz;  sf%sigma(3,1,i) = sxz
      end do
      close(unit)

      print '(A,I8,A,A)', "  ", n, " points lus depuis ", trim(filename)
   end subroutine read_stress_field

   ! Write the final damage field in a ASCII file (human readable format used for inspection)
   subroutine write_damage_field(filename, sf, D, N_fail)
      character(len=*), intent(in)             :: filename
      type(stress_field_t), intent(in)         :: sf
      real(dp), intent(in)                     :: D(:), N_fail(:)

      integer :: unit, i, n
      n = size(D)

      open(newunit=unit, file=filename, status='replace', action='write')
      write(unit, '(A)') &
         '# elem_id    node_id        X            Y            Z         &
         &Damage        N_failure'
      do i = 1, n
         write(unit, '(I9,I11,5ES15.6)') sf%elem_id(i), sf%node_id(i), sf%coords(:,i), D(i), N_fail(i)
      end do
      close(unit)

      print '(A,A)', "  écrit : ", trim(filename)
   end subroutine write_damage_field

   ! write the damage field as CSV in the format expected by the external data block of Workbench 
   subroutine write_damage_csv(filename, sf, D)

      character(len=*), intent(in)     :: filename
      type(stress_field_t), intent(in) :: sf 
      real(dp), intent(in)             :: D(:)
      integer                          :: unit, i

      open(newunit=unit, file=filename, status='replace', action='write')
      write(unit, '(A)') "X,Y,Z,Damage"
      do i = 1, size(D)
         write(unit, '(F12.4,",",F12.4,",",F12.4,",",ES14.6)') &
               sf%coords(1,i), sf%coords(2,i), sf%coords(3,i), D(i)
      end do
      close(unit)
      print '(A,A)', "  écrit : ", trim(filename)

   end subroutine write_damage_csv

   ! read the numpy binary cache produced by the Python driver 
   ! n points records, each of 11 float64 values
   ! the number of points if inferred from the file size (11*8=88 bytes per point, that's why there is a verification on fsize -> multiple of 88)
   subroutine read_stress_field_binary(filename,sf,stress_scale)
      character(len=*), intent(in)      :: filename 
      type(stress_field_t), intent(out) :: sf 
      real(dp), intent(in)              :: stress_scale

      integer         :: unit, ios, i, n 
      integer(kind=8) :: fsize 
      real(dp), dimension(:), allocatable :: buf(:)

      !Infer point count from file size (11 float64 = 88 bytes per point)
      inquire(file=filename, size=fsize)
      if(mod(fsize,11_8*8_8)/=0_8) then 
         print '(A,A)', "ERROR : binary file size not multiple of 88"
         stop 1 
      end if 
      n = int(fsize/(11_8*8_8))

      sf%n_points = n 
      allocate(sf%elem_id(n), sf%node_id(n), sf%coords(3,n), sf%sigma(3,3,n), buf(11*n))

      !read the whole file in one shot as a raw byte stream 
      open(newunit=unit, file=filename, access = 'stream', form = 'unformatted', status = 'old', iostat = ios)
      if(ios /= 0) then 
         print '(A,A)', "ERROR : opening binary file : ", trim(filename)
         stop 1 
      end if 
      read(unit) buf 
      close(unit)

      ! unpack the flat buffer into the stress_field_t structure
      do i =1,n 
         sf%elem_id(i) = nint(buf((i-1)*11+1))
         sf%node_id(i) = nint(buf((i-1)*11+2))
         sf%coords(1,i) = buf((i-1)*11+3)
         sf%coords(2,i) = buf((i-1)*11+4)
         sf%coords(3,i) = buf((i-1)*11+5)
         ! the .bin cache stores the raw exported values, so the unit conversion
         ! happens here exactly as it does for the ASCII reader
         sf%sigma(1,1,i) = buf((i-1)*11+6)*stress_scale
         sf%sigma(2,2,i) = buf((i-1)*11+7)*stress_scale
         sf%sigma(3,3,i) = buf((i-1)*11+8)*stress_scale
         sf%sigma(1,2,i) = buf((i-1)*11+9)*stress_scale
         sf%sigma(2,1,i) = sf%sigma(1,2,i)
         sf%sigma(2,3,i) = buf((i-1)*11+10)*stress_scale
         sf%sigma(3,2,i) = sf%sigma(2,3,i)
         sf%sigma(1,3,i) = buf((i-1)*11+11)*stress_scale
         sf%sigma(3,1,i) = sf%sigma(1,3,i)
      end do 

      deallocate(buf)
      print '(A,I8,A,A)', " ", n, " points (binary) in  ", trim(filename)

   end subroutine read_stress_field_binary

   ! auto-selected reader , chosose the .txt file (loadcase.txt and also the binary file to read the stress values, if it was created before)
   subroutine read_stress_field_auto(txtname, sf, stress_scale)
      character(len=*), intent(in)      :: txtname
      type(stress_field_t), intent(out) :: sf
      real(dp), intent(in)              :: stress_scale
      
      character(len=512) :: binname 
      integer            :: n 
      logical            :: bin_exists 

      ! build the .bin path from the .txt path 
      n = len_trim(txtname)
      if (n>=4 .and. txtname(n-3:n) == '.txt') then 
         binname = txtname(1:n-4) // '.bin'
      else 
         binname = trim(txtname) // '.bin'
      end if 

      inquire(file=trim(binname), exist=bin_exists)
      if (bin_exists) then 
         call read_stress_field_binary(trim(binname), sf, stress_scale)
      else 
         call read_stress_field(trim(txtname), sf, stress_scale)
      end if 
   end subroutine read_stress_field_auto

end module osr_io_mod

! shared state of the current run
! this module holds all the variables that need to be visible from the callback functions sigma_at_t / dsigma_at_t without appearing in their signatures
module pipeline_state_mod
   use osr_model_mod, only: dp
   implicit none

   ! Number of unit loading cases 
   integer, save :: n_cases = 0

   ! loading mode selector
   integer, parameter :: MODE_SINUS = 0, MODE_HISTORY = 1
   integer, save      :: loading_mode = MODE_SINUS 

   ! integrator selector 
   character(len=32), public :: integrator_mode = "AUTO"

   ! unit of the stresses exported by Ansys, and the factor that converts them
   ! into the MPa the material parameters are calibrated for.
   ! An Ansys model solved in MKS exports Pa, so it needs 1.0e-6 ; a model solved
   ! in the mm-kg-N unit system exports MPa directly, so it needs 1.0.
   character(len=32), save :: stress_unit  = "MPA"
   real(dp),          save :: stress_scale = 1.0_dp

   ! history mode storage : sample times and tabulated lambda values 
   integer, save      ::  n_history = 0 
   real(dp), dimension(:), allocatable, save   :: history_time   ! sample time 
   real(dp), dimension(:,:), allocatable, save :: history_lambda ! lambda_k at each sample

   ! stress tensor sigma^(k)  of the current material point for each loading case at current point : shape (3, 3, n_cases)
   real(dp), allocatable, save :: sigma_at_point(:,:,:)
   !$OMP THREADPRIVATE(sigma_at_point)

   ! temporal loading parameters (sinusoidal) :
   ! lambda_k(t) = mean_k + amp_k * sin(omega*t - phase_k)
   real(dp),              save :: omega = 0.0_dp
   real(dp), allocatable, save :: amp(:), phase(:), mean(:)
   integer, save               :: n_cycles_max, n_steps_per_cycle
end module pipeline_state_mod

! module -> sigma(t) reconstruction 
! provide the two callbacks passed to osr_integrate_rk4 : sigma_at_t(t) and dsigma_at_t(t)
! both are built by superposition of the unit stress fields with time coefficients lambda_k(t)
module pipeline_callbacks_mod
   use osr_model_mod, only: dp
   use pipeline_state_mod
   implicit none

   integer, save :: cached_history_interval = 1
   !$OMP THREADPRIVATE(cached_history_interval)

contains

   ! sigma(t) = sum_k lambda_k(t) * sigma^(k)
   function sigma_at_t(t) result(s)
      real(dp), intent(in) :: t
      real(dp) :: s(3,3)
      real(dp) :: lambda_k, tau, denom
      integer  :: k, j

      s = 0.0_dp
      
      if (loading_mode == MODE_HISTORY) then 
         j = find_history_interval(t)
         denom = history_time(j+1)-history_time(j)

         if (denom<=0.0_dp) then 
            tau = 0.0_dp 
         else 
            tau = (t-history_time(j))/denom
         end if 

         do k = 1, n_cases
            lambda_k = (1.0_dp-tau)*history_lambda(j,k) + tau*history_lambda(j+1,k)
            s = s + lambda_k*sigma_at_point(:,:,k)
         end do 
      
      else 
         
         do k = 1, n_cases
            lambda_k = mean(k) + amp(k) * sin(omega * t + phase(k))
            s = s + lambda_k * sigma_at_point(:,:,k)
         end do 

      end if 

   end function sigma_at_t

   ! dsigma/dt = sum_k dlambda_k/dt * sigma^(k)
   ! in HISTORY mode : piecewise constant on each interval, given by the slope of the current segment 
   function dsigma_at_t(t) result(s)
      real(dp), intent(in) :: t
      real(dp) :: s(3,3)
      real(dp) :: dlambda_k, dt_hist
      integer  :: k, j 

      s = 0.0_dp
      if (loading_mode == MODE_HISTORY) then 
         j = find_history_interval(t)
         dt_hist = history_time(j+1)-history_time(j)

         do k = 1, n_cases
 
            dlambda_k = (history_lambda(j+1,k)-history_lambda(j,k))/dt_hist
            s = s +dlambda_k*sigma_at_point(:,:,k)
         end do 

      else 
         do k = 1, n_cases 
            dlambda_k = amp(k)*omega*cos(omega*t+phase(k))
            s = s+dlambda_k*sigma_at_point(:,:,k)
         end do 
      end if 

   end function dsigma_at_t

   ! for history mode  
   ! find the interval containing time t 
   ! values outside the tabulated range are clamped to the first/last interval
   integer function find_history_interval(t) result(idx)
      real(dp), intent(in) :: t
      integer              :: j 

      ! bounds
      if (t<= history_time(1)) then 
         idx=1
         cached_history_interval = idx
         return 
      end if 

      if (t>= history_time(n_history)) then 
         idx = n_history-1 
         cached_history_interval = idx
         return 
      end if 

      !O(1) amortised path : try cached interval first 
      j = min(max(cached_history_interval,1), n_history-1)

      if (t>=history_time(j) .and. t<history_time(j+1)) then 
         idx = j 
         return 
      end if 

      ! try next interval 
      if (j<n_history-1) then 
         if (t>=history_time(j+1) .and. t<history_time(j+2)) then 
            idx = j+1 
            cached_history_interval = idx 
            return 
         end if 
      end if 

      !  try previous interval, useful when a new material points starts from t=0 
      if (j>1) then 
         if (t>=history_time(j-1) .and. t<history_time(j)) then 
            idx = j-1 
            cached_history_interval = idx 
            return 
         end if 
      end if 

      ! fallback: binary search, O(log N)
      idx = binary_search_history_interval(t)
      cached_history_interval = idx 

   end function find_history_interval

   integer function binary_search_history_interval(t) result(idx)
      real(dp), intent(in) :: t 
      integer              :: lo, hi, mid 

      lo = 1 
      hi = n_history-1

      do while (lo<=hi)
         mid = (lo+hi)/2

         if (t<history_time(mid)) then 
            hi = mid-1 
         else if (t>=history_time(mid+1)) then 
            lo = mid+1 
         else 
            idx = mid 
            return 
         end if 
      end do 

      ! safety fallback 
      idx = min(max(lo,1), n_history-1)

   end function binary_search_history_interval

   ! linear interpolation of lambda_k between two tabulated samples 
   real(dp) function lambda_history(k,t) result(lam)
      integer, intent(in)  :: k 
      real(dp), intent(in) :: t
      
      integer              :: j 
      real(dp)             :: tau, dt 

      j = find_history_interval(t)
      dt = history_time(j+1) - history_time(j)

      if (dt <= 0.0_dp) then 
         lam = history_lambda(j,k)
         return 
      end if 

      tau = (t-history_time(j))/dt 
      lam = (1.0_dp - tau) * history_lambda(j,k) + tau*history_lambda(j+1,k)
   end function lambda_history 

   ! time derivative of the interpolated lambda_k : slope of the current segment (constant on the interval)
   real(dp) function dlambda_history(k,t) result(dlam)
      integer, intent(in)  :: k 
      real(dp), intent(in) :: t  

      integer :: j 
      real(dp) :: dt 

      j = find_history_interval(t)
      dt = history_time(j+1)-history_time(j)

      if (dt<= 0.0_dp) then 
         dlam = 0.0_dp 
      else 
         dlam = (history_lambda(j+1,k) - history_lambda(j,k))/dt 
      end if 
   end function dlambda_history 

end module pipeline_callbacks_mod