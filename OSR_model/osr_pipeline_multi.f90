! driver programof the OSR pipeline
! reads the loading config and the per-case stress fields exported by Ansys, then integrates the OSR ODEs at every mesh point in parallel
! and writes the resulting damage field back to disk visualisation 

program osr_pipeline_multi

    use osr_io_mod
    use osr_model_mod
    use pipeline_callbacks_mod, only: sigma_at_t, dsigma_at_t
    use pipeline_state_mod
    use omp_lib  ! for parallelisation 

    implicit none 

    type(osr_params) :: p 
    
    ! params of the material 
    p%A        = 0.225_dp
    p%sigma_oe = 490.0_dp
    p%C        = 1.25_dp
    p%K        = 2.65e-5_dp
    p%L        = 14.4_dp

    ! run the whole pipeline end-to-end for these material parameters 
    call run_pipeline(p)

contains

    ! reads the configuration and stress fields, integrates the OSR ODEs at every mesh point in parallel, and writes the damage field 
    subroutine run_pipeline(p)
        type(osr_params), intent(in) :: p 

        character(len=256), dimension(:),  allocatable :: fnames
        type(stress_field_t),dimension(:), allocatable :: sf
        real(dp), dimension(:), allocatable            :: D_field, N_fail
        real(dp), dimension(3,3)                       :: alpha_local
        real(dp)                                       :: D_local, t_failure 
        logical                                        :: failed
        integer                                        :: i, k, n_skipped, N_pts, n_failed_pts, report_step, n_cycles_max, n_steps_per_cycle, n_substeps_history
        character(len=32)                              :: selected_integrator

        ! timing 
        integer(8) :: t_start_int, t_end_int, count_rate 
        real(dp)   :: t_wall_seconds

        ! 1. read the loading configuration and also fills the module variables of pipeline_state_mod 
        call read_load_config("loadcases.txt",fnames, n_cycles_max, n_steps_per_cycle, n_substeps_history, selected_integrator)

        if (trim(selected_integrator) == "AUTO") then 

            select case (loading_mode)

            case (MODE_SINUS)
                selected_integrator = "RK4"

            case (MODE_HISTORY)
                selected_integrator = "DOPRI_EVENTS"

            case default 
                print *, "ERROR: unknown loading mode"
                stop 1

            end select 
        end if 

        print *, "Selected integrator = ", trim(selected_integrator)

        if (loading_mode==MODE_HISTORY) then 
            if (n_history<2) then 
                print *, "ERROR : HISTORY mode but n_history<2 qfter reading history_file"
                stop 1
            end if 
        end if 

        print '(A,I0,A,ES12.5)',"", n_cases, "loading case, omega =",omega
        print '(A,A,A,ES12.5,A)', " stress unit of the export = ", trim(stress_unit), &
              " (x", stress_scale, " -> MPa)"
        print *,"n_cycles_max = ", n_cycles_max,"", "n_steps_per_cycle = ", n_steps_per_cycle
        print '(A,I0)', "OpenMP threads used : ", omp_get_max_threads()
        
        ! 2. read the stress field of every unit loading case  
        allocate(sf(n_cases),sigma_at_point(3,3,n_cases))
        print '(A)', "reading stresses fields"
        
        do k=1,n_cases 
            call read_stress_field_auto(trim(fnames(k)),sf(k),stress_scale)
        end do 

        ! 3. consistency check : all sf(k) must have the same number of points 
        N_pts = sf(1)%n_points
        do k = 2, n_cases
            if (sf(k)%n_points/=N_pts) then
                print*, "ERROR : not the same number of points between cases"
                stop 1 
            end  if 
        end do 

        ! 3b. the material parameters are calibrated in MPa : refuse to integrate a
        ! stress field that cannot be in that unit rather than returning D = 1 everywhere
        call check_stress_magnitude(p, sf, N_pts)

        ! 4. main loop on integration points, parallelised with OpenMP
        ! PRIVATE : alpha_local, D_local, t_failure, failed, k -> each thread carries its own copy of these local variables  
        ! SCHEDULE(dynamic,32) : dynamic scheduling matters here because the cost per point is highly non-uniform bc a safe point costs cheap call to point is safe 
        ! while a point that fails runs a full RK4 integration until D reaches 1 
        allocate(D_field(N_pts), N_fail(N_pts))
        D_field = 0.0_dp
        N_fail = huge(1.0_dp)
        report_step = max(1, N_pts/10) ! to show the progress 10% by 10%
        n_failed_pts = 0 ! number of failed points
        n_skipped = 0

        call system_clock(t_start_int, count_rate)

        !$OMP PARALLEL DO                          &
        !$OMP&    PRIVATE(alpha_local, D_local,    &
        !$OMP&            t_failure, failed, k)    &
        !$OMP&    REDUCTION(+: n_skipped, n_failed_pts) &
        !$OMP&    SCHEDULE(dynamic, 32)
        
        do i = 1, N_pts
            
            ! allocation per thread (one time per thread)
            if (.not. allocated(sigma_at_point)) then
                allocate(sigma_at_point(3,3,n_cases))
            end if

            ! copy the tensors of the current point (3x3 per case) into the thread-private sigma_at_point, which the callbacks 
            ! sigma_at_t and dsigma_at_t will then read at each RK4 sub-step
            do k = 1, n_cases
                sigma_at_point(:,:,k) = sf(k)%sigma(:,:,i)
            end do

            ! analytical pre-filter
            if (point_is_safe(p)) then
                D_field(i)  = 0.0_dp
                n_skipped   = n_skipped + 1
                cycle
            end if

            alpha_local = 0.0_dp
            D_local     = 0.0_dp

            ! integration : one call to the generic RK4 integrator per non-skipped point, with the time budget adapted to the loading mode  
            select case (loading_mode)

            case (MODE_SINUS)

                select case (trim(selected_integrator))

                case ("RK4")
                    
                    call osr_integrate_rk4(p, sigma_at_t, dsigma_at_t, &
                                    0.0_dp, real(n_cycles_max,dp), n_cycles_max*n_steps_per_cycle, &
                                    alpha_local, D_local, t_failure, failed)

                case default 

                    print *, "ERROR : integrator ", trim(selected_integrator)," is not supported for SINUS mode"
                    stop 1 
                end select 

            case (MODE_HISTORY)

                select case (trim(selected_integrator))

                case("RK4")
                ! fixed step RK4 integration interval by interval over the tabulated history 
                    call integrate_history_rk4_by_interval(p, n_substeps_history, &
                                        alpha_local, D_local, t_failure, failed)

                case("DOPRI")
                    call integrate_history_dopri_by_interval(p, n_substeps_history, .false., &
                                                             alpha_local, D_local, t_failure, failed)
                                                             
                case("DOPRI_EVENTS")
                    call integrate_history_dopri_by_interval(p, n_substeps_history, .true., &
                                                             alpha_local, D_local, t_failure, failed)

                case default 
                    print *, "ERROR: unknown integrator = ", trim(selected_integrator)
                    stop 1 
                end select 

            case default 

                print *, "ERROR: unknown loading mode"
                stop 1 

            end select 
                                
            D_field(i) = D_local
            if (failed) then
                N_fail(i)     = t_failure
                n_failed_pts  = n_failed_pts + 1
            end if

            ! progress : approximative in parallel (the order in which threads reach this line is arbitrary, so this is only an approximate indicator)
            if (mod(i, report_step) == 0) then

                !$OMP CRITICAL 
                print '(A,I8,A,I8)', " progress ~", i, " / ", N_pts
                !$OMP END CRITICAL
            end if
        end do
        
        !$OMP END PARALLEL DO 

        call system_clock(t_end_int)
        t_wall_seconds = real(t_end_int-t_start_int,dp)/real(count_rate,dp)
        print '(A,F10.3,A)', "Integration wall time : ", t_wall_seconds," s"
        if (N_pts-n_skipped>0) then 
            print '(A,ES12.4,A)', "Time per integrated point : ", t_wall_seconds/real(N_pts-n_skipped,dp), "s"
        end if 

        ! write the damage field 
        call write_damage_field("damage_field.txt",sf(1), D_field, N_fail)
        call write_damage_csv("damage.csv", sf(1), D_field)

        print '(A,I8)', " total points : ",N_pts
        print '(A,I8)', " skipped points : ", n_skipped
        print '(A,I8)', " failed points : ", n_failed_pts
        if (n_failed_pts>0) then
            print '(A,ES12.4)', "N failure min : ", minval(N_fail)
            print '(A,ES12.4)', "N failure max(finite) : ", maxval(N_fail, mask=(N_fail<huge(1.0_dp)))
        end if 
    end subroutine run_pipeline

    ! max over the loading of beta(t, alpha=0) for the material point currently
    ! loaded in sigma_at_point
    ! for SINUS mode : we only need the location of the maximum : the exact value is recomputed properly by the RK4 integrator anyway
    ! for HISTORY mode : use the tabulated instants themselves. The reconstructed signal is piecewise linear, so beta reaches its extrema at the samples 
    ! so sampling them exactly is also exact (not just economical)
    real(dp) function max_beta_at_point(p) result(maxbeta)
        type(osr_params), intent(in) :: p 

        integer, parameter :: nsamp = 64 ! time division 
        integer            :: m 
        real(dp)           :: s(3,3), seff, t, alpha0(3,3), dvm(3,3), period, beta

        alpha0 = 0.0_dp 
        maxbeta = -huge(1.0_dp)

        if (loading_mode == MODE_HISTORY) then 

            ! evaluate beta at each tabulated instant 
            do m = 1, n_history
                t = history_time(m)
                s= sigma_at_t(t)
                call osr_beta(p, s, alpha0, beta, seff, dvm)

                if (beta>maxbeta) then 
                    maxbeta = beta 
                end if 
            end do 

        else 

            ! evaluate beta on 64 equally-spaced samples of one period 
            period = 2.0_dp*acos(-1.0_dp)/omega
            do m = 0, nsamp-1
                t = period*real(m,dp)/real(nsamp,dp) ! equal step 
                s = sigma_at_t(t)
                call osr_beta(p,s,alpha0,beta,seff,dvm)
                if ( beta>maxbeta ) then
                    maxbeta = beta 
                end if 
            end do  
        
        end if    

    end function max_beta_at_point

    ! pre-filter : a point whose beta never reaches zero never damages, so it needs no integration 
    logical function point_is_safe(p) result(safe)
        type(osr_params), intent(in) :: p 

        safe = (max_beta_at_point(p) < 0.0_dp)

    end function point_is_safe

    ! resolve the stress unit declared in loadcases.txt into the factor that brings
    ! the exported stresses to the MPa the material parameters are calibrated for.
    ! A bare number is accepted as an explicit factor, for a unit not listed here.
    subroutine set_stress_unit(token)
        character(len=*), intent(in) :: token

        character(len=32) :: t
        real(dp)          :: factor
        integer           :: ios, j

        t = adjustl(token)
        do j = 1, len_trim(t)
            if (t(j:j) >= 'a' .and. t(j:j) <= 'z') t(j:j) = achar(iachar(t(j:j)) - 32)
        end do

        stress_unit = t

        select case (trim(t))
        case ("MPA", "N/MM2", "NMM")
            stress_scale = 1.0_dp
        case ("PA", "N/M2", "MKS", "SI")
            stress_scale = 1.0e-6_dp
        case ("KPA")
            stress_scale = 1.0e-3_dp
        case ("GPA")
            stress_scale = 1.0e3_dp
        case ("PSI")
            stress_scale = 6.894757e-3_dp
        case ("KSI")
            stress_scale = 6.894757_dp
        case default
            read(t, *, iostat=ios) factor
            if (ios /= 0 .or. factor <= 0.0_dp) then
                print *, "ERROR : unknown stress unit in loadcases.txt : ", trim(t)
                print *, "Expected PA, KPA, MPA, GPA, PSI, KSI, or a positive conversion factor to MPa."
                stop 1
            end if
            stress_scale = factor
            stress_unit  = "custom"
        end select

    end subroutine set_stress_unit

    ! guard against a stress field that is not in the declared unit
    ! the OSR damage rate is K*exp(L*beta), so a field exported in Pa but read as MPa
    ! gives beta of the order of 1.0e5 at every point : the exponential saturates, D
    ! reaches 1 during the very first step and the whole model "fails at zero cycle".
    ! That answer looks like a result but is only a unit mismatch, so stop instead.
    subroutine check_stress_magnitude(p, sf, N_pts)
        type(osr_params), intent(in)     :: p
        type(stress_field_t), intent(in) :: sf(:)
        integer, intent(in)              :: N_pts

        ! beta = 10 means an equivalent stress of 11*sigma_oe, ie ~5.4 GPa for steel :
        ! above any real fracture stress, so such a field cannot be in MPa
        real(dp), parameter :: BETA_ABSURD = 10.0_dp

        real(dp) :: beta_max, b, sigma_peak
        integer  :: i, k, n_active

        beta_max = -huge(1.0_dp)
        n_active = 0

        !$OMP PARALLEL DO PRIVATE(k, b) REDUCTION(max: beta_max) REDUCTION(+: n_active) SCHEDULE(static)
        do i = 1, N_pts
            if (.not. allocated(sigma_at_point)) then
                allocate(sigma_at_point(3,3,n_cases))
            end if
            do k = 1, n_cases
                sigma_at_point(:,:,k) = sf(k)%sigma(:,:,i)
            end do
            b = max_beta_at_point(p)
            beta_max = max(beta_max, b)
            if (b >= 0.0_dp) n_active = n_active + 1
        end do
        !$OMP END PARALLEL DO

        ! beta = (sigma_eff + A*I1)/sigma_oe - 1, so this is the peak of sigma_eff + A*I1
        sigma_peak = (1.0_dp + beta_max)*p%sigma_oe

        print '(A,ES12.4,A,I0,A,I0,A)', " max beta(alpha=0) over the mesh = ", beta_max, &
              "  (", n_active, " / ", N_pts, " points above the endurance surface)"

        if (beta_max > BETA_ABSURD) then
            print *, ""
            print *, "ERROR : the stress field is not consistent with the material parameters."
            print '(A,ES12.4,A)', " Peak of sigma_eff + A*I1 read by the solver : ", sigma_peak, " MPa"
            print '(A,F10.2,A)',  " Endurance limit sigma_oe                    : ", p%sigma_oe, " MPa"
            print '(A,A,A,ES12.4,A)', " Declared stress unit of the export          : ", trim(stress_unit), &
                  " (factor ", stress_scale, ")"
            print *, ""
            print *, "Every point would fail during the first step, which is a unit mismatch,"
            print *, "not a fatigue result. An Ansys model solved in MKS exports stresses in Pa:"
            print *, "declare the unit on the first data line of loadcases.txt, for instance"
            print *, ""
            print *, "    1 6.283185307 SINUS PA"
            print *, ""
            print *, "(PA, KPA, MPA, GPA, PSI, KSI, or an explicit factor to MPa are accepted.)"
            stop 1
        end if

    end subroutine check_stress_magnitude

    ! parser for loadcases.txt 
    subroutine read_load_config(filename, fnames, n_cycles_max, n_steps_per_cycle, n_substeps_history, selected_integrator)

        character(len=*), intent(in)                               :: filename 
        character(len=256), dimension(:), allocatable, intent(out) :: fnames 
        integer, intent(out)                                       :: n_cycles_max, n_steps_per_cycle, n_substeps_history
        character(len=32), intent(out)                             :: selected_integrator

        integer  :: ios, unit, k, nc, nc_max, nsteps, nsubhist 
        real(dp) :: om, a, ph, mn ! omega, amplitude, phase, mean 
        character(len=1024) :: line 
        character(len=256)  :: fn, history_file
        character(len=32)   :: mode_str, unit_str
        
        open (newunit = unit, file = filename, status = 'old', action = 'read', iostat = ios)
        if (ios /= 0 ) then 
            print '(A,A)', "error with the opening of the file :", trim(filename)
            stop 1 
        end if 

        !first useful line 
        call next_data_line(unit, line)
        read(line,*) nc, om, mode_str ! number of unit loading case and omega 
        ! the stress unit is an optional 4th token : a file written before it existed
        ! declares nothing, which means MPa (the unit the model has always assumed)
        call nth_token(line, 4, unit_str)
        if (len_trim(unit_str) == 0) unit_str = "MPA"
        omega = om 
        n_cases = nc 
        call set_stress_unit(unit_str)
        select case (trim(mode_str))

        case ("SINUS")
            loading_mode = MODE_SINUS
        case ("HISTORY")
            loading_mode = MODE_HISTORY
        case default
            print *, "ERROR : unknown loading mode in loadcases.txt"
            print *, "Expected SINUS or HISTORY"
            stop 1 

        end select 

        ! reallocate the SINUS parameter arrays for this run
        if(allocated(amp)) deallocate(amp)
        if(allocated(phase)) deallocate(phase)
        if(allocated(mean)) deallocate(mean)

        allocate (fnames(nc), amp(nc), phase(nc), mean(nc))

        ! n_cases lines 
        do k = 1, n_cases
            ! filling the tables 
            call next_data_line(unit,line)
            read(line,*) fn, a, ph, mn 
            fnames(k) = fn 
            amp(k) = a 
            phase(k) = ph 
            mean(k) = mn 
        end do 
        
        call next_data_line(unit,line)
        read(line,*) nc_max, nsteps, nsubhist
        n_cycles_max = nc_max
        n_steps_per_cycle = nsteps
        n_substeps_history = nsubhist

        ! history filename
        call next_data_line(unit,line)
        read(line,*) history_file

        ! integrator 
        call next_data_line(unit,line)
        read(line,*) selected_integrator
        selected_integrator = adjustl(selected_integrator)

        close(unit)

        if (loading_mode == MODE_HISTORY) then 
            if (trim(history_file)=="-") then 
                print *, "ERROR : HISTORY mode selected but history file is '-'"
                stop 1 
            end if 

            call read_history_file(trim(history_file))

            print *, "Loading mode = HISTORY"
            print *, "n_history = ", n_history 
            print *, "history first time = ", history_time(1)
            print *, "history last time = ", history_time(n_history)
        else 
            print *, "loading mode = 'SINUS'"
        end if 

    end subroutine read_load_config

    ! read the next non-blank, non-comment line from an open file 
    ! the caller passes a character variable of any length : declaring the dummy
    ! with an assumed length keeps the whole of it defined. With a fixed len=256
    ! dummy, a longer actual argument kept undefined bytes past column 256.
    subroutine next_data_line(unit, line)

        integer, intent(in)             :: unit 
        character(len=*), intent(out)   :: line 

        integer :: ios

        do 
            read (unit, '(A)', iostat = ios) line  
            if ( ios/=0 ) then 
                print *, " error : end of the file "
                stop 1 
            end if 
            line = adjustl(line)
            if (len_trim(line)==0) cycle
            if (line(1:1)=="#") cycle 
            return 
        end do 

    end subroutine next_data_line

    ! return the n-th blank separated token of a line, blank if the line is shorter
    subroutine nth_token(line, n, token)

        character(len=*), intent(in)  :: line
        integer,          intent(in)  :: n
        character(len=*), intent(out) :: token

        integer :: i, first, count_tok, len_line

        token = ""
        len_line = len_trim(line)
        i = 1
        count_tok = 0

        do while (i <= len_line)
            ! .and. is not short-circuiting in Fortran, so the bound is tested first
            do
                if (i > len_line) exit
                if (line(i:i) /= " ") exit
                i = i + 1
            end do
            if (i > len_line) exit

            first = i
            do
                if (i > len_line) exit
                if (line(i:i) == " ") exit
                i = i + 1
            end do

            count_tok = count_tok + 1
            if (count_tok == n) then
                token = line(first:i-1)
                return
            end if
        end do

    end subroutine nth_token

    ! parser for history.txt
    subroutine read_history_file(filename)

        character(len=*), intent(in) :: filename
        integer :: unit, ios, n, j, k 
        character(len=1024) :: line 

        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios/=0) then 
            print '(A,A)', "ERROR can not open history file : ", trim(filename)
            stop 1
        end if 

        ! count usable lines
        n=0 
        do
            read(unit,'(A)', iostat=ios) line 
            if (ios/=0) exit

            line = adjustl(line)
            if (len(trim(line))==0) cycle
            if (line(1:1)=='#') cycle 
            n=n+1
        end do 

        if (n<2) then 
            print *, "ERROR : history file must contain at least two data lines"
            stop 1
        end if 

        n_history = n 

        if(allocated(history_time)) deallocate(history_time)
        if(allocated(history_lambda)) deallocate(history_lambda)

        allocate(history_time(n_history))
        allocate(history_lambda(n_history,n_cases))

        ! rewind and actually read the samples
        rewind(unit)

        j=0
        do 
            read(unit,'(A)',iostat=ios) line 
            if (ios/=0) exit

            line = adjustl(line)
            if (len(trim(line))==0) cycle
            if (line(1:1)=='#') cycle 
            j=j+1

            read(line, *, iostat=ios) history_time(j), (history_lambda(j,k), k=1,n_cases)
            if (ios/=0) then
                print *, "ERROR reading history line :"
                print *, trim(line)
                stop 1
            end if 
        end do 

        close(unit)

        print *, "history file read : ", trim(filename)
        print *, "n_history = ", n_history
        print *, "t_start = ", history_time(1), "t_end = ", history_time(n_history)

    end subroutine read_history_file 

    subroutine integrate_history_dopri_by_interval(p, n_substeps_history, use_events, alpha, D, t_failure, failed)

        type(osr_params), intent(in) :: p 
        integer, intent(in)          :: n_substeps_history
        logical, intent (in)         :: use_events
        real(dp), intent(inout)      :: alpha(3,3), D
        real(dp), intent(out)        :: t_failure 
        logical, intent(out)         :: failed 

        integer  :: ih, n_used 
        real(dp) :: t_failure_loc 
        logical  :: failed_loc
        
        failed = .false. 
        t_failure = -1.0_dp 

        do ih = 1, n_history-1
            n_used = 0 
            call osr_integrate_dopri(p, sigma_at_t, dsigma_at_t, history_time(ih), history_time(ih+1), &
                                     max(1,n_substeps_history), alpha, D, t_failure_loc, failed_loc, &
                                     use_events = use_events, n_steps_used = n_used)

            if (failed_loc) then 
                failed = .true. 
                t_failure = t_failure_loc 
                return 
            end if 
        end do 

    end subroutine integrate_history_dopri_by_interval

    subroutine integrate_history_rk4_by_interval(p, n_substeps_history, alpha, D, t_failure, failed)

    type(osr_params), intent(in) :: p
    integer, intent(in)          :: n_substeps_history
    real(dp), intent(inout)      :: alpha(3,3), D
    real(dp), intent(out)        :: t_failure
    logical, intent(out)         :: failed

    integer  :: ih
    real(dp) :: t_failure_loc
    logical  :: failed_loc
        
    failed = .false.
    t_failure = -1.0_dp

    do ih = 1, n_history-1

        call osr_integrate_rk4(p, sigma_at_t, dsigma_at_t, &
                               history_time(ih), history_time(ih+1), &
                               max(1,n_substeps_history), &
                               alpha, D, t_failure_loc, failed_loc)

        if (failed_loc) then
            failed = .true.
            t_failure = t_failure_loc
            return
        end if

    end do

end subroutine integrate_history_rk4_by_interval

end program osr_pipeline_multi