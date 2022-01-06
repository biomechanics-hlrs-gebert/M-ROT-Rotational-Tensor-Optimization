
!------------------------------------------------------------------------------
! MODULE: auxiliaries_of_tensor_optimizer
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
! @Description:
!> Module containing additional routines for the main program.
!------------------------------------------------------------------------------
MODULE auxiliaries_of_tensor_optimizer

USE global_std
USE mechanical
USE vtk_meta_data
USE raw_binary
USE opt_stiffness

IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: write_criteria_space_to_vtk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Write a criteria space to a vtk file. MAKE SURE (!) that all the (global)
!> variables are not deallocated or "reallocated" (!) It simplifies interfaces,
!> but it is not a clean-code paradigm...
!
!> @param[in] file File to write to.
!------------------------------------------------------------------------------  
SUBROUTINE write_criteria_space_to_vtk(file)

    CHARACTER(LEN=*), INTENT(IN) :: file

    INTEGER(KIND=ik), DIMENSION(3) :: ttl_steps
    INTEGER(KIND=ik) :: fh
    LOGICAL :: fex

    fex = .FALSE.

    INQUIRE(FILE=TRIM(file), EXIST=fex)
    IF(fex) WRITE(std_out, FMT_WRN) TRIM(file)//" already exists."

    fh = give_new_unit()

    ttl_steps = (pm_steps*2_ik)+1_ik

    CALL write_vtk_struct_points_header(fh, file, TRIM('rk8'), &
        [1._rk, 1._rk, 1._rk], [0._rk, 0._rk, 0._rk], ttl_steps)

    CALL ser_write_raw(fh, file, crit)

    CALL write_vtk_struct_points_footer(fh, file)

END SUBROUTINE write_criteria_space_to_vtk

END MODULE auxiliaries_of_tensor_optimizer

!------------------------------------------------------------------------------
! PROGRAM: Rotational Tensor Optimization
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> \date 27.10.2020
!> \modified 06.01.2022
!------------------------------------------------------------------------------
PROGRAM tensor_optimizer

USE global_std
USE meta
USE user_interaction
USE formatted_plain
USE MPI
USE opt_stiffness
USE auxiliaries_of_tensor_optimizer

IMPLICIT NONE

! Parameter
INTEGER(KIND=ik), PARAMETER :: debug = 2   ! Choose an even integer!!
REAL   (KIND=rk), PARAMETER :: lower_thres_percentage = 0.1

! Variables
TYPE(tensor_2nd_rank_R66) :: dummy
TYPE(tensor_2nd_rank_R66), DIMENSION(:), ALLOCATABLE :: tglbl_in
TYPE(tensor_2nd_rank_R66), DIMENSION(:,:), ALLOCATABLE :: tglbl_res
TYPE(tensor_2nd_rank_R66), DIMENSION(:), ALLOCATABLE :: tlcl_res

TYPE(materialcard) :: bone

CHARACTER(LEN=mcl), DIMENSION(:), ALLOCATABLE :: m_rry      
CHARACTER(LEN=scl) :: restart, restart_cmd_arg, re_mono, re_orth, re_ani1, re_ani2
CHARACTER(LEN=scl) :: suf_covo, suf_mono, suf_orth, suf_ani1, suf_ani2, temp_suf    
CHARACTER(LEN=scl) :: binary, dmn_no, suffix
CHARACTER(LEN=  8) :: date
CHARACTER(LEN= 10) :: time

REAL(KIND=rk) :: start, end

INTEGER(KIND=ik) :: fh_covo, fh_mono, fh_orth, fh_ani1, fh_ani2, fhwcrit
INTEGER(KIND=ik) :: exp_dmn_crit, covo_amnt_lines, zero_matrix_counter = 0
INTEGER(KIND=ik) :: jj, kk, mm, zz

LOGICAL, DIMENSION(4) :: execute_optimization = .FALSE.
LOGICAL :: stp, print_criteria

INTEGER(KIND=mik) :: ierr, my_rank, size_mpi, mii
INTEGER(KIND=mik) :: active, feed_ranks, crs, crs_counter
INTEGER(KIND=mik), DIMENSION(MPI_STATUS_SIZE) :: stmpi
INTEGER(KIND=mik), DIMENSION(:,:), ALLOCATABLE :: statuses_mpi
INTEGER(KIND=mik), DIMENSION(:), ALLOCATABLE :: req_list

INTEGER(KIND=mik) :: MPI_tensor_2nd_rank_R66
INTEGER(KIND=mik), DIMENSION(6) :: blocklen, dtype 
INTEGER(KIND=MPI_ADDRESS_KIND)  :: disp(6), base

! Initialize MPI Environment
CALL MPI_INIT(ierr)
CALL mpi_err(ierr,"MPI_INIT failed.")

CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
CALL mpi_err(ierr,"MPI_COMM_RANK couldn't be retrieved")

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)
CALL mpi_err(ierr,"MPI_COMM_SIZE couldn't be retrieved")

IF (size_mpi < 2) CALL print_err_stop(std_out, "At least two ranks required to execute this program.", 1)

!------------------------------------------------------------------------------
! Redirect std_out into a file in case std_out is not useful by environment.
! Place these lines before handle_lock_file :-)
!------------------------------------------------------------------------------
CALL MPI_GET_ADDRESS(dummy%dmn, disp(1), ierr) 
CALL MPI_GET_ADDRESS(dummy%density, disp(2), ierr) 
CALL MPI_GET_ADDRESS(dummy%doa_zener, disp(3), ierr) 
CALL MPI_GET_ADDRESS(dummy%doa_gebert, disp(4), ierr) 
CALL MPI_GET_ADDRESS(dummy%pos, disp(5), ierr) 
CALL MPI_GET_ADDRESS(dummy%mat, disp(6), ierr) 
	
base = disp(1) 
disp = disp - base 

blocklen(1:4) = 1 
blocklen(5) = 3 
blocklen(6) = 36 

dtype(1) = MPI_INTEGER8 
dtype(2:6) = MPI_DOUBLE_PRECISION

CALL MPI_TYPE_CREATE_STRUCT(6_mik, blocklen, disp, dtype, MPI_tensor_2nd_rank_R66, ierr) 
CALL mpi_err(ierr,"MPI_tensor_2nd_rank_R66 couldn't be created.")

CALL MPI_TYPE_COMMIT(MPI_tensor_2nd_rank_R66, ierr)
CALL mpi_err(ierr,"MPI_tensor_2nd_rank_R66 couldn't be commited.")

!------------------------------------------------------------------------------
! Initialize program itself
!------------------------------------------------------------------------------
IF(my_rank == 0) THEN

    !------------------------------------------------------------------------------
    ! Redirect std_out into a file in case std_out is not useful by environment.
    ! Place these lines before handle_lock_file :-)
    !------------------------------------------------------------------------------
    std_out = determine_stout()

    IF(std_out/=6) CALL meta_start_ascii(std_out, '.std_out')

    CALL show_title()
 
    IF(debug >=0) WRITE(*,FMT_MSG) "Post mortem info probably in ./datasets/.temporary.std_out"

    CALL CPU_TIME(start)

    !------------------------------------------------------------------------------
    ! Parse the command arguments
    !------------------------------------------------------------------------------
    CALL get_cmd_args(binary, in%full, stp, restart, restart_cmd_arg)
    IF(stp) GOTO 1001

    !------------------------------------------------------------------------------
    ! Check and open the input file; Modify the Meta-Filename / Basename
    ! Define the new application name first
    !------------------------------------------------------------------------------
    global_meta_prgrm_mstr_app = 'ROTO' 
    global_meta_program_keyword = 'ROT_TENSOR_OPT'
    CALL meta_append(m_rry)

    !------------------------------------------------------------------------------
    ! Parse input
    !------------------------------------------------------------------------------
    CALL meta_read(std_out, 'REQUEST_MONO'   , m_rry, re_mono)
    CALL meta_read(std_out, 'REQUEST_ORTH'   , m_rry, re_orth)
    CALL meta_read(std_out, 'REQUEST_ANI_1'  , m_rry, re_ani1)
    CALL meta_read(std_out, 'REQUEST_ANI_2'  , m_rry, re_ani2)
    CALL meta_read(std_out, 'RESTART'        , m_rry, restart)
    CALL meta_read(std_out, 'EXPORT_DMN_CRIT', m_rry, exp_dmn_crit)
    CALL meta_read(std_out, 'YOUNG_MODULUS'  , m_rry, bone%E)
    CALL meta_read(std_out, 'POISSON_RATIO'  , m_rry, bone%nu)

    !------------------------------------------------------------------------------
    ! Restart handling
    ! Done after meta_io to decide based on keywords
    !------------------------------------------------------------------------------
    CALL meta_handle_lock_file(restart, restart_cmd_arg)

    !------------------------------------------------------------------------------
    ! Spawn a log and a monitoring file (to check the behavior of MPI)
    !------------------------------------------------------------------------------
    CALL meta_start_ascii(fh_log, log_suf)
    
    IF(debug >= 0) CALL meta_start_ascii(fh_mon, mon_suf)

    CALL DATE_AND_TIME(date, time)

    WRITE(fh_log, FMT_TXT_SEP)  
    WRITE(fh_log, FMT_TXT)      TRIM(ADJUSTL(longname))//" Results"
    WRITE(fh_log, FMT_TXT)     "Date: "//date//" [ccyymmdd]"
    WRITE(fh_log, FMT_TXT)     "Time: "//time//" [hhmmss.sss]"
    WRITE(fh_log, FMT_TXT_SEP)  
    WRITE(fh_log, FMT_MSG_AI0) "Processors:", size_mpi  
    
    !------------------------------------------------------------------------------
    ! Create/Open tensor files. Basically tuned csv data.
    !------------------------------------------------------------------------------
    fh_covo = give_new_unit()
    suf_covo = ".tcr.covo" ! control volume (in situ orientation)
    CALL meta_existing_ascii(fh_covo, suf_covo, covo_amnt_lines)

    crs = 0_mik
    IF(TRIM(re_mono) == "YES") THEN
        fh_mono = give_new_unit()
        suf_mono = ".tcr.mono" ! monotropic optimization 
        CALL meta_start_ascii(fh_mono, suf_mono) 
        execute_optimization(1) = .TRUE. 
        crs=crs+1_mik
    END IF

    IF(TRIM(re_orth) == "YES") THEN
        fh_orth = give_new_unit()
        suf_orth = ".tcr.orth" ! orthotropic optimization 
        CALL meta_start_ascii(fh_orth, suf_orth) 
        execute_optimization(2) = .TRUE. 
        crs=crs+1_mik
    END IF

    IF(TRIM(re_ani1) == "YES") THEN
        fh_ani1 = give_new_unit()
        suf_ani1 = ".tcr.an1" ! anisotropic optimization, version 1 
        CALL meta_start_ascii(fh_ani1, suf_ani1) 
        execute_optimization(3) = .TRUE. 
        crs=crs+1_mik
    END IF

    IF(TRIM(re_ani2) == "YES") THEN
        fh_ani2 = give_new_unit()
        suf_ani2 = ".tcr.an2" ! anisotropic optimization, version 2 
        CALL meta_start_ascii(fh_ani2, suf_ani2) 
        execute_optimization(4) = .TRUE. 
        crs=crs+1_mik
    END IF

    !------------------------------------------------------------------------------
    ! Ensure, the user is sane and report if so (or not).
    !------------------------------------------------------------------------------
    IF(crs == 0_mik) THEN
        mssg = "No optimization requested. Program aborts."
        CALL print_err_stop(std_out, mssg, 1)
    END IF

    IF((exp_dmn_crit < -0_ik) .AND. (covo_amnt_lines*crs > 100_ik)) THEN
        mssg = "Are you sure to export the criteria spaces of more than 100 &
            &optimizations (dmns*criterias) to the file system via *.vtk files?! &
            &If yes, change the source and recompile :-)"
        CALL print_err_stop(std_out, mssg, 1)
    END IF

    IF(debug >= 0) THEN
        WRITE(std_out, FMT_TXT_SEP)
        WRITE(std_out, FMT_MSG_AI0) "Debug Level:", debug
        WRITE(std_out, FMT_MSG_AI0) "Processors:", size_mpi  
        WRITE(std_out, FMT_MSG_AI0) "Amount of domains:", covo_amnt_lines-1_ik
        WRITE(std_out, FMT_TXT_SEP)
        FLUSH(std_out)
    END IF

    !------------------------------------------------------------------------------
    ! Parse the global input data
    ! Allocate amnt_lines -1_ik due to header
    !------------------------------------------------------------------------------
    ALLOCATE(tglbl_in(covo_amnt_lines-1_ik))
    CALL parse_tensor_2nd_rank_R66(fh_covo, in%p_n_bsnm//TRIM(suf_covo), &
        covo_amnt_lines, tglbl_in)

    !------------------------------------------------------------------------------
    ! Allocate global results array. Need to loop over execute_optimization 
    ! to properly assign the tglbl_res(entries, :)
    !------------------------------------------------------------------------------
    ALLOCATE(tglbl_res(crs, covo_amnt_lines-1_ik))

END IF ! (my_rank == 0)

CALL MPI_BCAST( in%p_n_bsnm, INT(meta_mcl, KIND=mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(out%p_n_bsnm, INT(meta_mcl, KIND=mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST(execute_optimization, 4_mik, MPI_LOGICAL, 0_mik, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST(exp_dmn_crit, 1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(crs, 1_mik, MPI_INTEGER, 0_mik, MPI_COMM_WORLD, ierr)

!------------------------------------------------------------------------------
! All Ranks -- Init MPI request and status lists
!------------------------------------------------------------------------------
ALLOCATE(req_list(size_mpi-1_mik))
req_list=0_mik

ALLOCATE(statuses_mpi(MPI_STATUS_SIZE, size_mpi-1_mik))
statuses_mpi=0_mik

!------------------------------------------------------------------------------
! Rank 0 -- Process master - Start working process
!------------------------------------------------------------------------------
IF (my_rank==0) THEN

    mii = 1_mik
    feed_ranks = 1_mik 

    !------------------------------------------------------------------------------
    ! Supply all worker masters with their first tensor
    ! Master does not compute tensors
    !
    ! Basically, the master loops over all workers to send and receive from them
    !
    ! Each processor -> 1 domain. Therefore, nn (domain) = mii (process)
    ! @Struct process -> 1 nn can have multiple processes
    !------------------------------------------------------------------------------
    DO WHILE (mii <= covo_amnt_lines-1_ik)

        !------------------------------------------------------------------------------
        ! Check whether it is a zero tensor:
        !------------------------------------------------------------------------------
        IF(SUM(tglbl_in(mii)%mat) <= lower_thres_percentage * bone%E) THEN
            DO mm=1, crs
                tglbl_res(mm, mii)%density = 0._rk
                tglbl_res(mm, mii)%doa_zener = 0._rk
                tglbl_res(mm, mii)%doa_gebert = 0._rk
                tglbl_res(mm, mii)%pos = 0._rk
                tglbl_res(mm, mii)%mat = 0._rk
                tglbl_res(mm, mii)%dmn = tglbl_in(mii)%dmn
            END DO
            zero_matrix_counter = zero_matrix_counter + 1_ik

            CYCLE
        END IF

        !------------------------------------------------------------------------------
        ! Give feedback to user
        !------------------------------------------------------------------------------
        IF(std_out == 6_ik) THEN
            CALL EXECUTE_COMMAND_LINE('printf "\033c"')
            
            WRITE(std_out, FMT_TXT_A3I0) "Processed domains: ", mii, " of ", covo_amnt_lines-1_ik
            WRITE(std_out, FMT_TXT) "Most current input to compute:"
            WRITE(std_out, FMT_TXT_A3I0) "Total 0-matrices: ", zero_matrix_counter
            WRITE(std_out, FMT_TXT) ""
            WRITE(dmn_no, '(I0)') tin%dmn
            CALL write_matrix(std_out, "Domain "//TRIM(dmn_no), 'spl', '', tglbl_in(mii)%mat)
        END IF

        !------------------------------------------------------------------------------
        ! Amount of entities to compute. 
        ! -1_ik due to header of file.
        !------------------------------------------------------------------------------
        IF (feed_ranks > (size_mpi-1_mik)) feed_ranks = 1_ik 

        !------------------------------------------------------------------------------
        ! Call the domain an active one.
        !------------------------------------------------------------------------------
        CALL MPI_SEND(1_mik, 1_mik, MPI_INTEGER, feed_ranks, -1_mik, MPI_COMM_WORLD,ierr)
        CALL print_err_stop(std_out, "MPI_SEND of activity didn't succeed", INT(ierr, KIND=ik))

        !------------------------------------------------------------------------------
        ! Send tensor
        !------------------------------------------------------------------------------
        CALL MPI_SEND(tglbl_in(mii), 1_mik, MPI_tensor_2nd_rank_R66, feed_ranks, &
            feed_ranks, MPI_COMM_WORLD, ierr)
        CALL print_err_stop(std_out, "MPI_SEND of tin failed.", INT(ierr, KIND=ik))            
         
        !------------------------------------------------------------------------------
        ! Log to monitor file
        !------------------------------------------------------------------------------
        IF(debug >= 0) THEN
            WRITE(fh_mon, FMT_MSG_2AI0)"MPI rank: ", feed_ranks, " Domain number: ", tglbl_in(mii)%dmn
            FLUSH(fh_mon)
        END IF

        !------------------------------------------------------------------------------
        ! Place a receive for the tensor. 
        ! Tensors are not sortded automatically.
        !------------------------------------------------------------------------------
        CALL MPI_IRECV(tglbl_res(:, mii), crs, MPI_tensor_2nd_rank_R66, feed_ranks, &
            INT(tglbl_in(mii)%dmn, KIND=mik), MPI_COMM_WORLD, req_list(mii), ierr)
        CALL print_err_stop(std_out, "MPI_IRECV of activity(mii) failed.", INT(ierr, KIND=ik))

        feed_ranks = feed_ranks + 1_ik  
    END DO
    ! Still my_rank == 0

    CALL MPI_WAITALL(size_mpi-1_mik, req_list, statuses_mpi, ierr)
    CALL print_err_stop(std_out, "MPI_WAITANY on req_list for IRECV of tglbl_res failed.", &
        INT(ierr, KIND=ik))

    !------------------------------------------------------------------------------
    ! Stop workers properly by sending an activity = -1 info
    !------------------------------------------------------------------------------
    DO zz = 1_mik, size_mpi-1_mik
        CALL MPI_SEND(-1_mik, 1_mik, MPI_INTEGER, INT(zz, mik), -1_mik, MPI_COMM_WORLD,ierr)
        CALL print_err_stop(std_out, "MPI_SEND of activity didn't succeed", INT(ierr, KIND=ik))
    END DO
    
!------------------------------------------------------------------------------
! Ranks > 0 -- Workers
! Worker: Must use my_rank as tag
!------------------------------------------------------------------------------
ELSE
    !------------------------------------------------------------------------------
    ! Allocate local results
    !------------------------------------------------------------------------------
    ALLOCATE(tlcl_res(crs))

    DO
        !------------------------------------------------------------------------------
        ! Stop (gracefully) if workers receive an active = -1 information.
        !------------------------------------------------------------------------------
        CALL MPI_RECV(active, 1_mik, MPI_INTEGER, 0_mik, -1_mik, MPI_COMM_WORLD, stmpi, ierr)
        CALL print_err_stop(std_out, "MPI_RECV on active failed.", INT(ierr, KIND=ik))

        IF(active == -1) EXIT

        !------------------------------------------------------------------------------
        ! Receive tensor (derived type of tensor_2nd_rank_R66)
        !------------------------------------------------------------------------------
        CALL MPI_RECV(tin, 1_mik, MPI_tensor_2nd_rank_R66, 0_mik, my_rank, MPI_COMM_WORLD, stmpi, ierr)
        CALL print_err_stop(std_out, "MPI_RECV on tin failed.", INT(ierr, KIND=ik))

        !------------------------------------------------------------------------------
        ! Check whether to write the requested domain(s) to *.vtk. At this point,
        ! the program does not care how many of all domains are for export.
        !------------------------------------------------------------------------------
        print_criteria = .FALSE.

        IF((exp_dmn_crit == tin%dmn) .OR. (exp_dmn_crit < 0_ik)) print_criteria = .TRUE.

        WRITE(dmn_no, '(I0)') tin%dmn

        !------------------------------------------------------------------------------
        ! Start computation of domains by workers.
        !------------------------------------------------------------------------------
        crs_counter = 1_mik

        DO jj = 1_ik, 4_ik
            IF(.NOT. execute_optimization(jj)) CYCLE

            !------------------------------------------------------------------------------
            ! Optimization is capable of accepting dig /= 0._rk. Tensors can be optimized
            ! several times and do not need to be restarted from the position of the
            ! control volume within the bone.
            !
            ! Assign a domain specific variable to a global, currently validone
            !------------------------------------------------------------------------------
            dig = tin%pos 

            DO kk = 1_ik, 2_ik
                !------------------------------------------------------------------------------
                ! Variables are global ones, which are valid for each specific rank (!)
                !
                ! The optimization is done in two steps. First, search the best orientation
                ! by 1° steps. Second, search at this positon with 0.025° steps.
                !------------------------------------------------------------------------------
                IF(kk == 1_ik) THEN
                    intervall = 1._rk
                    pm_steps = 182_ik
                ELSE
                    tin%mat = tout%mat
                    dig = tout%pos
                    intervall = 0.025_rk
                    pm_steps = 80_ik
                END IF

                SELECT CASE(jj)
                    CASE(1); CALL opt_stiff('monotropic'); temp_suf = "mono"
                    CASE(2); CALL opt_stiff('orthotropic'); temp_suf = "orth"
                    CASE(3); CALL opt_stiff('anisotropic1'); temp_suf = "an1"
                    CASE(4); CALL opt_stiff('anisotropic2'); temp_suf = "an2"
                END SELECT        

                !------------------------------------------------------------------------------
                ! Tilts until S11 > S22 > S33 
                !------------------------------------------------------------------------------
                CALL tilt_tensor(tout%mat)
                
                tout%dmn = tin%dmn
                tout%doa_zener = doa_zener(tout%mat)
                tout%doa_gebert = doa_gebert(tout%mat)
                tout%density = gebert_density_voigt(tout%mat, bone%E, bone%nu)

            END DO

            !------------------------------------------------------------------------------
            ! At the end of the second step, the results acutally get written to the 
            ! output variables that are then send back to my_rank=0 / main process.
            !
            ! Domain always gets retrieved from tin directly.
            !------------------------------------------------------------------------------
            tlcl_res(crs_counter) = tout

            CALL write_criteria_space_to_vtk(&
                out%p_n_bsnm//".tcr."//dmn_no//"."//TRIM(temp_suf)//"."//vtk_suf)

            crs_counter = crs_counter + 1_mik

        END DO

        CALL MPI_SEND(tlcl_res, INT(crs, KIND=mik), MPI_tensor_2nd_rank_R66, 0_mik, &
            INT(tout%dmn, KIND=mik), MPI_COMM_WORLD, ierr)
        CALL print_err_stop(std_out, "MPI_SEND on tlcl_res failed.", INT(ierr, KIND=ik))

        ! CALL MPI_SEND(tlcl_res%dmn, INT(crs, KIND=mik), MPI_INTEGER8, 0_mik, tin%dmn, &
        !     MPI_COMM_WORLD, ierr)
        ! CALL print_err_stop(std_out, "MPI_SEND on tlcl_res%dmn failed.", INT(ierr, KIND=ik))

        ! CALL MPI_SEND(tlcl_res%density, INT(crs, KIND=mik), MPI_DOUBLE_PRECISION, 0_mik, tin%dmn, &
        !     MPI_COMM_WORLD, ierr)
        ! CALL print_err_stop(std_out, "MPI_SEND on tlcl_res%density failed.", INT(ierr, KIND=ik))

        ! CALL MPI_SEND(tlcl_res%doa_zener, INT(crs, KIND=mik), MPI_DOUBLE_PRECISION, 0_mik, tin%dmn, &
        !     MPI_COMM_WORLD, ierr)
        ! CALL print_err_stop(std_out, "MPI_SEND on tlcl_res%doa_zener failed.", INT(ierr, KIND=ik))

        ! CALL MPI_SEND(tlcl_res%doa_gebert, INT(crs, KIND=mik), MPI_DOUBLE_PRECISION, 0_mik, tin%dmn, &
        !     MPI_COMM_WORLD, ierr)
        ! CALL print_err_stop(std_out, "MPI_SEND on tlcl_res%doa_gebert failed.", INT(ierr, KIND=ik))

        ! CALL MPI_SEND(tlcl_res%pos, INT(crs, KIND=mik)*3_mik, MPI_DOUBLE_PRECISION, 0_mik, tin%dmn, &
        !     MPI_COMM_WORLD, ierr)
        ! CALL print_err_stop(std_out, "MPI_SEND on tlcl_res%mat failed.", INT(ierr, KIND=ik))

        ! CALL MPI_SEND(tlcl_res%mat, INT(crs, KIND=mik)*36_mik, MPI_DOUBLE_PRECISION, 0_mik, tin%dmn, &
        !     MPI_COMM_WORLD, ierr)
        ! CALL print_err_stop(std_out, "MPI_SEND on tlcl_res%mat failed.", INT(ierr, KIND=ik))

    END DO

END IF ! Worker processes since "ELSE"

IF(my_rank == 0) THEN
    !------------------------------------------------------------------------------
    ! Write data to files.
    ! Crs counter mandatory! Otherwise, jj and tglbl_res will not match if not all
    ! kinds of optimizations are requested. 
    !------------------------------------------------------------------------------
    crs_counter = 0_mik
    
    DO jj = 1_ik, 4_ik
        IF(.NOT. execute_optimization(jj)) CYCLE

        crs_counter = 1_mik
    
        SELECT CASE(jj)
            CASE(1); fhwcrit = fh_mono; suffix = suf_mono
            CASE(2); fhwcrit = fh_orth; suffix = suf_orth
            CASE(3); fhwcrit = fh_ani1; suffix = suf_ani1
            CASE(4); fhwcrit = fh_ani2; suffix = suf_ani2
        END SELECT    

        CALL write_tensor_2nd_rank_R66(fhwcrit, covo_amnt_lines, tglbl_res(crs_counter, :))                    
        CALL meta_stop_ascii(fhwcrit, suffix)

    END DO

    CALL meta_write(fhmeo, 'ZERO_MATRICES', '(-)', zero_matrix_counter)

    CALL meta_signing(binary)
END IF ! (my_rank == 0)

!------------------------------------------------------------------------------
! Only used in specific cases to finish more gracefully. (grep -i "GOTO")
!------------------------------------------------------------------------------
1001 Continue

IF(my_rank == 0) THEN

    !------------------------------------------------------------------------------
    ! Finish the program
    !------------------------------------------------------------------------------
    CALL meta_close()

    CALL meta_stop_ascii(fh_log, log_suf)
    IF(debug >= 0) CALL meta_stop_ascii(fh_mon, mon_suf)


    IF(std_out/=6) CALL meta_stop_ascii(fh=std_out, suf='.std_out')

    CALL CPU_TIME(end)

    WRITE(std_out,FMT_MSG_SEP)
    WRITE(std_out,FMT_TXT_AF0A) 'Overall Time = ', (end-start) / 60._rk ,' Minutes'
    WRITE(std_out,FMT_TXT_AF0A) 'CPU time = ', (end-start) / 60._rk / 60._rk * size_mpi,' Hours'
    WRITE(std_out,FMT_MSG_SEP)

    WRITE(std_out,FMT_TXT_AF0A) 'Program finished successfully.'
    WRITE(std_out,FMT_MSG_SEP)

END IF ! (my_rank == 0)

CALL MPI_FINALIZE(ierr)
CALL print_err_stop(std_out, "MPI_FINALIZE failed.", INT(ierr, KIND=ik))

END PROGRAM tensor_optimizer

