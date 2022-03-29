
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
SUBROUTINE write_criteria_space_to_vtk(file, steps)

    CHARACTER(LEN=*), INTENT(IN) :: file
    INTEGER(KIND=ik), DIMENSION(3) :: steps

    INTEGER(KIND=ik), DIMENSION(3) :: ttl_steps
    INTEGER(KIND=ik) :: fh
    LOGICAL :: fex

    fex = .FALSE.

    INQUIRE(FILE=TRIM(file), EXIST=fex)
    IF(fex) WRITE(std_out, FMT_WRN) TRIM(file)//" already exists."

    fh = give_new_unit()

    ttl_steps = (steps*2_ik)+1_ik

    CALL write_vtk_struct_points_header(fh, file, TRIM('rk4'), &
        [1._rk, 1._rk, 1._rk], [0._rk, 0._rk, 0._rk], ttl_steps)

    CALL ser_write_raw(fh, file, REAL(crit, KIND=REAL32), 'BIG_ENDIAN')

    CALL write_vtk_struct_points_footer(fh, file)

END SUBROUTINE write_criteria_space_to_vtk

!------------------------------------------------------------------------------
! FUNCTION: stop_workers
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Stop all processs gracefully.
!
!> @param[in] size_mpi Amount of processes 
!> @return stat Status
!------------------------------------------------------------------------------  
FUNCTION stop_workers(size_mpi) RESULT (ierr)

    INTEGER(KIND=mik), INTENT(IN) :: size_mpi

    INTEGER(KIND=mik) :: ierr, zz

    !------------------------------------------------------------------------------
    ! Stop workers properly by sending an activity = -1 info
    !------------------------------------------------------------------------------
    DO zz = 1, size_mpi-1_mik
        CALL MPI_SEND(-1_mik, 1_mik, MPI_INTEGER, zz, zz, MPI_COMM_WORLD, ierr)
        CALL print_err_stop(std_out, "MPI_SEND of stopping workers failed.", INT(ierr, KIND=ik))
    END DO

END FUNCTION stop_workers



END MODULE auxiliaries_of_tensor_optimizer

!------------------------------------------------------------------------------
! PROGRAM: Rotational Tensor Optimization
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> \date 27.10.2020
!> \modified 07.01.2022
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
INTEGER(KIND=ik), PARAMETER :: debug = 0
REAL   (KIND=rk), PARAMETER :: lower_thres_percentage = 0.01

! Variables
TYPE(tensor_2nd_rank_R66) :: dummy
TYPE(tensor_2nd_rank_R66), DIMENSION(:), ALLOCATABLE :: tglbl_in
TYPE(tensor_2nd_rank_R66), DIMENSION(:,:), ALLOCATABLE :: tglbl_res
TYPE(tensor_2nd_rank_R66), DIMENSION(:), ALLOCATABLE :: tlcl_res

TYPE(materialcard) :: bone

CHARACTER(LEN=mcl), DIMENSION(:), ALLOCATABLE :: m_rry      
CHARACTER(LEN=mcl) :: crit_file, cmd_arg_history=''
CHARACTER(LEN=scl) :: restart, restart_cmd_arg, re_mono, re_orth, re_ani1, re_ani2
CHARACTER(LEN=scl) :: suf_covo, suf_mono, suf_orth, suf_ani1, suf_ani2, temp_suf    
CHARACTER(LEN=scl) :: binary, dmn_no, suffix
CHARACTER(LEN=  1) :: stg='1'
CHARACTER(LEN=  8) :: date
CHARACTER(LEN= 10) :: time

REAL(KIND=rk) :: start, end, sym
REAL(KIND=rk), DIMENSION(2,3) :: intervall 

INTEGER(KIND=ik) :: fh_covo, fh_mono, fh_orth, fh_ani1, fh_ani2, fhwcrit
INTEGER(KIND=ik) :: exp_dmn_crit, covo_amnt_lines, zero_matrix_counter = 0
INTEGER(KIND=ik) :: jj, kk, mm, stat, invalid_entries, iostat
INTEGER(KIND=ik), DIMENSION(2,3) :: steps

LOGICAL, DIMENSION(4) :: execute_optimization = .FALSE.
LOGICAL :: stp, print_criteria, fex, crit_exp_req = .FALSE., dmn_found = .FALSE.

INTEGER(KIND=mik) :: ierr, my_rank, size_mpi, mii
INTEGER(KIND=mik) :: active, feed_ranks, crs, crs_counter
INTEGER(KIND=mik), DIMENSION(MPI_STATUS_SIZE) :: stmpi
INTEGER(KIND=mik), DIMENSION(:,:), ALLOCATABLE :: statuses_mpi
INTEGER(KIND=mik), DIMENSION(:), ALLOCATABLE :: req_list

INTEGER(KIND=mik) :: MPI_tensor_2nd_rank_R66
INTEGER(KIND=mik), DIMENSION(7) :: blocklen, dtype 
INTEGER(KIND=MPI_ADDRESS_KIND) :: disp(7), base

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
CALL MPI_GET_ADDRESS(dummy%sym, disp(5), ierr) 
CALL MPI_GET_ADDRESS(dummy%pos, disp(6), ierr) 
CALL MPI_GET_ADDRESS(dummy%mat, disp(7), ierr) 
	
base = disp(1) 
disp = disp - base 

blocklen(1:5) = 1 
blocklen(6) = 3 
blocklen(7) = 36 

dtype(1) = MPI_INTEGER8 
dtype(2:7) = MPI_DOUBLE_PRECISION

CALL MPI_TYPE_CREATE_STRUCT(7_mik, blocklen, disp, dtype, MPI_tensor_2nd_rank_R66, ierr) 
CALL mpi_err(ierr,"MPI_tensor_2nd_rank_R66 couldn't be created.")

CALL MPI_TYPE_COMMIT(MPI_tensor_2nd_rank_R66, ierr)
CALL mpi_err(ierr,"MPI_tensor_2nd_rank_R66 couldn't be commited.")

!------------------------------------------------------------------------------
! Initialize program itself
!------------------------------------------------------------------------------
IF(my_rank == 0) THEN
    CALL CPU_TIME(start)

    !------------------------------------------------------------------------------
    ! Parse the command arguments
    !------------------------------------------------------------------------------
    CALL get_cmd_args(binary, in%full, stp, restart_cmd_arg, cmd_arg_history)
    IF(stp) GOTO 1001
    
    IF (in%full=='') THEN
        CALL usage(binary)    

        !------------------------------------------------------------------------------
        ! On std_out since file of std_out is not spawned
        !------------------------------------------------------------------------------
        CALL print_err_stop(6, "No input file given", 1)
    END IF

    !------------------------------------------------------------------------------
    ! Check and open the input file; Modify the Meta-Filename / Basename
    ! Define the new application name first
    !------------------------------------------------------------------------------
    global_meta_prgrm_mstr_app = 'rot' 
    global_meta_program_keyword = 'ROT_TENSOR_OPT'
    CALL meta_append(m_rry, size_mpi)
    
    !------------------------------------------------------------------------------
    ! Redirect std_out into a file in case std_out is not useful by environment.
    ! Place these lines before handle_lock_file :-)
    !------------------------------------------------------------------------------
    std_out = determine_stout()

    !------------------------------------------------------------------------------
    ! Spawn standard out after(!) the basename is known
    !------------------------------------------------------------------------------
    IF(std_out/=6) CALL meta_start_ascii(std_out, '.std_out')

    CALL show_title(["Johannes Gebert, M.Sc. (HLRS, NUM)"])

    IF(debug >=0) WRITE(std_out, FMT_MSG) "Post mortem info probably in ./datasets/temporary.std_out"

    !------------------------------------------------------------------------------
    ! Parse input
    !------------------------------------------------------------------------------
    CALL meta_read('REQUEST_MONO'   , m_rry, re_mono)
    CALL meta_read('REQUEST_ORTH'   , m_rry, re_orth)
    CALL meta_read('REQUEST_ANI_1'  , m_rry, re_ani1)
    CALL meta_read('REQUEST_ANI_2'  , m_rry, re_ani2)
    CALL meta_read('RESTART'        , m_rry, restart)
    CALL meta_read('EXPORT_DMN_CRIT', m_rry, exp_dmn_crit)
    CALL meta_read('YOUNG_MODULUS'  , m_rry, bone%E)
    CALL meta_read('POISSON_RATIO'  , m_rry, bone%nu)

    !------------------------------------------------------------------------------
    ! Restart handling
    ! Done after meta_io to decide based on keywords
    !------------------------------------------------------------------------------
    CALL meta_handle_lock_file(restart, restart_cmd_arg)

    !------------------------------------------------------------------------------
    ! Start monitoring of MPI
    !------------------------------------------------------------------------------
    CALL meta_start_ascii(fh_mon, mon_suf)

    CALL DATE_AND_TIME(date, time)

    WRITE(std_out, FMT_TXT_SEP)  
    WRITE(std_out, FMT_TXT) TRIM(ADJUSTL(longname))//" Results"
    WRITE(std_out, FMT_TXT) "Date: "//date//" [ccyymmdd]"
    WRITE(std_out, FMT_TXT) "Time: "//time//" [hhmmss.sss]"
    WRITE(std_out, FMT_TXT) "Program invocation:"//TRIM(cmd_arg_history)          
    WRITE(std_out, FMT_TXT_SEP)  
    
    !------------------------------------------------------------------------------
    ! Create/Open tensor files. Basically tuned csv data.
    !------------------------------------------------------------------------------
    fh_covo = give_new_unit()
    suf_covo = ".stage-0.covo" ! control volume (in situ orientation)
    CALL meta_existing_ascii(fh_covo, suf_covo, covo_amnt_lines)

    crs = 0_mik
    IF(TRIM(re_mono) == "YES") THEN
        fh_mono = give_new_unit()
        suf_mono = ".stage-2.mono" ! monotropic optimization 
        CALL meta_start_ascii(fh_mono, suf_mono) 
        execute_optimization(1) = .TRUE. 
        crs=crs+1_mik
    END IF

    IF(TRIM(re_orth) == "YES") THEN
        fh_orth = give_new_unit()
        suf_orth = ".stage-2.orth" ! orthotropic optimization 
        CALL meta_start_ascii(fh_orth, suf_orth) 
        execute_optimization(2) = .TRUE. 
        crs=crs+1_mik
    END IF

    IF(TRIM(re_ani1) == "YES") THEN
        fh_ani1 = give_new_unit()
        suf_ani1 = ".stage-2.an1" ! anisotropic optimization, version 1 
        CALL meta_start_ascii(fh_ani1, suf_ani1) 
        execute_optimization(3) = .TRUE. 
        crs=crs+1_mik
    END IF

    IF(TRIM(re_ani2) == "YES") THEN
        fh_ani2 = give_new_unit()
        suf_ani2 = ".stage-2.an2" ! anisotropic optimization, version 2 
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

    !------------------------------------------------------------------------------
    ! Export domain if requested number is positive (and matches an entry in file)
    ! Export all domains if number requested = -10 (and matrix /= 0._rk)
    ! Export nothing if number ((<0_ik) .AND. (/= 10_ik))
    !------------------------------------------------------------------------------
    IF((exp_dmn_crit == -10_ik) .OR. (exp_dmn_crit >= 0_ik))  crit_exp_req = .TRUE.

    IF((exp_dmn_crit == -10_ik) .AND. (covo_amnt_lines * crs > 100_ik)) THEN
        mssg = "Are you sure to export the criteria spaces of more than 100 &
            &optimizations (dmns*criterias) to the file system via *.vtk files?! &
            &If yes, change the source and recompile :-)"
        CALL print_err_stop(std_out, mssg, 1)
    END IF

    IF(debug >= 0) THEN
        WRITE(std_out, FMT_DBG_xAI0) "Debug Level:", debug
        WRITE(std_out, FMT_DBG_xAI0) "Processors:", size_mpi  
        WRITE(std_out, FMT_DBG_xAI0) "Amount of domains:", covo_amnt_lines-1_ik
        WRITE(std_out, FMT_TXT_SEP)
        FLUSH(std_out)
    END IF

    !------------------------------------------------------------------------------
    ! Parse the global input data
    ! Allocate amnt_lines -1_ik due to header
    !------------------------------------------------------------------------------
    ALLOCATE(tglbl_in(covo_amnt_lines-1_ik))
    CALL parse_tensor_2nd_rank_R66(fh_covo, in%p_n_bsnm//TRIM(suf_covo), &
        covo_amnt_lines, tglbl_in, invalid_entries)

    IF(covo_amnt_lines == invalid_entries+1_ik) THEN
        
        stat = stop_workers(size_mpi)
        
        mssg = "No valid data found. Probably an implementation issue or an &
            &invalid file format."
        CALL print_err_stop(std_out, mssg, 1)
        
    END IF
    !------------------------------------------------------------------------------
    ! Allocate global results array. Need to loop over execute_optimization 
    ! to properly assign the tglbl_res(entries, :)
    !------------------------------------------------------------------------------
    ALLOCATE(tglbl_res(crs, covo_amnt_lines-1_ik))

END IF ! (my_rank == 0)

!------------------------------------------------------------------------------
! MPI derived datatypes are another way of sending the meta information.
! But only the p_n_bsnm/path are used and they are used rarely. Therefore not
! considered as a major issue.
!------------------------------------------------------------------------------
CALL MPI_BCAST( in%p_n_bsnm, INT(meta_mcl, KIND=mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(out%path    , INT(meta_mcl, KIND=mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(out%p_n_bsnm, INT(meta_mcl, KIND=mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST(restart, 1_mik, MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)

!------------------------------------------------------------------------------
! Variables are global ones, which are valid for each specific rank (!)
!
! The optimization is done in two steps. First, search the best orientation
! by 1° steps. Second, search at this positon with 0.025° steps.
!------------------------------------------------------------------------------
! Shorten turnaround times
!------------------------------------------------------------------------------
intervall(1,:) = 1._rk
intervall(2,:) = 0.025_rk

steps(1,:) = 182_ik
steps(2,:) = 80_ik

IF(debug >= 3) THEN
    steps(1,:) = 30_ik
    steps(2,:) = 30_ik
END IF

!------------------------------------------------------------------------------
! Why 8? Only 4 entries in array. 
! MPI_LOGICAL = 4 Byte, Fortran_Logical = 8 Byte?
!------------------------------------------------------------------------------
CALL MPI_BCAST(execute_optimization, 8_mik, MPI_LOGICAL, 0_mik, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST(exp_dmn_crit, 1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(crs, 1_mik, MPI_INTEGER, 0_mik, MPI_COMM_WORLD, ierr)

!------------------------------------------------------------------------------
! All Ranks -- Init MPI request and status lists
!------------------------------------------------------------------------------
ALLOCATE(req_list(size_mpi-1_mik))
req_list=1_mik

ALLOCATE(statuses_mpi(MPI_STATUS_SIZE, size_mpi-1_mik))
statuses_mpi=0_mik

!------------------------------------------------------------------------------
! Rank 0 -- Process master - Start working process
!------------------------------------------------------------------------------
IF (my_rank==0) THEN

    IF(debug >= 0) THEN
        WRITE(std_out, FMT_DBG) "Starting worker supply."
        WRITE(std_out, FMT_DBG_xAF0) "Young modulus: ", bone%E
        WRITE(std_out, FMT_DBG_xAF0) "Poissions ratio: ", bone%nu
        WRITE(std_out, FMT_DBG_xAF0) "Mat out threshold: ", lower_thres_percentage * bone%E
        WRITE(std_out, FMT_DBG_SEP)
        FLUSH(std_out)
    END IF

    mii = 1_mik
    feed_ranks = 1_mik 
        
    IF(debug >= 2) THEN
        ! DO mii = 1, 9
        WRITE(std_out, FMT_DBG_AI0AxI0) "tglbl_in(" ,mii, ")%dmn: ", tglbl_in(mii)%dmn
        WRITE(std_out, FMT_DBG_AI0AxF0) "tglbl_in(" ,mii, ")%density: ", tglbl_in(mii)%density
        WRITE(std_out, FMT_DBG_AI0AxF0) "tglbl_in(" ,mii, ")%doa_zener: ", tglbl_in(mii)%doa_zener
        WRITE(std_out, FMT_DBG_AI0AxF0) "tglbl_in(" ,mii, ")%doa_gebert: ", tglbl_in(mii)%doa_gebert
        WRITE(std_out, FMT_DBG_AI0AxF0) "tglbl_in(" ,mii, ")%pos: ", tglbl_in(mii)%pos
        WRITE(std_out, FMT_DBG_AI0AxF0) "tglbl_in(" ,mii, ")%mat(1,1:3): ", tglbl_in(mii)%mat(1,1:3)
        WRITE(std_out, FMT_TXT_SEP)
        FLUSH(std_out)
        ! end do
    END IF

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
        ! Track whether the domain requested is contained in the input data.
        !------------------------------------------------------------------------------
        IF(exp_dmn_crit >= 0_ik) THEN
            IF(tglbl_in(mii)%dmn == exp_dmn_crit) dmn_found = .TRUE.
        END IF

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
 
            mii = mii + 1_mik
    
            CYCLE
        END IF

        !------------------------------------------------------------------------------
        ! Give feedback to user. Only active if the machine provides an interactive 
        ! terminal, which is determined by the environmemnt.
        !------------------------------------------------------------------------------)
        IF(std_out == 6_ik) THEN
            CALL EXECUTE_COMMAND_LINE("clear")

            CALL show_title(["Johannes Gebert, M.Sc. (HLRS, NUM)"])

            WRITE(std_out, FMT_TXT_xAI0) "Processed domains: ", mii, " of ", covo_amnt_lines-1_ik
            WRITE(std_out, FMT_TXT_xAI0) "Most current input to compute:", tglbl_in(mii)%dmn
            WRITE(std_out, FMT_TXT_xAI0) "Total 0-matrices: ", zero_matrix_counter
            WRITE(std_out, FMT_TXT) ""

            WRITE(dmn_no, '(I0)') tglbl_in(mii)%dmn ! Write the domain number to string (!)

            CALL write_matrix(std_out, "Domain "//TRIM(dmn_no), tglbl_in(mii)%mat, 'spl')

            WRITE (std_out, FMT_TXT) "Lots of optimization steps requested."
            WRITE (std_out, FMT_TXT) "The computation may take a long time."
                
        END IF

        !------------------------------------------------------------------------------
        ! Amount of entities to compute. 
        ! -1_ik due to header of file.
        !------------------------------------------------------------------------------
        IF (feed_ranks > (size_mpi-1_mik)) feed_ranks = 1_ik 

        !------------------------------------------------------------------------------
        ! Call the domain an active one.
        !------------------------------------------------------------------------------
        CALL MPI_SEND(1_mik, 1_mik, MPI_INTEGER, feed_ranks, feed_ranks, MPI_COMM_WORLD,ierr)
        CALL print_err_stop(std_out, "MPI_SEND of activity didn't succeed", INT(ierr, KIND=ik))

        !------------------------------------------------------------------------------
        ! Send tensor
        !------------------------------------------------------------------------------
        CALL MPI_SEND(tglbl_in(mii), 1_mik, MPI_tensor_2nd_rank_R66, feed_ranks, &
            feed_ranks, MPI_COMM_WORLD, ierr)
        CALL print_err_stop(std_out, "MPI_SEND of tin failed.", INT(ierr, KIND=ik))            
         
        !------------------------------------------------------------------------------
        ! Log to monitor file (first worker thread)
        !------------------------------------------------------------------------------
        IF(debug >= 0) THEN
            WRITE(fh_mon, DBG//"3(A,1x,"//FMT_INT//",1x),A,3(I5,1x),A,3(F8.4,1x))") &
                "MPI rank: ", feed_ranks, &
                " Domain number: ", tglbl_in(mii)%dmn
            FLUSH(fh_mon)
        END IF
        !------------------------------------------------------------------------------
        ! Place a receive for the tensor. 
        ! Tensors are not sortded automatically.
        !------------------------------------------------------------------------------
        CALL MPI_IRECV(tglbl_res(:, mii), crs, MPI_tensor_2nd_rank_R66, feed_ranks, &
            INT(tglbl_in(mii)%dmn, KIND=mik), MPI_COMM_WORLD, req_list(feed_ranks), ierr)
        CALL print_err_stop(std_out, "MPI_IRECV of activity(mii) failed.", INT(ierr, KIND=ik))

        feed_ranks = feed_ranks + 1_ik  

        mii = mii + 1_mik
    END DO
    ! Still my_rank == 0

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
        CALL MPI_RECV(active, 1_mik, MPI_INTEGER, 0_mik, my_rank, MPI_COMM_WORLD, stmpi, ierr)
        CALL print_err_stop(std_out, "MPI_RECV on active failed.", INT(ierr, KIND=ik))

        IF(active == -1) GOTO 1001

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

        IF((exp_dmn_crit == tin%dmn) .OR. (exp_dmn_crit == -10_ik)) print_criteria = .TRUE.

        WRITE(dmn_no, '(I0)') tin%dmn

        !------------------------------------------------------------------------------
        ! Start computation of domains by workers.
        !------------------------------------------------------------------------------
        crs_counter = 1_mik

        DO jj = 1_ik, 4_ik
            IF(.NOT. execute_optimization(jj)) CYCLE

            IF(debug >= 3) THEN
                WRITE(std_out, FMT_DBG_AI0AxI0) "Rank ", my_rank, " tin%dmn: ", tin%dmn
                WRITE(std_out, FMT_DBG_AI0AxF0) "Rank ", my_rank, " tin%density: ", tin%density
                WRITE(std_out, FMT_DBG_AI0AxF0) "Rank ", my_rank, " tin%doa_zener: ", tin%doa_zener
                WRITE(std_out, FMT_DBG_AI0AxF0) "Rank ", my_rank, " tin%doa_gebert: ", tin%doa_gebert
                WRITE(std_out, FMT_DBG_AI0AxF0) "Rank ", my_rank, " tin%pos: ", tin%pos
                WRITE(std_out, FMT_DBG_AI0AxF0) "Rank ", my_rank, " tin%mat(1,1:3): ", tin%mat(1,1:3)
                WRITE(std_out, FMT_TXT_SEP)
                FLUSH(std_out)
            END IF

            !------------------------------------------------------------------------------
            ! Optimization is capable of accepting dig /= 0._rk. Tensors can be optimized
            ! several times and do not need to be restarted from the position of the
            ! control volume within the bone.
            !
            ! Assign a domain specific variable to a global, currently validone
            !------------------------------------------------------------------------------
            dig = tin%pos 

            DO kk = 1_ik, 2_ik
                IF(debug >= 3) THEN
                    WRITE(std_out, FMT_DBG_AI0AxI0) "Rank ", my_rank, " tin%dmn: ", tin%dmn
                    WRITE(std_out, FMT_DBG_AI0AxI0) "Rank ", my_rank, " Optimization step: ", kk
                    WRITE(std_out, FMT_DBG_AI0AxI0) "Rank ", my_rank, " Optimization case: ", jj
                    WRITE(std_out, FMT_TXT_SEP)
                    FLUSH(std_out)
                END IF

                !------------------------------------------------------------------------------
                ! Reset for stage 2
                ! Dig = best position of stage 1. Dig is a global variable (!)
                ! Optimization always begins at dig - (intervall * steps / 2._rk)
                !------------------------------------------------------------------------------
                IF(kk == 2_ik) THEN
                    tin%mat = tout%mat
                    dig = tout%pos
                END IF
                
                !------------------------------------------------------------------------------
                ! Optimize tensors
                !------------------------------------------------------------------------------
                SELECT CASE(jj)
                    CASE(1)
                        CALL opt_stiff('monotropic'  , steps(kk,:), intervall(kk,:))
                        temp_suf = ".mono"
                    CASE(2)
                        CALL opt_stiff('orthotropic' , steps(kk,:), intervall(kk,:))
                        temp_suf = ".orth"
                    CASE(3)
                        CALL opt_stiff('anisotropic1', steps(kk,:), intervall(kk,:))
                        temp_suf = ".an1"
                    CASE(4)
                        CALL opt_stiff('anisotropic2', steps(kk,:), intervall(kk,:))
                        temp_suf = ".an2"
                END SELECT        

                !------------------------------------------------------------------------------
                ! Tilts until S11 > S22 > S33 
                !------------------------------------------------------------------------------
                CALL tilt_tensor(tout%mat)
                CALL check_sym(tout%mat, sym)

                tout%dmn = tin%dmn
                tout%doa_zener = doa_zener(tout%mat)
                tout%doa_gebert = doa_gebert(tout%mat)
                tout%density = gebert_density_voigt(tout%mat, bone%E, bone%nu)
                tout%sym = sym

                !------------------------------------------------------------------------------
                ! Print vtk files of criteria spaces
                !------------------------------------------------------------------------------
                IF(print_criteria) THEN 
                    WRITE(stg, '(I0)') kk

                    crit_file = TRIM(out%p_n_bsnm)//".stage-"//stg//"."//TRIM(dmn_no)//TRIM(temp_suf)//vtk_suf

                    INQUIRE(FILE=TRIM(crit_file), EXIST=fex)

                    IF(fex) THEN
                        IF(restart == 'Y') THEN
                            CALL EXECUTE_COMMAND_LINE("rm -f "//TRIM(crit_file), CMDSTAT=iostat)

                            CALL print_err_stop(std_out, &
                                "Removing the file "//TRIM(crit_file)//" failed.", iostat)
                        ELSE
                            CALL print_err_stop(std_out, &
                                "File "//TRIM(crit_file)//" exists and restart is 'N'.", iostat)
                        END IF
                    END IF

                    CALL write_criteria_space_to_vtk(TRIM(crit_file), steps(kk,:))
                END IF

            END DO
            !------------------------------------------------------------------------------
            ! At the end of the second step, the results acutally get written to the 
            ! output variables that are then send back to my_rank=0 / main process.
            !
            ! Domain always gets retrieved from tin directly.
            !------------------------------------------------------------------------------
            tlcl_res(crs_counter) = tout

            crs_counter = crs_counter + 1_mik

        END DO

        CALL MPI_SEND(tlcl_res, INT(crs, KIND=mik), MPI_tensor_2nd_rank_R66, 0_mik, &
            INT(tout%dmn, KIND=mik), MPI_COMM_WORLD, ierr)
        CALL print_err_stop(std_out, "MPI_SEND on tlcl_res failed.", INT(ierr, KIND=ik))

    END DO

END IF ! Worker processes since "ELSE"

IF(my_rank == 0) THEN

    !------------------------------------------------------------------------------
    ! Print last information
    !------------------------------------------------------------------------------
    WRITE(std_out, FMT_TXT_xAI0) "Processors:", size_mpi  
    WRITE(std_out, FMT_TXT_SEP)  
    WRITE(std_out, FMT_TXT) "Optimization parameters:"
    WRITE(std_out, FMT_TXT_AxI0) "Steps of stage 1:", steps(1,:)  
    WRITE(std_out, FMT_TXT_AxF0) "Intervall of stage 1:", intervall(1,:)  
    WRITE(std_out, FMT_TXT) ""
    WRITE(std_out, FMT_TXT_AxI0) "Steps of stage 2:", steps(2,:)  
    WRITE(std_out, FMT_TXT_AxF0) "Intervall of stage 2:", intervall(2,:)  
    WRITE(std_out, FMT_TXT_SEP)  

    !------------------------------------------------------------------------------
    ! Wait for all processes
    !------------------------------------------------------------------------------
    CALL MPI_WAITALL(size_mpi-1_mik, req_list, statuses_mpi, ierr)
    CALL print_err_stop(std_out, "MPI_WAITANY on req_list for IRECV of tglbl_res failed.", &
        INT(ierr, KIND=ik))

    stat = stop_workers(size_mpi)

    !------------------------------------------------------------------------------
    ! Write data to files. Crs -> criteria space counter
    ! Crs counter mandatory! Otherwise, jj and tglbl_res will not match if not all
    ! kinds of optimizations are requested. 
    !------------------------------------------------------------------------------
    crs_counter = 0_mik
    
    DO jj = 1_ik, 4_ik
        IF(.NOT. execute_optimization(jj)) CYCLE

        crs_counter = crs_counter + 1_mik
    
        SELECT CASE(jj)
            CASE(1); fhwcrit = fh_mono; suffix = suf_mono
            CASE(2); fhwcrit = fh_orth; suffix = suf_orth
            CASE(3); fhwcrit = fh_ani1; suffix = suf_ani1
            CASE(4); fhwcrit = fh_ani2; suffix = suf_ani2
        END SELECT    

        CALL write_tensor_2nd_rank_R66(fhwcrit, covo_amnt_lines, tglbl_res(crs_counter, :))                    
        CALL meta_stop_ascii(fhwcrit, suffix)

    END DO

    CALL meta_stop_ascii(fh_covo, suf_covo)

    CALL meta_write('ZERO_MATRICES', '(-)', zero_matrix_counter)

    CALL meta_signing(binary)
END IF ! (my_rank == 0)

!------------------------------------------------------------------------------
! Only used in specific cases to finish more gracefully. (grep -i "GOTO")
!------------------------------------------------------------------------------
1001 Continue

!------------------------------------------------------------------------------
! Stop monitoring of rank 1
!------------------------------------------------------------------------------
IF((debug >= 0) .AND. (my_rank == 1)) CALL meta_stop_ascii(fh_mon, mon_suf)

IF(my_rank == 0) THEN

    !------------------------------------------------------------------------------
    ! Finish the program
    !------------------------------------------------------------------------------
    CALL meta_close()

    CALL CPU_TIME(end)

    WRITE(std_out,FMT_TXT_xAF0) 'Overall Time = ', (end-start) / 60._rk ,' Minutes'
    WRITE(std_out,FMT_TXT_xAF0) 'CPU time = ', (end-start) / 60._rk / 60._rk * size_mpi,' Hours'
    WRITE(std_out,FMT_TXT_SEP)

    IF((.NOT. dmn_found) .AND. (crit_exp_req) .AND. (exp_dmn_crit /= -10_ik)) THEN
        WRITE(std_out,FMT_TXT) "Domain number for crit export requested, but no &
        &corresponding tensor found:"
        WRITE(std_out,FMT_TXT_xAI0) "Domain ", exp_dmn_crit
        WRITE(std_out,FMT_TXT) "If problem was on OSI-Layer 8, then give a &
            &*.stage-0.covo file with only the domain(s) requested."
        WRITE(std_out,FMT_TXT) ""      
    END IF
    
    IF(.NOT. crit_exp_req) THEN
        WRITE(std_out,FMT_TXT) "No domain criteria space export requested."
        WRITE(std_out,FMT_TXT) ""      
    END IF

    WRITE(std_out,FMT_TXT) 'Program finished.'
    WRITE(std_out,FMT_TXT_SEP)

    IF(std_out/=6) CALL meta_stop_ascii(fh=std_out, suf='.std_out')

END IF ! (my_rank == 0)

CALL MPI_FINALIZE(ierr)
CALL print_err_stop(std_out, "MPI_FINALIZE failed.", INT(ierr, KIND=ik))

END PROGRAM tensor_optimizer

