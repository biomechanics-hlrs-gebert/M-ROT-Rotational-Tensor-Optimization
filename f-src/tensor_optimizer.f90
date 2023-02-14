
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
USE mpi_binary
USE ser_binary
USE opt_stiffness
USE mpi_user_interaction

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

    CHARACTER(*), INTENT(IN) :: file
    INTEGER(ik), DIMENSION(3) :: steps

    INTEGER(ik), DIMENSION(3) :: ttl_steps
    INTEGER(ik) :: fh
    LOGICAL :: fex

    fex = .FALSE.

    INQUIRE(FILE=TRIM(file), EXIST=fex)
    IF(fex) WRITE(std_out, FMT_WRN) TRIM(file)//" already exists."

    fh = give_new_unit()

    ttl_steps = (steps*2_ik)+1_ik

    CALL write_vtk_struct_points_header(fh, file, TRIM('rk4'), &
        [1._rk, 1._rk, 1._rk], [0._rk, 0._rk, 0._rk], ttl_steps)

    CALL ser_write_binary(fh, file, REAL(crit, REAL32), 'BIG_ENDIAN')

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

    INTEGER(mik), INTENT(IN) :: size_mpi

    INTEGER(mik) :: ierr, zz

    !------------------------------------------------------------------------------
    ! Stop workers properly by sending an activity = -1 info
    !------------------------------------------------------------------------------
    DO zz = 1, size_mpi-1_mik
        CALL MPI_SEND(-1_mik, 1_mik, MPI_INTEGER, zz, zz, MPI_COMM_WORLD, ierr)
        CALL print_err_stop(std_out, "MPI_SEND of stopping workers failed.", INT(ierr, ik))
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

!------------------------------------------------------------------------------
! Parameter - They *MUST* fit to each other!
!------------------------------------------------------------------------------
CHARACTER(LEN=5), PARAMETER :: lower_thres_percentage_char = "<3% E"
REAL(rk), PARAMETER :: lower_thres_percentage = 0.03
REAL(rk), PARAMETER :: exec_thres = 0.0001

! Variables
TYPE(domain_data) :: dummy
TYPE(domain_data), DIMENSION(:), ALLOCATABLE :: tglbl_in, tlcl_res
TYPE(domain_data), DIMENSION(:,:), ALLOCATABLE :: tglbl_res

TYPE(materialcard) :: bone

CHARACTER(mcl), DIMENSION(:), ALLOCATABLE :: m_rry      
CHARACTER(mcl) :: crit_file, cmd_arg_history='', stat=""
CHARACTER(scl) :: restart, restart_cmd_arg, binary, dmn_no, suf_files(4), suf_covo
CHARACTER(  1) :: stg='1'
CHARACTER(  8) :: date
CHARACTER( 10) :: time

REAL(rk) :: start, end, res_mono, res_orth, res_ani1, res_ani2
REAL(rk), DIMENSION(2) :: step_width, swept_range
REAL(rk), DIMENSION(3) :: dmn_size, spcng
REAL(rk), DIMENSION(4) :: exec_opt

INTEGER(ik) :: fh_covo, no_stages, fh_files(4), &
    exp_dmn_crit, covo_no_lines, zero_matrix_counter = 0, &
    ii, jj, kk, mm, xx, invalid_entries, iostat, steps(2,3)
INTEGER(ik), DIMENSION(3) :: dims, grid

LOGICAL :: print_criteria, fex, crit_exp_req = .FALSE., dmn_found = .FALSE.

INTEGER(mik) :: ierr, my_rank, size_mpi, mii, active, feed_ranks, crs, crs_counter, ios
INTEGER(mik), DIMENSION(MPI_STATUS_SIZE) :: stmpi
INTEGER(mik), DIMENSION(:), ALLOCATABLE :: req_list, statInt

INTEGER(mik) :: MPI_tensor_2nd_rank_R66
INTEGER(mik), DIMENSION(24) :: blocklen, dtype 
INTEGER(MPI_ADDRESS_KIND) :: disp(24), base


LOGICAL :: abrt = .FALSE.

! Initialize MPI Environment
CALL MPI_INIT(ierr)
CALL mpi_err(ierr,"MPI_INIT failed.")

CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
CALL mpi_err(ierr,"MPI_COMM_RANK couldn't be retrieved")

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)
CALL mpi_err(ierr,"MPI_COMM_SIZE couldn't be retrieved")

IF (size_mpi < 2) THEN
    CALL print_err_stop(std_out, "At least two ranks required to execute this program.", 1)
END IF 

ALLOCATE(statInt(size_mpi))


!------------------------------------------------------------------------------
! Redirect std_out into a file in case std_out is not useful by environment.
! Place these lines before handle_lock_file :-)
!------------------------------------------------------------------------------
CALL MPI_GET_ADDRESS(dummy%section       , disp( 1), ierr)  ! Position within the CT image
CALL MPI_GET_ADDRESS(dummy%dmn           , disp( 2), ierr)  ! Number of the control volume
CALL MPI_GET_ADDRESS(dummy%no_elems      , disp( 3), ierr)  ! Number of elements of the domain
CALL MPI_GET_ADDRESS(dummy%no_nodes      , disp( 4), ierr)  ! Number of nodes of the domain
CALL MPI_GET_ADDRESS(dummy%collected_logs, disp( 5), ierr)  ! Timestamps during domain computation
CALL MPI_GET_ADDRESS(dummy%dmn_size      , disp( 6), ierr)  ! Size of the control volume
CALL MPI_GET_ADDRESS(dummy%t_start       , disp( 7), ierr)  ! Start of the domain after program start (s)
CALL MPI_GET_ADDRESS(dummy%t_duration    , disp( 8), ierr)  ! Duration of the computation of the domain (s)
CALL MPI_GET_ADDRESS(dummy%phy_dmn_bnds  , disp( 9), ierr)  ! Physical domain boundaries (x,y,z - lo,hi)
CALL MPI_GET_ADDRESS(dummy%opt_res       , disp(10), ierr)  ! Resolution the covo was optimized with
CALL MPI_GET_ADDRESS(dummy%pos           , disp(11), ierr)  ! Position (deg) of alpha, eta, phi
CALL MPI_GET_ADDRESS(dummy%sym           , disp(12), ierr)  ! Symmetry deviation (quotient)
CALL MPI_GET_ADDRESS(dummy%DA            , disp(13), ierr)  ! Degree of anisotropy - Bone gold standard
CALL MPI_GET_ADDRESS(dummy%bvtv          , disp(14), ierr)  ! Bone volume/total volume
CALL MPI_GET_ADDRESS(dummy%gray_density  , disp(15), ierr)  ! Density based on grayscale values.
CALL MPI_GET_ADDRESS(dummy%doa_zener     , disp(16), ierr)  ! Zener degree of anisotropy
CALL MPI_GET_ADDRESS(dummy%doa_gebert    , disp(17), ierr)  ! Gebert degree of anisotropy (modified Zener)
CALL MPI_GET_ADDRESS(dummy%mps           , disp(18), ierr)  ! Mean principal stiffness
CALL MPI_GET_ADDRESS(dummy%spec_norm     , disp(19), ierr)  ! Gebert degree of anisotropy (modified Zener)
CALL MPI_GET_ADDRESS(dummy%num           , disp(20), ierr)  ! An actual numerical stiffness tensor
CALL MPI_GET_ADDRESS(dummy%mat           , disp(21), ierr)  ! An actual stiffness tensor
CALL MPI_GET_ADDRESS(dummy%opt_crit      , disp(22), ierr)  ! Optimization information (e.g. criteria)
CALL MPI_GET_ADDRESS(dummy%ma_el_type    , disp(23), ierr)  ! Optimization information (e.g. criteria)
CALL MPI_GET_ADDRESS(dummy%mi_el_type    , disp(24), ierr)  ! Optimization information (e.g. criteria)

base = disp(1) 
disp = disp - base 

blocklen(1)  = 3_mik
blocklen(2)  = 1_mik
blocklen(3)  = 1_mik
blocklen(4)  = 1_mik
blocklen(5)  = 24_mik
blocklen(6)  = 1_mik
blocklen(7)  = 1_mik
blocklen(8)  = 1_mik
blocklen(9)  = 6_mik
blocklen(10) = 1_mik
blocklen(11) = 3_mik
blocklen(12) = 1_mik
blocklen(13) = 1_mik
blocklen(14) = 1_mik
blocklen(15) = 1_mik
blocklen(16) = 1_mik 
blocklen(17) = 1_mik 
blocklen(18) = 1_mik 
blocklen(19) = 1_mik 
blocklen(20) = 576_mik 
blocklen(21) = 36_mik 
blocklen(22) = INT(scl, mik)
blocklen(23) = INT(scl, mik)
blocklen(24) = INT(scl, mik)

dtype(1:5)   = MPI_INTEGER8 
dtype(6:21)  = MPI_DOUBLE_PRECISION
dtype(22:24) = MPI_CHARACTER 

CALL MPI_TYPE_CREATE_STRUCT(24_mik, blocklen, disp, dtype, MPI_tensor_2nd_rank_R66, ierr) 
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
    CALL get_cmd_args(binary, in%full, restart_cmd_arg, cmd_arg_history)
    
    IF (in%full=='') THEN
        CALL usage(binary)    

        !------------------------------------------------------------------------------
        ! On std_out since file of std_out is not spawned
        !------------------------------------------------------------------------------
        CALL print_err_stop(6, "No input file given", 1)
    END IF

    CALL show_title(["Johannes Gebert, M.Sc. (HLRS, NUM)"])

    !------------------------------------------------------------------------------
    ! Check and open the input file; Modify the Meta-Filename / Basename
    ! Define the new application name first
    !------------------------------------------------------------------------------
    global_meta_prgrm_mstr_app = 'rot'
    global_meta_program_keyword = 'ROT_TENSOR_OPT'
    CALL meta_append(m_rry, size_mpi, binary, stat)

    !------------------------------------------------------------------------------
    ! Redirect std_out into a file in case std_out is not useful by environment.
    ! Place these lines before handle_lock_file :-)
    !------------------------------------------------------------------------------
    ! CALL determine_std_fh(std_out, std_err)
    ! write(*,*) "std_out: ", std_out
    ! write(*,*) "std_err: ", std_err

    !------------------------------------------------------------------------------
    ! Spawn standard out after(!) the basename is known
    !------------------------------------------------------------------------------
    ! IF(std_out/=6) CALL meta_start_ascii(std_out, '.std_out')
    ! IF(std_out/=0) CALL meta_start_ascii(std_out, '.std_err')

    IF(out_amount=="DEBUG") THEN
        WRITE(std_out, FMT_MSG) "Post mortem info probably in ./datasets/temporary.std_out"
    END IF 

    !------------------------------------------------------------------------------
    ! Parse input
    !------------------------------------------------------------------------------
    CALL meta_read('RESOLUTION_MONO', m_rry, res_mono, stat); CALL mest(stat, abrt)
    CALL meta_read('RESOLUTION_ORTH', m_rry, res_orth, stat); CALL mest(stat, abrt)
    CALL meta_read('RESOLUTION_ANI1', m_rry, res_ani1, stat); CALL mest(stat, abrt)
    CALL meta_read('RESOLUTION_ANI2', m_rry, res_ani2, stat); CALL mest(stat, abrt)

    CALL meta_read('EXPORT_DMN_CRIT', m_rry, exp_dmn_crit, stat); CALL mest(stat, abrt)

    CALL meta_read('RESTART'        , m_rry, restart, stat); CALL mest(stat, abrt)

    CALL meta_read('YOUNG_MODULUS'  , m_rry, bone%E ,  stat); CALL mest(stat, abrt)
    CALL meta_read('POISSON_RATIO'  , m_rry, bone%nu,  stat); CALL mest(stat, abrt)
    CALL meta_read('SIZE_DOMAIN'    , m_rry, dmn_size, stat); CALL mest(stat, abrt)
    CALL meta_read('DIMENSIONS'     , m_rry, dims,     stat); CALL mest(stat, abrt)
    CALL meta_read('SPACING'        , m_rry, spcng,    stat); CALL mest(stat, abrt)

    !------------------------------------------------------------------------------
    ! Restart handling
    ! Done after meta_io to decide based on keywords
    !------------------------------------------------------------------------------
    CALL meta_handle_lock_file(restart, restart_cmd_arg)

    CALL DATE_AND_TIME(date, time)

    IF(out_amount=="DEBUG") THEN
        WRITE(std_out, FMT_TXT_SEP)  
        WRITE(std_out, FMT_TXT) TRIM(ADJUSTL(longname))//" Results"
        WRITE(std_out, FMT_TXT) "Date: "//date//" [ccyymmdd]"
        WRITE(std_out, FMT_TXT) "Time: "//time//" [hhmmss.sss]"
        WRITE(std_out, FMT_TXT) "Program invocation:"//TRIM(cmd_arg_history)          
        WRITE(std_out, FMT_TXT_SEP)  
    END IF

    !------------------------------------------------------------------------------
    ! Create/Open tensor files. Basically tuned csv data.
    !------------------------------------------------------------------------------
    fh_covo = give_new_unit()
    suf_covo = ".covo" ! control volume (in situ orientation)
    CALL meta_existing_ascii(fh_covo, suf_covo, covo_no_lines)

    !------------------------------------------------------------------------------
    ! Package user request on optimization
    !------------------------------------------------------------------------------
    exec_opt(1) = res_mono; suf_files(1) = "mono" 
    exec_opt(2) = res_orth; suf_files(2) = "orth" 
    exec_opt(3) = res_ani1; suf_files(3) = "an1" 
    exec_opt(4) = res_ani2; suf_files(4) = "an2" 

    crs = 0_mik
    DO ii=1, 4

        fh_files(ii) = give_new_unit()
        
        !------------------------------------------------------------------------------
        ! Exec thresh to check whether the angle ~= 0.
        !------------------------------------------------------------------------------
        IF(exec_opt(ii) > exec_thres) THEN
            crs=crs+1_mik
            
            CALL meta_start_ascii(fh_files(ii), "."//TRIM(suf_files(ii))) 
        END IF
    END DO

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

    IF((exp_dmn_crit == -10_ik) .AND. (covo_no_lines * crs > 100_ik)) THEN
        mssg = "Are you sure to export the criteria spaces of more than 100 &
            &optimizations (dmns*criterias) to the file system via *.vtk files?! &
            &If yes, change the source and recompile :-)"
        CALL print_err_stop(std_out, mssg, 1)
    END IF

    IF(out_amount=="DEBUG") THEN
        WRITE(std_out, FMT_DBG_xAI0) "Processors:", size_mpi  
        WRITE(std_out, FMT_DBG_xAI0) "Number of domains:", covo_no_lines-1_ik
        WRITE(std_out, FMT_TXT_SEP)
        FLUSH(std_out)
    END IF

    !------------------------------------------------------------------------------
    ! Parse the global input data
    ! Allocate amnt_lines -1_ik due to header
    !------------------------------------------------------------------------------
    ALLOCATE(tglbl_in(covo_no_lines-1_ik))

    CALL parse_tensor_2nd_rank_R66(fh_covo, TRIM(in%p_n_bsnm)//TRIM(suf_covo), &
        covo_no_lines, tglbl_in, invalid_entries)

    IF (invalid_entries > 0) THEN
        WRITE(std_out, FMT_WRN_AI0AxI0) "", invalid_entries, " invalid domains."
    END IF

    IF(covo_no_lines == invalid_entries+1_ik) THEN
        
        statInt = stop_workers(size_mpi)
        
        mssg = "No valid data found. Probably an implementation issue or an &
            &invalid file format."
        CALL print_err_stop(std_out, mssg, 1)
        
    END IF
    !------------------------------------------------------------------------------
    ! Allocate global results array. Need to loop over exec_opt 
    ! to properly assign the tglbl_res(entries, :)
    !------------------------------------------------------------------------------
    ALLOCATE(tglbl_res(crs, covo_no_lines-1_ik))

END IF ! (my_rank == 0)

!------------------------------------------------------------------------------
! MPI derived datatypes are another way of sending the meta information.
! But only the p_n_bsnm/path are used and they are used rarely. Therefore not
! considered as a major issue.
!------------------------------------------------------------------------------
CALL MPI_BCAST( in%p_n_bsnm, INT(meta_mcl, mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(out%path    , INT(meta_mcl, mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(out%p_n_bsnm, INT(meta_mcl, mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(suf_files   , INT(scl*4_ik, mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(restart, 1_mik, MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(exec_opt, 4_mik, MPI_DOUBLE_PRECISION, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(exp_dmn_crit, 1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(crs, 1_mik, MPI_INTEGER, 0_mik, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST(dmn_size, 3_mik, MPI_DOUBLE_PRECISION, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(spcng,    3_mik, MPI_DOUBLE_PRECISION, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(dims,     3_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)

!------------------------------------------------------------------------------
! All Ranks -- Init MPI request and status lists
!------------------------------------------------------------------------------
ALLOCATE(req_list(size_mpi-1_mik))
req_list=0_mik

!------------------------------------------------------------------------------
! Rank 0 -- Process master - Start working process
!------------------------------------------------------------------------------
IF (my_rank==0) THEN

    IF(out_amount=="DEBUG") THEN
        WRITE(std_out, FMT_DBG) "Starting worker supply."
        WRITE(std_out, FMT_DBG_xAF0) "Young modulus: ", bone%E
        WRITE(std_out, FMT_DBG_xAF0) "Poissions ratio: ", bone%nu
        WRITE(std_out, FMT_DBG_xAF0) "Mat out threshold: ", lower_thres_percentage * bone%E
        WRITE(std_out, FMT_DBG_SEP)
        FLUSH(std_out)
    END IF

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
    DO WHILE (mii <= covo_no_lines-1_ik)

        !------------------------------------------------------------------------------
        ! Track whether the domain requested is contained in the input data.
        !------------------------------------------------------------------------------
        IF(exp_dmn_crit >= 0_ik) THEN
            IF(tglbl_in(mii)%dmn == exp_dmn_crit) dmn_found = .TRUE.
        END IF

        IF(out_amount=="DEBUG") THEN

        END IF 

        !------------------------------------------------------------------------------
        ! Check whether it is a zero tensor:
        !------------------------------------------------------------------------------
        IF(MAXVAL(tglbl_in(mii)%mat) <= lower_thres_percentage * bone%E) THEN
            DO mm=1, crs
                tglbl_res(mm, mii) = tglbl_in(mii)

                tglbl_res(mm, mii)%opt_crit = lower_thres_percentage_char
                tglbl_res(mm, mii)%mat  = 0._rk

            END DO
            zero_matrix_counter = zero_matrix_counter + 1_ik
 
            mii = mii + 1_mik
    
            CYCLE
        ELSE ! Is not a zero tensor - "standard"
            tglbl_res(:, mii)%opt_crit = "std."
        END IF

        !------------------------------------------------------------------------------
        ! Give feedback to user. Only active if the machine provides an interactive 
        ! terminal, which is determined by the environmemnt.
        !------------------------------------------------------------------------------)
        IF((std_out == 6_ik) .AND. (out_amount == "DEBUG")) THEN
            CALL EXECUTE_COMMAND_LINE("clear")

            CALL show_title(["Johannes Gebert, M.Sc. (HLRS, NUM)"])

            WRITE(std_out, FMT_TXT_xAI0) "Processed domains: ", mii, " of ", covo_no_lines-1_ik
            WRITE(std_out, FMT_TXT_xAI0) "Most current input to compute:", tglbl_in(mii)%dmn
            WRITE(std_out, FMT_TXT_xAI0) "Total 0-matrices: ", zero_matrix_counter
            WRITE(std_out, FMT_TXT) ""

            WRITE(dmn_no, '(I0)') tglbl_in(mii)%dmn ! Write the domain number to string (!)

            write(std_out,FMT_TXT_AxI0) "tglbl_in(mii)%section(3)       : ", tglbl_in(mii)%section(:)       
            write(std_out,FMT_TXT_AxI0) "tglbl_in(mii)%dmn              : ", tglbl_in(mii)%dmn              
            write(std_out,FMT_TXT_AxF0) "tglbl_in(mii)%dmn_size         : ", tglbl_in(mii)%dmn_size         
            write(std_out,FMT_TXT_AxF0) "tglbl_in(mii)%phy_dmn_bnds(3,2): ", tglbl_in(mii)%phy_dmn_bnds(:,:)
            write(std_out,FMT_TXT_AxF0) "tglbl_in(mii)%opt_res          : ", tglbl_in(mii)%opt_res          
            write(std_out,FMT_TXT_AxF0) "tglbl_in(mii)%pos(3)           : ", tglbl_in(mii)%pos(:)           
            write(std_out,FMT_TXT_AxF0) "tglbl_in(mii)%sym              : ", tglbl_in(mii)%sym              
            write(std_out,FMT_TXT_AxF0) "tglbl_in(mii)%DA               : ", tglbl_in(mii)%DA               
            write(std_out,FMT_TXT_AxF0) "tglbl_in(mii)%bvtv             : ", tglbl_in(mii)%bvtv             
            write(std_out,FMT_TXT_AxF0) "tglbl_in(mii)%gray_density     : ", tglbl_in(mii)%gray_density     
            write(std_out,FMT_TXT_AxF0) "tglbl_in(mii)%doa_zener        : ", tglbl_in(mii)%doa_zener        
            write(std_out,FMT_TXT_AxF0) "tglbl_in(mii)%doa_gebert       : ", tglbl_in(mii)%doa_gebert     
            write(std_out,FMT_TXT_AxF0) "tglbl_in(mii)%mps              : ", tglbl_in(mii)%mps  
            write(std_out,FMT_TXT_AxF0) "TRIM(tglbl_in(mii)%opt_crit)   : "//TRIM(tglbl_in(mii)%opt_crit)
            
            CALL write_matrix(std_out, "Domain "//TRIM(dmn_no), tglbl_in(mii)%mat, 'spl')

            write(std_out,FMT_SEP)

            IF(out_amount /= "DEBUG") THEN
                WRITE (std_out, FMT_TXT) "Lots of optimization steps requested."
                WRITE (std_out, FMT_TXT) "The computation may take a long time."
            END IF 
        END IF

        !------------------------------------------------------------------------------
        ! Number of entities to compute. 
        ! -1_ik due to header of file.
        !------------------------------------------------------------------------------
        IF (feed_ranks > (size_mpi-1_mik)) feed_ranks = 1_ik 

        !------------------------------------------------------------------------------
        ! Call the domain an active one.
        !------------------------------------------------------------------------------
        CALL MPI_SEND(1_mik, 1_mik, MPI_INTEGER, feed_ranks, feed_ranks, MPI_COMM_WORLD,ierr)
        CALL print_err_stop(std_out, "MPI_SEND of activity didn't succeed", INT(ierr, ik))

        !------------------------------------------------------------------------------
        ! Send tensor
        !------------------------------------------------------------------------------
        CALL MPI_SEND(tglbl_in(mii), 1_mik, MPI_tensor_2nd_rank_R66, &
            feed_ranks, feed_ranks, MPI_COMM_WORLD, ierr)
        CALL print_err_stop(std_out, "MPI_SEND of tin failed.", INT(ierr, ik))            
         
        !------------------------------------------------------------------------------
        ! Place a receive for the tensor. 
        ! Tensors are not sortded automatically.
        !------------------------------------------------------------------------------
        CALL MPI_IRECV(tglbl_res(:, mii), crs, MPI_tensor_2nd_rank_R66, feed_ranks, &
            INT(tglbl_in(mii)%dmn, mik), MPI_COMM_WORLD, req_list(feed_ranks), ierr)
        CALL print_err_stop(std_out, "MPI_IRECV of activity(mii) failed.", INT(ierr, ik))

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
        ! Stop if all domains are computed
        ! Stop slaves in case of an error
        !------------------------------------------------------------------------------
        CALL MPI_RECV(active, 1_mik, MPI_INTEGER, 0_mik, my_rank, MPI_COMM_WORLD, stmpi, ierr)
        CALL print_err_stop(std_out, "MPI_RECV on active failed.", INT(ierr, ik))

        IF(active == -1) GOTO 1001
        
        !------------------------------------------------------------------------------
        ! Receive tensor (derived type of tensor_2nd_rank_R66)
        !------------------------------------------------------------------------------
        CALL MPI_RECV(tin, 1_mik, MPI_tensor_2nd_rank_R66, 0_mik, my_rank, MPI_COMM_WORLD, stmpi, ierr)
        CALL print_err_stop(std_out, "MPI_RECV on tin failed.", INT(ierr, ik))
   
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

        !------------------------------------------------------------------------------
        ! Loop over modes
        !------------------------------------------------------------------------------
        mm=0_ik
        DO jj = 1_ik, 4_ik

            !------------------------------------------------------------------------------
            ! Cycle if criteria is not requested. Is as good as any other if/else.
            !------------------------------------------------------------------------------
            IF(exec_opt(jj) < exec_thres) CYCLE

            !------------------------------------------------------------------------------
            ! The optimization is done in two steps if the chosen step size is too small.
            ! If the user chooses a large step size, then one stage is sufficient.
            ! Done this way to shorten turnaround times
            !------------------------------------------------------------------------------
            IF (exec_opt(jj) < 2._rk) THEN
                no_stages = 2_ik
                
                swept_range(1) = 180._rk
                swept_range(2) = 2._rk

                step_width(1) = 2._rk
                step_width(2) = exec_opt(jj)
            ELSE
                no_stages = 1_ik

                ! All entries assigned to help preventing unintended situations.
                swept_range(1) = 180._rk
                swept_range(2) = 180._rk

                step_width(1) = exec_opt(jj)
                step_width(2) = exec_opt(jj)
            END IF 

            !----------------------------------------------------------------------------------------------
            ! Copy input to output tensor. Important to keep domain data etc.
            !----------------------------------------------------------------------------------------------
            tout = tin      

            !------------------------------------------------------------------------------
            ! Optimization is capable of accepting dig /= 0._rk. Tensors can be optimized
            ! several times and do not need to be restarted from the position of the
            ! control volume within the bone.
            !
            ! Assign a domain specific variable to a global, currently validone
            !------------------------------------------------------------------------------
            dig = tin%pos 

            DO kk = 1_ik, no_stages

                !------------------------------------------------------------------------------
                ! Reset for stage 2
                ! Dig = best position of stage 1. Dig is a global variable (!)
                ! Optimization always begins at dig - (intervall * steps / 2._rk)
                !------------------------------------------------------------------------------
                IF(kk == 2_ik) THEN
                    !------------------------------------------------------------------------------
                    ! DO NOT assign tin%mat = tout%mat since dig=tout%pos already is an absolute
                    ! angle, refering to tin%mat (global coordinate system of the bone.)
                    !------------------------------------------------------------------------------
                    !!! tin%mat = tout%mat
                    dig = tout%pos
                END IF
                
                !------------------------------------------------------------------------------
                ! Optimize tensors
                ! Implementation should ensure positive deltas in exec_opt to step along 
                ! the range with. Otherwise, small ranges and moving backwards can cause
                ! unintended results/behavior.
                !------------------------------------------------------------------------------
                CALL opt_stiff(suf_files(jj), swept_range(kk), step_width(kk))

                !------------------------------------------------------------------------------
                ! Tilts until S11 > S22 > S33 
                !------------------------------------------------------------------------------
                CALL tilt_tensor(tout%mat)
                
                !------------------------------------------------------------------------------
                ! Print vtk files of criteria spaces
                !------------------------------------------------------------------------------
                IF((print_criteria) .AND. (exp_dmn_crit >= 0_ik))THEN 
                    WRITE(stg, '(I0)') kk

                    crit_file = TRIM(out%p_n_bsnm)//".stage-"//stg//"."//TRIM(dmn_no)//"."//TRIM(suf_files(jj))//vtk_suf

                    INQUIRE(FILE=TRIM(crit_file), EXIST=fex)

                    IF(fex) THEN
                        WRITE(std_out, '(A)') ""
                        WRITE(std_out, FMT_WRN) "Deleting existing *.vtk file."

                        CALL EXECUTE_COMMAND_LINE("rm -f "//TRIM(crit_file), CMDSTAT=iostat)

                        CALL print_err_stop(std_out, &
                            "Removing the file "//TRIM(crit_file)//" failed.", iostat)
                    END IF

                    CALL write_criteria_space_to_vtk(TRIM(crit_file), steps(kk,:))
                END IF

            END DO



            !------------------------------------------------------------------------------
            ! Compute additional parameters.
            !------------------------------------------------------------------------------
            grid = get_grid(dmn_size, dims, spcng)
            tout%section = domain_no_to_section(tout%dmn, grid)
          
            DO xx=1, 3
                tout%phy_dmn_bnds(xx,1) =  tout%section(xx)   *tout%dmn_size
                tout%phy_dmn_bnds(xx,2) = (tout%section(xx)+1)*tout%dmn_size
            END DO

            tout%sym = check_sym(tout%mat)
            tout%mps = mps(tout%mat)
            tout%doa_zener = doa_zener(tout%mat)
            tout%doa_gebert = doa_gebert(tout%mat)
            tout%opt_res = exec_opt(jj)

            ! CALL spectral_norm(tout%mat, 6_mik, tout%spec_norm, ios)

            IF(ios /= 0_ik) WRITE(std_out, FMT_ERR_AI0AxI0) &
                    "Computing the Eigenvalue for domain ", tout%dmn, " failed."

            !------------------------------------------------------------------------------
            ! At the end of the second step, the results acutally get written to the 
            ! output variables that are then send back to my_rank=0 / main process.
            !
            ! Domain always gets retrieved from tin directly.
            !------------------------------------------------------------------------------
            tlcl_res(crs_counter) = tout

            crs_counter = crs_counter + 1_mik

        END DO

        CALL MPI_SEND(tlcl_res, INT(crs, mik), MPI_tensor_2nd_rank_R66, 0_mik, &
            INT(tout%dmn, mik), MPI_COMM_WORLD, ierr)
        CALL print_err_stop(std_out, "MPI_SEND on tlcl_res failed.", INT(ierr, ik))

    END DO

END IF ! Worker processes since "ELSE"

IF(my_rank == 0) THEN

    !------------------------------------------------------------------------------
    ! Wait for all processes
    !------------------------------------------------------------------------------
    CALL MPI_WAITALL(size_mpi-1_mik, req_list, MPI_STATUSES_IGNORE, ierr)
    CALL print_err_stop(std_out, "MPI_WAITANY on req_list for IRECV of tglbl_res failed.", &
        INT(ierr, ik))
    
    statInt = stop_workers(size_mpi)

    !------------------------------------------------------------------------------
    ! Write data to files. Crs -> criteria space counter
    ! Crs counter mandatory! Otherwise, jj and tglbl_res will not match if not all
    ! kinds of optimizations are requested. 
    !------------------------------------------------------------------------------
    crs_counter = 0_mik
    
    DO jj = 1_ik, 4_ik
        
        IF(exec_opt(jj) <= exec_thres) CYCLE

        crs_counter = crs_counter + 1_mik
        
        CALL write_tensor_2nd_rank_R66(fh_files(jj), covo_no_lines, tglbl_res(crs_counter, :))                    
        CALL meta_stop_ascii(fh_files(jj), "."//TRIM(suf_files(jj)))

    END DO

    CALL meta_stop_ascii(fh_covo, suf_covo)

    CALL meta_write('ZERO_MATRICES', '(-)', zero_matrix_counter)
END IF ! (my_rank == 0)

!------------------------------------------------------------------------------
! Only used in specific cases to finish more gracefully. (grep -i "GOTO")
!------------------------------------------------------------------------------
1001 Continue


IF(my_rank == 0) THEN

    !------------------------------------------------------------------------------
    ! Finish the program
    ! ------------------------------------------------------------------------------
    CALL meta_close(m_rry)

    CALL CPU_TIME(end)

    WRITE(std_out, FMT_TXT_xAF0) 'Overall Time = ', (end-start) / 60._rk ,' Minutes'
    WRITE(std_out, FMT_TXT_xAF0) 'CPU time = ', (end-start) / 60._rk / 60._rk * size_mpi,' Hours'
    WRITE(std_out, FMT_TXT_SEP)

    IF((.NOT. dmn_found) .AND. (crit_exp_req) .AND. (exp_dmn_crit /= -10_ik)) THEN
        WRITE(std_out,FMT_TXT) "Domain number for crit export requested, but no &
            &corresponding tensor found:"
        WRITE(std_out,FMT_TXT_xAI0) "Domain ", exp_dmn_crit
        WRITE(std_out,FMT_TXT) "If problem was on OSI-Layer 8, then give a &
            &*.stage-0.covo file with only the domain(s) requested."
        WRITE(std_out,FMT_TXT) ""      
    END IF
    
    IF(.NOT. crit_exp_req) THEN
        WRITE(std_out, FMT_TXT) "No domain criteria space export requested."
        WRITE(std_out, FMT_TXT) ""      
    END IF

    WRITE(std_out, FMT_TXT) 'Program finished.'
    WRITE(std_out, FMT_TXT_SEP)

END IF ! (my_rank == 0)

CALL MPI_FINALIZE(ierr)
CALL print_err_stop(std_out, "MPI_FINALIZE failed.", INT(ierr, ik))

END PROGRAM tensor_optimizer

