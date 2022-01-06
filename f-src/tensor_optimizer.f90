
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
USE opt_stiffness
USE vtk_meta_data
USE raw_binary

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

	REAL(KIND=rk), DIMENSION(3) :: spcng = 0._rk, origin = 0._rk
	INTEGER(KIND=ik), DIMENSION(3) :: ttl_steps
	INTEGER(KIND=ik) :: fh_temp
	LOGICAL :: fex

	fex = .FALSE.

	INQUIRE(FILE=TRIM(file), EXIST=fex)
	IF(fex) WRITE(std_out, FMT_WRN) TRIM(file)//" already exists."

	fh_temp = give_new_unit()

	ttl_steps = (pm_steps*2_ik)+1_ik

	CALL write_vtk_struct_points_header(fh_temp, file, TRIM('rk8'), &
		spcng, origin, ttl_steps)

	CALL ser_write_raw(fh_temp, file, crit)

	CALL write_vtk_struct_points_footer(fh_temp, file)

END SUBROUTINE write_criteria_space_to_vtk

END MODULE auxiliaries_of_tensor_optimizer

!------------------------------------------------------------------------------
! PROGRAM: Rotational Tensor Optimization
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> \date 27.10.2020
!> \modified 05.01.2022
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

  INTEGER (KIND=ik) :: zero_matrix_counter, stat
  REAL    (KIND=rk) :: percent

		!-- Allgemeine Laufvariablen um numerischen Kernel zu steuern
		INTEGER     (KIND=ik) :: ii_lines, last_ii_lns_rlvnt, &
			ll, mm,  oo, pp, vali, zz
		REAL        (KIND=rk), DIMENSION(6,6) :: mono_opt, orth_opt, mono_sym, orth_sym, M_Null
		INTEGER     (KIND=ik), PARAMETER      :: header_row_header=1
		CHARACTER   (LEN=200)                 :: input_file, header
		CHARACTER   (LEN=1000)                :: line, CR_str
		CHARACTER   (LEN=12)                  :: string
		INTEGER     (KIND=ik), PARAMETER      :: n_cols=38, header_row=1
		INTEGER     (KIND=4)                  :: lines, lunit=10, how_many_lines=10
		CHARACTER   (LEN=10)                  :: hw_mny_lns
		INTEGER     (KIND=ik), DIMENSION(4)   :: un
		INTEGER  (KIND=ik)                    :: ios, ntokens
		CHARACTER(len=mcl)                    :: tokens(100)


		

		TYPE(csv_data), DIMENSION(:,:) , ALLOCATABLE :: t2nd                  ! Tensor of 2nd order R6x6
		CHARACTER(LEN=200), DIMENSION(4) :: flnm
		CHARACTER(LEN=*)  , PARAMETER :: int_fmt='(I5)'
		CHARACTER(LEN=*)  , PARAMETER :: REAL_fmt='(F9.3)'
		CHARACTER(LEN=*)  , PARAMETER :: REAL_int='(F9.0)'
		REAL     (KIND=rk)            :: start, finish
		INTEGER  (KIND=ik)            :: ETAt

		!-- Optimization variables
		REAL(KIND=rk) :: a,p,e
		REAL(KIND=rk) :: Sum_DoA_Mono=0, Sum_DoA_Orth=0, Sum_DoA_InO=0, Sum_DoA_InM=0
		REAL(KIND=rk) :: Avg_DoA_Mono, Avg_DoA_Orth, Avg_DoA_InO, Avg_DoA_InM
		REAL(KIND=rk) :: am1, pm1, em1, am2, pm2, em2, ao1, po1, eo1, ao2, po2, eo2

		!-- Parameters
		INTEGER     (KIND=ik), PARAMETER :: EY=5600
		REAL        (KIND=rk), PARAMETER :: factor_low_thres=0.001
		CHARACTER   (LEN=3)              :: calculat
		INTEGER     (KIND=ik)            :: calcfail=0
		LOGICAL                          :: rprt_mat, rprt

		!-- DEFAULTS ----- rprt=.TRUE. ----- rprt_mat=.FALSE. -----------------------------------------
		rprt=.TRUE.
		rprt_mat=.FALSE.
		!-- DEFAULTS ----------------------------------------------------------------------------------


TYPE(tensor_2nd_rank_R66), DIMENSION(:), ALLOCATABLE :: tglbl_in
TYPE(tensor_2nd_rank_R66), DIMENSION(:,:), ALLOCATABLE :: tglbl_res
TYPE(tensor_2nd_rank_R66), DIMENSION(:), ALLOCATABLE :: tlcl_res

TYPE(materialcard) :: bone

CHARACTER(LEN=mcl), DIMENSION(:), ALLOCATABLE :: m_rry      
CHARACTER(LEN=scl) :: restart, restart_cmd_arg, re_mono, re_orth, re_ani1, re_ani2
CHARACTER(LEN=scl) :: suf_covo, suf_mono, suf_orth, suf_ani1, suf_ani2, temp_suf	
CHARACTER(LEN=  8) :: date
CHARACTER(LEN= 10) :: time

INTEGER(KIND=mik) :: ierr, my_rank, size_mpi, mii, index, request, active
INTEGER(KIND=mik), DIMENSION(MPI_STATUS_SIZE) :: stmpi
INTEGER(KIND=mik), DIMENSION(:,:), ALLOCATABLE :: statuses_mpi
INTEGER(KIND=ik), DIMENSION(:), ALLOCATABLE :: activity, req_list

INTEGER(KIND=ik) :: fh_covo, fh_mono, fh_orth, fh_ani1, fh_ani2
INTEGER(KIND=ik) :: exp_dmn_crit, covo_amnt_lines, crs, crs_counter
INTEGER(KIND=ik) :: jj, kk, nn

LOGICAL, DIMENSION(4) :: execute_optimization = .FALSE.
LOGICAL :: stp, print_criteria



! Initialize MPI Environment
CALL MPI_INIT(ierr)
CALL MPI_ERR(ierr,"MPI_INIT failed.")

CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
CALL MPI_ERR(ierr,"MPI_COMM_RANK couldn't be retrieved")

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)
CALL MPI_ERR(ierr,"MPI_COMM_SIZE couldn't be retrieved")

IF (size_mpi < 2) CALL print_err_stop(std_out, "At least two ranks required to execute this program.", 1)

! Initialize program itself
IF(my_rank == 0) THEN

    !------------------------------------------------------------------------------
    ! Redirect std_out into a file in case std_out is not useful by environment.
    ! Place these lines before handle_lock_file :-)
    !------------------------------------------------------------------------------
    std_out = determine_stout()

    IF(std_out/=6) CALL meta_start_ascii(std_out, '.std_out')

    CALL show_title()
 
    IF(debug >=0) WRITE(*,FMT_MSG) "Post mortem info probably in ./datasets/.temporary.std_out"

    CALL CPU_TIME(global_start)

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
    global_meta_program_keyword = 'ROTATIONAL_TENSOR_OPTIMIZATION'
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

	!------------------------------------------------------------------------------
    ! Restart handling
    ! Done after meta_io to decide based on keywords
    !------------------------------------------------------------------------------
    CALL meta_handle_lock_file(restart, restart_cmd_arg)

	!------------------------------------------------------------------------------
	! Spawn a log and a monitoring file (to check the behavior of MPI)
	!------------------------------------------------------------------------------
    CALL meta_start_ascii(fh_log, log_suf)
    CALL meta_start_ascii(fh_mon, mon_suf)

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
    suf_covo = ".tcr.covo" ! control volume (in situ orientation)
	suf_mono = ".tcr.mono" ! monotropic optimization
	suf_orth = ".tcr.orth" ! orthotropic optimization
	suf_ani1 = ".tcr.an1" ! anisotropic optimization, version 1
	suf_ani2 = ".tcr.an2" ! anisotropic optimization, version 2

	fh_covo = give_new_unit()
	fh_mono = give_new_unit()
	fh_orth = give_new_unit()
	fh_ani1 = give_new_unit()
	fh_ani2 = give_new_unit()

	CALL meta_existing_ascii(fh_covo, suf_covo, covo_amnt_lines)

	crs = 0_ik
	IF(TRIM(re_mono) == "YES") execute_optimization(1) = .TRUE.; crs=crs+1_ik
	IF(TRIM(re_orth) == "YES") execute_optimization(2) = .TRUE.; crs=crs+1_ik
	IF(TRIM(re_ani1) == "YES") execute_optimization(3) = .TRUE.; crs=crs+1_ik
	IF(TRIM(re_ani2) == "YES") execute_optimization(4) = .TRUE.; crs=crs+1_ik

    !------------------------------------------------------------------------------
    ! Ensure, the user is sane and report if so (or not).
    !------------------------------------------------------------------------------
	IF(crs == 0_ik) THEN
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
	
	!------------------------------------------------------------------------------
    ! Generate internal logging of activities.
    ! Only worker processes involved (!)
    !------------------------------------------------------------------------------
	ALLOCATE(activity(size_mpi-1)
	activity=1

END IF ! (my_rank == 0)

CALL MPI_BCAST( in%p_n_bsnm, INT(meta_mcl, KIND=mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(out%p_n_bsnm, INT(meta_mcl, KIND=mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST(execute_optimization, 4_mik, MPI_LOGICAL, 0_mik, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST(exp_dmn_crit, 1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(crs, 1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)

!------------------------------------------------------------------------------
! All Ranks -- Init MPI request and status lists
!------------------------------------------------------------------------------
ALLOCATE(req_list(size_mpi-1_mik))
req_list=0_ik

ALLOCATE(statuses_mpi(MPI_STATUS_SIZE, size_mpi-1_mik))
statuses_mpi=0

activity = 1_ik

!------------------------------------------------------------------------------
! Rank 0 -- Process master - Start working process
!------------------------------------------------------------------------------
IF (my_rank==0) THEN

	nn = 1_ik
	mii = 1_mik

	!------------------------------------------------------------------------------
	! Supply all worker masters with their first tensor
	! Master does not compute tensors
	!
	! Basically, the master loops over all workers to send and receive from them
	!------------------------------------------------------------------------------
	DO WHILE (mii <= (size_mpi-1_mik))

		!------------------------------------------------------------------------------
		! Amount of entities to compute. 
		! -1_ik due to header of file.
		!------------------------------------------------------------------------------
		IF (nn > covo_amnt_lines-1_ik) EXIT 

		!------------------------------------------------------------------------------
		! Call the domain an active one.
		!------------------------------------------------------------------------------
		activity(mii) = 1_ik

		!------------------------------------------------------------------------------
		! Activity = 1 (Set above during init)
		!------------------------------------------------------------------------------
		CALL MPI_SEND(Activity(mii), 1_mik, MPI_INTEGER8, mii, mii, MPI_COMM_WORLD,ierr)
		CALL print_err_stop(std_out, "MPI_SEND of activity failed.", INT(ierr, KIND=ik))

		!------------------------------------------------------------------------------
		! MPI_Send(buf, count, datatype, dest, tag, comm, ierror)
		!------------------------------------------------------------------------------
		CALL MPI_SEND(tin%dmn, 1_mik, MPI_INTEGER8, mii, mii, MPI_COMM_WORLD, ierr)
		CALL print_err_stop(std_out, "MPI_SEND of domain number failed.", INT(ierr, KIND=ik))			
         
		!------------------------------------------------------------------------------
		! Log to monitor file
		!------------------------------------------------------------------------------
		WRITE(fh_mon, FMT_MSG_2AI0)"MPI rank: ",mii, " Domain number: ", tin%dmn
		FLUSH(fh_mon)

		nn = nn + 1_ik
		
		CALL MPI_IRECV(activity(mii), 1_mik, MPI_INTEGER8, mii, mii, MPI_COMM_WORLD, &
			req_list(mii), ierr)
		CALL print_err_stop(std_out, "MPI_IRECV of activity(mii) failed.", INT(ierr, KIND=ik))

	END DO

	!------------------------------------------------------------------------------
	! Waits for a process which completed its first asssignment.
	!------------------------------------------------------------------------------
	CALL MPI_WAITANY(size_mpi-1_mik, req_list(mii), index, stmpi, ierr)
	CALL print_err_stop(std_out, "MPI_WAITANY on req_list for IRECV of activity(mii) failed.", &
		INT(ierr, KIND=ik))

	mii = index

	!------------------------------------------------------------------------------
	! As long as domains are still left 
	!------------------------------------------------------------------------------
	DO WHILE (nn <= covo_amnt_lines-1_ik)

		activity(mii) = 1_ik
		
		CALL MPI_SEND(activity(mii), 1_mik, MPI_INTEGER8 , mii, mii, MPI_COMM_WORLD, ierr)
		CALL print_err_stop(std_out, "MPI_SEND of activity failed.", INT(ierr, KIND=ik))

		CALL MPI_SEND(tin%dmn, 1_mik, MPI_INTEGER8, mii, mii, MPI_COMM_WORLD, ierr)
		CALL print_err_stop(std_out, "MPI_SEND of domain number failed.", INT(ierr, KIND=ik))			

		!------------------------------------------------------------------------------
		! Log to monitor file
		!------------------------------------------------------------------------------
		WRITE(fh_mon, FMT_MSG_2AI0)"MPI rank: ",mii, " Domain number: ", tin%dmn
		FLUSH(fh_mon)
		
		nn = nn + 1_ik
		
		CALL MPI_IRECV(activity(mii), 1_mik, MPI_INTEGER8, mii, mii, MPI_COMM_WORLD, req_list(mii), ierr)
		CALL print_err_stop(std_out, "MPI_IRECV of Activity(ii) failed.", INT(ierr, KIND=ik))

		Call MPI_WAITANY(size_mpi-1_mik, req_list, index, stmpi, ierr)
		CALL print_err_stop(std_out, "MPI_WAITANY on req_list for IRECV of Activity(ii) failed.", INT(ierr, KIND=ik))

		mii = index

	END DO

	CALL MPI_WAITALL(size_mpi-1_mik, req_list, statuses_mpi, ierr)
	CALL print_err_stop(std_out, "MPI_WAITANY on req_list for IRECV of activity(ii) failed.", &
		INT(ierr, KIND=ik))

!------------------------------------------------------------------------------
! Ranks > 0 -- Workers
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

		IF(active == -1) Exit

		!------------------------------------------------------------------------------
		! Receive domain number 
		! Receive tensor 
		! Receive tensor position in respect to the control volume
		!------------------------------------------------------------------------------
		CALL MPI_RECV(tin%dmn, 1_mik, MPI_INTEGER8, 0_mik, my_rank, MPI_COMM_WORLD, stmpi, ierr)
		CALL print_err_stop(std_out, "MPI_RECV on tin%dmn failed.", INT(ierr, KIND=ik))

		CALL MPI_RECV(tin%mat, 36_mik, MPI_DOUBLE_PRECISION, 0_mik, my_rank, MPI_COMM_WORLD, stmpi, ierr)
		CALL print_err_stop(std_out, "MPI_RECV on tin%mat failed.", INT(ierr, KIND=ik))

		CALL MPI_RECV(tin%pos, 3_mik, MPI_DOUBLE_PRECISION, 0_mik, my_rank, MPI_COMM_WORLD, stmpi, ierr)
		CALL print_err_stop(std_out, "MPI_RECV on tin%pos failed.", INT(ierr, KIND=ik))


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
		crs_counter = 1_ik

		DO jj = 1_ik, 4_ik
			IF(execute_optimization(jj)) THEN

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

					SELECT CASE(jj))
						CASE(1); CALL opt_stiff('monotropic'); temp_suf = "mono"
						CASE(2); CALL opt_stiff('orthotropic'); temp_suf = "orth"
						CASE(3); CALL opt_stiff('anisotropic1'); temp_suf = "an1"
						CASE(4); CALL opt_stiff('anisotropic2'); temp_suf = "an2"
					END SELECT			

				END DO

				!------------------------------------------------------------------------------
				! At the end of the second step, the results acutally get written to the 
				! output variables that are then send back to my_rank=0 / main process.
				!
				! Domain always gets retrieved from tin directly.
				!------------------------------------------------------------------------------
				tlcl_res%dmn(crs_counter) = tin%dmn
				tlcl_res%pos(crs_counter) = tout%pos
				tlcl_res%mat(crs_counter) = tout%mat

				CALL write_criteria_space_to_vtk(&
					out%p_n_bsnm//".tcr."//dmn_no//"."//TRIM(temp_suf)//"."//vtk_suf)

				crs_counter = crs_counter + 1_ik

			END IF
		END DO
		!------------------------------------------------------------------------------
		! tin%dmn stays the same for all tlcl_res(:)
		!------------------------------------------------------------------------------
		CALL MPI_SEND(tlcl_res%dmn, INT(crs, KIND=mik)*1_mik, MPI_INTEGER8, 0_mik, tin%dmn, &
			MPI_COMM_WORLD, ierr)
		CALL print_err_stop(std_out, "MPI_SEND on tlcl_res%dmn failed.", INT(ierr, KIND=ik))

		CALL MPI_SEND(tlcl_res%mat, INT(crs, KIND=mik)*36_mik, MPI_DOUBLE_PRECISION, 0_mik, tin%dmn, &
			MPI_COMM_WORLD, ierr)
		CALL print_err_stop(std_out, "MPI_SEND on tlcl_res%mat failed.", INT(ierr, KIND=ik))

		CALL MPI_SEND(tlcl_res%pos, INT(crs, KIND=mik)*3_mik, MPI_DOUBLE_PRECISION, 0_mik, tin%dmn, &
			MPI_COMM_WORLD, ierr)
		CALL print_err_stop(std_out, "MPI_SEND on tlcl_res%mat failed.", INT(ierr, KIND=ik))

		CALL MPI_SEND(active, 1_mik, MPI_INTEGER, 0_mik, tin%dmn, MPI_COMM_WORLD, ierr)
		CALL print_err_stop(std_out, "MPI_SEND on active failed.", INT(ierr, KIND=ik))

	END DO

END IF









!-- Read the raw data; Sort the raw data (via array index); Thresholding on raw data
DO kk=2 , lines
   READ (un(1), '(A)') line

   CALL parse(str=line, delims=",", args=tokens, nargs=ntokens)
   CALL value_di(tokens(1), vali, ios=ios)

   vali = vali+1_ik        !-- vali+1_ik; otherwise index=0 due to dmnnr=0
   t2nd(:, vali)%dn = vali-1_ik

   IF ( vali .GT. lines ) THEN
      WRITE(*,'(A)')
      WRITE(*,'(A)') "Program aborted because the line count is less than Domains are numbered."
      WRITE(*,'(A)') "Either the file or simulation is corrupted or it was cleaned with another program before."
      WRITE(*,'(A)') "This program features a complete tensor processing for a struct-process chain."
      WRITE(*,'(A)')
      GOTO 9999            !-- jump to  end program - other ways to do that?
   END IF
   oo=2_ik                 !-- oo mind. 2! 1 = Domain number
   DO mm=1,6
      DO nn=1,6
         IF (LEN_TRIM(tokens(oo)) > 0_ik ) CALL value_dr(tokens(oo), t2nd(1, vali)%Smat(nn,mm), ios=ios)
         oo=oo+1_ik
      END DO
   END DO
   t2nd(1, vali)%thres=0_ik
   DO pp=1,6                                                                ! calculate the mean value of the Trace
      t2nd(1, vali)%thres = t2nd(1, vali)%thres+t2nd(1, vali)%Smat(pp,pp)
   END DO
   t2nd(1, vali)%thres=t2nd(1, vali)%thres/6._rk
   IF ( t2nd(1, vali)%thres .LT. EY*factor_low_thres ) THEN
      t2nd(1, vali)%thres_low = .FALSE.
      t2nd(1, vali)%zrmtx=1                                                 !-- Assignment criteria 0
   END IF
END DO
write(*,*) "How many lines: ", how_many_lines
write(*,*) "         lines: ",          lines


IF (how_many_lines > lines) THEN
   how_many_lines = lines
ELSE IF (how_many_lines /= 0_ik) THEN
   lines=how_many_lines
END IF

zero_matrix_counter=0

!-- Beginn processing of each individual domain specific tensor
WRITE(*, '(A)')

!last_ii_lns_rlvnt=50
DO ii_lines = 1, lines
   IF ( t2nd(1, ii_lines)%thres_low .EQV. .TRUE. .AND. t2nd(1, ii_lines)%thres_high .EQV. .TRUE. ) THEN
      CALL EXECUTE_COMMAND_LINE('printf "\033c"',CMDSTAT=stat)      !-- Clears the command line
      !-- User information
      WRITE(*,"(A)")
      WRITE(*,"(2A)") " Input file: ", input_file
      WRITE(*,'(A)') "-----------------------------------------------------------------"
      WRITE(*,'(A)') " Program cleans, sorts and optimises the tensors."
      WRITE(*,'(A)') " Entries with less then 0.1% of the true stiffness are ignored."
      WRITE(*,'(A)') " Entries with more then 100% of the true stiffness are ignored."
      WRITE(*,*)
      WRITE(*,'(A,I5,A)')    " Assumed true (monolithic) stiffness: ",EY," MPa"
      WRITE(*,*)
      WRITE(*,"(A)") " It can take several seconds until an update occurs."
      WRITE(*,'(A)') "-----------------------------------------------------------------"
      WRITE(*,"(A)")
      WRITE(*,'(A,I5,A,I5)') " Domain No.:", INT(t2nd(1, ii_lines)%dn, ik), " of", INT(how_many_lines, ik)
      WRITE(*,"(A)")
      !------------------------------------------------------------------------------
      WRITE(*,'(A,F6.1,A,F6.1,A,F6.1,2A)')" alpha: ",am1,"  phi: ",pm1,&
           & "  eta: ",em1, "  Min ","monotropic"
      WRITE(*,'(A,F6.1,A,F6.1,A,F6.1,2A)')" alpha: ",am2,"  phi: ",pm2,&
           & "  eta: ",em2, "  Min ","monotropic"
      WRITE(*,'(A,F6.1,A,F6.1,A,F6.1,2A)')" alpha: ",ao1,"  phi: ",po1,&
           & "  eta: ",eo1, "  Min ","orthotropic"
      WRITE(*,'(A,F6.1,A,F6.1,A,F6.1,2A)')" alpha: ",ao2,"  phi: ",po2,&
           & "  eta: ",eo2, "  Min ","orthotropic"
      WRITE(*,"(A)")

      ETAt=INT((finish-start)*lines-((finish-start)*ii_lines),ik)
      percent=(REAL(ii_lines,rk)/REAL(lines-1,rk))*100.0_rk

      !-- Print current status
      WRITE(*,"(A)",Advance="NO")  " "
      CALL etimea(ETAt)
      WRITE(*,'(A,F5.1,A)',advance='YES')" (",percent,"%)"
      WRITE(*,"(A)")  "----------------------------------------------------------------"

      IF (ii_lines-zero_matrix_counter .GT. 3_ik ) THEN
         WRITE(*,"(A)") " Previous Tensor input:"
         DO zz=1,6
            WRITE(*,"(A,6F10.3,A)")"[", t2nd(1, last_ii_lns_rlvnt)%Smat(:,zz) ,"  ]"
         END DO

         DO ii=3, 4
            WRITE(*,"(A)") "----------------------------------------------------------------"
            IF (ii == 3) WRITE(*, "(A)") " Previous Tensor optimized monotropic:"
            IF (ii == 4) WRITE(*, "(A)") " Previous Tensor optimized orthotropic:"
            DO zz=1,6
               WRITE(*,"(A,6F10.3,A)")"[", t2nd(ii, last_ii_lns_rlvnt)%Smat(:,zz) ,"  ]"
            END DO
            IF (t2nd(ii, last_ii_lns_rlvnt)%optimised .EQ. 1_ik ) THEN
               WRITE(*,"(A)") " Tensor optimized:     YES"
            ELSE
               WRITE(*,"(A)") " Tensor optimized:     NO"
            END IF
         END DO

         WRITE(*, "(A)")  "----------------------------------------------------------------"
!        WRITE(*, "(3A,I6)") " Calculat orth =  ", Calculat, "       Amount of non-123: ", calcfail
      END IF

      CALL CPU_TIME(start)

      !-- Calculate optimisation
      IF ( rprt_mat .EQV. .TRUE. ) THEN
         WRITE(*, '(A)')"Monotropic  - 1st stage"
         CALL write_matrix(t2nd(1, ii_lines)%Smat,1_ik,"smat","")
      END IF

      CALL opt_eff_stiff(1_ik, t2nd(1,ii_lines)%Smat, deg_a=-91_ik, stp_a=182_ik, deg_p=-91_ik, &
           stp_p=182_ik, deg_e=-91_ik, stp_e=182_ik, intervall=1._rk, opt_tensor=mono_opt, a=a, p=p, e=e, outp=rprt)
      am1=a     ! variables stored for use in user output
      pm1=p
      em1=e

      IF ( rprt_mat .EQV. .TRUE. ) THEN
         CALL write_matrix(mono_opt,1_ik,"Monotropic optimised stage 1","")
         WRITE(*, '(A)')"Monotropic  - 2nd stage"
      END IF

      !-- tune optimization, there can be a shitload of steps requested!!
      CALL opt_eff_stiff(1_ik, t2nd(1, ii_lines)%Smat, deg_a=INT(a-1,ik), &
           stp_a=80, deg_p=INT(p-1,ik), stp_p=80, deg_e=INT(e-1,ik), stp_e=80, &
           intervall=0.025_rk, opt_tensor=mono_opt, a=a, p=p, e=e, outp=rprt)
      am2=a
      pm2=p
      em2=e

      t2nd(3, ii_lines)%euclid_ang = (/ a, p, e /)

      IF ( rprt_mat .EQV. .TRUE. ) THEN
         CALL write_matrix(mono_opt, 1_ik, "Monotrop optimised stage 2", "")
      END IF

      CALL tilt_tensor (mono_opt, mono_opt)

      CALL checksym (mono_opt, mono_sym)

      IF ( DoAM(t2nd(1, ii_lines)%Smat) .GT. DoAM(mono_opt)) THEN
         t2nd(3, ii_lines)%optimised=1_ik
      END IF

      !-- Ortho Opt
      CALL opt_eff_stiff(2_ik,t2nd(1, ii_lines)%Smat, deg_a=-91_ik, stp_a=182_ik, deg_p=-91_ik,&
           stp_p=182_ik, deg_e=-91_ik, stp_e=182_ik, intervall=1._rk, &
           opt_tensor=orth_opt, a=a, p=p, e=e, outp=rprt)
      ao1=a
      po1=p
      eo1=e

      IF ( rprt_mat .EQV. .TRUE. ) THEN
         CALL write_matrix(orth_opt, 1_ik, "Orthotropic optimised stage 1", "")
         WRITE(*, '(A)') "Orthotropic - 1st stage"
      END IF

      !-- Tilting ------------------------------------------------------------------------------
      CALL tilt_tensor (orth_opt, orth_opt)

      !-- tune optimization, there can be a shitload of steps requested!!
      CALL opt_eff_stiff(2_ik, orth_opt, deg_a=-1_ik, &
           stp_a=80_ik, deg_p=-1_ik, stp_p=80_ik, deg_e=-1_ik, stp_e=80_ik, &
           intervall=0.025_rk, opt_tensor=orth_opt, a=a, p=p, e=e, outp=rprt)
      ao2=a
      po2=p
      eo2=e

      t2nd(4, ii_lines)%euclid_ang = (/ a, p, e /)

      IF ( rprt_mat .EQV. .TRUE. ) THEN
         CALL write_matrix(orth_opt, 1_ik, "Orthotropic optimised stage 2", "")
      END IF

      !-- Tilting ------------------------------------------------------------------------------
      CALL tilt_tensor (orth_opt, orth_opt)

      CALL checksym (orth_opt, orth_sym)

   IF ( DoAO(t2nd(1, ii_lines)%Smat) .GT. DoAO(orth_opt)  ) THEN
      t2nd(4, ii_lines)%optimised=1_ik
   END IF

   !-- Degree of Anisotropy - average 0 divided by average non-zero entry     ! All values are calculated again, because there 
   t2nd(1, ii_lines)%DoA = DoAM(t2nd(1, ii_lines)%Smat)                 ! Special Case!! Ortho DoA saved in file, Mono DoA within runtime variable only!!
   t2nd(3, ii_lines)%DoA = DoAM(mono_opt)                                ! Special Case!! Ortho DoA saved in file, Mono DoA within runtime variable only!!
   Sum_DoA_InM  = Sum_DoA_InM  + t2nd(1, ii_lines)%DoA                   ! All of this stuff may be written in arrays to perfome the operations once.
   Sum_DoA_Mono = Sum_DoA_Mono + t2nd(3, ii_lines)%DoA
   Avg_DoA_InM  = Sum_DoA_InM  / REAL(ii_lines, rk)
   Avg_DoA_Mono = Sum_DoA_Mono / REAL(ii_lines, rk)
   !-- Degree of Anisotropy - average 0 divided by average non-zero entry
   t2nd(2, ii_lines)%DoA = DoAO(t2nd(1, ii_lines)%Smat)
   t2nd(4, ii_lines)%DoA = DoAO(orth_opt)
   Sum_DoA_InO  = Sum_DoA_InO  + t2nd(2, ii_lines)%DoA
   Sum_DoA_Orth = Sum_DoA_Orth + t2nd(4, ii_lines)%DoA
   Avg_DoA_InO  = Sum_DoA_InO  / REAL(ii_lines, rk)                            ! Valid, because zeromatrices are considererd empty volume - therefore part of average
   Avg_DoA_Orth = Sum_DoA_Orth / REAL(ii_lines, rk)

   CALL CPU_TIME(finish)

   CALL checksym (t2nd(1, ii_lines)%Smat, t2nd(2, ii_lines)%Smat)     !-- Assignment criteria 0
   t2nd(3, ii_lines)%Smat = mono_sym                                  !-- Assignment criteria 1 - monotropic
   t2nd(4, ii_lines)%Smat = orth_sym                                  !-- Assignment criteria 2 - orthotropic
   last_ii_lns_rlvnt=ii_lines
  ELSE
     DO ii=2, 4
        t2nd(ii, ii_lines)%Smat      = M_Null                         !-- Assignment criteria 0
        t2nd(ii, ii_lines)%DoA       = 0.0_rk                         !-- Assignment criteria 0
        t2nd(ii, ii_lines)%optimised = 0_ik                           !-- Assignment criteria 0
     END DO
     zero_matrix_counter = zero_matrix_counter+1
  END IF
END DO

! copy optimised tensors to files
! write header
DO ii =2, 4
   OPEN( UNIT=un(ii), FILE=TRIM(flnm(ii)), STATUS='NEW' )      ! Open all output files
   CALL insertstr(header, "euclid_alpha, euclid_phi, euclid_eta", LEN_TRIM(header)+1_ik)
   WRITE(un(ii), '(A)') TRIM(header)
END DO

DO ii = 1, how_many_lines
   IF ( t2nd(1, ii)%thres_low      .EQV.  .TRUE. ) THEN
      IF (  t2nd(1, ii)%thres_high .EQV.  .TRUE. ) THEN
         IF ( t2nd(3, ii)%del      .EQV. .FALSE. ) THEN
            IF ( t2nd(4, ii)%del   .EQV. .FALSE. ) THEN
               DO jj=2, 4
                  CR_str = " "
                  WRITE( string, '(F12.3)' ) t2nd(1, ii)%dn
                  CALL insertstr(CR_str, TRIM(string), LEN_TRIM(CR_str)+1_ik)
                  DO mm=1, 6
                     DO nn=1, 6             ! Acts like a header sorting algorithm
                        WRITE( string, '(F12.3)' ) t2nd(jj, ii)%Smat(nn,mm)
                        CALL insertstr(CR_str, ", "//TRIM(string), LEN_TRIM(CR_str)+1_ik)
                     END DO
                  END DO

                  ! write euclidean angles of optimized tensors to file - spatial orientation
                  CALL insertstr(CR_str, ", ", LEN_TRIM(CR_str)+1_ik) ! hardcoded workaround..... for specific header
                  DO kk=1, 3
                     WRITE( string, '(F12.3)' ) t2nd(jj, ii)%euclid_ang(kk)
                     ! write(*,*)jj
                     ! write(*,*)kk
                     ! write(*,*) t2nd(jj, ii_lines)%euclid_ang(kk)
                     CALL insertstr(CR_str, ", "//TRIM(string), LEN_TRIM(CR_str)+1_ik)
                  END DO

                  ! write string to file
                  WRITE (un(jj), '(A)') TRIM(CR_str)
               END DO
            END IF
         END IF
      END IF
   END IF
END DO

! close files
DO ii =2, 4
   CLOSE( UNIT=un(ii))
END DO

ii=lines

! Print final output
WRITE(*, '(A)')
WRITE(*, '(A)'   , advance='NO' ) " Program index."
WRITE(*, '(I5,A)', advance='YES') zero_matrix_counter," lines contained thresholded tensors."
WRITE(*, '(A)')
WRITE(*, "(A)")  "----------------------------------------------------------------"


!------------------------------------------------------------------------------
! Log a few last information
!------------------------------------------------------------------------------
IF(my_rank == 0) THEN
    CALL CPU_TIME(global_finish)

    WRITE(fh_res,FMT_TXT_AF15A) 'Overall Time = ', (global_finish - global_start) / 60,' Minutes'
    WRITE(fh_res,FMT_TXT_SEP)  
    WRITE(fh_res,FMT_TXT_AF0A) 'CPU time = ', (global_finish - global_start) / 60 / 60 * size_mpi,' Hours'
    WRITE(fh_res,FMT_MSG_SEP)f
END IF


!------------------------------------------------------------------------------
! Only used in specific cases to finish more gracefully. (grep -i "GOTO")
!------------------------------------------------------------------------------
1001 Continue

IF(my_rank == 0) THEN
    !------------------------------------------------------------------------------
    ! Finish the program
    !------------------------------------------------------------------------------
    CALL meta_signing(binary)
    CALL meta_close()

    CALL meta_stop_ascii(fh_log, log_suf)
    CALL meta_stop_ascii(fh_mon, mon_suf)

    CALL meta_stop_ascii(fh_mono, suf_mono)
    CALL meta_stop_ascii(fh_orth, suf_orth)
    CALL meta_stop_ascii(fh_ani1, suf_ani1)
    CALL meta_stop_ascii(fh_ani2, suf_ani2)

    IF(std_out/=6) CALL meta_stop_ascii(fh=std_out, suf='.std_out')

    WRITE(std_out,FMT_TXT) 'Program index successfully.'
    WRITE(std_out,FMT_MSG_SEP)

END IF ! (my_rank == 0)

CALL MPI_FINALIZE(ierr)
CALL print_err_stop(std_out, "MPI_FINALIZE failed.", INT(ierr, KIND=ik))

END PROGRAM tensor_optimizer

