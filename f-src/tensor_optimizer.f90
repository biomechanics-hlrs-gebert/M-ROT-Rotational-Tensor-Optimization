!----------------------------------------------------------------------------------------------
!--
!-- Program to clean and optimise the entries of a csv file containing 2nd order R6x6
!-- stiffness tensors. It's intended to provide plausible data from bone tissue
!--
!-- Additionally, the tensors will be optimised and sorted
!--
!-- \author:                Johannes Gebert        (HLRS, NUM)
!--
!-- \date 27.10.2020
!-- \modified 09.03.2021
!--
!-- The program is the result of a couple of investigations in optimising, sorting and
!-- rotating R6x6 2nd order stiffness tensors.
!--
!-- Known deficits within this Program:
!--
!--    There still are invalid optimised tensors.
!--
!--    Crit_mono, Crit_orth and Degree of Anisotropy address the same issue, however it's not
!--       entirely clear, which solution will be the gold standard for ongoing work.
!--    Comparing Input DoA with optimised DoA may be unnecessary.
!--    Some Variables can be combined to an array to apply recurring operations to one op.
!--    The file handling and checking is at a very early stage of development and may need
!--       further improvements.
!--    Compilation is done via a bash script. A makefile may be nice :-)
!--
!-- Remarks:
!--
!--    The Tensors are read via a csv-file. On Linux, it needs no file name extension
!--    The Tensors are read within a loop in which they are filtered first.
!--    During reading the csv, the Tensors are placed within the array according to their
!--       domain number! Not according to their position within the csv file.
!--       This ensures sorting the Tensors without much computational effort.
!--    Optimising by two stages to 1° / 0.025° provides excellent results.
!--       However, for monotropic optimisations, there are some pretty bad optimised tensors...
!--    After optimising the tensors, it's necessary to tilt them to S11>S22>S33
!--       This will give a plausible sorting AND a measure for stiffness trajectories
!--
!----------------------------------------------------------------------------------------------
!-- Entries with less then 0.1% of the true stiffness are designated as remainder
!-- Entries with more then 100% of the true stiffness are designated as remainder
!-- The tensors are rotated during the optimisation to get the following criteria:
!-- S11 > S22 > S33
!-- It works with a hardcoded, parametrized Young-Modulus (E/EY)
!--
!-- Keeps the original file with has to be given as command argument
!-- CR0 - sym checked tensors
!-- CR1 - monotropic  optimised tensors
!-- CR2 - orthotropic optimised tensors
!-- ignores zero-tensors/domains
!----------------------------------------------------------------------------------------------
PROGRAM tensor_optimizer

USE, INTRINSIC :: iso_fortran_env, ONLY : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit
USE standards
USE strings
USE opt_stiffness
! USE operating_System


IMPLICIT NONE

  !----------------------------------------------------------------------------------------------
  !-- Ignore zero matrices by commenting in the according lines.
  !-- Handling der Nullmatrizen über diverse Kommentierungen im file.
  !--
  !-- Variablen zum Programm verwalten und Rechenzeit sparen
  INTEGER     (KIND=ik)                                      :: zero_matrix_counter, stat
  REAL        (KIND=rk)                                      :: percent

  !-- Allgemeine Laufvariablen um numerischen Kernel zu steuern
  INTEGER     (KIND=ik)                                      :: ii_lines, last_ii_lns_rlvnt, &
       ii, jj, kk, ll, mm, nn, oo, pp, vali, zz
  REAL        (KIND=rk), DIMENSION(6,6)                      :: mono_opt, orth_opt, mono_sym, orth_sym, M_Null
  INTEGER     (KIND=ik), PARAMETER                           :: header_row_header=1
  CHARACTER   (LEN=200)                                      :: input_file, header
  CHARACTER   (LEN=1000)                                     :: line, CR_str
  CHARACTER   (LEN=12)                                       :: string
  INTEGER     (KIND=ik), PARAMETER                           :: n_cols=38, header_row=1
  INTEGER     (KIND=4)                                       :: lines, lunit=10, how_many_lines=10
  CHARACTER   (LEN=10)                                       :: hw_mny_lns
  INTEGER     (KIND=ik), DIMENSION(4)                        :: un
  INTEGER  (KIND=ik)                                         :: ios, ntokens
  CHARACTER(len=mcl)                                         :: tokens(100)

  TYPE :: csv_data
      REAL    (KIND=rk)                                      :: dn, thres, Doa        ! REAL damit write bei csv geht!!
      REAL    (KIND=rk), DIMENSION(6,6)                      :: Smat
      REAL    (KIND=rk), DIMENSION(3)                        :: euclid_ang=(/ 0._rk, 0._rk, 0._rk /)
      INTEGER (KIND=ik)                                      :: zrmtx=0, optimised=0
      LOGICAL (KIND=rk)                                      :: thres_high=.TRUE., thres_low=.TRUE., del=.FALSE.   ! initialisation
  END TYPE csv_data

  TYPE(csv_data)       , DIMENSION(:,:) , ALLOCATABLE        :: t2nd                  ! Tensor of 2nd order R6x6
  CHARACTER   (LEN=200), DIMENSION(4)                        :: flnm
  CHARACTER   (LEN=*)  , PARAMETER                           :: int_fmt='(I5)'
  CHARACTER   (LEN=*)  , PARAMETER                           :: REAL_fmt='(F9.3)'
  CHARACTER   (LEN=*)  , PARAMETER                           :: REAL_int='(F9.0)'
  REAL        (KIND=rk)                                      :: start, finish
  INTEGER     (KIND=ik)                                      :: ETAt

  !-- Optimization variables
  REAL        (KIND=rk)                                      :: a,p,e
  REAL        (KIND=rk)                                      :: Sum_DoA_Mono=0, Sum_DoA_Orth=0, Sum_DoA_InO=0, Sum_DoA_InM=0
  REAL        (KIND=rk)                                      :: Avg_DoA_Mono, Avg_DoA_Orth, Avg_DoA_InO, Avg_DoA_InM
  REAL        (KIND=rk)                                      :: am1, pm1, em1, am2, pm2, em2, ao1, po1, eo1, ao2, po2, eo2

  !-- Parameters
  INTEGER     (KIND=ik), PARAMETER                           :: EY=5600
  REAL        (KIND=rk), PARAMETER                           :: factor_low_thres=0.001
  CHARACTER   (LEN=3)                                        :: calculat
  INTEGER     (KIND=ik)                                      :: calcfail=0
  LOGICAL                                                    :: rprt_mat, rprt

!-- DEFAULTS ----- rprt=.TRUE. ----- rprt_mat=.FALSE. -----------------------------------------
rprt=.TRUE.
rprt_mat=.FALSE.
!-- DEFAULTS ----------------------------------------------------------------------------------

CALL GET_COMMAND_ARGUMENT(1, input_file)
CALL GET_COMMAND_ARGUMENT(2, hw_mny_lns)
READ(hw_mny_lns,"(I10)") how_many_lines

!-- Abstand in cli
WRITE(*, '(A)')

!-- read lines
lines = 0
OPEN (lunit, file = input_file)
DO
    READ (lunit,*, END=10)
    lines = lines + 1_4
END DO
10 CLOSE (lunit)

!-- add flnm-extension
flnm(1) = TRIM(input_file)//".meta"
flnm(2) = TRIM(input_file(1:(LEN_TRIM(input_file)-4)))//"_CR0.csv" ! _c_ -> cleaned
flnm(3) = TRIM(input_file(1:(LEN_TRIM(input_file)-4)))//"_CR1.csv" ! MANDEL!
flnm(4) = TRIM(input_file(1:(LEN_TRIM(input_file)-4)))//"_CR2.csv"

ALLOCATE(t2nd(4, lines))

M_Null(:,:)=.0_rk

un(1) = 29_ik            ! input
un(2) = 30_ik            ! out CR0
un(3) = 31_ik            ! out CR1
un(4) = 32_ik            ! out CR2

!-- Initialize and open  the csv-files to write optimized data
OPEN( UNIT=un(1), FILE=TRIM(input_file ), STATUS='OLD' )
READ (un(1), '(A)') header

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
      !----------------------------------------------------------------------------------------------
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
WRITE(*, '(A)'   , advance='NO' ) " Program finished."
WRITE(*, '(I5,A)', advance='YES') zero_matrix_counter," lines contained thresholded tensors."
WRITE(*, '(A)')
WRITE(*, "(A)")  "----------------------------------------------------------------"

9999 CONTINUE

END PROGRAM tensor_optimizer

