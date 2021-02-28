!----------------------------------------------------------------------------------------------
!--
!-- Program to clean and optimise the entries of a csv file containing 2nd order R6x6
!-- stiffness tensors. It's intended to provide plausible data from bone tissue
!--
!-- Additionally, the tensors will be optimised and sorted
!--
!-- \author:                Johannes Gebert        (HLRS, NUM)
!--
!-- \date 27.10.2020        opt_eff_stiffness:::::Evo:3_x86_64
!--
!-- The program is the result of a couple of investigations in optimising, sorting and
!-- rotating R6x6 2nd order stiffness tensors.
!--
!-- Known deficits within this Program:
!--
!--    There still are invalid optimised tensors. 
!--
!--    The csv-header needs to be provided as an extra file within the directory
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
--
!----------------------------------------------------------------------------------------------
!-- Requires a header as follows:
!-- Domain No,S11,S12,S13,S14,S15,S16,S21,S22,S23,S24,S25,S26,S31,S32,S33,S34,S35,S36, ...
!-- ... S41,S42,S43,S44,S45,S46,S51,S52,S53,S54,S55,S56,S61,S62,S63,S64,S65,S66,
!--
!-- Entries with less then 0.1% of the true stiffness are designated as remainder
!-- Entries with more then 100% of the true stiffness are designated as remainder
!-- The tensors are rotated during the optimisation to get the following criteria:
!-- S11 > S22 > S33
!-- It works with a hardcoded Young-Modulus (E/EY)
!--
!-- Keeps the original file with has to be given as command argument
!-- CR0 - sym checked tensors
!-- CR1 - monotropically  optimised tensors
!-- CR2 - orthotropically optimised tensors
!-- ignores zero-tensors/domains
!----------------------------------------------------------------------------------------------
program opt_and_sort_eff_stiffness

use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit
#ifdef f2003
#else
#define stdin  5
#define stdout 6
#define stderr 0
#endif

use csv_module
use vtkio
use mat_matrices
use opt_stiffness
Use Operating_System


implicit none

  !----------------------------------------------------------------------------------------------
  !-- Ignore zero matrices by commenting in the according lines.
  !-- Handling der Nullmatrizen über diverse Kommentierungen im file.
  !--
  !-- Variablen zum Programm verwalten und Rechenzeit sparen
  integer(kind=ik)                                   :: zero_matrix_counter, stat
  real(kind=rk)                                      :: percent

  !-- Allgemeine Laufvariablen um numerischen Kernel zu steuern
  integer(kind=ik)                                   :: kk
  integer(kind=ik)                                   :: ii_lines, last_ii_lns_rlvnt, ii, ll, mm, nn, oo, pp, vali, zz
  real(kind=rk), dimension(6,6)                      :: mono_opt, orth_opt, mono_sym, orth_sym, M_Null

  !-- Initialisierung zum öffnen der Datei zum einlesen des Headers
  integer(kind=ik),parameter                         :: header_row_header=1

  !-- Objekte und Variablen zur internen Verwaltung der Daten
  type(csv_file)                                     :: input, input_header, output_CR0, output_CR1, output_CR2
  character(len=1), parameter                        :: quote=' ', delimiter=','
  logical(kind=rk)                                   :: status_ok
  logical(kind=rk), parameter                        :: verbose=.FALSE.
  real(kind=rk)                                      :: val                                  !! the value to add
  real(kind=rk), dimension(36)                       :: buffer_vector


  !-- Variablen zur Eingabe der Rohdaten
  character(len=200)                                 :: input_file
  type(csv_string),dimension(:),allocatable          :: header
  integer(kind=ik), parameter                        :: n_cols=38, header_row=1
  integer(kind=4)                                    :: lines, lunit=10, how_many_lines=10
  character(len=10)                                  :: hw_mny_lns

  !-- Variablen zur Ausgabe der Daten
  TYPE :: csv_data 
      Real(Kind=rk)                                  :: domain_no, threshold, Doa        ! Real damit write bei csv geht!!
      Real(Kind=rk), Dimension (6,6)                 :: Smat
      integer(kind=ik)                               :: zrmtx=0, optimised=0
      logical(kind=rk)                               :: threshold_high=.TRUE., threshold_low=.TRUE., del=.FALSE.   ! initialisation
  END TYPE csv_data

  type(csv_data), Dimension (:), Allocatable         :: Smat_dnr, Smat_dnr_CR0, Smat_dnr_CR1, Smat_dnr_CR2
  character(len=200)                                 :: op_CR0_flnm, op_CR1_flnm, op_CR2_flnm, meta,ip_header_flnm
  character(len=*), parameter                        :: int_fmt='(I5)'
  character(len=*), parameter                        :: real_fmt='(F9.3)'
  character(len=*), parameter                        :: real_int='(F9.0)'

  !-- Variablen für Berechnung der Restzeit - wär eigentlich was für type.
  real(kind=rk)                                      :: start, finish
  integer(kind=ik)                                   :: ETAt

  !-- Filehandling
  logical                                            :: me_L_ex, CR0_L_ex, CR1_L_ex, CR2_L_ex, in_L_ex, rprt, rprt_mat
  integer(kind=ik)                                   :: in_f_ex=1, CR0_f_ex, CR1_f_ex, CR2_f_ex, me_f_ex

  !-- hardcoded parameters
  integer(kind=ik), parameter                        :: EY=5600
  real(kind=rk), parameter                           :: factor_low_threshold=0.001

  !-- Optimization variables
  real(kind=rk)                                      :: a,p,e
  real(kind=rk)                                      :: DoAMonoOverall, DoAOrthOverall, DoAOverallInO, DoAOverallInM
  real(kind=rk)                                      :: Sum_DoA_Mono=0, Sum_DoA_Orth=0, Sum_DoA_InO=0, Sum_DoA_InM=0
  real(kind=rk)                                      :: Avg_DoA_Mono, Avg_DoA_Orth, Avg_DoA_InO, Avg_DoA_InM
  real(kind=rk)                                      :: am1, pm1, em1, am2, pm2, em2, ao1, po1, eo1, ao2, po2, eo2

  !-- additional transformational variables... double spending because it's implemented within mod_opt_stiffness.f90
  character(len=3)                                   :: calculat
  integer(kind=ik)                                   :: calcfail=0

!----------------------------------------------------------------------------------------------
!-- DEFAULTS ----- rprt=.TRUE. ----- rprt_mat=.FALSE. -----------------------------------------
rprt=.TRUE.
rprt_mat=.FALSE.
!-- DEFAULTS ----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------


call get_command_argument(1,input_file)
call get_command_argument(2,hw_mny_lns)
read(hw_mny_lns,"(I10)") how_many_lines 

!-- Abstand in cli
write(*,*) 

!-- read lines
lines = 0
open (lunit, file = input_file)
do
    read (lunit,*, END=10)
    lines = lines + 1_4
end do
10 close (lunit)

!-- add flnm-extension
op_CR0_flnm = trim(input_file(1:(len_trim(input_file)-4)))//"_MANDEL_CR0.csv" ! _c_ -> cleaned
op_CR1_flnm = trim(input_file(1:(len_trim(input_file)-4)))//"_MANDEL_CR1.csv"
op_CR2_flnm = trim(input_file(1:(len_trim(input_file)-4)))//"_MANDEL_CR2.csv"
meta        = trim(input_file)//".meta"
!-- Can be avoided with some precautions
ip_header_flnm="paraview-header.csv"

allocate(Smat_dnr(lines))
allocate(Smat_dnr_CR0(lines))
allocate(Smat_dnr_CR1(lines))
allocate(Smat_dnr_CR2(lines))

M_Null(:,:)=.0_rk

!-- Check whether file exist. Needs some rework.....
in_F_ex=1
CR0_F_ex=1
CR1_F_ex=1
CR2_F_ex=1
me_F_ex=1

inquire(file=input_file,  EXIST =  in_L_ex)
inquire(file=op_CR0_flnm, EXIST = CR0_L_ex)
inquire(file=op_CR1_flnm, EXIST = CR1_L_ex)
inquire(file=op_CR2_flnm, EXIST = CR2_L_ex)
inquire(file=meta,        EXIST =  me_L_ex)

if ( in_L_ex .eqv. .TRUE. ) then
in_f_ex=0
if ( CR0_L_ex .eqv. .FALSE. ) then
CR0_f_ex=0
if ( CR1_L_ex .eqv. .FALSE. ) then
CR1_f_ex=0
if ( CR2_L_ex .eqv. .FALSE. ) then
CR2_f_ex=0
if ( me_L_ex .eqv. .FALSE.  ) then
me_f_ex=0

!----------------------------------------------------------------------------------------------
!-- Initialize and open  the csv-files to write optimized data
call input%initalize(quote,delimiter,.FALSE.,.FALSE.,'T','F',100,verbose)
call output_CR0%initalize(quote,delimiter,.FALSE.,.FALSE.,'T','F',100,verbose)
call output_CR1%initalize(quote,delimiter,.FALSE.,.FALSE.,'T','F',100,verbose)
call output_CR2%initalize(quote,delimiter,.FALSE.,.FALSE.,'T','F',100,verbose)

write(*,*)"csv initialization done"
call input%read(input_file,header_row,status_ok=status_ok)
write(*,*)"csv input read"

call output_CR0%open(op_CR0_flnm,n_cols+2_ik,status_ok)
call output_CR1%open(op_CR1_flnm,n_cols+2_ik,status_ok)
call output_CR2%open(op_CR2_flnm,n_cols+2_ik,status_ok)


!----------------------------------------------------------------------------------------------
call input_header%initalize(quote,delimiter,.FALSE.,.FALSE.,'T','F',100,.FALSE.)
call input_header%read(ip_header_flnm,1,status_ok=status_ok)
call input_header%get_header(header,status_ok)

call output_CR0%add(header);    call output_CR0%next_row()
call output_CR1%add(header);    call output_CR1%next_row()
call output_CR2%add(header);    call output_CR2%next_row()
!----------------------------------------------------------------------------------------------
!-- Read the raw data
!-- Sort the raw data (via array index)
!-- Threshholding on raw data
!--
write(*,*)"csv header added"
! call sleep(3)
do kk=1 , lines
    call input%csv_get_value(kk,1,val,status_ok)
    vali=int(val,ik)+1
    if ( vali .gt. lines )then
       write(*,*) 
       write(*,'(A)')"Program aborted because the line count is less than Domains are numbered."
       write(*,'(A)')"Either the file or simulation is corrupted or it was cleaned with another program before."
       write(*,'(A)')"This program features a complete tensor processing for a struct-process chain."
       write(*,*) 
       goto 9999            !-- jump to  end program - other ways to do that?
    end if
    Smat_dnr(vali)%domain_no=real(val,rk)
    do ll=2,37
        call input%csv_get_value(kk,ll,val,status_ok)
        buffer_vector(ll-1)=real(val,rk)
        if ( buffer_vector(ll-1) .gt. EY ) then
           Smat_dnr(vali)%threshold_high=.FALSE.
        end if
    end do
    oo=1
    mm=1
    do mm=1,6
        nn=1
        do nn=1,6
            Smat_dnr(vali)%Smat(nn,mm)=buffer_vector(oo)
            oo=oo+1
        end do
     end do
     Smat_dnr(vali)%threshold=0
     do pp=1,6                                                     ! calculate the mean value of the Trace
        Smat_dnr(vali)%threshold = Smat_dnr(vali)%threshold+Smat_dnr(vali)%Smat(pp,pp)
     end do
     Smat_dnr(vali)%threshold=Smat_dnr(vali)%threshold/6._rk
     if ( Smat_dnr(vali)%threshold .lt. EY*factor_low_threshold ) then
        Smat_dnr(vali)%threshold_low=.FALSE.
        Smat_dnr(vali)%zrmtx=1                                              !-- Assignment criteria 0
     end if
end do

if (how_many_lines > lines) then
   how_many_lines = lines
else
   lines=how_many_lines
end if

zero_matrix_counter=0

!-- Beginn processing of each individual domain specific tensor
write(*,*)""

!last_ii_lns_rlvnt=50
Do ii_lines = 1, how_many_lines
   if ( Smat_dnr(ii_lines)%threshold_low .eqv. .TRUE. .and. Smat_dnr(ii_lines)%threshold_high .eqv. .TRUE. ) then
     Call execute_command_line('printf "\033c"',CMDSTAT=stat)
      !-- User information
      write(*,"(A)")
      write(*,"(2A)") " Input file: ",input_file
      write(*,'(A)') "-----------------------------------------------------------------"
      write(*,'(A)') " Program cleans, sorts and optimises the tensors."
      write(*,'(A)') " Entries with less then 0.1% of the true stiffness are ignored."
      write(*,'(A)') " Entries with more then 100% of the true stiffness are ignored."
      write(*,*)
      write(*,'(A,I5,A)')    " Assumed true (monolithic) stiffness: ",EY," MPa"
      write(*,*)
      write(*,"(A)") " It can take several seconds until an update occurs."
      write(*,'(A)') "-----------------------------------------------------------------"
      write(*,"(A)")
      write(*,'(A,I5,A,I5)') " Domain No.:", int(Smat_dnr(ii_lines)%domain_no,ik), " of", int(how_many_lines,ik)
      write(*,"(A)")
      !----------------------------------------------------------------------------------------------
      write(*,'(A,F6.1,A,F6.1,A,F6.1,2A)')" alpha: ",am1,"  phi: ",pm1,&
           & "  eta: ",em1, "  Min ","monotropic"
      write(*,'(A,F6.1,A,F6.1,A,F6.1,2A)')" alpha: ",am2,"  phi: ",pm2,&
           & "  eta: ",em2, "  Min ","monotropic"
      write(*,'(A,F6.1,A,F6.1,A,F6.1,2A)')" alpha: ",ao1,"  phi: ",po1,&
           & "  eta: ",eo1, "  Min ","orthotropic"
      write(*,'(A,F6.1,A,F6.1,A,F6.1,2A)')" alpha: ",ao2,"  phi: ",po2,&
           & "  eta: ",eo2, "  Min ","orthotropic"
     write(*,"(A)")
      ETAt=int((finish-start)*lines-((finish-start)*ii_lines),ik)
      percent=(real(ii_lines,rk)/real(lines-1,rk))*100.0_rk
      !-- Print current status
      write(*,"(A)",Advance="NO")  " "
      call etimea(ETAt)
      write(*,'(A,F5.1,A)',advance='YES')" (",percent,"%)"
      write(*,"(A)")  "----------------------------------------------------------------"
      if (ii_lines-zero_matrix_counter .gt. 3_ik ) then
         write(*,"(A)") " Previous Tensor input:"
         do zz=1,6
            write(*,"(A,6F10.3,A)")"[", Smat_dnr(last_ii_lns_rlvnt)%Smat(:,zz) ,"  ]"
         end do
         write(*,"(A)") "----------------------------------------------------------------"
         write(*,"(A)") " Previous Tensor optimised monotropic:"
         do zz=1,6
            write(*,"(A,6F10.3,A)")"[", Smat_dnr_CR1(last_ii_lns_rlvnt)%Smat(:,zz) ,"  ]"
         end do
         if (Smat_dnr_CR1(last_ii_lns_rlvnt)%optimised .eq. 1_ik ) then
            write(*,"(A)") " Tensor optimised:     YES"
         else
            write(*,"(A)") " Tensor optimised:     NO"
         end if
         write(*,"(A)") "----------------------------------------------------------------"
         write(*,"(A)") " Previous Tensor optimised orthotropic:"
         do zz=1,6
            write(*,"(A,6F10.3,A)")"[", Smat_dnr_CR2(last_ii_lns_rlvnt)%Smat(:,zz) ,"  ]"
         end do
         if (Smat_dnr_CR2(last_ii_lns_rlvnt)%optimised .eq. 1_ik ) then
            write(*,"(A)") " Tensor optimised:     YES"
         else
            write(*,"(A)") " Tensor optimised:     NO"
         end if
         write(*,"(A)") "----------------------------------------------------------------"
         write(*,"(A,F6.2,A)") " Previous DoA - input     monotropic:  ",abs(Smat_dnr    (last_ii_lns_rlvnt)%DoA),"%"
         write(*,"(A,F6.2,A)") " Previous DoA - optimised monotropic:  ",abs(Smat_dnr_CR1(last_ii_lns_rlvnt)%DoA),"%"
         write(*,"(A)")
         write(*,"(A,F6.2,A)") " Average  DoA - input     monotropic:  ",abs(Avg_DoA_InM),"%"
         write(*,"(A,F6.2,A)") " Average  DoA - optimised monotropic:  ",abs(Avg_DoA_Mono),"%"
         write(*,"(A)")
         write(*,"(A,F6.2,A)") " Maximum  DoA - input     monotropic:  ",abs(DoAOverallInM),"%"
         write(*,"(A,F6.2,A)") " Maximum  DoA - optimised monotropic:  ",abs(DoAMonoOverall),"%"
         write(*,"(A)") "------------------------------------"
         write(*,"(A,F6.2,A)") " Previous DoA - input     orthotropic: ",abs(Smat_dnr_CR0(last_ii_lns_rlvnt)%DoA),"%"
         write(*,"(A,F6.2,A)") " Previous DoA - optimised orthotropic: ",abs(Smat_dnr_CR2(last_ii_lns_rlvnt)%DoA),"%"
         write(*,"(A)")
         write(*,"(A,F6.2,A)") " Average  DoA - input     orthotropic: ",abs(Avg_DoA_InO),"%"
         write(*,"(A,F6.2,A)") " Average  DoA - optimised orthotropic: ",abs(Avg_DoA_Orth),"%"
         write(*,"(A)")
         write(*,"(A,F6.2,A)") " Maximum  DoA - input     orthotropic: ",abs(DoAOverallInO),"%"
         write(*,"(A,F6.2,A)") " Maximum  DoA - optimised orthotropic: ",abs(DoAOrthOverall),"%"
         write(*,"(A)")  "----------------------------------------------------------------"
         write(*,"(3A,I6)") " Calculat orth =  ", Calculat, "       Amount of non-123: ", calcfail
      end if
      call cpu_time(start)
       !----------------------------------------------------------------------------------------------
       !-- DEBUG -------------------------------------------------------------------------------------
       ! Smat_dnr(ii_lines) = Matrix_Iso(5600._rk,0.3_rk)
       ! Call inverse(EE, 6_ik,6)
       !-- DEBUG -------------------------------------------------------------------------------------
       !----------------------------------------------------------------------------------------------

        !-- Calculate optimisation
        if ( rprt_mat .eqv. .TRUE. ) then
            write(*,'(A)')"Monotropic  - 1st stage"
            call Write_real_matrix(6,Smat_dnr(ii_lines)%Smat,6_ik,6_ik,"smat","")
        endif

        call opt_eff_stiff(1_ik,Smat_dnr(ii_lines)%Smat,deg_a=-91_ik,stp_a=182_ik,deg_p=-91_ik,&
            &stp_p=182_ik,deg_e=-91_ik,stp_e=182_ik,intervall=1._rk,opt_tensor=mono_opt, a=a,p=p,e=e,outp=rprt)
        am1=a     ! variables stored for use in user output
        pm1=p
        em1=e
        if ( rprt_mat .eqv. .TRUE. ) then
            call Write_real_matrix(6,mono_opt,6_ik,6_ik,"Monotropic optimised stage 1","")
            write(*,'(A)')"Monotropic  - 2nd stage"
        endif

        !-- tune optimization, there can be a shitload of steps requested!!
        call opt_eff_stiff(1_ik ,Smat_dnr(ii_lines)%Smat,deg_a=int(a-1,ik), &
           & stp_a=80,deg_p=int(p-1,ik),stp_p=80,deg_e=int(e-1,ik),stp_e=80,&
           & intervall=0.025_rk,opt_tensor=mono_opt, a=a,p=p,e=e,outp=rprt)
        am2=a
        pm2=p
        em2=e
        if ( rprt_mat .eqv. .TRUE. ) then
          call Write_real_matrix(6,mono_opt,6_ik,6_ik,"Monotrop optimised stage 2","")
        endif

        !-- Tilting ------------------------------------------------------------------------------
        call tilt_tensor (mono_opt,mono_opt)

        !-- ------------------------------------------------------------------ RS ordering

        call checksym (mono_opt,mono_sym)

        if ( DoAM(Smat_dnr(ii_lines)%Smat) .gt. DoAM(mono_opt)) then
           Smat_dnr_CR1(ii_lines)%optimised=1_ik
        end if


        !-- Ortho Opt
        !-- Ortho Opt

        call opt_eff_stiff(2_ik,Smat_dnr(ii_lines)%Smat,deg_a=-91_ik,stp_a=182_ik,deg_p=-91_ik,&
             &stp_p=182_ik,deg_e=-91_ik,stp_e=182_ik,intervall=1._rk,opt_tensor=orth_opt, a=a,p=p,e=e,outp=rprt)
        ao1=a
        po1=p
        eo1=e
        if ( rprt_mat .eqv. .TRUE. ) then
           call Write_real_matrix(6,orth_opt,6_ik,6_ik,"Orthotropic optimised stage 1","")
           write(*,'(A)')"Orthotropic - 1st stage" 
        endif

        !-- Tilting ------------------------------------------------------------------------------
        call tilt_tensor (orth_opt,orth_opt)

        !-- tune optimization, there can be a shitload of steps requested!!
        call opt_eff_stiff(2_ik ,orth_opt,deg_a=-1_ik, &
             & stp_a=80_ik,deg_p=-1_ik,stp_p=80_ik,deg_e=-1_ik,stp_e=80_ik,&
             & intervall=0.025_rk,opt_tensor=orth_opt, a=a,p=p,e=e,outp=rprt)
        ao2=a
        po2=p
        eo2=e
        if ( rprt_mat .eqv. .TRUE. ) then
           call Write_real_matrix(6,orth_opt,6_ik,6_ik,"Orthotropic optimised stage 2","")
        endif

        !-- Tilting ------------------------------------------------------------------------------
        call tilt_tensor (orth_opt,orth_opt)

!---------- DEBUG --------------------------------------------------------------------------------
        If ( (orth_opt(1,1) > orth_opt(2,2)) .AND.  &
             (orth_opt(1,1) > orth_opt(3,3)) .AND.  (orth_opt(2,2) > orth_opt(3,3)) ) then
           calculat="123"
        Else If ( (orth_opt(1,1) > orth_opt(2,2)) .AND.  &
             (orth_opt(1,1) > orth_opt(3,3)) .AND.  (orth_opt(2,2) < orth_opt(3,3)) ) then
           calculat="132"
           calcfail=calcfail+1
        Else If ( (orth_opt(1,1) < orth_opt(2,2)) .AND.  &
             (orth_opt(1,1) > orth_opt(3,3)) .AND.  (orth_opt(2,2) > orth_opt(3,3)) ) then
           calculat="213"
           calcfail=calcfail+1
        Else If ( (orth_opt(1,1) < orth_opt(2,2)) .AND.  &
             (orth_opt(1,1) < orth_opt(3,3)) .AND.  (orth_opt(2,2) > orth_opt(3,3)) ) then
           calculat="231"
           calcfail=calcfail+1
        Else If ( (orth_opt(1,1) > orth_opt(2,2)) .AND.  &
             (orth_opt(1,1) < orth_opt(3,3)) .AND.  (orth_opt(2,2) < orth_opt(3,3)) ) then
           calculat="312"
           calcfail=calcfail+1
        Else If ( (orth_opt(1,1) < orth_opt(2,2)) .AND.  &
             (orth_opt(1,1) < orth_opt(3,3)) .AND.  (orth_opt(2,2) < orth_opt(3,3)) ) then
           calculat="321"
           calcfail=calcfail+1
        End If
!---------- DEBUG --------------------------------------------------------------------------------


        call checksym (orth_opt,orth_sym)

        if ( DoAO(Smat_dnr(ii_lines)%Smat) .gt. DoAO(orth_opt)  ) then
           Smat_dnr_CR2(ii_lines)%optimised=1_ik
        end if

        !-- Degree of Anisotropy - average 0 divided by average non-zero entry     ! All values are calculated again, because there 
        Smat_dnr    (ii_lines)%DoA = DoAM(Smat_dnr(ii_lines)%Smat)                 ! Special Case!! Ortho DoA saved in file, Mono DoA within runtime variable only!!
        Smat_dnr_CR1(ii_lines)%DoA = DoAM(mono_opt)                                ! Special Case!! Ortho DoA saved in file, Mono DoA within runtime variable only!!
        Sum_DoA_InM  = Sum_DoA_InM  + Smat_dnr    (ii_lines)%DoA                   ! All of this stuff may be written in arrays to perfome the operations once.
        Sum_DoA_Mono = Sum_DoA_Mono + Smat_dnr_CR1(ii_lines)%DoA
        Avg_DoA_InM  = Sum_DoA_InM  / real(ii_lines,rk)
        Avg_DoA_Mono = Sum_DoA_Mono / real(ii_lines,rk)
        !-- Degree of Anisotropy - average 0 divided by average non-zero entry
        Smat_dnr_CR0(ii_lines)%DoA = DoAO(Smat_dnr(ii_lines)%Smat)
        Smat_dnr_CR2(ii_lines)%DoA = DoAO(orth_opt)
        Sum_DoA_InO  = Sum_DoA_InO  + Smat_dnr_CR0(ii_lines)%DoA
        Sum_DoA_Orth = Sum_DoA_Orth + Smat_dnr_CR2(ii_lines)%DoA
        Avg_DoA_InO  = Sum_DoA_InO  / real(ii_lines,rk)                            ! Valid, because zeromatrices are considererd empty volume - therefore part of average
        Avg_DoA_Orth = Sum_DoA_Orth / real(ii_lines,rk)

        if ( Smat_dnr    (ii_lines)%DoA .gt. DoAOverallInM ) then
           DoAOverallInM = Smat_dnr    (ii_lines)%DoA
        end if
        if ( Smat_dnr_CR1(ii_lines)%DoA .gt. DoAMonoOverall ) then
           DoAMonoOverall= Smat_dnr_CR1(ii_lines)%DoA
        end if
        if ( Smat_dnr_CR0(ii_lines)%DoA .gt. DoAOverallInO ) then
           DoAOverallInO = Smat_dnr_CR0(ii_lines)%DoA
        end if
        if ( Smat_dnr_CR2(ii_lines)%DoA .gt. DoAOrthOverall ) then
           DoAOrthOverall= Smat_dnr_CR2(ii_lines)%DoA
        end if

        !----------------------------------------------------------------------------------------------
        call cpu_time(finish)
        !----------------------------------------------------------------------------------------------

        call checksym (Smat_dnr(ii_lines)%Smat,Smat_dnr_CR0(ii_lines)%Smat)     !-- Assignment criteria 0
        Smat_dnr_CR1(ii_lines)%Smat=mono_sym                                    !-- Assignment criteria 1 - monotropic
        Smat_dnr_CR2(ii_lines)%Smat=orth_sym                                    !-- Assignment criteria 2 - orthotropic
        last_ii_lns_rlvnt=ii_lines
    else
        Smat_dnr_CR0(ii_lines)%Smat=M_Null                                      !-- Assignment criteria 0
        Smat_dnr_CR1(ii_lines)%Smat=M_Null                                      !-- Assignment criteria 1 - monotropic
        Smat_dnr_CR2(ii_lines)%Smat=M_Null                                      !-- Assignment criteria 2 - orthotropic
        Smat_dnr_CR0(ii_lines)%DoA=0.0_rk                                       !-- Assignment criteria 0
        Smat_dnr_CR1(ii_lines)%DoA=0.0_rk                                       !-- Assignment criteria 1 - monotropic
        Smat_dnr_CR2(ii_lines)%DoA=0.0_rk                                       !-- Assignment criteria 2 - orthotropic
        Smat_dnr_CR0(ii_lines)%optimised=0_ik                                   !-- Assignment criteria 0
        Smat_dnr_CR1(ii_lines)%optimised=0_ik                                   !-- Assignment criteria 1 - monotropic
        Smat_dnr_CR2(ii_lines)%optimised=0_ik                                   !-- Assignment criteria 2 - orthotropic
        zero_matrix_counter = zero_matrix_counter+1
    end if
 end Do

! ----------------------------------------------------------------------------------------------
! -- copy optimised tensors to files

 Do ii = 1, how_many_lines
    if ( Smat_dnr(ii)%threshold_low .eqv. .FALSE. .and. Smat_dnr(ii)%threshold_high .eqv. .FALSE. ) then
        if ( Smat_dnr_CR1(ii)%del .eqv. .FALSE.) then
            call output_CR1%add(Smat_dnr(ii)%domain_no,real_fmt=real_int)
            do mm=1, 6
                do nn=1, 6             ! Acts like a header sorting algorithm
                    call output_CR1%add(Smat_dnr_CR1(ii)%Smat(nn,mm),real_fmt=real_fmt)
                end do
            end do
            call output_CR1%add(Smat_dnr_CR1(ii)%DoA,real_fmt=real_fmt)
            call output_CR1%add(Smat_dnr_CR1(ii)%optimised,real_fmt=int_fmt)
            call output_CR1%next_row()
        end if

        if ( Smat_dnr_CR2(ii)%del .eqv. .FALSE.) then
            call output_CR2%add(Smat_dnr(ii)%domain_no,real_fmt=real_int)
            do mm=1, 6
                do nn=1, 6             ! Acts like a header sorting algorithm!
                    call output_CR2%add(Smat_dnr_CR2(ii)%Smat(nn,mm),real_fmt=real_fmt)
                end do
            end do
            call output_CR2%add(Smat_dnr_CR2(ii)%DoA,real_fmt=real_fmt)
            call output_CR2%add(Smat_dnr_CR2(ii)%optimised,real_fmt=int_fmt)
            call output_CR2%next_row()
        end if

        call output_CR0%add(Smat_dnr(ii)%domain_no,real_fmt=real_int)
        do mm=1, 6
            do nn=1, 6             ! Acts like a header sorting algorithm
               call output_CR0%add(Smat_dnr_CR0(ii)%Smat(nn,mm),real_fmt=real_fmt)
            end do
        end do
        call output_CR0%add(Smat_dnr    (ii)%DoA,real_fmt=real_fmt)
        call output_CR0%add(Smat_dnr_CR0(ii)%DoA,real_fmt=real_fmt)
        call output_CR0%next_row()
    end if
 end Do



ii=lines      ! loop incremented to lines+1 which causes a runtime error 
!----------------------------------------------------------------------------------------------
!-- Print final output
!--
write(*,*)
write(*,'(A)',advance='NO')" Program finished."
write(*,'(I5,A)',advance='YES')zero_matrix_counter," lines contained thresholded tensors."
write(*,*)
write(*,"(A)")  "----------------------------------------------------------------"

end if
end if
end if
end if
end if

!-- lots of potential to simplify

if ( (CR0_F_ex+CR1_F_ex+CR2_F_ex+me_F_ex+in_F_ex ) .ne. 0 ) then
   write(*,'(A)',advance='YES')"---------------------------------------------------------------------------------"
   write(*,*) 
   if ( in_L_ex .eqv. .FALSE. ) then
      write(*,'(A)',advance='YES')" Input-file does not exist."
      write(*,*) 
   end if
   if ( CR0_L_ex .eqv. .TRUE. ) then
       write(*,'(A)',advance='YES')" File containing optimised eff stiffness matrices criteria 0 exists."  
   end if
   if ( CR1_L_ex .eqv. .TRUE. ) then
       write(*,'(A)',advance='YES')" File containing optimised eff stiffness matrices criteria 1 exists."  
   end if
   if ( CR2_L_ex .eqv. .TRUE. ) then
     write(*,'(A)',advance='YES')" File containing optimised eff stiffness matrices criteria 2 exists."
   end if
   if ( me_L_ex .eqv. .TRUE. ) then
       write(*,'(A)',advance='YES') "Meta-data-file exists."
   end if
   write(*,*)
   write(*,*) 
   write(*,'(A)',advance='YES')" Program aborted."
   write(*,'(A)',advance='YES')"---------------------------------------------------------------------------------"
end if
9999 continue

End program opt_and_sort_eff_stiffness

