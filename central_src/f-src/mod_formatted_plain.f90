!------------------------------------------------------------------------------
! MODULE: formatted_plain
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
! @description: 
!> "Formatted I/O for plain ascii files that are not related to PureDat, Meta, 
!> and log/mon messages. 
!------------------------------------------------------------------------------
MODULE formatted_plain

USE ISO_FORTRAN_ENV
USE global_std
USE mechanical
USE strings
USE user_interaction
USE math

IMPLICIT NONE

INTERFACE write_matrix
  MODULE PROCEDURE write_matrix_int
  MODULE PROCEDURE write_matrix_real 
END INTERFACE write_matrix

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: write_matrix_real
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
!
!> @brief
!> Subroutine to print regular tensors respectively matrices.
!
!> @Description
!> Please provide mat_real OR mat_in :-)
!> Automatically writes "sym" if hide_zeros .EQV. .TRUE. and left triangle = 0
!> Hide zeroes is set as default.
!> Accepted formats: 'std'/'standard' for scientific formatting and
!> 'spl'/'simple' for traditional formatting
!
!> @param[in] fh Handle of file to print to
!> @param[in] name Name of the object to print
!> @param[in] fmt Formatting of the data
!> @param[in] unit Physical unit of the information to print
!> @param[in] mat Actual matrix
!------------------------------------------------------------------------------
SUBROUTINE write_matrix_real(fh, name, fmt, unit, mat)

INTEGER(KIND=INT64), INTENT(IN) :: fh   
CHARACTER(LEN=*), INTENT(IN) :: name 
REAL   (KIND=rk), DIMENSION(:, :), INTENT(IN) :: mat    
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fmt 
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: unit 

! Internal variables 
INTEGER(KIND=ik)   :: prec , fw, nm_fmt_lngth, ii, jj, kk, dim1, dim2
CHARACTER(LEN=mcl) :: fmt_a, sep, nm_fmt
CHARACTER(LEN=mcl) :: text
CHARACTER(LEN=mcl) :: fmt_u
REAL(KIND=rk) :: sym_out
LOGICAL :: sym_u

!------------------------------------------------------------------------------
! Initialize and check for presence of the variables
!------------------------------------------------------------------------------
dim1 = SIZE(mat, 1)
dim2 = SIZE(mat, 2)
fmt_u = 'standard'
mssg='' 
text = ''

prec = PRECISION(mat)
fw = prec+8
sym_u = .FALSE.

IF (PRESENT(unit)) THEN
    IF (unit /= '') text = " Unit: ("//TRIM(unit)//")"
END IF

IF (dim1 == dim2) THEN
    CALL check_sym(INT(fh, KIND=ik), mat, name, sym_out=sym_out)
    sym_u = .TRUE.
END IF

!------------------------------------------------------------------------------
! Generate formats
!------------------------------------------------------------------------------
IF(PRESENT(fmt)) fmt_u = fmt

SELECT CASE (TRIM(fmt_u))
   CASE ('std', 'standard')
        WRITE(fmt_a, "(3(A,I0),A)") "(",dim2,"(E",fw,".",prec,"E2))"
        WRITE(sep  , "(A,I0,A)")    "(",fw*dim2,"('-'))"

        ! Calculate text and unit length. If name to long - overflow formaming to the right
        nm_fmt_lngth  = fw*dim2-4-LEN_TRIM(name)-LEN_TRIM(text)

   CASE ('spl', 'simple')
        WRITE(fmt_a,  "(3(A,I0),A)") "(",dim2,"(F10.3))"
        WRITE(sep  ,  "(A,I0,A)")    "(",dim2*10,"('-'))"        

        ! Calculate text and unit length. If name to long - overflow formaming to the right
        nm_fmt_lngth  = dim2*10-4-2-LEN_TRIM(name)-LEN_TRIM(text) 

   CASE('wxm', 'wxmaxima')
       WRITE(fmt_a, '(5(A,I0),A)')  "(' [',",dim2-1,"(E",fw,".",prec,"E2,','),E",fw,".",prec,"E2,'],' )"

       WRITE(fh,"(A,A)")TRIM(name),": matrix("

       DO kk = 1, dim1 - 1
          WRITE(fh, fmt_a) mat(kk,:)
       END DO

       WRITE(fmt_a,'(5(A,I0),A)')  "(' [',",dim2-1,"(E",fw,".",prec,"E2,','),E",fw,".",prec,"E2,']);' )"

       WRITE(fh, fmt_a) mat(dim1, :)
END SELECT

IF (nm_fmt_lngth .LT. 1_ik) nm_fmt_lngth = 1_ik
WRITE(nm_fmt, "(A,I0,A)")  "(2('-') ,3A,", nm_fmt_lngth ,"('-'), A)"    



!------------------------------------------------------------------------------
! Write output
!------------------------------------------------------------------------------
WRITE(fh, '(A)')
WRITE(fh, sep)                                    ! Separator
WRITE(fh, nm_fmt) ' ',TRIM(name), ' ', TRIM(text) ! Named separator

DO ii=1, dim1
DO jj=1, dim2

    IF ((sym_u) .AND. (ii==dim1) .AND. (jj==1)) THEN
        SELECT CASE(fmt_u)
        CASE('spl', 'simple')
            WRITE(fh, '(A)', ADVANCE='NO') "symmetric "
        CASE('std', 'standard')
            WRITE(fh, '(A)', ADVANCE='NO') "   symmetric           "
        END SELECT

    ELSE IF ((sym_u) .AND. (ii==dim1) .AND. (jj==2)) THEN
        IF (ABS(sym_out) <=  10E-08) sym_out = 0._rk
        WRITE(fh, fmt_a, ADVANCE='NO') sym_out          
    ELSE
        IF ((ABS(mat(ii,jj)) >=  10E-08) .AND. ((.NOT. sym_u) .OR. ((sym_u) .AND. (jj .GE. ii)))) THEN 
            WRITE(fh, fmt_a, ADVANCE='NO') mat (ii,jj)
        ELSE
            SELECT CASE(fmt_u)
            CASE('spl', 'simple')
                WRITE(fh, '(A)', ADVANCE='NO') "      .   "
            CASE('std', 'standard')
                WRITE(fh, '(A)', ADVANCE='NO') "   .                   "
            END SELECT
        END IF
    END IF

END DO
WRITE(fh,'(A)') ''
END DO        

WRITE(fh, '(A)') ''                               ! Newline & Carriage return
End Subroutine write_matrix_real


!------------------------------------------------------------------------------
! SUBROUTINE: write_matrix_int
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
!
!> @brief
!> Subroutine to print regular tensors respectively matrices.
!
!> @Description
!> Please provide mat_real OR mat_in :-)
!> Accepted formats: 'std'/'standard' for scientific formatting and
!> 'spl'/'simple' for traditional formatting
!
!> @param[in] fh Handle of file to print to
!> @param[in] name Name of the object to print
!> @param[in] fmt Formatting of the data
!> @param[in] unit Physical unit of the information to print
!> @param[in] mat Actual matrix
!------------------------------------------------------------------------------
SUBROUTINE write_matrix_int(fh, name, fmt, unit, mat)

INTEGER(KIND=INT64), INTENT(IN) :: fh   
CHARACTER(LEN=*), INTENT(IN) :: name 
INTEGER(KIND=INT64), DIMENSION(:, :), INTENT(IN) :: mat    
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fmt 
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: unit 

! Internal variables 
INTEGER(KIND=ik)   :: nm_fmt_lngth, ii, jj, dim1, dim2
CHARACTER(LEN=mcl) :: fmt_a, sep, nm_fmt
CHARACTER(LEN=mcl) :: text
CHARACTER(LEN=mcl) :: fmt_u
REAL(KIND=rk) :: sym_out
LOGICAL :: sym_u

!------------------------------------------------------------------------------
! Initialize and check for presence of the variables
!------------------------------------------------------------------------------
dim1 = SIZE(mat, 1)
dim2 = SIZE(mat, 2)
fmt_u = 'standard'
sym_u = .FALSE.
mssg='' 
text = ''

IF (PRESENT(unit)) THEN
    IF (unit /= '') text = " Unit: ("//TRIM(unit)//")"
END IF

IF (dim1 == dim2) THEN
    ! ik for file handle required to support different standard iks.
    CALL check_sym(INT(fh, KIND=ik), REAL(mat, KIND=rk), name, sym_out=sym_out)
    sym_u = .TRUE.
END IF

!------------------------------------------------------------------------------
! Generate formats
!------------------------------------------------------------------------------
IF(PRESENT(fmt)) fmt_u = fmt

SELECT CASE (TRIM(fmt_u))
   CASE ('std', 'standard')

        WRITE(sep, "(A,I0,A)") "(",dim2,"('-'))"

        ! Calculate text and unit length. If name to long - overflow formatting to the right
        nm_fmt_lngth = dim2-4-2-LEN_TRIM(name)-LEN_TRIM(text)

   CASE ('spl', 'simple')

        WRITE(sep, "(A,I0,A)") "(",dim2*10,"('-'))"        

        ! Calculate text and unit length. If name to long - overflow formaming to the right
        nm_fmt_lngth  = dim2*10-4-LEN_TRIM(name)-LEN_TRIM(text) 
END SELECT

IF (nm_fmt_lngth .LT. 1_ik) nm_fmt_lngth = 1_ik
WRITE(nm_fmt, "(A,I0,A)")  "(2('-') ,3A,", nm_fmt_lngth ,"('-'), A)"    

WRITE(fmt_a, "(3(A,I0),A)") "(",dim2,"(I10))"

!------------------------------------------------------------------------------
! Write output
!------------------------------------------------------------------------------
WRITE(fh, sep)                                    ! Separator
WRITE(fh, nm_fmt) ' ',TRIM(name), ' ', TRIM(text) ! Named separator

DO ii=1, dim1
DO jj=1, dim2

    IF ((sym_u) .AND. (ii==dim1) .AND. (jj==1)) THEN
        WRITE(fh, '(A)', ADVANCE='NO') " symmetric"
    ELSE IF ((sym_u) .AND. (ii==dim1) .AND. (jj==2)) THEN
        IF (ABS(sym_out) <=  10E-08) sym_out = 0._rk
        WRITE(fh, fmt_a, ADVANCE='NO') sym_out          
    ELSE
        IF ((mat(ii,jj) /= 0) .AND. ((.NOT. sym_u) .OR. ((sym_u) .AND. (jj .GE. ii)))) THEN 
            WRITE(fh, fmt_a, ADVANCE='NO') mat (ii,jj)
        ELSE
            WRITE(fh, '(A)', ADVANCE='NO') "         ."
        END IF
    END IF

END DO
WRITE(fh,'(A)') ''
END DO        

WRITE(fh, '(A)') ''
End Subroutine write_matrix_int


!------------------------------------------------------------------------------
! SUBROUTINE: underscore_to_blank
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Replaces underscores with blanks.
!
!> @param[in] instring Input string
!> @param[out] outstring Output string
!------------------------------------------------------------------------------  
SUBROUTINE underscore_to_blank (instring, outstring)
  ! This whole subroutine is a workaround :-)
  CHARACTER(LEN=*) :: instring
  CHARACTER(LEN=*) :: outstring
  INTEGER(KIND=ik) :: ii

  outstring=instring
  DO ii=1, LEN_TRIM(instring)
    IF (instring(ii:ii) == '_')  outstring(ii:ii) = ' '
  END DO

  outstring=ADJUSTL(TRIM(outstring))
END SUBROUTINE underscore_to_blank

!------------------------------------------------------------------------------
! SUBROUTINE: check_tensor_2nd_rank_R66_header
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Validate the header of a tensor_2nd_rank_R66 file
!
!> @param[in] header Header of the tensor_2nd_rank_R66 file
!> @param[out] abort Feedback whether to stop the program
!------------------------------------------------------------------------------  
SUBROUTINE check_tensor_2nd_rank_R66_header(header, abort)

CHARACTER(LEN=*), INTENT(IN) :: header
LOGICAL, INTENT(OUT) :: abort

CHARACTER(LEN=mcl) :: tokens(100)
INTEGER(KIND=ik) :: ntokens

abort = .FALSE.

CALL parse(TRIM(ADJUSTL(header)), ",", args=tokens, nargs=ntokens)

IF(ntokens /= 42_ik) abort = .TRUE.

!------------------------------------------------------------------------------
! Implementation might look silly. But calculating the indices during runtime, 
! to translate them into string and compare them with an if-function or 
! parsing the strings to numbers for comparison will consume way more time.
!------------------------------------------------------------------------------
IF(TRIM(ADJUSTL(tokens( 1))) /= "Domain")    abort = .TRUE.
IF(TRIM(ADJUSTL(tokens( 2))) /= "Density")   abort = .TRUE.
IF(TRIM(ADJUSTL(tokens( 3))) /= "DoA")       abort = .TRUE.
IF(TRIM(ADJUSTL(tokens( 4))) /= "pos_alpha") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens( 5))) /= "pos_eta")   abort = .TRUE.
IF(TRIM(ADJUSTL(tokens( 6))) /= "pos_phi")   abort = .TRUE.
IF(TRIM(ADJUSTL(tokens( 7))) /= "S11") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens( 8))) /= "S21") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens( 9))) /= "S31") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(10))) /= "S41") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(11))) /= "S51") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(12))) /= "S61") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(13))) /= "S12") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(14))) /= "S22") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(15))) /= "S32") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(16))) /= "S42") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(17))) /= "S52") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(18))) /= "S62") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(19))) /= "S13") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(20))) /= "S23") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(21))) /= "S33") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(22))) /= "S43") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(23))) /= "S53") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(24))) /= "S63") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(25))) /= "S14") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(26))) /= "S24") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(27))) /= "S34") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(28))) /= "S44") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(29))) /= "S54") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(30))) /= "S64") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(31))) /= "S15") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(32))) /= "S25") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(33))) /= "S35") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(34))) /= "S45") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(35))) /= "S55") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(36))) /= "S65") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(37))) /= "S16") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(38))) /= "S26") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(39))) /= "S36") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(40))) /= "S46") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(41))) /= "S56") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(42))) /= "S66") abort = .TRUE.

END SUBROUTINE check_tensor_2nd_rank_R66_header

!------------------------------------------------------------------------------
! SUBROUTINE: parse_tensor_2nd_rank_R66_row
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Parse the header of a tensor_2nd_rank_R66 file
!
!> @description
!> Function returns a negative first column (domain numer) if the data 
!> within a row does not meet the only valid amount of entries.
!
!> @param[in] row An unparsed row of the tensor_2nd_rank_R66 file
!> @param[out] tensor_of_row The resulting TYPE tensor_2nd_rank_R66 file
!------------------------------------------------------------------------------  
SUBROUTINE parse_tensor_2nd_rank_R66_row(row, tensor_of_row)

CHARACTER(LEN=*), INTENT(IN) :: row
TYPE(tensor_2nd_rank_R66), INTENT(OUT) :: tensor_of_row

CHARACTER(LEN=mcl) :: tokens(100)
INTEGER(KIND=ik) :: ntokens
LOGICAL :: invalid

invalid = .FALSE.

CALL parse(TRIM(ADJUSTL(row)), ",", args=tokens, nargs=ntokens)

IF(ntokens /= 42_ik) invalid = .TRUE.

IF(invalid) THEN
    tensor_of_row%dmn     = -INT(tokens(1), KIND=ik)
    tensor_of_row%density = 0._rk
    tensor_of_row%Doa     = 0._rk
    tensor_of_row%pos     = 0._rk
    tensor_of_row%mat     = 0._rk    
ELSE
    READ(tensor_of_row%dmn,      '(I0)'    ) tokens( 1)
    READ(tensor_of_row%density,  '(F30.15)') tokens( 2)
    READ(tensor_of_row%DoA,      '(F30.15)') tokens( 3)
    READ(tensor_of_row%pos(1),   '(F30.15)') tokens( 4)
    READ(tensor_of_row%pos(2),   '(F30.15)') tokens( 5)
    READ(tensor_of_row%pos(3),   '(F30.15)') tokens( 6)
    READ(tensor_of_row%mat(1,1), '(F30.15)') tokens( 7)
    READ(tensor_of_row%mat(2,1), '(F30.15)') tokens( 8)
    READ(tensor_of_row%mat(3,1), '(F30.15)') tokens( 9)
    READ(tensor_of_row%mat(4,1), '(F30.15)') tokens(10)
    READ(tensor_of_row%mat(5,1), '(F30.15)') tokens(11)
    READ(tensor_of_row%mat(6,1), '(F30.15)') tokens(12)
    READ(tensor_of_row%mat(1,2), '(F30.15)') tokens(13)
    READ(tensor_of_row%mat(2,2), '(F30.15)') tokens(14)
    READ(tensor_of_row%mat(3,2), '(F30.15)') tokens(15)
    READ(tensor_of_row%mat(4,2), '(F30.15)') tokens(16)
    READ(tensor_of_row%mat(5,2), '(F30.15)') tokens(17)
    READ(tensor_of_row%mat(6,2), '(F30.15)') tokens(18)
    READ(tensor_of_row%mat(1,3), '(F30.15)') tokens(19)
    READ(tensor_of_row%mat(2,3), '(F30.15)') tokens(20)
    READ(tensor_of_row%mat(3,3), '(F30.15)') tokens(21)
    READ(tensor_of_row%mat(4,3), '(F30.15)') tokens(22)
    READ(tensor_of_row%mat(5,3), '(F30.15)') tokens(23)
    READ(tensor_of_row%mat(6,3), '(F30.15)') tokens(24)
    READ(tensor_of_row%mat(1,4), '(F30.15)') tokens(25)
    READ(tensor_of_row%mat(2,4), '(F30.15)') tokens(26)
    READ(tensor_of_row%mat(3,4), '(F30.15)') tokens(27)
    READ(tensor_of_row%mat(4,4), '(F30.15)') tokens(28)
    READ(tensor_of_row%mat(5,4), '(F30.15)') tokens(29)
    READ(tensor_of_row%mat(6,4), '(F30.15)') tokens(30)
    READ(tensor_of_row%mat(1,5), '(F30.15)') tokens(31)
    READ(tensor_of_row%mat(2,5), '(F30.15)') tokens(32)
    READ(tensor_of_row%mat(3,5), '(F30.15)') tokens(33)
    READ(tensor_of_row%mat(4,5), '(F30.15)') tokens(34)
    READ(tensor_of_row%mat(5,5), '(F30.15)') tokens(35)
    READ(tensor_of_row%mat(6,5), '(F30.15)') tokens(36)
    READ(tensor_of_row%mat(1,6), '(F30.15)') tokens(37)
    READ(tensor_of_row%mat(2,6), '(F30.15)') tokens(38)
    READ(tensor_of_row%mat(3,6), '(F30.15)') tokens(39)
    READ(tensor_of_row%mat(4,6), '(F30.15)') tokens(40)
    READ(tensor_of_row%mat(5,6), '(F30.15)') tokens(41)
    READ(tensor_of_row%mat(6,6), '(F30.15)') tokens(42)
END IF

END SUBROUTINE parse_tensor_2nd_rank_R66_row

!------------------------------------------------------------------------------
! SUBROUTINE: parse_tensor_2nd_rank_R66
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Parse a "csv" file with the structure of a tensor_2nd_rank_R66 data type.
!
!> @param[in] fh Input file handle (file existane was checked before).
!> @param[in] filename Filename of the input file
!> @param[out] tensors_in Array of all input tensors
!------------------------------------------------------------------------------  
SUBROUTINE parse_tensor_2nd_rank_R66(fh, filename, amnt_lines, tensors_in)

INTEGER(KIND=ik), INTENT(IN) :: fh, amnt_lines
CHARACTER(LEN=*), INTENT(IN) :: filename
TYPE(tensor_2nd_rank_R66), DIMENSION(:), INTENT(OUT) :: tensors_in

!------------------------------------------------------------------------------
! Line of csv may be pretty long (>40 floats)
!------------------------------------------------------------------------------
CHARACTER(LEN=10_ik*mcl) :: header, line
CHARACTER(LEN=mcl) :: tokens(100)
INTEGER(KIND=ik) :: lines, ii, ntokens, entry_no
LOGICAL :: abort

!------------------------------------------------------------------------------
! Parse the header of the "csv" data
! Error if the header does not indicate 42 entries (determined by TYPE)
!------------------------------------------------------------------------------
READ(fh, '(A)') header
CALL check_tensor_2nd_rank_R66_header(header, abort)

IF(abort) THEN
	mssg = "The input file "//TRIM(ADJUSTL(filename))//&
		" is not a valid 'tensor_2nd_rank_R66' file."
	CALL print_err_stop(fh, mssg, 1)
END IF

DO ii = 2_ik, amnt_lines
	READ(fh, '(A)') line

	CALL parse_tensor_2nd_rank_R66_row(line, tensors_in(ii-1_ik))
END DO

END SUBROUTINE parse_tensor_2nd_rank_R66

END MODULE formatted_plain
