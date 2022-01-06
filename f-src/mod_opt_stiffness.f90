!------------------------------------------------------------------------------
! MODULE: opt_stiffness
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
! @description: 
!> Module for optimizing R6x6 anisotropic tensors.
!------------------------------------------------------------------------------
MODULE opt_stiffness

USE global_std
USE mechanical
USE math
USE user_interaction

IMPLICIT NONE

!------------------------------------------------------------------------------
! Following variables are declared globally, but mustn't be send via mpi(!)
! Declaring them globally makes calling functions and exporting the computed
! criteria space of a tensor easier. 
!------------------------------------------------------------------------------  
! tin            Input tensor
! dig            Input angles (Often = 0)
! pm_steps       +/- steps of the angle to check with.
! intervall      stps +/- intervall = range of optimization of this call.
! tout           Output tensor
! dog            Output angles
!------------------------------------------------------------------------------  
TYPE(tensor_2nd_rank_R66) :: tin, tout
INTEGER(KIND=ik), DIMENSION(3) :: pm_steps
REAL(KIND=rk)   , DIMENSION(3) :: intervall, dig, dog
REAL(KIND=rk)   , DIMENSION(:,:,:), ALLOCATABLE :: crit

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: opt_stiff
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @Brief:
!> Rotate around all 3 rot dof until objective function is min. 
!> Initial degrees must be integer!
!
!> @param[in] mode "monotropic" or "orthotropic"
!------------------------------------------------------------------------------ 
SUBROUTINE opt_stiff(mode)

CHARACTER(LEN=*), INTENT(IN) :: mode

INTEGER(KIND=ik), DIMENSION(3) :: ttl_steps
INTEGER(KIND=ik) :: ii, jj, kk

REAL(KIND=rk) :: min, alpha, eta, phi
REAL(KIND=rk), DIMENSION(6,6) :: tmp_r6x6, mask

!----------------------------------------------------------------------------------------------
! Initialize variables
!----------------------------------------------------------------------------------------------
ttl_steps = 0_ik
mask = 0_ik
min = 10E09_ik

!----------------------------------------------------------------------------------------------
! Total amount of steps: requested angle +/- steps --> 2*steps + value --> (pm_steps*2)+1
!----------------------------------------------------------------------------------------------
ttl_steps = (pm_steps*2_ik)+1_ik

!----------------------------------------------------------------------------------------------
! Copy input to output tensor. Important to keep domain data etc.
!----------------------------------------------------------------------------------------------
tout = tin

!----------------------------------------------------------------------------------------------
! Choose between monotropic or orthotropic optimization
!----------------------------------------------------------------------------------------------
SELECT CASE(TRIM(ADJUSTL(mode)))
    CASE('monotropic')
    !----------------------------------------------------------------------------------------------
    ! After multiplication with mask, the 2nd rank R6x6 tensor entries that are zero in an 
    ! monotropic scenario are given like the tensor was given to the call opt_stiff.
    !
    ! The values of the mask are 2 instead of 1, because only the first and the fourth quadrant
    ! of the matrix are taken into account.
    !
    ! The matrix mask needs no recurring initialization, as the values mustn't change during 
    ! optimization (program runtime in general).
    !----------------------------------------------------------------------------------------------
    mask(1:4,5:6) = 2_ik

CASE('orthotropic')
    ! Explanation @ opt_stiff_mono
    mask(1:3,4:6) = 2_ik
    mask(  4,5:6) = 2_ik
    mask(  5,  6) = 2_ik

CASE('anisotropic1')
    ! S11 available only once.
    mask(1,1) = 1_ik

CASE('anisotropic2')
    ! S11, S22, S44 available only once.
    mask(1,1) = 1_ik
    mask(2,2) = 1_ik
    mask(3,3) = 1_ik

CASE DEFAULT
    mssg = "No valid optimization chosen. Check your implementation!"
    CALL print_err_stop(std_out, mssg, 1)

END SELECT

!----------------------------------------------------------------------------------------------
! Allocate memory for angle +/- steps.
! Each call, the system checks whether the space of criteria is allocated.
! If not --> allocate
! If ttl_steps changes --> reallocate
! If Program closes --> Everything alright
! Nothing gets deallocated without explicit order. This feature enables writing the 
! space of criteria to a vtk (i.e., scalar field). 
! @MPI: With mpi, crit may have a different state across the ranks!
!----------------------------------------------------------------------------------------------
IF(ALLOCATED(crit)) THEN
    IF((SIZE(crit, DIM=1) /= ttl_steps(1)) .OR. &
       (SIZE(crit, DIM=2) /= ttl_steps(2)) .OR. &
       (SIZE(crit, DIM=3) /= ttl_steps(3))) THEN

        DEALLOCATE(crit)
    END IF
END IF

IF(.NOT. ALLOCATED(crit)) THEN
    ALLOCATE(crit(ttl_steps(1), ttl_steps(2), ttl_steps(3)))
END IF

alpha = dig(1) - (intervall(1) * pm_steps(1))
DO kk = 1_ik, ttl_steps(1)

    eta = dig(2) - (intervall(2) * pm_steps(2))
    DO jj = 1_ik, ttl_steps(2)

        phi = dig(3) - (intervall(3) * pm_steps(3))
        DO ii = 1_ik, ttl_steps(3)

            CALL transpose_mat (tin%mat, [ alpha, eta, phi ] , tmp_r6x6)

            !-------------------------------------------------------------------------------
            ! Criteria: Entries that are not zero after multiplication with mask are 
            ! squared and summed up. Integer power (**2) is quicker than real power (**2.0).
            ! https://twitter.com/fortrantip/status/1478765410405298176?s=24
            !-------------------------------------------------------------------------------
            crit(kk, jj, ii) = SUM((tmp_r6x6 * mask)**2_ik)  

            !-------------------------------------------------------------------------------
            ! Update the best position of the output matrix if the criteria is below the 
            ! current minimum.
            !-------------------------------------------------------------------------------
            IF ( min > crit(kk, jj, ii)) THEN

                min = crit(kk, jj, ii)

                tout%pos = [ alpha, eta, phi ]
            END IF

        phi = phi + intervall(3)
        END DO

    eta = eta + intervall(2)
    END DO

alpha = alpha + intervall(1)
END DO

!----------------------------------------------------------------------------------------------
! Recalculate the best combination of angles.
!----------------------------------------------------------------------------------------------
CALL transpose_mat (tin%mat, tout%pos, tout%mat)

END SUBROUTINE opt_stiff


!------------------------------------------------------------------------------
! FUNCTION: DoAM
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief:
!> Degree of anisotropy - monotropic
!
!> @param[in] mat 2nd rank R6x6 tensor
!> @return[out] res Result
!------------------------------------------------------------------------------
FUNCTION DoAM(mat) RESULT(res)
    REAL(KIND=rk), DIMENSION(6,6), INTENT(in) :: mat
    REAL(KIND=rk) :: res

    !------------------------------------------------------------------------------
    !-- Degree of Anisotropy - average 0 divided by average non-zero entry
    !------------------------------------------------------------------------------
    res= (( + ABS(mat(1,5)) +    ABS(mat(1,6)) &
            + ABS(mat(2,5)) +    ABS(mat(2,6)) &
            + ABS(mat(3,5)) +    ABS(mat(3,6)) &
            + ABS(mat(4,5)) +    ABS(mat(4,6)) )/8.0_rk ) / &
            ((ABS(mat(1,1)) + 2* ABS(mat(1,2)) + 2* ABS(mat(1,3)) + 2* ABS(mat(1,4)) &
            + ABS(mat(2,2)) + 2* ABS(mat(2,3)) + 2* ABS(mat(2,4)) &
            + ABS(mat(3,3)) + 2* ABS(mat(3,4)) &
            + ABS(mat(4,4)) &
            + ABS(mat(5,5)) + 2* ABS(mat(5,6)) &
            + ABS(mat(6,6)))/20.0_rk)*100.0_rk
END FUNCTION DoAM          

!------------------------------------------------------------------------------
! FUNCTION: DoAO
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief:
!> Degree of anisotropy - orthotropic
!
!> @param[in] mat 2nd rank R6x6 tensor
!> @return[out] res Result
!------------------------------------------------------------------------------
FUNCTION DoAO(mat) RESULT(res)
    REAL(KIND=rk), DIMENSION(6,6), INTENT(in) :: mat
    REAL(KIND=rk) :: res

    !------------------------------------------------------------------------------
    !-- Degree of Anisotropy - average 0 divided by average non-zero entry
    !------------------------------------------------------------------------------
    res=(( ABS(mat(1,4)) + ABS(mat(1,5)) + ABS(mat(1,6)) &
        + ABS(mat(2,4)) + ABS(mat(2,5)) + ABS(mat(2,6)) &
        + ABS(mat(3,4)) + ABS(mat(3,5)) + ABS(mat(3,6)) &
        + ABS(mat(4,5)) + ABS(mat(4,6)) &
        + ABS(mat(5,6)) )/12.0_rk ) / &
        ((ABS(mat(1,1)) +  2*ABS(mat(1,2)) +  2*ABS(mat(1,3)) &
        + ABS(mat(2,2)) +  2*ABS(mat(2,3)) &
        + ABS(mat(3,3)) &
        + ABS(mat(4,4)) + ABS(mat(5,5)) + ABS(mat(6,6)) )/12.0_rk)*100.0_rk
END FUNCTION DoAO

!------------------------------------------------------------------------------
! SUBROUTINE: tilt_tensor
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Tilts a tin by 90° around axis x, y or z - respectively 1, 2 or 3
!> other angles are not implemented. Changes S11, S22 and S33 to get 
!> an S11 > S22 > S33 ordering.
!
!> @param[in] mat Input tensor 
!------------------------------------------------------------------------------  
SUBROUTINE tilt_tensor(mat)

REAL(KIND=rk), DIMENSION(6,6), INTENT(INOUT)  :: mat ! input tensor

REAL(KIND=rk), DIMENSION(6,6) :: BB1, BB2, BB3, tmp6x6
REAL(KIND=rk), DIMENSION(3)   :: n1, n2, n3
REAL(KIND=rk), DIMENSION(3,3) :: aa1, aa2, aa3

n1 = [1,0,0]
n2 = [0,1,0]
n3 = [0,0,1]
aa1 = rot_alg (n1, pihalf)
aa2 = rot_alg (n2, pihalf)
aa3 = rot_alg (n3, pihalf)
BB1 = tra_R6 (aa1)
BB2 = tra_R6 (aa2)
BB3 = tra_R6 (aa3)

     IF((mat(1,1)>mat(2,2)) .AND. (mat(1,1)>mat(3,3)) .AND. (mat(2,2)>mat(3,3))) THEN
    mat = mat !-- 123
ELSE IF((mat(1,1)>mat(2,2)) .AND. (mat(1,1)>mat(3,3)) .AND. (mat(2,2)<mat(3,3))) THEN
    mat  = MATMUL(MATMUL(TRANSPOSE(BB1),mat),BB1) !-- 132 => 123
ELSE IF((mat(1,1)<mat(2,2)) .AND. (mat(1,1)<mat(3,3)) .AND. (mat(2,2)>mat(3,3))) THEN
    tmp6x6 = MATMUL(MATMUL(TRANSPOSE(BB3),mat),BB3) !-- 231 => 132
        
    mat  = MATMUL(MATMUL(TRANSPOSE(BB1),tmp6x6),BB1) !-- 132 => 123
ELSE IF((mat(1,1)<mat(2,2)) .AND. (mat(1,1)>mat(3,3)) .AND. (mat(2,2)>mat(3,3))) THEN
    mat  = MATMUL(MATMUL(TRANSPOSE(BB3),mat),BB3) !-- 213 => 123
ELSE IF((mat(1,1)>mat(2,2)) .AND. (mat(1,1)<mat(3,3)) .AND. (mat(2,2)<mat(3,3))) THEN
    tmp6x6 = MATMUL(MATMUL(TRANSPOSE(BB2),mat),BB2) !-- 312 => 132
        
    mat  = MATMUL(MATMUL(TRANSPOSE(BB1),tmp6x6),BB1) !-- 132 => 123
ELSE IF((mat(1,1)<mat(2,2)) .AND. (mat(1,1)<mat(3,3)) .AND. (mat(2,2)<mat(3,3))) THEN
    mat  = MATMUL(MATMUL(TRANSPOSE(BB2),mat),BB2) !-- 321 => 123
END IF

END SUBROUTINE tilt_tensor

END MODULE opt_stiffness
