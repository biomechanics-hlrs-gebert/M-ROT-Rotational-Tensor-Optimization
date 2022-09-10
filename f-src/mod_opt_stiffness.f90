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
! Following variables are declared globally, but mustn't be send via mpi,
! because this may alter the results of the other ranks!
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
REAL(rk), DIMENSION(3) :: dig, dog
REAL(rk), DIMENSION(:,:,:), ALLOCATABLE :: crit

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
!> @param[in] range Range to sweep over
!> @param[in] resolution Distance between steps of searching the optimum
!------------------------------------------------------------------------------ 
SUBROUTINE opt_stiff(mode, range, resolution)
    
CHARACTER(*), INTENT(IN) :: mode
REAL(rk), INTENT(IN):: range, resolution

INTEGER(ik) :: ii, jj, kk, best_position(3), steps

REAL(rk) :: min, alpha, eta, phi, offset
REAL(rk), DIMENSION(6,6) :: tmp_mat, mask

!----------------------------------------------------------------------------------------------
! Initialize variables
!----------------------------------------------------------------------------------------------
mask = 0_ik
min = 10E09_ik

!----------------------------------------------------------------------------------------------
! Copy input to output tensor. Important to keep domain data etc.
!----------------------------------------------------------------------------------------------
tout = tin

!----------------------------------------------------------------------------------------------
! Choose between monotropic or orthotropic optimization
!----------------------------------------------------------------------------------------------
SELECT CASE(TRIM(ADJUSTL(mode)))
    !----------------------------------------------------------------------------------------------
    !
    ! General stiffness matrix, S_ij = S_ji:             Convention: 
    ! 
    ! [ S_11, S_12, S_13, | S_14, S_15, S_16 ]           [            |            ] 
    ! [ S_21, S_22, S_23, | S_24, S_25, S_26 ]           [ Quadrant 2 | Quadrant 1 ] 
    ! [ S_31, S_32, S_33, | S_34, S_35, S_36 ]           [            |            ] 
    ! --------------------|-------------------           --------------------------- 
    ! [ S_41, S_42, S_43, | S_44, S_45, S_46 ]           [            |            ] 
    ! [ S_51, S_52, S_53, | S_54, S_55, S_56 ]           [ Quadrant 3 | Quadrant 4 ] 
    ! [ S_61, S_62, S_63, | S_64, S_65, S_66 ]           [            |            ] 
    !
    ! 
    ! Monotropic, case 1 of 3.                           Monotropic, case 2 of 3.                   
    !                                                                                                                                                                                
    ! [ S_11, S_12, S_13, | S_14,   0 ,   0  ]           [ S_11, S_12, S_13, |   0 , S_15,   0  ] 
    ! [ S_21, S_22, S_23, | S_24,   0 ,   0  ]           [ S_21, S_22, S_23, |   0 , S_25,   0  ] 
    ! [ S_31, S_32, S_33, | S_34,   0 ,   0  ]           [ S_31, S_32, S_33, |   0 , S_35,   0  ] 
    ! --------------------|-------------------           --------------------|------------------- 
    ! [ S_41, S_42, S_43, | S_44,   0 ,   0  ]           [   0 ,   0 ,   0 , | S_44,   0 , S_46 ] 
    ! [   0 ,   0 ,   0 , |   0 , S_55, S_56 ]           [ S_51, S_52, S_53, |   0 , S_55,   0  ] 
    ! [   0 ,   0 ,   0 , |   0 , S_65, S_66 ]           [   0 ,   0 ,   0 , | S_64,   0 , S_66 ] 
    !
    ! Monotropic, case 3 of 3.                           Orthotropic:  
    !                                                                                                 
    ! [ S_11, S_12, S_13, |   0 ,   0 , S_16 ]           [ S_11, S_12, S_13, |   0 ,   0 ,   0  ] 
    ! [ S_21, S_22, S_23, |   0 ,   0 , S_26 ]           [ S_21, S_22, S_23, |   0 ,   0 ,   0  ] 
    ! [ S_31, S_32, S_33, |   0 ,   0 , S_36 ]           [ S_31, S_32, S_33, |   0 ,   0 ,   0  ] 
    ! --------------------|-------------------           --------------------|------------------- 
    ! [   0 ,   0 ,   0 , | S_44, S_45,   0  ]           [   0 ,   0 ,   0 , | S_44,   0 ,   0  ] 
    ! [   0 ,   0 ,   0 , | S_54, S_55,   0  ]           [   0 ,   0 ,   0 , |   0 , S_55,   0  ] 
    ! [ S_61, S_62, S_63, |   0 ,   0 , S_66 ]           [   0 ,   0 ,   0 , |   0 ,   0 , S_66 ] 
    !
    !----------------------------------------------------------------------------------------------
    ! At the monotropic and the orthotropic optimization, the zero entries are minimized. 
    ! 
    ! Example orthotropic: The 1st and the 3rd quadrant must be 0.
    ! Additionally, the minor diagonals of the 4th quadrant must be 0. Therefore, each of the zero-
    ! entries occurs in the upper right and the lower left triangle of the R6x6 matrix. 
    ! Subsequently, the mask multiplies them by 2_ik.
    !
    ! The matrix mask needs no recurring initialization, as the values mustn't change during 
    ! optimization (program runtime in general).
    !----------------------------------------------------------------------------------------------
    CASE('mono')
        mask(1:4,5:6) = 2_ik

    CASE('orth')
        ! Explanation @ opt_stiff_mono
        mask(1:3,4:6) = 2_ik
        mask(  4,5:6) = 2_ik
        mask(  5,  6) = 2_ik

    CASE('an1')
        ! S11 available only once.
        mask(1,1) = 1_ik

    CASE('an2')
        ! S_11, S_22, S_33 available only once.
        mask(1,1) = 1_ik
        mask(2,2) = 1_ik
        mask(3,3) = 1_ik

    CASE DEFAULT
        mssg = "No valid optimization chosen. Check your implementation!"
        CALL print_err_stop(std_out, mssg, 1)

END SELECT

!----------------------------------------------------------------------------------------------
! alpha = input angle - 1/2 of the swept angles.
!----------------------------------------------------------------------------------------------
offset = (range / 2._rk)
steps  = CEILING(range/resolution) ! not to miss a step

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
    IF((SIZE(crit, DIM=1) /= steps) .OR. &
       (SIZE(crit, DIM=2) /= steps) .OR. &
       (SIZE(crit, DIM=3) /= steps)) THEN

        DEALLOCATE(crit)
    END IF
END IF

IF(.NOT. ALLOCATED(crit)) THEN
    ALLOCATE(crit(steps, steps, steps))
END IF 

phi = dig(3) - offset
DO kk = 1_ik, steps

    ! IF ((phi < -45)  .OR. (phi > 45)) CYCLE

    eta = dig(2) - offset
    DO jj = 1_ik, steps

        ! IF ((eta < -45)  .OR. (eta > 45)) CYCLE

        alpha = dig(1) - offset
        DO ii = 1_ik, steps

            ! IF ((alpha < -45)  .OR. (alpha > 45)) CYCLE

            CALL transpose_mat (tin%mat, [ alpha, eta, phi ] , tmp_mat)

            !-------------------------------------------------------------------------------
            ! Criteria: Entries that are not zero after multiplication with mask are 
            ! squared and summed up. Integer power (**2) is quicker than real power (**2.0).
            ! https://twitter.com/fortrantip/status/1478765410405298176?s=24
            !-------------------------------------------------------------------------------
            crit(ii, jj, kk) = SUM((tmp_mat * mask)**2_ik)  

        alpha = alpha + resolution
        END DO

    eta = eta + resolution
    END DO

phi = phi + resolution
END DO

!-------------------------------------------------------------------------------
! Find the best position (min or max) within the criteria space.
!-------------------------------------------------------------------------------
SELECT CASE(TRIM(ADJUSTL(mode)))
    CASE('mono', 'orth'); best_position = MINLOC(crit(:, :, :))
    CASE('an1', 'an2'); best_position = MAXLOC(crit(:, :, :))
END SELECT

!-------------------------------------------------------------------------------
! Calculate the angles of the best position again. 
! Alternatively the angles can be stored in an array like crit(:,:,:), but with 
! REAL64. But storing this array and the if/else within the nested loops of
! ii, jj, kk will consume lots of memory and compute time within the array.
!-------------------------------------------------------------------------------
tout%pos = dig(:) - (range / 2._rk) + (resolution * best_position(:))

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
    REAL(rk), DIMENSION(6,6), INTENT(in) :: mat
    REAL(rk) :: res

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
    REAL(rk), DIMENSION(6,6), INTENT(in) :: mat
    REAL(rk) :: res

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

REAL(rk), DIMENSION(6,6), INTENT(INOUT)  :: mat ! input tensor

REAL(rk), DIMENSION(6,6) :: BB1, BB2, BB3, tmp6x6
REAL(rk), DIMENSION(3)   :: n1, n2, n3
REAL(rk), DIMENSION(3,3) :: aa1, aa2, aa3

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
