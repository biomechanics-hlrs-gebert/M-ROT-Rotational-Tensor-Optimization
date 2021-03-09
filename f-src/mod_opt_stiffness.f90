!----------------------------------------------------------------------------------------------
!--
!-- Module for optimising R6x6 anisotropic tensors
!--
!-- author:                Johannes Gebert        (HLRS, NUM)
!-- numeric kernel by:     Ralf Schneider         (HLRS, NUM)
!--
!-- mod 27.06.2020         Previous programs ignored
!-- mod 06.03.2021
!--
!----------------------------------------------------------------------------------------------

MODULE opt_stiffness

USE standards
USE math_routines

IMPLICIT NONE

CONTAINS

  SUBROUTINE opt_eff_stiff (mode, EE, deg_a, stp_a, deg_p, stp_p, deg_e, stp_e, intervall, &
       opt_tensor, a, p, e, outp)

    !-- Choose monotropic OR otrhotropic optimization
    !-- Initial degrees must be integer!
    !-- Change mode to strings - way more convenient

    INTEGER(KIND=ik), INTENT(in)                                           :: mode                        ! mono=1 and orth=2
    REAL   (KIND=rk), INTENT(in)           , DIMENSION(6,6)                :: EE                          ! Tensor to opimise
    INTEGER(KIND=ik), INTENT(in)                                           :: deg_a, deg_e, deg_p         ! orientation to start with
    INTEGER(KIND=ik), INTENT(in)                                           :: stp_a, stp_e, stp_p         ! Amount of steps to increase
    REAL   (KIND=rk), INTENT(in) , OPTIONAL                                :: intervall                   ! stepwidth
    REAL   (KIND=rk), INTENT(out), OPTIONAL, DIMENSION(6,6)                :: opt_tensor                  ! optimised tensor
    REAL   (KIND=rk), INTENT(out), OPTIONAL                                :: a,p,e                       ! resulting angles
    LOGICAL         , INTENT(in) , OPTIONAL                                :: outp                        ! check whether print result

    !----------------------------------------------------------------------------------------------

    REAL   (KIND=rk)                                                       :: kk,   kk_phi,   kk_eta
    REAL   (KIND=rk)                                                       :: kk_s, kk_phi_s, kk_eta_s    ! stored variables
    INTEGER(KIND=ik)                                                       :: ii, ii_phi, ii_eta, sw1, sw2
    INTEGER(KIND=ik)                                                       :: ii_c, ii_phi_c, ii_eta_c
    REAL   (KIND=rk)                                                       :: alpha, phi, eta
    REAL   (KIND=rk)                                                       :: min, intval
    REAL   (KIND=rk)                       , DIMENSION(3)                  :: n
    REAL   (KIND=rk)                       , DIMENSION(3,3)                :: aa
    REAL   (KIND=rk)                       , DIMENSION(6,6)                :: BB, tmp_r6x6
    REAL   (KIND=rk)                       , DIMENSION(:,:,:), ALLOCATABLE :: crit
    LOGICAL                                                                :: op=.FALSE.
    CHARACTER(len=50)                                                      :: opt_str

    ALLOCATE(crit(stp_a, stp_p, stp_e))

    IF (PRESENT(intervall)) intval = intervall
    IF (PRESENT(outp))          op = outp
    min=10000000_ik        !-- arbitrary initial value
    !----------------------------------------------------------------------------------------------
    !-- DEBUG -------------------------------------------------------------------------------------
    ! write(*,'(A,I5)')"deg_a=", deg_a
    ! write(*,'(A,I5)')"deg_p=", deg_p
    ! write(*,'(A,I5)')"deg_e=", deg_e
    ! write(*,'(A,I5)')"stp_a=", stp_a
    ! write(*,'(A,I5)')"stp_p=", stp_p
    ! write(*,'(A,I5)')"stp_e=", stp_e
    ! write(*,'(A,F5.3)')"intervall=", intervall
    !-- DEBUG -------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------

    kk_eta = deg_e
    DO ii_eta_c = 1, stp_e
       kk_phi = deg_p
       DO ii_phi_c = 1, stp_p
          kk = deg_a
          DO ii_c = 1, stp_a
             CALL transpose_mat (EE, kk, kk_phi, kk_eta, tmp_r6x6,BB)
             IF ( mode .eq. 1_ik ) THEN
                CALL mono_crit (tmp_r6x6,res=crit(ii_c,ii_phi_c,ii_eta_c))
             ELSE IF ( mode .EQ. 2_ik ) THEN
                CALL ortho_crit (tmp_r6x6,crit(ii_c,ii_phi_c,ii_eta_c))
             END IF
             IF ( min .GT. crit(ii_c,ii_phi_c,ii_eta_c) ) THEN
                min = crit(ii_c,ii_phi_c,ii_eta_c)
                kk_s     = kk
                kk_phi_s = kk_phi
                kk_eta_s = kk_eta
             END IF
             kk = kk + intervall
          END DO
          kk_phi = kk_phi + intervall
       END DO
       kk_eta = kk_eta + intervall
    END DO

    !-- After minimising the whole field --> recalc the opt_tensor depending on opt angles
    CALL transpose_mat (EE, kk_s, kk_phi_s, kk_eta_s, opt_tensor,BB)

    IF (PRESENT(a)) a = kk_s
    IF (PRESENT(p)) p = kk_phi_s
    IF (PRESENT(e)) e = kk_eta_s

  END SUBROUTINE opt_eff_stiff

  subroutine opt_eff_stiff_r812 (mode,EE,deg_a,stp_a,deg_p,stp_p,deg_e,stp_e,intervall,opt_tensor,a,p,e,outp)

    !-- Choose monotropic OR otrhotropic optimization 
    !-- Initial degrees must be integer!
    !-- Change mode to strings - way more convenient

    integer(kind=ik), intent(in)                               :: mode                   ! Switch between mono=1 and orth=2
    real(kind=rk), dimension(6,6), intent(in)                  :: EE                     ! Tensor to opimise
    integer(kind=rk), intent(in)                               :: deg_a, deg_e, deg_p    ! orientation to start with
    integer(kind=rk), intent(in)                               :: stp_a, stp_e, stp_p    ! Amount of steps to increase
    real(kind=rk), intent(in), optional                        :: intervall              ! stepwidth
    real(kind=rk), dimension(6,6), intent(out), optional       :: opt_tensor             ! optimised tensor - should be a mandatory one
    real(kind=rk), intent(out), optional                       :: a,p,e                  ! resulting angles
    logical, intent(in), optional                              :: outp                   ! check whether print result

    !----------------------------------------------------------------------------------------------

    real(kind=rk)                                              :: kk,   kk_phi,   kk_eta
    real(kind=rk)                                              :: kk_s, kk_phi_s, kk_eta_s    ! stored variables
    integer(kind=ik)                                           :: ii, ii_phi, ii_eta, sw1, sw2
    integer(kind=ik)                                           :: ii_c, ii_phi_c, ii_eta_c
    real(kind=rk)                                              :: alpha, phi, eta
    real(kind=rk)                                              :: min, intval
    real(kind=rk), dimension(3)                                :: n
    real(kind=rk), dimension(3,3)                              :: aa
    real(kind=rk), dimension(6,6)                              :: BB, tmp_r6x6
    real(kind=rk), dimension(:,:,:), allocatable               :: crit
    logical                                                    :: op=.FALSE.
    character(len=50)                                          :: opt_str
    Real(kind=8)   , Dimension( 8)                             :: tmp_r8 
    Real(kind=8)   , Dimension(12)                             :: tmp_r12
  
    allocate(crit(stp_a,stp_p,stp_e))

   
    if (present(intervall)) intval = intervall
    if (present(outp)) op = outp
    if ( mode .eq. 1_ik ) then
        opt_str="monotrop:  "
    elseif ( mode .eq. 2_ik ) then
        opt_str="orthotrop:"
     end if
     min=10000000_ik        !-- arbitrary initial value
     !----------------------------------------------------------------------------------------------
     !-- DEBUG -------------------------------------------------------------------------------------
     ! write(*,'(A,I5)')"deg_a=", deg_a 
     ! write(*,'(A,I5)')"deg_p=", deg_p 
     ! write(*,'(A,I5)')"deg_e=", deg_e 
     ! write(*,'(A,I5)')"stp_a=", stp_a 
     ! write(*,'(A,I5)')"stp_p=", stp_p 
     ! write(*,'(A,I5)')"stp_e=", stp_e 
     ! write(*,'(A,F5.3)')"intervall=", intervall
     !-- DEBUG -------------------------------------------------------------------------------------
     !----------------------------------------------------------------------------------------------

    kk_eta = deg_e
    Do ii_eta_c = 1, stp_e
       kk_phi = deg_p
        Do ii_phi_c = 1, stp_p
            kk = deg_a
                Do ii_c = 1, stp_a
                   call transpose_mat (EE, kk, kk_phi, kk_eta, tmp_r6x6,BB)
                   if ( mode .eq. 1_ik ) then
                    tmp_r8(1) = &
                         BB(6,1) * &
                         (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) &
                         +BB(5,1) * &
                         (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) &
                         +BB(4,1) * &
                         (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) &
                         +BB(3,1) * &
                         (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) &
                         +BB(2,1) * &
                         (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) &
                         +BB(1,1) * &
                         (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
                    tmp_r8(2) = &
                         BB(6,1) * &
                         (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) &
                         +BB(5,1) * &
                         (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) &
                         +BB(4,1) * &
                         (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) &
                         +BB(3,1) * &
                         (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) &
                         +BB(2,1) * &
                         (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) &
                         +BB(1,1) * &
                         (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
                    tmp_r8( 3) = &
                         BB(6,2) * &
                         (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) &
                         +BB(5,2) * &
                         (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) &
                         +BB(4,2) * &
                         (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) &
                         +BB(3,2) * &
                         (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) &
                         +BB(2,2) * &
                         (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) &
                         +BB(1,2) * &
                         (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
                    tmp_r8( 4) = &
                         BB(6,2) * &
                         (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) &
                         +BB(5,2) * &
                         (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) &
                         +BB(4,2) * &
                         (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) &
                         +BB(3,2) * &
                         (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) &
                         +BB(2,2) * &
                         (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) &
                         +BB(1,2) * &
                         (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
                    tmp_r8( 5) = &
                         BB(6,3) * &
                         (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) &
                         +BB(5,3) * &
                         (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) &
                         +BB(4,3) * &
                         (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) &
                         +BB(3,3) * &
                         (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) &
                         +BB(2,3) * &
                         (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) &
                         +BB(1,3) * &
                         (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
                    tmp_r8( 6) = &
                         BB(6,3) * &
                         (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) &
                         +BB(5,3) * &
                         (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) &
                         +BB(4,3) * &
                         (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) &
                         +BB(3,3) * &
                         (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) &
                         +BB(2,3) * &
                         (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) &
                         +BB(1,3) * &
                         (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
                    tmp_r8( 7) = &
                         BB(6,4) * &
                         (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) &
                         +BB(5,4) * &
                         (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) &
                         +BB(4,4) * &
                         (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) &
                         +BB(3,4) * &
                         (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) &
                         +BB(2,4) * &
                         (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) &
                         +BB(1,4) * &
                         (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
                    tmp_r8( 8) = &
                         BB(6,4) * &
                         (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) &
                         +BB(5,4) * &
                         (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) &
                         +BB(4,4) * &
                         (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) &
                         +BB(3,4) * &
                         (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) &
                         +BB(2,4) * &
                         (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) &
                         +BB(1,4) * &
                         (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))

                    crit(ii_c,ii_phi_c,ii_eta_c) = ( &
                         tmp_r8( 1)*tmp_r8( 1) + tmp_r8( 2)*tmp_r8( 2) + &
                         tmp_r8( 3)*tmp_r8( 3) + tmp_r8( 4)*tmp_r8( 4) + &
                         tmp_r8( 5)*tmp_r8( 5) + tmp_r8( 6)*tmp_r8( 6) + &
                         tmp_r8( 7)*tmp_r8( 7) + tmp_r8( 8)*tmp_r8( 8)   &
                         )
                 elseif ( mode .eq. 2_ik ) then
                    tmp_r12(1) = &
                         BB(6,1) * &
                         (BB(6,4)*EE(6,6)+BB(5,4)*EE(6,5)+BB(4,4)*EE(6,4)+BB(3,4)*EE(6,3)+BB(2,4)*EE(6,2)+BB(1,4)*EE(6,1)) + &
                         BB(5,1) * &
                         (EE(5,6)*BB(6,4)+BB(5,4)*EE(5,5)+BB(4,4)*EE(5,4)+BB(3,4)*EE(5,3)+BB(2,4)*EE(5,2)+BB(1,4)*EE(5,1)) + &
                         BB(4,1) * &
                         (EE(4,6)*BB(6,4)+EE(4,5)*BB(5,4)+BB(4,4)*EE(4,4)+BB(3,4)*EE(4,3)+BB(2,4)*EE(4,2)+BB(1,4)*EE(4,1)) + &
                         BB(3,1) * &
                         (EE(3,6)*BB(6,4)+EE(3,5)*BB(5,4)+EE(3,4)*BB(4,4)+EE(3,3)*BB(3,4)+BB(2,4)*EE(3,2)+BB(1,4)*EE(3,1)) + &
                         BB(2,1) * &
                         (EE(2,6)*BB(6,4)+EE(2,5)*BB(5,4)+EE(2,4)*BB(4,4)+EE(2,3)*BB(3,4)+EE(2,2)*BB(2,4)+BB(1,4)*EE(2,1)) + &
                         BB(1,1) * &
                         (EE(1,6)*BB(6,4)+EE(1,5)*BB(5,4)+EE(1,4)*BB(4,4)+EE(1,3)*BB(3,4)+EE(1,2)*BB(2,4)+EE(1,1)*BB(1,4))
                    tmp_r12(2) =  &
                         BB(6,1) * &
                         (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) + &
                         BB(5,1) * &
                         (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) + &
                         BB(4,1) * &
                         (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) + &
                         BB(3,1) * &
                         (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) + &
                         BB(2,1) * &
                         (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) + &
                         BB(1,1) * &
                         (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
                    tmp_r12(3) = &
                         BB(6,1) * &
                         (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) + &
                         BB(5,1) * &
                         (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) + &
                         BB(4,1) * &
                         (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) + &
                         BB(3,1) * &
                         (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) + &
                         BB(2,1) * &
                         (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) + &
                         BB(1,1) * &
                         (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
                    tmp_r12(4) =  &
                         BB(6,2) * &
                         (BB(6,4)*EE(6,6)+BB(5,4)*EE(6,5)+BB(4,4)*EE(6,4)+BB(3,4)*EE(6,3)+BB(2,4)*EE(6,2)+BB(1,4)*EE(6,1)) + &
                         BB(5,2) * &
                         (EE(5,6)*BB(6,4)+BB(5,4)*EE(5,5)+BB(4,4)*EE(5,4)+BB(3,4)*EE(5,3)+BB(2,4)*EE(5,2)+BB(1,4)*EE(5,1)) + &
                         BB(4,2) * &
                         (EE(4,6)*BB(6,4)+EE(4,5)*BB(5,4)+BB(4,4)*EE(4,4)+BB(3,4)*EE(4,3)+BB(2,4)*EE(4,2)+BB(1,4)*EE(4,1)) + &
                         BB(3,2) * &
                         (EE(3,6)*BB(6,4)+EE(3,5)*BB(5,4)+EE(3,4)*BB(4,4)+EE(3,3)*BB(3,4)+BB(2,4)*EE(3,2)+BB(1,4)*EE(3,1)) + &
                         BB(2,2) * &
                         (EE(2,6)*BB(6,4)+EE(2,5)*BB(5,4)+EE(2,4)*BB(4,4)+EE(2,3)*BB(3,4)+EE(2,2)*BB(2,4)+BB(1,4)*EE(2,1)) + &
                         BB(1,2) * &
                         (EE(1,6)*BB(6,4)+EE(1,5)*BB(5,4)+EE(1,4)*BB(4,4)+EE(1,3)*BB(3,4)+EE(1,2)*BB(2,4)+EE(1,1)*BB(1,4))
                    tmp_r12( 5) = &
                         BB(6,2) * &
                         (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) + &
                         BB(5,2) * &
                         (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) + &
                         BB(4,2) * &
                         (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) + &
                         BB(3,2) * &
                         (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) + &
                         BB(2,2) * &
                         (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) + &
                         BB(1,2) * &
                         (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
                    tmp_r12( 6) = &
                         BB(6,2) * &
                         (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) + &
                         BB(5,2) * &
                         (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) + &
                         BB(4,2) * &
                         (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) + &
                         BB(3,2) * &
                         (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) + &
                         BB(2,2) * &
                         (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) + &
                         BB(1,2) * &
                         (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
                    tmp_r12( 7) = &
                         BB(6,3) * &
                         (BB(6,4)*EE(6,6)+BB(5,4)*EE(6,5)+BB(4,4)*EE(6,4)+BB(3,4)*EE(6,3)+BB(2,4)*EE(6,2)+BB(1,4)*EE(6,1)) + &
                         BB(5,3) * &
                         (EE(5,6)*BB(6,4)+BB(5,4)*EE(5,5)+BB(4,4)*EE(5,4)+BB(3,4)*EE(5,3)+BB(2,4)*EE(5,2)+BB(1,4)*EE(5,1)) + &
                         BB(4,3) * &
                         (EE(4,6)*BB(6,4)+EE(4,5)*BB(5,4)+BB(4,4)*EE(4,4)+BB(3,4)*EE(4,3)+BB(2,4)*EE(4,2)+BB(1,4)*EE(4,1)) + &
                         BB(3,3) * &
                         (EE(3,6)*BB(6,4)+EE(3,5)*BB(5,4)+EE(3,4)*BB(4,4)+EE(3,3)*BB(3,4)+BB(2,4)*EE(3,2)+BB(1,4)*EE(3,1)) + &
                         BB(2,3) * &
                         (EE(2,6)*BB(6,4)+EE(2,5)*BB(5,4)+EE(2,4)*BB(4,4)+EE(2,3)*BB(3,4)+EE(2,2)*BB(2,4)+BB(1,4)*EE(2,1)) + &
                         BB(1,3) * &
                         (EE(1,6)*BB(6,4)+EE(1,5)*BB(5,4)+EE(1,4)*BB(4,4)+EE(1,3)*BB(3,4)+EE(1,2)*BB(2,4)+EE(1,1)*BB(1,4))
                    tmp_r12( 8) = &
                         BB(6,3) * &
                         (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) + &
                         BB(5,3) * &
                         (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) + &
                         BB(4,3) * &
                         (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) + &
                         BB(3,3) * &
                         (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) + &
                         BB(2,3) * &
                         (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) + &
                         BB(1,3) * &
                         (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
                    tmp_r12( 9) = &
                         BB(6,3) * &
                         (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) + &
                         BB(5,3) * &
                         (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) + &
                         BB(4,3) * &
                         (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) + &
                         BB(3,3) * &
                         (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) + &
                         BB(2,3) * &
                         (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) + &
                         BB(1,3) * &
                         (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
                    tmp_r12(10) = &
                         BB(6,4) * &
                         (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) + &
                         BB(5,4) * &
                         (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) + &
                         BB(4,4) * &
                         (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) + &
                         BB(3,4) * &
                         (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) + &
                         BB(2,4) * &
                         (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) + &
                         BB(1,4) * &
                         (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
                    tmp_r12(11) = &
                         BB(6,4) * &
                         (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) + &
                         BB(5,4) * &
                         (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) + &
                         BB(4,4) * &
                         (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) + &
                         BB(3,4) * &
                         (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) + &
                         BB(2,4) * &
                         (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) + &
                         BB(1,4) * &
                         (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
                    tmp_r12(12) = &
                         BB(6,5) * &
                         (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) + &
                         BB(5,5) * &
                         (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) + &
                         BB(4,5) * &
                         (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) + &
                         BB(3,5) * &
                         (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) + &
                         BB(2,5) * &
                         (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) + &
                         BB(1,5) * &
                         (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))

                    crit(ii_c,ii_phi_c,ii_eta_c) = (&
                         tmp_r12( 1)*tmp_r12( 1) + tmp_r12( 2)*tmp_r12( 2) + tmp_r12( 3)*tmp_r12( 3) + &
                         tmp_r12( 4)*tmp_r12( 4) + tmp_r12( 5)*tmp_r12( 5) + tmp_r12( 6)*tmp_r12( 6) + &
                         tmp_r12( 7)*tmp_r12( 7) + tmp_r12( 8)*tmp_r12( 8) + tmp_r12( 9)*tmp_r12( 9) + &
                         tmp_r12(10)*tmp_r12(10) + tmp_r12(11)*tmp_r12(11) + &
                         tmp_r12(12)*tmp_r12(12) &
                         )
                   end if
                   if ( min .gt. crit(ii_c,ii_phi_c,ii_eta_c) ) then
                        min = crit(ii_c,ii_phi_c,ii_eta_c)
                        kk_s     = kk
                        kk_phi_s = kk_phi
                        kk_eta_s = kk_eta
                   end if
                   kk = kk + intervall
            End Do
            kk_phi = kk_phi + intervall
        end Do
        kk_eta = kk_eta + intervall
     end Do

    !-- After minimising the whole field --> recalc the opt_tensor depending on opt angles
    call transpose_mat (EE, kk_s, kk_phi_s, kk_eta_s, opt_tensor,BB)




    if ( outp .eqv. .TRUE. ) then
       write(*,'(A,F6.1,A,F10.1,A,F8.1,A,A,F10.1)',advance='YES')"alpha: ",kk_s,"      phi: ",kk_phi_s,&
            & "         eta: ",kk_eta_s, "        Min ",trim(opt_str),min
    end if

    if (present(a)) a = kk_s
    if (present(p)) p = kk_phi_s
    if (present(e)) e = kk_eta_s

end subroutine opt_eff_stiff_r812


SUBROUTINE etimea (sec)

    INTEGER(KIND=ik), INTENT(in)                       :: sec
    INTEGER(KIND=ik)                                   :: mins, hours, secs, seconds, mins_s, remainder

    mins_s    = MODULO(sec,  3600_ik)
    seconds   = MODULO(mins_s, 60_ik)
    remainder = MODULO(seconds, 1_ik)

    hours = (sec-mins_s)     / 3600_ik
    mins  = (mins_s-seconds) / 60_ik
    secs  = (seconds-remainder)

    write(*,'(A,I3,A,I2,A,I2,A)',advance='NO')"ETA: ",&
         hours,":",mins,":",secs,"     hhh:mm:ss"

END SUBROUTINE etimea


subroutine transpose_mat (EE, kk, kk_phi, kk_eta, tmp_r6x6,BB)

    real(kind=rk), intent(in)                                  :: kk,   kk_phi,   kk_eta
    real(kind=rk), dimension(6,6), intent(in)                  :: EE                     ! Tensor to opimise
    real(kind=rk), dimension(6,6), intent(out)                 :: tmp_r6x6
    real(kind=rk), dimension(6,6), intent(out),optional        :: BB
    real(kind=rk)                                              :: alpha, phi, eta
    real(kind=rk), dimension(3)                                :: n
    real(kind=rk), dimension(3,3)                              :: aa

    !-- Degrees as input, radian as output to sin/cos
    alpha = Real(kk)     * pi / 180._rk
    phi   = Real(kk_phi) * pi / 180._rk
    eta   = Real(kk_eta) * pi / 180._rk

    n = [cos(phi)*sin(eta) , sin(phi)*sin(eta) , cos(eta) ]
    n = n / sqrt(sum(n*n))

    aa = rot_alg(n,alpha)

    BB = tra_R6(aa)

    tmp_r6x6 = matmul(matmul(transpose(BB),EE),BB)
end subroutine transpose_mat

Subroutine mono_crit (EE,res)
    real(kind=rk), dimension(6,6), intent(in)         :: EE
    real(kind=rk), intent (out)                       :: res

    !-- Square based method
    res = sum(EE(1:4,5:6)*EE(1:4,5:6))
end subroutine mono_crit

subroutine ortho_crit (EE,res)
    real(kind=rk), dimension(6,6), intent(in)         :: EE
    real(kind=rk), intent (out)                       :: res

    !-- Square based method
    res = ( sum(EE(1:3,4:6) * EE(1:3,4:6)) + &
            sum(EE( 4 ,5:6) * EE( 4 ,5:6)) + &
                EE( 5 , 6 ) * EE( 5 , 6 ))
end subroutine ortho_crit

SUBROUTINE checksym (TT,sym)

  INTEGER(kind=ik)                                      :: i,j
  REAL   (kind=rk), DIMENSION(6,6), INTENT(in)          :: TT     ! Tensor
  REAL   (kind=rk), DIMENSION(6,6), INTENT(out)         :: sym    ! Tensor
  sym=TT

  sym(1  ,2)=   ABS(TT(1,2))-ABS(TT(2,1))
  sym(1:2,3)=(/ ABS(TT(1,3))-ABS(TT(3,1)), ABS(TT(2,3))-ABS(TT(3,2))  /)
  sym(1:3,4)=(/ ABS(TT(1,4))-ABS(TT(4,1)), ABS(TT(2,4))-ABS(TT(4,2)), ABS(TT(3,4))-ABS(TT(4,3))  /)
  sym(1:4,5)=(/ ABS(TT(1,5))-ABS(TT(5,1)), ABS(TT(2,5))-ABS(TT(5,2)), ABS(TT(3,5))-ABS(TT(5,3)), &
                ABS(TT(4,5))-ABS(TT(5,4))  /)
  sym(1:5,6)=(/ ABS(TT(1,6))-ABS(TT(6,1)), ABS(TT(2,6))-ABS(TT(6,2)), ABS(TT(3,6))-ABS(TT(6,3)), &
                ABS(TT(4,6))-ABS(TT(6,4)), ABS(TT(5,6))-ABS(TT(6,5)) /)

  !-- Arbitrary threshold!!
  !-- Needs some specifications
  DO i=1,5
     DO j=1,i
        IF ( ABS(sym(j,i+1)) .LT. 3._rk) THEN
           sym(j,i+1) = 0._rk
        END IF
     END DO
  END DO
! sym=transpose(sym)

END SUBROUTINE checksym



FUNCTION DoAO(EE) RESULT(OUT)
  !-- Degree of Orthotropic anisotropy
  REAL   (kind=rk), DIMENSION(6,6), INTENT(in) :: EE
  REAL   (kind=rk)                             :: OUT
  !-- Degree of Anisotropy - average 0 divided by average non-zero entry
  OUT = ((     ABS(EE(1,4)) +    ABS(EE(1,5)) +    ABS(EE(1,6)) &
          +    ABS(EE(2,4)) +    ABS(EE(2,5)) +    ABS(EE(2,6)) &
          +    ABS(EE(3,4)) +    ABS(EE(3,5)) +    ABS(EE(3,6)) &
                            +    ABS(EE(4,5)) +    ABS(EE(4,6)) &
                                              +    ABS(EE(5,6)) )/12.0_rk ) / &
        ((     ABS(EE(1,1)) +  2*ABS(EE(1,2)) +  2*ABS(EE(1,3)) &
                            +    ABS(EE(2,2)) +  2*ABS(EE(2,3)) &
                                              +    ABS(EE(3,3)) &
          +    ABS(EE(4,4)) +    ABS(EE(5,5)) +    ABS(EE(6,6)) )/12.0_rk)*100.0_rk
END FUNCTION DoAO

FUNCTION DoAM(EE) RESULT(OUT)
    !-- Degree of Orthotropic anisotropy
    REAL   (kind=rk), DIMENSION(6,6), INTENT(in) :: EE
    REAL   (kind=rk)                             :: OUT
    !-- Degree of Anisotropy - average 0 divided by average non-zero entry
    OUT = (( +    ABS(EE(1,5)) +    ABS(EE(1,6)) &
             +    ABS(EE(2,5)) +    ABS(EE(2,6)) &
             +    ABS(EE(3,5)) +    ABS(EE(3,6)) &
             +    ABS(EE(4,5)) +    ABS(EE(4,6)) )/8.0_rk ) / &
         ((       ABS(EE(1,1)) + 2* ABS(EE(1,2)) + 2* ABS(EE(1,3)) + 2* ABS(EE(1,4)) &
                               +    ABS(EE(2,2)) + 2* ABS(EE(2,3)) + 2* ABS(EE(2,4)) &
                                                 +    ABS(EE(3,3)) + 2* ABS(EE(3,4)) &
                                                                   +    ABS(EE(4,4)) &
                                                                   +    ABS(EE(5,5)) + 2* ABS(EE(5,6)) &
                                                                   +    ABS(EE(6,6)))/20.0_rk)*100.0_rk
END FUNCTION DoAM

SUBROUTINE tilt_tensor (TENSOR,T_OUT)

  !-- Tilts a tensor by 90° around axis x, y or z - respectively 1, 2 or 3
  !-- other angles are not implemented
  !-- intended to change S11, S22 and S33 to get an S11 .gt. S22 .gt. S33 ordering

  INTEGER(kind=ik)                                      :: axis
  REAL   (kind=rk), DIMENSION(6,6), INTENT(in)          :: TENSOR
  REAL   (kind=rk), DIMENSION(6,6), INTENT(out)         :: T_OUT
  REAL   (kind=rk), DIMENSION(6,6)                      :: BB1, BB2, BB3, tmp6x6
  REAL   (kind=rk), DIMENSION(3)                        :: n1, n2, n3
  REAL   (kind=rk), DIMENSION(3,3)                      :: aa1, aa2, aa3
  REAL   (kind=rk), PARAMETER                           :: alpha = 1.570796326794897 ! 90°*pi/180

  n1 = [1,0,0]
  n2 = [0,1,0]
  n3 = [0,0,1]
  aa1 = rot_alg(n1,alpha)
  aa2 = rot_alg(n2,alpha)
  aa3 = rot_alg(n3,alpha)
  BB1 = tra_R6 (aa1)
  BB2 = tra_R6 (aa2)
  BB3 = tra_R6 (aa3)

  IF ( (TENSOR(1,1) .GT. TENSOR(2,2)) .AND.  &
  (TENSOR(1,1) .GT. TENSOR(3,3)) .AND.  (TENSOR(2,2) .GT. TENSOR(3,3)) ) THEN
          !-- 123
          T_OUT = TENSOR
  Else IF ( (TENSOR(1,1) .GT. TENSOR(2,2)) .AND.  &
  (TENSOR(1,1) .GT. TENSOR(3,3)) .AND.  (TENSOR(2,2) .LT. TENSOR(3,3)) ) THEN
          !-- 132 => 123
          T_OUT  = MATMUL(MATMUL(TRANSPOSE(BB1),TENSOR),BB1)
  Else IF ( (TENSOR(1,1) .LT. TENSOR(2,2)) .AND.  &
  (TENSOR(1,1) .LT. TENSOR(3,3)) .AND.  (TENSOR(2,2) .GT. TENSOR(3,3)) ) THEN
          !-- 231 => 132
          tmp6x6 = MATMUL(MATMUL(TRANSPOSE(BB3),TENSOR),BB3)
          !-- 132 => 123
          T_OUT  = MATMUL(MATMUL(TRANSPOSE(BB1),tmp6x6),BB1)
  Else IF ( (TENSOR(1,1) .LT. TENSOR(2,2)) .AND.  &
  (TENSOR(1,1) .GT. TENSOR(3,3)) .AND.  (TENSOR(2,2) .GT. TENSOR(3,3)) ) THEN
          !-- 213 => 123
          T_OUT  = MATMUL(MATMUL(TRANSPOSE(BB3),TENSOR),BB3)
  Else IF ( (TENSOR(1,1) .GT. TENSOR(2,2)) .AND.  &
  (TENSOR(1,1) .LT. TENSOR(3,3)) .AND.  (TENSOR(2,2) .LT. TENSOR(3,3)) ) THEN
          !-- 312 => 132
          tmp6x6 = MATMUL(MATMUL(TRANSPOSE(BB2),TENSOR),BB2)
          !-- 132 => 123
          T_OUT  = MATMUL(MATMUL(TRANSPOSE(BB1),tmp6x6),BB1)
  Else IF ( (TENSOR(1,1) .LT. TENSOR(2,2)) .AND.  &
  (TENSOR(1,1) .LT. TENSOR(3,3)) .AND.  (TENSOR(2,2) .LT. TENSOR(3,3)) ) THEN
          !-- 321 => 123
          T_OUT  = MATMUL(MATMUL(TRANSPOSE(BB2),TENSOR),BB2)
  END IF

END SUBROUTINE tilt_tensor

END MODULE opt_stiffness
