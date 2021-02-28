!----------------------------------------------------------------------------------------------
!--
!-- Module for optimising R6x6 anisotropic tensors
!--
!-- author:                Johannes Gebert        (HLRS, NUM)
!-- numeric kernel by:     Ralf Schneider         (HLRS, NUM)
!--
!-- mod  27.06.2020        Previous programs ignored
!--
!----------------------------------------------------------------------------------------------

module opt_stiffness

use kinds
use mat_matrices

Implicit None

    contains

    subroutine opt_eff_stiff (mode,EE,deg_a,stp_a,deg_p,stp_p,deg_e,stp_e,intervall,opt_tensor,a,p,e,outp)

    !-- Choose monotropic OR otrhotropic optimization 
    !-- Initial degrees must be integer!
    !-- Change mode to strings - way more convenient

    integer(kind=ik), intent(in)                               :: mode                        ! Switch between mono=1 and orth=2
    real(kind=rk), dimension(6,6), intent(in)                  :: EE                          ! Tensor to opimise
    integer(kind=rk), intent(in)                               :: deg_a, deg_e, deg_p         ! orientation to start with
    integer(kind=rk), intent(in)                               :: stp_a, stp_e, stp_p         ! Amount of steps to increase
    real(kind=rk), intent(in), optional                        :: intervall                   ! stepwidth
    real(kind=rk), dimension(6,6), intent(out), optional       :: opt_tensor                  ! optimised tensor - should be a mandatory one
    real(kind=rk), intent(out), optional                       :: a,p,e                       ! resulting angles
    logical, intent(in), optional                              :: outp                        ! check whether print result

    !----------------------------------------------------------------------------------------------

    real(kind=ik)                                              :: kk,   kk_phi,   kk_eta
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

    allocate(crit(stp_a,stp_p,stp_e))


    if (present(intervall)) intval = intervall
    if (present(outp)) op = outp
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
                      call mono_crit (tmp_r6x6,res=crit(ii_c,ii_phi_c,ii_eta_c))
                   elseif ( mode .eq. 2_ik ) then
                      call ortho_crit (tmp_r6x6,crit(ii_c,ii_phi_c,ii_eta_c))
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

    if (present(a)) a = kk_s
    if (present(p)) p = kk_phi_s
    if (present(e)) e = kk_eta_s

  end subroutine opt_eff_stiff

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

    real(kind=ik)                                              :: kk,   kk_phi,   kk_eta
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


subroutine etimea (sec)

    integer(kind=ik), intent(in)                       :: sec
    integer(kind=ik)                                   :: mins,    hours,  secs, seconds, mins_s, rest

    mins_s=modulo(sec, 3600)
    seconds=modulo(mins_s, 60)
    rest=modulo(seconds, 1)

    hours=(sec-mins_s)/3600
    mins=(mins_s-seconds)/60
    secs=(seconds-rest)

    write(*,'(A,I3,A,I2,A,I2,A)',advance='NO')"ETA: ",&
         hours,":",mins,":",secs,"     hhh:mm:ss"

end subroutine etimea


subroutine transpose_mat (EE, kk, kk_phi, kk_eta, tmp_r6x6,BB)

    real(kind=ik), intent(in)                                  :: kk,   kk_phi,   kk_eta
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

subroutine checksym (TT,sym)

  integer(kind=ik)                                   :: i,j
  real(kind=rk), dimension(6,6), intent(in)          :: TT     ! Tensor
  real(kind=rk), dimension(6,6), intent(out)         :: sym    ! Tensor
  sym=TT

  sym(1  ,2)=   abs(TT(1,2))-abs(TT(2,1))
  sym(1:2,3)=(/ abs(TT(1,3))-abs(TT(3,1)), abs(TT(2,3))-abs(TT(3,2))  /)
  sym(1:3,4)=(/ abs(TT(1,4))-abs(TT(4,1)), abs(TT(2,4))-abs(TT(4,2)), abs(TT(3,4))-abs(TT(4,3))  /)
  sym(1:4,5)=(/ abs(TT(1,5))-abs(TT(5,1)), abs(TT(2,5))-abs(TT(5,2)), abs(TT(3,5))-abs(TT(5,3)), &
       & abs(TT(4,5))-abs(TT(5,4))  /)
  sym(1:5,6)=(/ abs(TT(1,6))-abs(TT(6,1)), abs(TT(2,6))-abs(TT(6,2)), abs(TT(3,6))-abs(TT(6,3)), &
       & abs(TT(4,6))-abs(TT(6,4)), abs(TT(5,6))-abs(TT(6,5)) /)

  !-- Arbitrary threshold!!
  !-- Needs some specifications
  do i=1,5
     do j=1,i
        if ( abs(sym(j,i+1)) .lt. 3._rk) then
           sym(j,i+1) = 0._rk
        endif
     end do
  end do
! sym=transpose(sym)

end subroutine checksym



Function DoAO(EE) Result(OUT)
  !-- Degree of Orthotropic anisotropy
  Real(kind=rk), Dimension(6,6), intent(in) :: EE
  Real(kind=rk)                             :: OUT
  !-- Degree of Anisotropy - average 0 divided by average non-zero entry
  OUT = ((     abs(EE(1,4)) +    abs(EE(1,5)) +    abs(EE(1,6)) &
          +    abs(EE(2,4)) +    abs(EE(2,5)) +    abs(EE(2,6)) &
          +    abs(EE(3,4)) +    abs(EE(3,5)) +    abs(EE(3,6)) &
                            +    abs(EE(4,5)) +    abs(EE(4,6)) &
                                              +    abs(EE(5,6)) )/12.0_rk ) / &
        ((     abs(EE(1,1)) +  2*abs(EE(1,2)) +  2*abs(EE(1,3)) &
                            +    abs(EE(2,2)) +  2*abs(EE(2,3)) &
                                              +    abs(EE(3,3)) &
          +    abs(EE(4,4)) +    abs(EE(5,5)) +    abs(EE(6,6)) )/12.0_rk)*100.0_rk
End Function DoAO

Function DoAM(EE) Result(OUT)
    !-- Degree of Orthotropic anisotropy
    Real(kind=rk), Dimension(6,6), intent(in) :: EE
    Real(kind=rk)                             :: OUT
    !-- Degree of Anisotropy - average 0 divided by average non-zero entry
    OUT = (( +    abs(EE(1,5)) +    abs(EE(1,6)) &
             +    abs(EE(2,5)) +    abs(EE(2,6)) &
             +    abs(EE(3,5)) +    abs(EE(3,6)) &
             +    abs(EE(4,5)) +    abs(EE(4,6)) )/8.0_rk ) / &
         ((       abs(EE(1,1)) + 2* abs(EE(1,2)) + 2* abs(EE(1,3)) + 2* abs(EE(1,4)) &
                               +    abs(EE(2,2)) + 2* abs(EE(2,3)) + 2* abs(EE(2,4)) &
                                                 +    abs(EE(3,3)) + 2* abs(EE(3,4)) &
                                                                   +    abs(EE(4,4)) &
                                                                   +    abs(EE(5,5)) + 2* abs(EE(5,6)) &
                                                                   +    abs(EE(6,6)))/20.0_rk)*100.0_rk
End Function DoAM

subroutine tilt_tensor (TENSOR,T_OUT)

  !-- Tilts a tensor by 90° around axis x, y or z - respectively 1, 2 or 3
  !-- other angles are not implemented
  !-- intended to change S11, S22 and S33 to get an S11 .gt. S22 .gt. S33 ordering

  integer(kind=ik)                                   :: axis
  real(kind=rk), dimension(6,6), intent(in)          :: TENSOR
  real(kind=rk), dimension(6,6), intent(out)         :: T_OUT
  real(kind=rk), dimension(6,6)                      :: BB1, BB2, BB3, tmp6x6
  real(kind=rk), dimension(3)                        :: n1, n2, n3
  real(kind=rk), dimension(3,3)                      :: aa1, aa2, aa3
  real(kind=rk), parameter                           :: alpha = 1.570796326794897 ! 90°*pi/180

  n1 = [1,0,0]
  n2 = [0,1,0]
  n3 = [0,0,1]
  aa1 = rot_alg(n1,alpha)
  aa2 = rot_alg(n2,alpha)
  aa3 = rot_alg(n3,alpha)
  BB1 = tra_R6(aa1)
  BB2 = tra_R6(aa2)
  BB3 = tra_R6(aa3)

  If ( (TENSOR(1,1) .gt. TENSOR(2,2)) .AND.  &
  (TENSOR(1,1) .gt. TENSOR(3,3)) .AND.  (TENSOR(2,2) .gt. TENSOR(3,3)) ) then
          !-- 123
          T_OUT = TENSOR
  Else If ( (TENSOR(1,1) .gt. TENSOR(2,2)) .AND.  &
  (TENSOR(1,1) .gt. TENSOR(3,3)) .AND.  (TENSOR(2,2) .lt. TENSOR(3,3)) ) then
          !-- 132 => 123
          T_OUT  = matmul(matmul(transpose(BB1),TENSOR),BB1)
  Else If ( (TENSOR(1,1) .lt. TENSOR(2,2)) .AND.  &
  (TENSOR(1,1) .lt. TENSOR(3,3)) .AND.  (TENSOR(2,2) .gt. TENSOR(3,3)) ) then
          !-- 231 => 132
          tmp6x6 = matmul(matmul(transpose(BB3),TENSOR),BB3)
          !-- 132 => 123
          T_OUT  = matmul(matmul(transpose(BB1),tmp6x6),BB1)
  Else If ( (TENSOR(1,1) .lt. TENSOR(2,2)) .AND.  &
  (TENSOR(1,1) .gt. TENSOR(3,3)) .AND.  (TENSOR(2,2) .gt. TENSOR(3,3)) ) then
          !-- 213 => 123
          T_OUT  = matmul(matmul(transpose(BB3),TENSOR),BB3)
  Else If ( (TENSOR(1,1) .gt. TENSOR(2,2)) .AND.  &
  (TENSOR(1,1) .lt. TENSOR(3,3)) .AND.  (TENSOR(2,2) .lt. TENSOR(3,3)) ) then
          !-- 312 => 132
          tmp6x6 = matmul(matmul(transpose(BB2),TENSOR),BB2)
          !-- 132 => 123
          T_OUT  = matmul(matmul(transpose(BB1),tmp6x6),BB1)
  Else If ( (TENSOR(1,1) .lt. TENSOR(2,2)) .AND.  &
  (TENSOR(1,1) .lt. TENSOR(3,3)) .AND.  (TENSOR(2,2) .lt. TENSOR(3,3)) ) then
          !-- 321 => 123
          T_OUT  = matmul(matmul(transpose(BB2),TENSOR),BB2)
  End If

end subroutine tilt_tensor



end module opt_stiffness
