!==============================================================================
!> \file mod_math.f90
!> Module with mathematical subroutines
!>
!> \author Ralf Schneider
!> \date 11.06.2012

Module math

  USE global_std
  use chain_routines
  USE auxiliaries

  Implicit None

  Real(kind=rk), Parameter :: num_zero   = 1.E-9_rk
  Real(kind=rk), Parameter :: sq2        = sqrt(2._rk)
  Real(kind=rk), Parameter :: pi         = acos(-1._rk)
  Real(kind=rk), Parameter :: inv180     = 1._rk/180._rk
  Real(kind=rk), Parameter :: pi_div_180 = acos(-1._rk)/180._rk
  

  !-- Higher dimensional numbers
  TYPE Quaternion
     REAL (KIND=rk)            :: w,x,y,z
  END TYPE Quaternion

  Logical, Parameter       :: mmdbg=.false.

contains

  Function rot_x(alpha) Result(aa)

    Real(kind=rk),                 intent(in) :: alpha
    Real(kind=rk), Dimension(3,3)             :: aa
    
    aa(1,:) = (/ 1._rk ,   0._rk    ,   0._rk     /)
    aa(2,:) = (/ 0._rk , cos(alpha) , -sin(alpha) /)
    aa(3,:) = (/ 0._rk , sin(alpha) ,  cos(alpha) /)

  End Function rot_x

  Function rot_y(alpha) Result(aa)

    Real(kind=rk),                 intent(in) :: alpha
    Real(kind=rk), Dimension(3,3)             :: aa
    
    aa(1,:) = (/ cos(alpha) ,   0._rk    ,  sin(alpha) /)
    aa(2,:) = (/   0._rk    ,   1._rk    ,   0._rk     /)
    aa(3,:) = (/-sin(alpha) ,   0._rk    ,  cos(alpha) /)

  End Function rot_y

  Function rot_z(alpha) Result(aa)

    Real(kind=rk),                 intent(in) :: alpha
    Real(kind=rk), Dimension(3,3)             :: aa
    
    aa(1,:) = (/ cos(alpha) , -sin(alpha) ,   0._rk     /)
    aa(2,:) = (/ sin(alpha) ,  cos(alpha) ,   0._rk     /)
    aa(3,:) = (/   0._rk    ,   0._rk     ,   1._rk     /)

  End Function rot_z

  !****************************************************************************
  !*                                                                         **
  !* Subroutine to perform a rotation on a list of coordinates               **
  !*                                                                         **
  !* last edited : on  05.03.2009                                            **
  !*               by  Ralf Schneider                                        **
  !*                                                                         **  
  Function rot_alg(axis, angle) Result(rr)

    real(kind=rk), dimension(3)           , Intent(In) :: axis
    !real(kind=rk), dimension(3)           , Intent(InOut) :: axis

    real(kind=rk)                         , Intent(In)    :: angle

    real(kind=rk), dimension(3, 3)                        :: rr

    !**************************************************************************

    !** normalize rotation axis ***********************************************
    !axis = axis / sqrt(sum(axis*axis))

    !** Setup transformation matrix *******************************************
    rr(1,1) = cos(angle) + axis(1)*axis(1)* (1._8 - cos(angle))
    rr(1,2) = axis(1)*axis(2)* (1._8 - cos(angle)) - axis(3) * sin(angle)
    rr(1,3) = axis(1)*axis(3)* (1._8 - cos(angle)) + axis(2) * sin(angle)

    rr(2,1) = axis(2)*axis(1)* (1._8 - cos(angle)) + axis(3) * sin(angle)
    rr(2,2) = cos(angle) + axis(2)*axis(2)* (1._8 - cos(angle))
    rr(2,3) = axis(2)*axis(3)* (1._8 - cos(angle)) - axis(1) * sin(angle)

    rr(3,1) = axis(3)*axis(1)* (1._8 - cos(angle)) - axis(2) * sin(angle)
    rr(3,2) = axis(3)*axis(2)* (1._8 - cos(angle)) + axis(1) * sin(angle)
    rr(3,3) = cos(angle) + axis(3)*axis(3)* (1._8 - cos(angle))

  End Function rot_alg

  Function tra_R6(aa) Result(BB)

    Real(kind=rk), Dimension(3,3), intent(in) :: aa
    Real(kind=rk), Dimension(6,6)             :: BB

    Real(kind=rk), Parameter :: sq2 = sqrt(2._rk)

    BB(1,:) = (/ aa(1,1)**2 , aa(1,2)**2 , aa(1,3)**2 , sq2*aa(1,1)*aa(1,2) , sq2*aa(1,1)*aa(1,3), sq2*aa(1,2)*aa(1,3)  /)
    BB(2,:) = (/ aa(2,1)**2 , aa(2,2)**2 , aa(2,3)**2 , sq2*aa(2,1)*aa(2,2) , sq2*aa(2,1)*aa(2,3), sq2*aa(2,2)*aa(2,3)  /)
    BB(3,:) = (/ aa(3,1)**2 , aa(3,2)**2 , aa(3,3)**2 , sq2*aa(3,1)*aa(3,2) , sq2*aa(3,1)*aa(3,3), sq2*aa(3,2)*aa(3,3)  /)
    
    BB(4,:) = (/ sq2*aa(2,1)*aa(1,1) , sq2*aa(2,2)*aa(1,2) , sq2*aa(2,3)*aa(1,3) , &
         aa(2,1)*aa(1,2)+aa(2,2)*aa(1,1) , aa(2,1)*aa(1,3)+aa(2,3)*aa(1,1) , aa(2,2)*aa(1,3)+aa(2,3)*aa(1,2) /)
    BB(5,:) = (/ sq2*aa(1,1)*aa(3,1) , sq2*aa(1,2)*aa(3,2) , sq2*aa(1,3)*aa(3,3) ,  &
         aa(1,1)*aa(3,2)+aa(1,2)*aa(3,1) , aa(1,1)*aa(3,3)+aa(1,3)*aa(3,1) , aa(1,2)*aa(3,3)+aa(1,3)*aa(3,2) /)
    BB(6,:) = (/ sq2*aa(2,1)*aa(3,1) , sq2*aa(2,2)*aa(3,2) , sq2*aa(2,3)*aa(3,3) ,  &
         aa(2,1)*aa(3,2)+aa(2,2)*aa(3,1) , aa(2,1)*aa(3,3)+aa(2,3)*aa(3,1) , aa(2,2)*aa(3,3)+aa(2,3)*aa(3,2) /)

  End Function tra_R6

  !============================================================================
  !> Subroutine which optimizes a 6x6 marix according to an monotropic
  !> criteria
  !>
  !> -- CR1 Monotropic ----------------------------------------
  !>  Sum over Aij*Aij with i in [1,4]**2 and j in [5,6]**2
  !>  +-                       -+
  !>  | a11 a12 a13 a14 A15 A16 |
  !>  | a12 a22 a23 a24 A25 A26 |
  !>  | a13 a12 a33 a34 A35 A36 |
  !>  | a14 a12 a13 a44 A45 A46 |
  !>  | a15 a12 a13 a14 a55 a56 |
  !>  | a16 a12 a13 a14 a15 a66 |
  !>  +-                       -+
  !>
  subroutine opt_mono(EE_orig, EE, fdir, fcrit)

    Real(kind=rk)    , Dimension(6,6), intent(in)      :: EE_orig
    Real(kind=rk)    , Dimension(6,6), intent(out)     :: EE     
    real(kind=rk)                    , intent(out)     :: fcrit
    Real(kind=rk)    , Dimension(3,3), intent(out)     :: fdir

    Integer(kind=ik)                                   :: kk_eta, kk_phi, kk
    Integer(kind=ik)                                   :: ii_eta, ii_phi, ii, jj

    Real(kind=rk)                                      :: alpha, phi, eta

    Real(kind=rk)    , Dimension(3)                    :: n
    Real(kind=rk)    , Dimension(3,3)                  :: aa
    Real(kind=rk)    , Dimension(6,6)                  :: BB, tmp_r6x6

    Integer          , Dimension(:,:,:,:), ALLOCATABLE :: ang
    Real(kind=rk)    , Dimension(:,:,:)  , ALLOCATABLE :: crit
    Real(kind=rk)    , Dimension(0 :16)                :: crit_min
    Integer          , Dimension(3)                    :: s_loop,e_loop, mlc

    Real(kind=rk), Parameter :: num_zero = 1.E-9_rk
    Real(kind=rk), Parameter :: pi  = acos(-1._rk)
    !**************************************************************************

    ALLOCATE(ang(3,0:180,0:180,0:90))
    ALLOCATE(crit( 0:180,0:180,0:90))
    
    ee = ee_orig

    kk_eta = 0
    kk_phi = 0
    kk = 0

    !** Initial parameter space sweep *****************************************
    Do ii_eta = 0 , 90 , 1

       kk_phi = 0

       Do ii_phi = 0 , 180 , 1

          kk = 0

          Do ii = 0 , 180 , 1

             alpha = Real(ii,rk)     * pi / 180._rk
             phi   = Real(ii_phi,rk) * pi / 180._rk
             eta   = Real(ii_eta,rk) * pi / 180._rk

             n = [cos(phi)*sin(eta) , sin(phi)*sin(eta) , cos(eta) ]
             n = n / sqrt(sum(n*n))

             aa = rot_alg(n,alpha)
             BB = tra_R6(aa)
             tmp_r6x6 = matmul(matmul(transpose(BB),EE),BB)

             ang(:,kk,kk_phi,kk_eta)  = [ii,ii_phi,ii_eta]

             !-- CR1 Monotropic ----------------------------------------
             crit(kk,kk_phi,kk_eta) = (&
                   sum((tmp_r6x6(1:4,5:6))*(tmp_r6x6(1:4,5:6))))

             kk = kk + 1
          End Do
          kk_phi = kk_phi + 1
       end Do
       kk_eta = kk_eta + 1
    end Do

    !==========================================================================
    !== Iteration of Crit =====================================================

    crit_min = 0._rk
    crit_min(0) = minval(crit)

    mlc = minloc(crit)-1
    If (mmdbg) then
       write(un_lf,FMT_MSG_A3I0)'Initial Minloc  CR_1 : ',mlc
       write(un_lf,FMT_MSG_AF0) 'Initial Minimum CR_1 : ',crit_min(0)
    End If

    jj = 1

    Do 

       mlc = minloc(crit(0:kk-1,0:kk_phi-1,0:kk_eta-1))-1

       s_loop = (ang(:,mlc(1),mlc(2),mlc(3))-1)*10
       e_loop = (ang(:,mlc(1),mlc(2),mlc(3))+1)*10

       kk_eta = 0
       kk_phi = 0
       kk = 0

       If (mmdbg) then
          write(un_lf,FMT_MSG_AI0) 'Iteration            : ',jj
          
          write(un_lf,FMT_MSG_A3I0)'Loop start           : ',s_loop
          write(un_lf,FMT_MSG_A3I0)'Loop end             : ',e_loop
       End If

       Do ii_eta = s_loop(3), e_loop(3)

          kk_phi = 0

          Do ii_phi = s_loop(2), e_loop(2)

             kk = 0

             Do ii = s_loop(1), e_loop(1)

                alpha = Real(ii,rk)     * pi / (180._rk*(10._rk**jj))
                phi   = Real(ii_phi,rk) * pi / (180._rk*(10._rk**jj))
                eta   = Real(ii_eta,rk) * pi / (180._rk*(10._rk**jj))

                n = [cos(phi)*sin(eta) , sin(phi)*sin(eta) , cos(eta) ] 
                n = n / sqrt(sum(n*n))

                aa = rot_alg(n,alpha)
                BB = tra_R6(aa)
                tmp_r6x6 = matmul(matmul(transpose(BB),EE),BB)

                ang(:,kk,kk_phi,kk_eta)  = [ii,ii_phi,ii_eta]

                !-- CR1 Monotropic ---------------------------------------
                crit(kk,kk_phi,kk_eta) = (&
                   sum((tmp_r6x6(1:4,5:6))*(tmp_r6x6(1:4,5:6)))) 
                
                kk = kk + 1

             End Do
             kk_phi = kk_phi + 1
          end Do
          kk_eta = kk_eta + 1
       end Do

       crit_min(jj) = minval(crit(0:kk-1,0:kk_phi-1,0:kk_eta-1))

       If (mmdbg) then
          write(un_lf,FMT_MSG_AF0)'Minimum CR_1         : ',crit_min(jj)
       End If

       If ( (abs(crit_min(jj-1)-crit_min(jj)) < num_zero) .OR. (jj > 16)) Exit

       jj = jj + 1

    End Do

    mlc = minloc(crit(0:kk-1,0:kk_phi-1,0:kk_eta-1))-1

    alpha = Real( ang(1,mlc(1),mlc(2),mlc(3)),rk ) * pi / & 
            (180._rk*(10._rk**jj-1))
    phi   = Real( ang(2,mlc(1),mlc(2),mlc(3)),rk ) * pi / &
            (180._rk*(10._rk**jj-1))
    eta   = Real( ang(3,mlc(1),mlc(2),mlc(3)),rk ) * pi / &
            (180._rk*(10._rk**jj-1))

    n = [cos(phi)*sin(eta) , sin(phi)*sin(eta) , cos(eta) ] 
    n = n / sqrt(sum(n*n))

    fcrit =  minval(crit(0:kk-1,0:kk_phi-1,0:kk_eta-1))

    If (mmdbg) then
       write(un_lf,*)
       Write(un_lf,"(A,I0,A) ")" Solution converged after :",jj," iterations"
       Write(un_lf,"(A,F17.6)")" With final citerion      :",fcrit
       
       Write(un_lf,"(A,F17.6)")" With final epsilon       :",&
            crit_min(jj-1)-crit_min(jj)
       Write(un_lf,"(A,3F17.9,A)")" Final rotation angle  is : ",alpha
       Write(un_lf,"(A,3F17.9,A)")" Final rotation vector is : ",n
       Write(un_lf,"(4(A,F17.9),A)")"B(rot([",n(1),",",n(2),",",n(3),"],",alpha,"))"
       Write(un_lf,*)
    End If

    !=========================================================================
    !== Inlining of EE =======================================================

    aa = rot_alg(n,alpha)
    BB = tra_R6(aa)
    EE = matmul(matmul(transpose(BB),EE),BB)

    If (mmdbg) then
       Call Write_matrix(un_lf, "Backrotated anisotropic stiffness CR_1", 'std', 'MPa', EE)
      !  call write_real_matrix(un_lf,EE_orig,6_ik,6_ik,&
      !       "eff_S","wxmaxima")
      !  Call Write_matrix(un_lf, "Backrotated anisotropic stiffness CR_1", 'std', 'MPa', EE_orig)
      !  call write_real_matrix(un_lf,EE,6_ik,6_ik,&
      !       "eff_S_CR_1","wxmaxima")
      !  Call Write_matrix(un_lf, "Backrotated anisotropic stiffness CR_1", 'std', 'MPa', EE)
    End If

    !=========================================================

    If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) < EE(3,3))         ) then

       If (mmdbg) write(un_lf,*)"123"

    Else If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) > EE(3,3)) ) then

       If (mmdbg) write(un_lf,*)"132"

       !** 132 => 123 *****************
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) > EE(3,3)) ) then

       If (mmdbg) write(un_lf,*)"231"

       !** 231 => 132 *****************
       n = aa(:,2)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)

       !** 132 => 123 *****************
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) < EE(3,3)) ) then

       If (mmdbg) write(un_lf,*)"213"

       !** 213 => 123 *****************
       n = aa(:,3)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) < EE(3,3)) ) then

       If (mmdbg) write(un_lf,*)"312"

       !** 312 => 132 *****************
       n = aa(:,3)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

       !** 132 => 123 *****************
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)  

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) > EE(3,3)) ) then

       If (mmdbg) write(un_lf,*)"321"

       !** 321 => 123 *****************
       n = aa(:,2)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)  

    End If

    BB = tra_R6(aa)
    EE = matmul(matmul(transpose(BB),EE_Orig),BB)

    If (mmdbg) then 
       Call Write_matrix(un_lf, "Final coordinate system CR_1", 'std', mat=aa)
      !  call write_real_matrix(un_lf,aa,3_ik,3_ik,&
      !       "R3_trafo","wxmaxima")
      !  call write_real_matrix(un_lf, EE,6_ik,6_ik,&
      !       "Inlined anisotropic stiffness CR_1")
      !  call check_sym(A=EE,name="Inlined anisotropic stiffness CR_1")
    End If

    fdir = aa

  End subroutine opt_mono

  !============================================================================
  !> Subroutine which optimizes a 6x6 marix according to an orthotropic
  !> criteria
  !>
  !> -- CR2 Orthotropic ----------------------------------------
  !>  Sum over Aij*Aij + A45**2 + A46**2 + A56**2
  !>  with i in [1,3] and j in [4,6]
  !>  +-                       -+
  !>  | a11 a12 a13 A14 A15 A16 |
  !>  | a12 a22 a23 A24 A25 A26 |
  !>  | a13 a12 a33 A34 A35 A36 |
  !>  | a14 a12 a13 a44 A45 A46 |
  !>  | a15 a12 a13 a14 a55 A56 |
  !>  | a16 a12 a13 a14 a15 a66 |
  !>  +-                       -+
  !>
  subroutine opt_ortho(EE_orig, EE, fdir, fcrit)

    Real(kind=rk)   , Dimension(6,6)    , intent(in)  :: EE_orig
    Real(kind=rk)   , Dimension(6,6)    , intent(out) :: EE     
    real(kind=rk)                       , intent(out) :: fcrit
    Real(kind=rk)   , Dimension(3,3)    , intent(out) :: fdir

    Integer(kind=ik)                                   :: kk_eta, kk_phi, kk
    Integer(kind=ik)                                   :: ii_eta, ii_phi, ii, jj
 
   Real(kind=rk)                                       :: alpha, phi, eta

    Real   (kind=rk), Dimension(3)                    :: n
    Real   (kind=rk), Dimension(3,3)                  :: aa
    Real   (kind=rk), Dimension(6,6)                  :: BB, tmp_r6x6

    Integer(kind=ik), Dimension(:,:,:,:), ALLOCATABLE :: ang
    Real   (kind=rk), Dimension(:,:,:)  , ALLOCATABLE :: crit
    Real   (kind=rk), Dimension(0 :16)                :: crit_min
    Integer         , Dimension(3)                    :: s_loop,e_loop, mlc

    Real   (kind=rk), Parameter                       :: num_zero = 1.E-9_rk
    Real   (kind=rk), Parameter                       :: pi  = acos(-1._rk)
    !**************************************************************************

    ALLOCATE(ang(3,0:180,0:180,0:90))
    ALLOCATE(crit( 0:180,0:180,0:90))

    ee = ee_orig

    kk_eta = 0
    kk_phi = 0
    kk = 0

    !** Initial parameter space sweep *****************************************
    Do ii_eta = 0 , 90 , 1

       kk_phi = 0

       Do ii_phi = 0 , 180 , 1

          kk = 0

          Do ii = 0 , 180 , 1

             alpha = Real(ii,rk)     * pi / 180._rk
             phi   = Real(ii_phi,rk) * pi / 180._rk
             eta   = Real(ii_eta,rk) * pi / 180._rk

             n = [cos(phi)*sin(eta) , sin(phi)*sin(eta) , cos(eta) ]
             n = n / sqrt(sum(n*n))

             aa = rot_alg(n,alpha)
             BB = tra_R6(aa)
             tmp_r6x6 = matmul(matmul(transpose(BB),EE),BB)

             ang(:,kk,kk_phi,kk_eta)  = [ii,ii_phi,ii_eta]

             !-- CR2 Orthotropic ---------------------------------------
             crit(kk,kk_phi,kk_eta) = (&
                  sum(tmp_r6x6(1:3,4:6) * tmp_r6x6(1:3,4:6)) + &
                  sum(tmp_r6x6( 4 ,5:6) * tmp_r6x6( 4 ,5:6)) + &
                  tmp_r6x6( 5 , 6 ) * tmp_r6x6( 5 , 6 )     )

             kk = kk + 1
          End Do
          kk_phi = kk_phi + 1
       end Do
       kk_eta = kk_eta + 1
    end Do

    !==========================================================================
    !== Iteration of Crit =====================================================

    crit_min = 0._rk
    crit_min(0) = minval(crit)

    mlc = minloc(crit)-1
    If (mmdbg) then
       write(un_lf,FMT_MSG_A3I0)'Initial Minloc  CR_1 : ',mlc
       write(un_lf,FMT_MSG_AF0) 'Initial Minimum CR_1 : ',crit_min(0)
    End If

    jj = 1

    Do 

       mlc = minloc(crit(0:kk-1,0:kk_phi-1,0:kk_eta-1))-1

       s_loop = (ang(:,mlc(1),mlc(2),mlc(3))-1)*10
       e_loop = (ang(:,mlc(1),mlc(2),mlc(3))+1)*10

       kk_eta = 0
       kk_phi = 0
       kk = 0

       If (mmdbg) then
          write(un_lf,FMT_MSG_AI0) 'Iteration            : ',jj
          
          write(un_lf,FMT_MSG_A3I0)'Loop start           : ',s_loop
          write(un_lf,FMT_MSG_A3I0)'Loop end             : ',e_loop
       End If

       Do ii_eta = s_loop(3), e_loop(3)

          kk_phi = 0

          Do ii_phi = s_loop(2), e_loop(2)

             kk = 0

             Do ii = s_loop(1), e_loop(1)

                alpha = Real(ii,rk)     * pi / (180._rk*(10._rk**jj))
                phi   = Real(ii_phi,rk) * pi / (180._rk*(10._rk**jj))
                eta   = Real(ii_eta,rk) * pi / (180._rk*(10._rk**jj))

                n = [cos(phi)*sin(eta) , sin(phi)*sin(eta) , cos(eta) ] 
                n = n / sqrt(sum(n*n))

                aa = rot_alg(n,alpha)
                BB = tra_R6(aa)
                tmp_r6x6 = matmul(matmul(transpose(BB),EE),BB)

                ang(:,kk,kk_phi,kk_eta)  = [ii,ii_phi,ii_eta]

                !-- CR2 Orthotropic ---------------------------------------
                crit(kk,kk_phi,kk_eta) =  (&
                     sum(tmp_r6x6(1:3,4:6) * tmp_r6x6(1:3,4:6)) + &
                     sum(tmp_r6x6( 4 ,5:6) * tmp_r6x6( 4 ,5:6)) + &
                     tmp_r6x6( 5 , 6 ) * tmp_r6x6( 5 , 6 )     )
                
                kk = kk + 1

             End Do
             kk_phi = kk_phi + 1
          end Do
          kk_eta = kk_eta + 1
       end Do

       crit_min(jj) = minval(crit(0:kk-1,0:kk_phi-1,0:kk_eta-1))

       If (mmdbg) then
          write(un_lf,FMT_MSG_AF0)'Minimum CR_1         : ',crit_min(jj)
       End If

       If ( (abs(crit_min(jj-1)-crit_min(jj)) < num_zero) .OR. (jj > 16)) Exit

       jj = jj + 1

    End Do

    mlc = minloc(crit(0:kk-1,0:kk_phi-1,0:kk_eta-1))-1

    alpha = Real( ang(1,mlc(1),mlc(2),mlc(3)),rk ) * pi / &
            (180._rk*(10._rk**jj-1))
    phi   = Real( ang(2,mlc(1),mlc(2),mlc(3)),rk ) * pi / &
            (180._rk*(10._rk**jj-1))
    eta   = Real( ang(3,mlc(1),mlc(2),mlc(3)),rk ) * pi / &
            (180._rk*(10._rk**jj-1))

    n = [cos(phi)*sin(eta) , sin(phi)*sin(eta) , cos(eta) ] 
    n = n / sqrt(sum(n*n))

    fcrit =  minval(crit(0:kk-1,0:kk_phi-1,0:kk_eta-1))

    If (mmdbg) then
       write(un_lf,*)
       Write(un_lf,"(A,I0,A) ")" Solution converged after :",jj," iterations"
       Write(un_lf,"(A,F17.6)")" With final citerion      :",fcrit
       
       Write(un_lf,"(A,F17.6)")" With final epsilon       :",&
            crit_min(jj-1)-crit_min(jj)
       Write(un_lf,"(A,3F17.9,A)")" Final rotation angle  is : ",alpha
       Write(un_lf,"(A,3F17.9,A)")" Final rotation vector is : ",n
       Write(un_lf,*)
    End If

    !=========================================================================
    !== Inlining of EE =======================================================

    aa = rot_alg(n,alpha)
    BB = tra_R6(aa)
    EE = matmul(matmul(transpose(BB),EE),BB)

    If (mmdbg) then
       Call Write_matrix(un_lf, "Backrotated anisotropic stiffness CR_1", 'std', 'MPa', EE)
    End If
    !=========================================================

    If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) < EE(3,3)) ) then

       If (mmdbg) write(un_lf,*)"123"

    Else If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) > EE(3,3)) ) then

       If (mmdbg) write(un_lf,*)"132"

       !** 132 => 123 *****************
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) > EE(3,3)) ) then

       If (mmdbg) write(un_lf,*)"231"

       !** 231 => 132 *****************
       n = aa(:,2)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)

       !** 132 => 123 *****************
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) < EE(3,3)) ) then

       If (mmdbg) write(un_lf,*)"213"

       !** 213 => 123 *****************
       n = aa(:,3)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) < EE(3,3)) ) then

       If (mmdbg) write(un_lf,*)"312"

       !** 312 => 132 *****************
       n = aa(:,3)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

       !** 132 => 123 *****************
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)  

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) > EE(3,3)) ) then

       If (mmdbg) write(un_lf,*)"321"

       !** 321 => 123 *****************
       n = aa(:,2)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)  

    End If

    BB = tra_R6(aa)
    EE = matmul(matmul(transpose(BB),EE_Orig),BB)

    If (mmdbg) then 
       Call Write_matrix(un_lf, "Final coordinate system CR_1", 'std', mat=aa)

       Call Write_matrix(un_lf, "Inlined anisotropic stiffness CR_1", 'std', 'MPa', EE)
       CALL check_sym(un_lf, EE, "Inlined anisotropic stiffness CR_1")
    End If

    fdir = aa

  End subroutine opt_ortho

  !============================================================================
  !> Subroutine which calculates a histogram
  subroutine calc_hist(phi, hist, csize)

    Real(kind=rk), dimension(:,:,:) , intent(in)              :: phi

    Integer(kind=ik) , dimension(:), Allocatable, intent(out) :: hist

    Integer(kind=ik) , optional                               :: csize

    Integer(kind=ik) , dimension(3)                           :: sphi

    Integer(kind=ik) :: hist_lbound, hist_ubound, ii,jj,kk

    Real(kind=ik)                                             :: ics_loc

    !--------------------------------------------------------------------------

    if (present(csize)) then
       ics_loc = 1._rk / Real(csize,rk)
    Else
       ics_loc = 1._rk
    End if


    hist_lbound = NInt(minval(phi) * ics_loc,ik)
    hist_ubound = NInt(maxval(phi) * ics_loc,ik)

    Write(*, *)
    Write(*, *) 'Minimal class value in phi = ',hist_lbound
    Write(*, *) 'Maximal class value in phi = ',hist_ubound

    Allocate(hist(hist_lbound:hist_ubound))
    hist = 0

    sphi = shape(phi)

    Do kk = 1, sphi(3)
       Do jj = 1, sphi(2)
          Do ii = 1, sphi(1)
             
             hist(NInt(phi(ii,jj,kk) * ics_loc,ik)) = &
             hist(NInt(phi(ii,jj,kk) * ics_loc,ik)) + 1
             
          End Do
       End Do
    End Do
    
  End Subroutine calc_hist

  !============================================================================
  !> Subroutine which calculates the area of a triangulated surface
  !>
  !> Subroutine which calculates the area of a triangulated surface given as
  !> node and element list. The nodelist is supposed to be ordered 
  !> consecutively starting at 1.
  Subroutine area(node_ls,cell_ls,arr)

    Real(kind=rk)   , Dimension(:,:), Allocatable, Intent(In)  :: node_ls
    Integer(kind=ik), Dimension(:,:), Allocatable, Intent(in)  :: cell_ls
    Real(Kind=rk)                                , Intent(Out) :: arr

    Real(Kind=rk), Dimension(3) :: c, a, cxa

    Integer(Kind=ik)            :: ii

    !--------------------------------------------------------------------------

    arr = 0._rk

    Do ii = 1, size(cell_ls(1,:))

       c = node_ls(:,cell_ls(2,ii))-node_ls(:,cell_ls(1,ii))
       a = node_ls(:,cell_ls(3,ii))-node_ls(:,cell_ls(1,ii))
     
       cxa(1) = c(2)*a(3)-c(3)*a(2)
       cxa(2) = c(3)*a(1)-c(1)*a(3)
       cxa(3) = c(1)*a(2)-c(2)*a(1)

       arr = arr + 0.5_rk * sqrt(cxa(1)*cxa(1)+cxa(2)*cxa(2)+cxa(3)*cxa(3))
 
    End Do

  End Subroutine area

  !============================================================================
  !> Subroutine which sorts a number of arrays acccording to a leading one
  Subroutine sort_arrays(lead,Real_1D)

    Integer(kind=ik), Dimension(:)      , Intent(inOut) :: lead

    Real(kind=rk), Dimension(:,:)  , optional, Intent(inOut) :: Real_1D

    Integer(kind=ik), Dimension(3)                      :: lb

    !**************************************************************************
    !* Variables declaration    
    Logical                                  :: unsorted
    Integer                                  :: ii,dim
    Real(Kind=rk)                            :: tmp_r

    !==========================================================================

    lb  = 1
    dim = size(lead)

    If (present(Real_1D)) then

       lb(1:2) = lbound(Real_1D)
write(*,*)lb(1:2)
!!$       if (dim /= size(Real_1D(lbound(Real_1D(:,:),1),:))) then
!!$          Write(*,*)"Only Arrays with equal sizes can be ordered"
!!$          Goto 1000
!!$       End if
    End If

!!$if ( present(Real_2D) )
!!$    if (dim /= size(Real_1D(1,:))) then
!!$       Write(*,*)"Only Arrays with equal sizes can be ordered"
!!$       Goto 1000
!!$    End if

    unsorted = .True.

    Do While (unsorted)

       unsorted = .false.

       Do ii = 2, dim
          
          If (lead(ii-1) > lead(ii)) Then

             tmp_r      = lead(ii-1)
             lead(ii-1) = lead(ii)
             lead(ii)   = tmp_r

!!$             tmp_r        = real_1(ii-1)
!!$             real_1(ii-1) = real_1(ii)
!!$             real_1(ii)   = tmp_r
!!$
!!$             tmp_i       = int_1(ii-1)
!!$             int_1(ii-1) = int_1(ii)
!!$             int_1(ii)   = tmp_i

             unsorted = .True.

          End IF

       End Do

    End Do

  End Subroutine sort_arrays

End Module math
