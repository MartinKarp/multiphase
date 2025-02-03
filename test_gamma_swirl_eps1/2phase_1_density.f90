module user
  use neko
  implicit none

  real(kind=rp) :: eps
  real(kind=rp) :: gamma, u_max = 0.0

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%scalar_user_ic => set_ic
    u%scalar_user_f_vector => forcing
    u%material_properties => set_material_properties
    u%user_check => usr_check
  end subroutine user_setup

  !> Dummy user check
  subroutine usr_check(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    integer :: i
    real(kind=rp) :: x, y, u_tmp(1)
    u_tmp = 0
    do i = 1, u%size()
       x = u%dof%x(i,1,1,1)
       y = u%dof%y(i,1,1,1)
       u%x(i,1,1,1) = -2*sin(pi * x) **2 * sin(pi *y )* cos(pi*y)*cos(pi*t/4.0_xp)
       v%x(i,1,1,1) = 2*sin(pi * x) *cos(pi*x) * sin(pi *y )**2*cos(pi*t/4.0_xp)
       u_tmp(1) = max(u_tmp(1),sqrt(u%x(i,1,1,1)**2+v%x(i,1,1,1)**2))
    end do
    u_max = glmax(u_tmp,1)
   
  end subroutine usr_check


  subroutine set_material_properties(t, tstep, rho, mu, cp, lambda, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp), intent(inout) :: rho, mu, cp, lambda
    type(json_file), intent(inout) :: params
    real(kind = rp) :: del
    type(field_t), pointer :: u

    call json_get(params, "case.scalar.epsilon", eps)
    call json_get(params, "case.scalar.gamma", gamma)

    del = gamma*u_max
    lambda = eps*del
    !print *, lambda,u_max, 'yuuu'
    mu = 1.0_rp 
    rho = 1.0_rp
    cp = 1.0_rp
  end subroutine set_material_properties
  !> User initial condition
  subroutine set_ic(s, params)
    type(field_t), intent(inout) :: s
    type(json_file), intent(inout) :: params
    real(kind=rp) ::  rad, eps
    integer :: i

 
    call json_get(params, "case.scalar.epsilon",eps)
    
    do i = 1, s%dof%size()
       rad = sqrt((s%dof%x(i,1,1,1)-0.5)**2 +(s%dof%y(i,1,1,1)-0.75)**2)
       s%x(i,1,1,1) = 0.5*(1+tanh((rad-0.15_xp)/(2*eps)))
    end do

    if ((NEKO_BCKND_DEVICE .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_memcpy(s%x, s%x_d, s%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
    end if


  end subroutine set_ic

  !> Forcing
  subroutine forcing(f, t)
    class(scalar_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    integer :: i
    type(field_t), pointer :: u, v, w, s
    real(kind=rp) :: absgrad
    integer :: ind(4)
    type(field_t), pointer :: work1, work2, work3, work4
    s => neko_field_registry%get_field('s')

    call neko_scratch_registry%request_field(work1, ind(1))
    call neko_scratch_registry%request_field(work2, ind(2))
    call neko_scratch_registry%request_field(work3, ind(3))
    call neko_scratch_registry%request_field(work4, ind(4))
      
    call grad(work1%x, work2%x, work3%x, s%x, f%coef)
    
    call f%coef%gs_h%op(work1,GS_OP_ADD)
    call f%coef%gs_h%op(work2,GS_OP_ADD)
    call f%coef%gs_h%op(work3,GS_OP_ADD)
    call col2(work1%x, f%coef%mult, work4%size())
    call col2(work2%x, f%coef%mult, work4%size())
    call col2(work3%x, f%coef%mult, work4%size())
    !call field_vdot3(work4,3work1, work2, work3,work1, work2, work3, work4%size())



    do i = 1, work4%size()
       absgrad = sqrt(work1%x(i,1,1,1)**2+work2%x(i,1,1,1)**2+work3%x(i,1,1,1)**2)
       if (absgrad == 0.0 ) then 
          print *, 'warning', absgrad
          absgrad = 1e21
       end if
          
       work1%x(i,1,1,1) = - s%x(i,1,1,1)*(1-s%x(i,1,1,1))*(work1%x(i,1,1,1)/absgrad)
       work2%x(i,1,1,1) = - s%x(i,1,1,1)*(1-s%x(i,1,1,1))*(work2%x(i,1,1,1)/absgrad)
       work3%x(i,1,1,1) = - s%x(i,1,1,1)*(1-s%x(i,1,1,1))*(work3%x(i,1,1,1)/absgrad)
    end do
    !call cdtp(work4%x, work1%x, f%coef%drdx, f%coef%dsdx, f%coef%dtdx, f%coef)
    !call add2(f%s, work4%x,work4%size())
    !call cdtp(work4%x, work2%x, f%coef%drdy, f%coef%dsdy, f%coef%dtdy, f%coef)
    !call add2(f%s, work4%x,work4%size())
    !call cdtp(work4%x, work3%x, f%coef%drdz, f%coef%dsdz, f%coef%dtdz, f%coef)
    !call add2(f%s, work4%x,work4%size())
    !call cmult(f%s, gamma,work4%size())
    !call f%coef%gs_h%op(f%s,GS_OP_ADD,work4%size())
    !call col2(f%s, f%coef%Binv, work4%size())
    call dudxyz(work4%x, work1%x, f%coef%drdx, f%coef%dsdx, f%coef%dtdx, f%coef)
    !call f%coef%gs_h%op(work4,GS_OP_ADD)
    !call col2(work4%x, f%coef%mult, work4%size())
    call copy(f%s, work4%x,work4%size())
    call dudxyz(work4%x, work2%x, f%coef%drdy, f%coef%dsdy, f%coef%dtdy, f%coef)
    !call f%coef%gs_h%op(work4,GS_OP_ADD)
    !call col2(work4%x, f%coef%mult, work4%size())
    call add2(f%s, work4%x,work4%size())
    call dudxyz(work4%x, work3%x, f%coef%drdz, f%coef%dsdz, f%coef%dtdz, f%coef)
    !call f%coef%gs_h%op(work4,GS_OP_ADD)
    !call col2(work4%x, f%coef%mult, work4%size())
    call add2(f%s, work4%x,work4%size())
    absgrad = gamma*u_max
    call cmult(f%s, absgrad,work4%size())
    !call col2(f%s, f%coef%B, work4%size())
    !call f%coef%gs_h%op(f%s,GS_OP_ADD,work4%size())
    !call col2(f%s, f%coef%Binv, work4%size())


    call neko_scratch_registry%relinquish_field(ind)
    !if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
    !   .or. (NEKO_BCKND_OPENCL .eq. 1)) then
    !   call device_cmult2(f%u_d,v%x_d,Ta2Pr,f%dm%size())
    !   call device_cmult2(f%v_d,u%x_d,Ta2Pr,f%dm%size())
    !   call device_cmult2(f%w_d,s%x_d,rapr,f%dm%size())
    !else
    !   call cmult2(f%u,v%x,Ta2Pr,f%dm%size())
    !   call cmult2(f%v,u%x,Ta2Pr,f%dm%size())
    !   call cmult2(f%w,s%x,rapr,f%dm%size())
    !end if
  end subroutine forcing
end module user
