module utility
   use iso_fortran_env, wp => real64
   implicit none
contains
   subroutine cross(lhs, rhs, res)
      implicit none
      real(wp), intent(in) :: lhs(3), rhs(3)
      real(wp), intent(out) :: res(3)
      res(1) = lhs(2) * rhs(3) - lhs(3) * rhs(2)
      res(2) = lhs(3) * rhs(1) - lhs(1) * rhs(3)
      res(3) = lhs(1) * rhs(2) - lhs(2) * rhs(1)
   end subroutine cross

   subroutine normalize(v)
      implicit none
      real(wp), intent(inout) :: v(3)
      v = v / norm2(v)
   end subroutine normalize

   subroutine localglobalrot(normal, m)
      implicit none
      real(wp), intent(in) :: normal(3)
      real(wp), intent(out) :: m(3,3)
      real(wp)             :: xvec(3), yvec(3)
      xvec = [1.0_wp, 0.0_wp, 0.0_wp]
      if (norm2(xvec - normal) < 0.00001_wp) then
         xvec = [1.0_wp, 1.0_wp, 0.0_wp] / sqrt(2.0_wp)
      end if
      xvec = xvec-dot_product(xvec, normal)/dot_product(normal,normal)*normal
      xvec = xvec/norm2(xvec)
      call cross(normal, xvec, yvec)
      yvec = yvec/norm2(yvec)
      m(:,1) = xvec
      m(:,2) = yvec
      m(:,3) = normal
   end subroutine localglobalrot
end module utility
