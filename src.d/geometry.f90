module geometry
   use iso_fortran_env, wp => real64
   implicit none

   type Ray
      real(wp) :: o(3), d(3)
      real(wp) :: t
   contains
      procedure, pass(arg) :: normalize => normalize_ray
   end type Ray

   type Material
      real(wp)     :: c(3)
      logical      :: is_diffuse = .false.
      logical      :: is_mirror  = .false.
      real(wp)     :: shininess  = 1.0_wp
      procedure(absorb), pointer, pass(m) :: shade
   end type Material

   type Intersection
      real(wp) :: pos(3), n(3), lum
      type(Material) :: mat
   end type Intersection

   type, abstract :: Object
      real(wp)       :: lum
      type(Material) :: m
   end type Object

   type, extends(Object) :: Sphere
      real(wp) :: p(3)
      real(wp) :: r
   contains
      procedure, pass(s) :: intersect => sphere_intersect
   end type Sphere

   type, extends(Object) :: Plane
      real(wp) :: n(3), p(3)
   contains
      procedure, pass(p) :: intersect => plane_intersect
   end type Plane

   type, extends(Object) :: Triangle
      real(wp) :: a(3), b(3), c(3)
   contains
      procedure, pass(t) :: intersect => triangle_intersect
   end type Triangle
contains
   pure subroutine cross(a, b, c)
      implicit none
      real(wp), intent(in) :: a(3), b(3)
      real(wp), intent(out) :: c(3)
      c(1) = a(2) * b(3) - a(3) * b(2)
      c(2) = a(3) * b(1) - a(1) * b(3)
      c(3) = a(1) * b(2) - a(2) * b(1)
   end subroutine

   subroutine point_on_ray(p, r)
      implicit none
      real(wp), intent(out) :: p(3)
      type(Ray), intent(in) :: r
      p = r%o + r%t * r%d
   end subroutine

   subroutine normalize_ray(arg)
      implicit none
      class(Ray), intent(inout) :: arg
      arg%d = arg%d / norm2(arg%d)
   end subroutine

   subroutine plane_intersect(p, r, hit, isect)
      implicit none
      class(Plane), intent(in)   :: p
      type(Ray), intent(inout)   :: r
      logical, intent(inout)     :: hit
      type(Intersection), intent(inout) :: isect
      real(wp)                   :: t

      if (dot_product(p%n, r%d) .le. epsilon(t)) then
      !if (abs(dot_product(p%n, r%d)) .le. epsilon(t)) then
         return ! ray nearly parallel to plane, no hit
      endif

      t = (dot_product(p%n,p%p - r%o) / dot_product(p%n, r%d))
      if (0 .lt. t .and. t .lt. r%t) then
         r%t = t
         hit = .true.
         call point_on_ray(isect%pos, r)
         isect%n= p%n
         isect%mat = p%m
         isect%lum = p%lum
      endif
   end subroutine

   subroutine sphere_intersect(s, r, hit, isect)
      implicit none
      class(Sphere), intent(in)  :: s
      type(Ray), intent(inout)   :: r
      logical, intent(inout)     :: hit
      type(Intersection), intent(inout) :: isect
      real(wp)                   :: t1, t2, dd, root, o(3), solution

      dd = dot_product(r%d, r%d) ! as the direction is normalized this is 1
      o = r%o - s%p
      root = (dot_product(r%d,o))**2 - (dot_product(o,o) - s%r**2)
      if (root .le. 0) then ! no real solution
         return
      endif
      t1 = dot_product(r%d,o) - root
      t2 = dot_product(r%d,o) + root
      solution = min(t1, t2)
      if (solution .le. 0) then
         solution = max(t1, t2)
      endif
      if (solution .le. 0) then
         return
      endif
      if (solution .le. r%t) then
         hit = .true.
         r%t = solution
         call point_on_ray(isect%pos, r)
         isect%n= isect%pos - s%p
         isect%n= isect%n/ norm2(isect%n)
         isect%mat = s%m
         isect%lum = s%lum
      endif
   end subroutine

   subroutine triangle_intersect(t, r, hit, isect)
      implicit none
      class(Triangle), intent(in):: t
      type(Ray), intent(inout)   :: r
      logical, intent(inout)     :: hit
      type(Intersection), intent(inout) :: isect
   end subroutine

   function make_plane(normal, pos, m) result(p)
      implicit none
      real(wp), intent(in)       :: normal(3), pos(3)
      type(Material), intent(in) :: m
      type(Plane)                :: p
      p%n = normal / norm2(normal)
      p%p = pos
      p%m = m
   end function make_plane

   function make_sphere(pos, r, m) result(s)
      implicit none
      real(wp), intent(in)       :: pos(3), r
      type(Material), intent(in) :: m
      type(Sphere)               :: s
      s%p = pos
      s%r = r
      s%m = m
   end function make_sphere

   function make_triangle(u, v, w) result(t)
      implicit none
      real(wp), intent(in)       :: u(3), v(3), w(3)
      type(Triangle)             :: t
      t%a = u
      t%b = v
      t%c = w
   end function make_triangle

   function make_material(color, is_diffuse, is_mirror, shininess) result(m)
      implicit none
      real(wp), intent(in) :: color(3), shininess
      logical, intent(in)  :: is_diffuse, is_mirror
      type(Material) :: m
      m%c = color
      m%shininess = shininess
      m%is_diffuse = is_diffuse
      m%is_mirror = is_mirror
      if (m%is_diffuse .and. m%is_mirror) then
         m%shade => absorb
      else if (m%is_diffuse) then
         m%shade => diffuse
      else if (m%is_mirror) then
         m%shade => mirror
      else
         m%shade => phong
      endif
   end function make_material

   function diffuse(m, normal, r_in, r_out) result(w)
      implicit none
      class(Material), intent(in) :: m
      real(wp), intent(in) :: normal(3), r_in(3), r_out(3)
      real(wp)             :: w
      w = dot_product(r_in, normal)
   end function diffuse

   function mirror(m, normal, r_in, r_out) result(w)
      implicit none
      class(Material), intent(in) :: m
      real(wp), intent(in) :: normal(3), r_in(3), r_out(3)
      real(wp)             :: w, v(3)
      ! o - <o,i>/<i,i>
      v = 2.0_wp * dot_product(normal, r_in) * normal - r_in
      if (norm2(v - r_out) .lt. 10 * epsilon(w)) then
         w = 1.0_wp
      else
         w = 0.0_wp
      endif
   end function mirror

   function phong(m, normal, r_in, r_out) result(w)
      implicit none
      class(Material), intent(in) :: m
      real(wp), intent(in) :: normal(3), r_in(3), r_out(3)
      real(wp)             :: w, v(3)
      v = 2.0_wp * dot_product(normal, r_in) * normal - r_in
      w = dot_product(v, r_out) ** m%shininess
   end function phong

   function absorb(m, normal, r_in, r_out) result(w)
      implicit none
      class(Material), intent(in) :: m
      real(wp), intent(in) :: normal(3), r_in(3), r_out(3)
      real(wp)             :: w
      w = 0.0_wp
   end function absorb

   subroutine shade(isect, c, spheres, planes, triangles, w_out, n)
      implicit none
      type(Intersection), intent(in) :: isect
      type(Intersection)             :: is
      real(wp), intent(out)          :: c(3)
      type(Plane), intent(in)        :: planes(:)
      type(Sphere), intent(in)       :: spheres(:)
      type(Triangle), intent(in)     :: triangles(:)
      real(wp), intent(in)           :: w_out(3)
      real(wp),parameter :: pi = 4.0_wp * atan(1.0_wp)
      integer, intent(in)            :: n
      integer                        :: i,k
      real(wp)                       :: phi,theta,v(3),w_in(3),m(3,3),w
      real(wp)                       :: acc(3)
      type(Ray)                      :: r
      logical                        :: hit
      acc = 0
      do i = 1, n
         call random_number(theta)
         call random_number(phi)
         theta = acos(theta)
         phi = phi * 2 * pi
         v(1) = sin(theta) * cos(phi)
         v(2) = sin(theta) * sin(phi)
         v(3) = cos(theta)
         r%o = isect%pos
         call localglobalrot(isect%n, m)
         w_in = reshape(matmul(m, reshape(v, [3, 1])), [3])
         w_in = w_in / norm2(w_in)
         r%d = w_in
         r%t = 100000.0_wp
         hit = .false.
         do k = 1, size(planes, 1)
            call planes(k)%intersect(r, hit, is)
         end do
         do k = 1, size(spheres, 1)
            call spheres(k)%intersect(r, hit, is)
         end do
         do k = 1, size(triangles, 1)
            call triangles(k)%intersect(r, hit, is)
         end do
         if (hit) then
            acc = acc + is%lum * cos(theta) * &! isect%mat%c * &
                  isect%mat%shade(isect%n, w_in, w_out)
         end if
      end do
      c = acc * pi**2 / n
   contains
      subroutine localglobalrot(normal, m)
         implicit none
         real(wp), intent(in) :: normal(3)
         real(wp), intent(out) :: m(3,3)
         real(wp)             :: xvec(3), yvec(3)
         xvec = [1.0_wp, 0.0_wp, 0.0_wp]
         if (norm2(xvec - normal) < 0.00001_wp) then
            xvec = [1.0_wp, 1.0_wp, 0.0_wp] / sqrt(2.0_wp)
         end if
         xvec = xvec - dot_product(xvec, normal) * normal
         xvec = xvec / norm2(xvec)
         call cross(normal, xvec, yvec)
         yvec = yvec / norm2(yvec)
         m(1:3,1) = xvec
         m(1:3,2) = yvec
         m(1:3,3) = normal
         m = transpose(m)
      end subroutine localglobalrot
   end subroutine shade
end module geometry
