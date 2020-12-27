module geometry
   use iso_fortran_env, wp => real64
   use ieee_arithmetic
   use utility
   implicit none

   type Ray
      real(wp) :: o(3), d(3)
      real(wp) :: t
   contains
      procedure, pass(r) :: at => point_on_ray
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

   type ObjGeometry
      real(wp)       :: lum
      type(Material) :: m
   end type ObjGeometry

   type, extends(ObjGeometry) :: Sphere
      real(wp) :: p(3)
      real(wp) :: r
   contains
      procedure, pass(s) :: intersect => sphere_intersect
   end type Sphere

   type, extends(ObjGeometry) :: Plane
      real(wp) :: n(3), p(3)
   contains
      procedure, pass(p) :: intersect => plane_intersect
   end type Plane

   type Object
      type(ObjGeometry) :: geom
   end type Object
contains
   subroutine point_on_ray(r, p)
      implicit none
      class(Ray), intent(in) :: r
      real(wp), intent(out) :: p(3)
      p = r%o + r%t * r%d
   end subroutine

   subroutine plane_intersect(p, r, hit, isect)
      implicit none
      class(Plane), intent(in)   :: p
      type(Ray), intent(inout)   :: r
      logical, intent(inout)     :: hit
      type(Intersection), intent(inout) :: isect
      real(wp)                   :: a,b,c,t

      c = dot_product(p%n, r%d)
      if (abs(c) .le. 10. * epsilon(t)) then
         return ! ray nearly parallel to plane, no hit
      endif

      a = dot_product(p%n, p%p - r%o)
      t = (a / c)
      if (epsilon(t) .lt. t .and. t .lt. r%t) then
         r%t = t
         hit = .true.
         call r%at(isect%pos)
         isect%n = p%n
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
      real(wp)                   :: t1, t2, dd, root, q(3), solution

      dd = dot_product(r%d,r%d) ! the direction is normalized, this should be 1
      q = s%p - r%o
      root = (dot_product(r%d,q)/dd)**2 - (dot_product(q,q) - s%r**2)/dd
      if (root .lt. 0) then ! no real solution
         return
      endif
      t1 = dot_product(r%d,q)/dd - root
      t2 = dot_product(r%d,q)/dd + root
      solution = min(t1, t2)
      if (solution .le. 0) then
         solution = max(t1, t2)
      endif
      if ((epsilon(0.0_wp) .lt. solution) .and. (solution .lt. r%t)) then
         hit = .true.
         r%t = solution
         call r%at(isect%pos)
         isect%n= isect%pos - s%p
         call normalize(isect%n)
         isect%mat = s%m
         isect%lum = s%lum
      endif
   end subroutine

   function make_plane(normal, pos, m) result(p)
      implicit none
      real(wp), intent(in)       :: normal(3), pos(3)
      type(Material), intent(in) :: m
      type(Plane)                :: p
      p%n = normal
      call normalize(p%n)
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
      real(wp), intent(in) :: normal(3),r_in(3),r_out(3)
      real(wp)             :: w,v(3)
      v = 2.0_wp*dot_product(r_in,normal)*normal-r_in
      if (norm2(v - r_out) .lt. 0.1) then
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
      v = 2.0_wp*dot_product(normal,r_in)*normal-r_in
      w = dot_product(r_out,v) ** m%shininess
   end function phong

   function absorb(m, normal, r_in, r_out) result(w)
      implicit none
      class(Material), intent(in) :: m
      real(wp), intent(in) :: normal(3), r_in(3), r_out(3)
      real(wp)             :: w
      w = 0.0_wp
   end function absorb

   subroutine intersect_all(r,hit,isect,planes,spheres)
      implicit none
      type(Ray), intent(inout)        :: r
      logical, intent(out)            :: hit
      type(Intersection), intent(out) :: isect
      type(Plane), intent(in)         :: planes(:)
      type(Sphere), intent(in)        :: spheres(:)
      integer :: k
      do k = 1, size(planes, 1)
         call planes(k)%intersect(r, hit, isect)
      end do
      do k = 1, size(spheres, 1)
         call spheres(k)%intersect(r, hit, isect)
      end do
   end subroutine

   recursive subroutine shade(isect,c,spheres,planes,w_out,n,depth,rands,qmc)
      implicit none
      type(Intersection), intent(in) :: isect
      type(Intersection)             :: is
      real(wp), intent(out)          :: c(3)
      type(Plane), intent(in)        :: planes(:)
      type(Sphere), intent(in)       :: spheres(:)
      real(wp), intent(in)           :: w_out(3)
      real(wp),parameter :: pi = 4.0_wp * atan(1.0_wp)
      real(wp),parameter :: square = pi * 0.5_wp * pi * 2.0_wp
      real(wp),parameter :: hemisphere = 4.0_wp/2.0_wp * pi
      integer, intent(in)            :: n, depth
      integer                        :: i,j,k,rep,max_rep
      real(wp)                       :: phi,theta,v(3),m(3,3),w
      real(wp)                       :: acc(3), racc(3),area
      real(wp)                       :: rands(:,:,:)
      type(Ray)                      :: r
      logical                        :: hit,qmc
      max_rep = size(rands,3)
      area = square
      acc = 0.0_wp
      do i = 1, n
         do rep = 1,max_rep

            theta = rands(1,i,rep)
            phi   = rands(2,i,rep)
            call eval_mc(real(i,wp),real(n,wp),theta,phi)
!           call eval_mc_H2(real(i,wp),real(n,wp),theta,phi)
!           call eval_mc_cos(real(i,wp),real(n,wp),theta,phi)
            v = [sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)]
            r%o = isect%pos
            call localglobalrot(isect%n, m)
            r%d = reshape(matmul(m, reshape(v, [3, 1])), [3])
            call normalize(r%d)
            r%t = 1000000.0_wp
            hit = .false.
            call intersect_all(r, hit, is, planes, spheres)
            if (hit) then
               if (depth .gt. 0) then
                  if (is%lum .eq. 0.0_wp) then
                     call shade(is,racc,spheres,planes,-r%d,n,depth-1,rands,qmc)
                  else
                     racc = is%lum * is%mat%c
                  end if
                  acc = acc + racc
               else
                  acc = acc + is%lum * cos(theta) * &
                        isect%mat%c * isect%mat%shade(isect%n, r%d, w_out)
               end if
            end if

         end do
      end do
      c = acc * area /(n*max_rep)
   contains
      subroutine eval_midpoint(n,i,m,j,theta,phi)
         implicit none
         real(wp),intent(in) :: n,i,m,j
         real(wp),intent(out) :: theta,phi
         theta = i / n * 0.5_wp * pi
         phi = j / m * 2.0_wp * pi
      end subroutine
      ! theta-phi, standard monte-carlo
      subroutine eval_mc(n,i,theta,phi)
         implicit none
         real(wp),intent(in) :: n,i
         real(wp),intent(inout) :: theta,phi
         theta = theta * 0.5_wp * pi
         phi = phi * 2.0_wp * pi
      end subroutine
      ! 
      subroutine eval_mc_H2(n,i,theta,phi)
         implicit none
         real(wp),intent(in) :: n,i
         real(wp),intent(inout) :: theta,phi
         theta = acos(theta)
         phi = phi * 2.0_wp * pi
      end subroutine
      ! important sampling
      subroutine eval_mc_cos(n,i,theta,phi)
         implicit none
         real(wp),intent(in) :: n,i
         real(wp),intent(inout) :: theta,phi
         theta = acos(sqrt(theta))
         phi = phi * 2.0_wp * pi
      end subroutine
   end subroutine shade
end module geometry
