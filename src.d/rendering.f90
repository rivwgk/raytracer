module rendering
   use geometry
   use iso_fortran_env, wp => real64
   implicit none
   type Camera
      real(wp) :: pos(3), dir(3), up(3), left(3), right(3), top(3), bot(3)
      real(wp) :: fov, length
      integer :: xdim, ydim
      real(wp), allocatable :: img(:,:,:)
   end type
contains
   function make_camera(pos, dir, up, fov, length, res_x, res_y) result(cam)
      implicit none
      real(wp), intent(in) :: pos(3), dir(3), up(3)
      real(wp), intent(in) :: fov, length
      integer, intent(in) :: res_x, res_y
      type(Camera) :: cam
      cam%pos = pos
      cam%dir = dir
      cam%up = up
      cam%fov = fov * atan(1.0_wp) / 90.0_wp
      cam%length = length
      cam%xdim = res_x
      cam%ydim = res_y
      allocate(cam%img(res_x, res_y, 3), source=0.0_wp)
   end function

   function generate_ray(cam, x, y) result(r)
      implicit none
      type(Camera), intent(in)   :: cam
      type(Ray)                  :: r
      integer, intent(in)        :: x, y
      real(wp)                   :: dx, dy, res_max 
      res_max = tan(cam%fov) * cam%length
      dx = 2 * res_max / cam%xdim
      dy = 2 * res_max / cam%ydim
      r%o = [ 0., 0., 0. ]
      r%d(1) = r%o(1) -res_max + x * dx
      r%d(2) = r%o(2) +res_max - y * dy
      r%d(3) = r%o(3) + cam%length
      r%t = 10000.0_wp
      call r%normalize()
   end function

   subroutine render_cam(cam, filename, spheres, planes, triangles)
      implicit none
      type(Camera), intent(inout):: cam
      character(*), intent(in)   :: filename
      type(Sphere), intent(in)   :: spheres(:)
      type(Plane), intent(in)    :: planes(:)
      type(Triangle), intent(in) :: triangles(:)
      type(Ray)                  :: r
      integer                    :: i, j, k
      integer, parameter         :: n = 64
      logical                    :: hit
      type(Intersection)         :: isect
      real(wp)                   :: c(3)
      do i = 1,cam%ydim
         do j = 1,cam%xdim
            hit = .false.
            r = generate_ray(cam, j, i)
            do k = 1, size(planes, 1)
               call planes(k)%intersect(r, hit, isect)
            end do
            do k = 1, size(spheres, 1)
               call spheres(k)%intersect(r, hit, isect)
            end do
            do k = 1, size(triangles, 1)
               call triangles(k)%intersect(r, hit, isect)
            end do
            if (hit) then
               if (isect%lum .eq. 0.0_wp) then
                  call shade(isect, c, spheres, planes, triangles, -r%d, n)
               else
                  c = isect%lum !isect%mat%c
               endif
               cam%img(j, i, 1:3) = c
            else
               cam%img(j, i, 1:3) = 0
            endif
         end do
      end do

      call save_to_file(cam, filename)
   end subroutine

   subroutine save_to_file(cam, filename)
      implicit none
      type(Camera), intent(in)   :: cam
      character(*), intent(in)   :: filename
      integer                    :: i,j,u=42
      integer,parameter :: byte = selected_int_kind(2)
      integer(byte),parameter :: zero = 0_byte
      integer(byte),parameter :: urgb = 2_byte
      integer(byte),allocatable :: img(:,:,:)
      real(wp)                  :: sup, inf
      logical,parameter :: greyscale = .true.
      sup = maxval(cam%img)
      inf = minval(cam%img)
      print*,'Writing to disk: ' // filename,sup,inf
      allocate(img(cam%xdim,cam%ydim,3))
      img = nint(255 * (cam%img - inf) / (sup - inf), byte)
      open(newunit=u, file=filename, access='stream')
      write(u) zero
      write(u) zero
      write(u) urgb
      write(u) zero
      write(u) zero
      write(u) zero
      write(u) zero
      write(u) zero
      write(u) zero ! X origin
      write(u) zero
      write(u) zero ! Y origin
      write(u) zero
      write(u) int(ibits(cam%xdim, 0,  7), byte)
      write(u) int(ibits(cam%xdim, 8, 16), byte)
      write(u) int(ibits(cam%ydim, 0,  7), byte)
      write(u) int(ibits(cam%ydim, 8, 16), byte)
      write(u) int(32_byte, byte)
      write(u) zero
      do i=1,cam%ydim
         do j=1,cam%xdim
            if (greyscale) then
               write(u) img(j, i, 3)  ! b
               write(u) img(j, i, 2)  ! g
               write(u) img(j, i, 1)  ! r
            else 
               write(u) nint(cam%img(j, i, 3), byte) !img(j, i, 3)  ! b
               write(u) nint(cam%img(j, i, 2), byte) !img(j, i, 2)  ! g
               write(u) nint(cam%img(j, i, 1), byte) !img(j, i, 1)  ! r
            endif
            write(u) int(Z'FF', byte)
         enddo
      enddo
      close(u)
      deallocate(img)
   contains
      pure elemental integer function clamp(val,minval,maxval)
         implicit none
         integer,intent(in) :: val,minval,maxval
         clamp = max(minval,min(maxval,val))
      end function clamp
      subroutine hsltorgb(hsl, rgb)
         implicit none
         real(wp),intent(in) :: hsl(3)
         real(wp),intent(out) :: rgb(3)
         rgb = [f(hsl(1), hsl(2), hsl(3), 2.0_wp), &
                f(hsl(1), hsl(2), hsl(3), 0.0_wp), &
                f(hsl(1), hsl(2), hsl(3), 4.0_wp)]
      end subroutine hsltorgb
      function f(H, S, L, n) result(p)
         implicit none
         real(wp),intent(in) :: H,S,L,n
         real(wp) :: a,k,p
         k = modulo((n + H / 60.0_wp), 6.0_wp)
         a = S * min(L, 1._wp - L)
         p = L - a + 2._wp * a * max(min(k, 4.0_wp - k, 1.0_wp), 0.0_wp)
      end function f
   end subroutine
end module rendering
