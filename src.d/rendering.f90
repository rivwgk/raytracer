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
      real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
      integer, intent(in) :: res_x, res_y
      type(Camera) :: cam
      cam%pos = pos
      cam%dir = dir
      cam%up = up
      cam%fov = fov * pi / 180.0_wp
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
      res_max = tan(cam%fov / 2.0_wp) * cam%length
      dx = 2 * res_max / cam%xdim
      dy = 2 * res_max / cam%ydim
      r%o = cam%pos
      r%d(1) = r%o(1) - res_max + x * dx
      r%d(2) = r%o(2) - res_max + y * dy
      r%d(3) = r%o(3) + cam%length
      r%t = huge(r%t)
      call normalize(r%d)
   end function

   subroutine render_cam(cam, filename, spheres, planes)
      implicit none
      type(Camera), intent(inout):: cam
      character(*), intent(in)   :: filename
      type(Sphere), intent(in)   :: spheres(:)
      type(Plane), intent(in)    :: planes(:)
      type(Ray)                  :: r
      integer                    :: i, j, k, dim, x, y
      integer, parameter         :: n = 2**8, depth = 0
      logical                    :: hit
      logical,parameter          :: qmc = .false.
      type(Intersection)         :: isect
      real(wp)                   :: c(3)
      real(wp),allocatable       :: rands(:,:,:)
      if (qmc) then
         dim = 1 
      else
         dim = 20
      endif
      allocate(rands(2,n,dim))
      if (qmc) then
         call halton_sequence(2,n,rands(:,:,1))
      else
         call random_number(rands)
      endif
      print*,'tracing rays...'
!$omp parallel do &
!$omp num_threads(8) &
!$omp shared(cam,spheres,planes,rands) &
!$omp private(r,isect,c,hit,i,j,k)
      do i = 1,cam%ydim
         do j = 1,cam%xdim
            hit = .false.
            r = generate_ray(cam, j, i)
            call intersect_all(r,hit,isect,planes,spheres)
            if (hit) then
               if (isect%lum .eq. 0.0_wp) then
                  call shade(isect,c,spheres,planes,-r%d,n,depth,rands(:,:,:),qmc)
               else
                  c = isect%lum * isect%mat%c
               endif
               cam%img(j, i, 1:3) = c
            else
               cam%img(j, i, 1:3) = 0
            endif
         end do
      end do
!$omp end parallel do
      deallocate(rands)
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
!     integer(byte),allocatable :: img(:,:,:)
      real(wp)                  :: sup, inf
      logical,parameter :: greyscale = .false.
      sup = maxval(cam%img)
      inf = minval(cam%img)
      print*,'Writing to disk: ' // filename,sup,inf
!     allocate(img(cam%xdim,cam%ydim,3))
!     img = nint(255 * (cam%img - inf) / (sup - inf), byte)
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
!           if (greyscale) then
!              write(u) img(j, i, 3)  ! b
!              write(u) img(j, i, 2)  ! g
!              write(u) img(j, i, 1)  ! r
!           else 
               write(u) nint(cam%img(j, i, 3)/sup * 255, byte) !img(j, i, 3)  ! b
               write(u) nint(cam%img(j, i, 2)/sup * 255, byte) !img(j, i, 2)  ! g
               write(u) nint(cam%img(j, i, 1)/sup * 255, byte) !img(j, i, 1)  ! r
!           endif
            write(u) int(Z'FF', byte)
         enddo
      enddo
      close(u)
!     deallocate(img)
   contains
      pure elemental integer function clamp(val,minval,maxval)
         implicit none
         integer,intent(in) :: val,minval,maxval
         clamp = max(minval,min(maxval,val))
      end function clamp
   end subroutine
end module rendering
