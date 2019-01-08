program renderer
   use iso_fortran_env, wp => real64
   use geometry
   use rendering
   implicit none
   type(Camera)   :: cam
   real(wp)       :: pos(3), dir(3), up(3)
   real(wp)       :: fov, length, colours(6,3)
   integer        :: res, argc, err
   character(:),allocatable :: fname
   type(Sphere)   :: spheres(3)
   type(Plane)    :: planes(2)
   type(Triangle) :: triangles(0)
   type(Material) :: materials(6)

!     -y -z
!      | /
!      |/
! x ---o
   pos = [ 0.0_wp, 0.0_wp, 0.0_wp ]
   dir = [ 0.0_wp, 0.0_wp,-1.0_wp ]
   up  = [ 0.0_wp,-1.0_wp, 0.0_wp ]
   fov = 45.0_wp ! degrees
   length = 1.0_wp
   res = 512

   argc = command_argument_count()
   call read_arg(1, fname, err)
   if (err /= 0) stop 'requires filename.tga as first argument'

   cam = make_camera(pos, dir, up, fov, length, res, res)

   colours(1,1:3) = [255.0_wp, 0.0_wp, 0.0_wp]
   colours(2,1:3) = [0.0_wp, 255.0_wp, 0.0_wp]
   colours(3,1:3) = [0.0_wp, 0.0_wp, 125.0_wp]
   colours(4,1:3) = [128.0_wp, 128.0_wp, 128.0_wp]
   colours(5,1:3) = [255.0_wp, 255.0_wp, 0.0_wp]
   colours(6,1:3) = [0.0_wp, 0.0_wp, 50.0_wp]
   materials(1) = make_material(reshape(colours(1,1:3),[3]), &
                                .true., .false., 5.0_wp)
   materials(2) = make_material(reshape(colours(2,1:3),[3]), &
                                .true., .false., 0.0_wp)
   materials(3) = make_material(reshape(colours(3,1:3),[3]), &
                                .true., .false., 0.0_wp)
   materials(4) = make_material(reshape(colours(4,1:3),[3]), &
                                .true., .false., 0.0_wp)
   materials(5) = make_material(reshape(colours(5,1:3),[3]), &
                                .false., .true., 0.0_wp)
   materials(6) = make_material(reshape(colours(6,1:3),[3]), &
                                .true., .false., 0.0_wp)
   planes(1) = make_plane(normal=[0.0_wp, 1.0_wp, 0.0_wp],&
                          pos   =[0.0_wp, 0.2_wp, 0.0_wp],&
                          m     =materials(1))
   planes(2) = make_plane(normal=[0.0_wp,-1.0_wp, 0.0_wp],&
                          pos   =[0.0_wp,-1.45_wp, 0.0_wp],&
                          m     =materials(6))
   planes(2)%lum = 20.0_wp
   spheres(1) = make_sphere(pos =[0.0_wp, 0.0_wp, -7.0_wp], r=1.0_wp,&
                            m   =materials(3))
   spheres(2) = make_sphere(pos =[0.2_wp, 0.1_wp, -1.0_wp], r=0.1_wp,&
                            m   =materials(4))
   spheres(3) = make_sphere(pos =[-0.5_wp, 0.8_wp, -2.0_wp], r=0.1_wp,&
                            m   =materials(5))
   spheres(3)%lum = 255.0_wp

   call render_cam(cam, fname, spheres, planes, triangles)
contains
   subroutine read_arg(i,arg, err)
      use iso_fortran_env
      implicit none
      integer,intent(in):: i
      integer,intent(inout) :: err
      character(len=:),allocatable,intent(inout) :: arg
      integer :: l
      intrinsic :: get_command_argument
      if (allocated(arg)) deallocate(arg)
      call get_command_argument(i,length=l,status=err)
      if (err /= 0) return
      allocate(character(len=l) :: arg)
      call get_command_argument(i,value=arg)
   end subroutine read_arg
end program


