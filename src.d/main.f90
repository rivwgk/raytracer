program renderer
   use iso_fortran_env, wp => real64
   use geometry
   use rendering
   implicit none
   type(Camera)   :: cam
   real(wp)       :: pos(3), dir(3), up(3)
   real(wp),dimension(3) :: red,blue,green,yellow,grey,lightblue,white
   integer        :: argc, err
   character(:),allocatable :: fname
   type(Sphere)   :: spheres(3)
   type(Plane)    :: planes(2)
   type(Material) :: materials(6)
   real(wp),parameter :: length = 1.0_wp, fov = 30.0_wp
   integer,parameter :: res = 512

! y  z
! | /
! |/
! o---x
   pos = [ 0.0_wp, 0.0_wp, 0.0_wp ]
   dir = [ 0.0_wp, 1.0_wp, 0.0_wp ]
   up  = [ 0.0_wp, 0.0_wp, 1.0_wp ]

   argc = command_argument_count()
   if (argc /= 1) stop'requires filename.tga as first and only argument'
   call read_arg(1, fname, err)
   if (err /= 0) stop'requires filename.tga as first and only argument'

   cam = make_camera(pos, dir, up, fov, length, res, res)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   red = [1.0_wp, 0.0_wp, 0.0_wp]
   blue = [0.0_wp, 0.0_wp, 1.0_wp]
   green = [0.0_wp, 1.0_wp, 0.0_wp]
   yellow = [1.0_wp, 1.0_wp, 0.0_wp]
   grey = [0.6_wp, 0.6_wp, 0.6_wp]
   lightblue = [0.3_wp, 0.728_wp, 0.993_wp]
   white = [1.0_wp, 1.0_wp, 1.0_wp]

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   materials(1) = make_material(red, .false., .false., 2.0_wp)
   materials(2) = make_material(green, .true., .false., 0.0_wp)
   materials(3) = make_material(blue, .false., .false., 8.0_wp)
   materials(4) = make_material(grey, .true., .false., 0.0_wp)
   materials(5) = make_material(yellow, .true., .false., 8.0_wp)
   materials(6) = make_material(lightblue, .true., .false., 0.0_wp)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   planes(1) = make_plane(normal=[0.0_wp, 1.0_wp, 0.0_wp],&
                          pos   =[0.0_wp,-0.15_wp, 1.0_wp],&
                          m     =materials(1))
   planes(2) = make_plane(normal=[0.0_wp,-1.0_wp,0.0_wp],&
                          pos   =[0.0_wp,10.0_wp,1.0_wp],&
                          m     =materials(6))
   planes(2)%lum = 0.5_wp
   spheres(1) = make_sphere(pos =[ 0.0_wp,-0.05_wp,  1.2_wp], r=0.1_wp,&
                            m   =materials(3))
   spheres(2) = make_sphere(pos =[ 0.2_wp, 0.13_wp,  1.25_wp], r=0.15_wp,&
                            m   =materials(4))
   spheres(3) = make_sphere(pos =[-0.1_wp, 0.1_wp, 1.1_wp], r=0.05_wp,&
                            m   =materials(5))
   spheres(3)%lum = 1.0_wp

   call render_cam(cam, fname, spheres, planes)
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
