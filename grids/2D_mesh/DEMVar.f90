!##############################################################       
!     DEMVar
!       + module for the geeral variables (DEM defintion and values and
!       YAMS parameters)
!
!
!  Author : F. Gillet-Chaulet; LGGE, Grenoble, France
!###########################################################     
      module DEMVar
         real,allocatable :: DEM(:,:)
         real :: xmin,ymin,dx,dy,noval,R, threshold
         real,dimension(3) :: Param, Param2  !hmin,hmax,err
         integer :: nx,ny
      end module
