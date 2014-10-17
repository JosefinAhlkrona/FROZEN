!##############################################################       
!     READ_DEM
!       + read the velocity data used to create the metric
!       + has to initialise : nx,ny,xmin,ymin,dx,dy,noval,dem(nx,ny)
!       +  Typically to be modified by the user in function of the data set format
!
! INPUT : "../Data/GreenlandVelocity.dem" Greenland velocity DEM produced from
!            I. Joughin data sets available in the nsidc website
! OUTPUT: nx,ny,xmin,ymin,dx,dy,noval,dem(nx,ny)
!
!  Author : F. Gillet-Chaulet; LGGE, Grenoble

! Changed by Evan Gowan to work with his files
!##############################################################       
      subroutine  READ_DEM
      use DEMVar
        ! real,allocatable :: DEM(:,:)
        ! real :: xmin,ymin,dx,dy,noval,R, threshold
        ! real,dimension(3) :: Param, Param2  !hmin,hmax,err
        ! integer :: nx,ny

      implicit none

      integer :: i,j
      real :: u,v,e

	character(len=255) :: filename

	call getarg(1,filename)

	open(unit=40, file=filename, access="sequential", form="formatted", status="old")

	read(40,*) nx, ny, xmin, ymin, dx, dy, noval



     Allocate(dem(nx,ny))

      Do i=1,nx
     		Do j=1,ny

           read(40,*) e
           if (e==noval) then
                   dem(i,j)=noval
           else
                    dem(i,j)=e
           endif

		write(678,*) real(i-1)*dx +xmin, real(j-1)*dy +ymin, dem(i,j)

        End do
     End do
     close(40)

    End
