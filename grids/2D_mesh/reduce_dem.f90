program reduce_dem

! this program reads in the DEM and outline from the previous steps, and outputs the parameters needed for the SIF file (assuming you want
! to use Grid2DInterpolator, which is a bilinear interpolation method

! input parameters include the contour file (in Cartesian coordinates), and the headered xyz large scale dem

! requires files READ_DEM.o and DEMVar.o

	use DEMVar
	implicit none
	character(len=255) :: filename_contour

	integer :: istat, start_x_counter, start_y_counter, end_x_counter, end_y_counter, x_counter, y_counter

	real :: min_contour_x, min_contour_y, max_contour_x, max_contour_y, x, y, min_output_x, min_output_y, max_output_x,         &
		  max_output_y, x_offset, y_offset


	character(len=255), parameter :: output_dem = "reduced_dem.xyz", sif_parameter_file = "sif_parameters.txt"

	! first, read in the DEM. the DEM file is the first command line argument

	call READ_DEM

	! second, read in the contour file, and find the minimum and maximum x and y points

	min_contour_x = 9.e10
	min_contour_y = 9.e10
	max_contour_x = -9.e10
	max_contour_y = -9.e10


	call getarg(2,filename_contour)

	open(unit=30, file=filename_contour, access="sequential", form="formatted", status="old")

	read_contour: do

		read(30,*,iostat=istat) x, y
		if(istat /=0) THEN
			exit read_contour
		end if

		if(x > max_contour_x) THEN
			max_contour_x = x
		endif

		if(x < min_contour_x) THEN
			min_contour_x = x
		endif

		if(y > max_contour_y) THEN
			max_contour_y = y
		endif

		if(y < min_contour_y) THEN
			min_contour_y = y
		endif


	end do read_contour

	close(unit=30)

	x_offset = xmin - real(floor((xmin / dx)))*dx
	y_offset = ymin - real(floor((ymin / dy)))*dy


	! find the range of interest

	min_output_x = floor((min_contour_x-x_offset) / real(dx)) * dx + x_offset
	min_output_y = floor((min_contour_y-y_offset) / real(dy)) * dy + y_offset
	max_output_x = ceiling((max_contour_x-x_offset) / real(dx)) * dx + x_offset
	max_output_y = ceiling((max_contour_y-y_offset) / real(dy)) * dy + y_offset


!	write(6,*) min_contour_x, min_output_x, int(min_contour_x / real(dx)), dx
!	write(6,*) min_contour_y, min_output_y, int(min_contour_y / real(dy)), dy

!	stop

	open(unit=50, file=output_dem, access="sequential", form="formatted", status="replace")

	start_x_counter = nint((min_output_x - xmin)/dx) + 1
	start_y_counter = nint((min_output_y - ymin)/dy) + 1
	end_x_counter = nint((max_output_x - xmin)/dx) + 1
	end_y_counter = nint((max_output_y - ymin)/dy) + 1


!	write(6,*) min_output_x,max_output_x, xmin, dx
!	write(6,*) min_output_y, max_output_y, ymin, dy

!	write(6,*) start_x_counter, end_x_counter, nx
!	write(6,*) start_y_counter, end_y_counter, ny

!	stop

	do x_counter = start_x_counter, end_x_counter, 1
		do y_counter = start_y_counter, end_y_counter, 1

			if(x_counter <= 0 .or. y_counter <= 0 .or. x_counter > nx .or. y_counter > ny ) THEN
				write(50,*) real(x_counter-1)*dx + xmin, real(y_counter-1)*dy+ymin, noval
			else

				write(50,*) real(x_counter-1)*dx + xmin, real(y_counter-1)*dy+ymin, DEM(x_counter,y_counter)

			endif

		end do
	end do

	close(unit=50)
	!write out the parameters that need to go into the SIF file

	open(unit=60, file=sif_parameter_file, access="sequential", form="formatted", status="replace")

	write(60,*) "! If you have multiple DEMs to read into your solver, you will have to change some of the parameters listed below."
	write(60,*) "! The primary thing will be the variable number. Secondly, you will have to change the first line "
	write(60,*) "! with the name of the DEM. Here, I have just made the generic name 'outDEM'"
	write(60,*) ""
	write(60,*) 'Variable 1 = String "outDEM" ! change this variable when you put it in a SIF file'
	write(60,*) 'Variable 1 data file = File "' // trim(adjustl(output_dem)) // '" ! adjust the path if needed later'
	write(60,*) "Variable 1 x0 = Real ", min_output_x
	write(60,*) "Variable 1 y0 = Real ", min_output_y 
	write(60,*) "Variable 1 lx = Real ", max_output_x - min_output_x 
	write(60,*) "Variable 1 ly = Real ", max_output_y - min_output_y 
	write(60,*) "Variable 1 Nx = Integer ", end_x_counter - start_x_counter+ 1
	write(60,*) "Variable 1 Ny = Integer",  end_y_counter - start_y_counter + 1
	write(60,*) "Variable 1 Invert = Logical True"
	write(60,*) "Variable 1 Fill = Logical False "
	write(60,*) "Variable 1 Position Tol = Real 1.0e-1"
	write(60,*) "Variable 1 No Data = Real ", noval
	write(60,*) "Variable 1 No Data Tol = Real 1.0 "

	close(unit=60)



end program reduce_dem
