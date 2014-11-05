program read_dem_austfonna

! this program converts the DEM files in the Austfonna format to a format that can be used in the programs set up
! to create an extruded mesh, as per the Wiki (https://github.com/JosefinAhlkrona/FROZEN/wiki/Create-Grids)

! run using:
! > read_dem_austfonna input_dem output_file

      use DEMVar
	implicit none

	character(len=255) infile, outfile

	real :: x1, y1, z1, x2, y2, z2

	integer :: x_counter, y_counter, x_start, y_start, x_end, y_end, x_increment, y_increment, istat



	logical :: found_dx

	noval = -9999999.

	call getarg(1,infile)
	call getarg(2,outfile)


	open(unit=10, file=infile, access="sequential", form="formatted", status="old")

	! determine dx and dy

	read(10,*) x1, y1, z1

	xmin = x1

	found_dx = .false.

	nx = 0

	determine_dx_dy: do

		! the files are structure so that it is in rows of x

		read(10,*) x2, y2, z2

		if(.not. found_dx) THEN
			dx = x2 - x1
			found_dx = .true.
		endif

		if(y2 /=y1) THEN

			dy = y2-y1

	

			nx = nx+1

			if(x2 < xmin) THEN
				xmin = x2
			endif


			exit determine_dx_dy

		endif

		! assuming this is a perfectly square grid...
		nx = nx+1

		! might as well find xmin as well

		if(x2 < xmin) THEN
			xmin = x2
		endif


	end do determine_dx_dy


	rewind (unit=10)


	! find the number of y

	read(10,*) x1, y1, z1
	ny = 1
	ymin = y1

	find_ny: do

		read(10,*, iostat=istat) x2, y2, z2
		if(istat /=0) THEN
			exit find_ny
		endif

		if(y2 /= y1) THEN
			ny = ny+1
			y1 = y2

			if(y2 < ymin) THEN
				ymin=y2
			endif

		endif


	end do find_ny

	! now, read in the entire dem to an array
	!
	! in one of the files, it appears that x is decreasing. Doing a generalization so that this is not a problem


	if(dx > 0) THEN

		x_start = 1
		x_end = nx
		x_increment = 1

	else

		x_start = nx
		x_end = 1
		x_increment = -1

	endif

	if(dy > 0) THEN

		y_start = 1
		y_end = ny
		y_increment = 1

	else

		y_start = ny
		y_end = 1
		y_increment = -1

	endif

	allocate(dem(nx,ny))

	rewind(unit=10)

	read_dem: do y_counter = y_start, y_end, y_increment
		do x_counter = x_start, x_end, x_increment

			read(10,*, iostat=istat) x1, y1, z1
			if(istat /=0) THEN
				write(6,*) "there is something wrong with the file structure in the dem"
				write(6,*) "aborting..."
				stop
			endif

			dem(x_counter,y_counter) = z1
	
		end do
	
	end do read_dem

	close(unit=10)

	! now write out the new dem file

	open(unit=40, file=outfile, access="sequential", form="formatted", status="replace")

	write(40,*) nx, ny, xmin, ymin, abs(dx), abs(dy), noval

	do x_counter = 1, nx

		do y_counter = 1, ny

			write(40,*) dem(x_counter,y_counter)

		end do
	end do

	deallocate(dem)

	close(unit=40)


end program read_dem_austfonna
