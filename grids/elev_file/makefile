
FC = gfortran

bicubic.o: bicubic.f90
	$(FC) -o bicubic.o -O2 bicubic.f90

nearest_int: nearest_int.f90 bicubic.f90
	$(FC) -o nearest_int -O2 nearest_int.f90 bicubic.f90
