elmerf90 -o fausto_temp fausto_temp.f90
elmerf90 -o geo_coord geo_coord.f90
elmerf90 -o IceFlowProperties IceFlowProperties.f90
elmerf90 -o initVeloTemp initVeloTemp.f90
elmerf90 -o interpolateByLayers interpolateByLayers_2.f90
elmerf90 -o lateralTemp lateralTemp.f90
elmerf90 -o limitVelo limitVelo.f90
elmerf90 -o pdd pdd.f90
elmerf90 -o VolumeSolver VolumeSolver.f90
elmerf90 -o TemperateIce TemperateIce.f90
#elmerf90 -I/home/josefin/PnetCDF/include -L/home/josefin/PnetCDF/lib output_netcdf.f90 -o output_netcdf2.so -lpnetcdf
elmerf90 -I/usr/include -L/usr/lib GridDataReader.f90 -o GridDataReader1.so -lnetcdf -lnetcdff
elmerf90 -I/usr/include -L/usr/lib GridDataReader.f90 -o GridDataReader2.so -lnetcdf -lnetcdff
elmerf90 -I/usr/include -L/usr/lib GridDataReader.f90 -o GridDataReader3.so -lnetcdf -lnetcdff
elmerf90 -I/usr/include -L/usr/lib GridDataReader.f90 -o GridDataReader4.so -lnetcdf -lnetcdff
