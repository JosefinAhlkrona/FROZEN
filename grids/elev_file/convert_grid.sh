#! /bin/bash

# Written by Evan Gowan

################################################################################
# MAKE SURE YOU RUN THIS IN BASH, NOT DASH (the default shell Ubuntu/Linux Mint)
################################################################################

# this script extracts a given area out of GEBCO and outputs an x-y grid centered at the given latitude and longitude

# get GEBCO here: http://www.gebco.net/


# Longitude region of interest
long_min=0
long_max=40

# latitude region of interest

lat_min=74
lat_max=85

# make a smaller grid, because calculating things with GEBCO requires a lot of memory

smaller_grid_file=smaller_grid_file.grd

grdcut gebco_08.nc -G${smaller_grid_file} -R${long_min}/${long_max}/${lat_min}/${lat_max}



# Alberts equal area conic projection. It probably doesn't matter what projection you use

center_longitude=18
center_latitude=78.75

southern_standard_parallel=76
northern_standard_parallel=82

map_width=20c # likely arbitrary

# output grid is in metres (change the -A option to -Ak for kilometers)

area_grid=area.grd

# note that the grid spacing must be an integer!

resolution=500

x_resolution=${resolution}
y_resolution=${resolution}


grdproject ${smaller_grid_file}  -R${long_min}/${long_max}/${lat_min}/${lat_max} -Jb${center_longitude}/${center_latitude}/${southern_standard_parallel}/${northern_standard_parallel}/${map_width} -G${area_grid} -D${x_resolution}/${y_resolution} -Ae -C -V # the -C option means that the grid is relative to the defined projection center (i.e. that is where the 0,0 coordinate is)


# the problem with grdproject is that the dx and dy are not exactly correction, probably due to a precision error. It also
# doesn't center things on the origin. I've created a program that will fix these issues

# get the x and y values. Note that this will not work in DASH shell (which is the default shell in Ubuntu and Linux Mint)

num_x=$( ncdump -h ${area_grid} | grep $'\tx = ' | sed -e 's/ //g' | sed 's/.*=//' | sed 's/;//g' )
num_y=$( ncdump -h ${area_grid} | grep $'\ty = ' | sed -e 's/ //g' | sed 's/.*=//' | sed 's/;//g' )

echo $num_x $num_y

ncdump  -v x area.grd  | awk '/data:/{y=1;next}y' | sed -e 's/[x={}; ]//g' | tr -d '\n' | sed -e 's/,/\n/g' > x_values.txt

ncdump  -v y area.grd  | awk '/data:/{y=1;next}y' | sed -e 's/[y={}; ]//g' | tr -d '\n' | sed -e 's/,/\n/g' > y_values.txt



# output to an ASCII file (warning, might be huge)

xyz_file="out.xyz"

# the -Cf option outputs it into array indices for reading into a Fortran program

grd2xyz ${area_grid} -V -sa -Cf > ${xyz_file}

# final step, since things are not in a nice grid starting at 0,0, this program converts it. You can likely change the resolution here if you want.

even_grid_file="even_grid.dat"

./nearest_int ${num_x} ${num_y} ${resolution} ${xyz_file} ${even_grid_file}







