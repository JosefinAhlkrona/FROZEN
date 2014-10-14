#! /bin/bash


# Longitude region of interest
long_min=0
long_max=40

# latitude region of interest

lat_min=74
lat_max=85


# Alberts equal area conic projection. It probably doesn't matter what projection you use, as long as it is consistent for every step

center_longitude=18
center_latitude=78.75

southern_standard_parallel=76
northern_standard_parallel=82

map_width=20c # likely arbitrary


outline_file=austfonna.csv

# the csv file is taken from the kml file, must convert to lat long only file, skip this step if it is a plain long-lat file

awk -F, '{print $1, $2 }' ${outline_file} > awk.out

# convert the longitude, latitude file to x-y using the given map projection

outline_xy="contour.dat"

mapproject awk.out  -R${long_min}/${long_max}/${lat_min}/${lat_max} -Jb${center_longitude}/${center_latitude}/${southern_standard_parallel}/${northern_standard_parallel}/${map_width} -C -Fe > ${outline_xy}

# create a geo file

geo_file="contour.geo"

./make_geo_file ${outline_xy} ${geo_file}

# create mesh file

mesh_file="contour.mesh"

./GeoToMesh ${geo_file} ${mesh_file}

