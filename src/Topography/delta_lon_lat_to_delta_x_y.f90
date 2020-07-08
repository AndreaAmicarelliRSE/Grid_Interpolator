!-------------------------------------------------------------------------------
! "Grid_Interpolator v.2.0" (DEM manager tool)
! Copyright 2016-2020 (RSE SpA)
! "Grid_Interpolator v.2.0" authors and email contact are provided on the 
! documentation file.
! This file is part of Grid_Interpolator v.2.0 .
! Grid_Interpolator v.2.0 is free software: you can redistribute it and/or 
! modify it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! Grid_Interpolator is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with Grid_Interpolator. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! This file was copied and pasted from DEM2xyz v.2.0 (RSE SpA), commit of 25 
! June 2020. The distribution of this file under the GNU-LGPL License is 
! authorized by the Copyright owner of DEM2xyz v.2.0 (RSE SpA).
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: delta_lon_lat_to_delta_x_y
! Description: Truncated Fourier's transforms for the x and y (cartographic 
!              conformal coordinates, in meters) intervals covered by the 
!              longitude and latitude input intervals. The conversion refers to 
!              the WGS84 ellipsoid (e.g., for the geodetic-cartographic system 
!              UTM-ED50). However, it can be also used for other 
!              geodetic-cartographic systems with larger approximations, 
!              provided that the input longitude and latitude intervals 
!              (in degrees) are small enough.        
!-------------------------------------------------------------------------------
subroutine delta_lon_lat_to_delta_x_y(delta_lon,delta_lat,abs_mean_latitude,   &
   delta_x,delta_y)
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(in) :: delta_lon,delta_lat,abs_mean_latitude
double precision,intent(out) :: delta_x,delta_y
double precision :: abs_mean_lat_rad
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
! Conversion degrees to radians
abs_mean_lat_rad = abs_mean_latitude / 180.d0 * 3.1415926
! Linear unit discretization: conversion (lon,lat) in degrees to (x,y) in meters
delta_x = delta_lon * (111412.84d0 * dcos(abs_mean_lat_rad) - 93.5d0 *         &
          dcos(3.d0 * abs_mean_lat_rad) + 0.118d0 * dcos(5.d0 *                &
          abs_mean_lat_rad))
delta_y = delta_lat * (111132.92d0 - 559.82d0 * dcos(2.d0 * abs_mean_lat_rad)  &
          + 1.175d0 * dcos(4.d0 * abs_mean_lat_rad) - 0.0023d0 * dcos(6.d0 *   &
          abs_mean_lat_rad))
!------------------------
! Deallocations
!------------------------
return
end subroutine delta_lon_lat_to_delta_x_y
