!-------------------------------------------------------------------------------
! "Grid_Interpolator v.2.0" (Grid manager tool)
! Copyright 2016-2021 (RSE SpA)
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
! Program unit: main_module
! Description: Fortran main module of Grid Interpolator
!-------------------------------------------------------------------------------
module main_module
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
! Number of output nodes (points) along x/y/z axis
integer(4) :: nx_out,ny_out,nz_out
! Error (IO) unit
integer(4) :: uerr = 0
! Total number of input vertices (considering all the input files)
integer(4) :: n_points_in
integer(4) :: distance_exponent,n_parts_out_x,n_parts_out_y,n_parts_out_z
double precision :: x_min_out,y_min_out,z_min_out,x_max_out,y_max_out,z_max_out
double precision :: dx_out,dy_out,dz_out,normalized_threshold_pos
double precision :: normalized_threshold_neg,normalized_influence_radius
double precision :: abs_mean_latitude,missing_data_value,delta_x,delta_y
double precision :: threshold_pos,threshold_neg,x_lon_trans,y_lat_trans
double precision :: delta_lat,delta_lon,dry_input_value
! Arrays of the edge node indices for the output files
integer(4),dimension(:),allocatable :: ix_out_min_file,ix_out_max_file
integer(4),dimension(:),allocatable :: iy_out_min_file,iy_out_max_file
integer(4),dimension(:),allocatable :: iz_out_min_file,iz_out_max_file
! Array of the input grid
double precision,dimension(:,:),allocatable :: field_in
! Arrays of the output grids
double precision,dimension(:,:,:,:),allocatable :: field_out
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
!------------------------
! Deallocations
!------------------------
end module main_module
