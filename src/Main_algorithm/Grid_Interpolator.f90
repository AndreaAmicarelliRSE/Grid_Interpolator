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
! Description. “Grid_Interpolator v.2.0” (RSE SpA) reads a 3D field of values 
!              from one or more input grid files (file format: 1 comment line 
!              plus records on 4 columns - x, y, z, variable-) and interpolates 
!              it on an output grid with any spatial resolution. The output 
!              field is available in both the file formats ".csv" and ".asc" 
!              and can be partitioned in multiple output files. Despiking is 
!              available. Output coordinates are also expressed as 
!              geographical coordinates (longitude, latitude) and are 
!              possibly translated on a global reference system.
!-------------------------------------------------------------------------------
program Grid_Interpolator
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
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
write(*,*) "Grid Interpolator: start"
write(*,*) "   Reading input data: start"
call reading_input_data
write(*,*) "   Reading input data: end"
write(*,*) "   Interpolation: start"
call interpolation
write(*,*) "   Interpolation: end"
write(*,*) "   Output grid partition: start"
call output_grid_partition
write(*,*) "   Output grid partition: end"
write(*,*) "   Writing output data: start"
call writing_output_data
write(*,*) "   Writing output data: end"
write(*,*) "Grid Interpolator: end"
!------------------------
! Deallocations
!------------------------
end program Grid_Interpolator
