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
! Program unit: output_grid_partition
! Description: Partition of the output grid in "n_parts_out_x" * 
!              "n_parts_out_y" * "n_parts_out_z" output files (for each type of 
!              output file)
!-------------------------------------------------------------------------------
subroutine output_grid_partition
!------------------------
! Modules
!------------------------
use main_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: iii,jjj,kkk,i_file_out
character(100) :: array_name
! Arrays of the node indices for the partition of the output grid
integer(4),dimension(:),allocatable :: ix_out,iy_out,iz_out
! Arrays of the file indices for the map of the output files
integer(4),dimension(:),allocatable :: ix_file_map,iy_file_map,iz_file_map
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine allocate_de_int4_r1(allocate_flag,array,extent_1,uerr,array_name)
      implicit none
      integer(4),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocate_flag
      integer(4),intent(in) :: extent_1
      integer(4),intent(in) :: uerr
      character(100),intent(in) :: array_name
   end subroutine allocate_de_int4_r1
end interface
!------------------------
! Allocations
!------------------------
array_name = "ix_out"
call allocate_de_int4_r1(.true.,ix_out,n_parts_out_x+1,uerr,array_name)
array_name = "iy_out"
call allocate_de_int4_r1(.true.,iy_out,n_parts_out_y+1,uerr,array_name)
array_name = "iz_out"
call allocate_de_int4_r1(.true.,iz_out,n_parts_out_z+1,uerr,array_name)
array_name = "ix_file_map"
call allocate_de_int4_r1(.true.,ix_file_map,                                   &
   n_parts_out_x*n_parts_out_y*n_parts_out_z,uerr,array_name)
array_name = "iy_file_map"
call allocate_de_int4_r1(.true.,iy_file_map,                                   &
   n_parts_out_x*n_parts_out_y*n_parts_out_z,uerr,array_name)
array_name = "iz_file_map"
call allocate_de_int4_r1(.true.,iz_file_map,                                   &
   n_parts_out_x*n_parts_out_y*n_parts_out_z,uerr,array_name)
array_name = "ix_out_min_file"
call allocate_de_int4_r1(.true.,ix_out_min_file,                               &
   n_parts_out_x*n_parts_out_y*n_parts_out_z,uerr,array_name)
array_name = "iy_out_min_file"
call allocate_de_int4_r1(.true.,iy_out_min_file,                               &
   n_parts_out_x*n_parts_out_y*n_parts_out_z,uerr,array_name)
array_name = "iz_out_min_file"
call allocate_de_int4_r1(.true.,iz_out_min_file,                               &
   n_parts_out_x*n_parts_out_y*n_parts_out_z,uerr,array_name)
array_name = "ix_out_max_file"
call allocate_de_int4_r1(.true.,ix_out_max_file,                               &
   n_parts_out_x*n_parts_out_y*n_parts_out_z,uerr,array_name)
array_name = "iy_out_max_file"
call allocate_de_int4_r1(.true.,iy_out_max_file,                               &
   n_parts_out_x*n_parts_out_y*n_parts_out_z,uerr,array_name)
array_name = "iz_out_max_file"
call allocate_de_int4_r1(.true.,iz_out_max_file,                               &
   n_parts_out_x*n_parts_out_y*n_parts_out_z,uerr,array_name)
!------------------------
! Initializations
!------------------------
ix_out(:) = 0
iy_out(:) = 0
iz_out(:) = 0
ix_file_map(:) = 0
iy_file_map(:) = 0
iz_file_map(:) = 0
ix_out_min_file(:) = 0
iy_out_min_file(:) = 0
iz_out_min_file(:) = 0
ix_out_max_file(:) = 0
iy_out_max_file(:) = 0
iz_out_max_file(:) = 0
!------------------------
! Statements
!------------------------
! Node indices for the partition of the output grid
do iii=1,n_parts_out_x
   ix_out(iii) = 1 + (iii - 1) * int(nx_out / n_parts_out_x)
enddo
ix_out(n_parts_out_x+1) = nx_out
do iii=1,n_parts_out_y
   iy_out(iii) = 1 + (iii - 1) * int(ny_out / n_parts_out_y)
enddo
iy_out(n_parts_out_y+1) = ny_out
do iii=1,n_parts_out_z
   iz_out(iii) = 1 + (iii - 1) * int(nz_out / n_parts_out_z)
enddo
iz_out(n_parts_out_z+1) = nz_out
! Output files are ordered as follows: x-values ascending, then y-values 
! descending, finally z-values ascending.
i_file_out = 1
do kkk=1,n_parts_out_z
   do jjj=n_parts_out_y,1,-1
      do iii=1,n_parts_out_x
         ix_file_map(i_file_out) = iii
         iy_file_map(i_file_out) = jjj
         iz_file_map(i_file_out) = kkk
         i_file_out = i_file_out + 1
      enddo
   enddo
enddo
! Minimum and maximum node indices for each file
do i_file_out=1,(n_parts_out_x*n_parts_out_y*n_parts_out_z)
   ix_out_min_file(i_file_out) = ix_out(ix_file_map(i_file_out))
   ix_out_max_file(i_file_out) = ix_out(ix_file_map(i_file_out)+1)
   iy_out_min_file(i_file_out) = iy_out(iy_file_map(i_file_out))
   iy_out_max_file(i_file_out) = iy_out(iy_file_map(i_file_out)+1)
   iz_out_min_file(i_file_out) = iz_out(iz_file_map(i_file_out))
   iz_out_max_file(i_file_out) = iz_out(iz_file_map(i_file_out)+1)
enddo
!------------------------
! Deallocations
!------------------------
array_name = "ix_out"
call allocate_de_int4_r1(.false.,ix_out,0,uerr,array_name)
array_name = "iy_out"
call allocate_de_int4_r1(.false.,iy_out,0,uerr,array_name)
array_name = "iz_out"
call allocate_de_int4_r1(.false.,iz_out,0,uerr,array_name)
array_name = "ix_file_map"
call allocate_de_int4_r1(.false.,ix_file_map,0,uerr,array_name)
array_name = "iy_file_map"
call allocate_de_int4_r1(.false.,iy_file_map,0,uerr,array_name)
array_name = "iz_file_map"
call allocate_de_int4_r1(.false.,iz_file_map,0,uerr,array_name)
return
end subroutine output_grid_partition
