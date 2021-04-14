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
! Program unit: reading_input_data
! Description: 
!-------------------------------------------------------------------------------
subroutine reading_input_data
!------------------------
! Modules
!------------------------
use main_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: n_input_grid_files,read_stat,i_igf,iii,i_record,n_missing_data
double precision :: mean,sigma,eps
character(100) :: array_name,file_name
integer(4),dimension(:),allocatable :: n_points_in_file
character(100),dimension(:),allocatable :: input_grid_file_name
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine allocate_de_ch100_r1(allocate_flag,array,extent_1,uerr,array_name)
      implicit none
      character(100),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocate_flag
      integer(4),intent(in) :: extent_1
      integer(4),intent(in) :: uerr
      character(100),intent(in) :: array_name
   end subroutine allocate_de_ch100_r1
   subroutine allocate_de_dp_r2(allocate_flag,array,extent_1,extent_2,uerr,    &
                             array_name)
      implicit none
      double precision,dimension(:,:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocate_flag
      integer(4),intent(in) :: extent_1
      integer(4),intent(in) :: extent_2
      integer(4),intent(in) :: uerr
      character(100),intent(in) :: array_name
   end subroutine allocate_de_dp_r2
   subroutine allocate_de_int4_r1(allocate_flag,array,extent_1,uerr,array_name)
      implicit none
      integer(4),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocate_flag
      integer(4),intent(in) :: extent_1
      integer(4),intent(in) :: uerr
      character(100),intent(in) :: array_name
   end subroutine allocate_de_int4_r1
   subroutine open_close_file(open_flag,file_unit,file_name,uerr)
      implicit none
      character(100),intent(inout) :: file_name
      logical,intent(in) :: open_flag
      integer(4),intent(in) :: file_unit,uerr
   end subroutine open_close_file
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
n_points_in = 0
n_missing_data = 0
mean = 0.d0
sigma = 0.d0
!------------------------
! Statements
!------------------------
file_name="Grid_Interpolator.inp"
call open_close_file(.true.,11,file_name,uerr)
read(11,*)
read(11,*) x_min_out,y_min_out,z_min_out,x_max_out,y_max_out,z_max_out,dx_out, &
   dy_out,dz_out
read(11,*)
read(11,*) normalized_threshold_pos,normalized_threshold_neg,                  &
   normalized_influence_radius,distance_exponent
read(11,*)
read(11,*) abs_mean_latitude
read(11,*)
read(11,*) missing_data_value,dry_input_value
read(11,*)
read(11,*) n_parts_out_x,n_parts_out_y,n_parts_out_z
read(11,*)
read(11,*) n_input_grid_files
array_name = "input_grid_file_name"
call allocate_de_ch100_r1(.true.,input_grid_file_name,n_input_grid_files,uerr, &
   array_name)
array_name = "n_points_in_file"
call allocate_de_int4_r1(.true.,n_points_in_file,n_input_grid_files,uerr,      &
   array_name)
do i_igf=1,n_input_grid_files
   read(11,*)
   read(11,*,IOSTAT=read_stat) input_grid_file_name(i_igf),                    &
      n_points_in_file(i_igf)
   if (read_stat/=0) then
      write(uerr,*) "Error in reading Grid_Interpolator.inp. The program stops."
      stop
   endif
   n_points_in = n_points_in + n_points_in_file(i_igf)
enddo
read(11,*)
read(11,*) x_lon_trans,y_lat_trans
file_name = "Grid_Interpolator.inp"
call open_close_file(.false.,11,file_name,uerr)
nx_out = int((x_max_out - x_min_out) / dx_out)
ny_out = int((y_max_out - y_min_out) / dy_out)
nz_out = int((z_max_out - z_min_out) / dz_out)
if ((n_parts_out_x>nx_out).or.(n_parts_out_y>ny_out).or.(n_parts_out_z>nz_out))&
   write(uerr,*) "The partition of the output grid is not correct. Please, ",  &
   "revise it. The execution stops here."
array_name = "field_in"
call allocate_de_dp_r2(.true.,field_in,n_points_in,4,uerr,array_name)
field_in(:,:) = -999.d0
i_record = 1
do i_igf=1,n_input_grid_files
   call open_close_file(.true.,12,input_grid_file_name(i_igf),uerr)
   read(12,*)
   do iii=1,n_points_in_file(i_igf)
     read(12,*) field_in(i_record,1),field_in(i_record,2),field_in(i_record,3),&
        field_in(i_record,4)
     eps = 1.d-9 * (-999.d0 ** 2) / dabs(-999.d0)
     if ((field_in(i_record,4)>(-999.d0+eps)).or.                              &
         (field_in(i_record,4)<(-999.d0-eps))) then
        mean = mean + field_in(i_record,4)
        sigma = sigma + field_in(i_record,4) ** 2
        else
           n_missing_data = n_missing_data + 1
     endif
     i_record = i_record + 1
   enddo
   call open_close_file(.false.,12,input_grid_file_name(i_igf),uerr)
enddo
if ((n_points_in-n_missing_data)>0) then
   mean = mean / (n_points_in - n_missing_data)
   sigma = sigma / (n_points_in - n_missing_data)
   sigma = dsqrt(sigma - mean ** 2)
   threshold_pos = mean + normalized_threshold_pos * sigma
   threshold_neg = mean - normalized_threshold_neg * sigma
endif
write(*,*) "      mean (input grid): ",mean
write(*,*) "      standard deviation (input grid): ",sigma
write(*,*) "      positive threshold (input grid): ",threshold_pos
write(*,*) "      negative threshold (input grid): ",threshold_neg
!------------------------
! Deallocations
!------------------------
array_name = "input_grid_file_name"
call allocate_de_ch100_r1(.false.,input_grid_file_name,0,uerr,array_name)
array_name = "n_points_in_file"
call allocate_de_int4_r1(.false.,n_points_in_file,0,uerr,array_name)
return
end subroutine reading_input_data
