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
! Program unit: writing_output_data
! Description: 
!-------------------------------------------------------------------------------
subroutine writing_output_data
!------------------------
! Modules
!------------------------
use main_module
!------------------------
! Declarations
!------------------------
implicit none
! Generic node indices for "do" constructs
integer(4) :: ix,iy,iz
integer(4) :: i_file_out
double precision :: x_min_out_file,y_min_out_file,dx_out_aux,dy_out_aux,eps
character(100) :: output_file_name,array_name
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine allocate_de_dp_r4(allocate_flag,array,extent_1,extent_2,         &
                                extent_3,extent_4,uerr,array_name)
      implicit none
      double precision,dimension(:,:,:,:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocate_flag
      integer(4),intent(in) :: extent_1
      integer(4),intent(in) :: extent_2
      integer(4),intent(in) :: extent_3
      integer(4),intent(in) :: extent_4
      integer(4),intent(in) :: uerr
      character(100),intent(in) :: array_name
   end subroutine allocate_de_dp_r4
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
!------------------------
! Statements
!------------------------
! .csv output grid in both cartographic and geographic coordinates
!$omp parallel do default(none)                                                &
!$omp shared(n_parts_out_x,n_parts_out_y,n_parts_out_z,uerr,nz_out,ny_out)     &
!$omp shared(nx_out,ix_out_min_file,ix_out_max_file,iy_out_min_file)           &
!$omp shared(iy_out_max_file,iz_out_min_file,iz_out_max_file,field_out)        &
!$omp shared(abs_mean_latitude)                                                &
!$omp private(i_file_out,output_file_name,iz,iy,ix)
do i_file_out=1,(n_parts_out_x*n_parts_out_y*n_parts_out_z)
   write(output_file_name,'(a,i4.4,a)') 'output_field_',i_file_out,'.csv'
   call open_close_file(.true.,10+i_file_out,output_file_name,uerr)
   if (abs_mean_latitude>-1.d-9) then
      write(1010+i_file_out,'(4a)') "lam;","phi;","z;","variable"
      else
         write(10+i_file_out,'(4a)') "x;","y;","z;","variable"
   endif
   do iz=1,nz_out
      do iy=1,ny_out
         do ix=1,nx_out
            if ((ix>=ix_out_min_file(i_file_out)).and.                         &
                (ix<=ix_out_max_file(i_file_out)).and.                         &
                (iy>=iy_out_min_file(i_file_out)).and.                         &
                (iy<=iy_out_max_file(i_file_out)).and.                         &
                (iz>=iz_out_min_file(i_file_out)).and.                         &
                (iz<=iz_out_max_file(i_file_out))) then
               write(10+i_file_out,'(3(ES18.9,a),ES18.9)')                     &
                  field_out(ix,iy,iz,1),";",field_out(ix,iy,iz,2),";",         &
                  field_out(ix,iy,iz,4),";",field_out(ix,iy,iz,4)
            endif
         enddo
      enddo
   enddo
   call open_close_file(.false.,10+i_file_out,output_file_name,uerr)
enddo
!$omp end parallel do
! ".asc" output grid in cartographic coordinates (only for the bottom layer)
! The vertices of every ".asc" ouput file are ordered as follows: x-values 
! ascending, then y-values descending (single z-layer).
!$omp parallel do default(none)                                                &
!$omp shared(n_parts_out_x,n_parts_out_y,n_parts_out_z,uerr,ix_out_min_file)   &
!$omp shared(iy_out_min_file,x_min_out,dx_out,y_min_out,dy_out,field_out)      &
!$omp shared(missing_data_value,ny_out,nx_out,ix_out_max_file,iy_out_max_file) &
!$omp shared(abs_mean_latitude,x_lon_trans,y_lat_trans,delta_lon,delta_lat)    &
!$omp shared(delta_x,delta_y,dry_input_value)                                  &
!$omp private(i_file_out,output_file_name,x_min_out_file,y_min_out_file,ix)    &
!$omp private(iy,iz,dx_out_aux,dy_out_aux,eps)
do i_file_out=1,(n_parts_out_x*n_parts_out_y*n_parts_out_z)
   write(output_file_name,'(a,i4.4,a)') 'output_field_',i_file_out,'.asc'
   call open_close_file(.true.,10+i_file_out,output_file_name,uerr)
   write(10+i_file_out,'(a,i15)') "ncols ",                                    &
      (ix_out_max_file(i_file_out)-ix_out_min_file(i_file_out)+1)
   write(10+i_file_out,'(a,i15)') "nrows ",                                    &
      (iy_out_max_file(i_file_out)-iy_out_min_file(i_file_out)+1)
   x_min_out_file = x_min_out + (ix_out_min_file(i_file_out) - 1) * dx_out
   y_min_out_file = y_min_out + (iy_out_min_file(i_file_out) - 1) * dy_out
   dx_out_aux = dx_out
   dy_out_aux = dy_out
   if (abs_mean_latitude>-1.d-9) then
! Reference system conversion: cartographic to geographic: (X,Y) in (m) to 
! (lon,lat) in (°)
      x_min_out_file = delta_lon * x_min_out_file / delta_x
      y_min_out_file = delta_lat * y_min_out_file / delta_y
! Reference system conversion: cartographic to geographic: (X,Y) in (m) to 
! (lon,lat) in (°)
! In the future, dx_out and dy_out might depend on the output file. Thus, their 
! conversion is repeated each file
      dx_out_aux = delta_lon * dx_out_aux / delta_x
      dy_out_aux = delta_lat * dy_out_aux / delta_y
   endif
! Reference system conversion: local (SPHERA) to global (georeferenced)
   x_min_out_file = x_min_out_file + x_lon_trans
   y_min_out_file = y_min_out_file + y_lat_trans
   write(10+i_file_out,'(a,ES18.10)') "xllcorner ",x_min_out_file
   write(10+i_file_out,'(a,ES18.10)') "yllcorner ",y_min_out_file
   eps = 1.d-9 * (dy_out_aux ** 2) / dabs(dy_out_aux)
   if ((dx_out_aux>(dy_out_aux+eps)).or.(dx_out_aux<(dy_out_aux-eps))) then
! Rectangular cells
      write(10+i_file_out,'(2(a,g15.5))') "dx ",dx_out_aux," dy ", &
         dy_out_aux
      else
! Squared cells
         write(10+i_file_out,'(a,g15.5)') "cellsize ",dx_out_aux
   endif
   write(10+i_file_out,*) "NODATA_value ",missing_data_value
   iz = 1
! In the future, the variable "dry_input_value" might depend on the output 
! file. Thus, the assessment of the variable "eps_aux" is repeated each file.
   eps = 1.d-9 * (dry_input_value ** 2) / dabs(dry_input_value)
   do iy=iy_out_max_file(i_file_out),iy_out_min_file(i_file_out),-1
      do ix=ix_out_min_file(i_file_out),ix_out_max_file(i_file_out)
! In the output ".asc" file, the input value asociated with dry cells is 
! replaced by the missing data value, as requested by some GIS software. 
         if ((field_out(ix,iy,1,4)<(dry_input_value+eps)).and.                 &
            (field_out(ix,iy,1,4)>(dry_input_value-eps))) then
            field_out(ix,iy,1,4) = missing_data_value
         endif
         write(10+i_file_out,'(1x,g15.5)',ADVANCE='NO') field_out(ix,iy,1,4)
      enddo
      write(10+i_file_out,*)
   enddo
   call open_close_file(.false.,10+i_file_out,output_file_name,uerr)
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
array_name = "ix_out_min_file"
call allocate_de_int4_r1(.false.,ix_out_min_file,0,uerr,array_name)
array_name = "iy_out_min_file"
call allocate_de_int4_r1(.false.,iy_out_min_file,0,uerr,array_name)
array_name = "iz_out_min_file"
call allocate_de_int4_r1(.false.,iz_out_min_file,0,uerr,array_name)
array_name = "ix_out_max_file"
call allocate_de_int4_r1(.false.,ix_out_max_file,0,uerr,array_name)
array_name = "iy_out_max_file"
call allocate_de_int4_r1(.false.,iy_out_max_file,0,uerr,array_name)
array_name = "iz_out_max_file"
call allocate_de_int4_r1(.false.,iz_out_max_file,0,uerr,array_name)
array_name = "field_out"
call allocate_de_dp_r4(.false.,field_out,0,0,0,0,uerr,array_name)
return
end subroutine writing_output_data
