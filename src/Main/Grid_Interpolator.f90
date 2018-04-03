!-------------------------------------------------------------------------------
! "Grid_Interpolator v.1.0" 
! Copyright 2016 (RSE SpA)
! "Grid_Interpolator v.1.0" authors and email contact are provided on 
! the documentation file.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Description. “Grid_Interpolator v.1.0” (RSE SpA) reads a 3D field of values 
!              from an input grid and interpolates them on an output grid with 
!              a different spatial resolution. The input file is a xyz file 
!              (with two additional ad-hoc lines at the beginning). The output 
!              field is available in both the file formats xyz and DEM. This 
!              tool is also useful to post-process the 2D fields of the maximum 
!              specific height and the maximum water depth as estimated by 
!              the SPH code SPHERA v.8.0 (RSE SpA; available on GitHub).
!-------------------------------------------------------------------------------
program Grid_Interpolator
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: i,j,k,n_points_in,n_points_out,i_record,nx_out,ny_out,nz_out
double precision :: denom,dx_out,dy_out,dz_out,distance,x_min,y_min,z_min,x_max
double precision :: y_max,z_max
double precision,dimension(:,:),allocatable :: field_in,field_out
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
!Subroutine info
write(*,*) "***Grid Interpolator is running:  "
write(*,*) "   This tool reads a 3D field of values from an origin grid"
write(*,*) "      and interpolates them on a destination grid"
write(*,*) "   Algorithm:"
write(*,*) "   1) Reading the input field"
write(*,*) "   2) Interpolation"
write(*,*) "   3) Writing the output field"
! 1) Reading the input field
write(*,*) "1)  Reading the input field"
open(1,file='input_field.prn')
read(1,*) 
read(1,'(i12,9(g12.5))') n_points_in,x_min,y_min,z_min,x_max,y_max,z_max,      &
   dx_out,dy_out,dz_out
read(1,*)
nx_out = int((x_max - x_min) / dx_out)
ny_out = int((y_max - y_min) / dy_out)
nz_out = int((z_max - z_min) / dz_out)
n_points_out = nx_out * ny_out * nz_out
allocate(field_in(n_points_in,4)) 
allocate(field_out(n_points_out,4)) 
field_in = -999.d0
field_out = 0.d0
! Reading 
do i=1,n_points_in
  read (1,'(4(g12.5))') field_in(i,1),field_in(i,2),field_in(i,3),field_in(i,4)
enddo
close(1)
write(*,*),"End Reading the origin field"
! End 1)
! 2) Interpolation (inverse of distance**2)
write(*,*) "2) Interpolation "
i_record = 0
do k=1,nz_out
   do j=1,ny_out
      do i=1,nx_out
         i_record = i_record + 1
         field_out(i_record,1) = x_min + dx_out * (i - 0.5)
         field_out(i_record,2) = y_min + dy_out * (j - 0.5)
         field_out(i_record,3) = z_min + dz_out * (k - 0.5)
      enddo
   enddo
enddo
!$omp parallel do default(none)                                                &
!$omp shared(n_points_out,n_points_in,field_in,field_out,dx_out)               &
!$omp private(j,denom,distance,i) 
do j=1,n_points_out
   denom = 0.d0
   do i=1,n_points_in
      distance = dsqrt((field_in(i,1) - field_out(j,1)) ** 2 + (field_in(i,2)  &
                 - field_out(j,2)) ** 2 +(field_in(i,3) - field_out(j,3)) ** 2)
      if ((distance<=(dx_out*10.d0)).and.(field_in(i,4)/=-999.d0)) then 
         field_out(j,4) = field_out(j,4) + field_in(i,4) / distance ** 2
         denom = denom + 1.d0 / distance ** 2
      endif
   enddo
   if (denom/=0.d0) field_out(j,4) = field_out(j,4) / denom
enddo
!$omp end parallel do 
! End 2)
! 3) Writing the output field (both file formats: xyz -colour=z- and DEM)
write(*,*) "3)  Writing the output field "
open (3,file='output_field_xyz.txt')
write(3,'(4a)') "x;","y;","z;","colour"
do i=1,n_points_out
   write(3,'(3(e12.6,a),e12.6)') field_out(i,1),";",field_out(i,2),";",       &
      field_out(i,4),";",field_out(i,4) 
enddo
close(3)
open(2,file='output_field_DEM.txt')
write(2,'(a,i15)') "ncols ",nx_out
write(2,'(a,i15)') "nrows ",ny_out
write(2,'(a,g15.5)') "xllcorner ",x_min
write(2,'(a,g15.5)') "yllcorner ",y_min
write(2,'(a,g15.5)') "cellsize ",dx_out
write(2,'(2a)') "NODATA_value ","-999."  
do i=ny_out,1,-1
   do j=1,nx_out
      k = (i - 1) * nx_out + j 
      if (field_out(k,4)==0.d0) field_out(k,4) = -999.d0
      write(2,'(1x,g12.5)',ADVANCE='NO') field_out(k,4) 
   enddo
   write(2,*)
enddo
close(2)
write(*,*) "End Writing the output field "
! End 3)
!------------------------
! Deallocations
!------------------------
deallocate(field_in)
deallocate(field_out)
write(*,*) "***  Grid Interpolator has terminated "
end program Grid_Interpolator
