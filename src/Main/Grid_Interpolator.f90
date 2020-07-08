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
! Description. “Grid_Interpolator v.2.0” (RSE SpA) reads a 3D field of values 
!              from an input grid (file format: 1 comment line plus records on 4 
!              columns - x, y, z, variable-) and interpolates it on an output 
!              grid with a different spatial resolution. The output field is 
!              available in both the file formats ".csv" and ".dem". Despiking. 
!              This tool is also useful to post-process the 2D fields of the 
!              maximum specific height and the maximum water depth as estimated 
!              by the CFD-SPH code SPHERA v.9.0.0 (RSE SpA; available on 
!              github.com).
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
integer(4) :: read_stat,open_stat,alloc_stat,distance_exponent
double precision :: denom,dx_out,dy_out,dz_out,distance,x_min_out,y_min_out
double precision :: y_max_out,z_max_out,threshold_pos,threshold_neg,mean,sigma
double precision :: normalized_threshold_pos,normalized_threshold_neg,z_min_out
double precision :: x_max_out,normalized_influence_radius
double precision :: abs_mean_latitude,lam_min,phi_min,L_x,L_y
double precision,dimension(:,:),allocatable :: field_in,field_out
double precision,dimension(:,:),allocatable :: field_out_lon_lat
character(100) :: input_grid_file_name
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
mean = 0.d0
sigma = 0.d0
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
open(11,file='Grid_Interpolator.inp')
read(11,*)
read(11,*) n_points_in,x_min_out,y_min_out,z_min_out,x_max_out,y_max_out,      &
   z_max_out,dx_out,dy_out,dz_out
read(11,*)
read(11,*,IOSTAT=read_stat) input_grid_file_name
if (read_stat/=0) then
   write(0,*) "Error in reading Grid_Interpolator.inp. The program stops here. "
   stop
endif
read(11,*)
read(11,*) normalized_threshold_pos,normalized_threshold_neg,                  &
   normalized_influence_radius,distance_exponent
read(11,*)
read(11,*) lam_min,phi_min,abs_mean_latitude
close(11)
nx_out = int((x_max_out - x_min_out) / dx_out)
ny_out = int((y_max_out - y_min_out) / dy_out)
nz_out = int((z_max_out - z_min_out) / dz_out)
n_points_out = nx_out * ny_out * nz_out
if(.not.allocated(field_in)) then
   allocate(field_in(n_points_in,4),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(*,*) "Allocation of field_in failed; the execution terminates here."
      stop
   endif
endif
if(.not.allocated(field_out)) then
   allocate(field_out(n_points_out,4),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(*,*) "Allocation of field_out failed; the execution terminates ",  &
         "here."
      stop
   endif
endif
if(.not.allocated(field_out_lon_lat)) then
   allocate(field_out_lon_lat(n_points_out,2),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(*,*) "Allocation of field_out_lon_lat failed; the execution ",     &
         "terminates here."
      stop
   endif
endif
field_in = -999.d0
field_out = 0.d0
open(12,file=trim(input_grid_file_name),IOSTAT=open_stat)
if (open_stat/=0) then
   write(0,*) "Error in opening .txt input file. The program stops here. "
   stop
endif
read(12,*)
do i=1,n_points_in
  read (12,*) field_in(i,1),field_in(i,2),field_in(i,3),field_in(i,4)
  mean = mean + field_in(i,4)
  sigma = sigma + field_in(i,4)**2
enddo
close(12)
if (n_points_in>0) then
   mean = mean / n_points_in
   sigma = sigma / n_points_in
   sigma = dsqrt(sigma - mean**2)
   threshold_pos = mean + normalized_threshold_pos * sigma
   threshold_neg = mean - normalized_threshold_neg * sigma
endif
write(*,*) "mean: ",mean
write(*,*) "standard deviation: ",sigma
write(*,*) "positive threshold: ",threshold_pos
write(*,*) "negative threshold: ",threshold_neg
write(*,*) "End Reading the origin field"
! End 1)
! 2) Interpolation (inverse of distance**n)
write(*,*) "2) Interpolation "
i_record = 0
do k=1,nz_out
   do j=1,ny_out
      do i=1,nx_out
         i_record = i_record + 1
         field_out(i_record,1) = x_min_out + dx_out * (i - 0.5)
         field_out(i_record,2) = y_min_out + dy_out * (j - 0.5)
         field_out(i_record,3) = z_min_out + dz_out * (k - 0.5)
      enddo
   enddo
enddo
!$omp parallel do default(none)                                                &
!$omp shared(n_points_out,n_points_in,field_in,field_out,dx_out,threshold_pos) &
!$omp shared(threshold_neg)                                                    &
!$omp private(j,denom,distance,i)
do j=1,n_points_out
   denom = 0.d0
   do i=1,n_points_in
      if ((field_in(i,4)<=threshold_pos).and.(field_in(i,4)>=threshold_neg))   &
         then
         distance = dsqrt((field_in(i,1) - field_out(j,1)) ** 2 + (            &
                    field_in(i,2) - field_out(j,2)) ** 2 +(field_in(i,3) -     &
                    field_out(j,3)) ** 2)
         if (distance<(1.d-6*dx_out)) then
            field_out(j,4) = field_in(i,4)
            denom = 0.d0
            exit
            else
               if ((distance<=(normalized_influence_radius*dx_out*dsqrt(3.d0)))&
                  .and.(field_in(i,4)/=-999.d0)) then 
                  field_out(j,4) = field_out(j,4) + field_in(i,4) / distance   &
                                   ** abs(distance_exponent)
                  denom = denom + 1.d0 / distance ** abs(distance_exponent)
               endif
         endif
         else
            cycle
      endif
   enddo
   if (denom/=0.d0) field_out(j,4) = field_out(j,4) / denom
enddo
!$omp end parallel do
write(*,*) "End Interpolation "
! Grid conversion from DEM2xyz v.2.0 (RSE SpA)”:  
! (X,Y) in (m) to (lon,lat) in (°)
abs_mean_latitude = abs_mean_latitude / 180.d0 * 3.1415926
L_x = 111412.84d0 * dcos(abs_mean_latitude) - 93.5d0 * dcos(3.d0 *             &
        abs_mean_latitude) + 0.118d0 * dcos(5.d0 * abs_mean_latitude)
L_y = (111132.92d0 - 559.82d0 * dcos(2.d0 * abs_mean_latitude) +               &
        1.175d0 * dcos(4.d0 * abs_mean_latitude) - 0.0023d0 * dcos(6.d0 *      &
        abs_mean_latitude))
field_out_lon_lat(:,1) = field_out(:,1)/L_x + lam_min
field_out_lon_lat(:,2) = field_out(:,2)/L_y + phi_min
! End 2)
! 3) Writing the output field
write(*,*) "3)  Writing the output field "
open(13,file='output_field.csv')
write(13,'(4a)') "x;","y;","z;","variable"
do i=1,n_points_out
   write(13,'(3(e15.6,a),e15.6)') field_out(i,1),";",field_out(i,2),";",       &
      field_out(i,4),";",field_out(i,4) 
enddo
close(13)
open (14,file='output_field_lon_lat.csv')
write(14,'(4a)') "lam;","phi;","z;","variable"
do i=1,n_points_out
   write(14,'(2(f15.12,a),e15.6,a,e15.6)') field_out_lon_lat(i,1),";",         &
      field_out_lon_lat(i,2),";", field_out(i,4),";",field_out(i,4) 
enddo
close(14)
open(12,file='output_field.dem')
write(12,'(a,i15)') "ncols ",nx_out
write(12,'(a,i15)') "nrows ",ny_out
write(12,'(a,g15.5)') "xllcorner ",x_min_out
write(12,'(a,g15.5)') "yllcorner ",y_min_out
write(12,'(a,g15.5)') "cellsize ",dx_out
write(12,'(2a)') "NODATA_value ","-999."  
do i=ny_out,1,-1
   do j=1,nx_out
      k = (i - 1) * nx_out + j 
      if (field_out(k,4)==0.d0) field_out(k,4) = -999.d0
      write(12,'(1x,g15.5)',ADVANCE='NO') field_out(k,4) 
   enddo
   write(12,*)
enddo
close(12)
write(*,*) "End Writing the output field "
! End 3)
!------------------------
! Deallocations
!------------------------
if(allocated(field_in)) then
   deallocate(field_in,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(*,*) "Deallocation of field_in failed; the execution terminates ", &
         "here."
      stop
   endif
endif
if(allocated(field_out)) then
   deallocate(field_out,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(*,*) "Deallocation of field_out failed; the execution terminates ",&
         "here."
      stop
   endif
endif
if(allocated(field_out_lon_lat)) then
   deallocate(field_out_lon_lat,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(*,*) "Deallocation of field_out_lon_lat failed; the execution    ",&
         "terminates here."
      stop
   endif
endif
write(*,*) "***  Grid Interpolator has terminated. "
end program Grid_Interpolator
