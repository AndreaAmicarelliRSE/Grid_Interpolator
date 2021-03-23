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
! Program unit: allocate_de_dp_r4
! Description: Allocation/deallocation of a generic allocatable array of type 
!              "double precision" and range (number of dimensions) 4
!-------------------------------------------------------------------------------
subroutine allocate_de_dp_r4(allocate_flag,array,extent_1,extent_2,extent_3, &
                          extent_4,uerr,array_name)
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
double precision,dimension(:,:,:,:),allocatable,intent(inout) :: array
logical,intent(in) :: allocate_flag
integer(4),intent(in) :: extent_1
integer(4),intent(in) :: extent_2
integer(4),intent(in) :: extent_3
integer(4),intent(in) :: extent_4
! error unit
integer(4),intent(in) :: uerr
character(100),intent(in) :: array_name
integer(4) :: alloc_stat
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
if((allocate_flag.eqv..true.).and.(.not.allocated(array))) then
   allocate(array(extent_1,extent_2,extent_3,extent_4),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) "Allocation of ",trim(array_name),                         &
         " failed; the execution stops here."
      stop
   endif
endif
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
!------------------------
! Deallocations
!------------------------
if((allocate_flag.eqv..false.).and.(allocated(array))) then
   deallocate(array,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) "Deallocation of ",trim(array_name),                       &
         " failed; the execution stops here."
      stop
   endif
endif
return
end subroutine allocate_de_dp_r4
