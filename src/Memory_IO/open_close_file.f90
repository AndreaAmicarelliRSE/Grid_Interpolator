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
! Program unit: open_close_file
! Description: File opening or closing.
!-------------------------------------------------------------------------------
subroutine open_close_file(open_flag,file_unit,file_name,uerr)
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
character(100),intent(inout) :: file_name
logical,intent(in) :: open_flag
integer(4),intent(in) :: file_unit,uerr
integer(4) :: open_stat
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
if (open_flag.eqv..true.) then
   open(file_unit,file=trim(file_name),IOSTAT=open_stat)
   if (open_stat/=0) then
      write(uerr,*) "Error in opening the file ",trim(file_name),              &
         ". The program stops."
      stop
   endif
   else
      close(file_unit)
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine open_close_file
