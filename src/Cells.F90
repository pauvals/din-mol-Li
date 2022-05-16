! Copyright (c) 2020  Sergio Alexis Paz
!
!  This file is part of GEMS. GEMS is an Extensible Molecular Simulator.
!   .
!  GEMS is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!  .
!  GEMS is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  .
!  You should have received a copy of the GNU General Public License
!  along with GEMS.  If not, see <https://www.gnu.org/licenses/>.


module gems_cells
use gems_program_types, only: boxed, box, mic, n1cells
use gems_groups,        only: vdistance, sys, igroup, atom, atom_dclist, useghost
use gems_constants,     only: dp,cdm,dm,ui_ev
use gems_errors

implicit none

! From 1 to 13 are the upper cells. From 14 to 26 the lower. 0 is the center
integer,parameter :: map(3,0:26) = &
         reshape( [ 0, 0, 0,  &
                    1, 0, 0,   1, 1, 0,   0, 1, 0,  -1, 1, 0,&
                    1, 0,-1,   1, 1,-1,   0, 1,-1,  -1, 1,-1,&
                    1, 0, 1,   1, 1, 1,   0, 1, 1,  -1, 1, 1,&
                    0, 0, 1,  -1, 0, 0,  -1,-1, 0,   0,-1, 0,&
                    1,-1, 0,  -1, 0, 1,  -1,-1, 1,   0,-1, 1,&
                    1,-1, 1,  -1, 0,-1,  -1,-1,-1,   0,-1,-1,&
                    1,-1,-1,   0, 0,-1],[3,27])

! Group split in cells
type, extends(igroup), public :: cgroup

  ! TODO: Might be the head can be a linked list.

  ! Linked list for atom in each cell (using internal ngroup index)
  integer,allocatable           :: head(:,:,:)
  integer,allocatable           :: next(:)

  ! Cells
  real(dp)                      :: cmax(3),cmin(3)
  integer                       :: ncell=1
  integer                       :: cellaux(3)=[1,1,1]
  real(dp)                      :: cell(3)

  ! Setting this to one force the first cell update
  integer                       :: ncells(3)=[1,1,1]

  ! Cut radius
  real(dp)                      :: rcut=1.e10_dp,rcut2=1.e10_dp

  ! Tessellation
  ! TODO: tessellation update if box change a lot
  ! TODO: tessellation without box using min and max positions.
  logical  :: tessellated=.false.
  procedure(cgroup_tessellate),pointer :: tessellate => cgroup_tessellate

  ! Turn false after first tessellation
  logical                         :: b_out = .true.

  contains

  ! Parent procedures that set a polymorphic pointer to an extension
  ! (see group type construct).
  procedure :: cgroup_attach_atom
  procedure :: cgroup_detach_atom

  procedure :: init => cgroup_construct
  procedure :: dest => cgroup_destroy
  procedure :: attach_atom => cgroup_attach_atom
  procedure :: detach_atom => cgroup_detach_atom

  procedure :: sort => cgroup_sort, cgroup_sort_atom
  procedure :: unsort => cgroup_unsort_atom
             
end type


contains

! cgroup events
! =============

subroutine cgroup_construct(g)
class(cgroup),target         :: g
call g%igroup_construct()
allocate(g%next(g%pad))
end subroutine cgroup_construct

subroutine cgroup_destroy(g)
class(cgroup),target         :: g
call g%igroup%dest()
if(allocated(g%head)) deallocate(g%head)
if(allocated(g%next)) deallocate(g%next)
end subroutine cgroup_destroy

subroutine cgroup_attach_atom(g,a,l_)
use gems_errors, only: silent,errf
class(cgroup),target  :: g
class(atom),target    :: a
integer               :: i
integer,allocatable   :: t_next(:)
integer               :: aux1(dm)
integer,intent(out),optional :: l_
integer                      :: l

! Attempt to attach
call g%igroup_attach_atom(a,l)
if(present(l_)) l_=l
if(l==0) return
 
! Continue reallocations
! TODO: Implement a "preserve" boolean?
! TODO: Skip this if not allocated head?
if(size(g%next)<size(g%a)) then
  allocate(t_next(size(g%a)))
  t_next(:size(g%next)) = g%next(:)
  call move_alloc(to=g%next,from=t_next)
endif
              
! Add atom to cells
if(g%tessellated) then
  i=a%gid(g)

  ! Let me handle the errors
  silent=.true.

  ! Attempt to sort atom
  call cgroup_sort_atom(g,i)

  ! Request full update if fail
  if(errf) g%tessellated=.false.
  silent=.false.

endif
              
end subroutine cgroup_attach_atom

subroutine cgroup_detach_atom(g,a)
! Detach atom `a` from cgroup `g`
class(cgroup)         :: g
class(atom),target    :: a
class(atom),pointer   :: aj
integer               :: i
integer               :: rc(3)

! Attempt to remove atom from cells
if(g%tessellated) then
  i=a%gid(g)
  if(i>0) call cgroup_unsort_atom(g,i)
endif

! Detach atom
call g%igroup_detach_atom(a)

! Clean null items
if(g%update) then
  deallocate(g%next)
  allocate(g%next(size(g%a)))
  g%tessellated=.false.
endif
 
end subroutine cgroup_detach_atom
                 
! Set properties
! ==============

subroutine cgroup_tessellate(g)
class(cgroup)     :: g
integer,parameter :: mincells=60 ! (4x4x4)

! If the box has not changed or if only a small adjustment of the cell size
! is sufficient.
if(g%tessellated) then

  ! Make sure that the box has not grown large enough to include new cells.
  if(all(g%ncells(:)+1>=box(:)/g%rcut)) then

    ! Make sure that the box has not dropped below cut-off radius
    if(all(g%rcut<box(:)/g%ncells(:))) then

      ! Update cell size (needed for NPT)
      call wlog('NHB','Update cell size.',g%b_out)
      g%cell(:)=box(:)/g%ncells(:)
      return

    endif

  endif
endif

! If not a box defined, do not use linked cells
if(.not.boxed) then
  call wlog('NHB','Not using linked cells since box is not defined.',g%b_out)
  g%b_out=.false.
  return
endif
!TODO: else...
!   call get_boundingbox(g%alist,box_max,box_min,g%n(2))
!   call get_boundingbox(g%blist,amax,amin,g%n(4)-g%n(2))
!   box_max(:)=max(box_max(:),amax(:)+pad)
!   box_min(:)=min(box_min(:),amin(:)-pad)
!   set box
! endif

!TODO: if boxed but not PBC? ignore box?
 
! If rcut not defined, do not use linked cells
if(g%rcut==1.e10_dp) then
  call wlog('NHB','Not using linked cells since interaction does not have rcut.',g%b_out)
  g%b_out=.false.
  return
endif

! Redimension of ncells
g%ncells(:)=int(box(:)/g%rcut) ! Number of cells in each direction

! Less than 4 have no sense to use cells
if(all(g%ncells(:)<4)) then
  if(g%b_out) then
    call wlog('NHB','Not using linked cells since it would involve less than 4 cells.')
    call wlog('NHB'); write(logunit,fmt="(a,f10.5)") ' -total cut radious: ', g%rcut
  endif
  g%b_out=.false.
  return
endif

! Cell size
g%ncell = g%ncells(1)*g%ncells(2)*g%ncells(3)
g%cell(:) = box(:)/g%ncells(:)

! Reallocate head array
if(allocated(g%head)) deallocate(g%head)

! Using extra cells around the box size for ghost atmos
allocate(g%head(0:g%ncells(1)+1,0:g%ncells(2)+1,0:g%ncells(3)+1))

! Used to compute cell index
g%cellaux(1) = 1
g%cellaux(2) = g%ncells(1)
g%cellaux(3) = g%ncells(1)*g%ncells(2)

if(g%b_out) then
  call wlog('NHB','Using Linked Cells.')
  call wlog('NHB'); write(logunit,fmt="(a,f10.5)") ' cut radious: ', g%rcut
  call wlog('NHB'); write(logunit,fmt="(a,"//cdm//"(f10.5,2x))") ' cell size: ', g%cell(1:dm)
  call wlog('NHB'); write(logunit,fmt="(a,"//cdm//"(i3,2x))") ' cell numbers: ', g%ncells(1:dm)
  g%b_out=.false.
end if

g%tessellated = .true.

end subroutine cgroup_tessellate

subroutine cgroup_sort(g)
! Sort g atoms into the g cells.
! TODO: Flag when sorting is needed, if not skip.
class(cgroup),intent(inout)  :: g
integer                      :: i,aux1(dm)
type(atom),pointer           :: a

g%head(:,:,:) = 0
do i = 1,g%amax
  if(.not.associated(g%a(i)%o)) cycle
  call cgroup_sort_atom(g,i)
enddo

end subroutine cgroup_sort

subroutine cgroup_sort_atom(g,i)
! Sort atom ith of g into the g cells.
class(cgroup),intent(inout)  :: g
integer                      :: i,ci(dm)
type(atom),pointer           :: a
        
a => g%a(i)%o

ci(:)=int(a%pos(:)/g%cell(:))+1
! j = 1 + dot_product(ci,cellaux)

! If atom pos is outside the cells, mark to update
call werr('Particle out of tessellation',any(ci(:)<0))
call werr('Particle out of tessellation',any(ci(:)>g%ncells(:)+1))

! Silent mode return.
if(errf) return 

g%next(i) = g%head(ci(1),ci(2),ci(3))
g%head(ci(1),ci(2),ci(3)) = i

end subroutine cgroup_sort_atom
                  
subroutine cgroup_unsort_atom(g,i)
! Revert sorting of atom ith and fix cell index.
class(cgroup),intent(inout)  :: g
integer                      :: i,j,k,ci(dm)
type(atom),pointer           :: a
logical                      :: removed

a => g%a(i)%o
  
! Get cell index
ci(:)=int(a%pos(:)/g%cell(:))+1

! TODO: Add to the cgroup a vector similar 
! to hold the cell id for each atom. This will be fast and independent of
! atom position. It may be better to migrate %head(:,:,:) to a %head(:)
! scheme. If atom pos is outside the cells, mark to update
if(any(ci(:)<0)) then
  g%tessellated=.false.
  return
endif  
if(any(ci(:)>g%ncells(:))) then
  g%tessellated=.false.
  return
endif

! Remove i from its cell
removed=.false.
if(g%head(ci(1),ci(2),ci(3))==i) then
  g%head(ci(1),ci(2),ci(3))=g%next(i)
  removed=.true.
else

  j = g%head(ci(1),ci(2),ci(3))
  do while( j>0 )

    k = g%next(j)
    if (k==i) then
      g%next(j)=g%next(i)
      removed=.true.
      exit
    endif
    j = k

  enddo

endif
call werr('Unsort cgroup particle fail',.not.removed)

end subroutine cgroup_unsort_atom
 
#ifdef DIM3

function cell_pbc(g,r,mic)
! Perform PBC on cells. If cell should not be considered return .false.
integer            :: r(3)
class(cgroup),intent(in)  :: g
logical,intent(in) :: mic
logical            :: cell_pbc

cell_pbc=.true.

if (mic) then
  ! May use pbc (it will further depends on atom%pbc)
  r(:)=r(:)-1
  r(:)= mod( r(:) + g%ncells(:), g%ncells(:) )
  r(:)=r(:)+1
else
  ! Can not use PBC
  if(useghost) then
    if(any(r(:)<0)) cell_pbc=.false.
    if(any(r(:)>g%ncells(:)+1)) cell_pbc=.false.
  else
    if(any(r(:)<1)) cell_pbc=.false.
    if(any(r(:)>g%ncells(:))) cell_pbc=.false.
  endif
endif


end function cell_pbc

function vcell2icell(g,i,j,k,mic) result(uid)
!Convert 3-index cell id into 1-index cell id.
!The cell 1,1,1 is map to the cell 1, the 2,1,1 to 2...
class(cgroup),intent(in)  :: g
integer,intent(in) :: i,j,k
integer            :: r(3)
integer            :: uid
logical,intent(in) :: mic

r(1)=i-1
r(2)=j-1
r(3)=k-1

if (mic) then
  ! May use pbc (it will further depends on atom%pbc)
  r(:) =  mod( r(:) + g%ncells(:), g%ncells(:) )
else
  ! Can not use PBC
  if(any(r(:)<1)) then
    uid=0
    return
  endif
  if(any(r(:)>g%ncells(:))) then
    uid=0
    return
  endif
endif

uid = 1 + dot_product(r,g%cellaux(:))

! XXX: Not sure why this (I guess never happend).
if (uid<0) uid=0

end function vcell2icell

function icell2vcell(g,ic) result(r)
!Convert 1-index cell id into 3-index cell id.
!ic shul be 1-based integer
!The cell 1,1,1 is map to the cell 1, the 2,1,1 to 2...
class(cgroup),intent(in)  :: g
integer,intent(in) :: ic
integer            :: r(3),c0

c0=ic-1

r(3)=int(c0/g%cellaux(3))

r(1)=mod(c0,g%cellaux(3))
r(2)=int(r(1)/g%cellaux(2))

r(1)=mod(r(1),g%cellaux(2))

r(:)=r(:)+1

end function icell2vcell

#endif

end module gems_cells
