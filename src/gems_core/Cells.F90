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

! Group split in cells.
! cgroup is a spatial execution view: it reuses group/igroup membership
! machinery, but owns the candidate-space cutoff and cell decomposition
! derived from source atoms.
! Execution-view contract:
! - Source state: inherited igroup index, atom positions, box, and candidate
!   cutoff.
! - Derived state: cell geometry plus head(:,:,:) and next(:) chains.
! - Validity: tessellated means cell geometry is usable; sort must refresh
!   membership chains before linked-cell neighbor search.
! - real_cells_only means the current cell chains intentionally exclude ghosts
!   during the pre-halo stage; ghost attach/detach must not update head(:,:,:)
!   until the following full candidate sort.
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
  logical  :: real_cells_only=.false.
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

  procedure :: sync_real_candidates => cgroup_sync_real_candidates
  procedure :: sync_candidates => cgroup_sync_candidates
  procedure :: sort => cgroup_sort, cgroup_sort_atom
  procedure :: sort_real_atoms => cgroup_sort_real_atoms
  procedure :: unsort => cgroup_unsort_atom
  procedure :: needs_halo => cgroup_needs_halo
  procedure :: halo_cutoff => cgroup_halo_cutoff
  procedure :: is_boundary_cell => cgroup_is_boundary_cell
  procedure :: can_request_boundary_halo => cgroup_can_request_boundary_halo
  procedure :: count_boundary_atoms => cgroup_count_boundary_atoms
  procedure :: boundary_contains_atom => cgroup_boundary_contains_atom
  procedure :: request_halo => cgroup_request_halo
  procedure :: request_boundary_halo => cgroup_request_boundary_halo
  procedure :: request_linear_halo => cgroup_request_linear_halo
             
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

function cgroup_halo_cutoff(g) result(rcut)
! Return the halo thickness required by this spatial candidate view.
! Current ghost storage is global, but the cutoff belongs conceptually to the
! candidate view because different interactions may use different cgroup
! spatial layouts.
class(cgroup),intent(in) :: g
real(dp)                 :: rcut

rcut=g%rcut

end function cgroup_halo_cutoff

logical function cgroup_needs_halo(g)
! Return whether this spatial candidate view requires periodic-image
! synchronization. Today this is simply a positive cutoff, but the method is
! the future hook for view-specific halo policy.
class(cgroup),intent(in) :: g

cgroup_needs_halo = g%halo_cutoff()>0._dp

end function cgroup_needs_halo

logical function cgroup_can_request_boundary_halo(g)
! Return whether this spatial candidate view has the real-atom cell
! decomposition needed to generate halo requests from boundary cells.
!
! A view can require a halo without being able to produce it from cells: small
! boxes or large cutoffs can intentionally disable linked cells. Those cases
! use this view's linear halo producer instead.
class(cgroup),intent(in) :: g

cgroup_can_request_boundary_halo = g%needs_halo() .and. g%tessellated .and. allocated(g%head)

end function cgroup_can_request_boundary_halo

logical function cgroup_is_boundary_cell(g,c)
! Return true when cell c belongs to the periodic halo boundary of this view.
! This is not used to generate ghosts yet. It defines the cell-level predicate
! needed by a future boundary-cell halo path.
class(cgroup),intent(in) :: g
integer,intent(in)       :: c(dm)
integer                  :: n(dm)

cgroup_is_boundary_cell=.false.

if(.not.g%needs_halo()) return
if(.not.g%tessellated) return

! Number of cells needed to cover the halo width in each direction. Current
! tessellation chooses cell size >= rcut, so this is normally one cell; keep the
! formula explicit for future per-view halo widths and nonuniform cells.
n(:)=ceiling(g%halo_cutoff()/g%cell(:))
n(:)=max(n(:),1)
n(:)=min(n(:),g%ncells(:))

cgroup_is_boundary_cell = any(c(:)<=n(:)) .or. any(c(:)>g%ncells(:)-n(:))

end function cgroup_is_boundary_cell

integer function cgroup_count_boundary_atoms(g) result(natoms)
! Count real atoms currently sorted into this view's halo boundary cells.
! Ghost images are derived execution candidates, so this count intentionally
! skips atoms whose prime pointer is associated.
class(cgroup),intent(in) :: g
integer                  :: c(dm),i,j,k,l
type(atom),pointer       :: a

natoms=0

if(.not.g%can_request_boundary_halo()) return

do k=1,g%ncells(3)
  c(3)=k
  do j=1,g%ncells(2)
    c(2)=j
    do i=1,g%ncells(1)
      c(1)=i
      if(.not.g%is_boundary_cell(c)) cycle

      l=g%head(i,j,k)
      do while(l>0)
        a => g%a(l)%o
        if(associated(a)) then
          if(.not.associated(a%prime)) natoms=natoms+1
        endif
        l=g%next(l)
      enddo
    enddo
  enddo
enddo

end function cgroup_count_boundary_atoms

logical function cgroup_boundary_contains_atom(g,a) result(found)
! Return whether real atom a is currently sorted in a boundary cell of this
! candidate view. Used by diagnostics before replacing global halo scans with
! boundary-cell traversal.
class(cgroup),intent(in) :: g
class(atom),target       :: a
integer                  :: c(dm),i,j,k,l

found=.false.

if(.not.g%can_request_boundary_halo()) return
if(associated(a%prime)) return

do k=1,g%ncells(3)
  c(3)=k
  do j=1,g%ncells(2)
    c(2)=j
    do i=1,g%ncells(1)
      c(1)=i
      if(.not.g%is_boundary_cell(c)) cycle

      l=g%head(i,j,k)
      do while(l>0)
        if(associated(g%a(l)%o,target=a)) then
          found=.true.
          return
        endif
        l=g%next(l)
      enddo
    enddo
  enddo
enddo

end function cgroup_boundary_contains_atom

subroutine cgroup_request_boundary_halo(g)
! Request compact periodic images for real atoms sorted into this candidate
! view's boundary cells. The halo layer still owns membership mutation; cgroup
! only owns the cell traversal that identifies real boundary atoms.
use gems_halo, only: clear_atom_global_periodic_halo, sync_atom_global_periodic_halo
class(cgroup),intent(inout) :: g
integer                     :: c(dm),i,j,k,l
type(atom),pointer          :: a
logical                     :: boundary_cell

if(.not.g%can_request_boundary_halo()) return

do k=1,g%ncells(3)
  c(3)=k
  do j=1,g%ncells(2)
    c(2)=j
    do i=1,g%ncells(1)
      c(1)=i
      boundary_cell=g%is_boundary_cell(c)

      l=g%head(i,j,k)
      do while(l>0)
        a => g%a(l)%o
        if(associated(a)) then
          if(.not.associated(a%prime)) then
            if(boundary_cell) then
              call sync_atom_global_periodic_halo(a,g%halo_cutoff())
            else
              call clear_atom_global_periodic_halo(a)
            endif
          endif
        endif
        l=g%next(l)
      enddo
    enddo
  enddo
enddo

end subroutine cgroup_request_boundary_halo

subroutine cgroup_request_halo(g)
! Request compact periodic images for this spatial candidate view. The view
! owns the producer choice: boundary-cell traversal when a tessellation exists,
! otherwise a linear scan over this view's own real candidate membership.
class(cgroup),intent(inout) :: g

if(.not.g%needs_halo()) return

if(g%can_request_boundary_halo()) then
  call g%request_boundary_halo()
else
  call g%request_linear_halo()
endif

end subroutine cgroup_request_halo

subroutine cgroup_request_linear_halo(g)
! Per-view fallback halo producer for candidate views that cannot use linked
! cells. This is intentionally still scoped to the cgroup membership; it is not
! a semantic return to a global sys-wide halo policy.
use gems_halo, only: sync_atom_global_periodic_halo
class(cgroup),intent(inout) :: g
type(atom_dclist),pointer   :: la
type(atom),pointer          :: a
integer                    :: i

if(.not.g%needs_halo()) return

la => g%alist
do i=1,g%nat
  la => la%next
  a => la%o
  if(.not.associated(a)) cycle
  if(associated(a%prime)) cycle
  call sync_atom_global_periodic_halo(a,g%halo_cutoff())
enddo

end subroutine cgroup_request_linear_halo

subroutine cgroup_sync_candidates(g)
! Prepare the spatial candidate view for consumers such as neighbor search.
! cgroup owns the tessellation/sorting order because it owns the cell
! decomposition and its validity flag.
class(cgroup),intent(inout) :: g

call g%tessellate()
if(g%tessellated) call g%sort()

end subroutine cgroup_sync_candidates

subroutine cgroup_sync_real_candidates(g)
! Prepare the spatial layout using only real atoms from this candidate view.
! This stage exists so halo generation can eventually inspect real boundary
! cells before ghost images are activated and the full candidate view is sorted.
class(cgroup),intent(inout) :: g

call g%tessellate()
if(g%tessellated) call g%sort_real_atoms()

end subroutine cgroup_sync_real_candidates

subroutine cgroup_attach_atom(g,a,l_)
! Attach one atom to the indexed spatial view and, when cells are active,
! insert it into the current tessellation without forcing a full rebuild.
use gems_errors, only: silent,errf
class(cgroup),target  :: g
class(atom),target    :: a
integer               :: i
integer,allocatable   :: t_next(:)
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

  if(g%real_cells_only .and. associated(a%prime)) return

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
integer               :: i

! Attempt to remove atom from cells
if(g%tessellated) then
  i=a%gid(g)
  if(i>0) then
    if(.not.(g%real_cells_only .and. associated(a%prime))) call cgroup_unsort_atom(g,i)
  endif
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
! Build or refresh the cell geometry for this spatial view.
! For NPT-style box changes, reuse the existing cell count when it remains
! compatible with the cutoff; otherwise rebuild the tessellation.
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

g%tessellated=.true.

end subroutine cgroup_tessellate

subroutine cgroup_sort(g)
! Sort g atoms into the g cells.
! TODO: Flag when sorting is needed, if not skip.
class(cgroup),intent(inout)  :: g
integer                      :: i

g%head(:,:,:) = 0
g%real_cells_only=.false.
do i = 1,g%amax
  if(.not.associated(g%a(i)%o)) cycle
  call cgroup_sort_atom(g,i)
enddo

end subroutine cgroup_sort

subroutine cgroup_sort_real_atoms(g)
! Sort only real atoms into the g cells. Ghost atoms remain candidate execution
! images and are intentionally excluded from this real-boundary preparation
! stage.
class(cgroup),intent(inout)  :: g
integer                      :: i

g%head(:,:,:) = 0
g%real_cells_only=.true.
do i = 1,g%amax
  if(.not.associated(g%a(i)%o)) cycle
  if(associated(g%a(i)%o%prime)) cycle
  call cgroup_sort_atom(g,i)
enddo

end subroutine cgroup_sort_real_atoms

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
