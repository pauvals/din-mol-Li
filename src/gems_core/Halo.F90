! Copyright (c) 2020  Sergio Alexis Paz
!
!  This file is part of GEMS. GEMS is an Extensible Molecular Simulator.
!
!  GEMS is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
!  GEMS is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with GEMS.  If not, see <https://www.gnu.org/licenses/>.


module gems_halo
use gems_constants, only: dp, dm
implicit none
private

public :: sync_global_periodic_halo, sync_atom_global_periodic_halo
public :: clear_global_periodic_halo, clear_atom_global_periodic_halo
public :: periodic_image_in_global_halo, position_in_global_halo
public :: global_halo_reference_image

contains

subroutine sync_global_periodic_halo(rcut)
! Synchronize the periodic-image execution layer before spatial/neighbor
! views are consumed.
!
! Current implementation still uses preallocated atom ghost images and group
! membership side effects. Architecturally, this routine is the boundary where
! a future cell/boundary based halo implementation should live: real atoms
! remain in sys, ghosts are derived periodic image candidates, and neighbor
! views consume the resulting candidate membership.
use gems_groups, only: useghost, do_pbc, sys
real(dp),intent(in) :: rcut

if(useghost) then
  call sync_global_cutoff_ghosts(rcut)
else
  call do_pbc(sys)
endif

end subroutine sync_global_periodic_halo

subroutine clear_global_periodic_halo()
! Deactivate all compact periodic ghost images while preserving the usual
! real-atom PBC wrapping. This is the reset stage needed by candidate-halo
! producers that rebuild active ghost membership from selected candidate views.
use gems_groups, only: atom, atom_dclist, atom_nghost_images, useghost, do_pbc, sys
type(atom),pointer        :: o
type(atom_dclist),pointer :: la
integer                  :: i,j

call do_pbc(sys)

if(.not.useghost) return

la => sys%alist
do i = 1,sys%nat
  la => la%next
  o => la%o

  if(.not.allocated(o%ghost)) cycle

  do j=1,atom_nghost_images
    call halo_set_ghost_active(o%ghost(j)%o,.false.)
  enddo
enddo

end subroutine clear_global_periodic_halo

subroutine clear_atom_global_periodic_halo(o)
! Deactivate all compact periodic images for one real atom.
use gems_groups, only: atom, atom_nghost_images
type(atom),target,intent(inout) :: o
integer                         :: i

if(.not.allocated(o%ghost)) return

do i=1,atom_nghost_images
  call halo_set_ghost_active(o%ghost(i)%o,.false.)
enddo

end subroutine clear_atom_global_periodic_halo

subroutine sync_global_cutoff_ghosts(rcut)
! Historical compact ghost update moved behind the halo boundary.
! It scans real atoms and activates at most the seven preallocated face/edge/
! corner images used by atom_setpbc.
use gems_groups, only: atom, atom_dclist, sys, do_pbc
type(atom),pointer           :: o
integer                      :: i
type(atom_dclist),pointer    :: la
real(dp),intent(in)          :: rcut

call do_pbc(sys)

la => sys%alist
do i = 1,sys%nat
  la => la%next
  o => la%o

  if(.not.allocated(o%ghost)) cycle

  call sync_atom_global_periodic_halo(o,rcut)

enddo

end subroutine sync_global_cutoff_ghosts

subroutine sync_atom_global_periodic_halo(o,rcut)
! Update the compact periodic images for one real atom using a cutoff width.
! This is used by both the historical global scan and per-view halo producers.
use gems_groups, only: atom, atom_nghost_images
type(atom),target,intent(inout) :: o
real(dp),intent(in)             :: rcut
logical                         :: active(atom_nghost_images)
logical                         :: touched(atom_nghost_images)

call prepare_atom_global_cutoff_requests(o,rcut,active,touched)
call apply_atom_ghost_requests(o,active,touched)

end subroutine sync_atom_global_periodic_halo

subroutine global_halo_reference_image(o,image_id,rcut,active,r)
! Compute the geometry for one compact image exactly as the global cutoff halo
! would, without mutating the preallocated ghost atom. This is used by
! diagnostics that compare candidate-view halo producers against the historical
! global-union reference.
use gems_program_types, only: box
use gems_groups, only: atom, atom_nghost_images
type(atom),intent(in) :: o
integer,intent(in)    :: image_id
real(dp),intent(in)   :: rcut
logical,intent(out)   :: active
real(dp),intent(out)  :: r(dm)
integer               :: reps(dm),boxcr(dm)
integer               :: i,j,k

active=.false.
r(:)=o%pos(:)
if(image_id<1 .or. image_id>atom_nghost_images) return

reps(:)=0
boxcr(:)=0

do k=1,dm
  if(.not.o%pbc(k)) cycle
  if (o%pos(k)>.5_dp*box(k)) then
    reps(k)=-1
  else
    reps(k)=1
  endif
enddo

do k=1,dm
  if(image_id/=k) cycle
  if(.not.o%pbc(k)) return
  boxcr(:)=unit_image_shift(k,reps(k))
  r(:)=o%pos(:)+box(:)*boxcr(:)
  active=halo_contains_image(r,rcut)
  return
enddo

if(count(o%pbc)>1) then
  do j=1,dm-1
    if(.not.o%pbc(j)) cycle
    do k=j+1,dm
      if(.not.o%pbc(k)) cycle
      i=halo_pair_image_index(j,k)
      if(image_id/=i) cycle
      boxcr(:)=pair_image_shift(j,k,reps)
      r(:)=o%pos(:)+box(:)*boxcr(:)
      active=halo_contains_image(r,rcut)
      return
    enddo
  enddo
endif

if(count(o%pbc)>2 .and. image_id==atom_nghost_images) then
  r(:)=o%pos(:)+box(:)*reps(:)
  active=halo_contains_image(r,rcut)
endif

end subroutine global_halo_reference_image

subroutine prepare_atom_global_cutoff_requests(o,rcut,active,touched)
use gems_program_types, only: box
use gems_groups, only: atom, atom_nghost_images
type(atom),target,intent(inout) :: o
real(dp),intent(in)             :: rcut
logical,intent(out)             :: active(atom_nghost_images)
logical,intent(out)             :: touched(atom_nghost_images)
integer                         :: reps(dm),pbc
integer                         :: i,j,k

active(:)=.false.
touched(:)=.false.

pbc=count(o%pbc)
if (pbc==0) return

do k=1,dm
  if(.not.o%pbc(k)) cycle
  if (o%pos(k)>.5_dp*box(k)) then
    reps(k)=-1
  else
    reps(k)=1
  endif
  call prepare_atom_image_request(o,k,unit_image_shift(k,reps(k)),rcut, &
                                  active,touched)
enddo

if (pbc>1) then
  do j=1,dm-1
    if (.not.o%pbc(j)) cycle
    do k=j+1,dm
      if(.not.o%pbc(k)) cycle
      i=halo_pair_image_index(j,k)
      call prepare_atom_image_request(o,i,pair_image_shift(j,k,reps),rcut, &
                                      active,touched)
    enddo
  enddo
endif

if (pbc>2) then
  call prepare_atom_image_request(o,atom_nghost_images,reps,rcut, &
                                  active,touched)
endif

end subroutine prepare_atom_global_cutoff_requests

subroutine prepare_atom_image_request(o,image_id,boxcr,rcut,active,touched)
use gems_program_types, only: box
use gems_groups, only: atom, atom_nghost_images
type(atom),target,intent(inout) :: o
integer,intent(in)              :: image_id
integer,intent(in)              :: boxcr(dm)
real(dp),intent(in)             :: rcut
logical,intent(inout)           :: active(atom_nghost_images)
logical,intent(inout)           :: touched(atom_nghost_images)
type(atom),pointer              :: og

og=>o%ghost(image_id)%o
og%boxcr(:)=boxcr(:)
og%pos(:)=o%pos(:)+box(:)*og%boxcr(:)
active(image_id)=halo_contains_image(og%pos(:),rcut)
touched(image_id)=.true.

end subroutine prepare_atom_image_request

function unit_image_shift(axis,rep) result(boxcr)
integer,intent(in) :: axis,rep
integer            :: boxcr(dm)

boxcr(:)=0
boxcr(axis)=rep

end function unit_image_shift

function pair_image_shift(i,j,reps) result(boxcr)
integer,intent(in) :: i,j
integer,intent(in) :: reps(dm)
integer            :: boxcr(dm)

boxcr(:)=0
boxcr(i)=reps(i)
boxcr(j)=reps(j)

end function pair_image_shift

subroutine apply_atom_ghost_requests(o,active,touched)
use gems_groups, only: atom, atom_nghost_images
type(atom),target,intent(inout) :: o
logical,intent(in)              :: active(atom_nghost_images)
logical,intent(in)              :: touched(atom_nghost_images)
integer                         :: i

do i=1,atom_nghost_images
  if(.not.touched(i)) cycle
  call halo_set_ghost_active(o%ghost(i)%o,active(i))
enddo

end subroutine apply_atom_ghost_requests

subroutine halo_set_ghost_active(og,active)
use gems_groups, only: atom, group, ghost
type(atom),pointer      :: og, o
class(group),pointer    :: g
integer                 :: i
logical,intent(in)      :: active

o=>og%prime

if (active) then

  if(og%ngr>0) return

  do i=1,o%ngr
    g => o%gro(i)
    if(g%ghost) call g%attach(og)
  enddo
  call ghost%attach(og)

else

  if(og%ngr==0) return

  do while (og%ngr/=0)
    g => og%gro(og%ngr)
    call g%detach(og)
  enddo

endif

end subroutine halo_set_ghost_active

logical function periodic_image_in_global_halo(og,rcut) result(active)
use gems_groups, only: atom
type(atom),intent(in) :: og
real(dp),intent(in)   :: rcut

active=halo_contains_image(og%pos(:),rcut)

end function periodic_image_in_global_halo

logical function position_in_global_halo(r,rcut) result(active)
real(dp),intent(in) :: r(dm),rcut

active=halo_contains_image(r,rcut)

end function position_in_global_halo

integer function halo_pair_image_index(i,j) result(idx)
integer,intent(in) :: i,j

idx=i+j+1

end function halo_pair_image_index

function halo_contains_image(r,rcut) result(active)
use gems_program_types, only: box
real(dp),intent(in) :: r(dm),rcut
logical             :: active
integer             :: i

active=.false.
do i=1,dm
  if(r(i)<-rcut) return
  if(r(i)>rcut+box(i)) return
enddo
active=.true.

end function halo_contains_image

end module gems_halo
