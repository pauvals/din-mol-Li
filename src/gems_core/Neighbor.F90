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


module gems_neighbor
use gems_program_types, only: boxed, box, mic
use gems_groups,        only: vdistance, sys, igroup, group, atom, atom_dclist, &
                              atom_nghost_images, useghost, do_pbc
use gems_constants,     only: dp,cdm,dm,ui_ev
use gems_cells,         only: cgroup, map, n1cells, cell_pbc
use gems_errors

implicit none

public :: sync_execution_views
public :: validate_candidate_policies

! Group with neighbors.
! Forces are computed in group `focus` given group `candidate` (without
! newton reaction).
! The interaction is directional at this layer: `focus` is the row/receiver
! side and `candidate` is the searched side. Do not swap them for performance
! unless a concrete kernel later declares a symmetric force-accumulation
! contract that makes the swap legal.
! ngroup gives a unique index to each atom in `focus` or `candidate`.
! ngroup is an execution view. It is class(group)-compatible through the
! indexed group layer, but it owns neighbor execution storage rather than the
! canonical semantic state of the atoms it references.
! Execution-view contract:
! - Source state: focus/candidate memberships, atom positions, box/PBC state,
!   and cutoff.
! - Derived state: inherited local ids plus nn(:) and list(:,:) neighbor ids.
! - Validity: listed means neighbor arrays match the last rebuild criteria.
! - Synchronization: sync_views rebuilds the local neighbor view.
type, extends(igroup), public :: ngroup

  ! The force-receiving group.
  type(group)   :: focus

  ! Logical source groups used to build focus and candidate memberships.
  ! These pointers do not own atoms or membership; they identify the semantic
  ! groups that future shared execution views can use as registry keys.
  class(group),pointer :: focus_source => null()
  class(group),pointer :: candidate_source => null()

  ! Spatial candidate view for possible neighbors. It starts as a private
  ! construction view and is promoted to a shared registry-owned cgroup when
  ! the source/cutoff policy is known.
  type(cgroup),pointer  :: candidate => null()
  logical               :: owns_candidate = .true.
                       
  ! Explicit cutoff policy. A true value means this interaction participates in
  ! spatial candidate synchronization, halo policy, neighbor-list rebuilds, and
  ! Verlet displacement checks. Field-like and bias-like interactions can still
  ! use the ngroup interaction interface without becoming neighbor views.
  real(dp)               :: rcut=1.e10_dp,rcut2=1.e10_dp
  logical                :: has_rcut=.false.

  ! Maximum number of neighbor per atom
  integer                :: mnb=10000

  ! Neighbor list build in this module and used in any forcefield routine.
  integer,allocatable    :: nn(:)        ! Numero de vecinos al atomo i
  integer,allocatable    :: list(:,:)    ! Vecino m-th del i-th atom of the first iteracting group (local and ghost a)
  integer,allocatable    :: focus_id(:)  ! Local ngroup rows that receive forces.
  integer,allocatable    :: real_id(:)   ! Local real-atom row for each local row or ghost.

  ! Position snapshot owned by this neighbor execution view. It is used only
  ! for the Verlet displacement rebuild criterion and is keyed by local atom id.
  real(dp),allocatable    :: pos_old(:,:)

  ! Automatic switch between search algorithm
  logical :: autoswitch=.true.

  logical :: listed=.false.
  procedure(ngroup_cells),pointer :: lista => ngroup_cells
  procedure(ngroup_cells_atom),pointer :: lista_atom => ngroup_cells_atom

  contains

  procedure :: ngroup_construct
  procedure :: ngroup_attach_atom
  procedure :: ngroup_detach_atom
  procedure :: ngroup_sync_views

                                 
    procedure :: init => ngroup_construct
    procedure :: dest => ngroup_destroy
    procedure :: attach_atom => ngroup_attach_atom
    procedure :: detach_atom => ngroup_detach_atom

    procedure :: bind_groups => ngroup_bind_groups
    procedure :: setrc => ngroup_setrc
    procedure :: sync_views => ngroup_sync_views

endtype

type candidate_record
  class(group),pointer :: source => null()
  type(cgroup),pointer :: view => null()
  real(dp)             :: rcut = 0._dp
  integer              :: refs = 0
end type candidate_record

type candidate_ptr
  type(cgroup),pointer :: view => null()
end type candidate_ptr

type ngroup_ptr
  class(ngroup),pointer :: view => null()
end type ngroup_ptr

! Global variables
real(dp),target     :: nb_dcut=1._dp        ! The shell length for verlet update criteria
integer             :: nupd_vlist = 0       ! Number of updates to print in the log
integer             :: nn_vlist = 0         ! Maximum neighbor count seen by dana
real(dp)            :: maxrcut=0.0_dp       ! Maximum cut ratio
logical             :: validate_candidate_policies=.false.
logical             :: b_candidate_policy_env=.true.
logical             :: validate_candidate_halo_equivalence=.false.
logical             :: b_candidate_halo_equivalence_env=.true.
logical             :: use_candidate_halo=.true.
logical             :: b_candidate_halo_env=.true.
type(candidate_record),allocatable :: candidates(:)
integer             :: n_candidates = 0
type(ngroup_ptr),allocatable :: neighbor_views(:)
integer             :: n_neighbor_views = 0

contains

! ngroup events
! =============

subroutine ngroup_construct(g)
class (ngroup),target             :: g

! Init internal groups
call g%focus%init()
call ngroup_allocate_candidate_view(g)

! Initialize the indexed execution view after its role groups. Atom group ids
! are globally sorted, so dynamic membership copying depends on this order.
call g%igroup_construct()
call ngroup_registry_register(g)

! Ghost atoms must have its own index in order to comunicate back
! to the real image the computed properties
g%ghost=.true.

end subroutine ngroup_construct

subroutine ngroup_destroy(g)
class (ngroup),target  :: g

! Destroy internal groups
call g%focus%dest()
g%focus_source=>null()
g%candidate_source=>null()
call ngroup_release_candidate_view(g)
if(allocated(g%nn)) deallocate(g%nn)
if(allocated(g%list)) deallocate(g%list)
if(allocated(g%focus_id)) deallocate(g%focus_id)
if(allocated(g%real_id)) deallocate(g%real_id)
if(allocated(g%pos_old)) deallocate(g%pos_old)
call ngroup_registry_release(g)
call g%igroup%dest()

end subroutine ngroup_destroy

subroutine ngroup_registry_register(g)
! Track live neighbor execution views without depending on concrete kernels.
class(ngroup),target :: g
type(ngroup_ptr),allocatable :: tmp(:)

if(.not.allocated(neighbor_views)) then
  allocate(neighbor_views(8))
elseif(n_neighbor_views==size(neighbor_views)) then
  allocate(tmp(n_neighbor_views+8))
  tmp(:n_neighbor_views)=neighbor_views(:n_neighbor_views)
  call move_alloc(from=tmp,to=neighbor_views)
endif

n_neighbor_views=n_neighbor_views+1
neighbor_views(n_neighbor_views)%view=>g

end subroutine ngroup_registry_register

subroutine ngroup_registry_release(g)
class(ngroup),target :: g
integer              :: i

do i=1,n_neighbor_views
  if(.not.associated(neighbor_views(i)%view,g)) cycle
  if(i<n_neighbor_views) neighbor_views(i:n_neighbor_views-1)=neighbor_views(i+1:n_neighbor_views)
  n_neighbor_views=n_neighbor_views-1
  return
enddo

end subroutine ngroup_registry_release

subroutine ngroup_sync_views(g)
! Rebuild a single neighbor execution view when it becomes stale.
class(ngroup),intent(inout) :: g
real(dp)                    :: dispmax1,dispmax2

if(.not.g%has_rcut) return

if (.not.g%listed) then
  call ngroup_rebuild_view(g)
  return
endif

call inq_dispmax(g,dispmax1,dispmax2)
if (sqrt(dispmax1)+sqrt(dispmax2)>nb_dcut) call ngroup_rebuild_view(g)

end subroutine ngroup_sync_views

subroutine ngroup_rebuild_view(g)
! Recompute the neighbor list for one execution view.
class(ngroup),intent(inout) :: g

nupd_vlist = nupd_vlist + 1
call ngroup_setlista(g)
if(associated(g%lista)) call g%lista()
call ngroup_record_positions(g)

end subroutine ngroup_rebuild_view

subroutine ngroup_allocate_candidate_view(g)
! Create a private candidate cgroup for this neighbor view.
! setrc may later replace it with a shared registry-owned candidate view.
class(ngroup),target :: g

call ngroup_release_candidate_view(g)
allocate(g%candidate)
call g%candidate%init()
g%candidate%ghost=.true.
g%owns_candidate=.true.

end subroutine ngroup_allocate_candidate_view

subroutine ngroup_release_candidate_view(g)
! Release the current candidate view according to its ownership mode.
! Private views are destroyed here; shared registry views decrement refs.
class(ngroup),target :: g

if(associated(g%candidate)) then
  if(g%owns_candidate) then
    call g%candidate%dest()
    deallocate(g%candidate)
  else
    call candidate_registry_release(g%candidate)
  endif
  g%candidate=>null()
endif
g%owns_candidate=.false.

end subroutine ngroup_release_candidate_view

subroutine ngroup_acquire_candidate_view(g,candidate_rcut)
! Select the cgroup that will produce spatial candidates for this neighbor view.
! Compatible source/cutoff policies share one registered candidate view.
class(ngroup),target :: g
real(dp),intent(in)  :: candidate_rcut
integer              :: i

if(.not.associated(g%candidate_source)) then
  g%candidate%rcut=candidate_rcut
  return
endif

i=candidate_registry_find(g%candidate_source,candidate_rcut)
if(i>0) then

  if(associated(g%candidate,candidates(i)%view)) then
    g%candidate%rcut=candidate_rcut
    return
  endif

  call ngroup_release_candidate_view(g)
  g%candidate=>candidates(i)%view
  g%owns_candidate=.false.
  candidates(i)%refs=candidates(i)%refs+1
  return
endif

if(associated(g%candidate)) then
  if(.not.g%owns_candidate) call ngroup_release_candidate_view(g)
endif

if(.not.associated(g%candidate)) then
  allocate(g%candidate)
  call g%candidate%init()
  g%candidate%ghost=.true.
endif

call g%candidate%attach(g%candidate_source)
g%candidate%rcut=candidate_rcut
call candidate_registry_register(g%candidate_source,g%candidate,candidate_rcut)
g%owns_candidate=.false.

end subroutine ngroup_acquire_candidate_view

subroutine ngroup_allocate_neighbor_storage(g,n)
! Allocate the neighbor-list arrays for a cutoff-based ngroup.
! Non-cutoff kernel objects never call this path.
class(ngroup),intent(inout) :: g
integer,intent(in)          :: n

if(allocated(g%nn)) deallocate(g%nn)
if(allocated(g%list)) deallocate(g%list)
if(allocated(g%focus_id)) deallocate(g%focus_id)
if(allocated(g%real_id)) deallocate(g%real_id)
if(allocated(g%pos_old)) deallocate(g%pos_old)

allocate(g%nn(n))
allocate(g%list(n,g%mnb))
allocate(g%focus_id(0))
allocate(g%real_id(n))
g%real_id(:)=0
allocate(g%pos_old(dm,n))
g%pos_old = 1.e8_dp

end subroutine ngroup_allocate_neighbor_storage

subroutine ngroup_compile_focus_ids(g)
! Compile focus membership to local ngroup rows. Kernels consume this dense
! list instead of resolving atom%gid(g) inside force/analysis loops.
class(ngroup),intent(inout) :: g
type(atom),pointer          :: a
type(atom_dclist),pointer   :: la
integer                     :: i, local_id

if(allocated(g%focus_id)) deallocate(g%focus_id)
allocate(g%focus_id(g%focus%nat))

la => g%focus%alist
do i = 1,g%focus%nat
  la => la%next
  a => la%o
  local_id = a%gid(g)
  call werr('Focus atom is missing from compiled neighbor view.',local_id<1)
  g%focus_id(i)=local_id
enddo

end subroutine ngroup_compile_focus_ids

subroutine ngroup_compile_real_ids(g)
! Compile the local real-atom row for every local row. Ghost rows point back to
! their prime atom row, allowing kernels to avoid atom%gid(g) in pair filters.
class(ngroup),intent(inout) :: g
type(atom),pointer          :: a
integer                     :: i, local_id

if(.not.allocated(g%real_id)) allocate(g%real_id(size(g%a)))
if(size(g%real_id)<size(g%a)) then
  deallocate(g%real_id)
  allocate(g%real_id(size(g%a)))
endif

g%real_id(:)=0
do i = 1,g%amax
  if(.not.associated(g%a(i)%o)) cycle
  a => g%a(i)%o
  if(associated(a%prime)) then
    local_id = a%prime%gid(g)
    call werr('Ghost atom prime is missing from compiled neighbor view.',local_id<1)
    g%real_id(i)=local_id
  else
    g%real_id(i)=i
  endif
enddo

end subroutine ngroup_compile_real_ids

integer function candidate_registry_find(source,rcut) result(idx)
! Locate a shared candidate view with the same logical source and halo cutoff.
class(group),intent(in) :: source
real(dp),intent(in)     :: rcut
integer                 :: i

idx=0
do i=1,n_candidates
  if(candidates(i)%source%id/=source%id) cycle
  if(abs(candidates(i)%rcut-rcut)>epsilon(max(abs(rcut),1._dp))) cycle
  idx=i
  return
enddo

end function candidate_registry_find

subroutine candidate_registry_register(source,view,rcut)
! Register a candidate cgroup as shared storage for a source/cutoff policy.
! The registry becomes the owner until the last ngroup releases the view.
class(group),target        :: source
type(cgroup),target        :: view
real(dp),intent(in)        :: rcut
type(candidate_record),allocatable :: tmp(:)

call werr('Candidate registry already has this source/cutoff policy.', &
          candidate_registry_find(source,rcut)>0)

if(.not.allocated(candidates)) then
  allocate(candidates(4))
elseif(n_candidates==size(candidates)) then
  allocate(tmp(n_candidates+4))
  tmp(:n_candidates)=candidates(:)
  call move_alloc(from=tmp,to=candidates)
endif

n_candidates=n_candidates+1
candidates(n_candidates)%source=>source
candidates(n_candidates)%view=>view
candidates(n_candidates)%rcut=rcut
candidates(n_candidates)%refs=1

end subroutine candidate_registry_register

subroutine candidate_registry_release(view)
! Drop one reference to a shared candidate view and destroy it at zero users.
type(cgroup),pointer :: view
integer             :: i

do i=1,n_candidates
  if(.not.associated(candidates(i)%view,view)) cycle

  candidates(i)%refs=candidates(i)%refs-1
  if(candidates(i)%refs<=0) then
    call candidates(i)%view%dest()
    deallocate(candidates(i)%view)
    if(i<n_candidates) candidates(i:n_candidates-1)=candidates(i+1:n_candidates)
    n_candidates=n_candidates-1
  endif
  return
enddo

call werr('Borrowed candidate view is missing from the registry during release.',.true.)

end subroutine candidate_registry_release

logical function candidate_registry_contains(view)
! Return whether a candidate cgroup pointer is owned by the shared registry.
type(cgroup),pointer :: view
integer             :: i

candidate_registry_contains=.false.
do i=1,n_candidates
  if(associated(candidates(i)%view,view)) then
    candidate_registry_contains=.true.
    return
  endif
enddo

end function candidate_registry_contains

integer function active_candidate_view_count() result(n)
! Count spatial candidate views that must be synchronized this step.
! Shared registry views are counted once; private cutoff ngroup views are added.
integer               :: i

n=n_candidates

do i=1,n_neighbor_views
  if(.not.associated(neighbor_views(i)%view)) cycle
  if(.not.neighbor_views(i)%view%has_rcut) cycle
  if(candidate_registry_contains(neighbor_views(i)%view%candidate)) cycle
  n=n+1
enddo

end function active_candidate_view_count

subroutine collect_active_candidates(active,n)
! Materialize the active candidate-view list used by halo and cell sync stages.
! This keeps shared views from being synchronized more than once.
type(candidate_ptr),allocatable,intent(out) :: active(:)
integer,intent(out)                              :: n
integer                                          :: i
integer                                          :: j

n=active_candidate_view_count()
allocate(active(n))

j=0

do i=1,n_candidates
  j=j+1
  active(j)%view => candidates(i)%view
enddo

do i=1,n_neighbor_views
  if(.not.associated(neighbor_views(i)%view)) cycle
  if(.not.neighbor_views(i)%view%has_rcut) cycle
  if(candidate_registry_contains(neighbor_views(i)%view%candidate)) cycle
  j=j+1
  active(j)%view => neighbor_views(i)%view%candidate
enddo

end subroutine collect_active_candidates

subroutine ngroup_bind_groups(g,focus,candidate)
! Bind the logical source groups used by this neighbor execution view.
! `focus` is the force-receiving side; `candidate` is the logical source of
! the spatial candidate view. Keeping this setup in ngroup makes the future
! shared-candidate registry a local change instead of a CLI-side concern.
! This captures the current membership into execution views. Later mutations
! of the logical source groups are not tracked automatically by this binding.
class(ngroup),target :: g
class(group),target  :: focus,candidate

g%focus_source=>focus
g%candidate_source=>candidate
call g%focus%attach(focus)
call g%candidate%attach(candidate)
call g%attach(focus)
if(focus%id/=candidate%id) call g%attach(candidate)

end subroutine ngroup_bind_groups

subroutine ngroup_attach_atom(g,a,l_)
class(ngroup),target   :: g
class(atom),target     :: a
integer, allocatable   :: t_nn(:), t_list(:,:)
integer, allocatable   :: t_real_id(:)
real(dp), allocatable  :: t_pos_old(:,:)
integer                :: n,i
integer,intent(out),optional :: l_
integer                      :: l
logical                      :: focus_atom

call ngroup_attach_atom_roles(g,a)

! Attempt to attach
call g%igroup_attach_atom(a,l)
if(present(l_)) l_=l
if(l==0) return
if(.not.g%has_rcut) return
 
n=size(g%a)
focus_atom = a%gid(g%focus)>0

if(g%listed) then
  if(size(g%nn)<n) then
    allocate(t_nn(n))
    t_nn(:)=0
    t_nn(:size(g%nn)) = g%nn(:)
    call move_alloc(to=g%nn,from=t_nn)

    allocate(t_list(n,g%mnb))
    t_list(:,:)=0
    t_list(:size(g%list,1),:) = g%list(:,:)
    call move_alloc(to=g%list,from=t_list)

    allocate(t_real_id(n))
    t_real_id(:)=0
    t_real_id(:size(g%real_id)) = g%real_id(:)
    call move_alloc(to=g%real_id,from=t_real_id)

    allocate(t_pos_old(dm,n))
    t_pos_old = 1.e8_dp
    t_pos_old(:,:size(g%pos_old,2)) = g%pos_old(:,:)
    call move_alloc(to=g%pos_old,from=t_pos_old)
  endif

  i=a%id(l)
  if(associated(a%prime)) then
    g%real_id(i)=a%prime%gid(g)
  else
    g%real_id(i)=i
  endif
  call ngroup_sort_atom(g,i)
  if(focus_atom .and. .not.associated(a%prime)) call ngroup_compile_focus_ids(g)
  g%pos_old(:,i)=a%pos(:)
else
  if(size(g%nn)<n) then
    call ngroup_allocate_neighbor_storage(g,n)
  endif
endif
        
end subroutine ngroup_attach_atom

subroutine ngroup_attach_atom_roles(g,a)
! Keep focus/candidate role membership ahead of local neighbor indexing.
! Dynamic paths such as GCMC and halo activation may encounter the ngroup view
! before its internal role groups because group ids are globally sorted.
class(ngroup),target :: g
class(atom),target   :: a
class(atom),pointer  :: source_atom

source_atom => a
if(associated(a%prime)) source_atom => a%prime

if(.not.associated(a%prime)) then
  if(ngroup_atom_in_source(g%focus_source,source_atom)) call g%focus%attach(a)
endif

if(ngroup_atom_in_source(g%candidate_source,source_atom)) call g%candidate%attach(a)

end subroutine ngroup_attach_atom_roles

logical function ngroup_atom_in_source(source,a) result(found)
class(group),pointer :: source
class(atom),target   :: a

found=.true.
if(associated(source)) found=a%gri(source)>0

end function ngroup_atom_in_source

subroutine ngroup_detach_atom(g,a)
class(ngroup)          :: g
class(atom),target     :: a
integer                :: i
logical                :: real_atom
logical                :: focus_atom

if(.not.g%has_rcut) then
  if(a%gid(g)<1) return
  call g%focus%detach(a)
  call g%candidate%detach(a)
  call g%igroup_detach_atom(a)
  return
endif

real_atom=.not.associated(a%prime)
focus_atom=a%gid(g%focus)>0
if(g%listed) then
  i=a%gid(g)
  if(i<1) return
  g%nn(i)=0
else
  if(a%gid(g)<1) return
endif

! Detach atom
call g%focus%detach(a)
call g%candidate%detach(a)
call g%igroup_detach_atom(a)
 
! Clean null items
if(g%update) then
  call ngroup_allocate_neighbor_storage(g,size(g%a))
  g%listed=.false.
endif

if(g%listed) then
  g%a(i)%o=>g%limbo
  if(allocated(g%real_id)) g%real_id(i)=0
  g%nlimbo=g%nlimbo+1
  g%b_limbo=.true.
  if(real_atom .and. focus_atom) call ngroup_compile_focus_ids(g)
endif

end subroutine ngroup_detach_atom

subroutine ngroup_sort_atom(g,i)
! Incrementally update neighbor execution arrays for one newly attached atom.
! Used when the full neighbor view is already listed, especially for ghost
! switches inside sync_execution_views.
class(ngroup),intent(inout)  :: g
type(atom),pointer           :: a,aj
integer                      :: i,j,jj,m
real(dp)                     :: rd,vd(dm)
real(dp)                     :: rcut
type(atom_dclist),pointer    :: la

g%nn(i)=0
a => g%a(i)%o

if(a%gid(g%focus)==0) call g%lista_atom(i)

if(a%gid(g%candidate)>0) then

  rcut=(g%rcut+nb_dcut)
  rcut=rcut*rcut

  la => g%focus%alist
  do jj = 1,g%focus%nat
    la => la%next
    aj => la%o
    j = aj%gid(g)
    if(j<1 .or. j>size(g%nn)) cycle

    if(associated(aj,target=a)) cycle

    call vdistance(vd,a,aj,mic)
    rd = dot_product(vd,vd)

    if (rd>rcut) cycle

    m = g%nn(j)+1
    g%list(j,m)=i
    g%nn(j)=m

  enddo

endif

end subroutine ngroup_sort_atom
     
! Set properties
! ==============

subroutine ngroup_setrc(g,rc)
! I could use the g%candidate%rcut directly.. but may be confuse...
class(ngroup)       :: g
real(dp),intent(in) :: rc

! New cut radious
g%rcut = rc
g%rcut2 = rc*rc
g%has_rcut=.true.
if(.not.allocated(g%nn)) call ngroup_allocate_neighbor_storage(g,size(g%a))
maxrcut=max(maxrcut,rc+nb_dcut)
call ngroup_acquire_candidate_view(g,rc+nb_dcut)

end subroutine ngroup_setrc

subroutine ngroup_setlista(g)
class(ngroup) :: g

! Update cells
! FIXME:DSAFR: tessellate must be here
! call g%candidate%tessellate()

if(.not.g%autoswitch) return

if(g%candidate%tessellated) then
  g%lista => ngroup_cells
  g%lista_atom => ngroup_cells_atom
else
  g%lista => ngroup_verlet
  g%lista_atom => ngroup_verlet_atom
endif

end subroutine ngroup_setlista

! Search algorithms
! =================

subroutine ngroup_verlet(g)
! Build neighbors verlet list.
class(ngroup),intent(inout)  :: g
type(atom),pointer           :: ai,aj
integer                      :: i,ii,j,m
real(dp)                     :: rd,vd(dm)
real(dp)                     :: rcut
type(atom_dclist),pointer    :: la

! Set ceros (TODO: I think this is not needed)
g%nn(:)=0

! Cut radious
rcut=(g%rcut+nb_dcut)
rcut=rcut*rcut

! Clean atoms in limbo
call g%clean()
call ngroup_compile_real_ids(g)
call ngroup_compile_focus_ids(g)

!OMP: Creo varios threads
!$OMP PARALLEL

!OMP: Un thread ejecuta, los demas quedan idle.
!$OMP SINGLE
la => g%focus%alist
do ii = 1,g%focus%nat
  la => la%next
  ai => la%o
  i = g%focus_id(ii)

  !OMP: Se reparte la tarea entre los idle threads
  !$OMP TASK DEFAULT(NONE)     &
  !$OMP& FIRSTPRIVATE(ai,i)    &
  !$OMP& SHARED(g,rcut,mic)    &
  !$OMP& PRIVATE(aj,j,vd,rd,m)

  ! Reset number of neighbors
  m=0

  do j = 1, g%candidate%amax
    aj => g%candidate%a(j)%o
    if(.not.associated(aj)) cycle

    ! Skip autointeraction
    if(associated(aj,target=ai)) cycle

    call vdistance(vd,ai,aj,mic)
    rd = dot_product(vd,vd)

    if (rd>rcut) cycle

    ! Add j as neighbor of i.
    m=m+1
    g%list(i,m)=aj%gid(g)

  enddo
  g%nn(i)=m

  !$OMP END TASK
enddo

!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL

g%listed=.true.
if(g%focus%nat>0) nn_vlist=max(nn_vlist,maxval(g%nn(g%focus_id(:))))

end subroutine ngroup_verlet

subroutine ngroup_verlet_atom(g,i)
! Search neighbors for atom i.
class(ngroup),intent(inout)  :: g
type(atom),pointer           :: ai,aj
integer                      :: i,j,m
real(dp)                     :: rd,vd(dm)
real(dp)                     :: rcut

g%nn(i)=0
ai => g%a(i)%o

rcut=(g%rcut+nb_dcut)
rcut=rcut*rcut

m=0
do j = 1, g%candidate%amax
  aj=>g%candidate%a(j)%o
  if(.not.associated(aj)) cycle

  if(associated(aj,target=ai)) cycle

  call vdistance(vd,ai,aj,mic)
  rd = dot_product(vd,vd)

  if (rd>rcut) cycle

  m=m+1
  g%list(i,m)=aj%gid(g)

enddo
g%nn(i)=m

end subroutine ngroup_verlet_atom

subroutine ngroup_cells(g)
! Build neighbors verlet list over linked cells.
class(ngroup),intent(inout)  :: g
type(atom), pointer          :: ai, aj
integer                      :: i,ii,j,k
integer                      :: nabor
real(dp)                     :: rd,vd(dm)
real(dp)                     :: rcut
integer                      :: rc(dm),nc(dm)
type(atom_dclist),pointer    :: la

! Set ceros
g%nn(:)=0
 
! Cut radious
rcut=(g%rcut+nb_dcut)
rcut=rcut*rcut

! Clean atoms in limbo
call g%clean
call ngroup_compile_real_ids(g)
call ngroup_compile_focus_ids(g)
   
! Sort in cells
! FIXME:DSAFR: sort must be here
! call g%candidate%sort()

!Bucle sobre los atomos de g
la => g%focus%alist
!$OMP PARALLEL
!$OMP SINGLE
do ii = 1,g%focus%nat
  la => la%next
  ai => la%o

  !$OMP TASK DEFAULT(NONE)             &
  !$OMP& FIRSTPRIVATE(ai,ii)           &
  !$OMP& SHARED(g,rcut,mic)            &
  !$OMP& PRIVATE(rc,i,nabor,nc,j,aj,k,vd,rd)

  i = g%focus_id(ii)

  ! Get cell index
  rc(:)=int(ai%pos(:)/g%candidate%cell(:))+1

  !Bucle sobre las celdas vecinas a rc(:)
  do nabor=0,26
    nc(:)=map(:,nabor)+rc(:)

    ! Check if cells should be considered and apply PBC
    if(.not.cell_pbc(g%candidate,nc,mic)) cycle

    !Bucle sobre los atomos de la celda vecina
    j = g%candidate%head(nc(1),nc(2),nc(3))
    do while( j>0 )
      aj=>g%candidate%a(j)%o
      k = aj%gid(g)

      ! Next here to allow cycle
      j = g%candidate%next(j)

      ! Skip autointeraction
      if(associated(aj,target=ai)) cycle

      call vdistance(vd,aj,ai,mic) ! respetar el orden
      rd =  dot_product(vd,vd)

      if (rd<rcut) then
        g%nn(i)=g%nn(i)+1
        g%list(i,g%nn(i))=k
      endif

    enddo

  enddo

  !$OMP END TASK

enddo

!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL

g%listed=.true.
if(g%focus%nat>0) nn_vlist=max(nn_vlist,maxval(g%nn(g%focus_id(:))))
 
end subroutine ngroup_cells

subroutine ngroup_cells_atom(g,i)
! Search neighbors for atom i. Assume candidate is already sorted in cells.
class(ngroup),intent(inout)  :: g
type(atom), pointer          :: ai, aj
integer                      :: i,j,k
integer                      :: nabor
real(dp)                     :: rd,vd(dm)
real(dp)                     :: rcut
integer                      :: rc(dm),nc(dm)

g%nn(i)=0
ai => g%a(i)%o

rcut=(g%rcut+nb_dcut)
rcut=rcut*rcut

rc(:)=int(ai%pos(:)/g%candidate%cell(:))+1

do nabor=0,26
  nc(:)=map(:,nabor)+rc(:)

  if(.not.cell_pbc(g%candidate,nc,mic)) cycle

  j = g%candidate%head(nc(1),nc(2),nc(3))
  do while( j>0 )
    aj=>g%candidate%a(j)%o
    if(.not.associated(aj)) print *, j
    k = aj%gid(g)

    j = g%candidate%next(j)

    if(associated(aj,target=ai)) cycle

    call vdistance(vd,aj,ai,mic)
    rd =  dot_product(vd,vd)

    if (rd<rcut) then
      g%nn(i)=g%nn(i)+1
      g%list(i,g%nn(i))=k
    endif

  enddo

enddo

end subroutine ngroup_cells_atom

! Update search
! =============

subroutine rebuild_neighbor_views()
! Rebuild all registered neighbor execution views.
integer                      :: i

! call system_clock(t1)
nupd_vlist = nupd_vlist +1

do i = 1,n_neighbor_views
  if(.not.associated(neighbor_views(i)%view)) cycle
  call neighbor_views(i)%view%sync_views()
enddo

end subroutine rebuild_neighbor_views

subroutine ngroup_record_positions(g)
class(ngroup),intent(inout) :: g
integer                     :: i

do i = 1,g%amax
  if(.not.associated(g%a(i)%o)) cycle
  if(associated(g%a(i)%o,target=g%limbo)) cycle
  g%pos_old(:,i)=g%a(i)%o%pos(:)
enddo

end subroutine ngroup_record_positions

subroutine inq_dispmax(g,dispmax1,dispmax2)
class(ngroup)             :: g
type(atom_dclist),pointer :: la
integer                   :: i,id
real(dp),intent(out)      :: dispmax1, dispmax2
real(dp)                  :: rd,vd(dm)

!vecinos con el propio sistema
dispmax1 = 1.e-16_dp
dispmax2 = 1.e-16_dp

la => g%alist
do i = 1,g%nat
  la => la%next

  if(associated(la%o%prime)) cycle

  id = la%o%gid(g)
  if(id<1 .or. id>size(g%pos_old,2)) then
    rd = huge(rd)
  else
    vd(:) = la%o%pos(:) - g%pos_old(:,id)
    rd = dot_product(vd,vd)
  endif

  if (rd>dispmax1) then
    dispmax2 = dispmax1
    dispmax1 = rd
  else
    if (rd>dispmax2) dispmax2 = rd
  endif

enddo

end subroutine inq_dispmax

subroutine sync_execution_views()
! Synchronize non-bonded execution views before kernels consume them.
integer                    :: i

call sync_spatial_real_candidates()
call sync_spatial_candidate_halos()
if(candidate_policy_diagnostics_active()) call validate_boundary_halo_coverage()
if(candidate_halo_equivalence_diagnostics_active()) call validate_global_ghost_candidate_halo_equivalence()
call sync_spatial_candidates()

call sync_candidate_policy_diagnostics()

do i = 1,n_neighbor_views
  if(.not.associated(neighbor_views(i)%view)) cycle
  call neighbor_views(i)%view%sync_views()
enddo

end subroutine sync_execution_views

subroutine sync_spatial_real_candidates()
! Prepare each spatial candidate view using only real atoms. This does not feed
! neighbor kernels directly; it creates a stable pre-halo spatial stage for
! future boundary-cell halo generation.
integer               :: i,n
type(candidate_ptr),allocatable :: active(:)

! Use the active candidate-view iterator so registered shared views and private
! ngroup views follow the same has-cutoff and duplicate-suppression contract.
call collect_active_candidates(active,n)
do i = 1,n
  if(.not.associated(active(i)%view)) cycle
  call active(i)%view%sync_real_candidates()
enddo

end subroutine sync_spatial_real_candidates

subroutine sync_spatial_candidates()
! Synchronize each spatial candidate view once.
! Several ngroups can share the same cgroup through the candidate registry; cell
! tessellation/sorting belongs to that shared execution view, not to each
! neighbor-list consumer.
integer               :: i,n
type(candidate_ptr),allocatable :: active(:)

call collect_active_candidates(active,n)
do i = 1,n
  if(.not.associated(active(i)%view)) cycle
  call active(i)%view%sync_candidates()
enddo

end subroutine sync_spatial_candidates

subroutine sync_spatial_candidate_halos()
! Synchronize periodic candidate ghosts before full spatial candidate views are
! rebuilt. The normal path asks each active cgroup to produce its own halo; the
! global-union path remains available as a conservative reference/debug mode.
use gems_halo, only: sync_global_periodic_halo
real(dp) :: halo_rcut

if(candidate_halo_active()) then
  call sync_candidate_view_halos()
  return
endif

halo_rcut = candidate_registry_union_halo_cutoff()
if(halo_rcut>0._dp) then
  call sync_global_periodic_halo(halo_rcut)
else
  call sync_global_periodic_halo(maxrcut)
endif

end subroutine sync_spatial_candidate_halos

subroutine sync_candidate_view_halos()
! Candidate-view halo synchronization. Each cgroup owns its producer:
! cell-capable views request ghosts from real boundary cells, and views without
! a usable tessellation scan their own real candidate membership.
use gems_halo, only: sync_global_periodic_halo
real(dp) :: halo_rcut

if(.not.useghost) then
  halo_rcut=candidate_registry_union_halo_cutoff()
  if(halo_rcut>0._dp) then
    call sync_global_periodic_halo(halo_rcut)
  else
    call sync_global_periodic_halo(maxrcut)
  endif
  return
endif

halo_rcut = candidate_registry_union_halo_cutoff()
if(halo_rcut<=0._dp .and. maxrcut>0._dp) then
  call sync_global_periodic_halo(maxrcut)
else
  call do_pbc(sys)
  call request_candidate_halos_by_cutoff()
  call clear_stale_candidate_halos()
endif

end subroutine sync_candidate_view_halos

subroutine clear_stale_candidate_halos()
! Remove active compact ghosts only for real atoms that no longer belong to any
! active cutoff-based candidate view. Atoms still owned by a candidate view were
! just synchronized by that view's halo producer.
use gems_halo, only: clear_atom_global_periodic_halo
type(atom_dclist),pointer :: la
type(atom),pointer        :: a
integer                  :: i
type(candidate_ptr),allocatable :: active(:)
integer                  :: n_active

call collect_active_candidates(active,n_active)

la => sys%alist
do i=1,sys%nat
  la => la%next
  a => la%o
  if(.not.associated(a)) cycle
  if(atom_in_candidate_halo_views(a,active,n_active)) cycle
  call clear_atom_global_periodic_halo(a)
enddo

end subroutine clear_stale_candidate_halos

logical function atom_in_candidate_halo_views(a,active,n) result(found)
type(atom),intent(in)   :: a
type(candidate_ptr),intent(in) :: active(:)
integer,intent(in)      :: n
integer                 :: i

found=.false.

do i = 1,n
  if(.not.associated(active(i)%view)) cycle
  if(.not.active(i)%view%needs_halo()) cycle
  if(a%gid(active(i)%view)>0) then
    found=.true.
    return
  endif
enddo

end function atom_in_candidate_halo_views

subroutine request_candidate_halos_by_cutoff()
! `sync_atom_global_periodic_halo` replaces the active image set for one atom.
! Request per-view halos in increasing cutoff order so larger policies win when
! the same real atom belongs to more than one candidate view.
real(dp) :: previous_rcut,next_rcut
real(dp),parameter :: initial_rcut=-1._dp
type(candidate_ptr),allocatable :: active(:)
integer :: n_active

call collect_active_candidates(active,n_active)

previous_rcut=initial_rcut
do
  next_rcut=next_candidate_halo_cutoff(previous_rcut,active,n_active)
  if(next_rcut<0._dp) exit
  call request_candidate_halos_at_cutoff(next_rcut,active,n_active)
  previous_rcut=next_rcut
enddo

end subroutine request_candidate_halos_by_cutoff

function next_candidate_halo_cutoff(previous_rcut,active,n) result(next_rcut)
real(dp),intent(in)   :: previous_rcut
type(candidate_ptr),intent(in) :: active(:)
integer,intent(in)    :: n
real(dp)              :: next_rcut,rcut
integer               :: i

next_rcut=-1._dp

do i = 1,n
  if(.not.associated(active(i)%view)) cycle
  if(.not.active(i)%view%needs_halo()) cycle
  rcut=active(i)%view%halo_cutoff()
  if(rcut<=previous_rcut) cycle
  if(next_rcut<0._dp .or. rcut<next_rcut) next_rcut=rcut
enddo

end function next_candidate_halo_cutoff

subroutine request_candidate_halos_at_cutoff(target_rcut,active,n)
real(dp),intent(in)   :: target_rcut
type(candidate_ptr),intent(in) :: active(:)
integer,intent(in)    :: n
real(dp)              :: rcut,tol
integer               :: i

tol=epsilon(max(1._dp,abs(target_rcut)))

do i = 1,n
  if(.not.associated(active(i)%view)) cycle
  if(.not.active(i)%view%needs_halo()) cycle
  rcut=active(i)%view%halo_cutoff()
  if(abs(rcut-target_rcut)>tol) cycle
  call active(i)%view%request_halo()
enddo

end subroutine request_candidate_halos_at_cutoff

function candidate_registry_union_halo_cutoff() result(halo_rcut)
type(candidate_ptr),allocatable :: active(:)
real(dp)                :: halo_rcut,rcut
integer                 :: i,n

halo_rcut=0._dp

call collect_active_candidates(active,n)
do i = 1,n
  if(.not.associated(active(i)%view)) cycle
  if(.not.active(i)%view%needs_halo()) cycle
  rcut=active(i)%view%halo_cutoff()
  if(rcut<=halo_rcut) cycle
  halo_rcut=rcut
enddo

end function candidate_registry_union_halo_cutoff

logical function candidate_halo_active() result(active)
character(16) :: value
integer       :: stat

if(b_candidate_halo_env) then
  call get_environment_variable('GEMS_USE_GLOBAL_UNION_HALO',value,status=stat)
  if(stat==0) then
    select case(trim(adjustl(value)))
    case('1','t','T','true','TRUE','.true.','.TRUE.')
      use_candidate_halo=.false.
    end select
  endif
  b_candidate_halo_env=.false.
endif

active=use_candidate_halo

end function candidate_halo_active

subroutine sync_candidate_policy_diagnostics()

if(candidate_policy_diagnostics_active()) then
  call validate_candidate_registry()
  call validate_candidate_policy_duplicates()
endif

end subroutine sync_candidate_policy_diagnostics

logical function candidate_policy_diagnostics_active() result(active)
character(16) :: value
integer       :: stat

if(b_candidate_policy_env) then
  call get_environment_variable('GEMS_VALIDATE_CANDIDATE_POLICIES',value,status=stat)
  if(stat==0) then
    select case(trim(adjustl(value)))
    case('1','t','T','true','TRUE','.true.','.TRUE.')
      validate_candidate_policies=.true.
    end select
  endif
  b_candidate_policy_env=.false.
endif

active=validate_candidate_policies

end function candidate_policy_diagnostics_active

logical function candidate_halo_equivalence_diagnostics_active() result(active)
character(16) :: value
integer       :: stat

if(b_candidate_halo_equivalence_env) then
  call get_environment_variable('GEMS_VALIDATE_CANDIDATE_HALO_EQUIVALENCE',value,status=stat)
  if(stat==0) then
    select case(trim(adjustl(value)))
    case('1','t','T','true','TRUE','.true.','.TRUE.')
      validate_candidate_halo_equivalence=.true.
    end select
  endif
  b_candidate_halo_equivalence_env=.false.
endif

active=validate_candidate_halo_equivalence

end function candidate_halo_equivalence_diagnostics_active

subroutine validate_global_ghost_candidate_halo_equivalence()
use gems_halo, only: global_halo_reference_image
type(atom_dclist),pointer :: la
type(atom),pointer        :: a
real(dp)                  :: halo_rcut
real(dp)                  :: reference_pos(dm)
type(candidate_ptr),allocatable :: active(:)
logical                   :: reference_active
integer                   :: i,j
integer                   :: n_active

halo_rcut = candidate_registry_union_halo_cutoff()
if(halo_rcut<=0._dp) halo_rcut=maxrcut
if(halo_rcut<=0._dp) return

call collect_active_candidates(active,n_active)

la => sys%alist
do i=1,sys%nat
  la => la%next
  a => la%o
  if(.not.associated(a)) cycle
  if(.not.allocated(a%ghost)) cycle
  if(.not.atom_in_candidate_halo_views(a,active,n_active)) cycle

  do j=1,atom_nghost_images
    call global_halo_reference_image(a,j,halo_rcut,reference_active,reference_pos)
    if(.not.reference_active) cycle
    call werr('Candidate-view halo producers would miss a global-union reference ghost.', &
              .not.boundary_views_cover_ghost(a,reference_pos,active,n_active))
  enddo
enddo

end subroutine validate_global_ghost_candidate_halo_equivalence

logical function boundary_views_cover_ghost(prime,ghost_pos,active,n) result(covered)
use gems_halo, only: position_in_global_halo
type(atom),target,intent(in) :: prime
real(dp),intent(in)      :: ghost_pos(dm)
type(candidate_ptr),intent(in) :: active(:)
integer,intent(in)       :: n
integer                  :: i

covered=.false.

do i = 1,n
  if(.not.associated(active(i)%view)) cycle
  if(candidate_view_covers_ghost(active(i)%view,prime,ghost_pos)) then
    covered=.true.
    return
  endif
  if(candidate_view_global_fallback_covers_ghost(active(i)%view,ghost_pos)) then
    covered=.true.
    return
  endif
enddo

contains

logical function candidate_view_covers_ghost(view,real_atom,r) result(covers)
type(cgroup),intent(in)      :: view
type(atom),target,intent(in) :: real_atom
real(dp),intent(in)          :: r(dm)

covers=.false.
if(.not.view%can_request_boundary_halo()) return
if(.not.view%boundary_contains_atom(real_atom)) return
if(.not.position_in_global_halo(r,view%halo_cutoff())) return

covers=.true.

end function candidate_view_covers_ghost

logical function candidate_view_global_fallback_covers_ghost(view,r) result(covers)
type(cgroup),intent(in) :: view
real(dp),intent(in)     :: r(dm)

covers=.false.
if(view%can_request_boundary_halo()) return
if(.not.view%needs_halo()) return
if(.not.position_in_global_halo(r,view%halo_cutoff())) return

covers=.true.

end function candidate_view_global_fallback_covers_ghost

end function boundary_views_cover_ghost

subroutine validate_boundary_halo_coverage()
integer               :: i,n
type(candidate_ptr),allocatable :: active(:)

call collect_active_candidates(active,n)
do i = 1,n
  if(.not.associated(active(i)%view)) cycle
  call validate_candidate_boundary_halo_coverage(active(i)%view)
enddo

end subroutine validate_boundary_halo_coverage

subroutine validate_candidate_boundary_halo_coverage(view)
type(cgroup),intent(in)  :: view
type(atom_dclist),pointer :: la
type(atom),pointer        :: a
integer                   :: i

if(.not.view%can_request_boundary_halo()) return

la => view%alist
do i=1,view%nat
  la => la%next
  a => la%o
  if(.not.associated(a%prime)) cycle
  call werr('Boundary halo traversal would miss an active candidate ghost.', &
            .not.view%boundary_contains_atom(a%prime))
enddo

end subroutine validate_candidate_boundary_halo_coverage

subroutine validate_candidate_registry()
integer                :: i,j,actual_refs

do i=1,n_candidates
  call werr('Candidate registry record has no source group.', &
            .not.associated(candidates(i)%source))
  call werr('Candidate registry record has no candidate view.', &
            .not.associated(candidates(i)%view))
  call werr('Candidate registry record has invalid reference count.', &
            candidates(i)%refs<=0)

  actual_refs=0
  do j=1,n_neighbor_views
    if(.not.associated(neighbor_views(j)%view)) cycle
    if(associated(neighbor_views(j)%view%candidate,candidates(i)%view)) actual_refs=actual_refs+1
  enddo

  call werr('Candidate registry reference count does not match ngroup users.', &
            candidates(i)%refs/=actual_refs)
  call validate_candidate_boundary_cells(candidates(i)%view)
enddo

do i=1,n_neighbor_views
  if(.not.associated(neighbor_views(i)%view)) cycle
  if(.not.neighbor_views(i)%view%has_rcut) cycle
  if(neighbor_views(i)%view%owns_candidate) cycle
  call werr('Borrowed candidate view is missing from the registry.', &
            .not.candidate_registry_contains(neighbor_views(i)%view%candidate))
enddo

end subroutine validate_candidate_registry

subroutine validate_candidate_boundary_cells(view)
type(cgroup),intent(in) :: view
integer                 :: c(dm),n_boundary,n_boundary_atoms
integer                 :: i,j,k
real(dp)                :: tol

if(.not.view%can_request_boundary_halo()) return

tol=100._dp*epsilon(max(1._dp,maxval(abs(box(:)))))

call werr('Candidate tessellation has invalid cell counts.',any(view%ncells(:)<=0))
call werr('Candidate tessellation has inconsistent total cell count.', &
          view%ncell/=product(view%ncells(:)))
call werr('Candidate tessellation cell size is inconsistent with current box.', &
          any(abs(view%cell(:)*view%ncells(:)-box(:))>tol))
call werr('Candidate tessellation cell size fell below cutoff.', &
          any(view%cell(:)<view%rcut-tol))

n_boundary=0
do k=1,view%ncells(3)
  c(3)=k
  do j=1,view%ncells(2)
    c(2)=j
    do i=1,view%ncells(1)
      c(1)=i
      if(view%is_boundary_cell(c)) n_boundary=n_boundary+1
    enddo
  enddo
enddo

call werr('Candidate halo boundary has no cells.',n_boundary<=0)
call werr('Candidate halo boundary exceeds tessellated cell count.',n_boundary>view%ncell)

n_boundary_atoms=view%count_boundary_atoms()
call werr('Candidate halo boundary atom count is invalid.', &
          n_boundary_atoms<0 .or. n_boundary_atoms>view%nat)

end subroutine validate_candidate_boundary_cells

subroutine validate_candidate_policy_duplicates()
integer                :: i,j

do i=1,n_neighbor_views
  if(.not.associated(neighbor_views(i)%view)) cycle
  do j=1,i-1
    if(.not.associated(neighbor_views(j)%view)) cycle
    if(.not.same_candidate_policy(neighbor_views(i)%view,neighbor_views(j)%view)) cycle
    call werr('Candidate policy duplicate is not sharing its candidate view.', &
              .not.associated(neighbor_views(i)%view%candidate,neighbor_views(j)%view%candidate))
    call validate_candidate_view_match(neighbor_views(i)%view,neighbor_views(j)%view)
  enddo
enddo

end subroutine validate_candidate_policy_duplicates

subroutine validate_candidate_view_match(a,candidate)
class(ngroup),intent(in) :: a,candidate
integer                  :: i

call werr('Candidate policy duplicate has different candidate membership count.',a%candidate%nat/=candidate%candidate%nat)
call werr('Candidate policy duplicate has different indexed capacity.',a%candidate%amax/=candidate%candidate%amax)
call werr('Candidate policy duplicate has different tessellation state.', &
          a%candidate%tessellated.neqv.candidate%candidate%tessellated)

do i=1,a%candidate%amax
  call werr('Candidate policy duplicate has different atom association state.', &
            associated(a%candidate%a(i)%o).neqv.associated(candidate%candidate%a(i)%o))
  if(.not.associated(a%candidate%a(i)%o)) cycle
  call werr('Candidate policy duplicate has different atom ordering.', &
            .not.associated(a%candidate%a(i)%o,target=candidate%candidate%a(i)%o))
enddo

if(.not.a%candidate%tessellated) return

call werr('Candidate policy duplicate has different cell counts.',any(a%candidate%ncells(:)/=candidate%candidate%ncells(:)))
call werr('Candidate policy duplicate has different cell sizes.', &
          any(abs(a%candidate%cell(:)-candidate%candidate%cell(:))>epsilon(1._dp)))
call werr('Candidate policy duplicate has different cell head shape.',any(shape(a%candidate%head)/=shape(candidate%candidate%head)))
call werr('Candidate policy duplicate has different cell heads.',any(a%candidate%head/=candidate%candidate%head))
call werr('Candidate policy duplicate has different cell next arrays.', &
          any(a%candidate%next(:a%candidate%amax)/= &
              candidate%candidate%next(:candidate%candidate%amax)))

end subroutine validate_candidate_view_match

logical function same_candidate_policy(a,candidate) result(same)
class(ngroup),intent(in) :: a,candidate
real(dp)                 :: ar,br

same=.false.

if(.not.a%has_rcut) return
if(.not.candidate%has_rcut) return
if(.not.associated(a%candidate_source)) return
if(.not.associated(candidate%candidate_source)) return
if(a%candidate_source%id/=candidate%candidate_source%id) return

ar=a%candidate%halo_cutoff()
br=candidate%candidate%halo_cutoff()
if(abs(ar-br)>epsilon(max(abs(ar),abs(br),1._dp))) return

same=.true.

end function same_candidate_policy

end module
