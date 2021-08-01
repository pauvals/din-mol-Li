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

#define SMALL 0.0001

! This module compute the neighbor list of a group of atoms. It also sort the
! atoms into cells to improve the search efficiency. A second group of atoms
! is used to perform the neighbor search. This group of neighbor is
! considered only for reading operations, not for writing. That is, using a
! single loop to change the atoms and its neighbors in the same iteration
! (i.e. add a reaction force to a neighbor), can save a computing factor
! near 2x but is not thread safe. Thus I will ignore this possibility and
! consider that only the atom properties of the first group may change and
! not the atoms of the neighbor group.

use gems_groups
use gems_constants,     only: dp,cdm,dm
use gems_errors
! use gems_algebra

implicit none
private

public  :: update, test_update

logical :: mic=.true.

type(group),target,public    :: ghost

! System/Local atoms, index to allow reference by `tag`.
type(igroup),target,public  :: sys
   
! Box
real(dp),public,dimension(dm),target    :: box    =1.0e6_dp  
real(dp),public,dimension(dm)    :: box_old=0.0_dp      ! For pbcghost. Set this small to force initial ghost inclusion
logical                          :: boxed=.true.

integer,parameter,public,dimension(26,3) :: n1cells = transpose(reshape( &
                   [ 1, 0, 0, -1, 0, 0,  0, 1, 0,&
                     1, 1, 0, -1, 1, 0,  0,-1, 0,&
                     1,-1, 0, -1,-1, 0,  0, 0, 1,&
                     1, 0, 1, -1, 0, 1,  0, 1, 1,&
                     1, 1, 1, -1, 1, 1,  0,-1, 1,&
                     1,-1, 1, -1,-1, 1,  0, 0,-1,&
                     1, 0,-1, -1, 0,-1,  0, 1,-1,&
                     1, 1,-1, -1, 1,-1,  0,-1,-1,&
                     1,-1,-1, -1,-1,-1 ],[3,26]))
! integer,public,dimension(26,3) :: n1test = 0.5_dp*(sign(n1cells(:,:))+1)-abs(n1cells(:,:))
                                              
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

  ! Flag if linked cells is made
  logical   :: cells=.false.

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

  ! Protocol to set cells
  ! TODO: Build different update protocols. e.g update cells if box change a
  ! lot or without box but using min and max coordinates of the group.
  procedure(cgroup_first),pointer :: update => cgroup_first

  contains

    procedure :: dest => cgroup_destroy
    procedure :: grow => cgroup_grow

    procedure :: setrc => cgroup_setrc
    procedure :: sort => cgroup_sort

end type

! Group with neighbors
! ngroup holds all atoms in a interaction.
! A general interaction is of type ngroup%a<ngroup%b, i.e. forces are
! computed in group a given group b. The intrinsic ngroup gives an index to
! give an index to each atom.
type, extends(igroup), public :: ngroup
                       
  ! The reference group
  type(group)   :: ref
                      
  ! The group of possible neighbors. It might be sorted in cells.
  type(cgroup)  :: b

  ! rcut
  real(dp)               :: rcut=1.e10_dp,rcut2=1.e10_dp

  ! Maximum number of neighbor per atom
  integer                :: mnb=1000

  ! Neighbor list build in this module and used in any forcefield routine.
  integer,allocatable    :: nn(:)        ! Numero de vecinos al atomo i
  integer,allocatable    :: list(:,:)    ! Vecino m-th del i-th atom of the first iteracting group (local and ghost a)

  ! This boolean allows to skip all the interaction.
  ! From the interact routine. This allows to compute
  ! forces on selected interactions
  logical :: disable=.false.

  ! Always trying to use linked cell over verlet list
  procedure(ngroup_cells),pointer :: lista=>ngroup_cells

  contains

    procedure :: init => ngroup_construct
    procedure :: dest => ngroup_destroy
    procedure :: grow => ngroup_grow

    ! Keep access to procedures of abstract parent
    ! see (https://fortran-lang.discourse.group/t/call-overridden-procedure-of-abstract-parent-type/590/23)
    procedure :: init_abstract => ngroup_construct
    procedure :: dest_abstract => ngroup_destroy
    procedure :: grow_abstract => ngroup_grow

    procedure :: setrc => ngroup_setrc

endtype

abstract interface
 subroutine ngroup0(g)
  import ngroup
  class(ngroup),intent(inout)  :: g
 end subroutine
end interface

! Double linked list for ngroups used for modules that require to keep track of their ngroups
! (e.g. TB or EAM to perform the preinteraction)
#define SOFT
#define _NODE ngroup_dl
#define _CLASS class(ngroup)
#include "dlist_header.inc"

#define _NODE ngroup_aop
#define _CLASS class(ngroup)
#include "arrayofptrs_header.inc"

#define _NODE ngroup_vop
#define _TYPE type(ngroup_aop)
#include "vector_header.inc"

! VOP to collect the ngroups declared inside modules
! The order of execution is important! See
! bias interaction.
type(ngroup_vop),public :: ngindex

! Global variables
real(dp),target     :: nb_dcut=1._dp        ! The shell length for verlet update criteria
integer             :: nupd_vlist = 0       ! Number of updates to print in the log
real(dp)            :: maxrcut=0.0_dp       ! Maximum cut ratio

! Ghost atoms. No need index.
logical,public              :: useghost=.false.
logical,public              :: fullghost=.false.

contains

#define SOFT
#define _NODE ngroup_dl
#define _CLASS class(ngroup)
#include "dlist_body.inc"

#define _NODE ngroup_vop
#define _TYPE type(ngroup_aop)
#include "vector_body.inc"

! cgroup events
! =============

subroutine cgroup_destroy(g)
class(cgroup)         :: g
call g%igroup%dest()
if(allocated(g%head)) deallocate(g%head)
if(allocated(g%next)) deallocate(g%next)
end subroutine cgroup_destroy

subroutine cgroup_grow(g)
class(cgroup)         :: g
integer               :: n

call g%igroup%grow()

if(allocated(g%next)) then
  if(g%nat<size(g%next)) return
  deallocate(g%next)
endif

n=g%nat+g%pad
allocate(g%next(n))

end subroutine cgroup_grow

! Set properties
! --------------

subroutine cgroup_setrc (g,rc)
class(cgroup)       :: g
real(dp),intent(in) :: rc

! New cut radious
g%rcut = rc
g%rcut2 = rc*rc
maxrcut=max(maxrcut,rc)

end subroutine cgroup_setrc

subroutine cgroup_first(g)
! FIXME: Instead of print, use error handling, I can merge it with
! cgroup_update()
class(cgroup) :: g
integer,parameter :: mincells=60 ! (4x4x4)

g%cells = .false.

! If not a box defined, do not use linked cells
if(.not.boxed) then
  call wlog('NHB','Not using linked cells since box is not defined.')
  g%update => cgroup_update
  return
endif

! If rcut not defined, do not use linked cells
if(g%rcut==1.e10_dp) then
  call wlog('NHB','Not using linked cells since interaction does not have rcut.')
  g%update => cgroup_update
  return
endif
   
! Redimension of ncells
g%ncells(:)=int(box(:)/(g%rcut+nb_dcut)) ! Number of cells in each direction

! Less than 4 have no sense to use cells
if(all(g%ncells(:)<4)) then
  call wlog('NHB','Not using linked cells since it would involve less than 4 cells.')
  call wlog('NHB'); write(logunit,fmt="(a,f10.5)") ' -cut radious: ', g%rcut
  g%update => cgroup_update
  return
endif

! Flag that cells where build
g%cells = .true.

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

call wlog('NHB','Using Linked Cells.')
call wlog('NHB'); write(logunit,fmt="(a,f10.5)") ' cut radious: ', g%rcut
call wlog('NHB'); write(logunit,fmt="(a,"//cdm//"(f10.5,2x))") ' cell size: ', g%cell(1:dm)
call wlog('NHB'); write(logunit,fmt="(a,"//cdm//"(i3,2x))") ' cell numbers: ', g%ncells(1:dm)

! This will run only once
g%update => cgroup_update

end subroutine cgroup_first

subroutine cgroup_update(g)
! TODO: Instead of print, use error handling.
class(cgroup) :: g
real(dp)          :: rcut
integer,parameter :: mincells=60 ! (4x4x4)

rcut=g%rcut+nb_dcut

! Si ya viene configurado, chequear si hace falta cambios y si no actualizar
! el tamaño de las celdas solamente
! Esto evita alocatear head t odo el tiempo 
if(g%cells) then
  ! Check if it is not possible to include more cells
  if(any(box(:)-g%ncells(:)*rcut<=rcut)) then
       
    ! Check if the cell size do not drop below cut radious
    if(any(g%cell(:)>rcut)) then
                
      ! Update cell size (needed for NPT)
      g%cell(:)=box(:)/g%ncells(:)
      return

    endif
         
  endif
endif
             
g%cells = .false.

! If not a box defined, do not use linked cells
if(.not.boxed) return
 
!TODO: Cuando el sistema no esta en una caja.. se podría hacer algo asi.
! Pero para ello hay que primero generalizar GEMS para que la caja no tenga
! el punto minimo en (0,0,0).
! if(.not.boxed) then
!   call get_boundingbox(g%alist,box_max,box_min,g%n(2))
!   call get_boundingbox(g%blist,amax,amin,g%n(4)-g%n(2))
!   box_max(:)=max(box_max(:),amax(:))
!   box_min(:)=min(box_min(:),amin(:))
!   set box
! endif
 
! Redimension of ncells
g%ncells(:)=int(box(:)/(g%rcut+nb_dcut)) ! Number of cells in each direction

! Less than 4 have no sense to use cells
if(all(g%ncells(:)<4))  return

! Flag that cells where build
g%cells = .true.

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

end subroutine cgroup_update

subroutine cgroup_sort(g)
! Sort g atoms into the linked list of each cell.
! Esta subrrutina debe llamarse siempre antes del calculo de las fuerzas Las
! variables extrañas son seteadas en la subrutina maps.
class(cgroup),intent(inout)  :: g
integer                      :: i,aux1(dm)
real(dp)                     :: one_cell(3)
type(atom),pointer           :: a

g%head(:,:,:) = 0

! Inverse of cell size
one_cell(:)=1._dp/g%cell(:)

! Atoms of list a
! !$OMP  PARALLEL DO DEFAULT(NONE) &
! !$OMP& PRIVATE(i,a,aux1) &
! !$OMP& SHARED(g,one_cell)
do i = 1,g%nat
  a => g%a(i)%o

  ! FIXME
  ! ! Si no hago pbc podria salirse alguna fuera y dar un segfull
  ! where(a%pbc(:))
  !   a%pos(:)=mod(a%pos(:)+box(:),box(:))
  ! endwhere

  ! Get cell index
  aux1(:)=int(a%pos(:)*one_cell(:))+1
  ! j = 1 + dot_product(aux1,cellaux)

  ! Build cell list
  ! !$OMP CRITICAL
  g%next(i) = g%head(aux1(1),aux1(2),aux1(3))
  g%head(aux1(1),aux1(2),aux1(3)) = i
  ! !$OMP END CRITICAL

enddo
! !$OMP END PARALLEL DO

end subroutine cgroup_sort

! ngroup events
! =============

subroutine ngroup_construct(g)
class (ngroup),target             :: g

! Initialize
call g%igroup%init()

! Ghost atoms must have its own index in order to comunicate back 
! to the real image the computed properties
g%ghost=.true.

! Index the group
call ngindex%append()
ngindex%o(ngindex%size)%o=>g

! Init internal groups
call g%ref%init()
call g%b%init()
g%b%ghost=.true.

end subroutine ngroup_construct

subroutine ngroup_destroy(g)
class (ngroup)  :: g
integer         :: i

! Update ngindex of `g`
do i=1,ngindex%size
  if(ngindex%o(i)%o%id==g%id) exit
enddo
call ngindex%del(i,1)

! Destroy internal groups  
call g%ref%dest()
call g%b%dest()

! Destroy
call g%dest()
deallocate(g%nn)
deallocate(g%list)

end subroutine ngroup_destroy

subroutine ngroup_grow(g)
class(ngroup)          :: g
integer                :: n

call g%igroup%grow()

! TODO Check if neighbor list is needed

if(allocated(g%nn)) then
  if(g%nat<size(g%nn)) return
  deallocate(g%nn)
  deallocate(g%list)
endif

n=g%nat+g%pad
allocate(g%nn(n))
allocate(g%list(n,g%mnb))

end subroutine ngroup_grow

! Set properties
! --------------

subroutine ngroup_setrc(g,rc)
! I could use the g%b%rcut directly.. but may be confuse...
class(ngroup)       :: g
real(dp),intent(in) :: rc

! New cut radious
g%rcut = rc
g%rcut2 = rc*rc
maxrcut=max(maxrcut,rc)
call g%b%setrc(rc)
 
end subroutine ngroup_setrc

subroutine ngroup_setlista(g)
class(ngroup) :: g

! Check if neighbor list is needed
if(.not.associated(g%lista)) return

if(g%b%cells) then
  g%lista => ngroup_cells
else
  g%lista => ngroup_verlet
endif

end subroutine ngroup_setlista

! Search algorithms
! -----------------

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

! TODO: Otra opcion para evaluar sería
! !$OMP PARALLEL DO SCHEDULE(STATIC,5) PRIVATE(m,vd,rd)
! do i = 1,g%na+g%nag-1
! m = 0
! do j = i+1,g%na+g%nag
!   vd = g%at(j)%o%pos-g%at(i)%o%pos

!$OMP PARALLEL
!$OMP SINGLE
la => g%ref%alist
do ii = 1,g%ref%nat
  la => la%next
  ai => la%o
  i = ai%gid(g%id)

  !OMP: Se reparte la tarea entre los idle threads
  !$OMP TASK DEFAULT(NONE)     &
  !$OMP& FIRSTPRIVATE(ai,i)    &
  !$OMP& SHARED(g,rcut,mic)    &
  !$OMP& PRIVATE(aj,j,vd,rd,m)

  ! Reset number of neighbors
  m=0

  do j = 1, g%b%nat
    aj=>g%b%a(j)%o

    ! Skip autointeraction
    if(associated(aj,target=ai)) cycle

    vd = vdistance(ai,aj,mic)
    rd = dot_product(vd,vd)

    if (rd>rcut) cycle

    ! Add j as neighbor of i.
    m=m+1
    g%list(i,m)=aj%gid(g%id)

  enddo
  g%nn(i)=m

  !$OMP END TASK
enddo

!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL

end subroutine ngroup_verlet

subroutine ngroup_cells(g)
! Build neighbors verlet list over linked cells.
class(ngroup),intent(inout)  :: g
type(atom), pointer          :: ai, aj
integer                      :: i,ii,j,ic,k
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

! Sort in cells
call g%b%sort()

la => g%ref%alist

!Bucle sobre los atomos de g
!$OMP PARALLEL
!$OMP SINGLE 
do ii = 1,g%ref%nat
  la => la%next
  ai => la%o

  !$OMP TASK DEFAULT(NONE)             &
  !$OMP& FIRSTPRIVATE(ai)              &
  !$OMP& SHARED(g,rcut,mic)            &
  !$OMP& PRIVATE(rc,i,nabor,nc,j,aj,k,vd,rd)   

  i = ai%gid(g%id)

  ! Get cell index
  rc(:)=int(ai%pos(:)/g%b%cell(:))+1

  !Bucle sobre las celdas vecinas a rc(:)
  do nabor=0,26
    nc(:)=map(:,nabor)+rc(:)

    ! Check if cells should be considered and apply PBC
    if(.not.cell_pbc(g%b,nc,mic)) cycle

    !Bucle sobre los atomos de la celda vecina
    j = g%b%head(nc(1),nc(2),nc(3))
    do while( j>0 )
      aj=>g%b%a(j)%o
      k = aj%gid(g%id)

      ! Next here to allow cycle
      j = g%b%next(j)

      ! Skip autointeraction
      if(associated(aj,target=ai)) cycle

      vd = vdistance(aj,ai,mic) ! respetar el orden
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


end subroutine ngroup_cells

! Update search
! =============

subroutine update()
! Update all neighbor lists
class(ngroup), pointer     :: g
integer                    :: i

! call system_clock(t1)
nupd_vlist = nupd_vlist +1

do i = 1, sys%nat
  sys%a(i)%o%pos_old = sys%a(i)%o%pos
enddo

! Needed for NPT... naaa
if(boxed) box_old(:)=box(:)

! Simple loop
do i = 1, ngindex%size
  g => ngindex%o(i)%o

  ! Check to skip
  if (g%disable) cycle

  call g%b%update()
  call ngroup_setlista(g)

  if(associated(g%lista)) call g%lista()
enddo

end subroutine update

subroutine test_update()
! Check if neighbor update is needed
real(dp)            :: rd,dispmax1,dispmax2,vd(dm)
integer             :: i

! Update the ghost positions
if(useghost) call pbcghost_move

!vecinos con el propio sistema
dispmax1 = 1.e-16_dp
dispmax2 = 1.e-16_dp

do i = 1, sys%nat

  vd = sys%a(i)%o%pos - sys%a(i)%o%pos_old

  rd = dot_product(vd,vd)

  if (rd>dispmax1) then
    dispmax2 = dispmax1
    dispmax1 = rd
  else
    if (rd>dispmax2) dispmax2 = rd
  endif

enddo

if (sqrt(dispmax1)+sqrt(dispmax2)>nb_dcut) then
  if(useghost) then
    if(fullghost) then
       call pbcfullghost()
    else
       call pbcghost()
    endif
  endif
  call update()
endif

end subroutine test_update

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

! BOX

function vdistance(i,j,mic) result(vd)
!calculates the distance of two atoms with or without minimum image convention
real(dp),dimension(dm)  :: vd
type(atom),intent(in)   :: i,j
logical,intent(in)      :: mic
logical                 :: pbc(dm)
integer                 :: l
  
! Distancia
vd=i%pos-j%pos

! Convencion de imagen minima
if(.not.mic) return

! Mas rapido usar idnint que un if
pbc=i%pbc.or.j%pbc
do l = 1,dm
  if (pbc(l)) vd(l)=vd(l)-box(l)*idnint(vd(l)/box(l))
enddo

end function vdistance

subroutine set_pbc(g,pbc)
! rota un grupo un angulo r sobre el punto v
class(group),intent(inout)  :: g
class(atom_dclist),pointer  :: la
logical,intent(in)         :: pbc(dm)
integer                    :: i

la => g%alist
do i = 1,g%nat
  la => la%next
  la%o%pbc = pbc
enddo

end subroutine set_pbc

subroutine do_pbc(g)
class(group),intent(inout)  :: g
class(atom_dclist),pointer :: la
integer                    :: i,k

la => g%alist
do i = 1,g%nat
  la => la%next

  do k=1,dm
    if (la%o%pbc(k)) then
      if(la%o%pos(k)>=box(k)) then
        la%o%pos(k)=la%o%pos(k)-box(k)
        la%o%boxcr(k)=la%o%boxcr(k)+1
      elseif(la%o%pos(k)<0.0_dp) then
        la%o%pos(k)=la%o%pos(k)+box(k)
        la%o%boxcr(k)=la%o%boxcr(k)-1
      endif
    endif
  enddo
enddo
     
end subroutine do_pbc
         

! GHOSTS

subroutine new_ghost(o2,r)
class(atom),target,intent(in)  :: o2
type(atom),pointer             :: o
class(group),pointer           :: g
integer                        :: j,i
integer                        :: r(:)

! Allocate new ghost
allocate(o)
call o%init()
call ghost%attach(o)

! Assign image properties
o%ghost=>o2
call atom_asign(o,o2)
o%id(:)=o2%id(:)
o%pos(:)=o2%pos(:)+r(:)*box(:)
o%boxcr(:)=r(:)

! Add ghost to the image ngroups
do j=1,o2%ngr
  i = o2%gr(j)
  g => gindex%o(i)%o
  if(g%ghost) call g%attach(o)
enddo

end subroutine

function islocal(r)
! Return .true. if r is inside the box
! It is equivalent to `isghost(r,0.)`
real(dp),intent(in)  :: r(dm)
logical              :: islocal
integer              :: i

islocal=.false.
do i=1,dm
  if(r(i)<0._dp) then
    return
  elseif (r(i)>=box(i)) then
    return
  endif
enddo
islocal=.true.

end function islocal

function isghost(r,rc)
! Return .true. if r inside the box+rcut. So actually it says if some position
! belong to a ghost region but also if is a local region
real(dp),intent(in)  :: r(dm),rc
logical              :: isghost
integer              :: i

isghost=.false.
do i=1,dm
  ! s=sign(1,r(i)-box2(i))
  ! if (s*r(i)>rc+box*0.5_dp*(s+1)) return

  if(r(i)<-rc) then
    return
  elseif (r(i)>rc+box(i)) then
    return
  endif
enddo
isghost=.true.

end function isghost

function isoldghost(r,rc)
! Return .true. if r inside the box+rcut. So actually it says if some position
! belong to a ghost region but also if is a local region
real(dp),intent(in)  :: r(dm),rc
logical              :: isoldghost
integer              :: i

isoldghost=.false.
do i=1,dm
  ! s=sign(1,r(i)-box2(i))
  ! if (s*r(i)>rc+box*0.5_dp*(s+1)) return

  if(r(i)<-rc) then
    return
  elseif (r(i)>rc+box_old(i)) then
    return
  endif
enddo
isoldghost=.true.

end function isoldghost

function replicaisghost(n,r,rc)
! This expresion si general for any n(:) (neigh cells or not) and rc (bigger than
! box or not)
! This is the same that `call isghost(r+n*h,rc)`
real(dp),intent(in)  :: r(dm),rc
integer,intent(in)   :: n(dm)
logical              :: replicaisghost
integer              :: i,s

replicaisghost=.false.
do i=1,dm
  s=sign(1,n(i))
  if (s*r(i)>rc+box(i)*(0.5_dp*(s+1)-abs(n(i)))) return
enddo
replicaisghost=.true.

end function replicaisghost

subroutine pbcghost()
! Update ghost positions in a serial pbc simulation for only certain cut ratios.
! This should be call each time the neighbor list will be updated. This
! routine assume that the local atoms moves "slowly" between different calls. In
! other words,if there is a chance that local configurations used in two
! consecutive calls to pbcghost are uncorrelated, it is safer to use
! pbcfullghost.
use gems_groups, only: atom_dclist
real(dp)                   :: rcut,r(dm),rold(dm)
type(atom), pointer        :: o,o2
type(atom_dclist), pointer :: la
integer                    :: i,j,m
logical                    :: updatebcr

! if(ghost%nat==0) return

rcut=maxrcut+nb_dcut

! Check ghost status
updatebcr=.false.
la => ghost%alist
do i = 1,ghost%nat
  la => la%next
  o => la%o

  ! Leaved box+rcut, not a ghost any more
  if(.not.isghost(o%pos,rcut)) then
    call o%dest()
    cycle
  endif

  ! Ghost becomes local, swaps with its image.
  if(islocal(o%pos)) then

    updatebcr=.true.

    ! Fold its local image inside the box
    o2=>o%ghost
    do j=1,dm
      if(o2%pos(j)>=box(j)) then
        o2%pos(j)=o2%pos(j)-box(j)
        o2%pos_old(j)=o2%pos_old(j)-box(j)
        o2%boxcr(j)=o2%boxcr(j)+1
      elseif(o2%pos(j)<0.0_dp) then
        o2%pos(j)=o2%pos(j)+box(j)
        o2%pos_old(j)=o2%pos_old(j)+box(j)
        o2%boxcr(j)=o2%boxcr(j)-1
      endif
    enddo

    ! Pass the ghost atom to the opposite cell
    o%pos(:)=o2%pos(:)-o%boxcr(:)*box(:)

  endif

enddo

! Find new ghosts
! !$OMP PARALLEL DO PRIVATE(m,i,j,la,o,k,o2,g,r,rold)
do m =1,26

  do i = 1,sys%nat
    o2 => sys%a(i)%o

    ! Image position..
    ! TODO: it would be easy to check proximity to the border?
    r(:)=o2%pos(:)+n1cells(m,:)*box(:)

    ! if (.not.all(la%o%pbc(:)*n1cells(m,:))) cycle

    if (isghost(r(:),rcut)) then

      ! Check if this was a ghost before
      rold(:)=o2%pos_old(:)+n1cells(m,:)*box_old(:)
      if (.not.isoldghost(rold(:),rcut)) then

        ! !$OMP CRITICAL
        call new_ghost(o2,n1cells(m,:))
        ! !$OMP END CRITICAL

      endif
    end if
  enddo

enddo
! !$OMP END PARALLEL DO

! Actualizar boxcr ya que si un local es plegado hacia adentro de la caja,
! muchos fantasmas van a tener el o%boxcr(:) incorrecto. Abria que ver el
! virial... FIXME
la => ghost%alist
do i = 1,ghost%nat
  la => la%next
  o => la%o
  if (updatebcr) o%boxcr(:)=floor(o%pos(:)/box(:))
  ! Save positions
  ! o%pos_old(:)=o%pos(:)
enddo

end subroutine

subroutine pbcghost_move()
! Move ghost atoms to reflect the motion of their local images
use gems_groups, only: atom_dclist
type(atom_dclist), pointer :: la
type(atom), pointer        :: o
integer                    :: i

la => ghost%alist
do i = 1,ghost%nat
  la => la%next
  o  => la%o%ghost
  la%o%pos(:)=o%pos(:)+box(:)*la%o%boxcr(:)
enddo

end subroutine


! subroutine ghost_pbc
!  use gems_program_types, only:box,atom_dclist,nlocal,nghost,aghost,alocal
!  integer                 :: ii,i,j
!  type ( atom_dclist ), pointer :: la,lb,lc
!  type ( atom ), pointer :: o
!
!   la => aghost
!   do ii = 1,nghost
!     la => la%next
!
!     if(all(la%o%pos(:)>0.0_dp)) then
!     if(all(la%o%pos(:)<box(:))) then
!       j = la%o%tag
!       i = la%o%id
!       lb=>a(j)%o%link
!
!       ! Ajusto alist y aghost
!       o=>la%o
!       la%o=>lb%o
!       lb%o=>o
!
!       ! ! Ajusto a
!       ! a(i)%o=>a(j)%o
!       ! a(j)%o=>o
!       ! a(i)%id=i
!       ! a(j)%id=j
!
!       ! Ajusto los links
!       lc=>lb%o%link
!       lb%o%link=>la%o%link
!       la%o%link=>lc
!
!       ! With half list, changeing the id is require a full update
!     endif
!     endif
!   enddo
!
! end subroutine

subroutine pbcfullghost()
! Update ghost positions in a serial pbc simulation for only certain cut ratios.
! This should be call each time the neighbor list will be updated. I think this
! routine is prefered over pbcghost when it can not be assume that the local
! atoms moves "slowly" between different calls. In other words, if there is a
! chance that local configurations used in two consecutive calls to pbcghost are
! uncorrelated, it is safer to use pbcfullghost.
use gems_groups, only: atom_dclist
real(dp)            :: rcut,r(dm)
type(atom),pointer  :: o
integer             :: i,m

rcut=maxrcut+nb_dcut

! Destroy all ghost atoms
! FIXME: avoid this deallocate.
call ghost%destroy_all()

! Make local atoms pbc
call do_pbc(sys)

! Find new ghost atoms
! !$OMP PARALLEL DO PRIVATE(m,i,j,la,o,k,o,g,r,rold)
do m =1,26

  do i = 1,sys%nat
    o => sys%a(i)%o

    ! if (.not.all(la%o%pbc(:)*n1cells(m,:))) cycle

    r(:)=o%pos(:)+n1cells(m,:)*box(:)
    if (isghost(r(:),rcut)) then

      ! !$OMP CRITICAL
      call new_ghost(o,n1cells(m,:))
      ! !$OMP END CRITICAL

    end if
  enddo

enddo
! !$OMP END PARALLEL DO

end subroutine

function idoit(i,j,itag,jtag,vd)
! Use precomputed midpoint criterion to decide if interaction is owned.
logical                :: idoit
integer,intent(in)     :: i, j, itag, jtag
real(dp),intent(in)    :: vd(3)

idoit = .false.

if (i>sys%nat) return

if (j<sys%nat) then
  idoit = .true.
else if (itag < jtag) then
  idoit = .true.
else if (itag == jtag) then
  if (vd(3) > SMALL) then
    idoit = .true.
  else if (abs(vd(3)) < SMALL) then
    if (vd(2) > SMALL) then
      idoit = .true.
    else if (abs(vd(2)) < SMALL .and. vd(1) > SMALL) then
      idoit = .true.
    endif
  endif
endif

end function idoit

end module
