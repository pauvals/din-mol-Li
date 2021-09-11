! Copyright (c) 2020  Sergio Alexis Paz
!
!  This file is part of GEMS. GEMS is an Extensible Molecular Simulator.
!	 .
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


module gems_groups
use gems_constants,only:sp,dp,dm

implicit none
private

! Atom type
! =========

type, public :: atom

  ! Group membership
  ! ----------------
              
  ! The number of the groups that holds the atom
  integer                :: ngr=0
              
  ! The ids of the groups that holds the atom
  integer,allocatable    :: gr(:)

  ! The atom id for each igroup (regular groups has a 0)
  integer,allocatable    :: id(:)

  ! If the aotm is a ghost, point to the real image
  class(atom),pointer    :: ghost=>null()
                         
  ! Element properties
  ! ------------------

  ! Propiedades que defino afuera de e para que se mas rapidamente accedida
  ! (en general la 1/masa esta en los cuellos de botella de los algoritmos)
  integer           :: tipo=0
  real(dp)          :: m=1.0_dp
  character(2)      :: sym

  ! ----- Propiedades mecanicas
  real(dp) :: pos(3),   &!propieties of atom. [a][..][m/s][..]
              force(3), &
              acel(3),  & !aceleracion
              vel(3)

  !In a local atom it has the info to unwrap coordinates. In the ghost atom,
  !it has the info of the subdomain/processor it belongs.
  integer                :: boxcr(dm)=0
  logical                :: pbc(dm)=.false. !PBC para ese atomo

  !para ver el desplazamiento en la lista de vecinos. Esto lo establezco bien
  !grande para forzar la primera actualizacion del verlet
  real(dp),dimension(dm) :: pos_old =1.e8_dp, old_cg = 1.e8_dp

  real(dp)               :: energy=0.d0
         
  contains

  procedure :: init => atom_allocate
  procedure :: dest => atom_destroy

  procedure :: setz => atom_setelmnt

  !procedure :: del => atom_atom_del
  !procedure :: delall => atom_allatom_del
  !procedure :: addatom => atom_include
  !procedure :: atom_asign
  !generic   :: assignment(=) => atom_asign

  ! Group membership
  ! TODO: This should be an atom extenison defined in Group module. 
  ! See del_atom in Group module.
  procedure :: addgr => atom_addgr
  procedure :: delgr => atom_delgr
  procedure :: index => atom_index
  procedure :: gid => atom_id
  procedure :: try_dest => atom_destroy_attempt

end type atom

! Double Circular Linked List of atoms
#define _NODE atom_dclist
#define _CLASS class(atom)
#include "cdlist_header.inc"
public :: atom_dclist
                                  
! Array of Pointers to atoms
#define _NODE atom_ap
#define _CLASS class(atom)
#include "arrayofptrs_header.inc"
public :: atom_ap


! Module procedures 

public :: atom_asign
                   

! Group type
! ==========

! All groups must be declared with target atribute so they can be indexed in
! gindex

type, public :: group

  ! Group ID
  ! ---------

  ! All groups are indexed using the `group_vop` `gindex`.
  ! This ID is particularly useful to quickly know the group memberships of
  ! an atom. It can also be made using pointers... but I think an ID is
  ! clearer, and is the same approach taken with atom IDs, which has further
  ! advantages. See `id` properti inside `atom` type. 
  integer                :: id=0
 
  ! Flag this group to include ghost atoms. If a ghost is created, it will be
  ! included in any group that hold its image and have this flag true.
  logical                :: ghost=.false.
                       
  ! Atom members
  ! ------------
  integer                     :: nat=0
  type(atom_dclist),pointer   :: alist=>null()
 
  ! TODO: DELETE
  ! hypervectores: vectores en el espacio de las fases
  !------------------------------------------------------------------------------
  ! Ciertas subrutinas (i.e. lbfgs) necesitan de un vector del espacio de las
  ! fases (velocidad, posicion, etc). Otras (i.e. neb, mad, claculo de
  ! correlaciones ) necesitan en una acumulacion en la memoria de estos vectores.
  ! De alli estas declaraciones. La forma de asociar estas variables con los
  ! objetos correspondientes, es recurrir a las subrutinas de asociacion para
  ! estos vectores al final de este modulo :
  real(dp),pointer            :: pp(:)=>null(),pv(:)=>null(),pf(:)=>null(),pa(:)=>null() ! vectores actuales

  !Propiedades Mecanicas
  !=====================

  ! Estas propiedades deben calcularse usando el modulo de inquire_properties.
  ! Sus definiciones estan estrechamente relacionadas con este modulo


  ! -----temperatura
  real(dp)                    :: temp    =0.0_dp  ,&
                                 tempvib =0.0_dp  ,&
                                 temprot =0.0_dp  ,&
                                 mass    =0.0_dp
  ! -----energies
  real(dp)                    :: ekin    =0.0_dp  ,&
                                 epot    =0.0_dp
  ! -----center of mass
  real(dp),dimension(dm)      :: cm_pos  =0.0_dp  ,&
                                 cg_pos  =0.0_dp  ,&
                                 cm_vel  =0.0_dp
  ! -----rg
  real(dp)                    :: rg_pos  =0.0_dp
  ! -----angular properties
  real(dp)                    :: erot    =0.0_dp  ,&
                                 evib    =0.0_dp
  real(dp),dimension(3)       :: ang_mom =0.0_dp  ,&
                                 ang_vel =0.0_dp
  real(dp),dimension(3,3)     :: inercia =0.0_dp

  real(dp),dimension(3,3)     :: covar= 0.0_dp               ! Matriz de covarianza

  ! Reservado para las contribuciones a la energia potencial
  real(dp)                      :: epot1=0.0_dp, &
                                   epot2=0.0_dp, &
                                   epot3=0.0_dp, &
                                   epot4=0.0_dp, &
                                   epot5=0.0_dp


  ! ----- pressure
  real(dp),dimension(3,3)     :: virial=0.0_dp,  &
                                 pressure=0.0_dp


  ! -----geometria y morfologia
  real(dp),dimension(dm)      :: maxpos  =0.0_dp  ,&
                                 minpos  =0.0_dp
  real(sp),dimension(3)       :: mainaxis =0.0_sp

  ! Para cada propiedad le asigno un booleano que me indica si esta necesita
  ! ser calculada. Por ejemplo, me ayudaria para saber si tengo que recalcular
  ! la masa cuando estoy calculando el centro de masa o no. La idea es que el
  ! inquire sea mas eficientes

  ! -----temperatura
  logical                     :: b_temp    =.false. ,&
                                 b_tempvib =.false. ,&
                                 b_temprot =.false. ,&
                                 b_mass    =.false.
  ! -----energies
  logical                     :: b_ekin    =.false. ,&
                                 b_epot    =.false.
  ! ----- pressure
  logical                     :: b_virial   =.false. ,&
                                 b_pressure =.false.

  ! -----center of mass
  logical                     :: b_cm_pos  =.false. ,&
                                 b_cg_pos  =.false. ,&
                                 b_cm_vel  =.false.
  ! -----radio de giro
  logical                     :: b_rg_pos  =.false.

  ! -----angular properties
  logical                     :: b_erot    =.false. ,&
                                 b_evib    =.false.
  logical                     :: b_ang_mom =.false. ,&
                                 b_ang_vel =.false.
  logical                     :: b_inercia =.false. ,&
                                 b_covar   =.false.
  ! -----geometria y morfologia
  logical                     :: b_maxpos  =.false. ,&
                                 b_minpos  =.false. ,&
                                 b_mainaxis=.false.

  ! Esta indice permite saber si el grupo es usado en variables colectivas para
  ! evitar que sea borrado dejando colgado una variable colectiva.
  integer                     :: cvs=0

  contains

  procedure :: init => group_construct
  
  ! TODO: Make this FINAL
  procedure :: dest => group_destroy

  !procedure :: delall => group_all_detach
  procedure :: detach_atom => group_atom_detach
  procedure :: detach_all => group_all_detach
  generic   :: detach => detach_all, detach_atom

  procedure :: try_destroy_all => group_all_destroy_attempt
  procedure :: destroy_all => group_all_destroy

  procedure :: attach_atom => group_attach_atom
  procedure :: attach_group => group_attach_group
  generic   :: attach => attach_group,attach_atom

  ! devuelve un puntero correspondiendo a un indice desde el head
  ! XXX: Note that igroup might be a better resource
  procedure :: atom => group_atombyindex 

end type group
 
! An indexed group with atom pointers sorted in an array. A second index
! given by `id` is found useful to establish a link to other indexed groups.
type, extends(group), public :: igroup

  ! The array of pointers.
  type(atom_ap),allocatable :: a(:)
    
  ! The growth speed for reallocation
  integer                 :: pad=100

  contains
  procedure :: init => igroup_construct
  procedure :: dest => igroup_destroy

  procedure :: grow => igroup_grow

  procedure :: detach_atom => igroup_atom_detach
  procedure :: attach_atom => igroup_attach_atom
                                 
end type igroup

! Array of Pointers to groups
#define _NODE group_ap
#define _CLASS class(group)
#include "arrayofptrs_header.inc"
public :: group_ap
  
#define _NODE group_vop
#define _TYPE type(group_ap)
#include "vector_header.inc"


! Basic Groups
! ------------

! Index of groups (see `id` in group type)
type(group_vop),target,public   :: gindex

contains
      
! atom events
! ===========
 
#define _NODE atom_dclist
#define _CLASS class(atom)
#include "cdlist_body.inc"
 
subroutine atom_allocate(a)
! inicializo los punteros y allocateables que no
! se pueden inicializar en la declaración
class(atom),intent(inout)    :: a

! crear un atomo
a%acel(:)=0._dp
a%pos(:)=0._dp
a%vel(:)=0._dp
a%force(:)=0._dp
call a%setz('F')

end subroutine atom_allocate

subroutine atom_destroy(a)
class(atom)         :: a
integer             :: i,j

! Dettach the atom from all the groups
do i=1,a%ngr
  j=a%gr(i)
  call gindex%o(j)%o%detach(a)
enddo
if (allocated(a%gr)) deallocate(a%gr,a%id)

end subroutine atom_destroy
 
subroutine atom_destroy_attempt(a)
class(atom)         :: a

! Skip detroy if there is a group reference to this atom
if(a%ngr>0) return

! Destroy
call a%dest

end subroutine atom_destroy_attempt
 
subroutine atom_asign(a1,a2)
! Copia la informacion de a2 en a1. Esto se hace sin importar como estan !
! conformados los objetos, pudiendo por ejemplo a2%pos ser un puntero slice
! o un arreglo duro, no importa.
type(atom) :: a1,a2

a1%m        = a2%m
a1%pos(:)   = a2%pos(:)
a1%vel(:)   = a2%vel(:)
a1%force(:) = a2%force(:)
a1%acel(:)  = a2%acel(:)
call a1%setz(a2%sym)

a1 % pos_old(:) = a2 % pos_old(:)

a1 % energy    = a2 % energy

end subroutine atom_asign
  
! group membership
! ----------------
 
subroutine atom_addgr(a,uid)
! Register group uid into internal atom records  
class(atom)         :: a
integer,intent(in)  :: uid
integer,allocatable :: t_id(:),t_gr(:)
integer             :: n  

n=a%ngr

if(allocated(a%gr)) then
  if(n<size(a%gr)) then
    a%ngr=n+1 
    a%gr(n+1)=uid
    a%id(n+1)=0
    return
  endif
endif

allocate(t_id(n+5))
allocate(t_gr(n+5))
t_gr(1:n) = a%gr(1:n)
t_id(1:n) = a%id(1:n)
call move_alloc(to=a%id,from=t_id)
call move_alloc(to=a%gr,from=t_gr)

a%ngr=n+1
a%gr(n+1)=uid
a%id(n+1)=0
 
end subroutine atom_addgr

subroutine atom_delgr(a,uid)
! Unregister group from atom records
class(atom)         :: a
integer,intent(in)  :: uid
integer             :: i,j

j=0
do i=1,a%ngr
  if (j/=0) then
    a%gr(i-j)=a%gr(i)
    a%id(i-j)=a%id(i)
  endif
  if (a%gr(i)==uid) j=j+1
enddo
a%ngr=a%ngr-j
 
end subroutine atom_delgr

function atom_index(a,uid) result(id)
! Return 0 if atom do not belong to the group (given by uid)
!        i if atom belongs and group is registred in the ith element
class(atom)         :: a
integer,intent(in)  :: uid
integer             :: i,id

id=0 
do i=1,a%ngr
  if (a%gr(i)==uid) then
    id=i
    return
  endif
enddo
 
end function atom_index
               
function atom_id(a,uid) result(id)
! Return -1 if atom do not belong to this group (given by uid)
!        0  if atom do not have an id for this group
!        id if atom belongs and has id in the group
class(atom)         :: a
integer,intent(in)  :: uid
integer             :: i,id

id=0 
do i=1,a%ngr
  if (a%gr(i)==uid) then
    id=a%id(i)
    return
  endif
enddo
 
end function atom_id
               
! others
! -------

subroutine atom_setelmnt(a,z)
! Establece las propiedades relacionadas al elemento en un atomo
class(atom)               :: a
character(*),intent(in)   :: z


select case(z)
case('Li')
 a%sym='Li'
 a%tipo=1
case('CG')
 a%sym='CG'
 a%tipo=2
case('F')
 a%sym='F'
 a%tipo=3
case default
end select

end subroutine

! group events
! ============
               
#define _NODE group_vop
#define _TYPE type(group_ap)
#include "vector_body.inc"
       
subroutine group_construct(g)
class(group),target :: g

! Init atom list
allocate(g%alist)
call g%alist%init()
g%nat = 0

! Init head para acciones agrupadas
allocate(g%alist%o)
call g%alist%o%init()

! Index the group
call gindex%append()
g%id=gindex%size
gindex%o(g%id)%o=>g

end subroutine group_construct

subroutine group_destroy ( g )
class(group)    :: g
integer         :: i

! Deindex the group
do i=g%id,gindex%size-1
  gindex%o(i)%o=>gindex%o(i+1)%o
  gindex%o(i)%o%id=gindex%o(i)%o%id-1
enddo
 
! Remove atoms
call g%detach_all()

end subroutine group_destroy

! Include atoms
! -------------
 
subroutine group_attach_atom(g,a)
! Add a `soft atom`, i.e. a new link to atom a
class(group)          :: g
class(atom),target    :: a

! Check if this atom is already in group g
if(a%index(g%id)>0) return

! Add atom to group `alist`
call g%alist%add_before()
call g%alist%prev%point(a)
g%nat = g%nat + 1 ! numero de particulas

! Add group id to atom `gr` 
call a%addgr(g%id)

! propiedades basicas para modificar
g%mass = g%mass + a % m ! masa

! si esta vectorial, ahora no tiene sentido
if(associated(g%pp)) then
  deallocate(g%pp)
  deallocate(g%pv)
  deallocate(g%pa)
  deallocate(g%pf)
endif

end subroutine group_attach_atom

subroutine group_attach_group(g,g1)
! agrega los atomos del g1 al g2
class(group)              :: g
class(group)              :: g1
type(atom_dclist),pointer :: la
integer                   :: i

la => g1%alist
do i = 1,g1%nat
  la => la%next
  call g%attach(la%o)
enddo

end subroutine group_attach_group

! Remove atoms
! ------------

subroutine group_atom_detach(g,a)
! Detach atom from group `alist`
class(group)               :: g
type(atom_dclist), pointer :: la
class(atom),target         :: a
integer                    :: i

la => g%alist
do i =1,g%nat
  la => la%next

  ! Find the link to the atom
  if (.not.associated(la%o,target=a)) cycle
 
  ! Delete group from atom register
  call la%o%delgr(g%id)
                    
  ! Delete group id from atom `gr` list
  call la%deattach()
  deallocate(la)
  g%nat = g%nat - 1

  ! TODO: propiedades extras para modificar?
  g%mass = g%mass - a%m

  return
enddo

end subroutine group_atom_detach

subroutine group_all_detach(g)
! Detach all atoms from group
class(group)               :: g
type(atom_dclist), pointer :: la,next

! Circulo por la lista hasta que la vacio
la => g%alist%next
do while(g%nat/=0)
  next => la%next
    
  ! Delete group id from atom `gr` list
  call la%o%delgr(g%id)
    
  ! Deattach the link from the group
  call la%deattach()
  deallocate(la)
  g%nat = g%nat - 1

  la=>next
enddo

! TODO: propiedades extras para modificar?
g%nat = 0
g%mass = 0._dp

end subroutine group_all_detach

subroutine group_all_destroy_attempt(g)
! Detach all atoms from group
class(group)               :: g
type(atom_dclist), pointer :: la,next

! Circulo por la lista hasta que la vacio
la => g%alist
do while(g%nat/=0)
  next => la%next
    
  if (la%o%ngr==1) then
    ! Destroy atom, since `g` is its only group
    call la%o%try_dest()
  else
    ! Delete group id from atom `gr` list
    call la%o%delgr(g%id)
  endif
           
  ! Deattach the link from the group
  call la%deattach()
  deallocate(la)

  la=>next
enddo

! TODO: propiedades extras para modificar?
g%nat = 0
g%mass = 0._dp

end subroutine group_all_destroy_attempt

subroutine group_all_destroy(g)
! Detach all atoms from group
class(group)               :: g
type(atom_dclist), pointer :: la,next

! Circulo por la lista hasta que la vacio
la => g%alist
do while(g%nat/=0)
  next => la%next
    
  ! Destroy atom
  call la%o%dest()
         
  ! Deattach the link from the group
  call la%deattach()
  deallocate(la)

  la=>next
enddo

! TODO: propiedades extras para modificar?
g%nat = 0
g%mass = 0._dp

end subroutine group_all_destroy

! Select atoms
! ------------

function group_atombyindex(this,i) result(at)
class(group)            :: this
type(atom),pointer      :: at
integer,intent(in)      :: i
integer                 :: j
type ( atom_dclist ), pointer :: la

if (i>this%nat .or. i<0) at => null()

la => this%alist
do j=1,i
  la => la%next
enddo
at => la%o

endfunction

! igroup events (indexed atoms)
! =============================

subroutine igroup_construct(g)
class(igroup),target    :: g
call group_construct(g)
allocate(g%a(g%pad))
end subroutine igroup_construct

subroutine igroup_destroy (g)
class(igroup) :: g
call group_destroy(g)
! TODO: No needed, allocatable components are deallocated at type deallocation
deallocate(g%a) 
end subroutine igroup_destroy

! Include atoms
! -------------

subroutine igroup_grow(g)
! Increase size of group array
class(igroup)              :: g
type(atom_ap),allocatable  :: t_a(:)
integer                    :: n  

if(allocated(g%a)) then
  n=size(g%a)
  if(g%nat<n) return
endif

allocate(t_a(g%nat+g%pad))
t_a(1:n) = g%a(1:n)
call move_alloc(to=g%a,from=t_a)
       
end subroutine igroup_grow
  
subroutine igroup_attach_atom(g,a)
class(igroup)        :: g
class(atom),target   :: a
integer              :: i,n  

! Save current atom number
n=g%nat

! Attach
call g%group%attach(a)
call g%grow()

! Index new atom
n=g%nat-n
if(n==1) then
  g%a(g%nat)%o=>a
  i=a%index(g%id)
  a%id(i)=g%nat
endif

end subroutine igroup_attach_atom
           
! Remove atoms
! ------------

subroutine igroup_atom_detach(g,a)
! Remove soft atom (i.e. detach) from `alist` and `a`
class(igroup)              :: g
class(atom),target         :: a
class(atom),pointer        :: aj
integer                    :: i,j,n
 
! Search index of `a`
i=a%index(g%id)
if(i==0) return  
n=a%id(i)

! Update index
do j=n,g%nat-1
  aj=>g%a(j+1)%o
  aj%id(i)=j
  g%a(j)%o=>aj
enddo
g%a(j)%o=>null()
              
! Detach atom
call g%group%detach(a)

end subroutine igroup_atom_detach

end module gems_groups

