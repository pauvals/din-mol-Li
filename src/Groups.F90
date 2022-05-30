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
use gems_errors, only:werr
use gems_constants,only:sp,dp,dm,find_io
use gems_elements,only:elements,element,inq_z

implicit none
private


! Group type
! ==========

! All groups must be declared with target atribute so they can be indexed in
! gindex

type, public :: group

  ! An unique ID to identify a group.
  ! Allows to transform syntax like `associated(o,target=g)` into
  ! `o%id==g%id`, which works regardless the dynamic type of g, useful to
  ! avoid issues like in https://stackoverflow.com/a/69138397/1342186.
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

  ! https://stackoverflow.com/a/69138397/1342186.
  ! Parent procedures that set a polymorphic pointer to an extension
  ! must be call avoiding syntax like `g%group%xxx()` since `g%group` 
  ! has declared and dynamic type group (not polymorphic).
  ! Thus, below I keep accesible those procedures by direct call (i.e.
  ! `g%group_xxx()`) which can pass an extended type as dynamic type.
  procedure :: group_construct    ! Set gindex pointer
  procedure :: group_attach_atom  ! Set a%gr pointer via addgr

  procedure :: init => group_construct
  procedure :: dest => group_destroy ! TODO: Make this FINAL
  procedure :: try_destroy_all => group_all_destroy_attempt
  procedure :: destroy_all => group_all_destroy
    
  procedure :: attach_atom => group_attach_atom
  procedure :: attach_group => group_attach_group
  generic   :: attach => attach_group,attach_atom
   
  procedure :: detach_link => group_detach_link
  procedure :: detach_atom => group_detach_atom
  procedure :: detach_all => group_detach_all
  generic   :: detach => detach_all, detach_atom

  procedure :: epot_changed
  procedure :: pos_changed
  procedure :: vel_changed
  procedure :: all_changed

  ! devuelve un puntero correspondiendo a un indice desde el head
  ! XXX: Note that igroup might be a better resource
  procedure :: atom => group_atombyindex 

  procedure :: write => group_write
  generic :: write(formatted) => write
  ! procedure :: write => group_write
 
end type group
 
! An indexed group with atom pointers sorted in an array. A second index
! given by `id` is found useful to establish a link to other indexed groups.
type, extends(group), public :: igroup

  ! The array of pointers.
  type(atom_ap),allocatable :: a(:)

  ! Efective dimension of `a` (nat<amax<size(a)).
  integer             :: amax=0 
    
  ! A dummy atom used as a flag to skip detached atoms. 
  type(atom),pointer  :: limbo
  logical             :: b_limbo=.false.
  integer             :: nlimbo=0

  ! Update index if 10% of the index is null. 
  ! Set to 0 to avoid updates.
  real(dp)                :: aupd=0.5

  ! Signal out to rebuild internal arrays in extended types.
  ! It becomes `.true.` if there is an index update.
  logical                 :: update=.false. 
    
  ! The growth speed for reallocation
  integer                 :: pad=100

  contains
   
  ! Parent procedures that set a polymorphic pointer to an extension
  ! (see group type construct).
  procedure :: igroup_construct   
  procedure :: igroup_attach_atom 
  procedure :: igroup_detach_atom 
           
  procedure :: init => igroup_construct
  procedure :: dest => igroup_destroy
  procedure :: clean => igroup_clean

  procedure :: attach_atom  => igroup_attach_atom
  procedure :: detach_atom  => igroup_detach_atom
  procedure :: update_index => igroup_update_index
                              
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
public :: gindex_epot_changed, &
          gindex_all_changed, &
          gindex_vel_changed, &
          gindex_pos_changed 

! TODO
public :: group_switch_vectorial
public :: group_switch_objeto
 
type(igroup),target,public  :: sys
                                 
! Atom type
! =========

type, public :: atom

  ! Group membership
  ! ----------------
              
  ! The number of the groups that holds the atom
  integer                :: ngr=0
              
  ! IDs (i.e. index in `gindex) of the groups that contain the atom.
  integer,allocatable    :: gr(:)

  ! The atom ID for each igroup class it belongs (0 for regular groups)
  integer,allocatable    :: id(:)
  
  ! Ghost
  ! -----
       
  ! If the atom is a ghost, prime point to the real image
  type(atom),pointer     :: prime=>null()

  ! If the atom is real, ghost(i) point to its ith ghost (if active).
  type(atom_ap),allocatable :: ghost(:) ! Use allocatable since declaration comes later (F2008)
                         
  ! Element properties
  ! ------------------

  ! Propiedades que defino afuera de e para que se mas rapidamente accedida
  ! (en general la 1/masa esta en los cuellos de botella de los algoritmos)
  integer                  :: z=119 ! The generic element
  real(dp)                 :: mass=1.0_dp
  real(dp)                 :: q=0.0_dp  ! Carga
  real(dp)                 :: s=1.0_dp  ! sigma
  real(dp)                 :: e=0.0_dp  ! epsilon
  character(:),allocatable :: sym
  integer                  :: sp=0      ! Hybridization

  ! !  Enlaces y moleculas.... TOFIX
  ! integer      :: abondid(20)=0  ! el indicie dentro de la molecula de los asociados
  ! integer      :: abonds=0  ! el numero de asociados
  ! integer      :: molid=0   ! el indice de la molecula
  ! integer      :: amolid=0  ! el indice dentro de la molecula

  ! ----- Propiedades mecanicas
  real(dp),pointer       :: pos(:)=>null(),   &!propieties of atom. [a][..][m/s][..]
                            force(:)=>null(), &
                            acel(:)=>null(),  & !aceleracion
                            vel(:)=>null()

  !In a local atom it has the info to unwrap coordinates. In the ghost atom,
  !it has the info of the subdomain/processor it belongs.
  integer                :: boxcr(dm)=0
  logical                :: pbc(dm)=.false. !PBC para ese atomo

  real(dp),dimension(dm) :: acel2  =0._dp,& !derivada primera de la aceleración
                            acel3  =0._dp,& !derivada segunda de la aceleración
                            acel4  =0._dp,& !derivada tercera de la aceleración
                            pos_v  =0._dp,& !posicion relativa al punto v
                            vel_v  =0._dp

  !para ver el desplazamiento en la lista de vecinos. Esto lo establezco bien
  !grande para forzar la primera actualizacion del verlet
  real(dp),dimension(dm) :: pos_old =1.e8_dp, old_cg=1.e8_dp

  real(dp)               :: epot=0.d0,                & !energia potencial total[ev]
                            border=0.d0                 !Orden de Enlace

  contains

  procedure :: init => atom_allocate
  procedure :: dest => atom_destroy

  procedure :: setz => atom_setelmnt_byz
  procedure :: setsym => atom_setelmnt_bysym

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
  procedure :: gid => atom_id
  procedure :: gri => atom_gri
  procedure :: gro => atom_gro
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

interface atom_setelmnt
  module procedure atom_setelmnt_bysym,atom_setelmnt_byz
end interface
public :: atom_setelmnt,atom_asign, vdistance
                   
! Ghosts
! ======
public :: ghost_from_atom, pbcghost_move, pbcghost, pbchalfghost
public :: set_pbc, do_pbc
   
! Groups for selection
type(group),target,public    :: ghost

! Ghost atoms. No need index.
logical,public              :: useghost=.false.
    
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
allocate(a%gr(10))
allocate(a%id(10))
allocate(a%pos(dm))
allocate(a%vel(dm))
allocate(a%force(dm))
allocate(a%acel(dm))
a%acel(:)=0._dp
a%pos(:)=0._dp
a%vel(:)=0._dp
a%force(:)=0._dp
call atom_setelmnt(a,119)
end subroutine atom_allocate

subroutine atom_setpbc(a,pbc)
! Set atom pbc and allocate ghosts  
class(atom),target,intent(inout) :: a
logical,intent(in)               :: pbc(dm)
type(atom),pointer               :: og
integer                          :: i

a%pbc(:)=pbc(:)

! Already allocated, skip
if(allocated(a%ghost)) then

  ! Also deallocate if all false
  if(.not.any(pbc)) then
    do i=1,7  
      og => a%ghost(i)%o
      call og%dest()
      deallocate(og)
    enddo
    deallocate(a%ghost)
  endif  
 
endif  

if(.not.any(pbc)) return

! Allocate ghost 
! TODO: if count(pbc)<3 do not use 7
allocate(a%ghost(7))
do i=1,7  

  allocate(a%ghost(i)%o)
  og => a%ghost(i)%o

  call og%init()
  og%prime=>a
  call atom_asign(og,a)
  og%q=a%q
  og%pbc(:)=.false.
  
enddo

end subroutine atom_setpbc
    
subroutine atom_destroy(a)
class(atom)           :: a
type(atom),pointer    :: og
class(group),pointer  :: g
integer               :: i

! Detach the atom from all the groups using a LIFO scheme. Since some
! attach/deatach precedures may depend on other atom memberships, remeber
! to build these group dependencies considering the present LIFO
! destruction.
do while (a%ngr/=0)
  g => a%gro(a%ngr)
  call g%detach(a)
enddo
deallocate(a%gr,a%id)
 
deallocate(a%pos)
deallocate(a%vel)
deallocate(a%force)
deallocate(a%acel)

a%prime => null()  
if(.not.allocated(a%ghost)) return

do i=1,7  
  og=>a%ghost(i)%o
  do while (og%ngr/=0)
    g => og%gro(og%ngr)
    call g%detach(og)
  enddo
  deallocate(a%ghost(i)%o)
enddo
deallocate(a%ghost)
   
end subroutine atom_destroy
 
function atom_destroy_attempt(a) result(r)
class(atom)   :: a
logical       :: r

r=.false.

! Skip destroy if there is a group reference to this atom
if(a%ngr>0) return

! Destroy
call a%dest()
r=.true.  

end function atom_destroy_attempt
 
subroutine atom_asign(a1,a2)
! Copia la informacion de a2 en a1. Esto se hace sin importar como estan !
! conformados los objetos, pudiendo por ejemplo a2%pos ser un puntero slice
! o un arreglo duro, no importa.
type(atom) :: a1,a2

a1 % pos(:)   = a2 % pos(:)
a1 % vel(:)   = a2 % vel(:)
a1 % force(:) = a2 % force(:)
a1 % acel(:)  = a2 % acel(:)
a1 % pbc(:)   = a2 % pbc(:)
call atom_setelmnt(a1,a2%z)

a1 % acel(:)    = a2 % acel(:)
a1 % acel2(:)   = a2 % acel2(:)
a1 % acel3(:)   = a2 % acel3(:)
a1 % acel4(:)   = a2 % acel4(:)
a1 % pos_old(:) = a2 % pos_old(:)

a1 % epot    = a2 % epot

end subroutine atom_asign

! group membership
! ----------------
 
subroutine atom_addgr(a,g,l,found)
! Register group ID into internal atom gr(:) array.
! Return the index l where a%gr(l)=g%id.
class(atom)                  :: a
class(group),target          :: g
integer,intent(out)          :: l
logical,intent(out)          :: found
integer,allocatable          :: t_id(:),t_gr(:)
integer                      :: n,i

! Exit if group is already added
call binleft(a%gr(:a%ngr),g%id,l,found)
if(found) return

n=a%ngr
if(n==size(a%gr)) then
  allocate(t_gr(n+5))
  allocate(t_id(n+5))
  t_id(1:n) = a%id(1:n)
  t_gr(1:n) = a%gr(1:n)
  call move_alloc(to=a%id,from=t_id)
  call move_alloc(to=a%gr,from=t_gr)
  a%id(n+1) = 0
  a%gr(n+1) = g%id
endif

! Append
! r=a%ngr+1 

! Sorted insertion
l=l+1
do i=a%ngr,l,-1
  a%gr(i+1)=a%gr(i)
  a%id(i+1)=a%id(i)
enddo

! Insert
a%gr(l)=g%id
a%id(l)=0
a%ngr=a%ngr+1 

end subroutine atom_addgr

subroutine atom_delgr(a,g,found)
! Unregister group from atom records
class(atom)         :: a
class(group)        :: g
integer             :: i,j
logical,intent(out) :: found

! j=findloc(a%gr(:a%ngr),g%id,1)
call binleft(a%gr(:a%ngr),g%id,j,found)
if(.not.found) return  

a%ngr=a%ngr-1
do i=j,a%ngr
  a%gr(i)=a%gr(i+1)
  a%id(i)=a%id(i+1)
enddo

end subroutine atom_delgr
 
function atom_gri(a,g) result(i)
! Return  0 if atom do not belong to g
!         i as the index of a%gr vector associated with g
class(atom)         :: a
class(group),target :: g
integer             :: i
logical             :: found

! Exit if group is already added
call binleft(a%gr(:a%ngr),g%id,i,found)
if(.not.found) i=0
                                
! Linear search  
! i=findloc(a%gr(1:a%ngr),g%id,1)

end function atom_gri
 
function atom_gro(a,i) result(g)
class(atom)          :: a
class(group),pointer :: g
integer,intent(in)   :: i
g=>null()
call werr('Atom does not belong to that group',i>a%ngr)
g=>gindex%o(a%gr(i))%o
end function atom_gro
  
function atom_id(a,g) result(i)
! Return -1 if atom do not belong to this group (given by uid)
!        0  if atom do not have an id for this group
!        id if atom belongs and has id in the group
class(atom)         :: a
class(group),target :: g
integer             :: i

i=a%gri(g)
if(i==0) then
  i=-1
else  
  i=a%id(i)
endif
 
end function atom_id
               
! properties
! ----------

subroutine atom_setelmnt_byz(a,z)
! Establece las propiedades relacionadas al elemento en un atomo
class(atom)               :: a
integer                   :: z

! FIXME: What if is not there???

a%z = z
a%mass = elements%o(z)%mass
a%sym = elements%o(z)%sym

end subroutine

subroutine atom_setelmnt_bysym(a,sym)
use gems_elements,only:add_z, inq_z
class(atom)    :: a
character(*)   :: sym

! Adding z just in case is not there
call add_z(sym)

call atom_setelmnt_byz(a,inq_z(sym))
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

! Set defaults  
g%nat = 0

! Init head para acciones agrupadas
! allocate(g%alist%o)
! call g%alist%o%init()

! Index the group
call gindex%append()
g%id=gindex%size
gindex%o(g%id)%o=>g

end subroutine group_construct

subroutine group_destroy ( g )
class(group),target   :: g
integer               :: i

! Deindex the group (does not change group ids)
do i=1,gindex%size
  if(gindex%o(i)%o%id==g%id) exit
enddo
call gindex%del(i,1)
! Another way:
! j=0
! do i=1,gindex%size
!   if (j/=0) gindex%o(i-j)%o=>gindex%o(i)%o
!   if (associated(gindex%o(i)%o,target=g)) j=j+1
! enddo
! gindex%size=gindex%size-j
  
 
! Remove atoms
call g%detach_all()

! call g%alist%o%dest()  
! deallocate(g%alist%o)
deallocate(g%alist)  

end subroutine group_destroy
 
subroutine group_write(g, unit, iotype, v_list, iostat, iomsg)
use gems_strings, only: str
class(group),intent(in)   :: g
integer,intent(in)         :: unit
character(*),intent(in)    :: iotype
integer,intent(in)         :: v_list(:)
integer,intent(out)        :: iostat
character(*),intent(inout) :: iomsg  
character(:),allocatable   :: wfmt
 
! Format descriptor
select case(size(v_list))
case(0) ! default
  wfmt = '(i0)'
case(2)
  wfmt = '(i'//str(v_list(1))//')'
case default
  iostat = 1
  iomsg = 'wrong number of format descriptors'
  return
end select
           
iostat=0
select case(iotype)
case('DT','DTnat')
  write(unit,wfmt,iostat=iostat,iomsg=iomsg)  g%nat
case default
  iostat = 1
  iomsg = 'unexpected iotype'
end select
 
end subroutine
 
! Include atoms
! -------------
 
subroutine group_attach_atom(g,a,l_)
! Attach atom `a` to group `g`.
! Optional: return index l where a%gr(l)=g%id. 
class(group),target          :: g
class(atom),target           :: a
logical                      :: found
integer,intent(out),optional :: l_
integer                      :: l

! Add group to atom or skip if is already there
call a%addgr(g,l,found)
if(found) then
  if(present(l_)) l_=0
  return
else
  if(present(l_)) l_=l
endif
 
! Add atom to group `alist`
call g%alist%add_before()
call g%alist%prev%point(a)
g%nat = g%nat + 1 ! numero de particulas

! propiedades basicas para modificar
g%mass = g%mass + a % mass ! masa

! si esta vectorial, ahora no tiene sentido
if(associated(g%pp)) then
  deallocate(g%pp)
  deallocate(g%pv)
  deallocate(g%pa)
  deallocate(g%pf)
endif

call all_changed(g)

end subroutine group_attach_atom

subroutine group_attach_group(g,g1)
! agrega los atomos del g1 al g2
class(group),target       :: g
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
           
subroutine group_detach_link(g,la)
! Detach link `la` from group `alist`
! It return previous link in `la`
class(group)               :: g
class(atom), pointer       :: a
type(atom_dclist), pointer :: la, prev
logical                    :: found

! Delete group from atom register
a=>la%o
call a%delgr(g,found)
if(.not.found) return  
                  
! Delete group id from atom `gr` list
prev=>la%prev
call la%deattach()
deallocate(la)
la=>prev

! TODO: propiedades extras para modificar?
g%nat = g%nat - 1
g%mass = g%mass - a%mass

call all_changed(g)

end subroutine group_detach_link
           
subroutine group_detach_atom(g,a)
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
  call g%detach_link(la)
                    
  return
enddo

! TODO: this might be handle in a higher level?
! Also detach ghosts
if(allocated(a%ghost)) then
  do i=1,7  
    call g%detach(a%ghost(i)%o)
  enddo
endif  
 
end subroutine group_detach_atom

subroutine group_detach_all(g)
! Detach all atoms from `g`
class(group)        :: g
type(atom), pointer :: a

do while(g%nat/=0)
  a => g%alist%next%o
  call g%detach_atom(a)
enddo

end subroutine group_detach_all

subroutine group_all_destroy_attempt(g)
! Detach all atoms from `g` and destroy free atoms.
class(group)        :: g
type(atom), pointer :: a

do while(g%nat/=0)
  a => g%alist%next%o
  call g%detach_atom(a)
 
  if (a%ngr==0) then
    call a%dest()
    deallocate(a)
  endif
   
enddo

end subroutine group_all_destroy_attempt

subroutine group_all_destroy(g)
! Destroy all atoms from group
class(group)               :: g
type(atom), pointer        :: a

do while(g%nat/=0)
  a => g%alist%next%o
  ! call g%detach_atom(a)
  call a%dest()
  deallocate(a)
enddo
   
end subroutine group_all_destroy

! Group Update Properties 
! -----------------------
! A change in 1 group implies a potential change in all of them: if an atom
! changes, all groups might include that atom. Thus, all this update system
! is only useful while nothing change.
! TODO: I should build/check connections between groups to selectively update?
                 
subroutine gindex_all_changed()
integer :: i  
do i=1,gindex%size
  call gindex%o(i)%o%all_changed()
enddo
end subroutine gindex_all_changed
                  
subroutine gindex_pos_changed()
integer :: i  
do i=1,gindex%size
  call gindex%o(i)%o%pos_changed()
enddo
end subroutine gindex_pos_changed
                   
subroutine gindex_vel_changed()
integer :: i  
do i=1,gindex%size
  call gindex%o(i)%o%vel_changed()
enddo
end subroutine gindex_vel_changed
                    
subroutine gindex_epot_changed()
integer :: i  
do i=1,gindex%size
  call gindex%o(i)%o%epot_changed()
enddo
end subroutine gindex_epot_changed
                  
subroutine all_changed(g)
class(group) :: g
g%b_mass     =.false.
g%b_erot     =.false.
g%b_evib     =.false.
g%b_cm_vel   =.false.
g%b_ang_mom  =.false.
g%b_ang_vel  =.false.
g%b_maxpos   =.false.
g%b_minpos   =.false.
g%b_mainaxis =.false.
g%b_cm_pos   =.false.
g%b_rg_pos   =.false.
g%b_covar    =.false.
g%b_inercia  =.false.
g%b_virial   =.false.
g%b_pressure =.false.
g%b_ekin    =.false.
g%b_temp    =.false.
g%b_tempvib =.false.
g%b_temprot =.false.
g%b_mass     =.false.
call epot_changed(g)
end subroutine all_changed
            
subroutine pos_changed(g)
class(group) :: g
g%b_erot     =.false.
g%b_evib     =.false.
g%b_cm_vel   =.false.
g%b_ang_mom  =.false.
g%b_ang_vel  =.false.
g%b_maxpos   =.false.
g%b_minpos   =.false.
g%b_mainaxis =.false.
g%b_cm_pos   =.false.
g%b_rg_pos   =.false.
g%b_covar    =.false.
g%b_inercia  =.false.
g%b_virial   =.false.
g%b_pressure =.false.
call epot_changed(g)
end subroutine pos_changed

subroutine epot_changed(g)
class(group) :: g
g%b_epot  =.false.
g%b_virial =.false.
end subroutine epot_changed
                          
subroutine vel_changed(g)
class(group) :: g
g%b_ekin    =.false.
g%b_temp    =.false.
g%b_tempvib =.false.
g%b_temprot =.false.
g%b_erot    =.false.
g%b_evib    =.false.
g%b_cm_vel  =.false.
g%b_ang_mom =.false.
g%b_ang_vel =.false.
g%b_pressure =.false.
end subroutine vel_changed
                    
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

! Properties
! ----------

function vdistance(i,j,mic) result(vd)
!calculates the distance of two atoms with or without minimum image convention
use gems_program_types, only: box, one_box
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
  if (pbc(l)) vd(l)=vd(l)-box(l)*idnint(vd(l)*one_box(l))
enddo

end function vdistance
         
 
! igroup events (indexed atoms)
! =============================

subroutine igroup_construct(g)
class(igroup),target    :: g
call group_construct(g)
allocate(g%a(g%pad))
allocate(g%limbo)  
end subroutine igroup_construct

subroutine igroup_destroy(g)
class(igroup),target :: g
call group_destroy(g)
deallocate(g%a) 
deallocate(g%limbo)  
end subroutine igroup_destroy

subroutine igroup_clean(g)
class (igroup),target  :: g
integer                :: i

if(.not.g%b_limbo) return

do i =1,g%amax
  if(associated(g%a(i)%o,target=g%limbo)) then
    g%a(i)%o=>null()
  endif
enddo

g%nlimbo=0
g%b_limbo=.false.

! Update index might be done here.

end subroutine igroup_clean

! Include atoms
! -------------

subroutine igroup_attach_atom(g,a,l_)
! Attach atom `a` to igroup `g`.
! Optional: return l where a%gr(l) is the group ID and a%id(l) is the atom ID
! in g.
class(igroup),target         :: g
class(atom),target           :: a
type(atom_ap),allocatable    :: t_a(:)
integer                      :: n, m
integer,intent(out),optional :: l_
integer                      :: l

! Attempt to attach
call g%group_attach_atom(a,l)
if(present(l_)) l_=l
if(l==0) return

! Reallocate if needed
n=size(g%a)
m=g%nat+g%nlimbo
if(n<m) then
  allocate(t_a(m+g%pad))
  t_a(1:n) = g%a(1:n)
  call move_alloc(to=g%a,from=t_a)
endif
   
if(g%amax>=g%nat+g%nlimbo) then 
  ! Find index of previous deattached atom
  do n=1,g%amax
    if(.not.associated(g%a(n)%o)) exit
  enddo
  call werr('Index inconsistency',n>g%amax)
else
  ! Append
  g%amax=g%amax+1
  n=g%amax
endif
 
! Set atom index
g%a(n)%o=>a
! a%id(a%gri(g))=n
a%id(l)=n

end subroutine igroup_attach_atom
           
! Remove atoms
! ------------

subroutine igroup_detach_atom(g,a)
! Detach atom from `alist` and `a`
class(igroup)              :: g
class(atom),target         :: a
integer                    :: i
 
! Search index of `a`
i=a%gid(g)
if(i==-1) return  

! Nullify index
g%a(i)%o=>null()
              
! Detach atom
call group_detach_atom(g,a)

! FIXME: Seg fault with this  
! ! Update index if null count is above a fraction of the array size.  
! if ((g%amax-g%nat)>size(g%a)*g%aupd) then
!   call igroup_update_index(g)
! else
!   g%update=.false.
! endif

end subroutine igroup_detach_atom
     
subroutine igroup_update_index(g)
! Update index
use gems_strings, only: str
use gems_errors, only: wlog, werr
class(igroup)              :: g
class(atom),pointer        :: a
integer                    :: i,j,k

! Skip if there is not null components
if(g%amax==g%nat) return

i=0
do j=1,g%amax
  a=>g%a(j)%o
 
  ! Skip detached atoms
  if(.not.associated(a)) cycle
   
  ! Skip atoms in limbo
  if(g%b_limbo) then
    if(associated(a,target=g%limbo)) cycle
  endif
                          
  ! Skip if update is not needed
  i=i+1
  if(i==j) cycle

  ! Update atom index 
  g%a(i)%o=>a

  ! Update atom id
  k=a%gri(g)
  a%id(k)=i

enddo

call werr('Index inconsistency while updating',g%nat/=i)
g%amax=g%nat

! Signal out to rebuild internal arrays in extended types.
g%update=.true. 

! Write into log file 
! TODO: Build a warning if several calls are made to this subroutine
! and set up a CLI to control/deactivate `aupd`
call wlog('GRP','Index of group '//str(g%id)//' updated')

end subroutine igroup_update_index

! Ghosts
! ======

function new_ghost(o2,r) result(o)
use gems_program_types, only: box
class(atom),target,intent(in)  :: o2
type(atom),pointer             :: o
class(group),pointer           :: g
integer                        :: j
integer                        :: r(:)

! Allocate new ghost
allocate(o)
call o%init()
call ghost%attach(o)

! Assign image properties
o%prime=>o2
call atom_asign(o,o2)
o%q=o2%q
o%pbc(:)=.false.
o%pos(:)=o2%pos(:)+r(:)*box(:)
o%boxcr(:)=r(:)

! Add ghost to the image ngroups
do j=1,o2%ngr
  g => o2%gro(j)
  if(g%ghost) call g%attach(o)
enddo

end function

function isghost(r,rc)
! Return .true. if r inside the box+rcut. So actually it says if some position
! belong to a ghost region but also if is a local region
use gems_program_types, only: box
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

subroutine ghost_from_atom(a,rcut)
! Create ghost images for an atom.
use gems_program_types, only: box, n1cells
real(dp)                     :: r(dm)
real(dp),intent(in)          :: rcut
type(atom)                   :: a
integer                      :: m
 
do m =1,26
  r(:)=a%pos(:)+n1cells(m,:)*box(:)
  if (isghost(r(:),rcut)) a%ghost(m)%o => new_ghost(a,n1cells(m,:))
enddo
 
end subroutine ghost_from_atom
     
subroutine pbcghost_move()
! Move ghost atoms to reflect the motion of their local images
use gems_program_types, only: box
type(atom_dclist), pointer :: la
integer                    :: i

la => ghost%alist
do i = 1,ghost%nat
  la => la%next
  la%o%pos(:)=la%o%prime%pos(:)+box(:)*la%o%boxcr(:)
enddo

end subroutine

subroutine pbchalfghost(rcut)
! Create and destroy ghost for a pbc box.
use gems_program_types, only: box, n1cells
type(atom),pointer           :: o,og
integer                      :: i,j,k
type(atom_dclist),pointer    :: la
real(dp),intent(in)          :: rcut
integer                      :: reps(dm), pbc

! Make local atoms pbc
call do_pbc(sys)

la => sys%alist
do i = 1,sys%nat
  la => la%next
  o => la%o

  ! Check pbc
  pbc=count(o%pbc)

  ! TODO: when rcut is grater than half box add more images
  ! TODO: when rcut is less than half box, may have no sense to add/destroy
  ! ghost atoms

  ! Set ghost(1:3) at the box sides
  if (pbc==0) cycle
  do k=1,dm
    if (o%pos(k)>.5*box(k)) then
      reps(k)=-1
    else
      reps(k)=1
    endif
    og=>o%ghost(k)%o
    og%boxcr(:)=0
    og%boxcr(k)=reps(k)
    og%pos(:)=o%pos(:)+box(:)*og%boxcr(:)
    call ghost_switch(og,rcut)
  enddo 
      
  ! Set ghosts(4:6) at the box edges
  if (pbc==1) cycle
  do j=1,dm-1
    if (.not.o%pbc(j)) cycle
    do k=j+1,dm
      if(.not.o%pbc(k)) cycle
      og=>o%ghost(j+k+1)%o
      og%boxcr(:)=0
      og%boxcr(j)=reps(j)
      og%boxcr(k)=reps(k)
      og%pos(:)=o%pos(:)+box(:)*og%boxcr(:)
      call ghost_switch(og,rcut)
    enddo 
  enddo 
  
  ! Set ghosts(7) at the box korner
  if (pbc==2) cycle
  og=>o%ghost(7)%o
  og%boxcr(:)=reps(:)
  og%pos(:)=o%pos(:)+box(:)*og%boxcr(:)
  call ghost_switch(og,rcut)

enddo


end subroutine
     

subroutine ghost_switch(og,rcut)
type(atom),pointer      :: og, o
class(group),pointer    :: g
integer                 :: i
real(dp)                :: rcut
 
o=>og%prime

! Check if ghost atoms is needed
if (isghost(og%pos(:),rcut)) then

  ! If ghost already activated, exit
  if(og%ngr>0) return

  ! Ensure ghosts are in the same groups than prime atom
  ! FIXME: detach if ghost is in another group
  do i=1,o%ngr
    g => o%gro(i)
    if(g%ghost) call g%attach(og)
  enddo
  call ghost%attach(og)
            
else

  ! If ghost already deactivated, exit
  if(og%ngr==0) return

  do while (og%ngr/=0)
    g => og%gro(og%ngr)
    call g%detach(og)
  enddo   

end if

end subroutine    
                        

subroutine pbcghost(rcut)
! Create and destroy ghost for a pbc box.
use gems_program_types, only: box, n1cells
real(dp)                     :: r(dm)
type(atom),pointer           :: o,og
integer                      :: i,m
type(atom_dclist),pointer    :: la
real(dp),intent(in)          :: rcut

! Make local atoms pbc
call do_pbc(sys)

la => sys%alist
do i = 1,sys%nat
  la => la%next
  o => la%o

  ! Loop trough box images
  ! Note: seems to check proximity to the 6 faces can be fast?
  do m = 1,26
    og => la%o%ghost(m)%o

    ! Check if ghost atoms is needed
    r(:)=o%pos(:)+n1cells(m,:)*box(:)
    if (isghost(r(:),rcut)) then

      ! If ghost already exist, skip
      if(associated(og)) cycle
    
      ! Add new ghost
      og => new_ghost(o,n1cells(m,:))

    else
 
      ! If ghost exist, remove it
      if(associated(og)) call og%dest()
      og => null()
         
    end if

    o%ghost(m)%o => og

  enddo

enddo


end subroutine

! PBC
! ---

subroutine set_pbc(g,pbc)
class(group),intent(inout)  :: g
class(atom_dclist),pointer  :: la
logical,intent(in)         :: pbc(dm)
integer                    :: i

la => g%alist
do i = 1,g%nat
  la => la%next
  call atom_setpbc(la%o,pbc)
enddo

end subroutine set_pbc
 
subroutine do_pbc(g)
use gems_program_types, only: box
class(group),intent(inout)  :: g
class(atom_dclist),pointer  :: la
type(atom),pointer          :: o
integer                     :: i,k

la => g%alist
do i = 1,g%nat
  la => la%next
  o => la%o

  do k=1,dm
    if (o%pbc(k)) then
      if(o%pos(k)>=box(k)) then
        o%pos(k)=o%pos(k)-box(k)
        o%pos_old(k)=o%pos_old(k)-box(k)
        o%boxcr(k)=o%boxcr(k)+1
      elseif(o%pos(k)<0.0_dp) then
        o%pos(k)=o%pos(k)+box(k)
        o%pos_old(k)=o%pos_old(k)+box(k)
        o%boxcr(k)=o%boxcr(k)-1
      endif
    endif
  enddo
enddo
     
end subroutine do_pbc
                     

! Hyper vector Constructors
! =========================

! El echo de que el grupo sea una lista linkeada esta indicando que es una
! selecciona caprichosa de atomos. Esto permite aplicar disitntas subrutinas a
! distintas selecciones caprichosas.

! En muchas subrutinas, y para muchas cosas, es util constar con vectores que
! representen el estado de este grupo.

! Las subrutinas que siguen a continuacion son un parche para lograr esto. Las
! variables atomicas, se encuentran inicialmente apuntando a vectores en el
! systema. Luego de la seleccion caprichos por un grupo, se puede copiar
! las variables atomicas a un vector y asociar los atomos a ese vector,
! logrando asi un ordenamiento del grupo de forma vectorial. Si este
! ordenamiento no existe, los punteros del grupo tendran la condicion null. Si
! por el contrario este ordenamiento existe, los punteros del grupo apuntaran
! al vector determinado. Esto involucra ciertos manejos que se deben realizar
! con las siguientes subrutinas:

subroutine group_switch_vectorial(g,switched)
class(group),intent(inout)     :: g
logical,optional,intent(out)   :: switched
integer                        :: i
type (atom_dclist),pointer     :: la

if(present(switched)) switched=.false.
if(associated(g%pp)) return
if(present(switched)) switched=.true.

allocate(g%pp(g%nat*dm))
allocate(g%pv(g%nat*dm))
allocate(g%pa(g%nat*dm))
allocate(g%pf(g%nat*dm))

la => g%alist
do i = 1,g%nat
  la => la%next
  g%pp((i-1)*dm+1:i*dm) =  la%o%pos
  g%pv((i-1)*dm+1:i*dm) =  la%o%vel
  g%pf((i-1)*dm+1:i*dm) =  la%o%force
  g%pa((i-1)*dm+1:i*dm) =  la%o%acel
  deallocate(la%o%pos  )
  deallocate(la%o%vel  )
  deallocate(la%o%force)
  deallocate(la%o%acel )
  la%o%pos    => g%pp((i-1)*dm+1:i*dm)
  la%o%vel    => g%pv((i-1)*dm+1:i*dm)
  la%o%force  => g%pf((i-1)*dm+1:i*dm)
  la%o%acel   => g%pa((i-1)*dm+1:i*dm)
enddo

end subroutine

subroutine group_switch_objeto(g,switched)
class(group),intent(inout)     :: g
logical,optional,intent(out)   :: switched
integer                        :: i
type (atom_dclist),pointer     :: la

if(present(switched)) switched=.false.
if(.not.associated(g%pp)) return
if(present(switched)) switched=.true.


la => g%alist
do i = 1,g%nat
  la => la%next
  allocate(la%o%pos  (dm))
  allocate(la%o%vel  (dm))
  allocate(la%o%force(dm))
  allocate(la%o%acel (dm))
  la%o%pos    = g%pp((i-1)*dm+1:i*dm)
  la%o%vel    = g%pv((i-1)*dm+1:i*dm)
  la%o%force  = g%pf((i-1)*dm+1:i*dm)
  la%o%acel   = g%pa((i-1)*dm+1:i*dm)
enddo

deallocate(g%pp); g%pp=>null()
deallocate(g%pv); g%pv=>null()
deallocate(g%pa); g%pa=>null()
deallocate(g%pf); g%pf=>null()

end subroutine

subroutine binleft(a,val,m,found)
! Find the index m where a(m)<=val<a(m+1) using a binary search.
! The function return true if a(m)=val and false otherwise.
! a(:) must be sorted.
integer,intent(in)  :: a(:),val
integer             :: r,l,j
integer,intent(out) :: m
logical,intent(out) :: found

! Left and right bounds
l=lbound(a,1)-1
r=ubound(a,1)+1

! Shrink bounds
do while (r>l+1)
  m=(r+l)/2
  j=a(m)
  
  if(val==j) then
    found=.true.
    return
  endif

  if(val>j) then
    l=m
  else
    r=m
  endif
enddo  

! Worst case
if(r<=size(a)) then
  if(a(r)==val) then
    m=r
    found=.true.
    return
  endif
endif

! Not found
found=.false.
m=l 

end subroutine binleft
  
end module gems_groups

