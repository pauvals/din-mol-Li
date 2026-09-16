program dana
  use gems_groups, only: group, atom, atom_dclist, sys, gindex
  use gems_neighbor, only: gems_sync_execution_views => sync_execution_views, ngroup, nn_vlist, nupd_vlist, nb_dcut
  use gems_elements, only: set_z, elements_init,elements
  use gems_errors, only: timer_dump, timer_start, sclock_t1, sclock_t2, wstd, sclock_max, sclock_rate,logunit
  use gems_constants, only: time1, dp, sp

  implicit none
 
  ! Cosas del sistema  #what a quilombo
  integer             :: n                    ! Num de particulas
  real(dp), parameter :: o=0._dp               ! Origen de la caja
  real(dp)            :: z0, zmax              ! Limites verticales
  real(dp)            :: xi, yi               ! Tamaño en x e y de caja sim.
  real(dp)            :: domain_lo(3)[*], domain_hi(3)[*]
  real(dp)            :: domain_global_lo(3), domain_global_hi(3)
  integer             :: domain_dims(3)=[1,1,1], domain_coord(3)=[0,0,0]
  logical             :: domain_periodic(3)=[.true.,.true.,.false.]
  character(8)        :: domain_topology='z'
  real(dp)            :: dif, dif_sei, dif_sc ! Difusion
 
  ! Variables dinámicas, agrupadas en átomo
  type(atom), pointer             :: pa=>null()
  real(dp)                        :: prob,max_vel=0._dp, msd_u= 0._dp, msd_t= 0._dp, msd_max= 0._dp

  ! Contador de choques (superposiciones)
  integer             :: choques=0, choques2=0, choques3=0

  ! Parametros de integración
  real(dp)            :: h          ! paso de tiempo
  integer             :: nst        ! nro de pasos
  integer             :: nwr        ! paso de escritura
  real(dp)            :: t=0.0_dp   ! tiempo en ps

  logical             :: new_sim, electrodo

  ! Varios
  integer             :: idum         ! Semilla
  integer             :: i,j,k        ! Enteros
  integer             :: try, depo    ! Parts. cercanas a electrodo y deposic.
  integer             :: prof_integrate_tick, prof_sync_pos_tick, prof_sync_views_tick
  integer             :: prof_overlap_tick, prof_output_tick
  real(dp)            :: prof_integrate=0._dp, prof_sync_pos=0._dp, prof_sync_views=0._dp
  real(dp)            :: prof_overlap=0._dp, prof_output=0._dp
  integer             :: prof_n_integrate=0, prof_n_sync_pos=0, prof_n_sync_views=0
  integer             :: prof_n_overlap=0, prof_n_output=0
  integer,allocatable :: overlap_gid_lo(:), overlap_gid_hi(:), overlap_pair_owner(:)
  integer,allocatable :: overlap_pairs_by_owner(:)
  logical,allocatable :: overlap_cross_domain(:), overlap_deposition_pair(:)
  integer             :: n_overlap_pairs=0, overlap_conflict_scans=0
  integer             :: overlap_pair_observations=0, overlap_local_observations=0
  integer             :: overlap_cross_observations=0, overlap_deposition_observations=0
  integer             :: overlap_duplicate_observations=0

  real(dp),allocatable:: co_pos(:,:)[:]
  real(dp),allocatable:: co_old_cg(:,:)[:], co_vel(:,:)[:]
  real(dp),allocatable:: recv_pos(:,:), recv_old_cg(:,:), recv_vel(:,:)
  integer,allocatable :: atom_owner(:)
  integer,allocatable :: co_uid(:)[:], co_element(:)[:], co_new_owner(:)[:]
  integer,allocatable :: recv_uid(:), recv_element(:), recv_owner(:), uid_to_sys(:)
  logical,allocatable :: atom_visible(:), co_skip(:)[:]
  logical,allocatable :: recv_skip(:)
  integer,allocatable :: co_gcmc_kind(:)[:], co_gcmc_uid(:)[:]
  real(dp),allocatable:: co_gcmc_pos(:,:)[:]
  logical             :: focus_dirty = .false.
  integer             :: next_atom_uid=0
  integer             :: co_owned_count[*]=0
  integer             :: co_gcmc_count[*]=0, co_gcmc_creations[*]=0
  real(dp)            :: co_gcmc_active_volume[*]=0._dp
  real(dp)            :: gcmc_patch_volume=0._dp, gcmc_accessible_volume=0._dp
  integer             :: gcmc_local_attempts=0
  integer             :: overlap_round_actions[*]=0, overlap_round_unresolved[*]=0
  integer             :: output_try[*]=0, output_depo[*]=0
  integer             :: overlap_parallel_rounds=0, overlap_parallel_residuals=0
   
  ! GCMC
  type(group)         :: gcmc
  real(dp)            :: act
  integer             :: nadj
  type(atom_dclist), pointer :: la, lb 
  
  ! Group used to compute neighbors list
  type(ngroup)        :: hs
  type(group), target :: hs_focus, hs_candidate
  type(atom),pointer  :: o1
   
  ! Init wall time
  call system_clock(count_rate=sclock_rate)
  call system_clock(count_max=sclock_max)
  call system_clock(sclock_t1)

  sync all

  if (this_image() == 1) then
    open(25,File='version')
    write(25,*) PACKAGE_VERSION
    close(25)
  endif

  ! Chemical elements used
  call elements_init()
  call set_z(1,sym='Li',mass=6.94_dp)
  call set_z(2,sym='CG',mass=6.94_dp)
  call set_z(3,sym='F',mass=6.94_dp)
            
  ! Para contar partículas que intentan deposic.
  try= 0
  depo= 0
                      
  ! Init general del sistema de grupos
  call gindex%init()
  call sys%init()
 
  ! Leer run parameters:  semilla, probabilidad, nro. pasos
  ! y otros valores iniciales
  call entrada()
  if(num_images()>1) idum=idum-104729*(this_image()-1)
        
  ! Init lista para choques de esferas duras
  call hs_focus%init()
  call hs_candidate%init()
  call hs%init()
  call hs%bind_groups(hs_focus, hs_candidate)
  ! Halo churn reuses detached slots, so avoid compacting the candidate index
  ! while its membership is being updated incrementally.
  hs%candidate%aupd=1._dp
  call hs%setrc(3.2_dp) ! Maximo radio de corte

  ! Configure the GCMC reservoir.
  call config_run()

  ! Generar posiciones inic.
  if (new_sim .and. this_image() == 1) call gen_pos_inic()
  sync all

  ! Define tamaño en z del sistema (sin contar el reservorio), y lee posic. inic.
  ! de partículas 
  ! Read initial positions (optionally from a previous run).
  call config_inic()
  call sync_owned_positions()

  call gcmc%init()
  la=>sys%alist
  do j=1,sys%nat
    la=>la%next
    if(la%o%z==1) call gcmc%attach(la%o)
  enddo

  !Search neighbors
  call sync_execution_views()

  call sync_owned_positions()


  ! Output 
  if (this_image() == 1) then
    open(11,File='Li.xyz')
    open(15,File='try.dat')
    open(16,File='depo.dat')
  endif

  call salida() ! config. inic. en el primer paso ;)

  call timer_start()

  do i=1,nst
    ! Da un paso #wii

    call prof_start(prof_integrate_tick)
    call cbrownian_hs(hs,h)
    call prof_stop(prof_integrate_tick, prof_integrate, prof_n_integrate)
    call sync_owned_positions()

    ! Update neighbors for serial overlap solving. Multi-image overlap rebuilds
    ! the all-atom view on image 1 inside solve_overlaps_domain().
    if(num_images() == 1) call sync_execution_views()
 
    ! Solve overlaps
    call solve_overlaps_domain()
       
    ! Update neighbors (give error if not used again)
    call sync_execution_views()
         
    msd_t= msd_t/hs%focus%nat
    msd_max= max(msd_max,msd_t)
 
    ! Add new CG to the ref group
    ! NOTE: Should we allow chain reaction instead?
    la=> hs%focus%alist
    do j = 1, hs%focus%nat
      la=> la%next
      o1=> la%o
      if (o1%sym/='F') cycle
      call o1%setz(2)
      lb => la%prev
      call hs_focus%detach(o1)
      call hs%focus%detach(o1)
      la => lb
      call gcmc%detach(o1)
    enddo
    
    call gcmc_run(gcmc)
    call sync_owned_positions()

    ! Output 
    if (mod(i,nwr)==0) then
      call prof_start(prof_output_tick)
      call salida()
      call prof_stop(prof_output_tick, prof_output, prof_n_output)
      ! print*, 't', t
      ! print*, 'msd_t', msd_t
      ! print*, 'msd_max', msd_max
    endif
   
    ! Timing information
    call timer_dump(i,nup=nupd_vlist)
  
    t=t+h
   
  enddo
 
  if (this_image() == 1) then
    close(11)
    close(15)
    close(16)
  endif

  
  ! FIN
  call cpu_time(time1)
  if (time1<60) then
    call wstd(); write(logunit,'("cpu time: ",f10.3," s")') time1
  elseif(time1<3600) then
    call wstd(); write(logunit,'("cpu time: ",i0," m ",f10.3," s")') int(time1/60.0_dp),mod(time1,60.d0)
  else
    call wstd(); write(logunit,'("cpu time: ",i0," h ",i0," m")') int(time1/3600.0_dp),int(mod(time1,3600.d0)/60.0_dp)
  endif 

  call system_clock(sclock_t2)
  time1=(sclock_t2-sclock_t1)/real(sclock_rate)
  if (time1<60) then
    call wstd(); write(logunit,'("wall time: ",f10.3," s")') time1
  elseif(time1<3600) then
    call wstd(); write(logunit,'("wall time: ",i0," m ",f10.3," s")') int(time1/60.0_dp),mod(time1,60.d0)
  else
    call wstd(); write(logunit,'("wall time: ",i0," h ",i0," m")') int(time1/3600.0_dp),int(mod(time1,3600.d0)/60.0_dp)
  endif 
 
!  call wstd(); write(logunit,*) 'with ',rupdate,' neighbour list updates'
  call wstd(); write(logunit,'("coarray images: ",i0,"; active domain image: ",i0)') num_images(), this_image()
  call wstd(); write(logunit,'("domain topology: ",a,"; grid: ",i0," x ",i0," x ",i0, &
    & "; coordinates: ",i0,",",i0,",",i0)') trim(domain_topology), domain_dims, domain_coord
  call wstd(); write(logunit,'("domain x range: [",f10.3,",",f10.3,")")') domain_lo(1), domain_hi(1)
  call wstd(); write(logunit,'("domain y range: [",f10.3,",",f10.3,")")') domain_lo(2), domain_hi(2)
  call wstd(); write(logunit,'("domain z range: [",f10.3,",",f10.3,")")') domain_lo(3), domain_hi(3)
  call wstd(); write(logunit,'("domain atoms: focus=",i0," candidate=",i0)') hs_focus%nat, hs_candidate%nat
  if(allocated(atom_visible)) then
    call wstd(); write(logunit,'("domain sync visible atoms: ",i0)') count(atom_visible)
  endif
  call wstd(); write(logunit,'("vecinos actualizados: ",i0," veces")') nupd_vlist
  call wstd(); write(logunit,'("maximo numero de vecinos en algun paso: ",i0)') nn_vlist
  call wstd(); write(logunit,'("maximo desplazamiento en algun paso: ",e10.3)') sqrt(max_vel)*h
  call wstd(); write(logunit,'("numero total de choques: ",i0)') choques
  call wstd(); write(logunit,'("numero total de choques en 2da vuelta: ",i0)') choques2
  call wstd(); write(logunit,'("numero total de choques sin solucion: ",i0)') choques3
  call wstd(); write(logunit,'("MSD maximo en x-y: ",e10.3)') msd_max
  call wstd(); write(logunit,'("GCMC patch/access volume: ",f12.3," / ",f12.3, &
    & "; local attempts: ",i0)') gcmc_patch_volume, gcmc_accessible_volume, gcmc_local_attempts
  call report_profile()


  call elements%destroy()

  sync all



contains

subroutine prof_start(tick)
integer,intent(out) :: tick

call system_clock(tick)

end subroutine prof_start


subroutine prof_stop(tick_start, total, ncall)
integer,intent(in)    :: tick_start
real(dp),intent(inout) :: total
integer,intent(inout)  :: ncall
integer               :: tick

call system_clock(tick)
total = total + (tick-tick_start)/real(sclock_rate,dp)
ncall = ncall + 1

end subroutine prof_stop


subroutine report_profile()

  call wstd(); write(logunit,'("profile integrate: ",f10.3," s calls=",i0)') prof_integrate, prof_n_integrate
  call wstd(); write(logunit,'("profile sync_positions: ",f10.3," s calls=",i0)') prof_sync_pos, prof_n_sync_pos
  call wstd(); write(logunit,'("profile sync_neighbors: ",f10.3," s calls=",i0)') prof_sync_views, prof_n_sync_views
  call wstd(); write(logunit,'("profile overlap: ",f10.3," s calls=",i0)') prof_overlap, prof_n_overlap
  call wstd(); write(logunit,'("profile output: ",f10.3," s calls=",i0)') prof_output, prof_n_output
if(num_images()>1 .and. this_image()==1) then
  call wstd(); write(logunit,'("parallel overlap rounds: ",i0,"; residual observations: ",i0)') &
    overlap_parallel_rounds, overlap_parallel_residuals
  if(overlap_conflict_scans>0) then
    call wstd(); write(logunit,'("overlap conflict list: current=",i0," scans=",i0)') &
      n_overlap_pairs, overlap_conflict_scans
    call wstd(); write(logunit,'("overlap conflict observations: total=",i0," local=",i0," cross=",i0, &
      & " deposition=",i0," duplicates=",i0)') overlap_pair_observations, overlap_local_observations, &
      overlap_cross_observations, overlap_deposition_observations, overlap_duplicate_observations
    if(allocated(overlap_pairs_by_owner)) then
      call wstd(); write(logunit,'("overlap pair-owner observations: min=",i0," max=",i0)') &
        minval(overlap_pairs_by_owner), maxval(overlap_pairs_by_owner)
    endif
  endif
endif

end subroutine report_profile

subroutine entrada()
! Some run conditions from entrada.ini
use gems_neighbor, only: nb_dcut  
open(15,File='entrada.ini')
read(15,*) idum
read(15,*) prob
read(15,*) h 
read(15,*) nst
read(15,*) nwr
read(15,*) xi 
read(15,*) yi
read(15,*) z0 
read(15,*) zmax 
read(15,*) dif_sc 
read(15,*) dif_sei 
read(15,*) nb_dcut
read(15,*) new_sim
read(15,*) electrodo
close(15)
end subroutine entrada


subroutine gen_pos_inic()
! Create initial positions. For 1 M [Li+]
use gems_program_types, only: distance 
integer             :: n 
real(dp),allocatable:: r(:,:) !posic.
real(dp)            :: Mol, v1(3),v2(3),dif_pos(3), alto, dif2 ! molaridad, diferencia entre vectores posic., alto caja sim.
real(dp),parameter  :: r0=3.2_dp, mLi= 6.94_dp
integer             :: i,j,k,n_CG
logical,parameter   :: pbc(3)=[.true.,.true.,.false.]
character(10)       :: sym

! Calculo nro. de partículas
! n= Molaridad * Volumen(A^3) * 1e-27 L/A^3 * 6.022e23 (Nro Avogadro)
! Molaridad usada = 1 M
Mol = 1._dp 
alto = zmax - o
n = Mol * xi * yi * alto * 6.022e-4 

open(12,File='pos_inic.xyz')

! Leer primero las posic. inic. de CG.
if (electrodo) then
  open(13,File='CG_inic.xyz')
  read(13,*) n_CG
  read(13,*)

  ! nro. total de átomos
  allocate(r(n+n_CG,3))
  write(12,*) n+n_CG  
  write(12,*)

  do i=1, n_CG
    read(13,*) sym,r(i,:)
    ! Escribo CGs existentes. Explicito que sym va a ser CG.
    write(12,'(a,4(x,f25.12))') 'CG',r(i,1),r(i,2),r(i,3),mLi
  enddo

else
  n_CG=0
  ! nro. total de átomos
  allocate(r(n,3))
  write(12,*) n  
  write(12,*)
endif

! Recorro todas las particulas a crear
do i=n_CG+1, n+n_CG

  ! Recorro los intentos 
  intento: do k=1,10000
    !Posic. de Li
    r(i,1)=ran(idum)*xi
    r(i,2)=ran(idum)*yi
    r(i,3)=(ran(idum)*alto) + o ! en z 

    ! Recorro las particulas ya creadas 
    do j=1,i-1 

      !Veo distancia Li-Li
      v1(:)=r(i,:)
      v2(:)=r(j,:)
      dif_pos(:)=distance(v2,v1,pbc)
      dif2=dot_product(dif_pos,dif_pos)

      if (dif2<r0*r0) cycle intento

    enddo

    exit
  enddo intento

  if(k==10001) then
      print *, 'Maximo numero de intentos alcanzado'
      stop
  endif

write(12,'(a,4(x,f25.12))') 'Li',r(i,1),r(i,2),r(i,3),mLi

enddo
close(12)
if(electrodo) close(13)


end subroutine gen_pos_inic


subroutine config_run()
use gems_program_types, only: tbox, box_setvars  

! Set the box
tbox(:,:)= 0._dp
tbox(1,1)= xi
tbox(2,2)= yi
tbox(3,3)= zmax
call box_setvars()

open(15,File='movedor.ini')
read(15,*) act, nadj
close(15)

allocate(co_gcmc_kind(max(1,nadj))[*])
allocate(co_gcmc_uid(max(1,nadj))[*])
allocate(co_gcmc_pos(3,max(1,nadj))[*])
co_gcmc_kind(:)=0
co_gcmc_uid(:)=0
co_gcmc_pos(:,:)=0._dp
end subroutine config_run


subroutine config_inic()
! configuración inicial
use gems_program_types, only: tbox, box_setvars
integer :: i
!type(atom_dclist), pointer :: la
character(10)     :: sym

if(num_images() > 1) then
  tbox(3,3)=zmax
  call box_setvars()
endif
call set_domain_bounds()

open(11,File='pos_inic.xyz')
read(11,*) n
read(11,*)

if(allocated(atom_owner)) deallocate(atom_owner)
if(allocated(co_element)) deallocate(co_element)
if(allocated(co_uid)) deallocate(co_uid)
if(allocated(co_new_owner)) deallocate(co_new_owner)
if(allocated(atom_visible)) deallocate(atom_visible)
if(allocated(co_skip)) deallocate(co_skip)
if(allocated(co_pos)) deallocate(co_pos)
if(allocated(co_old_cg)) deallocate(co_old_cg)
if(allocated(co_vel)) deallocate(co_vel)
if(allocated(recv_pos)) deallocate(recv_pos)
if(allocated(recv_old_cg)) deallocate(recv_old_cg)
if(allocated(recv_vel)) deallocate(recv_vel)
if(allocated(recv_element)) deallocate(recv_element)
if(allocated(recv_uid)) deallocate(recv_uid)
if(allocated(recv_owner)) deallocate(recv_owner)
if(allocated(uid_to_sys)) deallocate(uid_to_sys)
if(allocated(recv_skip)) deallocate(recv_skip)
allocate(atom_owner(n))
allocate(co_element(n)[*])
allocate(co_uid(n)[*])
allocate(co_new_owner(n)[*])
allocate(atom_visible(n))
allocate(co_skip(n)[*])
allocate(co_pos(3,n)[*])
allocate(co_old_cg(3,n)[*],co_vel(3,n)[*])
allocate(recv_pos(3,n),recv_old_cg(3,n),recv_vel(3,n))
allocate(recv_uid(n),recv_element(n),recv_owner(n),recv_skip(n))
  atom_owner(:)=1
  co_element(:)=119
  co_uid(:)=0
  co_new_owner(:)=1
  atom_visible(:)=.true.
  co_skip(:)=.false.
  co_pos(:,:)=0._dp
  co_old_cg(:,:)=0._dp
  co_vel(:,:)=0._dp
  focus_dirty=.true.

next_atom_uid=0
do i=1,n

  ! Allocatea el atomo-como un puntero del tipo 'atom'
  allocate(pa)
  call pa%init()
  next_atom_uid=next_atom_uid+1
  pa%uid=next_atom_uid

  read(11,*) sym,pa%pos(:)
  call pa%setsym(sym) ! Asigna algunos valores según tipo de átomo

  ! Agrega atom al sistema (sys)
  call sys%attach(pa)

  ! Agrego atomos a los grupos
  ! NOTA: hs%attach tiene que estar al final
  ! para que ande la lista de vecinos
  atom_owner(pa%gid(sys)) = owner_image_for_position(pa%pos)

  if (candidate_atom(pa)) call hs_candidate%attach(pa)
  if (pa%sym/='CG') then
    if (owned_atom(pa)) call hs_focus%attach(pa)
  endif
  call hs%attach(pa)

  ! Set pbc
  pa%pbc(:)=.true.
  pa%pbc(3)=.false.

  ! Libero puntero para siguiente allocate
  pa=>null()

end do
close(11)

end subroutine config_inic


subroutine set_domain_bounds()
! Cartesian domain metadata. Keep z slabs as the validation topology while
! ownership and halo predicates are generalized for later x-y columns.
real(dp) :: width
integer  :: axis

call select_domain_grid()
call image_coordinates(this_image(),domain_coord)
domain_global_lo=[o,o,o]
domain_global_hi=[xi,yi,zmax]

do axis=1,3
  width=(domain_global_hi(axis)-domain_global_lo(axis))/real(domain_dims(axis),kind=dp)
  domain_lo(axis)=domain_global_lo(axis)+real(domain_coord(axis),kind=dp)*width
  if(domain_coord(axis)==domain_dims(axis)-1) then
    domain_hi(axis)=domain_global_hi(axis)
  else
    domain_hi(axis)=domain_lo(axis)+width
  endif
enddo
sync all

end subroutine set_domain_bounds


subroutine select_domain_grid()
use gems_errors, only: werr
integer  :: candidate_x, candidate_y, best_x, status
real(dp) :: aspect_error, best_error

call get_environment_variable('DANA_DOMAIN_TOPOLOGY',domain_topology,status=status)
if(status/=0 .or. len_trim(domain_topology)==0) domain_topology='z'

select case(trim(domain_topology))
case('z')
  domain_dims=[1,1,num_images()]
case('xy')
  best_x=1
  best_error=huge(best_error)
  do candidate_x=1,num_images()
    if(mod(num_images(),candidate_x)/=0) cycle
    candidate_y=num_images()/candidate_x
    aspect_error=abs(log((xi/real(candidate_x,dp))/(yi/real(candidate_y,dp))))
    if(aspect_error<best_error) then
      best_error=aspect_error
      best_x=candidate_x
    endif
  enddo
  domain_dims=[best_x,num_images()/best_x,1]
case default
  call werr('DANA_DOMAIN_TOPOLOGY must be z or xy',.true.)
end select

end subroutine select_domain_grid


subroutine image_coordinates(image,coord)
integer,intent(in)  :: image
integer,intent(out) :: coord(3)
integer             :: rank

rank=image-1
coord(1)=mod(rank,domain_dims(1))
rank=rank/domain_dims(1)
coord(2)=mod(rank,domain_dims(2))
coord(3)=rank/domain_dims(2)

end subroutine image_coordinates


integer function image_from_coordinates(coord) result(image)
integer,intent(in) :: coord(3)

image=1+coord(1)+domain_dims(1)*(coord(2)+domain_dims(2)*coord(3))

end function image_from_coordinates


logical function owned_atom(a) result(owned)
class(atom),intent(in) :: a

owned=owner_image_for_position(a%pos)==this_image()

end function owned_atom


logical function candidate_atom(a) result(candidate)
class(atom),intent(in) :: a

candidate=candidate_position(a%pos)

end function candidate_atom


logical function candidate_position(r) result(candidate)
real(dp),intent(in) :: r(3)
real(dp)            :: halo
integer             :: axis

if(num_images() == 1) then
  candidate = .true.
  return
endif

halo = hs%rcut + nb_dcut
candidate=owner_image_for_position(r)==this_image()
if(candidate) return

candidate=.true.
do axis=1,3
  if(.not.coordinate_in_domain_halo(r(axis),axis,halo)) then
    candidate=.false.
    return
  endif
enddo

end function candidate_position


logical function coordinate_in_domain_halo(value,axis,halo) result(visible)
real(dp),intent(in) :: value, halo
integer,intent(in)  :: axis
real(dp)            :: coordinate, length

visible=.true.
if(domain_dims(axis)==1) return

coordinate=value
length=global_domain_length(axis)
if(domain_periodic(axis)) coordinate=modulo(coordinate-domain_global_lo(axis),length)+domain_global_lo(axis)

visible=coordinate>=domain_lo(axis)-halo .and. coordinate<domain_hi(axis)+halo
if(visible .or. .not.domain_periodic(axis)) return
visible=coordinate+length>=domain_lo(axis)-halo .and. coordinate+length<domain_hi(axis)+halo
if(visible) return
visible=coordinate-length>=domain_lo(axis)-halo .and. coordinate-length<domain_hi(axis)+halo

end function coordinate_in_domain_halo


integer function owner_image_for_position(r) result(owner)
real(dp),intent(in) :: r(3)
integer             :: axis, coord(3), c
real(dp)            :: coordinate, lower, length, width

do axis=1,3
  lower=domain_global_lo(axis)
  length=global_domain_length(axis)
  coordinate=r(axis)
  if(domain_periodic(axis)) coordinate=modulo(coordinate-lower,length)+lower

  coord(axis)=domain_dims(axis)-1
  width=length/real(domain_dims(axis),kind=dp)
  do c=0,domain_dims(axis)-1
    if(coordinate>=lower+real(c,kind=dp)*width .and. &
      coordinate<lower+real(c+1,kind=dp)*width) then
      coord(axis)=c
      exit
    endif
  enddo
enddo

owner=image_from_coordinates(coord)

end function owner_image_for_position


real(dp) function global_domain_length(axis) result(length)
integer,intent(in) :: axis

select case(axis)
case(1)
  length=domain_global_hi(1)-domain_global_lo(1)
case(2)
  length=domain_global_hi(2)-domain_global_lo(2)
case(3)
  length=domain_global_hi(3)-domain_global_lo(3)
case default
  length=0._dp
end select

end function global_domain_length


subroutine sync_owned_positions(global_snapshot)
! Pack each owner's records contiguously so remote state arrives in a few array
! transfers instead of many scalar coindexed reads.
use gems_errors, only: werr
logical,intent(in),optional :: global_snapshot
type(atom_dclist), pointer :: la
type(atom), pointer        :: received_atom
integer                   :: i, gid, owner, nrecv, slot, new_owner
logical                   :: gather_all

call prof_start(prof_sync_pos_tick)
call ensure_sync_storage()
call rebuild_uid_index()
gather_all=.false.
if(present(global_snapshot)) gather_all=global_snapshot
if(num_images()==1) then
  call refresh_domain_ownership()
  call prof_stop(prof_sync_pos_tick,prof_sync_pos,prof_n_sync_pos)
  return
endif

atom_visible(:)=.false.
la=>sys%alist
do i=1,sys%nat
  la=>la%next
  gid=la%o%gid(sys)
  atom_visible(gid)=atom_owner(gid)==this_image()
enddo
co_owned_count=0
la => sys%alist
do i=1,sys%nat
  la => la%next
  gid=la%o%gid(sys)
  if(atom_owner(gid) /= this_image()) cycle
  new_owner=owner_image_for_position(la%o%pos)
  if(.not.gather_all .and. new_owner==this_image()) then
    if(.not.position_near_owner_boundary(la%o%pos,this_image())) cycle
  endif
  co_owned_count=co_owned_count+1
  slot=co_owned_count
  co_uid(slot) = la%o%uid
  co_new_owner(slot) = new_owner
  co_pos(:,slot) = la%o%pos(:)
  co_old_cg(:,slot) = la%o%old_cg(:)
  co_vel(:,slot) = la%o%vel(:)
  co_element(slot) = la%o%z
  co_skip(slot) = la%o%skip
enddo

sync all

do owner=1,num_images()
  nrecv=co_owned_count[owner]
  if(nrecv==0) cycle
  recv_uid(:nrecv)=co_uid(:nrecv)[owner]
  recv_owner(:nrecv)=co_new_owner(:nrecv)[owner]
  recv_pos(:,:nrecv)=co_pos(:,:nrecv)[owner]
  recv_old_cg(:,:nrecv)=co_old_cg(:,:nrecv)[owner]
  recv_vel(:,:nrecv)=co_vel(:,:nrecv)[owner]
  recv_element(:nrecv)=co_element(:nrecv)[owner]
  recv_skip(:nrecv)=co_skip(:nrecv)[owner]

  do slot=1,nrecv
    call werr('Invalid atom UID in coarray snapshot',recv_uid(slot)<1 .or. &
      recv_uid(slot)>size(uid_to_sys))
    i=uid_to_sys(recv_uid(slot))
    call werr('Unknown atom UID in coarray snapshot',i<1)
    received_atom=>sys%a(i)%o
    if(atom_owner(i)/=owner) then
      call wstd(); write(logunit,'("GCMC/domain ownership mismatch: image=",i0," uid=",i0, &
        & " sender=",i0," local_owner=",i0," sys_gid=",i0)') &
        this_image(),recv_uid(slot),owner,atom_owner(i),i
      call werr('Coarray atom ownership mismatch',.true.)
    endif
    if(received_atom%z/=recv_element(slot)) call received_atom%setz(recv_element(slot))
    received_atom%pos(:)=recv_pos(:,slot)
    received_atom%old_cg(:)=recv_old_cg(:,slot)
    received_atom%vel(:)=recv_vel(:,slot)
    received_atom%skip=recv_skip(slot)
    atom_owner(i)=recv_owner(slot)
    if(recv_owner(slot)==this_image() .or. candidate_atom(received_atom)) atom_visible(i)=.true.
  enddo
enddo

sync all

call refresh_domain_ownership()
call prof_stop(prof_sync_pos_tick, prof_sync_pos, prof_n_sync_pos)

end subroutine sync_owned_positions


logical function position_near_owner_boundary(r,owner) result(near_boundary)
real(dp),intent(in) :: r(3)
integer,intent(in)  :: owner
real(dp)            :: coordinate, halo, length, lower, upper, width
integer             :: axis, coord(3)

near_boundary=.false.
halo=hs%rcut+nb_dcut
call image_coordinates(owner,coord)
do axis=1,3
  if(domain_dims(axis)==1) cycle
  length=global_domain_length(axis)
  coordinate=r(axis)
  if(domain_periodic(axis)) coordinate=modulo(coordinate-domain_global_lo(axis),length)+domain_global_lo(axis)
  width=length/real(domain_dims(axis),dp)
  lower=domain_global_lo(axis)+real(coord(axis),dp)*width
  upper=lower+width
  if((coordinate>=lower-halo .and. coordinate<lower+halo) .or. &
    (coordinate>=upper-halo .and. coordinate<upper+halo)) then
    near_boundary=.true.
    return
  endif
enddo

end function position_near_owner_boundary


subroutine rebuild_uid_index()
use gems_errors, only: werr
type(atom_dclist),pointer :: list_atom
integer                   :: i, sys_id, uid

if(allocated(uid_to_sys)) deallocate(uid_to_sys)
allocate(uid_to_sys(max(1,next_atom_uid)))
uid_to_sys(:)=0

list_atom=>sys%alist
do i=1,sys%nat
  list_atom=>list_atom%next
  uid=list_atom%o%uid
  call werr('Physical atom has no valid UID',uid<1 .or. uid>size(uid_to_sys))
  sys_id=list_atom%o%gid(sys)
  call werr('Duplicate physical atom UID',uid_to_sys(uid)/=0)
  uid_to_sys(uid)=sys_id
enddo

end subroutine rebuild_uid_index


subroutine refresh_domain_ownership()
! Recompute Cartesian ownership from current replicated positions and rebuild
! the local focus group consumed by the neighbor view.
type(atom_dclist), pointer :: la
integer                   :: i, gid
integer                   :: new_owner
logical                   :: rebuild_focus, desired_candidate, desired_focus

if(.not.allocated(atom_owner)) return

rebuild_focus = focus_dirty .or. num_images() > 1

if(num_images()==1) then
  la => sys%alist
  do i=1,sys%nat
    la => la%next
    gid=la%o%gid(sys)
    new_owner = owner_image_for_position(la%o%pos)
    if(atom_owner(gid) /= new_owner) rebuild_focus = .true.
    atom_owner(gid) = new_owner
  enddo
endif

if(.not.rebuild_focus) then
  return
endif

! Positions may have changed since the previous cell sort. Invalidate the cell
! chains before removing the comparatively few atoms whose role changed.
hs%candidate%tessellated=.false.

if(num_images()==1 .and. allocated(atom_visible)) atom_visible(:) = .true.
la => sys%alist
do i=1,sys%nat
  la => la%next
  gid=la%o%gid(sys)
  desired_candidate=atom_visible(gid) .and. candidate_atom(la%o)
  desired_focus=la%o%sym /= 'CG' .and. atom_owner(gid) == this_image()

  if(desired_candidate) then
    if(la%o%gri(hs_candidate)<1) call hs_candidate%attach(la%o)
    if(la%o%gri(hs%candidate)<1) call hs%candidate%attach(la%o)
  else
    if(la%o%gri(hs%candidate)>0) call hs%candidate%detach(la%o)
    if(la%o%gri(hs_candidate)>0) call hs_candidate%detach(la%o)
  endif

  if(desired_focus) then
    if(la%o%gri(hs_focus)<1) call hs_focus%attach(la%o)
    if(la%o%gri(hs%focus)<1) call hs%focus%attach(la%o)
  else
    if(la%o%gri(hs%focus)>0) call hs%focus%detach(la%o)
    if(la%o%gri(hs_focus)>0) call hs_focus%detach(la%o)
  endif
enddo

hs%listed = .false.
focus_dirty = .false.

sync all

end subroutine refresh_domain_ownership


subroutine solve_overlaps_domain()
! Preserve the historical recursive solver in serial. Multiple images use
! bounded rounds in which each image can mutate only the atoms it owns.
call prof_start(prof_overlap_tick)
if(num_images() == 1) then
  call overlap_moveback(hs)
  call sync_owned_positions()
  call prof_stop(prof_overlap_tick, prof_overlap, prof_n_overlap)
  return
endif

call solve_overlaps_parallel(hs)
call prof_stop(prof_overlap_tick, prof_overlap, prof_n_overlap)

end subroutine solve_overlaps_domain


subroutine sync_execution_views()
call prof_start(prof_sync_views_tick)
call gems_sync_execution_views()
call prof_stop(prof_sync_views_tick, prof_sync_views, prof_n_sync_views)

end subroutine sync_execution_views


subroutine attach_all_focus_atoms()
type(atom_dclist), pointer :: la
integer                   :: i

call hs_focus%detach_all()
call hs_candidate%detach_all()

la => sys%alist
do i=1,sys%nat
  la => la%next
  call hs_candidate%attach(la%o)
  if(la%o%sym == 'CG') cycle
  call hs_focus%attach(la%o)
enddo

hs%listed = .false.

end subroutine attach_all_focus_atoms


subroutine solve_overlaps_parallel(g)
use gems_errors, only: werr
class(ngroup),intent(inout) :: g
logical,allocatable         :: rollback(:), deposition_contact(:)
integer                     :: round, img, local_actions, total_actions
integer                     :: local_unresolved, total_unresolved

allocate(rollback(sys%nat),deposition_contact(sys%nat))
total_actions=0
total_unresolved=0

do round=1,max(1,sys%nat)
  call sync_execution_views()
  call plan_overlap_round(g,rollback,deposition_contact,local_unresolved)
  call apply_overlap_round(rollback,deposition_contact,local_actions)

  overlap_round_actions=local_actions
  overlap_round_unresolved=local_unresolved
  sync all

  total_actions=0
  total_unresolved=0
  do img=1,num_images()
    total_actions=total_actions+overlap_round_actions[img]
    total_unresolved=total_unresolved+overlap_round_unresolved[img]
  enddo

  overlap_parallel_rounds=overlap_parallel_rounds+1
  if(total_actions==0) exit
  call sync_owned_positions()
enddo

call werr('Bounded parallel overlap solver did not converge',total_actions>0)
overlap_parallel_residuals=overlap_parallel_residuals+total_unresolved
choques3=choques3+local_unresolved
deallocate(rollback,deposition_contact)

end subroutine solve_overlaps_parallel


subroutine plan_overlap_round(g,rollback,deposition_contact,unresolved)
use gems_groups, only: vdistance
class(ngroup),intent(in)  :: g
logical,intent(out)       :: rollback(:), deposition_contact(:)
integer,intent(out)       :: unresolved
type(atom_dclist),pointer :: la
type(atom),pointer        :: o1, o2, o2_image
integer,allocatable       :: row_gid(:)
integer                   :: i, ii, j, jj, k, gid1, gid2, real_row, target_gid
real(dp)                  :: vd(3), dr, moved1, moved2, scale
logical                   :: movable1, movable2

rollback(:)=.false.
deposition_contact(:)=.false.
unresolved=0
allocate(row_gid(size(g%real_id)))
row_gid(:)=0

la=>sys%alist
do k=1,sys%nat
  la=>la%next
  i=la%o%gid(g)
  if(i>0) row_gid(i)=k
enddo

la=>g%focus%alist
do ii=1,g%focus%nat
  la=>la%next
  o1=>la%o
  i=o1%gid(g)
  if(i<1) cycle
  real_row=g%real_id(i)
  if(real_row<1 .or. real_row>size(row_gid)) cycle
  gid1=row_gid(real_row)
  if(gid1<1) cycle

  do jj=1,g%nn(i)
    j=g%list(i,jj)
    o2_image=>g%a(j)%o
    if(g%b_limbo) then
      if(associated(o2_image,target=g%limbo)) cycle
    endif
    call vdistance(vd,o1,o2_image,.true.)
    dr=dot_product(vd,vd)
    if(dr>g%rcut2) cycle

    real_row=g%real_id(j)
    if(real_row<1 .or. real_row>size(row_gid)) cycle
    gid2=row_gid(real_row)
    if(gid2<1 .or. gid1==gid2) cycle
    if(associated(o2_image%prime)) then
      o2=>o2_image%prime
    else
      o2=>o2_image
    endif

    if(o2%z==2) then
      if(o1%z==1) deposition_contact(gid1)=.true.
      cycle
    endif

    moved1=dot_product(o1%pos-o1%old_cg,o1%pos-o1%old_cg)
    moved2=dot_product(o2%pos-o2%old_cg,o2%pos-o2%old_cg)
    scale=max(1._dp,dot_product(o1%pos,o1%pos),dot_product(o2%pos,o2%pos))
    movable1=o1%z==1 .and. moved1>epsilon(scale)*scale
    movable2=o2%z==1 .and. moved2>epsilon(scale)*scale

    target_gid=0
    if(movable1 .and. .not.movable2) then
      target_gid=gid1
    elseif(movable2 .and. .not.movable1) then
      target_gid=gid2
    elseif(movable1 .and. movable2) then
      if(moved1>moved2) then
        target_gid=gid1
      elseif(moved2>moved1) then
        target_gid=gid2
      else
        target_gid=max(gid1,gid2)
      endif
    endif

    if(target_gid==gid1) then
      rollback(gid1)=.true.
    elseif(target_gid==0 .and. gid1<gid2) then
      unresolved=unresolved+1
    endif
  enddo
enddo

deallocate(row_gid)

end subroutine plan_overlap_round


subroutine apply_overlap_round(rollback,deposition_contact,actions)
logical,intent(in)        :: rollback(:), deposition_contact(:)
integer,intent(out)       :: actions
type(atom_dclist),pointer :: la
type(atom),pointer        :: a
integer                   :: gid
logical                   :: do_rollback

actions=0
la=>sys%alist
do gid=1,sys%nat
  la=>la%next
  if(atom_owner(gid)/=this_image()) cycle
  a=>la%o
  do_rollback=rollback(gid)

  if(deposition_contact(gid) .and. a%z==1) then
    try=try+1
    if(ran(idum)<prob) then
      depo=depo+1
      call a%setz(3)
      do_rollback=.false.
      actions=actions+1
    else
      do_rollback=.true.
    endif
  endif

  if(.not.do_rollback) cycle
  a%pos(:)=a%old_cg(:)
  a%vel(:)=0._dp
  a%skip=.false.
  choques=choques+1
  actions=actions+1
enddo

end subroutine apply_overlap_round


subroutine classify_overlap_pairs(g)
! Build a mutation-free list of unique physical overlaps. Neighbor rows for
! periodic ghosts are mapped back to their real atom before assigning a pair.
use gems_groups, only: vdistance
class(ngroup),intent(in)   :: g
type(atom_dclist),pointer  :: la
type(atom),pointer         :: o1, o2
integer,allocatable        :: row_gid(:)
integer                   :: i, ii, j, jj, k, gid1, gid2, gid_lo, gid_hi
integer                   :: max_pairs, pair_owner, real_row
real(dp)                  :: vd(3), dr
logical                   :: duplicate_pair

max_pairs=1
if(allocated(g%focus_id)) then
  if(size(g%focus_id)>0) max_pairs=max(1,sum(g%nn(g%focus_id)))
endif

if(allocated(overlap_gid_lo)) deallocate(overlap_gid_lo, overlap_gid_hi, overlap_pair_owner)
if(allocated(overlap_cross_domain)) deallocate(overlap_cross_domain, overlap_deposition_pair)
allocate(overlap_gid_lo(max_pairs), overlap_gid_hi(max_pairs), overlap_pair_owner(max_pairs))
allocate(overlap_cross_domain(max_pairs), overlap_deposition_pair(max_pairs))
allocate(row_gid(size(g%real_id)))
row_gid(:)=0
n_overlap_pairs=0

la => sys%alist
do k=1,sys%nat
  la => la%next
  i=la%o%gid(g)
  if(i>0) row_gid(i)=k
enddo

la => g%focus%alist
do ii=1,g%focus%nat
  la => la%next
  o1 => la%o
  i=o1%gid(g)
  if(i<1) cycle
  real_row=g%real_id(i)
  if(real_row<1 .or. real_row>size(row_gid)) cycle
  gid1=row_gid(real_row)
  if(gid1<1) cycle

  do jj=1,g%nn(i)
    j=g%list(i,jj)
    o2 => g%a(j)%o
    if(g%b_limbo) then
      if(associated(o2,target=g%limbo)) cycle
    endif

    call vdistance(vd,o1,o2,.true.)
    dr=dot_product(vd,vd)
    if(dr>g%rcut2) cycle

    real_row=g%real_id(j)
    if(real_row<1 .or. real_row>size(row_gid)) cycle
    gid2=row_gid(real_row)
    if(gid2<1 .or. gid1==gid2) cycle
    gid_lo=min(gid1,gid2)
    gid_hi=max(gid1,gid2)

    duplicate_pair=.false.
    do k=1,n_overlap_pairs
      if(overlap_gid_lo(k)==gid_lo .and. overlap_gid_hi(k)==gid_hi) then
        duplicate_pair=.true.
        exit
      endif
    enddo
    if(duplicate_pair) then
      overlap_duplicate_observations=overlap_duplicate_observations+1
      cycle
    endif

    n_overlap_pairs=n_overlap_pairs+1
    pair_owner=atom_owner(gid_lo)
    overlap_gid_lo(n_overlap_pairs)=gid_lo
    overlap_gid_hi(n_overlap_pairs)=gid_hi
    overlap_pair_owner(n_overlap_pairs)=pair_owner
    overlap_cross_domain(n_overlap_pairs)=atom_owner(gid_lo)/=atom_owner(gid_hi)
    overlap_deposition_pair(n_overlap_pairs)=o1%z==2 .or. o2%z==2
  enddo
enddo

if(.not.allocated(overlap_pairs_by_owner)) then
  allocate(overlap_pairs_by_owner(num_images()))
  overlap_pairs_by_owner(:)=0
endif
overlap_conflict_scans=overlap_conflict_scans+1
overlap_pair_observations=overlap_pair_observations+n_overlap_pairs
overlap_local_observations=overlap_local_observations+count(.not.overlap_cross_domain(:n_overlap_pairs))
overlap_cross_observations=overlap_cross_observations+count(overlap_cross_domain(:n_overlap_pairs))
overlap_deposition_observations=overlap_deposition_observations+count(overlap_deposition_pair(:n_overlap_pairs))
do k=1,n_overlap_pairs
  pair_owner=overlap_pair_owner(k)
  overlap_pairs_by_owner(pair_owner)=overlap_pairs_by_owner(pair_owner)+1
enddo

deallocate(row_gid)

end subroutine classify_overlap_pairs


subroutine ensure_sync_storage()
type(atom_dclist), pointer :: la
integer, allocatable       :: old_owner(:)
integer                    :: i, nold

if(allocated(atom_owner) .and. allocated(co_uid) .and. allocated(co_element) .and. &
  allocated(atom_visible) .and. allocated(co_skip) .and. allocated(co_pos) .and. &
  allocated(co_old_cg) .and. allocated(co_vel)) then
  if(size(atom_owner) >= sys%amax .and. size(atom_visible) >= sys%amax .and. &
    size(co_pos,2) == sys%nat) return
endif

focus_dirty = .true.

nold = 0
if(allocated(atom_owner)) then
  nold = size(atom_owner)
  allocate(old_owner(nold))
  old_owner(:) = atom_owner(:)
  deallocate(atom_owner)
endif
if(allocated(atom_visible)) deallocate(atom_visible)

sync all
if(allocated(co_element)) deallocate(co_element)
if(allocated(co_uid)) deallocate(co_uid)
if(allocated(co_new_owner)) deallocate(co_new_owner)
if(allocated(co_skip)) deallocate(co_skip)
if(allocated(co_pos)) deallocate(co_pos)
if(allocated(co_old_cg)) deallocate(co_old_cg)
if(allocated(co_vel)) deallocate(co_vel)
if(allocated(recv_pos)) deallocate(recv_pos)
if(allocated(recv_old_cg)) deallocate(recv_old_cg)
if(allocated(recv_vel)) deallocate(recv_vel)
if(allocated(recv_element)) deallocate(recv_element)
if(allocated(recv_uid)) deallocate(recv_uid)
if(allocated(recv_owner)) deallocate(recv_owner)
if(allocated(uid_to_sys)) deallocate(uid_to_sys)
if(allocated(recv_skip)) deallocate(recv_skip)

allocate(atom_owner(size(sys%a)))
allocate(co_element(sys%nat)[*])
allocate(co_uid(sys%nat)[*])
allocate(co_new_owner(sys%nat)[*])
allocate(atom_visible(size(sys%a)))
allocate(co_skip(sys%nat)[*])
allocate(co_pos(3,sys%nat)[*])
allocate(co_old_cg(3,sys%nat)[*],co_vel(3,sys%nat)[*])
allocate(recv_pos(3,sys%nat),recv_old_cg(3,sys%nat),recv_vel(3,sys%nat))
allocate(recv_uid(sys%nat),recv_element(sys%nat),recv_owner(sys%nat),recv_skip(sys%nat))
co_pos(:,:) = 0._dp
co_old_cg(:,:) = 0._dp
co_vel(:,:) = 0._dp
co_element(:) = 119
co_uid(:) = 0
co_new_owner(:) = 1
co_skip(:) = .false.
atom_visible(:) = .false.

atom_owner(:) = 1
if(nold > 0) atom_owner(:min(nold,size(atom_owner))) = old_owner(:min(nold,size(atom_owner)))

if(sys%nat > nold) then
  la => sys%alist
  do i=1,sys%nat
    la => la%next
    if(la%o%gid(sys) > nold) atom_owner(la%o%gid(sys)) = owner_image_for_position(la%o%pos)
  enddo
endif

if(allocated(old_owner)) deallocate(old_owner)
sync all

end subroutine ensure_sync_storage


! Reservorios
subroutine gcmc_run(g)
! Domain-local GCMC. Proposals are generated against a virtual local state,
! then accepted create/delete records are replayed in the same global order on
! every image. Safety margins keep concurrent trial insertions noninteracting.
use gems_program_types, only: distance
use gems_errors, only: werr
class(group)               :: g
real(dp)                   :: r(3), vd(3), dr, rc
real(dp)                   :: active_lo(3), active_hi(3), active_volume
real(dp)                   :: total_active_volume, raw_attempts, my_fraction
real(dp)                   :: other_raw, other_fraction, eps_fraction
real(dp),allocatable       :: eligible_pos(:,:)
type(atom_dclist), pointer :: la
type(atom),pointer         :: o, ref
integer,allocatable        :: eligible_uid(:)
integer                    :: i,j,n,m,owner,kind,uid
integer                    :: base_attempts,total_base,remaining,fraction_rank
integer                    :: uid_base,uid_prefix,total_creations
integer                    :: max_eligible
logical                    :: overlaps
 
rc=hs%rcut

call gcmc_accessible_bounds(active_lo,active_hi,gcmc_patch_volume,active_volume)
gcmc_accessible_volume=active_volume
co_gcmc_active_volume=active_volume
sync all

total_active_volume=0._dp
do owner=1,num_images()
  total_active_volume=total_active_volume+co_gcmc_active_volume[owner]
enddo
call werr('GCMC safety margins leave no accessible reservoir volume',total_active_volume<=0._dp)

! Keep exactly nadj global attempts while distributing them in proportion to
! accessible volume. Largest-remainder assignment avoids stochastic rounding.
raw_attempts=real(nadj,dp)*active_volume/total_active_volume
base_attempts=floor(raw_attempts)
total_base=0
do owner=1,num_images()
  total_base=total_base+floor(real(nadj,dp)*co_gcmc_active_volume[owner]/total_active_volume)
enddo
remaining=nadj-total_base
my_fraction=raw_attempts-real(base_attempts,dp)
fraction_rank=1
eps_fraction=32._dp*epsilon(1._dp)
do owner=1,num_images()
  if(owner==this_image()) cycle
  other_raw=real(nadj,dp)*co_gcmc_active_volume[owner]/total_active_volume
  other_fraction=other_raw-real(floor(other_raw),dp)
  if(other_fraction>my_fraction+eps_fraction .or. &
    (abs(other_fraction-my_fraction)<=eps_fraction .and. owner<this_image())) fraction_rank=fraction_rank+1
enddo
gcmc_local_attempts=base_attempts
if(fraction_rank<=remaining) gcmc_local_attempts=gcmc_local_attempts+1

co_gcmc_count=0
co_gcmc_creations=0
co_gcmc_kind(:)=0
co_gcmc_uid(:)=0
co_gcmc_pos(:,:)=0._dp

max_eligible=max(1,g%nat+gcmc_local_attempts)
allocate(eligible_uid(max_eligible),eligible_pos(3,max_eligible))

! Build the particles eligible for the reverse deletion proposal. Particles
! outside the accessible region remain overlap obstacles but cannot be deleted.
n=0
la=>g%alist
do j=1,g%nat
  la=>la%next
  if(owner_image_for_position(la%o%pos)/=this_image()) cycle
  if(.not.position_in_box(la%o%pos,active_lo,active_hi)) cycle
  n=n+1
  eligible_uid(n)=la%o%uid
  eligible_pos(:,n)=la%o%pos(:)
enddo

! Generate local operations without changing group/list order. This lets every
! image replay the complete operation stream identically after one exchange.
adj: do i=1,gcmc_local_attempts
  ref => g%alist%next%o
  call werr('No more particles',.not.associated(ref))

  ! Creation attempt
  if (ran(idum)<0.5) then
              
    ! Metropolis acceptance
    if(act*active_volume/(n+1)<ran(idum)) cycle
       
    ! Random coordinates
    do j=1,3
      r(j)=active_lo(j)+ran(idum)*(active_hi(j)-active_lo(j))
    enddo
                   
    ! Check overlap against the current group, excluding virtual deletions.
    overlaps=.false.
    la=>g%alist
    do j=1,g%nat
      la=>la%next
      o => la%o
      if(gcmc_uid_deleted(o%uid,co_gcmc_count)) cycle
      vd(:) = distance(o%pos,r,o%pbc)
      dr = dot_product(vd,vd)
      if(dr<rc*rc) then
        overlaps=.true.
        exit
      endif
    enddo
    if(overlaps) cycle

    ! Provisional creations are not in the physical groups yet.
    do j=1,n
      if(eligible_uid(j)>=0) cycle
      vd(:)=distance(eligible_pos(:,j),r,[.true.,.true.,.false.])
      if(dot_product(vd,vd)<rc*rc) then
        overlaps=.true.
        exit
      endif
    enddo
    if(overlaps) cycle

    co_gcmc_count=co_gcmc_count+1
    co_gcmc_creations=co_gcmc_creations+1
    co_gcmc_kind(co_gcmc_count)=1
    co_gcmc_uid(co_gcmc_count)=-co_gcmc_creations
    co_gcmc_pos(:,co_gcmc_count)=r(:)
    n=n+1
    eligible_uid(n)=-co_gcmc_creations
    eligible_pos(:,n)=r(:)
         
  ! Destruction attempt
  else 
            
    ! Metropolis acceptance
    if(real(n,dp)/(active_volume*act)<ran(idum)) cycle
                
    ! Choose a particle
    m=floor(ran(idum)*n)+1
    if(m>n) m=n
    call werr('Chosen particle does not exist',m<1 .or. m>n)

    co_gcmc_count=co_gcmc_count+1
    co_gcmc_kind(co_gcmc_count)=-1
    co_gcmc_uid(co_gcmc_count)=eligible_uid(m)
    co_gcmc_pos(:,co_gcmc_count)=eligible_pos(:,m)
    n=n-1
    if(m<=n) then
      eligible_uid(m:n)=eligible_uid(m+1:n+1)
      eligible_pos(:,m:n)=eligible_pos(:,m+1:n+1)
    endif
       
  endif 
enddo adj

! Assign compact, globally unique UIDs by image-prefix after proposal counts
! are known. Temporary negative UIDs also let a later local deletion refer to a
! creation from the same batch.
sync all
uid_base=next_atom_uid
uid_prefix=0
do owner=1,this_image()-1
  uid_prefix=uid_prefix+co_gcmc_creations[owner]
enddo
total_creations=uid_prefix
do owner=this_image(),num_images()
  total_creations=total_creations+co_gcmc_creations[owner]
enddo
do i=1,co_gcmc_count
  if(co_gcmc_uid(i)<0) co_gcmc_uid(i)=uid_base+uid_prefix-co_gcmc_uid(i)
enddo
next_atom_uid=uid_base+total_creations
sync all

do owner=1,num_images()
  do i=1,co_gcmc_count[owner]
    kind=co_gcmc_kind(i)[owner]
    uid=co_gcmc_uid(i)[owner]
    r(:)=co_gcmc_pos(:,i)[owner]
    call werr('Invalid GCMC operation kind',kind/=-1 .and. kind/=1)
    if(kind==1) then
      call gcmc_apply_creation(g,uid,r)
    else
      call gcmc_apply_deletion(uid)
    endif
  enddo
enddo

call ensure_sync_storage()
la=>sys%alist
do i=1,sys%nat
  la=>la%next
  atom_owner(la%o%gid(sys))=owner_image_for_position(la%o%pos)
enddo
focus_dirty=.true.

deallocate(eligible_uid,eligible_pos)

end subroutine gcmc_run


subroutine gcmc_accessible_bounds(active_lo,active_hi,patch_volume,active_volume)
real(dp),intent(out) :: active_lo(3), active_hi(3), patch_volume, active_volume
real(dp)             :: reservoir_lo(3), reservoir_hi(3), margin, scale
integer              :: axis

reservoir_lo=[o,o,z0]
reservoir_hi=[xi,yi,zmax]
active_lo=max(domain_lo,reservoir_lo)
active_hi=min(domain_hi,reservoir_hi)
patch_volume=product(max(0._dp,active_hi-active_lo))
margin=0.5_dp*hs%rcut
scale=max(1._dp,maxval(abs(domain_global_hi-domain_global_lo)))

do axis=1,3
  if(active_hi(axis)<=active_lo(axis)) cycle
  if(domain_dims(axis)==1) cycle
  if(abs(active_lo(axis)-domain_lo(axis))<=epsilon(scale)*scale) then
    if(domain_periodic(axis) .or. domain_coord(axis)>0) active_lo(axis)=active_lo(axis)+margin
  endif
  if(abs(active_hi(axis)-domain_hi(axis))<=epsilon(scale)*scale) then
    if(domain_periodic(axis) .or. domain_coord(axis)<domain_dims(axis)-1) active_hi(axis)=active_hi(axis)-margin
  endif
enddo

active_volume=product(max(0._dp,active_hi-active_lo))

end subroutine gcmc_accessible_bounds


logical function position_in_box(r,lo,hi) result(inside)
real(dp),intent(in) :: r(3),lo(3),hi(3)

inside=all(r>=lo) .and. all(r<=hi)

end function position_in_box


logical function gcmc_uid_deleted(uid,nops) result(deleted)
integer,intent(in) :: uid,nops
integer            :: op

deleted=.false.
do op=1,nops
  if(co_gcmc_kind(op)==-1 .and. co_gcmc_uid(op)==uid) then
    deleted=.true.
    return
  endif
enddo

end function gcmc_uid_deleted


subroutine gcmc_apply_creation(g,uid,r)
use gems_groups, only: atom_asign
use gems_errors, only: werr
class(group)               :: g
integer,intent(in)         :: uid
real(dp),intent(in)        :: r(3)
type(atom),pointer         :: created,ref
class(group),pointer       :: gp
integer                    :: membership

ref=>g%alist%next%o
call werr('No GCMC template particle',.not.associated(ref))
allocate(created)
call created%init()
created%uid=uid
call atom_asign(created,ref)
created%pos(:)=r(:)
do membership=1,ref%ngr
  gp=>ref%gro(membership)
  call gp%attach(created)
enddo
created=>null()

end subroutine gcmc_apply_creation


subroutine gcmc_apply_deletion(uid)
use gems_errors, only: werr
integer,intent(in)         :: uid
type(atom_dclist),pointer  :: list_atom
type(atom),pointer         :: deleted
integer                    :: item

deleted=>null()
list_atom=>sys%alist
do item=1,sys%nat
  list_atom=>list_atom%next
  if(list_atom%o%uid/=uid) cycle
  deleted=>list_atom%o
  exit
enddo
call werr('GCMC deletion UID does not exist',.not.associated(deleted))
call deleted%dest()
deallocate(deleted)

end subroutine gcmc_apply_deletion

 
! Integradores de posics.
subroutine cbrownian_hs(g,h)
! Integrac. browniana
class(ngroup)              :: g
class(atom), pointer        :: o1 
type(atom_dclist), pointer :: la 
real(dp),intent(in)        :: h
real(dp)                   :: fac1, r1, posold, vd(3)
integer                    :: i,j,ii
logical                    :: depos

la => g%focus%alist
do ii = 1,g%focus%nat
  la => la%next
  o1 => la%o 

  !Para luego ver congelam.
  o1%old_cg(:)=o1%pos(:)

  ! Para cambio en coef. difusión SC-SEI
  if (o1%pos(3)>80._dp) then
     dif= dif_sc
  else
     dif= dif_sei
  endif

  fac1 = sqrt(2._dp*dif*h) 

  do j = 1,3
    r1=gasdev()
  
    posold = o1%pos(j)
    o1%pos(j) = posold + r1*fac1

    ! Velocidad derivada de euler para atras
    o1%vel(j) = (o1%pos(j)-posold)/h
  enddo

  call atom_pbc(o1, depos)
  ! Para no calcular choques si depositó
  if (depos) cycle

  max_vel=max(max_vel,dot_product(o1%vel,o1%vel))

  ! Mark atom to check for colisions
  o1%skip=.false.
enddo

            
end subroutine cbrownian_hs


recursive subroutine overlap_moveback(g)
! Search for colisions and try to solve them by a sequence of moving back
! the particles to its previous positions.
! NOTE: When using piston, overlaps due to reservoir compresion may remain
! unsolved.
use gems_groups, only: vdistance
use gems_errors, only: werr
class(ngroup)              :: g
class(atom),pointer        :: o1,o2
real(dp)                   :: ne, vd(3), dr
integer                    :: i, ii, j, jj
type(atom_dclist), pointer :: la 
logical                    :: again
integer, save              :: count=0

again=.false.

! Choque con las demas particulas
la => g%focus%alist
do ii = 1,g%focus%nat
  la => la%next
  o1 => la%o 
     
  if(o1%skip) cycle
  o1%skip=.true.

  i = o1%gid(g)
  do jj = 1, g%nn(i)  ! sobre los vecinos

    j = g%list(i,jj)
    o2 => g%a(j)%o
               
    ! Skip atoms in limbo
    if(g%b_limbo) then
      if(associated(o2,target=g%limbo)) cycle
    endif
               
    call vdistance(vd,o1,o2,.true.)
    dr = dot_product(vd,vd)

    if(dr>g%rcut2) cycle !Sí es necesario :B

    ! Deposicion por contacto con otra particula metalica
    ! 1.19 es el radio de Mayers

    if (o2%z==2) then
      ! if(dr> 1.4161_dp) cycle
      
      try= try + 1
      ne=ran(idum)
      if(ne<prob) then
        depo= depo + 1 ! para contar CG
        call o1%setz(3)   ! F. Es inerte en este paso de tiempo

        ! Permite deposicion en cadena pero depende del atom id.
        ! Si queremos deposicion en cadnea sería mejor programarla
        ! para que no dependa de el orden en que se ejecuta el do.
        if (o1%pos(3)>z0) then
           print *,'superó z0', o1%pos(3)
           stop ! raro 
        endif
      else
        o1%pos(:)= o1%old_cg(:) 
        ! TODO: resort?
        o1%skip=.false.
      endif

      exit

    endif

    ! Colision con otra particula
    o2%pos(:)=o2%old_cg(:) 
    ! TODO: resort?
    o2%vel(:)=0._dp
    o2%skip=.false.
    choques=choques+1
    again=.true.
  enddo
    
enddo

i=choques
if(again) then
  count=count+1
  call werr('Recursion limit in HS overlap solving',count==100)
  call overlap_moveback(g)
else
  count=0
endif
choques2=max(choques2,choques-i)

end subroutine overlap_moveback




!Output
subroutine salida()  ! Escribe los datos calc. :P
integer  :: j, image, total_try, total_depo
type(atom_dclist), pointer :: la 

call sync_owned_positions(.true.)
output_try=try
output_depo=depo
sync all

total_try=0
total_depo=0
do image=1,num_images()
  total_try=total_try+output_try[image]
  total_depo=total_depo+output_depo[image]
enddo
try=0
depo=0

if (this_image() /= 1) return

! Coords. de partíc.
write(11,*) sys%nat ! n
write(11,*) "info:",zmax,sys%nat
la=>sys%alist
do j=1,sys%nat
 la=>la%next 
 pa=>la%o
 write(11,*) pa%sym,pa%pos(:),pa%z

enddo

write(15,*)t,total_try
write(16,*)t,total_depo

flush(16)
flush(15)
flush(11)

end subroutine salida


! Mas subrutinas:
subroutine atom_pbc(o1, depos)
! Para PBC, e intento deposic. sobre electrodo
use gems_program_types, only: box
class(atom),pointer :: o1
integer     :: j
real(dp)    :: ne
logical     :: depos

depos= .false.

do j=1, 3
  if(j<3) then
    ! PBC en x e y
    if (o1%pos(j)>box(j)) then
       o1%pos(j)=o1%pos(j)-box(j)
    endif  
    if (o1%pos(j)<o) then
       o1%pos(j)=o1%pos(j)+box(j)
    endif

  else 
    ! Rebote en zmax
    if(o1%pos(j)>zmax) then
      o1%pos(:)= o1%old_cg(:)
    endif

  endif
end do

! En una componente (xt - x0)**2 -- mean square displacement
msd_u= o1%vel(1) * o1%vel(1) * h * h
msd_t= msd_t + msd_u


! Si toca el electrodo implicito ¿se congela? (probabilidad ne)
if(o1%pos(3)<=0._dp) then 
  try= try + 1
  ne=ran(idum)

  if (electrodo.eqv..false.) then
    if(ne<prob) then
      depo= depo + 1
      call o1%setz(3)   ! Es inerte en este paso de tiempo
      depos= .true.
    endif
  endif

  ! XXX: atomos depositados también rebotan. Sería + consistente
  ! si no lo hicieran (ver cbrownian_hs)
  o1%pos(:) = o1%old_cg(:)
  ! FIXME: con el XXX de arriba se debería chequear colisiones tb. en deposic.

endif
end subroutine atom_pbc




function gasdev() !Nro aleat.
  ! XXX: We use sp in rsq because log(rsq) gives different results in different machines when rsq is dp.
  real(sp)                  :: rsq 
  real(dp)                  :: v1,v2
  real(dp), save            :: g
  real(dp)                  :: gasdev
  logical, save             :: gaus_stored=.false.
  if (gaus_stored) then
    gasdev=g
    gaus_stored=.false.
  else
    do
      v1=2.0_dp*ran(idum)-1.0_dp
      v2=2.0_dp*ran(idum)-1.0_dp
      rsq=v1**2+v2**2
      if (rsq > 0._dp .and. rsq < 1._dp) then
             exit
      end if
    end do
    rsq=sqrt(-2.0_dp*log(rsq)/rsq)
    gasdev=v1*real(rsq,kind=dp)
    g=v2*rsq
    gaus_stored=.true.
  end if

end function gasdev
 

function ran(idum)
implicit none
integer, parameter     :: k4b=selected_int_kind(9)
integer, intent(inout) :: idum
real(dp)               :: ran
integer, parameter     :: ia=16807,im=2147483647,iq=127773,ir=2836
real(dp), parameter    :: am=nearest(1.0,-1.0)/real(im,kind=dp)
integer(k4b)           :: k
integer(k4b), save    :: ix=-1,iy=-1
if (idum <= 0 .or. iy < 0) then
  iy=ior(ieor(888889999,abs(idum)),1)
  ix=ieor(777755555,abs(idum))
  idum=abs(idum)+1
end if
ix=ieor(ix,ishft(ix,13))
ix=ieor(ix,ishft(ix,-17))
ix=ieor(ix,ishft(ix,5))
k=iy/iq
iy=ia*(iy-k*iq)-ir*k
if (iy < 0) iy=iy+im
ran=am*ior(iand(im,ieor(ix,iy)),1)   
end function ran 

end program
