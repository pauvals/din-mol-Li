program dana
  use gems_groups, only: group, atom, atom_dclist, sys, gindex
  use gems_neighbor, only: gems_sync_execution_views => sync_execution_views, ngroup, nn_vlist, nupd_vlist, nb_dcut
  use gems_elements, only: set_z, elements_init,elements
  use gems_errors, only: timer_dump, timer_start, sclock_t1, sclock_t2, wstd, sclock_max, sclock_rate,logunit
  use gems_constants, only: time1, dp, sp

  implicit none
 
  ! Cosas del sistema  #what a quilombo
  integer             :: n, nx, nchunk        ! Num de particulas
  real(dp), parameter :: o=0._dp, tau=0.1_dp  ! Origen de la caja, tau p/ rho
  real(dp)            :: z0, z1, zmax         ! Para ajuste del reservorio
  real(dp)            :: xi, yi               ! Tamaño en x e y de caja sim.
  real(dp)            :: domain_lo(3)[*], domain_hi(3)[*]
  real(dp)            :: domain_global_lo(3), domain_global_hi(3)
  integer             :: domain_dims(3)=[1,1,1], domain_coord(3)=[0,0,0]
  logical             :: domain_periodic(3)=[.true.,.true.,.false.]
  character(8)        :: domain_topology='z'
  real(dp), parameter :: gama=1._dp           ! Para fricción
  real(dp)            :: dif, dif_sei, dif_sc ! Difusion
  real(dp), parameter :: Tsist=300._dp        ! T= 300 K
  real(dp)            :: eps(3,3),r0(3,3)     ! Lennard-Jones epsilon, r0 

  ! Constantes Físicas y
  ! factores de conversión de unidades
  real(dp), parameter :: kB_eVK=8.617330350e-5_dp  ! cte. Boltzmann en eV
  real(dp), parameter :: kJm_ui=100._dp, eV_kJm=96.485_dp,eV_ui=eV_kJm*kJm_ui
  real(dp), parameter :: kB_ui=kB_eVK*eV_ui
  real(dp), parameter :: kT_ui=kB_ui*300._dp  !en ui, la T= 300 K
 
  ! Variables dinámicas, agrupadas en átomo
  type(atom), allocatable, target :: a(:) ! átomo
  type(atom), pointer             :: pa=>null(), pb=>null()
  real(dp)                        :: prob,max_vel=0._dp, msd_u= 0._dp, msd_t= 0._dp, msd_max= 0._dp

  ! Contador de choques (superposiciones)
  integer             :: choques=0, choques2=0, choques3=0

  ! Parametros de integración
  real(dp)            :: h          ! paso de tiempo
  integer             :: nst        ! nro de pasos
  integer             :: nwr        ! paso de escritura
  real(dp)            :: t=0.0_dp   ! tiempo en ps

  ! To configure the run: positions integrator and reservoir type 
  logical             :: integrador, s_piston=.false.,s_chunk=.false.,&
                         new_sim, electrodo, ideal_dep

  ! Observables
  real(dp)            :: temp,rho,rho0

  ! Varios
  integer             :: idum         ! Semilla
  integer             :: i,j,k        ! Enteros
  integer             :: try, depo    ! Parts. cercanas a electrodo y deposic.
  real(dp)            :: dist, rhomedia, cstdev, factor, factor2
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

  ! Ermak
  logical             :: ermak 
  real(dp)            :: cc0,cc1,cc2, cc0_sei,cc1_sei,cc2_sei,&
                         cc0_sc,cc1_sc,cc2_sc
  real(dp)            :: sdr,sdv, sdr_sei,sdv_sei, sdr_sc,sdv_sc, skt
  real(dp)            :: crv1,crv2, crv1_sei,crv2_sei, crv1_sc,crv2_sc
  real(dp),allocatable:: ranv(:,:)
  real(dp),allocatable:: co_pos(:,:)[:]
  real(dp),allocatable:: co_old_cg(:,:)[:], co_pos_old(:,:)[:]
  real(dp),allocatable:: co_vel(:,:)[:], co_acel(:,:)[:], co_force(:,:)[:]
  integer,allocatable :: atom_owner(:)
  integer,allocatable :: co_owner(:)[:], co_element(:)[:]
  logical,allocatable :: atom_visible(:), co_skip(:)[:]
  logical             :: focus_dirty = .false.
  integer             :: overlap_round_actions[*]=0, overlap_round_unresolved[*]=0
  integer             :: output_try[*]=0, output_depo[*]=0
  integer             :: overlap_parallel_rounds=0, overlap_parallel_residuals=0
   
  ! GCMC
  logical             :: s_gcmc=.false.
  type(group)         :: gcmc
  real(dp)            :: act
  integer             :: nadj
  type(atom_dclist), pointer :: la, lb 
  
  ! Group to add particles
  type(group)         :: chunk

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
            
  ! Lennard-Jones: epsilon,  r0 :P
  eps(:,2) = 0
  r0(:,2)  = 0
  eps(2,:) = 0
  r0(2,:)  = 0
  r0(2,1)  = 3.5_dp  
  r0(1,2)  = 3.5_dp
  eps(1,1) = 2313.6_dp 
  r0(1,1)  = 3.2_dp
  eps(3,3) = 121._dp
  r0(3,3)  = 3.61_dp
  eps(1,3) = 529.1_dp
  eps(3,1) = eps(1,3)
  r0(1,3)  = 1.564_dp
  r0(3,1)  = r0(1,3)
 
  ! Para contar partículas que intentan deposic.
  try= 0
  depo= 0
                      
  ! Init general del sistema de grupos
  call gindex%init()
  call sys%init()
 
  ! Leer run parameters:  semilla, probabilidad, nro. pasos
  ! y otros valores iniciales
  call entrada()
        
  ! Init lista para choques de esferas duras
  call hs_focus%init()
  call hs_candidate%init()
  call hs%init()
  call hs%bind_groups(hs_focus, hs_candidate)
  call hs%setrc(3.2_dp) ! Maximo radio de corte

  ! Configure integrator and reservoir for the run.
  call config_run()

  ! Generar posiciones inic.
  if (new_sim .and. this_image() == 1) call gen_pos_inic()
  sync all

  ! Define tamaño en z del sistema (sin contar el reservorio), y lee posic. inic.
  ! de partículas 
  ! Read initial positions (optionally from a previous run).
  call config_inic()
  call sync_owned_positions()

  ! Add to GCMC
  if(s_gcmc) then
    call gcmc%init()
    la=>sys%alist
    do j=1,sys%nat
     la=>la%next 
     if(la%o%z==1) call gcmc%attach(la%o)
    enddo
  endif

  !Search neighbors
  call sync_execution_views()

  if (integrador) call fuerza(hs,eps,r0)
  call sync_owned_positions()

  ! Calcula rho - densidad reservorio inicial
  call calc_rho(rho)
  rho0=rho
  

  ! Leer chunk de atoms (Reservoir 2):
  if(s_chunk)  call config_chunk()


  ! Output 
  if (this_image() == 1) then
    open(11,File='Li.xyz')
    open(12,File='E.dat')
    open(13,File='T.dat')
    open(14,File='rho.dat')
    open(15,File='try.dat')
    open(16,File='depo.dat')
  endif

  call salida() ! config. inic. en el primer paso ;)

  call timer_start()

  do i=1,nst
    ! Da un paso #wii

    call prof_start(prof_integrate_tick)
    if (integrador) then
      ! Para trabajar con din. Langevin
      call ermak_a(hs,ranv)
      ! FIXME: sync_execution_views must be always before force calculation
      ! TODO: ver de sacar ya knock
      ! Ve si congela o rebota
      ! call knock(hs)
      call fuerza(hs,eps,r0)
      call ermak_b(hs,ranv)
      call prof_stop(prof_integrate_tick, prof_integrate, prof_n_integrate)
      call sync_owned_positions()

    else
      ! Para usar din. Browniana, y ~esferas duras
      call cbrownian_hs(hs,h)
      call prof_stop(prof_integrate_tick, prof_integrate, prof_n_integrate)
      call sync_owned_positions()

    endif

    ! Update neighbors for serial overlap solving. Multi-image overlap rebuilds
    ! the all-atom view on image 1 inside solve_overlaps_domain().
    if(num_images() == 1) call sync_execution_views()
 
    ! Solve overlaps
    call solve_overlaps_domain()
       
    ! Update neighbors (give error if not used again)
    call sync_execution_views()
         
    msd_t= msd_t/hs%focus%nat
    msd_max= max(msd_max,msd_t)
 
    ! Para estudiar desposic. ideal (Cottrell):
    if (ideal_dep) then
      la=> hs%focus%alist
      do j = 1, hs%focus%nat
        la=> la%next
        o1=> la%o
        if (o1%sym/='F') cycle
        lb => la%prev
        if(s_gcmc) call gcmc%detach(o1)
        call o1%dest()
        deallocate(o1)
        la => lb
      ! Fin chequeo
      enddo
      la=> hs%focus%alist
      do j = 1, hs%focus%nat
        la=> la%next
        o1=> la%o
        if (o1%sym/='F') cycle 
        print *, 'error'
        stop
      enddo
    endif

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
      if(s_gcmc) call gcmc%detach(o1)
    enddo
    
    ! Reservoir 3 - GCMC:
    if(s_gcmc) then
      call gcmc_run(gcmc)
      call sync_owned_positions()
    endif
      
    ! Cant. de partículas en reservorio o sensor
    call calc_rho(rho)

    ! Reservoir 2 - blocks of atoms: 
    if (s_chunk) call bloques(sys, chunk, nx, rho)

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
  
    ! Reservoir 1 - piston:
    if (s_piston) call maxz(zmax)

    t=t+h
   
  enddo
 
  if (this_image() == 1) then
    close(11)
    close(12)
    close(13)
    close(14)
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
read(15,*) dist 
read(15,*) z0 
read(15,*) zmax 
read(15,*) dif_sc 
read(15,*) dif_sei 
read(15,*) nb_dcut
read(15,*) new_sim
read(15,*) electrodo
read(15,*) ideal_dep
close(15)
end subroutine entrada


subroutine gen_pos_inic()
! Create initial positions. For 1 M [Li+]
use gems_program_types, only: distance 
integer             :: n 
real(dp),allocatable:: r(:,:) !posic.
real(dp)            :: Mol, v1(3),v2(3),dif_pos(3), alto, dif2 ! molaridad, diferencia entre vectores posic., alto caja sim.
real(dp),parameter  :: r0=3.2_dp, mLi= 6.94_dp
integer             :: i,j,l,k,n_CG,n_Li
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
close(13)


end subroutine gen_pos_inic


subroutine config_run()
! Integrator and reservoir configuration 
use gems_errors, only: werr
use gems_program_types, only: tbox, box_setvars  
character(100)  :: a

! Set the box
tbox(:,:)= 0._dp
tbox(1,1)= xi
tbox(2,2)= yi
tbox(3,3)= zmax
call box_setvars()

open(15,File='movedor.ini') ! by Nohe :)
read(15,*) integrador 

! tipo de reservorio
read(15,*) a
select case(trim(a))
! Reservoir type 1
case('piston')
  s_piston=.true.

! Reservoir type 2
case('chunks')
  s_chunk=.true.
  call chunk%init()

! Reservoir type 3
case('gcmc')
  s_gcmc=.true.
  read(15,*) act, nadj
case default
  call werr('Unknown reservoir type',.true.)
end select

close(15)
end subroutine config_run


subroutine config_inic()
! configuración inicial
use gems_program_types, only: tbox, box_setvars
integer :: i
!type(atom_dclist), pointer :: la
character(10)     :: sym

! Tamaños iniciales de "reservorio"
if(s_chunk) then
  ! Al usar chunk con sensor
  dist= dist + hs%rcut
  z1 = z0 + dist
  zmax = z1 + dist 
endif

if(num_images() > 1) then
  tbox(3,3)=zmax
  call box_setvars()
endif
call set_domain_bounds()

! Para luego actualizar reservorio
rhomedia= 5.775329e-4_dp
cstdev= 1.0754306e-4_dp
factor= rhomedia - cstdev
factor2= rhomedia + cstdev
! 18.6211140525501
! 0.000577532941264052= STATS_mean
! 0.000107543058373559= 4*STATS_stddev


open(11,File='pos_inic.xyz')
read(11,*) n
read(11,*)

if(allocated(atom_owner)) deallocate(atom_owner)
if(allocated(co_owner)) deallocate(co_owner)
if(allocated(co_element)) deallocate(co_element)
if(allocated(atom_visible)) deallocate(atom_visible)
if(allocated(co_skip)) deallocate(co_skip)
if(allocated(co_pos)) deallocate(co_pos)
if(allocated(co_old_cg)) deallocate(co_old_cg)
if(allocated(co_pos_old)) deallocate(co_pos_old)
if(allocated(co_vel)) deallocate(co_vel)
if(allocated(co_acel)) deallocate(co_acel)
if(allocated(co_force)) deallocate(co_force)
allocate(atom_owner(n))
allocate(co_owner(n)[*])
allocate(co_element(n)[*])
allocate(atom_visible(n))
allocate(co_skip(n)[*])
allocate(co_pos(3,n)[*])
allocate(co_old_cg(3,n)[*],co_pos_old(3,n)[*])
allocate(co_vel(3,n)[*],co_acel(3,n)[*],co_force(3,n)[*])
  atom_owner(:)=1
  co_owner(:)=1
  co_element(:)=119
  atom_visible(:)=.true.
  co_skip(:)=.false.
  co_pos(:,:)=0._dp
  co_old_cg(:,:)=0._dp
  co_pos_old(:,:)=0._dp
  co_vel(:,:)=0._dp
  co_acel(:,:)=0._dp
  co_force(:,:)=0._dp
  focus_dirty=.true.

do i=1,n

  ! Allocatea el atomo-como un puntero del tipo 'atom'
  allocate(pa)
  call pa%init()

  read(11,*) sym,pa%pos(:)
  call pa%setsym(sym) ! Asigna algunos valores según tipo de átomo
  pa%force(:)=0._dp
  pa%pos_old(:)=pa%pos(:)

  ! Agrega atom al sistema (sys)
  call sys%attach(pa)

  ! Agrego atomos a los grupos
  ! NOTA: hs%attach tiene que estar al final
  ! para que ande la lista de vecinos
  atom_owner(i) = owner_image_for_position(pa%pos)

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

allocate(ranv(n,3))

! Calculo la velocidad neta del sistema/Sino como que se trasladaría todo el sist. en el espacio... Así trabajo c/ coords.
! internas ;)
! do k=1,3
!   vm(k)=sum(pa%pos(k)-pa%pos_old(k))/n
! enddo

! Sustraer la velocidad neta
! do i=1,n
!   pa%pos_old(i)=pa%pos_old(i)-vm(:)
! end do

!calculo vel. inic.
! pa%vel(:)=(pa%pos(:)-pa%pos_old(:))/h

! Para trabajar con din. Langevin
if (integrador) then
  ! Valores inic. de las ctes. de Ermak
  ! call set_ermak(h,gama_sc,Tsist,cc0_sc, cc1_sc, cc2_sc, sdr_sc, sdv_sc, crv1_sc, crv2_sc)
  ! call set_ermak(h,gama,Tsist,cc0_sei, cc1_sei, cc2_sei, sdr_sei, sdv_sei, crv1_sei, crv2_sei)
  call set_ermak(h,gama,Tsist,cc0, cc1, cc2, sdr, sdv, crv1, crv2)

endif

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


subroutine sync_owned_positions()
! Replicated-state bridge: each owner publishes the atom state needed by a new
! owner after migration and by halo consumers during overlap correction.
type(atom_dclist), pointer :: la
integer                   :: i, owner

call prof_start(prof_sync_pos_tick)
call ensure_sync_storage()
call mark_visible_atoms()

co_pos(:,:) = 0._dp
co_old_cg(:,:) = 0._dp
co_pos_old(:,:) = 0._dp
co_vel(:,:) = 0._dp
co_acel(:,:) = 0._dp
co_force(:,:) = 0._dp
co_element(:) = 119
co_skip(:) = .false.
la => sys%alist
do i=1,sys%nat
  la => la%next
  if(atom_owner(i) /= this_image()) cycle
  co_pos(:,i) = la%o%pos(:)
  co_old_cg(:,i) = la%o%old_cg(:)
  co_pos_old(:,i) = la%o%pos_old(:)
  co_vel(:,i) = la%o%vel(:)
  co_acel(:,i) = la%o%acel(:)
  co_force(:,i) = la%o%force(:)
  co_element(i) = la%o%z
  co_skip(i) = la%o%skip
enddo

sync all

la => sys%alist
do i=1,sys%nat
  la => la%next
  if(.not.atom_visible(i)) cycle
  owner = atom_owner(i)
  if(la%o%z /= co_element(i)[owner]) call la%o%setz(co_element(i)[owner])
  la%o%pos(:) = co_pos(:,i)[owner]
  la%o%old_cg(:) = co_old_cg(:,i)[owner]
  la%o%pos_old(:) = co_pos_old(:,i)[owner]
  la%o%vel(:) = co_vel(:,i)[owner]
  la%o%acel(:) = co_acel(:,i)[owner]
  la%o%force(:) = co_force(:,i)[owner]
  la%o%skip = co_skip(i)[owner]
enddo

sync all

call refresh_domain_ownership()
call prof_stop(prof_sync_pos_tick, prof_sync_pos, prof_n_sync_pos)

end subroutine sync_owned_positions


subroutine mark_visible_atoms()
type(atom_dclist), pointer :: la
integer                   :: i

if(.not.allocated(atom_visible)) return

if(this_image() == 1) then
  atom_visible(:) = .true.
  return
endif

la => sys%alist
do i=1,sys%nat
  la => la%next
  atom_visible(i) = atom_owner(i) == this_image() .or. candidate_atom(la%o)
enddo

end subroutine mark_visible_atoms


subroutine refresh_domain_ownership()
! Recompute Cartesian ownership from current replicated positions and rebuild
! the local focus group consumed by the neighbor view.
type(atom_dclist), pointer :: la
integer, allocatable      :: source_owner(:)
integer                   :: i
integer                   :: new_owner, published_owner
real(dp)                  :: published_pos(3)
logical                   :: rebuild_focus

if(.not.allocated(atom_owner) .or. .not.allocated(co_owner)) return

allocate(source_owner(sys%nat))
source_owner(:) = atom_owner(:)

rebuild_focus = focus_dirty .or. num_images() > 1

if(num_images() == 1) then
  la => sys%alist
  do i=1,sys%nat
    la => la%next
    new_owner = owner_image_for_position(la%o%pos)
    if(atom_owner(i) /= new_owner) rebuild_focus = .true.
    atom_owner(i) = new_owner
  enddo
else
  if(this_image() == 1) then
    la => sys%alist
    do i=1,sys%nat
      la => la%next
      co_owner(i) = owner_image_for_position(la%o%pos)
    enddo
  endif

  sync all

  do i=1,sys%nat
    new_owner = co_owner(i)[1]
    if(atom_owner(i) /= new_owner) rebuild_focus = .true.
    atom_owner(i) = new_owner
  enddo
endif

if(.not.rebuild_focus) then
  deallocate(source_owner)
  return
endif

call hs_focus%detach_all()
call hs_candidate%detach_all()

if(allocated(atom_visible)) atom_visible(:) = .false.
la => sys%alist
do i=1,sys%nat
  la => la%next
  published_owner = source_owner(i)
  published_pos = co_pos(:,i)[published_owner]
  if(candidate_position(published_pos)) then
    if(allocated(atom_visible)) atom_visible(i) = .true.
    if(la%o%z /= co_element(i)[published_owner]) call la%o%setz(co_element(i)[published_owner])
    la%o%pos(:) = co_pos(:,i)[published_owner]
    la%o%old_cg(:) = co_old_cg(:,i)[published_owner]
    la%o%pos_old(:) = co_pos_old(:,i)[published_owner]
    la%o%vel(:) = co_vel(:,i)[published_owner]
    la%o%acel(:) = co_acel(:,i)[published_owner]
    la%o%force(:) = co_force(:,i)[published_owner]
    la%o%skip = co_skip(i)[published_owner]
    call hs_candidate%attach(la%o)
  endif
  if(this_image() == 1) then
    if(allocated(atom_visible)) atom_visible(i) = .true.
    if(la%o%z /= co_element(i)[published_owner]) call la%o%setz(co_element(i)[published_owner])
    la%o%pos(:) = co_pos(:,i)[published_owner]
    la%o%old_cg(:) = co_old_cg(:,i)[published_owner]
    la%o%pos_old(:) = co_pos_old(:,i)[published_owner]
    la%o%vel(:) = co_vel(:,i)[published_owner]
    la%o%acel(:) = co_acel(:,i)[published_owner]
    la%o%force(:) = co_force(:,i)[published_owner]
    la%o%skip = co_skip(i)[published_owner]
  endif
  if(la%o%sym /= 'CG' .and. atom_owner(i) == this_image()) call hs_focus%attach(la%o)
enddo

hs%listed = .false.
focus_dirty = .false.

deallocate(source_owner)
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
  a%acel(:)=0._dp
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

if(allocated(atom_owner) .and. allocated(co_owner) .and. allocated(co_element) .and. &
  allocated(atom_visible) .and. allocated(co_skip) .and. allocated(co_pos) .and. &
  allocated(co_old_cg) .and. allocated(co_pos_old) .and. allocated(co_vel) .and. &
  allocated(co_acel) .and. allocated(co_force)) then
  if(size(atom_owner) == sys%nat .and. size(co_owner) == sys%nat .and. &
    size(atom_visible) == sys%nat .and. size(co_pos,2) == sys%nat) return
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
if(allocated(co_owner)) deallocate(co_owner)
if(allocated(co_element)) deallocate(co_element)
if(allocated(co_skip)) deallocate(co_skip)
if(allocated(co_pos)) deallocate(co_pos)
if(allocated(co_old_cg)) deallocate(co_old_cg)
if(allocated(co_pos_old)) deallocate(co_pos_old)
if(allocated(co_vel)) deallocate(co_vel)
if(allocated(co_acel)) deallocate(co_acel)
if(allocated(co_force)) deallocate(co_force)

allocate(atom_owner(sys%nat))
allocate(co_owner(sys%nat)[*])
allocate(co_element(sys%nat)[*])
allocate(atom_visible(sys%nat))
allocate(co_skip(sys%nat)[*])
allocate(co_pos(3,sys%nat)[*])
allocate(co_old_cg(3,sys%nat)[*],co_pos_old(3,sys%nat)[*])
allocate(co_vel(3,sys%nat)[*],co_acel(3,sys%nat)[*],co_force(3,sys%nat)[*])
co_pos(:,:) = 0._dp
co_old_cg(:,:) = 0._dp
co_pos_old(:,:) = 0._dp
co_vel(:,:) = 0._dp
co_acel(:,:) = 0._dp
co_force(:,:) = 0._dp
co_owner(:) = 1
co_element(:) = 119
co_skip(:) = .false.
atom_visible(:) = .true.

atom_owner(:) = 1
if(nold > 0) atom_owner(:min(nold,sys%nat)) = old_owner(:min(nold,sys%nat))

if(sys%nat > nold) then
  la => sys%alist
  do i=1,sys%nat
    la => la%next
    if(i > nold) atom_owner(i) = owner_image_for_position(la%o%pos)
  enddo
endif

if(allocated(old_owner)) deallocate(old_owner)
sync all

end subroutine ensure_sync_storage


! Reservorios
subroutine calc_rho(rho) 
! Densidad/concentrac. reservorio
use gems_program_types, only: box
! Dos opciones: pistón,
! o con sensor: "debajo" está el sistema, y "encima" está un volumen extra de átomos.
real(dp)::vol,min_vol,z
real(dp),intent(out)::rho
type(atom_dclist), pointer :: la 
integer::i,g

g=0
min_vol=box(1)*box(2)*2.5_dp

if(s_chunk) then
  z=z1
else
  z=zmax
endif

la => sys%alist
do i=1, sys%nat ! Para contar las partícs. por encima de z0
la => la%next
 pa=> la%o
 if (pa%pos(3)>z0 .and. pa%pos(3)<z) g= g+1 ! piston
enddo

vol=box(1)*box(2)*(z-z0) ! piston
rho=g/vol
end subroutine calc_rho


subroutine config_chunk()
! Leer chunk de atoms
integer :: j
character(10)  :: sym

open(11,File='chunk.xyz')
read(11,*) nchunk
read(11,*)

do j=1,nchunk
  allocate(pb)
  call pb%init()

  ! Asigna propiedades al bloque de partículas
  read(11,*) sym,pb%pos(:)
  call pb%setsym(sym)
  pb%force(:)=0._dp
  pb%pos_old(:)=pb%pos(:)

  ! Como sus posic. en z empiezan en z= 0, las subo en zmax+rcut para cuando
  ! haya que agregarlas
  pb%pos(3)=pb%pos(3) + zmax 
  pb%pos_old(3)=pb%pos(3) + zmax

  ! Agrego un átomo a un grupo chunk
  call chunk%attach(pb)

  ! Set pbc
  pb%pbc(:)=.true.
  pb%pbc(3)=.false.

  pb=>null()
enddo
close (11)

end subroutine config_chunk


subroutine gcmc_run(g)
! Reservoir 3 - GCMC
use gems_groups, only: atom_asign
use gems_program_types, only: box, distance 
use gems_constants, only: kB_ui
use gems_errors, only: werr
class(group)               :: g
real(dp)                   :: r(3), vd(3), dr, v, rc, temp, beta
type(atom_dclist), pointer :: la
type(atom),pointer         :: o, ref
class(group),pointer       :: gp
integer                    :: i,j,n,m
 
rc=hs%rcut
           
! Compute volume
v=box(1)*box(2)*(zmax-z0)
     
! Count particles in the control volume
n=0
la=>g%alist
do j=1,g%nat
  la=>la%next
  if(la%o%pos(3)<z0.or.la%o%pos(3)>zmax) cycle
  n=n+1
enddo
    
! Point to an atom that will work as template
! in order to add new atoms into groups.
       
! Attempted adjustments 
adj: do i=1,nadj
  ref => g%alist%next%o
  beta = sqrt(kB_ui*Tsist/ref%mass)
  call werr('No more particles',.not.associated(ref))

  ! Creation attempt
 if (ran(idum)<0.5) then
              
    ! Metropolis acceptance
    if(act*v/(n+1)<ran(idum)) cycle
       
    ! Random coordinates
    r(1)=ran(idum)*box(1)
    r(2)=ran(idum)*box(2)
    r(3)=ran(idum)*(zmax-z0)+z0
                   
    ! Check overlap
    la=>g%alist
    do j=1,g%nat
      la=>la%next
      o => la%o

      ! ! Skip particles outside the control volume.
      ! FIXME: consider PBC
      ! if(o%pos(3)<z0-rc) cycle
      ! if(o%pos(3)>zmax+rc) cycle

      ! Skip if overlapping
      vd(:) = distance(o%pos,r,o%pbc)
      dr = dot_product(vd,vd)
      if(dr<rc*rc) cycle adj

    enddo

    ! Add particle
    n=n+1

    ! Initialize particle from template.
    allocate(o)
    call o%init()
    call atom_asign(o,ref)
    o%pos(:)=r(:)
    o%pos_old(:)=r(:)
     
    ! Give a velocity from maxwell-boltzman distribution
    do j = 1,3
      la%o%vel(j) = beta*gasdev()
    enddo

    ! Add to the same groups of the template.
    do j=1,ref%ngr
      gp => ref%gro(j)
      call gp%attach(o)
    enddo
        
    ! Free pointer
    o=>null()
         
  ! Destruction attempt
  else 
            
    ! Metropolis acceptance
    if(n/(v*act)<ran(idum)) cycle
                
    ! Choose a particle
    m=floor(ran(idum)*n)+1
    if(m>n) m=n

    la=>g%alist
    do j=1,g%nat
      la=>la%next
      o => la%o

      ! Skip particles outside the control volume.
      if(o%pos(3)<z0) cycle
      if(o%pos(3)>zmax) cycle

      m=m-1  
      if(m==0) exit
    enddo
    call werr('Chosen particle does not exists',m>0)
            
    ! Remove particle
    n=n-1
    call o%dest()
    deallocate(o)
       
  endif 
   
enddo adj
      
           
end subroutine gcmc_run

 
subroutine bloques(g1, g2, nx, rho) 
! Reservoir 2 - blocks of atoms added to keep bulk density
use gems_program_types, only: tbox, box_setvars  
use gems_groups, only: group, atom, atom_asign, atom_dclist
class(group)    :: g1, g2
type(atom_dclist), pointer :: la
type(atom),pointer    :: o1,o2
real(dp),intent(in)   ::rho
real(dp)              :: drho
integer               :: j
integer, intent(in)   :: nx

! Criterio para actualizar. En base a 4*sigma de desviac. estándar
drho= rho - rhomedia
if(abs(drho)<(rhomedia*0.186)) return ! volver a 0.25 ¿?

z0 = z0 + dist
z1 = z1 + dist
zmax = zmax + dist

! Deactivate neighboor list to speed up multiple addition
hs%listed=.false.

! copiamos los atomos que estaban en chunk a sys, con allocate
la=>g2%alist ! apunta al chunk
do j=1,g2%nat
  la=>la%next
  o1=> la%o
  allocate(o2)
  call o2%init()
  call g1%attach(o2) ! Sys

  ! Lista de vecinos
  call hs_candidate%attach(o2)
  call hs_focus%attach(o2)
  call hs%attach(o2)

  ! Si vos cambias el rcut, esto te va a hacer acordar de crear un nuevo chunk:
  ! TODO: call wwan("Posibilidad de colision",o1%pos(3)<hs%rcut)

  ! Nuevas partícs. reciben props. de otras ya existentes 
  call atom_asign(o2, o1)
  o2%pbc(:) = o1%pbc(:)
  o2=>null()

  ! y actualiza altura del chunk
  o1%pos(3)=o1%pos(3)+ dist 
  o1%pos_old(3)=o1%pos_old(3)+ dist
enddo

n= n + nx 

! Agrega nuevos vecinos para los atomos agregados y los cercanos
tbox(3,3)=zmax 
call box_setvars()
call sync_owned_positions()
call sync_execution_views()
           
end subroutine bloques


subroutine maxz(zmax) 
! Reservoir 1 - piston
! Ajusta el "máximo absoluto" en z de la caja de simulación con el mov. de partícs.
! Hoang et al.
real(dp)::lohi !Like a valley/bird in the sky
real(dp),intent(inout)::zmax
integer::i

!El factor para corregir zmax y las pos(3) de las que estén sobre z0
lohi= ((h/tau)*((rho0-rho)/rho)) 

do i=1, sys%nat
 pa=>sys%a(i)%o
   if (pa%pos(3)>z0) pa%pos(3)=pa%pos(3)-lohi*(pa%pos(3)-z0)
enddo

zmax=zmax-lohi*(zmax-z0)

end subroutine maxz


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
    if(s_piston) then
      if(all(o2%pos(:)==o2%old_cg(:))) then
        if(all(o1%pos(:)==o1%old_cg(:))) then
          choques3=choques3+1
          cycle
        endif
      endif
    endif
    o2%pos(:)=o2%old_cg(:) 
    ! TODO: resort?
    o2%acel(:)=0._dp
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


! Integrac. Langevin
subroutine set_ermak(h,gama,Tsist,cc0, cc1, cc2, sdr, sdv, crv1, crv2) 
! Constantes p/ Ermak
real(dp),intent(in)  :: h,gama,Tsist
real(dp),intent(out) :: cc0, cc1, cc2, sdr, sdv, crv1, crv2

!En libro: xi=kB*T/m*D=gama
! Para cambio en coef. difusión SC-SEI

!Calcula las ctes.
cc0= exp(-h*gama)
cc1= (1._dp-cc0)/gama
cc2= (1._dp-cc1/h)/gama

!Desv. estándar
sdr= sqrt(h/gama*(2._dp-(3._dp-4._dp*cc0+cc0*cc0)/(h*gama)))
sdv= sqrt(1._dp-cc0*cc0)

!Acá calcula el coef. de correlac. posic.-vel.
crv1= (1._dp-cc0)*(1._dp-cc0)/(gama*sdr*sdv)
crv2= sqrt(1._dp-(crv1*crv1))

!Un factor útil :P/para la rutina que sigue
skt=sqrt(kB_ui*Tsist)

end subroutine set_ermak


subroutine ermak_a(g,ranv) 
!Actualiza posic., din. Langevin
! Algoritmo sacado del libro de "Computer simulation of liquids" de Allen Chap 9, "Brownian Dnamics", Pag 263. Ed. vieja
class(ngroup)    :: g
class(atom), pointer        :: o1 
type(atom_dclist), pointer :: la 
real(dp),intent(out)       :: ranv(:,:)
real(dp)                   :: r1 ,r2, ranr
integer                    :: i,j
logical                    :: depos

la => g%focus%alist
do i = 1,g%focus%nat
  la => la%next
  o1 => la%o 

  !Para luego ver congelam.
  o1%old_cg(:)=o1%pos(:)
 
  ! Para cambio SC-SEI
  ! if (o1%pos(3)>80._dp) then
  !    cc1= cc1_sc
  !    cc2= cc2_sc
  !    crv1= crv1_sc
  !    crv2= crv2_sc
  ! else
  !    cc1= cc1_sei
  !    cc2= cc2_sei
  !    crv1= crv1_sei
  !    crv2= crv2_sei
  ! endif

  do j = 1, 3
    r1=gasdev()

    ranr = skt/sqrt(o1%mass)*sdr*r1
    o1%pos(j) = o1%pos(j) + cc1*o1%vel(j) + cc2*h*o1%acel(j) + ranr
  
    ! Me guardo un nro random para la veloc manteniendo la correlac con la posición.                      
    r2=gasdev()
    ranv(i,j) = skt/sqrt(o1%mass)*sdv*(crv1*r1+crv2*r2)
  end do

  call atom_pbc(o1, depos) ! ver de seleccionar forma de CG sobre electrodo

  ! Para no calcular choques o deposic. sobre Li
  ! si depositó sobre electrodo
  if (depos) cycle
  ! probar, luego ver sino si usar knock
 
  ! Mark atom to check for colisions
  o1%skip=.false. 
enddo

end subroutine ermak_a


subroutine ermak_b(g,ranv)
! calcula veloc 
class(ngroup)              :: g
type(atom), pointer        :: o1 
type(atom_dclist), pointer :: la 
real(dp),intent(in)        :: ranv(:,:)
integer                    :: i

la => g%focus%alist

do i = 1,g%focus%nat
  la => la%next
  o1 => la%o 
 
  if(o1%sym=='CG') cycle !Así se ahorra un cálculo

  o1%vel(:) = cc0*o1%vel(:) + (cc1-cc2)*o1%acel(:) + cc2*o1%force(:)/o1%mass + ranv(i,:)
  o1%acel(:)=o1%force(:)/o1%mass

enddo

end subroutine ermak_b


subroutine fuerza(g,eps,r0) 
! Fuerzas - potencial LJ
use gems_program_types, only:box  
class(ngroup)    :: g
type(atom), pointer        :: o1, o2 
type(atom_dclist), pointer :: la!, lb 
real(dp)            :: vd(3),dr,aux,b,c
real(dp),intent(in) :: eps(3,3),r0(3,3)
integer         :: i,j,ii,jj,l,k,m

la => g%focus%alist

do i=1, g%focus%nat !n
  la => la%next
  o1 => la%o
  o1%force(:) = 0._dp
  o1%epot = 0._dp  

enddo

! Calcula interacc. (F) con las partícs. vecinas
la => g%focus%alist
do ii = 1, g%focus%nat
  la => la%next 
  o1 => la%o
  i = o1%gid(g)

  k=o1%z

  do jj = 1, g%nn(i) !sobre lo vecinos

    j = g%list(i,jj)
    o2 => g%a(j)%o

    ! Skip atoms in limbo
    if(g%b_limbo) then
      if(associated(o2,target=g%limbo)) cycle
    endif
         
    m=o2%z   ! para luego poder elegir los valores de eps y r0
    vd(:) = o1%pos(:)-o2%pos(:)
   
    ! PBC en x e y
    do l=1,2       !Sin contar en z ;)
      if (vd(l)>box(l)*.5_dp) then
        vd(l)=vd(l)-box(l)
      else if (vd(l)<-box(l)*.5_dp) then
        vd(l)=vd(l)+box(l)
      endif
      
      ! if (abs(vd(l))>box*.5_dp)  vd(l)=vd(l)-sign(vd(l))*box
    enddo

    if(o2%sym=='CG'.and.o1%sym=='CG') cycle ! REVISAR

    dr = dot_product(vd,vd)

    if(dr>r0(k,m)**2) cycle

    dr=sqrt(dr)

    !factores para f
    b=r0(k,m)**6
    c=eps(k,m)*12._dp*b 
    c=c/(dr**7)
    b=b/(dr**6)

    !derivado el pot de LJ
    aux=c*(b-1)    

    o1%force(:)=o1%force(:)+aux*vd(:)/dr
    o2%force(:)=o2%force(:)-aux*vd(:)/dr

    aux = eps(k,m)*b*(b-2)

    !el pot de LJ+epsilon
    aux=aux+eps(k,m)  

    o1%epot = o1%epot + aux*.5_dp
    o2%epot = o2%epot + aux*.5_dp

  enddo
enddo

end subroutine fuerza


!Output
subroutine salida()  ! Escribe los datos calc. :P
integer  :: j, image, total_try, total_depo
type(atom_dclist), pointer :: la 
real(dp) :: energia

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

energia= 0._dp

! Coords. de partíc.
write(11,*) sys%nat ! n
write(11,*) "info:",zmax,sys%nat
la=>sys%alist
do j=1,sys%nat
 la=>la%next 
 pa=>la%o
 write(11,*) pa%sym,pa%pos(:),pa%z

 energia= energia + pa%epot 
enddo

! t, Epot+Ecin 
write(12,*) t, energia !sum(sys%a(:)%o%epot)

call kion(sys,temp) 
write(13,*)t,temp

write(14,*)t,rho

write(15,*)t,total_try
write(16,*)t,total_depo

flush(16)
flush(15)
flush(14)
flush(13)
flush(12)
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
       o1%pos_old(j)=o1%pos_old(j)-box(j)
    endif  
    if (o1%pos(j)<o) then
       o1%pos(j)=o1%pos(j)+box(j)
       o1%pos_old(j)=o1%pos_old(j)+box(j)
    endif

  else 
    ! Rebote en zmax
    if(o1%pos(j)>zmax) then
      if (integrador) then
        ! Con Ermak
        o1%pos(3) = o1%pos(3) - 2*(o1%pos(3) - zmax)
        o1%vel(3) = -o1%vel(3)

      else
        ! Con esferas duras
        o1%pos(:)= o1%old_cg(:)
      endif

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


subroutine knock(list) 
! Deposic. sobre Li ya depositado, sino, rebote brusco
use gems_program_types, only: box
type(ngroup)              :: list
real(dp)                  :: ne,vd(3),dr
integer                   :: i,ii,j,jj,l 
type(atom),pointer        :: o1,o2
type(atom_dclist),pointer :: la
 
! Por todos los pares de particulas

la => list%focus%alist
do ii = 1,list%focus%nat
  la => la%next
  o1 => la%o ! o1 es el único que puede ser CG

  i = o1%gid(list)
  ! tp = o1%gid(sys)

  do jj = 1, list%nn(i)  ! sobre los vecinos

    j = list%list(i,jj)
    o2 => list%a(j)%o

    vd(:) = o1%pos(:)-o2%pos(:)

    !Condicion de imagen minima
    do l=1,2       !Sin contar en z ;)
       if (vd(l)>box(l)*.5_dp) then
         vd(l)=vd(l)-box(l)
       else if (vd(l)<-box(l)*.5_dp) then
         vd(l)=vd(l)+box(l)
       endif
    enddo

    dr = dot_product(vd,vd)

    if(dr>list%rcut2) cycle !Sí es necesario :B

    ne=ran(idum) ! nro aleatorio para decidir si congelar o no.
    if(ne<prob) then
       call o2%setz(3)

       ! Terminac. brusca del programa si la dendrita toca el z0
       if (o2%pos(3)>z0) stop 
       exit

    else

      ! FIXME
      ! Retorno la partícula a la solución-esto es una falla u.u
      ! Podrían superponerse partícs. al retroceder...
      ! Igual que Mayers
      o2%pos(:)= o2%old_cg(:) 

    endif

  enddo

enddo


! Add new CG to the ref group
la=> list%candidate%alist
do i = 1, list%candidate%nat
  la=> la%next
  o1=> la%o
  if (o1%sym/='F') cycle
  call o1%setz(2)
  call list%focus%attach(o1)
  call list%attach(o1)
enddo

! ESTO ES DE LA RUTINA VIEJA - REVISAR POR LAS DUDAS ANTES DE TIRAR

! Lo que sigue iba por si la posic. vieja del Li+ comparada con la del CG es gde.,
! pasa a decidir si congela.
! Si están "cerca" es probable que en un paso de sim.
! anterior ya haya intentado depositar, y no vuelve a probar
! if(m==1) then
!    vd(:)=o2%pos_old(:)-o1%pos(:)  
! else
!    vd(:)=o1%pos_old(:)-o2%pos(:)  
! endif

! vd(:)=a(lit)%pos_old(:)-a(cng)%pos(:)  
end subroutine knock

 
subroutine kion(g,temp)
! Calcula T
class(group)    :: g
type(atom), pointer        :: o1 
type(atom_dclist), pointer :: la
real(dp),intent(out)::temp
real(dp)::vd,vdn,vdac !módulo de la vel., y donde acumulo
integer::i,j

la => g%alist

vdac=0. !inicia la cuenta
j=0

do i=1,g%nat

la => la%next
o1 => la%o

if(o1%sym=='CG') cycle !No considera Li metálico
j=j+1 !Cuenta los iones en mov.

vd=dot_product(o1%vel(:),o1%vel(:))
vd=vd*o1%mass !vel. al cuadrado ;) *porq. |v|=sqrt vd...


vdn=vdac+vd !acumula m*vel**2
vdac=vdn !acá como que guardo en vdac el valor de vdn para la próx.

!vdac=vdac+vd !empiezo a acumular m*vel**2
enddo

temp=vdac/(j*3._dp*kB_ui)

end subroutine kion


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
