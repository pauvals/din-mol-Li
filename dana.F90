program din_mol_Li
  use gems_groups
  use gems_neighbor
  use gems_errors, only: timer_dump, timer_start, sclock_t1, sclock_t2, wstd, sclock_max, sclock_rate,logunit
  use gems_constants, only: time1

  implicit none
 
  integer,parameter     :: dp=8
 
  !Cosas del sistema  #what a quilombo
  integer               :: n, nx, nchunk        ! Nros de particulas
  real(dp), parameter   :: o=0._dp, tau=0.1_dp  ! Origen de la caja, tau p/ rho
  real(dp)              :: z0, z1, zmax         ! Para ajuste del reservorio
  real(dp), parameter   :: gama=1._dp           ! Para fricción (no usado)
  real(dp)              :: dif, dif_sei, dif_sc ! Difusion (USADO)
  real(dp), parameter   :: Tsist=300._dp        ! Temp. del sist., 300 K
  real(dp)              :: eps(3,3),r0(3,3)     ! eps y r0 definidas como matrices para c/ tipo

  ! Constantes Físicas y
  ! factores de conversión de unidades
  real(dp), parameter :: kB_eVK=8.617330350e-5_dp  ! cte. Boltzmann en eV
  real(dp), parameter :: kJm_ui=100._dp, eV_kJm=96.485_dp,eV_ui=eV_kJm*kJm_ui
  real(dp), parameter :: kB_ui=kB_eVK*eV_ui
  real(dp), parameter :: kT_ui=kB_ui*300._dp  !en ui, la T= 300 K
 
  ! Variables dinámicas, agrupadas en átomo
  type(atom), allocatable, target :: a(:) ! átomo
  type(atom), pointer             :: pa=>null(), pb=>null()
  real(dp)            :: prob,max_vel=0._dp, msd_u= 0._dp, msd_t= 0._dp, msd_max= 0._dp

  ! Parametros de integración
  real(dp)            :: h          ! paso de tiempo
  integer             :: nst        ! nro de pasos
  integer             :: nwr        ! paso de escritura
  real(dp)            :: t=0.0_dp   ! tiempo en ps

  ! Modelo (?) 
  logical             :: integrador, reserva    ! algoritmo y estilo de reservorio usados 

  !Observables
  real(dp)            :: temp,rho,rho0

  ! Varios
  integer             :: idum         ! Semilla
  integer             :: i,j,k        ! Enteros
  real(dp)            :: vm(3)        ! Vector 3D auxiliar
  real(dp)            :: dist, rhomedia, cstdev, factor, factor2

  ! Esto es para Ermak
  logical             :: ermak 
  real(dp)            :: cc0,cc1,cc2, cc0_sei,cc1_sei,cc2_sei,&
                         cc0_sc,cc1_sc,cc2_sc
  real(dp)            :: sdr,sdv, sdr_sei,sdv_sei, sdr_sc,sdv_sc, skt
  real(dp)            :: crv1,crv2, crv1_sei,crv2_sei, crv1_sc,crv2_sc
  real(dp),allocatable:: ranv(:,:)

  ! Grupos varios
  type(group)         :: chunk

  ! Group used to compute neighbors list
  type(ngroup)        :: hs
   
  ! Init wall time
  call system_clock(count_rate=sclock_rate)
  call system_clock(count_max=sclock_max)
  call system_clock(sclock_t1)
               
  ! Valores de epsilon y r0 :P
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
      
  ! Init general del sistema de grupos
  call gindex%init()
  call ngindex%init()
  call sys%init()

  ! Init lista para choques de esferas duras
  call hs%init()
  call hs%setrc(3.2_dp) ! Maximo radio de corte
  nb_dcut=10._dp        ! The shell length for verlet update criteria
 
  call config()

  ! Grupo chunk
  if (reserva.eqv..false.) call chunk%init()

  open(25,File='version')
  write(25,*) PACKAGE_VERSION
  close(25)

  ! Leer semilla, probabilidad, nro. pasos
  ! y otros valores iniciales
  call entrada()

  ! Inicializa tamaño en z de caja simulac., y lee posic. inic.
  ! de partículas 
  call config_inic()

  !Search neighbors
  call update()

  if (integrador.eqv..true.) call fuerza(hs,eps,r0)

  ! Calcula rho-densidad reservorio inicial
  call calc_rho(rho)
  rho0=rho
  

  ! Leer chunk de atoms:
  if (reserva.eqv..false.) call config_chunk()


  ! Abro archivos de salida
  open(11,File='Li.xyz')
  open(12,File='E.dat')
  open(13,File='T.dat')
  open(14,File='rho.dat')

  call salida() ! Escribe la config. inic. en el primer paso ;)

  call timer_start()

  do i=1,nst

    ! Da un paso #wii

    if (integrador.eqv..true.) then
      ! Para trabajar con din. Langevin
      call ermak_a(hs,ranv)
      ! TODO: ver de sacar ya knock
      ! Ve si congela o rebota
      ! call knock(hs)
      call fuerza(hs,eps,r0)
      call ermak_b(hs,ranv)

    else
      ! Para usar din. Browniana, y "esferas duras"
      call cbrownian_hs(hs,h)
    endif
 
    call test_update()
    ! call update()
    
    ! Calcula cant. de partículas en reservorio o sensor
    call calc_rho(rho)

    ! Agrega bloques de partículas
    if (reserva.eqv..false.) call bloques(sys, chunk, nx, rho)

    ! Salida
    if (mod(i,nwr)==0) then
      call salida()
      ! print*, 't', t
      ! print*, 'msd_t', msd_t
      ! print*, 'msd_max', msd_max
    endif
   
    ! Timing information
    call timer_dump(i,nup=nupd_vlist)
  
    ! Para usar reservorio-pistón:
    if (reserva.eqv..true.) call maxz(zmax)

    t=t+h
   
  enddo
 
  close(11)              
  close(12)              
  close(13)
  close(14)

  
  !FIN
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
  call wstd(); write(logunit,'("vecinos actualizados: ",i0," veces")') nupd_vlist
  call wstd(); write(logunit,'("maximo numero de vecinos en algun paso: ",i0)') nn_vlist
  call wstd(); write(logunit,'("maximo desplazamiento en algun paso: ",e10.3)') sqrt(max_vel)*h
  call wstd(); write(logunit,'("MSD maximo en x-y: ",e10.3)') msd_max
 
contains

subroutine config()
  open(15,File='movedor.ini') ! by Nohe :)
  read(15,*) integrador 
  read(15,*) reserva
  close(15)
end subroutine config 

subroutine entrada()
  open(15,File='entrada.ini')
  read(15,*) idum
  read(15,*) prob
  read(15,*) h 
  read(15,*) nst
  read(15,*) nwr
  read(15,*) dist 
  read(15,*) z0 
  read(15,*) dif_sc 
  read(15,*) dif_sei 
  close(15)
end subroutine entrada

! Leer configuración inicial
subroutine config_inic()
integer :: i
type(atom_dclist), pointer :: la

! Tamaños iniciales de "reservorio"
if (reserva.eqv..true.) then
  ! Al usar pistón
  zmax= 200._dp

else
  ! Al usar chunk con sensor
  dist= dist + hs%rcut
  z1 = z0 + dist
  zmax = z1 + dist 
endif

! Para luego actualizar reservorio
rhomedia= 5.775329e-4_dp
cstdev= 1.0754306e-4_dp
factor= rhomedia - cstdev
factor2= rhomedia + cstdev
! 18.6211140525501
! 0.000577532941264052= STATS_mean
! 0.000107543058373559= 4*STATS_stddev

open(11,File='posic_inic.xyz')
read(11,*) n
read(11,*)

do i=1,n

  ! Allocatea el atomo-como un puntero del tipo 'atom'
  allocate(pa)
  call pa%init()

  read(11,*) pa%sym,pa%pos(:),pa%m
  call pa%setz(pa%sym) ! Asigna algunos valores según tipo de átomo
  pa%force(:)=0._dp
  pa%pos_old(:)=pa%pos(:)

  ! Agrega atom al sistema (sys)
  call sys%attach(pa)

  ! Agrego atomos a los grupos 
  call hs%attach(pa)
  if (pa%sym=='CG') then
    call hs%b%attach(pa)
  else
    call hs%ref%attach(pa)
    call hs%b%attach(pa)
  endif

  ! Set pbc
  pa%pbc(:)=.true.
  pa%pbc(3)=.false.

  ! Libero puntero para siguiente allocate
  pa=>null()

end do
close(11)

allocate(ranv(n,3))

! Set the box
box(:)=100._dp
box(3)=zmax

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
if (integrador.eqv..true.) then
  ! Valores inic. de las ctes. de Ermak
  ! call set_ermak(h,gama_sc,Tsist,cc0_sc, cc1_sc, cc2_sc, sdr_sc, sdv_sc, crv1_sc, crv2_sc)
  ! call set_ermak(h,gama,Tsist,cc0_sei, cc1_sei, cc2_sei, sdr_sei, sdv_sei, crv1_sei, crv2_sei)
  call set_ermak(h,gama,Tsist,cc0, cc1, cc2, sdr, sdv, crv1, crv2)

endif

end subroutine config_inic

! Leer chunk de atoms:
subroutine config_chunk()
integer :: j

open(11,File='chunk.xyz')
read(11,*) nchunk
read(11,*)

do j=1,nchunk
  allocate(pb)
  call pb%init()

  ! Asigna propiedades al bloque de partículas
  read(11,*) pb%sym,pb%pos(:),pb%m
  call pb%setz(pb%sym)
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

subroutine salida()  ! Escribe los datos calc. :P
  integer  :: j !Siempre hay que definirlos =O siempre privados
  real(dp) :: energia
 
  energia= 0._dp

  ! Coords. de partíc.
  write(11,*) sys%nat ! n
  write(11,*) "info:",zmax,sys%nat
  do j =1, sys%nat
   pa=>sys%a(j)%o
   write(11,*) pa%sym,pa%pos(:),pa%tipo

   ! Para el otro archivo de salida
   energia= energia + pa%energy 
  enddo

  ! t, suma Epot+Ecin 
  write(12,*) t, energia !sum(sys%a(:)%o%energy)

  call kion(sys,temp) 
  write(13,*)t,temp
 
  write(14,*)t,rho
  flush(14)
  flush(13)
  flush(12)

end subroutine salida

subroutine calc_rho(rho) !Densidad/concentrac.
  ! Calculada para un volumen considerado reservorio
  ! Dos opciones: pistón,
  ! o con sensor: "debajo" está el sistema, y "encima" está 
  ! un volumen extra de átomos.

  integer::i,g
  real(dp)::vol,min_vol
  real(dp),intent(out)::rho

  g=0
  min_vol=box(1)*box(2)*2.5_dp

  do i=1, sys%nat ! Para contar las partícs. por encima de z0
   pa=>sys%a(i)%o
   if (reserva.eqv..true.) then 
     if (pa%pos(3)>z0 .and. pa%pos(3)<zmax) g= g+1 ! piston
   else
     if (pa%pos(3)>z0 .and. pa%pos(3)<z1) g= g+1   ! sensor+chunk
   endif
  enddo

  if (reserva.eqv..true.) then
    vol=box(1)*box(2)*(zmax-z0) ! piston
  else
    vol=box(1)*box(2)*(z1-z0)   ! sensor+chunk
  end if
  rho= g/vol
end subroutine calc_rho

subroutine bloques(g1, g2, nx, rho) ! Crece reservorio y agrega partículas
  class(group)    :: g1, g2
  type(atom_dclist), pointer :: la
  type(atom),pointer    :: o1,o2
  real(dp),intent(in)   ::rho
  real(dp)              :: drho
  integer               :: j
  integer, intent(in)   :: nx

  ! Criterio para actualizar. En base a 4*sigma de desviac. estándar
  drho= rho - rhomedia
  if(abs(drho)>(rhomedia*0.186)) return ! volver a 0.25 ¿?

  z0 = z0 + dist
  z1 = z1 + dist
  zmax = zmax + dist

  ! copiamos los atomos que estaban en chunk a sys, con allocate
  la=>g2%alist ! apunta al chunk
  do j=1,g2%nat
     la=>la%next
     o1=> la%o
     allocate(o2)
     call o2%init()
     call g1%attach(o2) ! Sys

     ! Lista de vecinos
     call hs%attach(o2)
     call hs%b%attach(o2)
     call hs%ref%attach(o2)

     ! TODO: call wwan("Posibilidad de colision",o1%pos(3)<hs%rcut)
     ! Si vos cambias el rcut, esto te va a hacer acordar de crear un nuevo chunk.

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
  box(3)=zmax 
  call update()
                   
end subroutine bloques


! Para trabajar con pistón

! Ajusta el "máximo absoluto" en z de la caja de simulación con el mov. de partícs.
! Hoang et al.
subroutine maxz(zmax) 
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

subroutine set_sym(a,z) !Asigna tipo a partíc.
  character(*),intent(in)     :: z
  class(atom)                 :: a

  if(i>n) then ! ¿o sys%nat?
    print *, '¡Error! partíc. no existe'
    stop
  endif

  !pa=>sys%a(i) !Lee el i (intent(in)) y lo asigna luego
  select case(z) !z=símbolo de átomo
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

end subroutine set_sym


! Integrac. browniana

subroutine cbrownian_hs(g,h)
class(ngroup)              :: g
type(atom), pointer        :: o1 
type(atom_dclist), pointer :: la 
real(dp),intent(in)        :: h
real(dp)                   :: fac1, r1, posold, vd(3)
integer                    :: i,j,ii
logical                    :: depos

la => g%ref%alist
do ii = 1,g%ref%nat
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

  ! TODO: chequear sea el algoritmo de mayers- ¿si? 24.2.22
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

  call atom_hs_choque(o1, g)
 
enddo

msd_t= msd_t/g%nat
msd_max= max(msd_max,msd_t)

! Add new CG to the ref group
la=> g%ref%alist
do i = 1, g%ref%nat
  la=> la%next
  o1=> la%o
  if (o1%sym/='F') cycle
  call o1%setz('CG')
  call g%ref%detach(o1)
enddo
               
end subroutine cbrownian_hs

subroutine atom_hs_choque(o1, g)
class(atom)     :: o1
class(atom),pointer :: o2
class(ngroup)   :: g
real(dp)          :: ne, vd(3), dr
integer         :: i, j, jj

! Choque con las demas particulas
i = o1%gid(g)
do jj = 1, g%nn(i)  ! sobre los vecinos

  j = g%list(i,jj)
  o2 => g%a(j)%o

  vd(:) = vdistance(o1,o2,.true.)
  dr = dot_product(vd,vd)

  if(dr>g%rcut2) cycle !Sí es necesario :B

  ! Deposicion por contacto con otra particula metalica
  ! 1.19 es el radio de Mayers
  if (o2%tipo==2) then
    ! if(dr> 1.4161_dp) cycle
    
    ne=ran(idum)
    if(ne<prob) then
      call o1%setz('F')   ! Es inerte en este paso de tiempo

      ! Permite deposicion en cadena pero depende del atom id.
      ! Si queremos deposicion en cadnea sería mejor programarla
      ! para que no dependa de el orden en que se ejecuta el do.
      if (o1%pos(3)>z0) then
         print*,'supero z0', o1%pos(3)
         stop ! raro 
      endif
    else
      o1%pos(:)= o1%old_cg(:) 
    endif

    exit

  endif

  ! Esto pasa si es tipo==1 y están bajo rcut
  o1%pos(:)= o1%old_cg(:) 
  exit
enddo
end subroutine atom_hs_choque

! Para PBC, e intento deposic. sobre electrodo
subroutine atom_pbc(o1, depos)
class(atom) :: o1
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
      if (integrador.eqv..true.) then
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
  ne=ran(idum)

  if(ne<prob) then
    call o1%setz('F')   ! Es inerte en este paso de tiempo
    depos= .true.

    ! Para chequeo Cottrell
    ! o1%pos(:)=[0.,0.,-1.e3] 

    !!!!!
    ! Vers. vieja (Langevin 2019)
    !
    ! call o1%setz('CG') !Le dice que se congele ;)
    ! call g%ref%detach(o1)
    ! o1%pos(3)=0._dp 
    ! cycle !Cicla el do más cercano
    !
    ! Vers. vieja (Langevin 2019)
    !!!!!


! else !Acá rechazo el congelamiento y rebota
! 
!   o1%pos(3)=o1%pos(3)+2*(1._dp-o1%pos(3))
!   o1%vel(3)=-o1%vel(3)
  endif

  o1%pos(:) = o1%old_cg(:)
  ! depos= .true. ! ¿No va dentro del anterior if?
endif
end subroutine atom_pbc

! Integrac. Langevin
! Constantes p/ Ermak
subroutine set_ermak(h,gama,Tsist,cc0, cc1, cc2, sdr, sdv, crv1, crv2) 
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

!Actualiza posic., din. Langevin

subroutine ermak_a(g,ranv) 

! Algoritmo sacado del libro de "Computer simulation of liquids" de Allen Chap 9, "Brownian Dnamics", Pag 263. Ed. vieja
class(ngroup)    :: g
type(atom), pointer        :: o1 
type(atom_dclist), pointer :: la 
real(dp),intent(out)       :: ranv(:,:)
real(dp)                   :: r1 ,r2, ranr
integer                    :: i,j
logical                    :: depos

la => g%ref%alist
do i = 1,g%ref%nat
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

    ranr = skt/sqrt(o1%m)*sdr*r1
    o1%pos(j) = o1%pos(j) + cc1*o1%vel(j) + cc2*h*o1%acel(j) + ranr
  
    ! Me guardo un nro random para la veloc manteniendo la correlac con la posición.                      
    r2=gasdev()
    ranv(i,j) = skt/sqrt(o1%m)*sdv*(crv1*r1+crv2*r2)
  end do

  call atom_pbc(o1, depos) ! ver de seleccionar forma de CG sobre electrodo

  ! Para no calcular choques o deposic. sobre Li
  ! si depositó sobre electrodo
  if (depos) cycle
  call atom_hs_choque(o1, g)
  ! probar, luego ver sino si usar knock
enddo
                                                                                            
! Add new CG to the ref group
la=> g%ref%alist
do i = 1, g%ref%nat
  la=> la%next
  o1=> la%o
  if (o1%sym/='F') cycle
  call o1%setz('CG')
  call g%ref%detach(o1)
enddo

end subroutine ermak_a


! calcula veloc 
subroutine ermak_b(g,ranv)
class(ngroup)              :: g
type(atom), pointer        :: o1 
type(atom_dclist), pointer :: la 
real(dp),intent(in)        :: ranv(:,:)
integer                    :: i

la => g%ref%alist

do i = 1,g%ref%nat
  la => la%next
  o1 => la%o 
 
  if(o1%sym=='CG') cycle !Así se ahorra un cálculo

  o1%vel(:) = cc0*o1%vel(:) + (cc1-cc2)*o1%acel(:) + cc2*o1%force(:)/o1%m + ranv(i,:)
  o1%acel(:)=o1%force(:)/o1%m

enddo

end subroutine ermak_b


subroutine kion(g,temp)
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
  vd=vd*o1%m !vel. al cuadrado ;) *porq. |v|=sqrt vd...


  vdn=vdac+vd !acumula m*vel**2
  vdac=vdn !acá como que guardo en vdac el valor de vdn para la próx.

  !vdac=vdac+vd !empiezo a acumular m*vel**2
  enddo
 
  temp=vdac/(j*3._dp*kB_ui)

end subroutine kion

! Fuerzas - potencial LJ
subroutine fuerza(g,eps,r0) 
class(ngroup)    :: g
type(atom), pointer        :: o1, o2 
type(atom_dclist), pointer :: la!, lb 
real(dp)            :: vd(3),dr,aux,b,c
real(dp),intent(in) :: eps(3,3),r0(3,3)
integer         :: i,j,ii,jj,l,k,m

la => g%ref%alist

do i=1, g%ref%nat !n
  la => la%next
  o1 => la%o
  o1%force(:) = 0._dp
  o1%energy = 0._dp  

enddo
! Calcula fuerzas con las partícs. vecinas
la => g%ref%alist
do ii = 1, g%ref%nat
  la => la%next 
  o1 => la%o
  i = o1%gid(g)

  k=o1%tipo

  do jj = 1, g%nn(i) !sobre lo vecinos

    j = g%list(i,jj)
    o2 => g%a(j)%o
    m=o2%tipo   ! para luego poder elegir los valores de eps y r0
   
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

    o1%energy = o1%energy + aux*.5_dp
    o2%energy = o2%energy + aux*.5_dp

  enddo
enddo

end subroutine fuerza

! Deposic. sobre Li ya depositado, sino, rebote brusco
subroutine knock(list) 
type(ngroup)              :: list
real(dp)                  :: ne,vd(3),dr
integer                   :: i,ii,j,jj,l 
type(atom),pointer        :: o1,o2
type(atom_dclist),pointer :: la
 
! Por todos los pares de particulas

la => list%ref%alist
do ii = 1,list%ref%nat
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
       call o2%setz('F') !call set_sym(o2,'F')

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
la=> list%b%alist
do i = 1, list%b%nat
  la=> la%next
  o1=> la%o
  if (o1%sym/='F') cycle
  call o1%setz('CG')
  call list%addref(o1)
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

! subroutine hspheres(list) ! Rebote brusco en solucion
! type(ngroup)              :: list
! real(dp)                  :: ne,vd(3),dr
! integer                   :: i,ii,k,j,jj,m,l,ts,tp 
! type(atom),pointer        :: o1,o2
! type(atom_dclist),pointer :: la
!  
! ! Por todos los pares de particulas
! 
! la => list%ref%alist
! do ii = 1,list%ref%nat
!   la => la%next
!   o1 => la%o ! o1 es el único que puede ser CG
! 
!   i = o1%gid(list)
!   k= o1%tipo
!   tp = o1%id(ii)
! 
!   do jj = 1, list%nn(i)  ! sobre los vecinos
! 
!     j = list%list(i,jj)
!     o2 => list%a(j)%o
! 
!     m = o2%tipo
!     ts = o2%id(jj)
! 
!     vd(:) = o1%pos(:)-o2%pos(:)
! 
!     !Condicion de imagen minima
!     do l=1,2       !Sin contar en z ;)
!        if (vd(l)>box(l)*.5_dp) then
!          vd(l)=vd(l)-box(l)
!        else if (vd(l)<-box(l)*.5_dp) then
!          vd(l)=vd(l)+box(l)
!        endif
!     enddo
! 
!     dr = dot_product(vd,vd)
! 
!     if(dr>list%rcut2) cycle !Sí es necesario :B
! 
!     ! Retorno la partícula a la solución
!     o2%pos(:)= o2%old_cg(:) 
! 
!   enddo
! 
! enddo
! 
! end subroutine hspheres

function gasdev() !Nro aleat.
  real(dp)                  :: rsq,v1,v2
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
    gasdev=v1*rsq
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

