program din_mol_Li
  use gems_groups
  use gems_neighbor

  implicit none
 
  integer,parameter     :: dp=8
 
  !Cosas del sistema  #what a quilombo
  integer               :: n, nx                ! Nro de particulas
  real(dp), parameter   :: o=0._dp, tau=0.1_dp  ! Largo de la caja, tau p/ rho
  real(dp)              :: z0, z1, zmax         ! Para ajuste del reservorio
  real(dp), parameter   :: gama=1._dp           ! Para fricción
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
  real(dp)            :: prob

  ! Parametros de integración
  real(dp), parameter :: h=1.e-2_dp ! paso de tiempo
  integer             :: nst        ! nro de pasos
  integer             :: nwr        ! paso de escritura
  real(dp)            :: t=0.0_dp   ! tiempo en ps

  !Observables
  real(dp)            :: temp,rho,rho0

  ! Varios
  integer             :: idum         ! Semilla
  integer             :: i,j,k        ! Enteros
  real(dp)            :: vm(3)        ! Vector 3D auxiliar
  real(dp)            :: dist

  ! Esto es para Ermak
  real(dp)            :: cc0,cc1,cc2
  real(dp)            :: sdr,sdv,skt
  real(dp)            :: crv1,crv2
  ! real(dp)            :: ranv(n,3),a(n)%acel(3) !de usar esto también hay que hacer allocatable

  ! Grupos varios
  type(group)         :: chunk

  ! Group used to compute neighbors list
  type(ngroup)        :: hs
  
  ! Valores de epsilon y r0 :P
  eps(:,2) = 0
  r0(:,2)  = 0
  eps(2,:) = 0
  r0(2,:)  = 0
  r0(2,1)  = 3.5_dp
  r0(1,2)  = 3.5_dp
  eps(1,1) = 2313.6_dp
  r0(1,1)  = 3.5_dp
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
  call hs%setrc(2.5_dp) ! Maximo radio de corte
 
  ! Grupo chunk
  call chunk%init()
 
  open(25,File='version')
  write(25,*) PACKAGE_VERSION
  close(25)

  ! Leer semilla, probabilidad, nro. pasos
  ! y otros valores iniciales
  open(15,File='entrada.ini')
  read(15,*) idum
  read(15,*) prob
  read(15,*) nst
  read(15,*) nwr
  read(15,*) dist 
  read(15,*) z0 
  close(15)

  ! Tamaños iniciales de "reservorio"
  z1 = z0 +dist
  zmax = z1+dist 

  ! Leer configuración inicial
  open(11,File='posic_inic.xyz')
  read(11,*) n
  read(11,*)

  do i=1,n

    ! Allocatea el atomo-como un puntero del tipo 'atom'
    allocate(pa)

    read(11,*) pa%sym,pa%pos(:),pa%m
    call pa%setz(pa%sym) ! set_sym(pa,pa%sym) Asigna algunos valores según tipo de átomo
    pa%force(:)=0._dp
    pa%pos_old(:)=pa%pos(:)

    ! Attach atom to the system
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


  ! Set the box
  box(:)=100._dp
  box(3)=zmax

  !Search neighbors
  call update()

  ! Calcula rho-densidad reservorio inicial
  call set_rho(rho)
  rho0=rho
                      
  ! Crear chunk de atoms:
  k=0
  do j=1,sys%nat
    pa=>sys%a(j)%o
    if (pa%pos(3) > z1 .and. pa%pos(3) < zmax) then
      allocate(pb)
      k=k+1

      ! Asigna propiedades al bloque de partículas
      call atom_asign(pb, pa)

      ! Agrego un átomo a un grupo chunk
      call chunk%attach(pb)
      pb=>null()
    endif
  enddo
  nx = k !guardo el valor de k

  ! Calculo la velocidad neta del sistema
  ! Sino se trasladaría todo el sist. en el espacio... Así trabajo c/ coords.
  ! internas ;)
  do k=1,sys%nat 
    vm(:)= sum(sys%a(k)%o%pos(:) - sys%a(k)%o%pos_old(:))/n
  enddo

  ! Sustraer la velocidad neta
  do i=1,sys%nat
    pa=>sys%a(i)%o
    pa%pos_old(:) = pa%pos_old(:)-vm(:)
  end do

  ! Valores inic. de las ctes. de Ermak
  ! call set_ermak(h,gama,Tsist)

  ! calculo vel. y acelerac. iniciales
  do i=1, sys%nat
    pa=>sys%a(i)%o
    pa%vel(:) = (pa%pos(:)-pa%pos_old(:))/h
    pa%acel(:) = pa%force(:)/pa%m
  enddo
 
  ! Abro archivos de salida
  open(11,File='Li.xyz')
  open(12,File='E.dat')
  open(13,File='T.dat')
  open(14,File='rho.dat')

  call salida() !Para que escriba la config. inic. en el primer paso ;)

  do i=1,nst

    box(3)=zmax

    ! Da un paso #wii
    ! call ermak_a(a,ranv) ! revisar declaración de ranv
    ! call fuerza(a,r0)
    ! call ermak_b(a,ranv)
    ! call cbrownian(sys,h,gama,Tsist)

    call cbrownian_hs(hs,h,gama,prob,Tsist)
 
    call test_update()
    ! call update()
    
    ! !Ve si congela o rebota
    ! call knock(list)          
    
    ! Ajuste de tamaño, y cant. de partículas en reservorio
    call set_rho(rho)
    call mv_reserva(sys, chunk, nx)

    ! Salida
    if (mod(i,nwr)==0) then
      call salida()
    endif

    t=t+h
   
  enddo
 
  close(11)              
  close(12)              
  close(13)
  close(14)


contains

  subroutine mv_reserva(g1, g2, nx) ! Muevo reservorio y agrego partículas
    class(group)    :: g1, g2 
    type(atom_dclist), pointer :: la
    type(atom),pointer        :: o1,o2
    integer :: j,l,k         ! pierde acceso a la k global
    integer, intent(in) :: nx

    ! Criterio para actualizar
    if(abs(rho0-rho)<0.2*rho0) return

    z0 = z0 + dist
    z1 = z1 + dist
    zmax = zmax + dist

    ! copiamos los atomos que estaban en chunk a sys, pero es una copia, por eso el allocate
    la=>g2%alist ! apunta al chunk
    do j=1,g2%nat
       la=>la%next
       o1=> la%o
       allocate(o2)
       call g1%attach(o2) ! Sys

       ! Lista de vecinos
       call hs%attach(o2)
       call hs%b%attach(o2)
       call hs%ref%attach(o2)

       ! Nuevas partícs. reciben props. de otras ya existentes
       o1%pos(3)=o1%pos(3)+ dist 
       o1%pos_old(3)=o1%pos_old(3)+ dist
       call atom_asign(o2, o1)
       o2=>null()
    enddo

    n= n + nx 
 
    ! Agrega nuevos vecinos para los atomos agregados y los cercanos
    call update()
                     
  endsubroutine

  subroutine salida()  !Donde escribe los datos calc. :P
    integer  :: j !Siempre hay que definirlos =O siempre privados
    real(dp) :: energia
   
    energia= 0._dp

    ! Coords. de partíc.
    write(11,*) sys%nat ! n
    write(11,*)
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

  endsubroutine


  subroutine set_rho(rho) !Densidad/concentrac.
    ! Calculada para un volumen considerado reservorio, "debajo" está el sistema, y "encima" está un volumen extra
    ! de átomos.

    integer::i,g
    real(dp)::vol,min_vol
    real(dp),intent(out)::rho

    g=0
    min_vol=box(1)*box(2)*2.5_dp

    do i=1, sys%nat !Esto debería contar las partícs. por encima de z0
     pa=>sys%a(i)%o
    if (pa%pos(3)>z0.and.pa%pos(3)<z1) then
      g=g+1
    endif
    enddo

    vol=box(1)*box(2)*(z1-z0)
    rho= g/vol
  endsubroutine

  subroutine maxz(zmax) !Ajusta caja con el mov. de partícs.
    ! Ajusta el "máximo absoluto" en z de la caja de simulación
    real(dp)::lohi !Like a valley/bird in the sky
    real(dp),intent(inout)::zmax
    integer::i

    lohi= ((h/tau)*((rho0-rho)/rho)) !El factor para corregir zmax y las a(i)%pos(3) de las que estén sobre z0
    ! print *, rho0-rho
   
    do i=1, sys%nat
     pa=>sys%a(i)%o
       if (pa%pos(3)>z0) pa%pos(3)=pa%pos(3)-lohi*(pa%pos(3)-z0)
    enddo

    zmax=zmax-lohi*(zmax-z0)

  endsubroutine


  subroutine set_sym(a,z) !Asigna tipo a partíc.
    !integer,intent(in)         :: i
    character(*),intent(in)    :: z
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


  end subroutine

! Integrac. browniana

subroutine cbrownian_hs(g,h,gama,prob,Tsist)
class(ngroup)              :: g
type(atom), pointer        :: o1, o2 
type(atom_dclist), pointer :: la 
real(dp),intent(in)        :: h,gama,Tsist,prob
real(dp)                   :: fac1, fac2, r1, posold, ne, vd(3), dr
integer                    :: i,j,ii,jj

fac1 = h/gama                     
fac2 = sqrt(2._dp*kB_ui*Tsist*fac1) 

la => g%ref%alist
do ii = 1,g%ref%nat
  la => la%next
  o1 => la%o ! o1 es el único que puede ser CG

  !Para luego ver congelam. en knock
  o1%old_cg(:)=o1%pos(:)

  do j = 1,3

    r1=gasdev()
  
    posold = o1%pos(j)
    o1%pos(j) = posold + fac1*o1%force(j)/o1%m + r1*fac2/sqrt(o1%m)

    ! Velocidad derivada de euler para atras
    o1%vel(j) = (o1%pos(j)-posold)/h

    if(j<3) then
      ! PBC en x e y
      if (o1%pos(j)>box(j)) o1%pos(j)=o1%pos(j)-box(j)
      if (o1%pos(j)<o) o1%pos(j)=o1%pos(j)+box(j)
    else 
      !Con el else veo en z;es para que en z rebote en zmax, y no atraviese capa CG
      ! Rebote en zmax
      ! if(o1%pos(j)>zmax) then
      !   o1%pos(3)=o1%pos(3)-2*(o1%pos(3)-zmax)
      !   o1%vel(3)=-o1%vel(3)   !También cambio el signo de la componente de la vel. ;)
      ! endif
      o1%pos(:)= o1%old_cg(:) 
    endif

  enddo
 
  ! Si toca el electrodo implicito ¿se congela? (probabilidad ne)
  if(o1%pos(3)<1._dp) then !No importa si < o <=
    ne=ran(idum)
    if(ne<prob) then
      call o1%setz('F')   ! Es inerte en este paso de tiempo
      ! o1%pos(:)=[0.,0.,-1.e3] ! Chequeo Cottrell
    endif

    o1%pos(:) = o1%old_cg(:) 
    cycle
  endif


  ! Choque con las demas particulas
  i = o1%gid(g)
  do jj = 1, g%nn(i)  ! sobre los vecinos

    j = g%list(i,jj)
    o2 => g%a(j)%o
 
    vd(:) = vdistance(o1,o2,.true.)
    dr = dot_product(vd,vd)

    if(dr>g%rcut2) cycle !Sí es necesario :B
      
    ! Deposicion por contacto con otra particula metalica
    if (o2%tipo==2) then
      
      ne=ran(idum)
      if(ne<prob) then
        call o1%setz('F')   ! Es inerte en este paso de tiempo
        !call o2%setz('CG') 
        ! Permite deposicion en cadena pero depende del atom id.
        ! Si queremos deposicion en cadnea sería mejor programarla
        ! para que no dependa de el orden en que se ejecuta el do.
        if (o2%pos(3)>z0) stop 
      endif
    endif

    o1%pos(:)= o1%old_cg(:) 
    exit
  enddo

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
               
end subroutine

! Integrac. Langevin

subroutine set_ermak(h,gama,Tsist) !Constantes p/ Ermak
!En libro: xi=kB*T/m*D=gama
real(dp),intent(in)::h,gama,Tsist !Acá irían dist. gama

!Calcula las ctes.
cc0=exp(-h*gama)
cc1=(1._dp-cc0)/gama
cc2=(1._dp-cc1/h)/gama

!Desv. estándar
sdr=sqrt(h/gama*(2._dp-(3._dp-4._dp*cc0+cc0*cc0)/(h*gama)))
sdv=sqrt(1._dp-cc0*cc0)

!Acá calcula el coef. de correlac. posic.-vel.
crv1=(1._dp-cc0)*(1._dp-cc0)/(gama*sdr*sdv)
crv2=sqrt(1._dp-(crv1*crv1))

!Un factor útil :P/para la rutina que sigue
skt=sqrt(kB_ui*Tsist)

endsubroutine

! subroutine ermak_a(g,ranv) !Actualiza posic.
!
! ! Algoritmo sacado del libro de "Computer simulation of liquids" de Allen Chap 9, "Brownian Dnamics", Pag 263.
! class(group)    :: g
! type(atom_dclist), pointer :: la 
! type(atom), pointer        :: o1 
! real(dp),intent(out)          :: ranv(n,3)
! real(dp)                      :: r1 ,r2, ranr, ne
! integer                       :: i,j
!
! la => g%alist
!
! do i = 1,g%nat
!   la => la%next
!   o1 => la%o 
!
!   !Si ya está congelada, ni le calcula una nueva posic ;)
!   if (o1%sym=='CG') cycle 
!
!   !Para luego ver congelam.
!   o1%old_cg(:)=o1%pos(:)
!  
!   do j = 1, 3
!     r1=gasdev()
!     ranr = skt/sqrt(o1%m)*sdr*r1
!     o1%pos(j) = o1%pos(j) + cc1*o1%vel(j) + cc2*h*o1%acel(j) + ranr
!
!     !P/ que Li pueda congelarse, o que siga al rebote
!     if (o1%sym/='CG')then 
!       if(o1%pos(3)<1._dp) then !No importa si < o <=
!     
!          ne=ran(idum)
!          write (15,*) 'llamada a ran ermak a'
!          if(ne<prob) then
!              if(o1%sym=='Li') then
!                call o1%setz('CG') !Le dice que se congele ;) La sub-r. lee el CG-set_sym(o1,'CG')
!                o1%pos(3)=1._dp
!                call list%addref(o1) 
!                cycle !Cicla el do más cercano
!              endif
!
!          else !Acá rechazo el congelamiento y rebota
!      
!              o1%pos(3)=o1%pos(3)+2*(1._dp-o1%pos(3))
!              o1%vel(3)=-o1%vel(3)
!
!          endif
!       endif
!     endif
!
!     ! Me guardo un nro random para la veloc manteniendo la correlac con la posición.                      
!     r2=gasdev()
!     ranv(i,j) = skt/sqrt(o1%m)*sdv*(crv1*r1+crv2*r2)
!
!     !Pongo condiciones de caja
!     if(j<3) then
!       if (o1%pos(j)>box(j)) o1%pos(j)= o1%pos(j)-box(j)
!       if (o1%pos(j)<o) o1%pos(j)= o1%pos(j)+box(j)
!       
!       !Con el else veo en z;es para que en z rebote en zmax, y no atraviese capa CG
!     else  
!
!       if(o1%pos(j)>zmax) then
!         o1%pos(3) = o1%pos(3) - 2*(o1%pos(3) - zmax)
!         o1%vel(3) = -o1%vel(3)
!       endif
!
!       if(o1%pos(j)<1._dp) then
!         o1%pos(3) = o1%pos(3) + 2*(1._dp - o1%pos(3))
!         o1%vel(3) = -o1%vel(3)
!       endif
!
!     endif
!
!   enddo
!
! enddo
!                                                                                             
! end subroutine   


! calcula veloc 
! subroutine ermak_b(g,ranv)
! class(group)               :: g
! type(atom), pointer        :: o1 
! type(atom_dclist), pointer :: la 
! real(dp),intent(in)        :: ranv(n,3)
! integer                    :: i
!
! la => g%alist
!
! do i = 1,g%nat
!   la => la%next
!   o1 => la%o 
!  
!   if(o1%sym=='CG') cycle !Así se ahorra un cálculo
!
!   o1%vel(:) = cc0*o1%vel(:) + (cc1-cc2)*o1%acel(:) + cc2*o1%force(:)/o1%m + ranv(i,:)
!   o1%acel(:)=o1%force(:)/o1%m
!
! enddo
!
! end subroutine


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

endsubroutine


subroutine fuerza(g,eps,r0) ! según LJ
class(group)    :: g
type(atom), pointer        :: o1, o2 
type(atom_dclist), pointer :: la 
real(dp)            :: vd(3),dr,aux,b,c
real(dp),intent(in) :: eps(3,3),r0(3,3)
integer         :: i,j,l,k,m

la => g%alist

do i=1, g%nat !n
  la => la%next
  o1 => la%o
  o1%force(:) = 0._dp
  o1%energy = 0._dp  

enddo


do i = 1, g%nat-1 
  la => la%next 
  o1 => la%o

  k=o1%tipo

  do j = i+1, g%nat
    la => la%next
    o2 => la%o

    m=o2%tipo   !Determina esto para luego poder elegir los valores de eps y r0
   
    vd(:) = o1%pos(:)-o2%pos(:)
   
    !Armar la caja
    do l=1,2       !Sin contar en z ;)
      if (vd(l)>box(l)*.5_dp) then
        vd(l)=vd(l)-box(l)
      else if (vd(l)<-box(l)*.5_dp) then
        vd(l)=vd(l)+box(l)
      endif
      
      ! if (abs(vd(l))>box*.5_dp)  vd(l)=vd(l)-sign(vd(l))*box
    enddo

    if(o2%sym=='CG'.and.o1%sym=='CG') cycle

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

subroutine knock(list) ! Rebote brusco
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

      ! Retorno la partícula a la solución
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
endsubroutine

subroutine hspheres(list) ! Rebote brusco en solucion
type(ngroup)              :: list
real(dp)                  :: ne,vd(3),dr
integer                   :: i,ii,k,j,jj,m,l,ts,tp 
type(atom),pointer        :: o1,o2
type(atom_dclist),pointer :: la
 
! Por todos los pares de particulas

la => list%ref%alist
do ii = 1,list%ref%nat
  la => la%next
  o1 => la%o ! o1 es el único que puede ser CG

  i = o1%gid(list)
  k= o1%tipo
  tp = o1%id(ii)

  do jj = 1, list%nn(i)  ! sobre los vecinos

    j = list%list(i,jj)
    o2 => list%a(j)%o

    m = o2%tipo
    ts = o2%id(jj)

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

    ! Retorno la partícula a la solución
    o2%pos(:)= o2%old_cg(:) 

  enddo

enddo

endsubroutine

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
integer, intent(inout) :: idum
integer, parameter     :: k4b=selected_int_kind(9)
real                   :: ran
integer, parameter     :: ia=16807,im=2147483647,iq=127773,ir=2836
real, save             :: am
integer(k4b), save     :: ix=-1,iy=-1,k
if (idum <= 0 .or. iy < 0) then
  am=nearest(1.0,-1.0)/im
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

