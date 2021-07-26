program din_mol_Li
  
  implicit none
  
  integer,parameter     :: dp=8
  
  !Cosas del sistema  #what a quilombo
  integer               :: n, nx                !Nro de particulas
  real(dp), parameter   :: rmax=100._dp, o=0._dp, box=rmax-o, tau=0.1_dp !Largo de la caja, tau p/ rho
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
  type :: atom ! :o
          real(dp), dimension(3) :: r, rold, rnew, v, f, acel
          real(dp) :: m, energy
          character(2) :: sym
          integer :: tipo
  endtype
  type(atom), allocatable :: a(:) !,xa(:) ! átomo y el auxiliar
  type(atom), allocatable :: chunk(:)
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
  real(dp),parameter  :: dist=50._dp
  integer             :: PACKAGE_VERSION ! REVISAR

  ! Esto es para Ermak
  real(dp)            :: cc0,cc1,cc2
  real(dp)            :: sdr,sdv,skt
  real(dp)            :: crv1,crv2
  ! real(dp)            :: ranv(n,3),a(n)%acel(3) !de usar esto también hay que hacer allocatable

  open(25,File='version')
  write(25,*) PACKAGE_VERSION
  close(25)

  ! Leer semilla, probabilidad, nro. pasos
  open(15,File='entrada.ini')
  read(15,*) idum 
  read(15,*) prob
  read(15,*) nst
  read(15,*) nwr
  close(15)

  ! Leer configuración inicial
  open(11,File='posic_inic.xyz')
  read(11,*) n
  read(11,*)
  allocate(a(n))

  do i=1,n
    read(11,*) a(i)%sym,a(i)%r(:),a(i)%m
    call set_sym(i,a(i)%sym)  ! Asigna algunos valores según tipo de átomo 
    a(i)%f(:)=0._dp
    a(i)%rold(:)=a(i)%r(:)

  end do
  close(11)

 
  ! Valores iniciales
  z0=100._dp
  z1 = z0 +dist 
  zmax = z1+dist ! 200 A

  ! Calcula rho-densidad reservorio inicial
  call set_rho(rho) 
  rho0=rho
                       
  ! Crear chunk of atoms:
  ! Cuenta partículas en volumen a modificar
  k = 0
  do j=1,n
    if (a(j)%r(3) > z1 .and. a(j)%r(3) < zmax) then 
        k = k + 1
    endif
  enddo
  allocate(chunk(k))
  nx = size(chunk) !guardo el valor de k

  ! Asigna propiedades al bloque de partículas
  k=0
  do j=1,n
    if (a(j)%r(3) > z1 .and. a(j)%r(3) < zmax) then 
      k=k+1
      chunk(k)=a(j)
    endif
  enddo
 
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

  ! Calculo la velocidad neta del sistema
  ! Sino se trasladaría todo el sist. en el espacio... Así trabajo c/ coords.
  ! internas ;)
  do k=1,3
    vm(k)=sum(a(:)%r(k)-a(:)%rold(k))/n
  enddo

  ! Sustraer la velocidad neta
  do i=1,n
    a(i)%rold(:) = a(i)%rold(:)-vm(:)
    ! a(i)%v(:) = (a(:)%r(:)-a(:)%rold(:))/h ! calcula vel. inic.
  end do

  ! calculo vel. inic.
  ! a(:)%v(:)=(a(:)%r(:)-a(:)%rold(:))/h

  ! Valores inic. de las ctes. de Ermak
  ! call set_ermak(h,gama,Tsist) 

  ! calculo vel. y acelerac. iniciales
  do i=1,n
     a(i)%v(:) = (a(i)%r(:)-a(i)%rold(:))/h
     a(i)%acel(:) = a(i)%f(:)/a(i)%m
 enddo
  
  ! Abro archivos de salida
  open(11,File='Li.xyz')
  open(12,File='E.dat')
  open(13,File='T.dat')
  open(14,File='rho.dat')

  call salida() !Para que escriba la config. inic. en el primer paso ;)


  do i=1,nst

    ! Da un paso #wii
    ! call ermak_a(a,ranv) ! revisar declaración de ranv
    ! call fuerza(a,r0)
    ! call ermak_b(a,ranv)
    call cbrownian(a,h,gama,Tsist)
    call knock2(a)          !Ve si congela o no

    ! Ajuste de tamaño, y cant. de partículas en reservorio
    call set_rho(rho)
    call mv_reserva(a, nx)

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

  subroutine mv_reserva(a, nx) ! Muevo reservorio y agrego partículas
    integer :: j,l,k         ! pierde acceso a la k global 
    integer, intent(in) :: nx 
    type(atom),intent(inout), allocatable   :: a(:)
    type(atom), allocatable        :: xa(:)

    if(abs(rho0-rho)>0.2*rho0) then ! de esta forma ya cumple tb. con agrandar reservorio
       z0 = z0 + dist
       z1 = z1 + dist
       zmax = zmax + dist

       ! Hacer deallocate hace que pierda la info ya guardada en las variables... : haremos move alloc
       if (size(a) < n+nx) then 
         allocate (xa(n+size(chunk)))
         xa(1:n)= a(1:n)        ! copiado de los datos
         call move_alloc(from= xa, to= a)
       endif

       ! Nuevas partícs. reciben props. de otras ya existentes 
       chunk(:)%r(3)=chunk(:)%r(3)+ dist
       chunk(:)%rold(3)=chunk(:)%rold(3)+ dist
       chunk(:)%rnew(3)=chunk(:)%rnew(3)+ dist
       a(n+1:n+size(chunk))=chunk(:)

       n= n + size(chunk)
    endif

  endsubroutine

  subroutine salida()  !Donde escribe los datos calc. :P
    integer :: j !Siempre hay que definirlos =O siempre privados
    
    ! t, suma Epot+Ecin
    write(12,*) t,sum(a(:)%energy)

    ! Coords. de partíc.
    write(11,*) n
    write(11,*) 
    
    do j =1,n
     write(11,*) a(j)%sym,a(j)%r(:),a(j)%tipo
    enddo

    call kion(a,temp) 
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
    min_vol=rmax**2*2.5_dp

    do i=1,n !Esto debería contar las partícs. por encima de z0
    if (a(i)%r(3)>z0.and.a(i)%r(3)<z1) then
      g=g+1
    endif
    enddo

    vol=rmax**2*(z1-z0) 
    rho= g/vol
  endsubroutine

  subroutine maxz(zmax) !Ajusta caja con el mov. de partícs.
    ! Ajusta el "máximo absoluto" en z de la caja de simulación
    real(dp)::lohi !Like a valley/bird in the sky
    real(dp),intent(inout)::zmax
    integer::i

    lohi= ((h/tau)*((rho0-rho)/rho)) !El factor para corregir zmax y las a(i)%r(3) de las que estén sobre z0
    print *, rho0-rho
    
    do i=1,n
       if (a(i)%r(3)>z0) a(i)%r(3)=a(i)%r(3)-lohi*(a(i)%r(3)-z0)
    enddo

    zmax=zmax-lohi*(zmax-z0)

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
        if (rsq > 0._dp .and. rsq < 1._dp) exit
      end do
      rsq=sqrt(-2.0_dp*log(rsq)/rsq)
      gasdev=v1*rsq
      g=v2*rsq
      gaus_stored=.true.
    end if
  
  end function gasdev

  subroutine set_sym(i,z) !Asigna tipo a partíc.
    integer,intent(in)            :: i
    character(2),intent(in)       :: z

    if(i>n) then
      print *, '¡Error! partíc. no existe'
      stop
    endif

    select case(z) !z=a(i)%sym
    case('Li')
     a(i)%sym='Li' !Lee el i (intent(in)) y lo asigna
     a(i)%tipo=1
    case('CG')
     a(i)%sym='CG'
     a(i)%tipo=2
    case('F')
     a(i)%sym='F'
     a(i)%tipo=3
    case default

    end select


  end subroutine

! Integrac. browniana

subroutine cbrownian(a,h,gama,Tsist)
type(atom),intent(inout)        :: a(n)
real(dp),intent(in)     :: h,gama,Tsist !Acá irían dist. gama
real(dp)                :: fac1,fac2,r1,posold,g
integer                 :: i,j

fac1 = h/gama                      
fac2 = sqrt(2._dp*kB_ui*Tsist*fac1)  

do i = 1,n

  !Si ya está congelada, ni le calcula una nueva posic ;)
  if (a(i)%sym=='CG') cycle  
 
  !Para luego ver congelam. en knock2
  a(i)%rold(:)=a(i)%r(:) 

  do j = 1,3

    r1=gasdev()
   
    posold = a(i)%r(j)
    a(i)%r(j) = posold + fac1*a(i)%f(j)/a(i)%m + r1*fac2/sqrt(a(i)%m)

    ! Velocidad derivada de euler para atras
    a(i)%v(j) = (a(i)%r(j)-posold)/h

    if(j<3) then
      ! PBC en x e y
      if (a(i)%r(j)>rmax) a(i)%r(j)=a(i)%r(j)-box
      if (a(i)%r(j)<o) a(i)%r(j)=a(i)%r(j)+box
    else  
      !Con el else veo en z;es para que en z rebote en zmax, y no atraviese capa CG
      ! Rebote en zmax
      if(a(i)%r(j)>zmax) then
        a(i)%r(3)=a(i)%r(3)-2*(a(i)%r(3)-zmax)
        a(i)%v(3)=-a(i)%v(3)   !También cambio el signo de la componente de la vel. ;)
      endif
    endif

  enddo
  
  ! Si toca el electrodo implicito
  ! se congela con probabilidad g
  g=ran(idum)
  if(a(i)%r(3)<1._dp) then !No importa si < o <=
    if(g<prob) then
      call set_sym(i,'CG') !Le dice que se congele ;) La sub-r. lee el CG
      a(i)%r(3)=1._dp
      ! a(i)%r(:)=[0.,0.,-1.e3] ! Chequeo Cottrell
      cycle !Cicla el do más cercano
    else !Acá rechazo el congelamiento y rebota
      a(i)%r(3)=a(i)%r(3)+2*(1._dp-a(i)%r(3))
      a(i)%v(3)=-a(i)%v(3)   !También cambio el signo de la componente de la vel. ;)
    endif 
  endif

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

subroutine ermak_a(a,ranv) !Actualiza posic.

! Algoritmo sacado del libro de "Computer simulation of liquids" de Allen Chap 9, "Brownian Dnamics", Pag 263.
type(atom),intent(inout)        :: a(n)
real(dp),intent(out)          :: ranv(n,3)                                                               
real(dp)                      :: r1 ,r2, ranr,g                                                               
integer                       :: i,j                                                                   
                                                                                         
do i = 1,n

  if (a(i)%sym=='CG') cycle  !Si ya está congelada, ni le calcula una nueva posic ;)
  a(i)%rold(:)=a(i)%r(:) !Para luego ver congelam.
                                                                                       
  do j = 1, 3
    r1=gasdev()
    ranr = skt/sqrt(a(i)%m)*sdr*r1 
    a(i)%r(j) = a(i)%r(j) + cc1*a(i)%v(j) + cc2*h*a(i)%acel(j) + ranr

    !P/ que Li pueda congelarse, o que siga al rebote
    if (a(i)%sym/='CG')then  
      if(a(i)%r(3)<1._dp) then !No importa si < o <=
     
         g=ran(idum)
     
         if(g<prob) then
             if(a(i)%sym=='Li') then
               call set_sym(i,'CG') !Le dice que se congele ;) La sub-r. lee el CG
               a(i)%r(3)=1._dp
               cycle !Cicla el do más cercano
             endif

         else !Acá rechazo el congelamiento y rebota
      
             a(i)%r(3)=a(i)%r(3)+2*(1._dp-a(i)%r(3))
             a(i)%v(3)=-a(i)%v(3)

         endif
      endif
    endif

    ! Me guardo un nro random para la veloc manteniendo la correlac con la posición.                       
    r2=gasdev()
    ranv(i,j) = skt/sqrt(a(i)%m)*sdv*(crv1*r1+crv2*r2)

    !Pongo condiciones de caja
    if(j<3) then
      if (a(i)%r(j)>rmax) a(i)%r(j)=a(i)%r(j)-box
      if (a(i)%r(j)<o) a(i)%r(j)=a(i)%r(j)+box
    else  !Con el else veo en z;es para que en z rebote en zmax, y no atraviese capa CG

      if(a(i)%r(j)>zmax) then
        a(i)%r(3)=a(i)%r(3)-2*(a(i)%r(3)-zmax)
        a(i)%v(3)=-a(i)%v(3)   !También cambio el signo de la componente de la vel. ;)
      endif

      if(a(i)%r(j)<1._dp) then
        a(i)%r(3)=a(i)%r(3)+2*(1._dp-a(i)%r(3))
        a(i)%v(3)=-a(i)%v(3)
      endif

    endif

  enddo

enddo
                                                                                             
end subroutine    
                                                                                             
subroutine ermak_b(a,ranv)                                                                                  
type(atom),intent(inout)        :: a(n)
real(dp),intent(in)          :: ranv(n,3)
integer                       :: i

! calcula veloc                                                                           
do i = 1,n
  
  if(a(i)%sym=='CG') cycle !Así se ahorra un cálculo

  a(i)%v(:) = cc0*a(i)%v(:) + (cc1-cc2)*a(i)%acel(:) + cc2*a(i)%f(:)/a(i)%m + ranv(i,:)
  a(i)%acel(:)=a(i)%f(:)/a(i)%m

enddo

end subroutine


subroutine kion(a,temp)
  type(atom),intent(in)::a(n)
  real(dp),intent(out)::temp
  real(dp)::vd,vdn,vdac !módulo de la vel., y donde acumulo
  integer::i,j

  vdac=0. !acá inicio la cuenta (?)
  j=0

  do i=1,n
 
  if(a(i)%sym=='CG') cycle !No considera Li metálico
  j=j+1 !Cuenta los iones en mov.

  vd=dot_product(a(i)%v(:),a(i)%v(:))
  vd=vd*a(i)%m !vel. al cuadrado ;) *porq. |v|=sqrt vd...


  vdn=vdac+vd !acumula m*vel**2
  vdac=vdn !acá como que guardo en vdac el valor de vdn para la próx.

  !vdac=vdac+vd !empiezo a acumular m*vel**2
  enddo
  
  temp=vdac/(j*3._dp*kB_ui)

endsubroutine


subroutine fuerza(a,eps,r0) ! según LJ
type(atom),intent(inout) :: a(n) ! esto era intent (in), más lo de abajo. ¿Ahora lo paso 
                                                    ! a intent(inout)?
!real(dp),intent(out):: a(n)%energy,a(n)%f(3)
real(dp)            :: vd(3),dr,aux,b,c
real(dp),intent(in) :: eps(3,3),r0(3,3)
integer         :: i,j,l,k,m

do i=1, n
  a(i)%f(:) = 0._dp
enddo

a(:)%energy = 0._dp

do i = 1, n-1

  k=a(i)%tipo

  do j = i+1,n
    
    m=a(j)%tipo   !Determina esto para luego poder elegir los valores de eps y r0
    
    vd(:) = a(i)%r(:)-a(j)%r(:) 
    
    !Armar la caja 
    do l=1,2       !Sin contar en z ;)
      if (vd(l)>box*.5_dp) then
        vd(l)=vd(l)-box
      else if (vd(l)<-box*.5_dp) then
        vd(l)=vd(l)+box             
      endif
       
      ! if (abs(vd(l))>box*.5_dp)  vd(l)=vd(l)-sign(vd(l))*box
    enddo

    if(a(j)%sym=='CG'.and.a(i)%sym=='CG') cycle

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

    a(i)%f(:)=a(i)%f(:)+aux*vd(:)/dr
    a(j)%f(:)=a(j)%f(:)-aux*vd(:)/dr

    aux = eps(k,m)*b*(b-2)

    !el pot de LJ+epsilon
    aux=aux+eps(k,m)   

    a(i)%energy = a(i)%energy + aux*.5_dp
    a(j)%energy = a(j)%energy + aux*.5_dp 

  enddo
enddo

end subroutine fuerza

subroutine knock(a) !Congela o no
!real(dp),intent(in)::a(n)%rold(3)
type(atom),intent(inout)::a(n)
real(dp)::g,vd(3),dr
integer::i,k,j,m,lit,cng,l !'lit' y 'cng' para identificar a Li y a la cong.

do i = 1, n-1

  k=a(i)%tipo

  do j = i+1,n

    m=a(j)%tipo

    vd(:) = a(i)%r(:)-a(j)%r(:)

     !Pruebo armar la caja 
     do l=1,2       !Sin contar en z ;)
       if (vd(l)>box*.5_dp) then
       vd(l)=vd(l)-box
       else if (vd(l)<-box*.5_dp) then
       vd(l)=vd(l)+box             
       endif
     
     ! if (abs(vd(l))>box*.5_dp)  vd(l)=vd(l)-sign(vd(l))*box

     enddo

    dr = dot_product(vd,vd)
    if(dr>r0(k,m)**2) cycle !Sí es necesario :B

    if(k*m==2) then !Esto se cumple sólo si son 1 y 2 :P (Li y CG)

      ! Identifica los id del Li y del CG
      if(m==1) then
         lit=j
         cng=i
      else
         lit=i
         cng=j
      endif

      vd(:)=a(lit)%rold(:)-a(cng)%r(:)  !Si la posic. vieja del Li comparada con la del CG es gde., pasa a decidir si congela.

      !Pruebo armar la caja 
      do l=1,2       !Sin contar en z ;)
        if (vd(l)>box*.5_dp) then
          vd(l)=vd(l)-box
        else if (vd(l)<-box*.5_dp) then
          vd(l)=vd(l)+box             
        endif
       
      ! if (abs(vd(l))>box*.5_dp)  vd(l)=vd(l)-sign(vd(l))*box
      enddo


      dr = dot_product(vd,vd)

      if(dr>r0(k,m)**2) then
        g=ran(idum) !nro aleatorio para decidir si congelar o no.
        
        if(g<prob) then
          call set_sym(lit,'CG') 
          if (a(lit)%r(3)>z0) stop  !Dentro de la subr. sigue siendo lit ;)
          cycle
        endif
      endif
    
    endif

  enddo

enddo

endsubroutine

subroutine knock2(a) ! Rebote brusco
!real(dp),intent(in)::a(n)%rold(3)
type(atom),intent(inout)::a(n)
real(dp)::g,vd(3),dr
integer::i,k,j,m,lit,cng,l !'lit' y 'cng' para identificar a Li y a la cong.
  
! Por todos los pares de particulas
do i = 1, n-1
  k=a(i)%tipo
  do j = i+1,n
    m=a(j)%tipo

    !Esto se cumple sólo si son 1 y 2 :P (Li y CG)
    if(k*m/=2) cycle

    vd(:) = a(i)%r(:)-a(j)%r(:)

    !Condicion de imagen minima
    do l=1,2       !Sin contar en z ;)
       if (vd(l)>box*.5_dp) then
         vd(l)=vd(l)-box
       else if (vd(l)<-box*.5_dp) then
         vd(l)=vd(l)+box             
       endif
    enddo

    dr = dot_product(vd,vd)

    if(dr>r0(k,m)**2) cycle !Sí es necesario :B

    ! Identifica los id del Li y del CG
    if(m==1) then
       lit=j
       cng=i
    else
       lit=i
       cng=j
    endif

    g=ran(idum) ! nro aleatorio para decidir si congelar o no.
    if(g<prob) then
       call set_sym(lit,'CG') 

       ! Terminac. brusca del programa si la dendrita toca el z0
       if (a(lit)%r(3)>z0) stop  

       cycle

    else

      ! Retorno la partícula a la solución
      ! Igual que Mayers
      a(lit)%r(:)=a(lit)%rold(:)

    endif

  enddo

enddo

endsubroutine

  
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

