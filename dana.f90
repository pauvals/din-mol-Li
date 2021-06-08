program din_mol_Li
  
  implicit none
  
  integer,parameter     :: dp=8
  
  !Cosas del sistema  #what a quilombo
  integer,  parameter :: n=2000         !Nro de particulas
  real(dp),parameter  :: rmax=100._dp,o=0._dp,box=rmax-o,z0=100._dp,tau=0.1_dp !Largo de la caja, máx. z fijo, tau p/ rho
  real(dp)::zmax !El máx. en z que va variando
  real(dp), parameter :: gama=1._dp !Un gamma general, vale 1/ps
  real(dp)::m(n) !gama(n)  !masa y el gamma(a implementar)
  real(dp),parameter::Tsist=300._dp !Temp. del sist., 300 K
  character(2)        :: sym(n),z
  integer::tipo(n)  !Tipo de partíc., asignado a c/una
  real(dp)::eps(3,3),r0(3,3) !eps y r0 definidas como matrices para c/ tipo

  ! Factores de conversión de unidades
  real(dp), parameter :: kJm_ui=100._dp, eV_kJm=96.485_dp,eV_ui=eV_kJm*kJm_ui

  ! Constantes Fisicas
  real(dp), parameter :: kB_eVK=8.617330350e-5_dp  !en eV 
  real(dp), parameter :: kB_ui=kB_eVK*eV_ui
  real(dp), parameter :: kT_ui=kB_ui*300._dp  !en ui, la T= 300 K
  
  ! Variables dinámicas
  real(dp)            :: r(n,3)       ! Posiciones
  real(dp)            :: rold(n,3)    ! Posiciones anteriores
  real(dp)            :: rnew(n,3)    ! Posiciones nuevas
  real(dp)::v(n,3) !veloc.
  real(dp)            :: f(n,3)       ! Fuerza
  real(dp)::prob

  ! Parametros de integración
  real(dp), parameter :: h=1.e-2_dp !paso de tiempo en ps
  integer,  parameter :: nst=100000 !numero de pasos
  integer,  parameter :: nwr=100    !paso de escritura
  real(dp)            :: t=0.0_dp   !tiempo en ps

  !Observables
  real(dp)            :: energy(n)
  real(dp)            :: ecin,epot
  real(dp)::temp,rho,rho0

  ! Varios
  integer             :: idum         ! Semilla
  integer             :: i,j,k        ! Enteros 
  real(dp)            :: vm(3)        ! Vector 3D auxiliar
  
  !Esto es para Ermak
  real(dp)::ranv(n,3),acel(n,3)
  real(dp)::cc0,cc1,cc2
  real(dp)::sdr,sdv,skt
  real(dp)::crv1,crv2

  ! Leer semilla y probabilidad
  open(15,File='entrada.ini')
  read(15,*) idum 
  read(15,*) prob
  close(15)


  ! Leer configuración inicial
  open(11,File='posic_inic.xyz')
  read(11,*)
  read(11,*)
  do i=1,n
    read(11,*) sym(i),r(i,:),m(i)
    !Y esto llama a la subrr. que asigna valores de q,eps,r0...
    call set_sym(i,sym(i))   
  end do
  close(11)


  ! Valores iniciales
  f(:,:)=0._dp
  rold(:,:)=r(:,:)
  zmax=400._dp 
  !Acá hace un call a la sub-r. que calcula rho, y establece el primer valor=rho0
  call set_rho(rho) 
  rho0=rho

  !Valores de epsilon y r0 :P
  eps(:,2)=0
  r0(:,2)=0
  eps(2,:)=0
  r0(2,:)=0
  r0(2,1)=3.5_dp
  r0(1,2)=3.5_dp
  eps(1,1)=2313.6_dp
  r0(1,1)=3.2_dp
  eps(3,3)=121._dp
  r0(3,3)=3.61_dp
  eps(1,3)=529.1_dp
  eps(3,1)=eps(1,3)
  r0(1,3)=1.564_dp
  r0(3,1)=r0(1,3)

  ! Calculo la velocidad neta del sistema/Sino como que se trasladaría todo el sist. en el espacio... Así trabajo c/ coords.
  ! internas ;)
  do k=1,3
    vm(k)=sum(r(:,k)-rold(:,k))/n
  enddo

  ! Sustraer la velocidad neta
  do i=1,n
    rold(i,:)=rold(i,:)-vm(:)
  end do

  !calculo vel. inic.
  v(:,:)=(r(:,:)-rold(:,:))/h

  !Valores inic. de las ctes. de Ermak
  call set_ermak(h,gama,Tsist) 

  call fuerza(r,f,eps,r0,energy)
  do i=1,n
     acel(i,:)=f(i,:)/m(i)
  enddo

  
  ! Abro archivos de salida
  open(11,File='Li.xyz')
  open(12,File='E.dat')
  open(13,File='T.dat')
  open(14,File='rho.dat')

  call salida() !Para que escriba la config. inic. en el primer paso ;)


  do i=1,nst

    ! Da un paso #wii
    call ermak_a(r,v,acel,ranv)
    call knock(r,rold)          !Ve si congela o no
    call fuerza(r,f,eps,r0,energy)
    call ermak_b(v,acel,f,ranv)
    call set_rho(rho)

    !Salida 
    if (mod(i,nwr)==0) then

      call salida()
       
    endif

    call maxz(zmax)

    t=t+h
    
  enddo


  close(11)               
  close(12)               
  close(13)
  close(14)


contains

        subroutine salida()  !Donde escribe los datos calc. :P

                !real(dp),intent(inout)::sym(n),r(n,3),t//le aclaro el intent al arg. declarado para la sub-r
                integer::j !Siempre hay que definirlos =O
                
                ! t, suma Epot+Ecin
                write(12,*) t,sum(energy)

                ! Coords. de partíc.
                write(11,*) n
                write(11,*) 
                
                do j =1,n
                 write(11,*) sym(j),r(j,:)
                enddo

                call kion(v,temp) 
                write(13,*)t,temp
                
                write(14,*)t,rho
                flush(14)
                flush(13)
                flush(12)

        endsubroutine


        subroutine set_rho(rho) !Densidad/concentrac.

                integer::i,g
                real(dp)::vol,min_vol
                real(dp),intent(out)::rho

                g=0
                min_vol=rmax**2*2.5_dp

                do i=1,n !Esto debería contar las partícs. por encima de z0
                if (r(i,3)>z0) then
                        g=g+1
                endif
                enddo

                vol=rmax**2*(zmax-z0) !zmax es lo que iré recalculando, pa la próx. iteración
                
                if (vol<min_vol) stop !Si el V del reservorio se hace pequeño, corta todo, sin escribir este último
                                      !resultado.

                rho= g/vol


        endsubroutine

        subroutine maxz(zmax) !Ajusta caja con el mov. de partícs.
                real(dp)::lohi !Like a valley/bird in the sky
                real(dp),intent(inout)::zmax
                integer::i

                lohi= ((h/tau)*((rho0-rho)/rho)) !El factor para corregir zmax y las r(i,3) de las que estén sobre z0
                
                do i=1,n
                   if (r(i,3)>z0) r(i,3)=r(i,3)-lohi*(r(i,3)-z0)
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

                select case(z) !z=sym(i)
                case('Li')
                 sym(i)='Li' !Lee el i (intent(in)) y lo asigna
                 tipo(i)=1
                case('CG')
                 sym(i)='CG'
                 tipo(i)=2
                case('F')
                 sym(i)='F'
                 tipo(i)=3
                case default

                end select


        end subroutine


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


  subroutine ermak_a(r,v,acel,ranv) !Actualiza posic.

 ! Algoritmo sacado del libro de "Computer simulation of liquids" de Allen Chap 9, "Brownian Dnamics", Pag 263.
 real(dp),intent(inout)        :: r(n,3),v(n,3)
 real(dp),intent(in)           :: acel(n,3)                                                              
 real(dp),intent(out)          :: ranv(n,3)                                                               
 real(dp)                      :: r1 ,r2, ranr,g                                                               
 integer                       :: i,j                                                                   
                                                                                                                  
 do i = 1,n

   if (sym(i)=='CG') cycle  !Si ya está congelada, ni le calcula una nueva posic ;)
   rold(i,:)=r(i,:) !Para luego ver congelam.
                                                                                                                
   do j = 1, 3
     r1=gasdev()
     ranr = skt/sqrt(m(i))*sdr*r1 
     r(i,j) = r(i,j) + cc1*v(i,j) + cc2*h*acel(i,j) + ranr

     !P/ que Li pueda congelarse, o que siga al rebote
     if (sym(i)/='CG')then  
          if(r(i,3)<1._dp) then !No importa si < o <=
                  
                  g=ran(idum)
                  
                  if(g<prob) then
                          if(sym(i)=='Li') then
                                  call set_sym(i,'CG') !Le dice que se congele ;) La sub-r. lee el CG
                                  r(i,3)=1._dp
                                  cycle !Cicla el do más cercano
                          endif

                  else !Acá rechazo el congelamiento y rebota
                         
                          r(i,3)=r(i,3)+2*(1._dp-r(i,3))
                          v(i,3)=-v(i,3)

                  endif
          endif
     endif

     ! Me guardo un nro random para la veloc manteniendo la correlac con la posición.                       
     r2=gasdev()
     ranv(i,j) = skt/sqrt(m(i))*sdv*(crv1*r1+crv2*r2)

     !Pongo condiciones de caja
     if(j<3) then
             if (r(i,j)>rmax) r(i,j)=r(i,j)-box
             if (r(i,j)<o) r(i,j)=r(i,j)+box
     else  !Con el else veo en z;es para que en z rebote en zmax, y no atraviese capa CG

             if(r(i,j)>zmax) then
                     r(i,3)=r(i,3)-2*(r(i,3)-zmax)
                     v(i,3)=-v(i,3)   !También cambio el signo de la componente de la vel. ;)
             endif

             if(r(i,j)<1._dp) then
                     r(i,3)=r(i,3)+2*(1._dp-r(i,3))
                     v(i,3)=-v(i,3)
             endif

     endif

   enddo



 enddo
                                                                                                                      
 end subroutine    
                                                                                                                      
 subroutine ermak_b(v,acel,f,ranv)                                                                                  
 real(dp),intent(inout)        :: v(n,3),acel(n,3)
 real(dp),intent(in)           :: f(n,3)                                                             
 real(dp),intent(in)          :: ranv(n,3)                                                                  
 integer                       :: i                                                                                  
 ! calcula veloc                                                                           
 do i = 1,n
   
   if(sym(i)=='CG') cycle !Así se ahorra un cálculo

   v(i,:) = cc0*v(i,:) + (cc1-cc2)*acel(i,:) + cc2*f(i,:)/m(i) + ranv(i,:)
   acel(i,:)=f(i,:)/m(i)

  enddo

 end subroutine
                                                                                                          
 subroutine kion(v,temp)
          real(dp),intent(in)::v(n,3)
          real(dp),intent(out)::temp
          real(dp)::vd,vdn,vdac,otrak !módulo de la vel., y donde acumulo
          integer::i,j

          vdac=0. !acá inicio la cuenta (?)
          j=0

          do i=1,n
        
          if(sym(i)=='CG') cycle !No considera Li metálico
          j=j+1 !Cuenta los iones en mov.

          vd=dot_product(v(i,:),v(i,:))
          vd=vd*m(i) !vel. al cuadrado ;) *porq. |v|=sqrt vd...


          vdn=vdac+vd !acumula m*vel**2
          vdac=vdn !acá como que guardo en vdac el valor de vdn para la próx.

          !vdac=vdac+vd !empiezo a acumular m*vel**2
          enddo
          
          temp=vdac/(j*3._dp*kB_ui)

  endsubroutine


  subroutine fuerza(r,f,eps,r0,energy)
    real(dp)            :: vd(3),dr,aux,a,b,g,rcg(3)
    real(dp),intent(in) :: r(n,3),eps(3,3),r0(3,3)
    real(dp),intent(out):: energy(n),f(n,3)
    integer         :: i,j,l,k,m

    f(:,:) = 0._dp
    energy(:) = 0._dp

    do i = 1, n-1

      k=tipo(i)

      do j = i+1,n

        m=tipo(j)   !Determina esto para luego poder elegir los valores de eps y r0

        vd(:) = r(i,:)-r(j,:) 

         !Armar la caja 
         do l=1,2       !Sin contar en z ;)
            if (vd(l)>box*.5_dp) then
                  vd(l)=vd(l)-box
            else if (vd(l)<-box*.5_dp) then
                  vd(l)=vd(l)+box             
            endif
          
         ! if (abs(vd(l))>box*.5_dp)  vd(l)=vd(l)-sign(vd(l))*box

         enddo

        if(sym(j)=='CG'.and.sym(i)=='CG') cycle

        dr = dot_product(vd,vd)

        if(dr>r0(k,m)**2) cycle

        dr=sqrt(dr)

        !factores para f
        b=r0(k,m)**6
        a=eps(k,m)*12._dp*b  
        a=a/(dr**7)
        b=b/(dr**6)

        !derivado el pot de LJ
        aux=a*(b-1)     

        f(i,:)=f(i,:)+aux*vd(:)/dr
        f(j,:)=f(j,:)-aux*vd(:)/dr
 
        aux = eps(k,m)*b*(b-2)

        !el pot de LJ+epsilon
        aux=aux+eps(k,m)   

        energy(i) = energy(i) + aux*.5_dp
        energy(j) = energy(j) + aux*.5_dp 

  
      enddo
    enddo

  end subroutine fuerza

  subroutine knock(r,rold) !Congela o no

          real(dp),intent(in)::rold(n,3)
          real(dp),intent(inout)::r(n,3)
          real(dp)::g,a(3),b(3),c,d,e,vd(3),dr
          integer::i,k,j,m,lit,cng,l !'lit' y 'cng' para identificar a Li y a la cong.
          
    do i = 1, n-1

      k=tipo(i)

      do j = i+1,n

        m=tipo(j)

        vd(:) = r(i,:)-r(j,:)


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
        if(dr>r0(k,m)**2) cycle !Sí es necesario, pavota :B


        if(k*m==2) then !Esto se cumple sólo si son 1 y 2 :P (Li y CG)

                ! Identifica los id del Li y del CG
                if(m==1) then
                   lit=j
                   cng=i
                else
                   lit=i
                   cng=j
                endif

                vd(:)=rold(lit,:)-r(cng,:)  !Si la posic. vieja del Li comparada con la del CG es gde., pasa a decidir si congela.

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
                                call set_sym(lit,'CG') !Independientemente de si es 'i' o 'j' la que es Li, ya sabe que es la 'l'
                                if (r(lit,3)>z0) stop  !Dentro de la subr. sigue siendo lit ;)
                                cycle
                        endif
                endif
        
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

