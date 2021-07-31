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

 
module gems_constants

implicit none
public 

!                                                           dimension variables
!------------------------------------------------------------------------------

! el ancho de linea para el interprete de comandos
integer,parameter       :: linewidth = 200


! Elijo el numero de dimensiones a tiempo de compilacion con el precompilador
#ifdef DIM2
  integer,parameter       :: dm=2    ! La dimension de la dinamica. Maximo 3
#else
#ifdef DIM1
  integer,parameter       :: dm=1    ! La dimension de la dinamica. Maximo 3                    
#else
  integer,parameter       :: dm=3    ! La dimension de la dinamica. Maximo 3                    
#endif
#endif
 
! Esta variable es necesaria para no tener que convertir
! cada vez que se necesite dm a tipo character
character(1)  :: cdm

! La dimension de proyeccion. Debe ser igual o menor a dm
! sirve para ignorar las otras dimensiones aunque esten
! compiladas. Es una forma de restringir la dimension en tiempo
! de ejecucion
integer       :: dd=dm    



! La precision del programa
! Evitar long integers por los indices de los arrays
integer, parameter    :: isp = 4 !normal=> 8 cifras
integer, parameter    :: idp  = selected_int_kind(18) !18 cifras
integer, parameter    :: sp   = kind(1.0)
integer, parameter    :: dp   = kind(1.0d0)
integer, parameter    :: ie1  = selected_int_kind(1)
! integer, parameter    :: dp = selected_real_kind(12)

! para medición de eficiencia, tiempo real
real(dp),public    :: time1,time2     


! Numero pi
!real(dp),parameter  :: pi =3.141592653589793238462643383279502884197_dp    ,&
!                       pio2=1.57079632679489661923132169163975144209858_dp ,&
!                       twopi=6.283185307179586476925286766559005768394_dp  ,&
real(dp),parameter  :: pi = 4.0_dp*atan(1.0_dp)  ,&
                       pio2=pi*0.5_dp            ,&
                       twopi=2.0_dp*pi           ,&
                       twoopi=2.0_dp/pi          ,& 
                       oneopi=1.0_dp/pi          ,& 
                       sqrtpi=sqrt(pi)           ,& 
                       oneospi=1.0_dp/sqrtpi     ,& 
                       twoospi=2.0_dp*oneospi    

! Numero de avogadro
real(dp),parameter  :: navog = 6.0221412927e23_dp

! Unidades
! *unit1_unit2 say "change the unit1 for unit2" (this is the operator *unit2/unit1)

! Sistema de Unidades Internas :
!   longitud ---> A
!   masa ---> uma
!   tiempo ---> ps
!   Esto deriva a 
!   energia --- uma*a**2/ps**2 ---> ui
!   fuerza  en uma*a/ps**2 = ui/a

! Otro Sistema de Unidades Internas (ue):
!   longitud en A
!   masa en uma
!   tiempo en ps
!   energia en eV


! Escala
real(dp),parameter  :: ps_us = 1.0e-6_dp,                           & ![micros/ps]
                       us_ps = 1/ps_us,                             & ![ps/micros]
                       ps_ms = 1.0e-9_dp,                           & ![milis/ps]
                       ms_ps = 1/ps_ms,                             & ![ps/milis]
                       ps_s = 1.0e-12_dp,                           & ![s/ps]
                       s_ps = 1/ps_s,                               & ![ps/s]
                       axps_mxs = 1.0e2_dp,                         & ![m/s*ps/a]
                       mxs_axps = 1.0e-2_dp                         

real(dp),parameter  :: m_a = 1.0e10_dp                                ![m/a]

! Constantes fisicas en SI 
real(dp),parameter  :: qe_si=1.60219e-19_dp                           ![C] carga del electron

real(dp),parameter  :: ev_joule = qe_si,                            & ![j/ev]
                       joule_ev = 1.0_dp/ev_joule,                  & !
                       joulemol_ev = joule_ev/navog,                & !
                       ev_joulemol = 1.0_dp/joulemol_ev,            & !
                       cal_ev = 2.61144768e19,                      & !
                       ev_cal = 1.0_dp/cal_ev,                      & !
                       calmol_ev = cal_ev/navog,                    & !
                       ev_calmol = 1.0_dp/calmol_ev,                & ![kj/ev]
                       uma_kg = 1.6605402e-27_dp,                   & ![kg/amu]
                       kg_uma = 0.60221366e27_dp,                   & !
                       ui_ev = axps_mxs*axps_mxs*uma_kg*joule_ev,   & !(1,0364190264575364 d-4)[(ps/a)**2/uma*ev]
                       ev_ui = 1.0_dp/ui_ev,                        & 
                       ui_joule = ui_ev*ev_joule,                   & !1.66054e-23
                       joule_ui = 1.0_dp/ui_joule,                  & 
                       kjm_ui = 1000.0_dp*joule_ui/navog,           & !
                       ui_kjm = 1.0_dp/kjm_ui,                      & !
                       e_coulomb = qe_si,                           & ![C/e]
                       coulomb_e = 1.0_dp/e_coulomb                   ![e/C]


real(dp),parameter  :: eV_kcm = 23.06035, kcm_eV=1._dp/eV_kcm,      &
                       ui_kcm = ui_eV*eV_kcm, kcm_ui=1._dp/ui_kcm

! Unidades atomicas 
real(dp),parameter  :: bohr_ang = 0.5291772108_dp,                  & !
                       ang_bohr = 1._dp/bohr_ang,                   & !
                       hartree_ev = 27.211383411_dp,                & !
                       ev_hartree = 1.0_dp/hartree_ev,              & !
                       hartree_ui = hartree_ev*ev_ui,               & !
                       ui_hartree = 1.0_dp/hartree_ui,              & !
                       hartreebohr_eVang=hartree_ev/bohr_ang,       & !
                       eVang_hartreebohr=1._dp/hartreebohr_eVang,   & !0.0194469064593167_dp
                       hartreebohr_ui=hartreebohr_eVang*ev_ui,      & !
                       ui_hartreebohr=1._dp/hartreebohr_ui
                                

! Radianes
 real(dp),parameter  :: rad120 =twopi/3.0_dp           ,& 
                        rad90 = pio2                   ,& 
                        rad60 =pi/3.0_dp               ,& 
                        rad109 =109.0_dp*pi/180.0_dp   ,& 
                        grad_rad = pi/180.0_dp         ,& 
                        rad_grad = 180.0_dp/pi

! Trigonometricos
 real(dp),parameter  :: cos120 =cos(rad120)             ,&
                        sin120=sin(rad120)              ,&
                        cos109=cos(rad109)              ,&
                        sin109=sin(rad109)              

! Otros numeros                       
real(dp),parameter  :: sqrt2=1.41421356237309504880168872420969807856967_dp,&  
                       sqrt3=sqrt(3.0_dp),&  
                       oneos2=1.0_dp/sqrt2                                 ,&  
                       euler=0.5772156649015328606065120900824024310422_dp

! Boltzman                       
real(dp),parameter  :: kB_ev = 8.617385e-05_dp                      ,& !ev/k
                       kB_j = kB_ev*ev_joule                  ,& !j/k
                       kB_ui = kB_ev*ev_ui                       !ui/k

! Permitividad del vacio
real(dp),parameter  :: permvac_si = 8.8541878176e-12                     ,& ![F/m]=[C*C/J/m]
                       permvac_ui = coulomb_e*coulomb_e/(joule_ui*m_a)      ![1/(ui A)]

! Coulomb
real(dp),parameter  :: kcoulomb_si = 1.0_dp/(4*pi*permvac_si)            ,& ![m/F]
                       kcoulomb_ui = kcoulomb_si/permvac_ui                 ![ui A]
                   
! Presion
real(dp),parameter  :: atm_bar = 1.01325                  ,&
                       bar_atm = 1._dp/atm_bar            ,& 
                       bar_pa = 1.e5_dp                   ,&
                       pa_bar = 1._dp/bar_pa              ,& 
                       pa_ui = joule_ui*1e-30             ,&
                       ui_pa = 1.0_dp/pa_ui               ,&
                       atm_ui = atm_bar*bar_pa*pa_ui      ,&
                       ui_atm = 1.0_dp/atm_ui
                       
                                                       
! Otras....                      
real(dp),parameter  :: planck_js = 6.6260755e-34_dp                          !js 


! Variables auxiliares (no constantes a compilacion)
real(dp)       :: beta,kt,temp

contains


integer function find_io(start)

!  find an unused unit number for input or output. unit n=start is used
!  if available; otherwise n is incremented until an unused unit is found.
!  unit numbers are limited to the range 1-100; if n reaches 100 the
!  search starts again at 1.

integer, intent(in) :: start
logical :: in_use, exists
integer :: n, n0
integer, parameter :: max_unit=99

n0=start
if (n0 <= 1 .or. n0 > max_unit) n0=1
n=n0
in_use=.true.
do while (in_use)
  inquire(n,opened=in_use,exist=exists)
  if (exists) then
    if (.not. in_use) exit
  else
    !FIXME write (unit=string,fmt="(a,i3,a)") "unit number", n, " out of range"
    !call report (string)
  endif
  n=n+1
  if (n > max_unit) n=1
  if (n == n0) then
    !FIXME call report ("no i/o unit available")
  end if
end do
find_io=n

end function find_io


! function uid(a,b)
! !http://stackoverflow.com/a/13871379/1342186
! !Cantor pairing function is really one of the better ones out there
! !considering its simple, fast and space efficient, but there is
! !something even better published at Wolfram by Matthew Szudzik,
! !here. The limitation of Cantor pairing function (relatively) is
! !that the range of encoded results doesn't always stay within the
! !limits of a 2N bit integer if the inputs are two N bit integers.
! !That is, if my inputs are two 16 bit integers ranging from 0 to
! !2^16 -1, then there are 2^16 * (2^16 -1) combinations of inputs
! !possible, so by the obvious Pigeonhole Principle, we need an output
! !of size at least 2^16 * (2^16 -1), which is equal to 2^32 - 2^16,
! !or in other words, a map of 32 bit numbers should be feasible
! !ideally. This may not be of little practical importance in
! !programming world.
!   integer           :: a,b !(must be positive)
!   integer           :: uid
!   if(a >= b) then 
!     uid = a*a+a+b 
!   else
!     uid = a+b*b
!   endif
! end function
!
! function uid(a,b)
! !http://stackoverflow.com/a/13871379/1342186
!   integer,intent(in) :: a,b !(can be negative)
!   integer            :: uid
!   integer            :: a2,b2
! 
!   if(a>=0) then
!     a2=2*a
!   else
!     a2=-2*a-1
!   endif
!       
!   if(b>=0) then
!     b2=2*a
!   else
!     b2=-2*a-1
!   endif
!       
!   if(a2 >= b2) then 
!     uid=a2*a2+a2+b2 
!   else
!     uid=a2+b2*b2
!   endif
!
!   uid=uid/2
!
!   if (a<0.and.b>0) uid=-uid-1
!   if (a>0.and.b<0) uid=-uid-1
!
! end function


function same_proc (a,b)
! Associated(procedure_pointer,procedure_target) do not work, this is a
! solution given in https://stackoverflow.com/a/48977783/1342186
use, intrinsic :: iso_c_binding
logical same_proc
type(c_funptr), intent(in) :: a,b

same_proc = transfer(a,0_C_INTPTR_T) == transfer(b,0_C_INTPTR_T)
end function same_proc

end module gems_constants

