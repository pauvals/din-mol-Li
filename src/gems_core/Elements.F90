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

module gems_elements
use gems_constants, only:dp

implicit none

private
public   :: inq_z, add_z, set_z, elements_init
 
type, public :: element
   character(:),allocatable  :: sym
   real(dp)                  :: carga=0._dp
   real(dp)                  :: mass=1._dp
   
  ! Radio de corte covalente, lo uso para definir el nro de enlaces covalentes.
  ! Los valores son toqueteados para que los algoritmos que los necesiten funcionen.
  ! En definitiva, van a funcionar si este numero distingue bien entre un enlace
  ! covalente y un atomo no enlazado vecino cercano. Por eso lo hago por especie.
 
   real(dp)       :: rad=1._dp
 
end type

#define _NODE element_v
#define _TYPE type(element)
#include "vector_header.inc"
   
type(element_v),target,public   :: elements

contains

#define _NODE element_v
#define _TYPE type(element)
#include "vector_body.inc"

subroutine elements_init
type(element),pointer :: e(:)

! There are 118 elements in the periodic table
! I use one more for the generic element A
call elements%init(119)
elements%size=119
e=>elements%o(:)

e(1  )%sym='H'  ; e(1  )%mass=1.008    
e(2  )%sym='He' ; e(2  )%mass=4.0032   
e(3  )%sym='Li' ; e(3  )%mass=6.941    
e(4  )%sym='Be' ; e(4  )%mass=9.012    
e(5  )%sym='B'  ; e(5  )%mass=10.81    
e(6  )%sym='C'  ; e(6  )%mass=12.01; e(6)%rad=1.74   
e(7  )%sym='N'  ; e(7  )%mass=14.01    
e(8  )%sym='O'  ; e(8  )%mass=16.00    
e(9  )%sym='F'  ; e(9  )%mass=19.00    
e(10 )%sym='Ne' ; e(10 )%mass=20.183   
e(11 )%sym='Na' ; e(11 )%mass=22.99    
e(12 )%sym='Mg' ; e(12 )%mass=24.31    
e(13 )%sym='Al' ; e(13 )%mass=26.98    
e(14 )%sym='Si' ; e(14 )%mass=28.09    
e(15 )%sym='P'  ; e(15 )%mass=30.97    
e(16 )%sym='S'  ; e(16 )%mass=32.07    
e(17 )%sym='Cl' ; e(17 )%mass=35.45    
e(18 )%sym='Ar' ; e(18 )%mass=39.954   
e(19 )%sym='K'  ; e(19 )%mass=39.10    
e(20 )%sym='Ca' ; e(20 )%mass=40.08    
e(21 )%sym='Sc' ; e(21 )%mass=44.96    
e(22 )%sym='Ti' ; e(22 )%mass=47.87    
e(23 )%sym='V'  ; e(23 )%mass=50.94    
e(24 )%sym='Cr' ; e(24 )%mass=52.00    
e(25 )%sym='Mn' ; e(25 )%mass=54.94    
e(26 )%sym='Fe' ; e(26 )%mass=55.84    
e(27 )%sym='Co' ; e(27 )%mass=58.93    
e(28 )%sym='Ni' ; e(28 )%mass=58.69    
e(29 )%sym='Cu' ; e(29 )%mass=63.55    
e(30 )%sym='Zn' ; e(30 )%mass=65.39    
e(31 )%sym='Ga' ; e(31 )%mass=69.72    
e(32 )%sym='Ge' ; e(32 )%mass=72.61    
e(33 )%sym='As' ; e(33 )%mass=74.92    
e(34 )%sym='Se' ; e(34 )%mass=78.96    
e(35 )%sym='Br' ; e(35 )%mass=79.90    
e(36 )%sym='Kr' ; e(36 )%mass=83.805   
e(37 )%sym='Rb' ; e(37 )%mass=85.47    
e(38 )%sym='Sr' ; e(38 )%mass=87.62    
e(39 )%sym='Y'  ; e(39 )%mass=88.91    
e(40 )%sym='Zr' ; e(40 )%mass=91.22    
e(41 )%sym='Nb' ; e(41 )%mass=92.91    
e(42 )%sym='Mo' ; e(42 )%mass=95.94    
e(43 )%sym='Tc' ; e(43 )%mass=99.      
e(44 )%sym='Ru' ; e(44 )%mass=101.07   
e(45 )%sym='Rh' ; e(45 )%mass=102.91   
e(46 )%sym='Pd' ; e(46 )%mass=106.42   
e(47 )%sym='Ag' ; e(47 )%mass=107.87   
e(48 )%sym='Cd' ; e(48 )%mass=112.41   
e(49 )%sym='In' ; e(49 )%mass=114.82   
e(50 )%sym='Sn' ; e(50 )%mass=118.71   
e(51 )%sym='Sb' ; e(51 )%mass=121.76   
e(52 )%sym='Te' ; e(52 )%mass=127.60   
e(53 )%sym='I'  ; e(53 )%mass=126.90   
e(54 )%sym='Xe' ; e(54 )%mass=131.296  
e(55 )%sym='Cs' ; e(55 )%mass=132.91   
e(56 )%sym='Ba' ; e(56 )%mass=137.33   
e(57 )%sym='La' ; e(57 )%mass=0.       
e(58 )%sym='Ce' ; e(58 )%mass=0.       
e(59 )%sym='Pr' ; e(59 )%mass=0.       
e(60 )%sym='Nd' ; e(60 )%mass=0.       
e(61 )%sym='Pm' ; e(61 )%mass=0.       
e(62 )%sym='Sm' ; e(62 )%mass=0.       
e(63 )%sym='Eu' ; e(63 )%mass=0.       
e(64 )%sym='Gd' ; e(64 )%mass=0.       
e(65 )%sym='Tb' ; e(65 )%mass=0.       
e(66 )%sym='Dy' ; e(66 )%mass=0.       
e(67 )%sym='Ho' ; e(67 )%mass=0.       
e(68 )%sym='Er' ; e(68 )%mass=0.       
e(69 )%sym='Tm' ; e(69 )%mass=0.       
e(70 )%sym='Yb' ; e(70 )%mass=0.       
e(71 )%sym='Lu' ; e(71 )%mass=0.       
e(72 )%sym='Hf' ; e(72 )%mass=178.49   
e(73 )%sym='Ta' ; e(73 )%mass=180.95   
e(74 )%sym='W'  ; e(74 )%mass=183.84   
e(75 )%sym='Re' ; e(75 )%mass=186.21   
e(76 )%sym='Os' ; e(76 )%mass=190.23   
e(77 )%sym='Ir' ; e(77 )%mass=192.22   
e(78 )%sym='Pt' ; e(78 )%mass=195.08   
e(79 )%sym='Au' ; e(79 )%mass=196.97   
e(80 )%sym='Hg' ; e(80 )%mass=200.59   
e(81 )%sym='Tl' ; e(81 )%mass=204.38   
e(82 )%sym='Pb' ; e(82 )%mass=207.2    
e(83 )%sym='Bi' ; e(83 )%mass=208.98   
e(84 )%sym='Po' ; e(84 )%mass=209.     
e(85 )%sym='At' ; e(85 )%mass=210.     
e(86 )%sym='Rn' ; e(86 )%mass=222.     
e(87 )%sym='Fr' ; e(87 )%mass=223.     
e(88 )%sym='Ra' ; e(88 )%mass=226.     
e(89 )%sym='Ac' ; e(89 )%mass=0.       
e(90 )%sym='Th' ; e(90 )%mass=0.       
e(91 )%sym='Pa' ; e(91 )%mass=0.       
e(92 )%sym='U'  ; e(92 )%mass=0.       
e(93 )%sym='Np' ; e(93 )%mass=0.       
e(94 )%sym='Pu' ; e(94 )%mass=0.       
e(95 )%sym='Am' ; e(95 )%mass=0.       
e(96 )%sym='Cm' ; e(96 )%mass=0.       
e(97 )%sym='Bk' ; e(97 )%mass=0.       
e(98 )%sym='Cf' ; e(98 )%mass=0.       
e(99 )%sym='Es' ; e(99 )%mass=0.       
e(100)%sym='Fm' ; e(100)%mass=0.       
e(101)%sym='Md' ; e(101)%mass=0.       
e(102)%sym='No' ; e(102)%mass=0.       
e(103)%sym='Lr' ; e(103)%mass=0.       
e(104)%sym='Rf' ; e(104)%mass=263.     
e(105)%sym='Db' ; e(105)%mass=262.     
e(106)%sym='Sg' ; e(106)%mass=266.     
e(107)%sym='Bh' ; e(107)%mass=264.     
e(108)%sym='Hs' ; e(108)%mass=269.     
e(109)%sym='Mt' ; e(109)%mass=268.     
e(110)%sym='Ds' ; e(110)%mass=272.     
e(111)%sym='Rg' ; e(111)%mass=272.     
e(112)%sym='Cn' ; e(112)%mass=277.     
e(113)%sym='Ut' ; e(113)%mass=284.     
e(114)%sym='Uq' ; e(114)%mass=289.     
e(115)%sym='Up' ; e(115)%mass=288.     
e(116)%sym='Uh' ; e(116)%mass=292.     
e(117)%sym='Us' ; e(117)%mass=291.     
e(118)%sym='Uo' ; e(118)%mass=293.     

! The generix element
e(119)%sym='A' ; e(119)%mass=1.

end subroutine
  
function inq_z(sym)
! Devuelve el id del elemento correspondiente al simbolo sym
use gems_strings, only:locase
use gems_errors, only: werr
integer                   :: inq_z
character(*), intent(in)  :: sym
character(:),allocatable  :: l1,l2

l1 = sym
call locase(l1)
do inq_z = 1, elements%size
  l2 = elements%o(inq_z)%sym
  call locase(l2)
  if (l1==l2)  return
enddo  

inq_z=0

end function inq_z
   
subroutine add_z(sym,mass,carga,rad)
use gems_errors, only:werr
integer                  :: i
character(*), intent(in) :: sym
real(dp),optional        :: mass,carga,rad

! Check already exisitng elements and the string lenght
i=inq_z(sym)
if (i/=0) return

! Add a new element
call elements%append()
i=elements%size
!elements%o(i)%sym=trim(adjustl(sym))
elements%o(i)%sym=trim(sym)

if(present(mass)) elements%o(i)%mass=mass
if(present(mass)) elements%o(i)%carga=carga
if(present(rad)) elements%o(i)%rad=rad

end subroutine add_z

subroutine set_z(z,sym,mass,carga,rad)
integer,intent(in)       :: z
character(*),optional, intent(in) :: sym
real(dp),optional,intent(in)      :: mass,carga,rad

if(present(sym)) elements%o(z)%sym=adjustl(sym)
if(present(mass)) elements%o(z)%mass=mass
if(present(carga)) elements%o(z)%carga=carga
if(present(rad)) elements%o(z)%rad=rad

end subroutine set_z
              
end module gems_elements
