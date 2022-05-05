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


module gems_program_types
! TODO: Module dedicated to BOX only
use gems_constants,only:sp,dp,dm,find_io
use gems_elements,only:elements,element,inq_z

implicit none
save
private

! The program and subprograms trace this variable to get the term order
! After the term_signal is read, is set to .false. too allow subprograms end
! without quit gems.
logical,public :: term_signal=.false., cycle_signal=.false.

! Program variables
real(dp),target,public   :: dm_steps=0._dp
real(dp),public,target   :: time=0.0_dp

integer,public           :: frame=0 ! Cuando escribe: el frame actual del outunit
integer,public           :: nstep=100 ! number of integration step

real(dp),target,public   :: dt=0.002

public :: box_setvars,box_expand,boxed
real(dp),public,dimension(dm,dm)     :: tbox   =0.0_dp
real(dp),public,dimension(dm),target :: box    =1.0e6_dp  
real(dp),public,dimension(dm)        :: box2   =5.0e5_dp  , &
                                        one_box=1.0e-6_dp , &
                                        one_box2=2.5e-5_dp
real(dp),public                      :: box_vol=1.0e18_dp
logical                              :: boxed=.false.

real(dp),public                      :: virial(dm,dm)=0._dp
logical,public                       :: b_gvirial=.false.

logical,public                       :: mic=.true.

public :: distance
public :: distance2line
 
! Neigh boxes to a box
! ====================

integer,parameter,public,dimension(26,3) :: n1cells = transpose(reshape( &
                   [ 1, 0, 0, -1, 0, 0,  0, 1, 0,&
                     1, 1, 0, -1, 1, 0,  0,-1, 0,&
                     1,-1, 0, -1,-1, 0,  0, 0, 1,&
                     1, 0, 1, -1, 0, 1,  0, 1, 1,&
                     1, 1, 1, -1, 1, 1,  0,-1, 1,&
                     1,-1, 1, -1,-1, 1,  0, 0,-1,&
                     1, 0,-1, -1, 0,-1,  0, 1,-1,&
                     1, 1,-1, -1, 1,-1,  0,-1,-1,&
                     1,-1,-1, -1,-1,-1 ],[3,26]))
 
! Others
! ======

contains
                      
! Box
! ===
           
subroutine box_setvars()

boxed=.true.

! For non-cubic boxes this should be box(:)=matmul(tbox(i,i),[1,1,1])
box(1)=tbox(1,1)
box(2)=tbox(2,2)
box(3)=tbox(3,3)

one_box(:) = 1.0_dp/box(:)

! For non-cubic boxes this should be the determinant of tbox
box_vol = box(1)*box(2)*box(3)

end subroutine box_setvars

subroutine box_expand(a,b,c)
real(dp)       :: a,b,c

tbox(1,1)=tbox(1,1)*a
tbox(2,2)=tbox(2,2)*b
tbox(3,3)=tbox(3,3)*c

call box_setvars()

end subroutine box_expand

! PBC
! ---

function itrans(r,pbc) result(cr)
! Gives integer with translation needs
! Note: It can not be `elemental` because it reads box(:)
real(dp),intent(in)     :: r(dm)
logical,intent(in)      :: pbc(dm)
real(dp)                :: cr(dm)
integer                 :: i

cr(:)=0
do i=1,dm
  if (pbc(i)) then
    if(r(i)>=box(i)) then
      cr(i)=cr(i)+1
    elseif(r(i)<0._dp) then
      cr(i)=cr(i)-1
    endif
  endif
enddo
     
end function itrans
 
function rtrans(r,pbc)
real(dp),intent(in)     :: r(dm)
logical,intent(in)      :: pbc(dm)
real(dp)                :: rtrans(dm)
rtrans(:)=r(:)-box(:)*itrans(r,pbc)
end function rtrans

function distance(r1,r2,pbc) result(vd)
!calculates the distance between two points considering pbc.
! TODO: elemental? check speed
real(dp),intent(in)     :: r1(dm),r2(dm)
logical,intent(in)      :: pbc(dm)
real(dp)                :: vd(dm)
integer                 :: l

! Mas rapido usar idnint que un if
vd(:)=r2(:)-r1(:)
do l = 1,dm
  if (pbc(l)) vd(l)=vd(l)-box(l)*idnint(vd(l)*one_box(l))
enddo

end function distance

function distance2line(r,p,v,pbc) result(vd)
!calculates the distance between a line and a points considering pbc.
real(dp),intent(in)     :: r(dm),p(dm),v(dm)
logical,intent(in)      :: pbc(dm)
real(dp)                :: vd(dm),aux(dm)

aux=distance(r,p,pbc)
vd(:)=aux(:)-(dot_product(aux,v(:)))*v(:)

end function distance2line
                  
end module gems_program_types
 
