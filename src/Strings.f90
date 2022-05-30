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

 
module gems_strings
use gems_constants, only:dp

implicit none

private

interface str
  module procedure int2char,float2char
end interface

public    :: upcase,locase,str,int2char0
public    :: contar_cifras

character(len=*),parameter,public ::    &
  chset_u="ABCDEFGHIJKLMNOPQRSTUVWXYZ", & ! upper
  chset_l="abcdefghijklmnopqrstuvwxyz", & ! lower
  chset_n="0123456789",                 & ! number
  chset_ln=chset_l//chset_n,            & ! lower and numbers
  chset_lun=chset_l//chset_u//chset_n,  & ! upper lower and numbers
  chset_var=chset_ln//"_[]:"             ! allowed in a variable name

contains

function contar_cifras(i) result(c)
! Count digits of integers.
! Negative integers returns digits+1
integer,intent(in)          :: i
integer                     :: a,c

a=abs(i)
if(a==0) then
  c=1
  return
endif

c=floor(log10(real(a))+1)
if(i<0) c=c+1

end function  

function int2char(i,fmt_)
! Convert integer to string
! https://stackoverflow.com/a/30103608/1342186
character(:),allocatable           :: int2char
integer,intent(in)                 :: i
character(*),intent(in),optional   :: fmt_
character(:),allocatable           :: wfmt
character(range(i)+2)              :: tmp
  
if(present(fmt_)) then
  wfmt=fmt_
else
  wfmt = '(i0)'
endif
                      
write(tmp,wfmt) i
int2char = trim(tmp)
end function  

function float2char(i,fmt_)
! Convert float to string
character(:),allocatable           :: float2char
real(dp),intent(in)                :: i
character(*),intent(in),optional   :: fmt_
character(:),allocatable           :: wfmt
character(50)                      :: tmp

if(present(fmt_)) then
  wfmt=fmt_
else
  wfmt = '(e25.12)'
endif

write(tmp,wfmt) i
float2char = trim(tmp)
end function   

subroutine upcase(word)
! Change the word to upercase
character(len=*), intent(inout) :: word
integer :: i,k

do i=1,len(word)
  k=index(chset_l,word(i:i))
  if (k .ne. 0) word(i:i)=chset_u(k:k)
end do

end subroutine upcase

subroutine locase(word)
! Change the word to lowercase
character(len=*), intent(inout) :: word
integer :: i,k

do i=1,len(word)
  k=index(chset_u,word(i:i))
  if (k .ne. 0) word(i:i)=chset_l(k:k)
end do
end subroutine locase


function int2char0(i1,i2) result(c)
! Like int2char but add ceros on the left until 
! fill the same digit number of num2.
integer,intent(in)          :: i1,i2
integer                     :: a,i
character(:),allocatable    :: c
a=contar_cifras(i2)
allocate(character(a) :: c)
write(c,'(i0)') i1
a=a-contar_cifras(i1)
if(a<=0) return
do i=1,a
  c(i:i) = '0' 
enddo
end function  



end module gems_strings 


