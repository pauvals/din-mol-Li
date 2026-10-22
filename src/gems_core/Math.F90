! Minimal gems_math subset required by the vendored GEMS neighbor core.
module gems_math
implicit none
private
public :: binleft

contains

subroutine binleft(a,val,m,found)
! Find the index m where a(m)<=val<a(m+1) using a binary search.
! The function returns true if a(m)=val and false otherwise.
! a(:) must be sorted.
integer,intent(in)  :: a(:),val
integer             :: r,l,j
integer,intent(out) :: m
logical,intent(out) :: found

l=lbound(a,1)-1
r=ubound(a,1)+1

do while (r>l+1)
  m=(r+l)/2
  j=a(m)

  if(val==j) then
    found=.true.
    return
  endif

  if(val>j) then
    l=m
  else
    r=m
  endif
enddo

if(r<=size(a)) then
  if(a(r)==val) then
    m=r
    found=.true.
    return
  endif
endif

found=.false.
m=l

end subroutine binleft

end module gems_math
