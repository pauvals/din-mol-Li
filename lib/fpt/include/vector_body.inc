#include "prep.inc"

#define _Create   _PASTE(_NODE,_create)
#define _Append   _PASTE(_NODE,_append)
#define _Put      _PASTE(_NODE,_put)
#define _Delete   _PASTE(_NODE,_delete)
#define _Destroy  _PASTE(_NODE,_destroy)
#define _Increase _PASTE(_NODE,_increase)
#define _Assign   _PASTE(_NODE,_assign)

! Dummy might be necesary in some cases. For example, if the object is a
! character with certain length and we want to use the add subrroutine with a
! shorter string we can set _DUMMY as character(*)
#ifndef _DUMMY
#define _DUMMY _TYPE
#endif
        
subroutine _Assign( a, b )
  class(_NODE),intent(out)      :: a
  class(_NODE),intent(in)       :: b
  integer                       :: i

  ! check that the vector does not have any data left
  if (a%size<b%size) then
    ! call a%destroy()
    call a%init(b%size)
  endif

  do i = 1,b%size
    a%o(i)=b%o(i)
  enddo

  a%size=b%size
  a%growth_rate=b%growth_rate

end subroutine _Assign

subroutine _Create( vec, capacity )
  ! Create a new vector with certain initial capacity
  class(_NODE)       :: vec
  integer, optional  :: capacity
  integer            :: cap

  ! WARNING: Assuming that _TYPE pointer is always allocated
  ! and never associated.
  if (associated(vec%o)) call _Destroy(vec)

  if (present(capacity)) then
      cap=max(1,capacity)
  else
      cap=10
  endif

  allocate(vec%o(1:cap))
  vec%size=0

end subroutine _Create

subroutine _Destroy( vec )
    class(_NODE)     :: vec

    ! WARNING: Assuming that _TYPE pointer is always allocated
    ! and never associated.
    if (associated(vec%o))  deallocate(vec%o)

    vec%o=>null()

    vec%size = 0
end subroutine _Destroy


!subroutine _insert_empty( vec, pos, number )
!    ! insert 'number' empty elements started from position 'pos'
!    class(_NODE)       :: vec
!    integer, intent(in)  :: pos
!    integer, intent(in)  :: number
!
!    integer              :: i
!
!    if ( number<1 .or. pos<1 .or. pos>vec%size ) return
!
!    if ( vec%size+number>=size(vec%o) ) then
!        call _increase_capacity( vec, vec%size+number )
!    endif
!
!    do i = vec%size,pos,-1
!        vec%o(i+number) = vec%o(i)
!    enddo
!
!    do i = 1,number
!        vec%o(pos+i-1) = empty__data
!    enddo
!
!    vec%size = vec%size + number
!end subroutine _insert_empty

subroutine _Delete( vec, pos, n )
    ! delete 'n' element begining with position 'pos'
    class(_NODE)         :: vec
    integer, intent(in)  :: pos
    integer, intent(in)  :: n
    integer              :: i

    if ( n<1 .or. pos<1 .or. pos>vec%size ) return

    do i = pos,vec%size-n
        vec%o(i) = vec%o(i+n)
    enddo

    vec%size = vec%size - n
end subroutine _Delete


subroutine _Append( vec, data )
  ! Append 'data' to vector 'vec'
  class(_NODE)                :: vec
  _DUMMY,optional,intent(in)   :: data

  call _Increase( vec, vec%size+1 )
  vec%size = vec%size + 1

  if(present(data)) vec%o(vec%size) = data

end subroutine _Append

subroutine _Put(vec,n,data)
  class(_NODE)     :: vec
  integer          :: n
  _DUMMY,optional,intent(in)   :: data

  if (n<1) return
  call _Increase(vec,n)
  vec%size = max(vec%size,n)

  if(present(data)) vec%o(n) = data

end subroutine _Put

subroutine _Increase( vec, capacity )
  !Expand the array holding the data
  !Only for internal use (note that vec%size does not change)
  !     capacity      minimum capacity
  class(_NODE),intent(inout)     :: vec
  integer                        :: capacity
  integer                        :: new_cap
  _TYPE, allocatable, target     :: aux(:)
  ! integer                        :: i

  if ( capacity <= size(vec%o) ) return

  !new_cap = max(capacity,nint(vec%growth_rate*size(vec%o)))
  new_cap = int(capacity*vec%growth_rate)
  
  allocate( aux(1:vec%size) )
  aux(1:vec%size) = vec%o(1:vec%size)
  ! do i=1,vec%size
  !   aux(i) = vec%o(i)
  ! enddo
  deallocate(vec%o)
  allocate(vec%o(new_cap))
  vec%o(1:vec%size) = aux(1:vec%size)
  ! do i=1,vec%size
  !   vec%o(i) = aux(i)
  ! enddo

  ! It seams that move alloc works if vec%o is allocatable (not pointer)
  ! call move_alloc(to=vec%o,from=aux)

end subroutine _Increase

#undef _NODE
#undef _DUMMY
#undef _TYPE
#undef _Create  
#undef _Append  
#undef _Put     
#undef _Delete  
#undef _Destroy 
#undef _Assign
#undef _Increase
 
