! For more information read list_header.inc.

#include "prep.inc"

#define _Alloc      _PASTE(_NODE,_alloc)
#define _Source     _PASTE(_NODE,_source)
#define _Point      _PASTE(_NODE,_point)
#define _Destroy    _PASTE(_NODE,_destroy)
      
#define _AddAfter   _PASTE(_NODE,_add_after)
#define _AddHead    _PASTE(_NODE,_add_head)
#define _Append     _PASTE(_NODE,_append)
#define _DestroyAll _PASTE(_NODE,_destroy_all)
#define _Deattach   _PASTE(_NODE,_deattach)
      
#include "object.inc"
    
! Hard allocation of a type migth required the allocation of its internal variables. If this allocation are slow (huge internal
! variables) it might be wanted to allocate the object and not its internal variables. May be it is possible to declare an
! allocation procedure for the type that requires minima allocation size and then to increase the size a reallocation should take
! place. However this is not considered here.
! 
! Soft and hard booleans allow for massive destroy procedures like _DestroyAll However, If _CLASS correspond to a polymorfphic
! object, the hard allocations without source will fail (i.e. `class(*)`) or only allow the allocation the parent class (e.g.
! class(animals). Thus, procedures that creates hard nodes do not allocate the object by themself unless a source is included.
! Therefore, the hard allocation without source should be something like:
!
! new=_AppendHard(node)
! allocate(node%o)
!
! Thus, hard boolean is not the only thing to check at destroy procedures. Associated should be used also if a robust deallocation is
! required.

subroutine _AddAfter(node, opt)
!Insert a new node after `node`
class(_NODE),target                          :: node
type(_NODE), pointer                         :: new
type(_NODE), pointer, optional, intent(out)  :: opt

allocate(new)

new%next  => node%next
node%next => new

if(present(opt)) opt=>new

end subroutine _AddAfter

function _AddHead(node) result(new)
!Insert a new head before node
!On output, node is the new head
class(_NODE),target       :: node
type(_NODE), pointer      :: new

allocate(new)
new%next => node

end function _AddHead

subroutine _Append(node, opt)
class(_NODE),target             :: node
type(_NODE),pointer             :: last
type(_NODE), pointer, optional  :: opt

! Find the last node
last => node
do while( associated(last%next) )
  last => last%next
enddo    

if(present(opt)) then
  call _AddAfter(last, opt)
else
  call _AddAfter(last)
endif
            
end subroutine _Append
         
subroutine _Deattach( node, prev )
class(_NODE)            :: node
type(_NODE)             :: prev
prev%next => node%next
end subroutine _Deattach
 
subroutine _DestroyAll( head )
! Destroy the entire list.
! WARNING: Do not forget to deallocate the node
class(_NODE),target     :: head
class(_NODE),pointer   :: aux

do while(associated(head%next))
  aux => head%next
  head%next => aux%next
  deallocate(aux)
enddo

call _Destroy(head)
! deallocate(node) ! final

end subroutine _DestroyAll

#undef _NODE
#undef _CLASS

#undef _Alloc
#undef _Source
#undef _Point
#undef _Destroy

#undef _AddAfter
#undef _AddHead
#undef _Append
#undef _DestroyAll 
#undef _Deattach   
