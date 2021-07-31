! For more information read cdlist_header.inc.

#include "prep.inc"

#define _Alloc      _PASTE(_NODE,_alloc)
#define _Source     _PASTE(_NODE,_source)
#define _Point      _PASTE(_NODE,_point)
#define _Destroy    _PASTE(_NODE,_destroy)
       
#define _Constructor      _PASTE(_NODE,_constructor)
#define _AddAfter   _PASTE(_NODE,_add_after)
#define _AddBefore  _PASTE(_NODE,_add_before)
#define _DestroyAll _PASTE(_NODE,_destroy_all)
#define _Deattach   _PASTE(_NODE,_deattach)
      
#include "object.inc"

subroutine _Constructor(node)
class(_NODE),target        :: node
node%next => node
node%prev => node
node%o  => null()
node%hard=.false.
end subroutine _Constructor

subroutine _AddAfter(node, opt)
!Insert a new node after `node`
class(_NODE),target                          :: node
type(_NODE), pointer                         :: new
type(_NODE), pointer, optional, intent(out)  :: opt

allocate(new)

new%next  => node%next
new%prev  => node
node%next%prev => new
node%next => new

if(present(opt)) opt=>new

end subroutine _AddAfter
            
subroutine _AddBefore(node, opt)
!Insert a new node bedore `node`
class(_NODE),target                          :: node
type(_NODE), pointer                         :: new
type(_NODE), pointer, optional, intent(out)  :: opt

allocate(new)

new%next => node
new%prev => node%prev
node%prev%next => new
node%prev => new

if(present(opt)) opt=>new

end subroutine _AddBefore
 
subroutine _Deattach( node )
class(_NODE)   :: node
node%prev%next => node%next
node%next%prev => node%prev
end subroutine _Deattach
                      
subroutine _DestroyAll(node)
! Destroy the entire list.
! WARNING: Do not forget to deallocate the node
class(_NODE),target        :: node
class(_NODE),pointer       :: aux


do while(.not.associated(node%next,target=node))
  aux => node%next
  call _Deattach(aux)
  deallocate(aux)
enddo

call _Destroy(node)
! deallocate(node) ! final
             
end subroutine _DestroyAll
        
#undef _NODE
#undef _CLASS
#undef _Constructor

#undef _Alloc
#undef _Source
#undef _Point
#undef _Destroy
 
#undef _AddAfter
#undef _AddBefore
#undef _DestroyAll 
#undef _Deattach   
