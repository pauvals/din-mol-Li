! Subroutines that are common to all structures

#ifndef SOFT
#ifndef HARD
#define SOFT
#define HARD
#endif
#endif
      
#ifdef HARD
  subroutine _Alloc(node)
  class(_NODE),target     :: node
  allocate(node%o)
  node%hard = .true.
  end subroutine _Alloc
                    
  subroutine _Source(node, o)
  class(_NODE),target         :: node
  _CLASS, target, intent(in)  :: o
  allocate(node%o,source=o)
  node%hard = .true.
  end subroutine _Source
#endif
                    
#ifdef SOFT
  subroutine _Point(node, o)
  class(_NODE),target         :: node
  _CLASS, target, intent(in)  :: o
  node%o   =>  o
  node%hard = .false.
  end subroutine _Point
#endif

subroutine _Destroy( node )
! Destroy the node. 
! WARNING: Do not forget to _Deattach first to keep structure integrity.
! WARNING: Do not forget to deallocate the node
type(_NODE)                     :: node

#ifdef HARD
if(node%hard) deallocate(node%o)
#endif

node%o => null()

end subroutine _Destroy

#undef _SOFT
#undef _HARD
