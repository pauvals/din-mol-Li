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

#ifndef SOFT
#ifndef HARD
#define SOFT
#define HARD
#endif
#endif
      
! ## Linked List (LL) template.
!
! See `../README.md` for information in how to use the template.


! Declaring the node
type :: _NODE
  class(_NODE), pointer :: next => null()
  _CLASS, pointer        :: o   => null()
  logical               :: hard=.false.
  contains
    procedure :: add_after     => _AddAfter
    procedure :: add_head      => _AddHead

#ifdef HARD    
    procedure :: alloc         => _Alloc
    procedure :: source        => _Source
#endif    

#ifdef SOFT    
    procedure :: point         => _Point
#endif

    procedure :: append        => _Append
    procedure :: deattach      => _Deattach
    procedure :: destroy_all   => _DestroyAll
    final     :: _Destroy

    ! procedure :: deattach           => _Deattach
    ! Note that deattach is not used because it would consist only of 1 line:
    !   prev => node%next
    ! Where prev is the previous pointer and node is the node to deattach.

end type _NODE

private _AddAfter    
private _AddHead
#ifdef HARD
  private _Alloc       
  private _Source      
#endif
#ifdef SOFT
  private _Point       
#endif
private _Append      
private _DestroyAll  
private _Destroy
private _Deattach    

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

#undef _SOFT
#undef _HARD
