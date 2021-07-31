#include "prep.inc"

#define _Alloc      _PASTE(_NODE,_alloc)
#define _Source     _PASTE(_NODE,_source)
#define _Point      _PASTE(_NODE,_point)
#define _Destroy    _PASTE(_NODE,_destroy)

#define _AddAfter   _PASTE(_NODE,_add_after)
#define _AddBefore  _PASTE(_NODE,_add_before)
#define _Append     _PASTE(_NODE,_append)
#define _DestroyAll _PASTE(_NODE,_destroy_all)
#define _Deattach   _PASTE(_NODE,_deattach)

#ifndef SOFT
#ifndef HARD
#define SOFT
#define HARD
#endif
#endif
      
! ## Double Linked List (DLL) template.
!
! See `../README.md` for information in how to use the template.


! Declaring the node
type :: _NODE
  class(_NODE), pointer :: prev => null()
  class(_NODE), pointer :: next => null()
  _CLASS, pointer       :: o   => null()
  logical               :: hard=.false.
  contains
    procedure :: add_after     => _AddAfter
    procedure :: add_before    => _AddBefore

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
end type _NODE

private _AddAfter    
private _AddBefore
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
#undef _AddBefore
#undef _Append
#undef _DestroyAll 
#undef _Deattach     

#undef _SOFT
#undef _HARD
