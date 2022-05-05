#include "prep.inc"

#define _Create   _PASTE(_NODE,_create)
#define _Append   _PASTE(_NODE,_append)
#define _Put      _PASTE(_NODE,_put)
#define _Delete   _PASTE(_NODE,_delete)
#define _Destroy  _PASTE(_NODE,_destroy)
#define _Increase _PASTE(_NODE,_increase)
#define _Assign   _PASTE(_NODE,_assign)


! If the _TYPE is polymorhpic, it must have an assigment procedure or
! it will give an error like:
!
! Nonallocatable variable must not be polymorphic in intrinsic assignment at (1)
! - check that there is a matching specific subroutine for '=' operator
!

type, public :: _NODE
  integer               :: size = 0

  ! I want _TYPE to be pointer to allow point to 
  ! a slice of the vector. However, this implies that
  ! the `allocated` intrinsic function can not be used.
  ! usually this will require to declare a boolean to know if the
  ! pointer is hard or soft linked (allocated or associated). However
  ! here the association does not have sense, so it is important
  ! to know that:
  ! WARNING: Assuming that _TYPE pointer is always allocated
  ! and never associated.
  ! _TYPE, allocatable    :: o(:) 
  _TYPE, pointer        :: o(:)=>null()
  real(4)               :: growth_rate = 1.1
  contains
    procedure :: init     => _Create
    procedure :: append   => _Append
    procedure :: put      => _Put
    procedure :: del      => _Delete
    procedure :: destroy  => _Destroy
    procedure :: ensure   => _Increase
    procedure :: _Assign
    generic   :: assignment(=) => _Assign
end type _NODE

private  _Create   
private  _Append   
private  _Put      
private  _Delete   
private  _Destroy  
private  _Increase
private  _Assign

#undef _NODE
#undef _TYPE
#undef _Create  
#undef _Append  
#undef _Put     
#undef _Delete  
#undef _Destroy 
#undef _Assign
#undef _Increase

