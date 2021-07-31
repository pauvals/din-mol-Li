! Defining some objects
module lion_class

public
type :: lion
  character(20) :: color
  contains
    procedure   :: printme=>printlion
end type

type, extends(lion) :: lionking
  character(20) :: name
  contains
    procedure   :: printme=>printlionking
end type    
  
! Creating a lion list
#define _NODE lion_list
#define _CLASS class(lion)
#include "list_header.inc"
 
! Creating a lion double list
#define _NODE lion_dlist
#define _CLASS class(lion)
#include "dlist_header.inc"
 
! Creating a lion list
#define _NODE lion_cdlist
#define _CLASS class(lion)
#include "cdlist_header.inc"

! Creating a lion circular double list
#define _NODE lion_vector
#define _TYPE type(lion)
#include "vector_header.inc"
                
! Creating a lion circular double list
#define _NODE lion_aop
#define _CLASS class(lion)
#include "arrayofptrs_header.inc"
                
contains
  
! Creating a lion list
#define _NODE lion_list
#define _CLASS class(lion)
#include "list_body.inc"
 
! Creating a lion double list
#define _NODE lion_dlist
#define _CLASS class(lion)
#include "dlist_body.inc"
 
! Creating a lion list
#define _NODE lion_cdlist
#define _CLASS class(lion)
#include "cdlist_body.inc"

! Creating a lion circular double list
#define _NODE lion_vector
#define _TYPE type(lion)
#include "vector_body.inc"

subroutine printlion(l)
  class(lion)      :: l
  print *, 'simple lion ', l%color
end subroutine

subroutine printlionking(l)
  class(lionking)      :: l
  print *, l%name, l%color
end subroutine

end module lion_class
 
