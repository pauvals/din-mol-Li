! The concept of "array of pointers" does not exist in Fortran. The way you have
! to do it is have an array of derived type with a pointer component (Steve
! Lionel, https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/280765).

! Declaring the node
type :: _NODE
  _CLASS, pointer        :: o   => null()
end type _NODE

#undef _NODE
#undef _CLASS
 
