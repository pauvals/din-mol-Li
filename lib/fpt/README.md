# Fortran Preprocessor Templates for Dynamic Data Structures (FPT-DDS)

Here I present my way to implement _fortran templates_ to create **dynamic data
structures using preprocessor directives**.

A data structure for an arbitrary data type can be easily constructed by
defining a few preprocessor variables and including the corresponding files
distributed in the `include` folder as it is shown in [the simple example
below](#-A-simple-example).

Another examples of use can be found in the `src` folder. To compile this
example just type:

    ./configure
    make

# A simple example

To compile the following example, do not forget to add
`-I$(path_to_include_files)` and use `.F90` extension to activate preprocessing.

    module intrinsic_class

    ! Here we declare the type of the dynamic data structure (DDS)
    ! 1. _NODE define the name of the data type of the DDS
    ! 2. _CLASS or _TYPE is the data type of the objects handle by the DDS. 
    !   If the DDS use _CLASS, it can be polymorphic. If not, it use _TYPE
    ! 3. include file name allow to select between the different DDS, e.g linked list, double linked list, etc. 
    #define _NODE integer_dlist
    #define _CLASS integer
    #include "dlist_header.inc"

    contains

    ! Here we declare the procedures of the new integer_dlist DDS.
    ! Definition should match the one used in the header
    #define _NODE integer_dlist
    #define _CLASS integer
    #include "dlist_body.inc"

    end module intrinsic_class


    program example
    use intrinsic_class
    type(integer_dlist),target		:: rlist
    type(integer_dlist),pointer		:: ptr
    integer,target                :: i
    integer                       :: j

    i=3

    ! Add, allocate and copy and item (i.e. allocated in the list)
    ! rlist: head -> 2  
    call rlist%add_after()
    call rlist%next%source(2)

    ! Add second item
    ! rlist: head -> 3 -> 2  
    call rlist%add_after()
    call rlist%next%source(i)

    ! Add third item as a pointer (i.e. allcoated elswhere)  
    ! rlist: head -> o -> 3 -> 2   with o -> i
    call rlist%add_after()
    call rlist%next%point(i)

    ! Print the list (output: 3 3 2)
    ptr=>rlist
    do j=1,3
      ptr=>ptr%next
      print *, ptr%o
    end do  

    ! Change the target of the first node
    i=4

    ! Print the list (output: 4 3 2)
    ptr=>rlist
    do j=1,3
      ptr=>ptr%next
      print *, ptr%o
    end do  
     
    end program example


# HARD and SOFT modes

`HARD` and `SOFT` modes can be activated by a preprocesor variable. For example

    #SOFT
    #define _NODE integer_dlist
    #define _CLASS integer
    #include "dlist_body.inc"

A `SOFT` declaration means that objects in the structure can only be pointers.

A `HARD` declaration means that objects will be allocated in each node.

If neither `SOFT` or `HARD` is used, then the strucutre can hold pointers or
allocations.


# Supported Dynamic Data Structures

## Linked list template (LL, polymorphic)

      #declare _NODE <LL name>
      #declare _CLASS <data type>
      #include list_header.inc

Any other type or object can be used. A common way to iterate might be:

    type(_NODE),pointer   :: head
    class(_NODE),pointer  :: node
   
    ...

    node => head
    do while(associated(node))
      ... !work with node%obj ("cycle" do not advance the node)
      node => node%next
    enddo

The problem with this iteration is that a `cycle` fortran keyword will skip the
advance of the node pointer to the next element. To fix that, I just skip
(avoid to use) the data object associated with first node of the list (the
`head` node) and start the iteration in this way:
  
    node => head
    do while(associated(node%next))
      node => node%next
      ... !work with node%obj ("cycle" advance the list)   
    enddo

## Double linked list template (DLL, polymorphic).

Same that LL but it also allows to iterate backwards:
   
    node => tail
    do while(associated(node))
      ... !work with node%obj
      node => node%prev
    enddo

To make `cycle` keyword work, see [LL template](##-Linked-list-template-(LL,-polymorphic))


## Circular double linked list template (CDLL, polymorphic).

CDLL might not have any particular head, so to iterate it requires to keep a
total count of the items to avoid infinite loops. Leaving the first node with a
null obj might be a way to avoid tracking the number of items: 
    
    node => head%next
    do while(associated(node%obj))
      ... !work with node%obj
      node => node%next
    enddo
 
CDLL simplify the procedures of adding and removing nodes avoiding the
association check needed for the beginning and final nodes of a DLL. 

CDLL requires a constructor, which is in contrast to LL and DLL. 

Other way to iterate that does not require empty head

    node => head%prev
    do 
      node => node%next

      ... !work with node%obj

      if(associated(node%next,alist)) exit
    enddo


## Vector

An array that automatically reallocate when needed.

## Array of Pointers (AOP, polymorphic)

An array of a user defined type that contains a pointer. With this structure
each element of the array might point to a different target.

# License

Fortran Preprocessor Templates for Dynamic Data Structures (FPT-DDS) is hosted
in [github](https://github.com/alexispaz/FortranTemplates) with a BSD 3-Clause License. 

