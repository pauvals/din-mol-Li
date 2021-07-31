
program pp
use intrinsics_classes
use lion_class

real(8)                   :: a=3.14

type(real_dlist), target   :: list
class(real_dlist), pointer  :: node

type(lion), target        :: lion1, lion2, lion3
type(lion), pointer       :: lionptr
type(lionking)            :: simba

type(lion_list), target   :: lions
class(lion_list), pointer  :: thislion

type(cdlion_list), target   :: clions
class(cdlion_list), pointer  :: thisclion


print *, 'Real double list:'
print *, 'First 3 nodes have real numbers'
print *, '4th node have a pointer to a real'

! Real list
call list%add_hardcpy( 1.0d0 )
call list%add_hardcpy( 2.6d0 )
call list%add_hardcpy( 3.5d0 )
call list%add_soft( a )

! Simple loop  
node => list
do while( associated(node%next) )
    node => node%next
    print *, node%obj
enddo
print *, ''

print *, 'Deleting second element'
node=>list%next
node=>node%next
call node%del_soft

! WARNING: Do not forget to deallocate the node
deallocate(node)

node => list
do while( associated(node%next) )
    node => node%next
    print *, node%obj
enddo
print *, ''
print *, ''
       


print *, 'Lion list:'
print *, 'First 2 nodes are pointers to the lion objects'
print *, '3rd node has allocated a lion object'
print *, '4th node is a pointer to a lionking object'

! Add lion 1 after the head
lion1%color = "orange"
call lions%add_soft( lion1 )

! Add lion 2 after the head
lion2%color = "red"
call lions%add_soft( lion2 )

! Add lion 3 after the head
lion3%color = "brown"
call lions%add_hardcpy( lion3 )

! Add lionking after the head
simba%color = "gold"
simba%name  = "simba"
call lions%add_soft( simba )

print *, 'Access common properties in the class'
thislion => lions
do while( associated(thislion%next) )
    thislion => thislion%next

    print *, 'a lion ', thislion%obj%color

enddo
print *, ''
              
print *, 'Access all properties using select type.'
thislion => lions
do while( associated(thislion%next) )
    thislion => thislion%next

    ! Print the lion name if is a lionking
    select type (a=>thislion%obj)
    type is (lionking)
      print *, a%name, a%color
    class default  
      ! Print the lion color
      print *, 'simple lion ', thislion%obj%color
    end select

enddo
print *, ''


print *, 'Access all properties using internal procedures (cleaner in my opinion)'
thislion => lions
do while( associated(thislion%next) )
    thislion => thislion%next
    call thislion%obj%printme
enddo
print *, ''

print *, 'Removing pointer to the red lion' 
call lions%del_soft(lion2)
thislion => lions
do while( associated(thislion%next) )
    thislion => thislion%next
    call thislion%obj%printme
enddo
print *, ''

print *, 'Deallocating the brown lion' 
thislion => lions
do while( associated(thislion%next) )
    thislion => thislion%next
    if (thislion%obj%color=='brown') then
      call lions%del_hard(thislion)
      exit
    endif
enddo

thislion => lions
do while( associated(thislion%next) )
    thislion => thislion%next
    call thislion%obj%printme
enddo
print *, ''
print *, ''
       


print *, 'Lion circular double list:'
print *, 'First 2 nodes are pointers to the lion objects'
print *, '3rd node has allocated a lion object'
print *, '4th node is a pointer to a lionking object'

! Initialize the list
call clions%init()

! Add lions
call clions%add_soft( lion1 )
call clions%add_soft( lion2 )
call clions%add_hardcpy( lion3 )
call clions%add_soft( simba )

print *, 'Access common properties in the class'
thisclion => clions
do while( associated(thisclion%next%obj) )
    thisclion => thisclion%next

    print *, 'a lion ', thisclion%obj%color

enddo
print *, ''
              
print *, 'Access all properties using select type.'
thisclion => clions
do while( associated(thisclion%next%obj) )
    thisclion => thisclion%next

    ! Print the lion name if is a lionking
    select type (a=>thisclion%obj)
    type is (lionking)
      print *, a%name, a%color
    class default  
      ! Print the lion color
      print *, 'simple lion ', thisclion%obj%color
    end select

enddo
print *, ''


print *, 'Access all properties using internal procedures (cleaner in my opinion)'
thisclion => clions
do while( associated(thisclion%next%obj) )
    thisclion => thisclion%next
    call thisclion%obj%printme
enddo
print *, ''

print *, 'Deallocating the brown lion' 
thisclion => clions
do while( associated(thisclion%next%obj) )
    thisclion => thisclion%next
    if (thisclion%obj%color=='brown') then
      call thisclion%del_hard()
      exit
    endif
enddo

thisclion => clions
do while( associated(thisclion%next%obj) )
    thisclion => thisclion%next
    call thisclion%obj%printme
enddo
print *, '' 
          

endprogram

