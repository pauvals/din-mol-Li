#define _IDENTITY(x) x
#define _PASTE(a,b) _IDENTITY(a)_IDENTITY(b)
#define _STR(name) "name"

! #define _PASTE(a,b) a ## b
! #define _STR(name) #name
      
