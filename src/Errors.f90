! Copyright (c) 2020  Sergio Alexis Paz
!
!  This file is part of GEMS. GEMS is an Extensible Molecular Simulator.
!	 .
!  GEMS is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!  .
!  GEMS is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  .
!  You should have received a copy of the GNU General Public License
!  along with GEMS.  If not, see <https://www.gnu.org/licenses/>.

 
module gems_errors
use gems_constants,     only: dp,dm,linewidth, time1
use, intrinsic :: iso_fortran_env, only: input_unit, output_unit, error_unit
implicit none

private 
public  :: msn,wlog,wwan,werr,wstd,silent,wprt,wref
public  :: timer_start, timer_dump

integer,public,target    :: logunit=output_unit
integer,public           :: inunit=input_unit,truelogunit=output_unit
integer,public,pointer   :: printunit=>null()

character(linewidth)     :: msn=''

! timer variables
real(dp)       :: vclock,sclock
integer,public :: sclock_rate, sclock_max,sclock_t1,sclock_t2
 
! Error message and error flag 
character(:),allocatable    :: errm
logical,public              :: errf

! Silent switch makes errors or messages to be supressed. When silent is on,
! errors will not halt the execution of the program. This provide a simple way to handle
! errors. Only first message will be save in errm variable.
! TODO Buffer messages in FIFO pile.
logical        :: silent=.false.

! True right after a message is set, false if not.
! TODO: switch stats to warm and stdm
logical,public :: stats

! True after a wstd/wwan occurs. Manually set to .false. is required.
! This is usefull for group a block of wwan calls.
! TODO: switch wstats and sstats to warm and stdm
logical,public :: wstats, sstats

contains
     
subroutine werr(message,logic)
! Test logic, report error and halt the program
logical,intent(in)           :: logic
character(*),intent(in)      :: message

errf=.false.
if(.not.logic) return

errf=.true.
errm='#-ERR-> '//message
           
if(silent) return

write(unit=error_unit, fmt='(a)') errm
write(logunit,'(a)') errm

! TODO:
! integer, parameter :: ALLOCATION_ERR = 1, &
!                       FILE_OPEN_ERR = 2, &
!                       CMD_LINE_NR_ARGS_ERR = 3, &
!                       CMD_LINE_ARG_VALUE_ERR = 4, &
!                       FILE_VALUE_ERR = 5
! stop ALLOCATION_ERR
stop 1

! TODO:
! Using compile macros 'at __FILE__ __LINE__'
 
endsubroutine werr

subroutine wstd(message,logic)
logical,intent(in),optional      :: logic
character(*),intent(in),optional :: message

stats=.false.
if(present(logic)) then
  if (.not.logic) return
endif
stats=.true.
sstats=.true.

if(silent) return

write(logunit,'(a)',advance='no') '#-STD-> '
if (present(message))  write(logunit,'(a)') trim(message)

call flush(logunit)

endsubroutine wstd

subroutine wwan(message,logic)
logical,intent(in),optional      :: logic
character(*),intent(in),optional :: message

stats=.false.
if(present(logic)) then
  if (.not.logic) return
endif
stats=.true.
wstats=.true.

if(silent) return

write(logunit,'(a)',advance='no') '#-WAN-> '
if (present(message))  write(logunit,'(a)') trim(message)

call flush(logunit)


endsubroutine wwan

subroutine wlog(flag,message,logic)
logical,intent(in),optional      :: logic
character(*),intent(in),optional :: message
character(*)                     :: flag

if(silent) return

if(present(logic)) then
  if (.not.logic) return
endif

write(logunit,'(a)',advance='no') '# '//flag//' '
if (present(message))  write(logunit,'(a)') trim(message)

call flush(logunit)

endsubroutine wlog
         
subroutine wprt(message)
character(*),intent(in),optional :: message

if(associated(printunit,target=logunit)) write(printunit,'(a)',advance='no') '# '
if (present(message))  write(printunit,'(a)') trim(message)
call flush(printunit)

endsubroutine wprt
            
subroutine timer_start(stime)
real(dp),intent(in),optional    :: stime
if(present(stime)) vclock=stime
call system_clock(sclock_t1)
end subroutine

subroutine timer_dump(ns,stime,nup)
real(dp),intent(in),optional    :: stime ! Simulated time
integer,intent(in)    :: ns    ! time step
integer,intent(in),optional    :: nup   ! number of nieghboor updates
real(dp)              :: a

! Only print when ns is 10,100,1000,...
a=log10(real(ns,dp))
if (a/=int(a)) return

! Print simulation information
call wstd(); write(logunit,'(i0," steps reached!")') ns
if(present(nup)) call wstd(); write(logunit,'(i0, "neigh searches")') nup
 
! Get cpu time
call system_clock(sclock_t2)
sclock=(sclock_t2-sclock_t1)/real(sclock_rate)/ns

if (sclock<60) then
  call wstd(); write(logunit,'("average wall time per step: ",f10.3," s")') sclock
else
  call wstd(); write(logunit,'("average wall time per step: ",i0," m ",f10.3," s")') &
                        int(sclock/60.0_dp),mod(sclock,60.d0)
endif 

! Simulated time information
if(.not.present(stime)) return

call wstd(); write(logunit,'(f10.2,"ps simulated!")') stime

sclock=sclock*ns/(stime-vclock)*1.e3

if (sclock<60) then
  call wstd(); write(logunit,'("average wall time per ns: ",f10.3," s")') sclock
elseif(sclock<3600) then
  call wstd(); write(logunit,'("average wall time per ns: ",i0," m ",f10.3," s")') &
                        int(sclock/60.0_dp),mod(sclock,60.d0)
else
  call wstd(); write(logunit,'("average wall time per ns: ",i0," h ",i0," m")') &
                        int(sclock/3600.0_dp),int(mod(sclock,3600.d0)/60.0_dp)
endif  
         
call flush(logunit)

end subroutine timer_dump

subroutine wref(cite)
character(*),intent(in),optional :: cite

call wlog('REF','This feature use the following references:')

select case(cite)
case('ermak')
  call wlog('REF','- D. L. Ermak and H. Buckholz, “Numerical integration of the Langevin equation:   ')
  call wlog('REF','  Monte Carlo simulation”, J. Comput. Phys. 35, 169 (1980)                        ')
  call wlog('REF','- M. P. Allen and D. J. Tildesley, "Computer simulation of liquids" Chap 9, "Brownian Dinamics", Pag 263  ')
  call wlog('REF','- I. Snook, "The ermak and Generalised Langevin Apprroach to the Dynamics of     ')
  call wlog('REF','  Atomic, Polymeric and Colloidal System". Elsevier, Sec 6.2.4 "A                 ')
  call wlog('REF','  third first-order BD algorithm", Pag 118.                                       ')
endselect

endsubroutine wref
           
end module gems_errors
 
