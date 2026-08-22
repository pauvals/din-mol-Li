module dana_parallel
  implicit none

  private
  public :: parallel_init, parallel_finalize, parallel_barrier
  public :: parallel_image, parallel_images, parallel_is_master
  public :: parallel_enabled, parallel_running_multi_image

  logical, save :: initialized = .false.
  logical, save :: parallel_enabled = .false.
  integer, save :: parallel_image = 1
  integer, save :: parallel_images = 1

contains

subroutine parallel_init()
#ifdef DANA_COARRAY
  parallel_enabled = .true.
  parallel_image = this_image()
  parallel_images = num_images()
  sync all
#else
  parallel_enabled = .false.
  parallel_image = 1
  parallel_images = 1
#endif
  initialized = .true.
end subroutine parallel_init

subroutine parallel_finalize()
  if (.not.initialized) return
#ifdef DANA_COARRAY
  sync all
#endif
  initialized = .false.
end subroutine parallel_finalize

subroutine parallel_barrier()
#ifdef DANA_COARRAY
  sync all
#endif
end subroutine parallel_barrier

logical function parallel_is_master()
  parallel_is_master = parallel_image == 1
end function parallel_is_master

logical function parallel_running_multi_image()
  parallel_running_multi_image = parallel_images > 1
end function parallel_running_multi_image

end module dana_parallel
