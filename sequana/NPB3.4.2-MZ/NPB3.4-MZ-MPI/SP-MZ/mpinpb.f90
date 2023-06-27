!---------------------------------------------------------------------
!---------------------------------------------------------------------
!
!  mpinpb module
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      module mpinpb

      use sp_data, only : max_zones, kind2
!
      include 'mpif.h'
!
!     zone_proc_id(MZ)     - process id each zone assigned to
!     proc_zone_id(MZ)     - list of zones assigned to this process
!     proc_num_zones       - number of zones assigned to this process
!     proc_zone_count(NP)  - number of zones assigned to each process
!     proc_num_threads(NP) - number of threads assigned to each process
!     proc_group(NP)       - group id each process assigned to
!
      integer   zone_proc_id(max_zones), proc_zone_id(max_zones),  &
     &          proc_num_zones
      integer, allocatable :: proc_zone_count(:),  &
     &          proc_num_threads(:), proc_group(:)
      double precision, allocatable :: proc_zone_size(:)
!
      integer   myid, root, comm_setup, ierror, dp_type
      integer   num_threads, max_threads, mz_bload, mz_bload_erank
      integer   num_procs, num_procs2
      logical   active
!
! ... Two adjustable parameters for MPI communication
!     max_reqs  -- max. number of async message requests
!     MSG_SIZE  -- optimal message size (in words) for communication
      integer   max_reqs, MSG_SIZE
      parameter (max_reqs=32, MSG_SIZE=400000)
!
      integer   requests(max_reqs), statuses(MPI_STATUS_SIZE,max_reqs)
!
      integer, allocatable :: pcomm_group(:)
      integer(kind2), allocatable :: qcomm_size(:)
      integer(kind2) ::  &
     &          qstart2_west (max_zones), qstart2_east (max_zones),  &
     &          qstart2_south(max_zones), qstart2_north(max_zones)

      end module mpinpb


!---------------------------------------------------------------------
!---------------------------------------------------------------------
!
!  Allocate process-based working arrays
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      subroutine alloc_proc_space

      use mpinpb
      implicit none

      integer ios

      allocate(  &
     &   proc_zone_count (num_procs ),  &
     &   proc_num_threads(num_procs ),  &
     &   proc_group      (num_procs ),  &
     &   proc_zone_size  (num_procs ),  &
     &   pcomm_group     (num_procs2),  &
     &   qcomm_size      (num_procs ), stat=ios )

      if (ios .ne. 0) call error_cond( 1, 'proc arrays' )

      return
      end
