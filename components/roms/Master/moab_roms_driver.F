#include "cppdefs.h"

MODULE moab_roms_driver_mod
!==================================================== John C. Warner ===
!                                                                      !
!  Master library interface to couple ROMS to MPAS-O model using the   !
!  Mesh Oriented DataBase (MOAB) library.                              !
!                                                                      !
!  The following models are coupled to ROMS:                           !
!                                                                      !
#ifdef HAVE_MPASO
!  MPAS-Ocean (MPAS-O) model:                                          !
!       https://mpas-dev.github.io/ocean/ocean.html                    !
!                                                                      !
#endif
!=======================================================================
#include "moab/MOABConfig.h"
   !USE iso_c_binding
   USE iMOAB
   !USE mod_param
#ifdef DISTRIBUTE
   !USE mod_parallel
#endif
   !USE mod_grid
   !USE mod_ocean
   !USE mod_scalars
   !USE mod_iounits
   !USE mod_scalars
   !
   !USE mod_coupler_kinds

   ! Include the ROMS run hooks
   USE ocean_control_mod, ONLY: ROMS_initialize
   USE ocean_control_mod, ONLY: ROMS_run
   USE ocean_control_mod, ONLY: ROMS_finalize

   implicit none
   include 'mpif.h'

   !  TYPE MOABCouplerData
   !     integer :: roms_pid     ! MOAB ROMS App ID
   !     integer :: roms_compid  ! MOAB ROMS Component ID
   !  END TYPE MOABCouplerData

   !  integer   :: ierr                   ! MOAB error code
   !  integer   :: mpi_size, mpi_rank

!
!  Local variable declarations.
!
   logical, save :: first

   integer :: MyColor, MyCOMM, MyError, MyKey, Nnodes
   integer :: MyRank, pelast
   integer :: MPASOcncolor
   integer :: ng, iw, io, ia, ih, icc, roms_exit
   integer :: Nocn_grids, Nmpas_grids
   integer :: N_moabmodels, NnodesOCN
   real(8) :: lcm

   real(8), dimension(:), pointer :: dtocn, dtmpas
   real(4) :: CouplingTime             ! single precision

   integer, dimension(:), pointer :: ocnids
   integer :: OCNid
!
!  Parallel nodes assined to the ocean model.
!
   integer :: peOCN_frst          ! first ocean parallel node
   integer :: peOCN_last          ! last  ocean parallel node
   integer, dimension(:), pointer :: roms_fmcoup ! roms-mpas coupling
   real(8) :: TI_MPAS2OCN           ! MPAS to ROMS ocean coupling interval
   real(8) :: TI_OCN2MPAS           ! ROMS ocean to MPAS coupling interval

   integer, dimension(:,:), pointer :: nOCN2MPAS
   integer, dimension(:,:), pointer :: nMPAS2OCN
   integer, dimension(:,:), pointer :: nOCNFMPAS
   integer, dimension(:,:), pointer :: nMPASFOCN

!
CONTAINS
!
!-----------------------------------------------------------------------
   subroutine errorout(ierr, message)
!-----------------------------------------------------------------------
      implicit none

      integer ierr
      character*(*) message

      if (ierr .ne. 0) then
         print *, message
         call exit(1)
      end if
      return
   end subroutine errorout
!
!
!-----------------------------------------------------------------------
   subroutine driver_initialize()
!-----------------------------------------------------------------------
      implicit none

      real :: gcdlcm_dt
!
!     This is roms exit flag if blows up.
      roms_exit = 0
      Nocn_grids = 1
      Nmpas_grids = 1
      N_moabmodels = 1

      allocate(dtmpas(Nmpas_grids))
      allocate(dtocn(Nocn_grids))
      dtmpas(1) = 1.0
      dtocn(1) = 1.0
      TI_MPAS2OCN = 1.0
      TI_OCN2MPAS = 1.0

      allocate(nOCN2MPAS(Nocn_grids, Nmpas_grids))
      allocate(nMPAS2OCN(Nmpas_grids, Nocn_grids))
      allocate(nOCNFMPAS(Nocn_grids, Nmpas_grids))
      allocate(nMPASFOCN(Nmpas_grids, Nocn_grids))
!
!-----------------------------------------------------------------------
!  Initialize distributed-memory (MPI) configuration
!-----------------------------------------------------------------------
!
!     Initialize MPI execution environment.
!
      CALL mpi_init(MyError)
!
!     Get rank of the local process in the group associated with the
!     communicator.
!
      CALL mpi_comm_size(MPI_COMM_WORLD, Nnodes, MyError)
      CALL mpi_comm_rank(MPI_COMM_WORLD, MyRank, MyError)
      NnodesOCN = Nnodes
!
!  Read in coupled model parameters from standard input.
!
!   CALL read_coawst_par(1)
!
!  Now that we know the input file names and locations for each model,
!  for each model read in the number of grids and the grid time steps.
!
!   CALL read_model_inputs
!
!   CALL allocate_coupler_params
!
#if defined MCT_INTERP_OC2WV || defined MCT_INTERP_OC2AT || \
      defined MCT_INTERP_WV2AT||defined MCT_INTERP_OC2HY
!
!  Read coupled model sparse matrix file names from standard input.
!
      CALL allocate_coupler_iounits
      CALL read_coawst_par(2)
#endif
!
!  Compute the mct send and recv instances.
!
!  For each model grid, determine the number of steps it should
!  compute before it sends data out.
!  For example, nMPAS2OCN(1,2) is the number of steps the MPAS model
!  grid 1 should take before it sends data to the ocn grid 2.
!
      DO iw = 1, Nmpas_grids
         DO io = 1, Nocn_grids
            lcm = gcdlcm_dt(dtmpas(iw), dtocn(io))
            IF (MOD(TI_MPAS2OCN, lcm) .eq. 0) THEN
               nMPAS2OCN(iw, io) = INT(TI_MPAS2OCN/dtmpas(iw))
            ELSE
               lcm = gcdlcm_dt(TI_MPAS2OCN, lcm)
               nMPAS2OCN(iw, io) = INT(lcm/dtmpas(iw))
            END IF
         END DO
      END DO
!
      DO io = 1, Nocn_grids
         DO iw = 1, Nmpas_grids
            lcm = gcdlcm_dt(dtmpas(iw), dtocn(io))
            IF (MOD(TI_OCN2MPAS, lcm) .eq. 0) THEN
               nOCN2MPAS(io, iw) = INT(TI_OCN2MPAS/dtocn(io))
            ELSE
               lcm = gcdlcm_dt(TI_OCN2MPAS, lcm)
               nOCN2MPAS(io, iw) = INT(lcm/dtocn(io))
            END IF
         END DO
      END DO
!
!  Similarly, for each model grid, determine the number of steps
!  it should compute before it recvs data from MPAS-O.
!  For example, nMPASFOCN(1,2) is the number of steps the MPAS model
!  grid 1 should take before it gets data from ocn grid 2.
!
      DO iw = 1, Nmpas_grids
         DO io = 1, Nocn_grids
            lcm = gcdlcm_dt(dtmpas(iw), dtocn(io))
            IF (MOD(TI_OCN2MPAS, lcm) .eq. 0) THEN
               nMPASFOCN(iw, io) = INT(TI_OCN2MPAS/dtmpas(iw))
            ELSE
               lcm = gcdlcm_dt(TI_OCN2MPAS, lcm)
               nMPASFOCN(iw, io) = INT(lcm/dtmpas(iw))
            END IF
         END DO
      END DO
!
      DO io = 1, Nocn_grids
         DO iw = 1, Nmpas_grids
            lcm = gcdlcm_dt(dtmpas(iw), dtocn(io))
            IF (MOD(TI_MPAS2OCN, lcm) .eq. 0) THEN
               nOCNFMPAS(io, iw) = INT(TI_MPAS2OCN/dtocn(io))
            ELSE
               lcm = gcdlcm_dt(TI_MPAS2OCN, lcm)
               nOCNFMPAS(io, iw) = INT(lcm/dtocn(io))
            END IF
         END DO
      END DO
!
!  Allocate coupling variables.
!
      allocate (ocnids(Nocn_grids))
!
      N_moabmodels = 0
      DO ng = 1, Nocn_grids
         N_moabmodels = N_moabmodels + 1
         ocnids(ng) = N_moabmodels
      END DO
!
!  Assign processors to the models.
!
      pelast = -1
      peOCN_frst = pelast + 1
      peOCN_last = peOCN_frst + NnodesOCN - 1
      pelast = peOCN_last
      IF (pelast .ne. Nnodes - 1) THEN
         IF (MyRank .eq. 0) THEN
            WRITE (stdout, 10) pelast + 1, Nnodes
10          FORMAT(/, ' moab_coupler - Number assigned processors: '       &
                                                 &            , i3.3, /, 15x, 'not equal to spawned MPI nodes: ', i3.3)
         END IF
         STOP
      ELSE
         IF (MyRank .eq. 0) THEN
            WRITE (stdout, 19)
19          FORMAT(/, ' Model Coupling: ',/)
            WRITE (stdout, 20) peOCN_frst, peOCN_last
20          FORMAT(/, 7x, 'Ocean Model MPI nodes: ', i3.3, ' - ', i3.3)
!
!  Write out some coupled model info.
!
            DO iw = 1, Nmpas_grids
               DO io = 1, Nocn_grids
                  WRITE (stdout, 25) iw, dtmpas(iw), io, dtocn(io),              &
            &                        TI_MPAS2OCN, nMPAS2OCN(iw, io)
25                FORMAT(/, 7x, 'MPASgrid ', i2.2, ' dt= ', f5.1, ' -to- OCNgrid ', &
                          i2.2, ' dt= ', f5.1, ', CplInt: ', f7.1, ' Steps: ', i3.3)
                  WRITE (stdout, 26) io, dtocn(io), iw, dtmpas(iw),              &
            &                        TI_OCN2MPAS, nOCN2MPAS(io, iw)
26                FORMAT(/, 7x, 'OCNgrid ', i2.2, ' dt= ', f5.1, ' -to- MPASgrid ', &
                          i2.2, ' dt= ', f5.1, ', CplInt: ', f7.1, ' Steps: ', i3.3)
               END DO
            END DO
         END IF
      END IF
!     CALL flush_coawst (stdout)
!
!  Split the communicator into SWAN or WW3, WRF, and ROMS subgroups based
!  on color and key.
!
      MPASOcncolor = 1
      MyKey = 0
      IF ((peOCN_frst .le. MyRank) .and. (MyRank .le. peOCN_last)) THEN
         MyColor = MPASOcncolor
      END IF
      CALL mpi_comm_split(MPI_COMM_WORLD, MyColor, MyKey, MyCOMM,      &
      &                     MyError)
!
!-----------------------------------------------------------------------
!  Run coupled models according to the processor rank.
!-----------------------------------------------------------------------
!
      IF (MyColor .eq. MPASOcncolor) THEN
         first = .TRUE.
         Nrun = 1
         IF (exit_flag .eq. NoError) THEN
            CALL ROMS_initialize(first, MyCOMM)
         END IF
      END IF

      return

   end subroutine driver_initialize

!
!-----------------------------------------------------------------------
   subroutine roms_run_coupled()
!-----------------------------------------------------------------------
!
      implicit none
!
      IF (MyColor .eq. MPASOcncolor) THEN

         IF (exit_flag .eq. NoError) THEN
            run_time = 0.0
            DO ng = 1, Ngrids
               run_time = MAX(run_time, dt(ng)*ntimes(ng))
            END DO

            CALL ROMS_run(run_time)
            roms_exit = exit_flag
         END IF

      END IF

      return

   end subroutine roms_run_coupled

!
!-----------------------------------------------------------------------
   subroutine driver_finalize()
!-----------------------------------------------------------------------
!
      implicit none
!
      IF (MyColor .eq. MPASOcncolor) THEN
         CALL ROMS_finalize
         roms_exit = roms_exit + exit_flag
      END IF

      !!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
      ! Terminates all the mpi - processing and coupling.
      !-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -!
      CALL mpi_barrier(MPI_COMM_WORLD, MyError)
      !CALL MCTWorld_clean()
      CALL mpi_finalize(MyError)

      return

   end subroutine driver_finalize

END MODULE moab_roms_driver_mod

FUNCTION gcdlcm_dt(dtAin, dtBin) result(dtGcdLcm)
!
!=======================================================================
!                                                                      !
!  This function computes the greatest common denominator              !
!  and lowest common multiple.                                         !
!                                                                      !
!  On Input:                                                           !
!     dtA        time step of model A                                  !
!     dtB        time step of model B                                  !
!                                                                      !
!  On Output:                                                          !
!     lcm        least common multiple                                 !
!                                                                      !
!=======================================================================
!
   implicit none
!
!  Imported variable declarations.
!
   real(8), intent(in) :: dtAin, dtBin
   real(8) :: dtGcdLcm
!
!  Local variable declarations.
!
   logical :: stayin
   real(8) :: r, m, n, p, gcd, dtA, dtB, scale
!  Initialize
   scale = 1000.0
!
!-----------------------------------------------------------------------
!  Compute greatest common denominator and least common multiplier.
!-----------------------------------------------------------------------
   dtA = dtAin*scale
   dtB = dtBin*scale
   m = dtA
   n = dtB
   IF (dtA .gt. dtB) THEN
      p = dtA
      dtA = dtB
      dtB = p
   END IF
   stayin = .true.
   DO WHILE (stayin)
      r = mod(dtB, dtA)
      IF (r .eq. 0) THEN
         gcd = dtA
         stayin = .false.
      ELSE
         dtB = dtA
         dtA = r
      END IF
   END DO
   dtGcdLcm = m*n/dtA/scale

   RETURN
END FUNCTION gcdlcm_dt

! -- ROMS is coupled with external MPAS --
! This is an example to show the simplified driver interface
! to couple ROMS with MPAS-Ocean model using the MOAB library.
!PROGRAM master
!    use moab_roms_driver_mod
!    implicit none
!
!    call driver_initialize
!    call roms_run_coupled
!    call driver_finalize
!
!END PROGRAM master

