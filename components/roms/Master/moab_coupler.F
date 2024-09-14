#include "cppdefs.h"
#ifdef HAVE_MOAB
MODULE moab_coupler_mod

#include "moab/MOABConfig.h"
   USE iso_c_binding
   USE iMOAB
   USE mod_param
#ifdef DISTRIBUTE
   USE mod_parallel
#endif
   USE mod_grid
   USE mod_ocean
   USE mod_scalars

   implicit none

   TYPE MOABMeshData
      integer :: roms_pid     ! MOAB ROMS App ID
      integer :: roms_compid  ! MOAB ROMS Component ID
      integer :: global_x, global_y, global_z
      integer :: LBi, UBi, LBj, UBj
      integer :: nghostlayers, Itile, Jtile

      integer, pointer :: vlocal_ids(:, :)
      integer, pointer :: vglobal_ids(:, :)
      integer, pointer :: elocal_ids(:, :)
      integer, pointer :: eglobal_ids(:, :)

   END TYPE MOABMeshData

   integer   :: ierr                   ! MOAB error code
   integer   :: mpi_size, mpi_rank

   integer   :: numgridlevels
   type(MOABMeshData), pointer :: griddata(:)

!
CONTAINS
!
!-----------------------------------------------------------------------
   SUBROUTINE errorout(ierr, message)
!-----------------------------------------------------------------------
      implicit none

      integer ierr
      character*(*) message

      if (ierr .ne. 0) then
         print *, message
         call exit(1)
      end if
      return
   end

!-----------------------------------------------------------------------
   SUBROUTINE MOAB_initialize() !(Iarg, Carg)
!-----------------------------------------------------------------------
!
      implicit none
!
!  Imported variable declarations.
!
!      integer, intent(in) :: Iarg
!      character (len=*), intent(inout) :: Carg
!
!  Local variable declarations.
!
      ierr = iMOAB_Initialize()
      CALL errorout(ierr, 'failed to initialize iMOAB')

      RETURN
   END SUBROUTINE MOAB_initialize

!-----------------------------------------------------------------------
   SUBROUTINE MOAB_finalize()
!-----------------------------------------------------------------------
!
      ierr = iMOAB_Finalize()
      CALL errorout(ierr, 'failed to finalize iMOAB')

      RETURN
   END SUBROUTINE MOAB_finalize

!-----------------------------------------------------------------------
   SUBROUTINE MOAB_RegisterROMSApplications(Ngrids)
!-----------------------------------------------------------------------
!
      implicit none
!
!  Local variable declarations.
!
      integer, intent(in) :: Ngrids
      character(len=8)  :: fmt             ! format descriptor
      character(len=4) :: fmtng    ! MOAB ROMS App name
      character(len=32) :: roms_appname    ! MOAB ROMS App name
      type(MOABMeshData) :: currentgrid
      integer :: ng

      numgridlevels = Ngrids
      fmt = '(I2.2)' ! an integer of width 5 with zeros at the left

      allocate (griddata(numgridlevels))

      DO ng = 1, numgridlevels
         write (fmtng, fmt) ng ! converting integer to string using a 'internal file'

         roms_appname = 'ROMS'//trim(fmtng)//C_NULL_CHAR

         ! specify a unique external id;
         griddata(ng)%roms_compid = 225 + ng

#ifdef DISTRIBUTE
#ifdef MOAB_HAVE_MPI
         ierr = iMOAB_RegisterApplication(roms_appname, OCN_COMM_WORLD, &
                                          griddata(ng)%roms_compid, &
                                          griddata(ng)%roms_pid)
         CALL errorout(ierr, 'failed to register application')
#else
         ierr = iMOAB_RegisterApplication(roms_appname, &
                                          griddata(ng)%roms_compid, &
                                          griddata(ng)%roms_pid)
         CALL errorout(ierr, 'failed to register application')
#endif
#else
         ierr = iMOAB_RegisterApplication(roms_appname, &
                                          griddata(ng)%roms_compid, &
                                          griddata(ng)%roms_pid)
         CALL errorout(ierr, 'failed to register application')
#endif

      END DO

      ! set defaults
      mpi_size = 1
      mpi_rank = 0

#ifdef DISTRIBUTE
#ifdef MOAB_HAVE_MPI
      call MPI_Comm_Size(OCN_COMM_WORLD, mpi_size, ierr)
      CALL errorout(ierr, 'failed to get size')
      call MPI_Comm_Rank(OCN_COMM_WORLD, mpi_rank, ierr)
      CALL errorout(ierr, 'failed to get rank')
#endif
#endif

      RETURN
   END SUBROUTINE MOAB_RegisterROMSApplications

!-----------------------------------------------------------------------
   SUBROUTINE MOAB_DeregisterROMSApplications()
!-----------------------------------------------------------------------
!
      integer :: ng
      DO ng = 1, numgridlevels
         ! destroy the application instance
         ierr = iMOAB_DeregisterApplication(griddata(ng)%roms_pid)
         CALL errorout(ierr, 'failed to de-register application')
      END DO

      RETURN
   END SUBROUTINE MOAB_DeregisterROMSApplications

!-----------------------------------------------------------------------
   SUBROUTINE MOAB_get_bounds(ng, tile, gtype)
!-----------------------------------------------------------------------
!
! Subroutine copied from get_bounds.F
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine compute grid bounds in the I- and J-directions.        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     tile       Domain partition.                                     !
!     gtype      C-grid type. If zero, compute array allocation bounds.!
!                  Otherwise, compute bounds for IO processing.        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Itile      I-tile coordinate (a value from 0 to NtileI(ng)).     !
!     Jtile      J-tile coordinate (a value from 0 to NtileJ(ng)).     !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
#ifdef NESTING
      USE mod_nesting
#endif
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, gtype
      integer :: mItile, mJtile, mLBi, mUBi, mLBj, mUBj
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax, Nghost
#ifdef DISTRIBUTE
      integer :: Iend, Istr, Jend, Jstr
      integer :: IstrM, IstrR, IstrU, IendR
      integer :: JstrM, JstrR, JstrV, JendR
      integer :: IstrB, IendB, IstrP, IendP, IstrT, IendT
      integer :: JstrB, JendB, JstrP, JendP, JstrT, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
      integer :: MyType
!
!-----------------------------------------------------------------------
!  Set array bounds in the I- and J-direction for distributed-memory
!  configurations.
!-----------------------------------------------------------------------
!
!  Set first and last grid-points according to staggered C-grid
!  classification.  If gtype = 0, it returns the values needed for
!  array allocation. Otherwise, it returns the values needed for IO
!  processing.
!
      Nghost = griddata(ng)%nghostlayers
      MyType = ABS(gtype)
      IF (MyType .eq. 0) THEN
# ifdef NESTING
         IF (RefinedGrid(ng) .and. (RefineScale(ng) .gt. 0)) THEN
            Imin = -NghostPoints
            Imax = Lm(ng) + NghostPoints
            Jmin = -NghostPoints
            Jmax = Mm(ng) + NghostPoints
         ELSE
            IF (CompositeGrid(iwest, ng) .or. EWperiodic(ng)) THEN
               Imin = -NghostPoints
            ELSE
               Imin = 0
            END IF
            IF (CompositeGrid(ieast, ng) .or. EWperiodic(ng)) THEN
               Imax = Lm(ng) + NghostPoints
            ELSE
               Imax = Im(ng) + 1
            END IF
            IF (CompositeGrid(isouth, ng) .or. NSperiodic(ng)) THEN
               Jmin = -NghostPoints
            ELSE
               Jmin = 0
            END IF
            IF (CompositeGrid(inorth, ng) .or. NSperiodic(ng)) THEN
               Jmax = Mm(ng) + NghostPoints
            ELSE
               Jmax = Jm(ng) + 1
            END IF
         END IF
# else
         IF (EWperiodic(ng)) THEN
            IF (NSperiodic(ng)) THEN
               Imin = -NghostPoints
               Imax = Im(ng) + NghostPoints
               Jmin = -NghostPoints
               Jmax = Jm(ng) + NghostPoints
            ELSE
               Imin = -NghostPoints
               Imax = Im(ng) + NghostPoints
               Jmin = 0
               Jmax = Jm(ng) + 1
            END IF
         ELSE
            IF (NSperiodic(ng)) THEN
               Imin = 0
               Imax = Im(ng) + 1
               Jmin = -NghostPoints
               Jmax = Jm(ng) + NghostPoints
            ELSE
               Imin = 0
               Imax = Im(ng) + 1
               Jmin = 0
               Jmax = Jm(ng) + 1
            END IF
         END IF
# endif
      ELSE
         IF ((MyType .eq. p2dvar) .or. (MyType .eq. u2dvar) .or.                &
      &      (MyType .eq. p3dvar) .or. (MyType .eq. u3dvar)) THEN
            Imin = 1
         ELSE
            Imin = 0
         END IF
         Imax = Lm(ng) + 1
         IF ((MyType .eq. p2dvar) .or. (MyType .eq. v2dvar) .or.                &
      &      (MyType .eq. p3dvar) .or. (MyType .eq. v3dvar)) THEN
            Jmin = 1
         ELSE
            Jmin = 0
         END IF
         Jmax = Mm(ng) + 1
      END IF
!
!  Set physical, overlapping (Nghost>0) or non-overlapping (Nghost=0)
!  grid bounds according to tile rank.
!
      CALL get_tile(ng, tile,                                          &
     &               mItile, mJtile,                                      &
     &               Istr, Iend, Jstr, Jend,                            &
     &               IstrM, IstrR, IstrU, IendR,                        &
     &               JstrM, JstrR, JstrV, JendR,                        &
     &               IstrB, IendB, IstrP, IendP, IstrT, IendT,          &
     &               JstrB, JendB, JstrP, JendP, JstrT, JendT,          &
     &               Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1,          &
     &               Iendp1, Iendp2, Iendp2i, Iendp3,                   &
     &               Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1,          &
     &               Jendp1, Jendp2, Jendp2i, Jendp3)
!
      IF ((mItile .eq. -1) .or. (mItile .eq. 0)) THEN
         mLBi = Imin
      ELSE
         mLBi = Istr - Nghost - 1
      END IF
      IF ((mItile .eq. -1) .or. (mItile .eq. (NtileI(ng) - 1))) THEN
         mUBi = Imax - 1
      ELSE
         mUBi = Iend + Nghost
      END IF
      IF ((mJtile .eq. -1) .or. (mJtile .eq. 0)) THEN
         mLBj = Jmin
      ELSE
         mLBj = Jstr - Nghost - 1
      END IF
      IF ((mJtile .eq. -1) .or. (mJtile .eq. (NtileJ(ng) - 1))) THEN
         mUBj = Jmax - 1
      ELSE
         mUBj = Jend + Nghost
      END IF
#else
!
!-----------------------------------------------------------------------
!  Set array allocation bounds in the I- and J-direction for serial and
!  shared-memory configurations.
!-----------------------------------------------------------------------
!
      IF (tile .eq. -1) THEN
         mItile(ng) = -1
         mJtile(ng) = -1
      ELSE
         mJtile(ng) = tile/NtileI(ng)
         mItile(ng) = tile - mJtile*NtileI(ng)
      END IF
!
# ifdef NESTING
      IF (RefinedGrid(ng) .and. (RefineScale(ng) .gt. 0)) THEN
         mLBi = -NghostPoints
         mUBi = Lm(ng) + NghostPoints
         mLBj = -NghostPoints
         mUBj = Mm(ng) + NghostPoints
      ELSE
         IF (CompositeGrid(iwest, ng) .or. EWperiodic(ng)) THEN
            mLBi = -NghostPoints
         ELSE
            mLBi = 0
         END IF
         IF (CompositeGrid(ieast, ng) .or. EWperiodic(ng)) THEN
            mUBi = Lm(ng) + NghostPoints
         ELSE
            mUBi = Im(ng) + 1
         END IF
         IF (CompositeGrid(isouth, ng) .or. NSperiodic(ng)) THEN
            mLBj = -NghostPoints
         ELSE
            mLBj = 0
         END IF
         IF (CompositeGrid(inorth, ng) .or. NSperiodic(ng)) THEN
            mUBj = Mm(ng) + NghostPoints
         ELSE
            mUBj = Jm(ng) + 1
         END IF
      END IF
# else

      IF (EWperiodic(ng)) THEN
         IF (NSperiodic(ng)) THEN
            mLBi = -NghostPoints
            mUBi = Im(ng) + NghostPoints
            mLBj = -NghostPoints
            mUBj = Jm(ng) + NghostPoints
         ELSE
            mLBi = -NghostPoints
            mUBi = Im(ng) + NghostPoints
            mLBj = 0
            mUBj = Jm(ng) + 1
         END IF
      ELSE
         IF (NSperiodic(ng)) THEN
            mLBi = 0
            mUBi = Im(ng) + 1
            mLBj = -NghostPoints
            mUBj = Jm(ng) + NghostPoints
         ELSE
            mLBi = 0
            mUBi = Im(ng) + 1
            mLBj = 0
            mUBj = Jm(ng) + 1
         END IF
      END IF
# endif
#endif

      ! if running on a single process, use standard bounds
      if (mpi_size == 1) then
         mLBi = BOUNDS(ng)%Istr(tile)
         mUBi = BOUNDS(ng)%Iend(tile)
         mLBj = BOUNDS(ng)%Jstr(tile)
         mUBj = BOUNDS(ng)%Jend(tile)
      end if

      griddata(ng)%Itile = mItile
      griddata(ng)%Jtile = mJtile
      griddata(ng)%LBi = mLBi
      griddata(ng)%UBi = mUBi
      griddata(ng)%LBj = mLBj
      griddata(ng)%UBj = mUBj

      RETURN
   END SUBROUTINE MOAB_get_bounds

!-----------------------------------------------------------------------
   SUBROUTINE MOAB_AllocateGridData(ng, tile)
!-----------------------------------------------------------------------

      USE mod_param
      USE mod_coupling
      USE mod_grid
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile

!     Nghost     Number of ghost-points in the halo region:            !
!                  Nghost = 0,  compute non-overlapping bounds.        !
!                  Nghost > 0,  compute overlapping bounds.            !
      griddata(ng)%nghostlayers = 0
      CALL MOAB_get_bounds(ng, tile, 0)

      griddata(ng)%global_x = Lm(ng)+3
      griddata(ng)%global_y = Mm(ng)+3
      griddata(ng)%global_z = N(ng)

      !print *, "Global Lm, Mm: ", Lm(ng), Mm(ng)
      print *, "Global data: ", griddata(ng)%global_x, &
         griddata(ng)%global_y, griddata(ng)%global_z
      print *, "Detailed info: ", griddata(ng)%LBi, griddata(ng)%UBi, &
         griddata(ng)%LBj, griddata(ng)%UBj

      RETURN
   END SUBROUTINE MOAB_AllocateGridData

!-----------------------------------------------------------------------
   SUBROUTINE MOAB_DeallocateGridData(ng)
!-----------------------------------------------------------------------
!  Imported variable declarations.
!
      integer, intent(in) :: ng

      ! de-allocate local data for grid
      deallocate (griddata(ng)%vlocal_ids)
      deallocate (griddata(ng)%vglobal_ids)
      deallocate (griddata(ng)%elocal_ids)
      deallocate (griddata(ng)%eglobal_ids)

      RETURN
   END SUBROUTINE MOAB_DeallocateGridData

!-----------------------------------------------------------------------
   SUBROUTINE MOAB_SetGrid(ng, tile)
!-----------------------------------------------------------------------

      USE mod_param
      USE mod_scalars
      ! USE mod_coupling
      USE mod_grid
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer :: mLBi, mUBi, mLBj, mUBj

      ! local data
      integer :: ierr, num_verts_in_cells, num_elements
      integer :: i, j, i1, offset
      real(r8), pointer :: moab_vert_coords(:)
      integer, pointer :: connectivity(:), connectivityg(:)
      character*100 tagname, lnum, filename, fileopts
      integer ndimcoord, ndimension, nmax_vpere, element_type
      integer tagtype, numco, tag_sto_len, ent_type, tagindex, block_id
      real(r8), pointer :: datasets(:, :), datasets3d(:, :, :)
      real(r8) :: factor
!
      ! allocate the z-levels for all tiles
      ! call MOAB_AllocateGridData()
      ! call MOAB_AllocateZLayers()
      ! call MOAB_SetZLayers (ng)

! #ifdef SOLVE3D
!       ! copy the z-levels for all tiles
!       zlevels = GRID(ng)%Hz(LBi:UBi, LBj:UBj, :)
! #endif

      ndimension = 3 ! dimension of coordinates
      element_type = 3 ! MBQUAD
      nmax_vpere = 4 ! MBQUAD

      mLBi = griddata(ng)%LBi-1
      mUBi = griddata(ng)%UBi+2
      mLBj = griddata(ng)%LBj-1
      mUBj = griddata(ng)%UBj+2
      !print *, "MOAB_SetGrid info: ", mLBi, mUBi, &
      !   mLBj, mUBj, ng, tile, griddata(ng)%global_x, griddata(ng)%global_y

      ! Eg.,
      ! LBi, UBi = [1, M]
      ! LBj, UBj = [1, N]
      ! num_verts_in_cells = M*N
      ! num_elems_in_cells = (M-1)*(N-1)
      num_verts_in_cells = (mUBi - mLBi + 1)*(mUBj - mLBj + 1)
      ndimcoord = 3*num_verts_in_cells
      allocate (moab_vert_coords(ndimcoord))
      ! allocate vertex ID arrays
      allocate (griddata(ng)%vlocal_ids(mLBi:mUBi, mLBj:mUBj))
      allocate (griddata(ng)%vglobal_ids(mLBi:mUBi, mLBj:mUBj))

      i1 = 1
      do j = mLBj, mUBj
         do i = mLBi, mUBi
            griddata(ng)%vlocal_ids(i, j) = i1
            griddata(ng)%vglobal_ids(i, j) = 1 + i + j*griddata(ng)%global_x
            !     globalIds(i,j) = j + (i-1) * griddata(ng)%global_y

#ifdef ANA_GRID
               factor = 1E3
               ! moab_vert_coords(3*i1 - 2) = GRID(ng)%xp(i, j)/1E3
               ! moab_vert_coords(3*i1 - 1) = GRID(ng)%yp(i, j)/1E3
               ! moab_vert_coords(3*i1)     = GRID(ng)%z_r(i, j, griddata(ng)%global_z)

               if (i .eq. mLBi .or. i .eq. mUBi .or. j .eq. mLBj .or. j .eq. mUBj) then
               ! Now set the coordinates for halo regions
                  if (i .eq. mLBi .and. j .eq. mLBj) then
                     moab_vert_coords(3*i1 - 2) = (2*GRID(ng)%xp(i+1, j+1)-GRID(ng)%xp(i+2, j+1))/factor
                     moab_vert_coords(3*i1 - 1) = (2*GRID(ng)%yp(i+1, j+1)-GRID(ng)%yp(i+1, j+2))/factor
                     moab_vert_coords(3*i1)     = GRID(ng)%z_r(i+1, j+1, griddata(ng)%global_z)
                  else if (i .eq. mLBi .and. j .eq. mUBj) then
                     moab_vert_coords(3*i1 - 2) = (2*GRID(ng)%xp(i+1, j-1)-GRID(ng)%xp(i+2, j-1))/factor
                     moab_vert_coords(3*i1 - 1) = (2*GRID(ng)%yp(i+1, j-1)-GRID(ng)%yp(i+1, j-2))/factor
                     moab_vert_coords(3*i1)     = GRID(ng)%z_r(i+1, j-1, griddata(ng)%global_z)
                  else if (i .eq. mUBi .and. j .eq. mLBj) then
                     ! print *, "Left corner: ", i, j, i1, GRID(ng)%xp(i-1, j+1),GRID(ng)%xp(i-2, j+1), GRID(ng)%yp(i-1, j+1),GRID(ng)%yp(i-1, j+2)
                     moab_vert_coords(3*i1 - 2) = (2*GRID(ng)%xp(i-1, j+1)-GRID(ng)%xp(i-2, j+1))/factor
                     moab_vert_coords(3*i1 - 1) = (2*GRID(ng)%yp(i-1, j+1)-GRID(ng)%yp(i-1, j+2))/factor
                     moab_vert_coords(3*i1)     = GRID(ng)%z_r(i-1, j+1, griddata(ng)%global_z)
                  else if (i .eq. mUBi .and. j .eq. mUBj) then
                     ! print *, "Right corner: ", i, j, i1, GRID(ng)%xp(i-1, j-1), GRID(ng)%xp(i-2, j-1), GRID(ng)%yp(i-1, j-1), GRID(ng)%yp(i-1, j-2)
                     moab_vert_coords(3*i1 - 2) = (2*GRID(ng)%xp(i-1, j-1)-GRID(ng)%xp(i-2, j-1))/factor
                     moab_vert_coords(3*i1 - 1) = (2*GRID(ng)%yp(i-1, j-1)-GRID(ng)%yp(i-1, j-2))/factor
                     moab_vert_coords(3*i1)     = GRID(ng)%z_r(i-1, j-1, griddata(ng)%global_z)
                  else if (i .eq. mLBi) then
                     moab_vert_coords(3*i1 - 2) = (2*GRID(ng)%xp(i+1, j)-GRID(ng)%xp(i+2, j))/factor
                     moab_vert_coords(3*i1 - 1) = GRID(ng)%yp(i+1, j)/factor
                     moab_vert_coords(3*i1)     = GRID(ng)%z_r(i+1, j, griddata(ng)%global_z)
                  else if (i .eq. mUBi) then
                     moab_vert_coords(3*i1 - 2) = (2*GRID(ng)%xp(i-1, j)-GRID(ng)%xp(i-2, j))/factor
                     moab_vert_coords(3*i1 - 1) = GRID(ng)%yp(i-1, j)/factor
                     moab_vert_coords(3*i1)     = GRID(ng)%z_r(i-1, j, griddata(ng)%global_z)
                  else if (j .eq. mLBj) then
                     moab_vert_coords(3*i1 - 2) = GRID(ng)%xp(i, j+1)/factor
                     moab_vert_coords(3*i1 - 1) = (2*GRID(ng)%yp(i, j+1)-GRID(ng)%yp(i, j+2))/factor
                     moab_vert_coords(3*i1)     = GRID(ng)%z_r(i, j+1, griddata(ng)%global_z)
                  else if (j .eq. mUBj) then
                     moab_vert_coords(3*i1 - 2) = GRID(ng)%xp(i, j-1)/factor
                     moab_vert_coords(3*i1 - 1) = (2*GRID(ng)%yp(i, j-1)-GRID(ng)%yp(i, j-2))/factor
                     moab_vert_coords(3*i1)     = GRID(ng)%z_r(i, j-1, griddata(ng)%global_z)
                  end if
               else
                  moab_vert_coords(3*i1 - 2) = GRID(ng)%xp(i, j)/factor
                  moab_vert_coords(3*i1 - 1) = GRID(ng)%yp(i, j)/factor
                  moab_vert_coords(3*i1)     = GRID(ng)%z_r(i, j, griddata(ng)%global_z)
               end if
#else
               if (i .eq. mLBi .or. i .eq. mUBi .or. j .eq. mLBj .or. j .eq. mUBj) then
               ! Now set the coordinates for halo regions
                  if (i .eq. mLBi .and. j .eq. mLBj) then ! Lower Left corner
                     moab_vert_coords(3*i1 - 2) = 2*GRID(ng)%lonp(i+1, j+1)-GRID(ng)%lonp(i+2, j+1)
                     moab_vert_coords(3*i1 - 1) = 2*GRID(ng)%latp(i+1, j+1)-GRID(ng)%latp(i+1, j+2)
                     moab_vert_coords(3*i1)     = GRID(ng)%z_r(i+1, j+1, griddata(ng)%global_z)
                  else if (i .eq. mUBi .and. j .eq. mLBj) then ! Lower Right corner
                     ! print *, "Lower Right corner: ", i, j, i1, GRID(ng)%lonp(i-1, j+1),GRID(ng)%lonp(i-2, j+1), GRID(ng)%latp(i-1, j+1),GRID(ng)%latp(i-1, j+2)
                     moab_vert_coords(3*i1 - 2) = 2*GRID(ng)%lonp(i-1, j+1)-GRID(ng)%lonp(i-2, j+1)
                     moab_vert_coords(3*i1 - 1) = 2*GRID(ng)%latp(i-1, j+1)-GRID(ng)%latp(i-1, j+2)
                     moab_vert_coords(3*i1)     = GRID(ng)%z_r(i-1, j+1, griddata(ng)%global_z)
                  else if (i .eq. mUBi .and. j .eq. mUBj) then ! Upper Left corner
                     ! print *, "Lower Left corner: ", i, j, i1, GRID(ng)%lonp(i-1, j-1), GRID(ng)%lonp(i-2, j-1), GRID(ng)%latp(i-1, j-1), GRID(ng)%latp(i-1, j-2)
                     moab_vert_coords(3*i1 - 2) = 2*GRID(ng)%lonp(i-1, j-1)-GRID(ng)%lonp(i-2, j-1)
                     moab_vert_coords(3*i1 - 1) = 2*GRID(ng)%latp(i-1, j-1)-GRID(ng)%latp(i-1, j-2)
                     moab_vert_coords(3*i1)     = GRID(ng)%z_r(i-1, j-1, griddata(ng)%global_z)
                  else if (i .eq. mLBi .and. j .eq. mUBj) then ! Upper Right corner
                     moab_vert_coords(3*i1 - 2) = 2*GRID(ng)%lonp(i+1, j-1)-GRID(ng)%lonp(i+2, j-1)
                     moab_vert_coords(3*i1 - 1) = 2*GRID(ng)%latp(i+1, j-1)-GRID(ng)%latp(i+1, j-2)
                     moab_vert_coords(3*i1)     = GRID(ng)%z_r(i+1, j-1, griddata(ng)%global_z)
                  else if (i .eq. mLBi) then
                     moab_vert_coords(3*i1 - 2) = 2*GRID(ng)%lonp(i+1, j)-GRID(ng)%lonp(i+2, j)
                     moab_vert_coords(3*i1 - 1) = GRID(ng)%latp(i+1, j)
                     moab_vert_coords(3*i1)     = GRID(ng)%z_r(i+1, j, griddata(ng)%global_z)
                  else if (i .eq. mUBi) then
                     moab_vert_coords(3*i1 - 2) = 2*GRID(ng)%lonp(i-1, j)-GRID(ng)%lonp(i-2, j)
                     moab_vert_coords(3*i1 - 1) = GRID(ng)%latp(i-1, j)
                     moab_vert_coords(3*i1)     = GRID(ng)%z_r(i-1, j, griddata(ng)%global_z)
                  else if (j .eq. mLBj) then
                     moab_vert_coords(3*i1 - 2) = GRID(ng)%lonp(i, j+1)
                     moab_vert_coords(3*i1 - 1) = 2*GRID(ng)%latp(i, j+1)-GRID(ng)%latp(i, j+2)
                     moab_vert_coords(3*i1)     = GRID(ng)%z_r(i, j+1, griddata(ng)%global_z)
                  else if (j .eq. mUBj) then
                     moab_vert_coords(3*i1 - 2) = GRID(ng)%lonp(i, j-1)
                     moab_vert_coords(3*i1 - 1) = 2*GRID(ng)%latp(i, j-1)-GRID(ng)%latp(i, j-2)
                     moab_vert_coords(3*i1)     = GRID(ng)%z_r(i, j-1, griddata(ng)%global_z)
                  end if
               else
                  moab_vert_coords(3*i1 - 2) = GRID(ng)%lonp(i, j) !xp(i,j)
                  moab_vert_coords(3*i1 - 1) = GRID(ng)%latp(i, j) !yp(i,j)
                  moab_vert_coords(3*i1)     = GRID(ng)%z_r(i, j, griddata(ng)%global_z)
               end if
#endif

            i1 = i1 + 1
         end do
      end do

      print *, "Rank: ", mpi_rank, " Total vertices: ", num_verts_in_cells, griddata(ng)%vglobal_ids(mUBi, mUBj)

      ! now let us create the vertices in MOAB
      ierr = iMOAB_CreateVertices(griddata(ng)%roms_pid, ndimcoord, ndimension, moab_vert_coords)
      call errorout(ierr, 'failed to create vertices')

      ! Eg.,
      ! UBi-LBi = 4
      ! UBj-LBj = 4
      ! num_elements = 16
      !
      !     5   10  15   20   25
      !       D   H    L    P
      !     4   9   14   19   24
      !       C   G    K    O
      !     3   8   13   18   23
      !       B   F    J    N
      !     2   7   12   17   22
      !       A   E    I    M
      !     1   6   11   16   21
      !
      ! eloffset = mpi_rank * num_elements
      num_elements = (mUBi - mLBi)*(mUBj - mLBj)
      allocate (connectivity(num_elements*4))
      allocate (connectivityg(num_elements*4))
      ! allocate element ID arrays
      allocate (griddata(ng)%elocal_ids(mLBi:mUBi - 1, mLBj:mUBj - 1))
      allocate (griddata(ng)%eglobal_ids(mLBi:mUBi - 1, mLBj:mUBj - 1))

      offset = 1
      i1 = 1
      do j = mLBj, mUBj - 1
         do i = mLBi, mUBi - 1
            griddata(ng)%elocal_ids(i, j) = i1
            griddata(ng)%eglobal_ids(i, j) = i + j*(griddata(ng)%global_x-1) + 1
            connectivity(offset) = griddata(ng)%vlocal_ids(i, j)
            connectivity(offset + 1) = griddata(ng)%vlocal_ids(i + 1, j)
            connectivity(offset + 2) = griddata(ng)%vlocal_ids(i + 1, j + 1)
            connectivity(offset + 3) = griddata(ng)%vlocal_ids(i, j + 1)
#ifdef MOAB_DEBUG
            connectivityg(offset) = griddata(ng)%vglobal_ids(i, j)
            connectivityg(offset + 1) = griddata(ng)%vglobal_ids(i + 1, j)
            connectivityg(offset + 2) = griddata(ng)%vglobal_ids(i + 1, j + 1)
            connectivityg(offset + 3) = griddata(ng)%vglobal_ids(i, j + 1)
#endif
            offset = offset + 4
            i1 = i1 + 1
         end do
      end do
#ifdef MOAB_DEBUG
      print *, "Rank: ", mpi_rank, " Consistency check for elements: ", num_elements, i1, OCN_COMM_WORLD
      ! call MPI_Barrier(OCN_COMM_WORLD, ierr)
      ! if ( mpi_rank .eq. 0 ) then
      !   print *, "Now printing local indices"
      !   offset = 1
      !   do j = 1, num_elements
      !     print *, "[0]: ", j, connectivity(offset), connectivity(offset+1), connectivity(offset+2), connectivity(offset+3)
      !     offset = offset + 4
      !   enddo
      ! endif
      ! call MPI_Barrier(OCN_COMM_WORLD, ierr)
      ! if ( mpi_rank .eq. 1 ) then
      !   offset = 1
      !   do j = 1, num_elements
      !     print *, "[1]: ", j, connectivity(offset), connectivity(offset+1), connectivity(offset+2), connectivity(offset+3)
      !     offset = offset + 4
      !   enddo
      ! endif
      ! call MPI_Barrier(OCN_COMM_WORLD, ierr)

      if (mpi_rank .eq. 0) then
         print *, "Now printing global indices"
         offset = 1
         do j = 1, num_elements
            print *, "[0]: ", j, connectivityg(offset), connectivityg(offset + 1), &
               connectivityg(offset + 2), connectivityg(offset + 3)
            offset = offset + 4
         end do
         print *, "  "
         print *, "  "
      end if
      call MPI_Barrier(OCN_COMM_WORLD, ierr)
      if (mpi_rank .eq. 1) then
         offset = 1
         do j = 1, num_elements
            print *, "[1]: ", j, connectivityg(offset), connectivityg(offset + 1), &
               connectivityg(offset + 2), connectivityg(offset + 3)
            offset = offset + 4
         end do
         print *, "  "
      end if
      call MPI_Barrier(OCN_COMM_WORLD, ierr)
#endif

      print *, "Rank: ", mpi_rank, " Total elements: ", num_elements, griddata(ng)%eglobal_ids(mUBi - 1, mUBj - 1)

      block_id = 100*griddata(ng)%roms_compid + mpi_rank ! we should have only one block right now

      ! next let us create the QUAD elements in MOAB with the right connectivity
      ierr = iMOAB_CreateElements(griddata(ng)%roms_pid, num_elements, element_type, nmax_vpere, connectivity, block_id);
      call errorout(ierr, 'failed to create quads')

      ! set the global id for vertices
      ! first, retrieve the tag
      tagname = 'GLOBAL_ID'//C_NULL_CHAR
      tagtype = 0  ! dense, integer
      numco = 1
      ierr = iMOAB_DefineTagStorage(griddata(ng)%roms_pid, tagname, tagtype, numco, tagindex)
      call errorout(ierr, 'failed to get global id tag')

      ! now set the values
      ent_type = 0 ! vertex type
      ierr = iMOAB_SetIntTagStorage(griddata(ng)%roms_pid, tagname, num_verts_in_cells, ent_type, &
                                    griddata(ng)%vglobal_ids)
      call errorout(ierr, 'failed to set global id tag for vertices')

      ! set global id tag for elements
      ent_type = 1 ! now set the global id tag on elements
      ierr = iMOAB_SetIntTagStorage(griddata(ng)%roms_pid, tagname, num_elements, ent_type, &
                                    griddata(ng)%eglobal_ids)
      call errorout(ierr, 'failed to set global id tag for quads')

#ifdef MOAB_HAVE_MPI
      if (mpi_size > 1) then
         ! get next block
         ierr = iMOAB_ResolveSharedEntities(griddata(ng)%roms_pid, num_verts_in_cells, &
                                            griddata(ng)%vglobal_ids)
         call errorout(ierr, 'failed to resolve shared entities')

         !   ierr = iMOAB_MergeVertices( griddata(ng)%roms_pid )
         !   call errorout(ierr, 'failed to parallel merge mesh')

         !   ierr = iMOAB_DetermineGhostEntities ( griddata(ng)%roms_pid, 2, 2, 1 )
         !   call errorout(ierr, 'failed to create ghost entities')
      end if
#endif

      allocate (datasets(mLBi:mUBi-1, mLBj:mUBj-1))
      allocate (datasets3d(mLBi:mUBi-1, mLBj:mUBj-1, 1:griddata(ng)%global_z))

      ! create a tag for bathymetry data and set in a tag
      tagname = 'bathymetry'//C_NULL_CHAR
      tagtype = 1  ! dense, double
      numco = 1
      ! first, create the tag
      ierr = iMOAB_DefineTagStorage(griddata(ng)%roms_pid, tagname, tagtype, numco, tagindex)
      call errorout(ierr, 'failed to get bathymetry tag')

      ! now set the values
      ent_type = 1 ! element type

      datasets = GRID(ng)%h(mLBi:mUBi-1, mLBj:mUBj-1)
      ierr = iMOAB_SetDoubleTagStorage(griddata(ng)%roms_pid, tagname, num_elements, ent_type, &
                                       datasets)
      call errorout(ierr, 'failed to set bathymetry tag on elements')

      ! create a tag for zlevels data and set in a tag
      tagname = 'coriolis'//C_NULL_CHAR
      tagtype = 1  ! dense, double
      numco = 1
      ! first, create the tag
      ierr = iMOAB_DefineTagStorage(griddata(ng)%roms_pid, tagname, tagtype, numco, tagindex)
      call errorout(ierr, 'failed to get coriolis tag')

      ! now set the values
      datasets = GRID(ng)%f(mLBi:mUBi-1, mLBj:mUBj-1)
      ierr = iMOAB_SetDoubleTagStorage(griddata(ng)%roms_pid, tagname, num_elements, ent_type, &
                                       datasets)
      call errorout(ierr, 'failed to set coriolis tag on elements')

#ifdef MASKING

      ! create a tag for zlevels data and set in a tag
      tagname = 'GRID_IMASK'//C_NULL_CHAR
      ! first, create the tag
      ierr = iMOAB_DefineTagStorage(griddata(ng)%roms_pid, tagname, tagtype, numco, tagindex)
      call errorout(ierr, 'failed to get GRID_IMASK tag')

      ! now set the values
      datasets = GRID(ng)%rmask(mLBi:mUBi - 1, mLBj:mUBj - 1)
      ierr = iMOAB_SetDoubleTagStorage(griddata(ng)%roms_pid, tagname, num_elements, ent_type, &
                                       datasets)
      call errorout(ierr, 'failed to set GRID_IMASK tag on elements')

#endif

      ! ! set the 3d column depth data for each element
      ! tagname = 'depth'//C_NULL_CHAR
      ! tagtype = 1  ! dense, double
      ! numco = griddata(ng)%global_z
      ! ! first, create the tag
      ! ierr = iMOAB_DefineTagStorage(griddata(ng)%roms_pid, tagname, tagtype, numco, tagindex)
      ! call errorout(ierr, 'failed to get depth tag')

      ! ent_type = 0 ! vertex type
      ! datasets3d = GRID(ng)%z_r(mLBi:mUBi, mLBj:mUBj, 1:griddata(ng)%global_z)
      ! ierr = iMOAB_SetDoubleTagStorage(griddata(ng)%roms_pid, tagname, num_verts_in_cells*numco, ent_type, &
      !                                  datasets3d)
      ! call errorout(ierr, 'failed to set depth tag on vertices')

      deallocate (datasets)
      deallocate (datasets3d)
      deallocate (moab_vert_coords)
      deallocate (connectivity)

      ierr = iMOAB_UpdateMeshInfo(griddata(ng)%roms_pid)
      call errorout(ierr, 'failed to update mesh info')

      filename = 'roms_2d_mesh.h5m'//C_NULL_CHAR
      fileopts = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
      ierr = iMOAB_WriteMesh(griddata(ng)%roms_pid, filename, fileopts)
      call errorout(ierr, 'failed to write mesh to disk')

      RETURN
   END SUBROUTINE MOAB_SetGrid

#ifdef SOLVE3D
!-----------------------------------------------------------------------
   SUBROUTINE MOAB_SetGrid3D(ng, tile)
!-----------------------------------------------------------------------

      USE mod_param
      USE mod_scalars
      ! USE mod_coupling
      USE mod_grid
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer :: mLBi, mUBi, mLBj, mUBj

      integer :: roms3d_pid     ! MOAB ROMS App ID
      integer :: roms3d_compid  ! MOAB ROMS Component ID

      ! local data
      integer :: ierr, num_verts_in_cells, num_elements
      integer :: i, j, k, i1, offset
      real(r8), pointer :: moab_vert_coords(:)
      integer, pointer :: connectivity(:), connectivityg(:)
      character*100 tagname, lnum, filename, fileopts
      integer ndimcoord, ndimension, nmax_vpere, element_type, block_id
      integer tagtype, numco, tag_sto_len, ent_type, tagindex
      real(r8), pointer :: datasets(:, :, :)
      character(len=32) :: roms_appname    ! MOAB ROMS App name

      integer, pointer :: v3dlocal_ids(:, :, :)
      integer, pointer :: v3dglobal_ids(:, :, :)
      integer, pointer :: e3dglobal_ids(:, :, :)
!

         ! specify a unique external id;
         roms3d_compid = 625 + ng
         roms_appname = 'ROMS3D'//C_NULL_CHAR
#ifdef DISTRIBUTE
#ifdef MOAB_HAVE_MPI
         ierr = iMOAB_RegisterApplication(roms_appname, OCN_COMM_WORLD, &
                                          roms3d_compid, &
                                          roms3d_pid)
         CALL errorout(ierr, 'failed to register application')
#else
         ierr = iMOAB_RegisterApplication(roms_appname, &
                                          roms3d_compid, &
                                          roms3d_pid)
         CALL errorout(ierr, 'failed to register application')
#endif
#else
         ierr = iMOAB_RegisterApplication(roms_appname, &
                                          roms3d_compid, &
                                          roms3d_pid)
         CALL errorout(ierr, 'failed to register application')
#endif

      ndimension = 3 ! dimension of coordinates
      element_type = 9 ! MBHEX
      nmax_vpere = 8 ! MBHEX

      mLBi = griddata(ng)%LBi
      mUBi = griddata(ng)%UBi
      mLBj = griddata(ng)%LBj
      mUBj = griddata(ng)%UBj

      ! Eg.,
      ! LBi, UBi = [1, 5]
      ! LBj, UBj = [1, 5]
      ! num_verts_in_cells = 25
      num_verts_in_cells = (mUBi - mLBi + 1)*(mUBj - mLBj + 1)*(griddata(ng)%global_z)
      ndimcoord = 3*num_verts_in_cells
      allocate (moab_vert_coords(ndimcoord))
      ! allocate vertex ID arrays
      allocate (v3dlocal_ids(mLBi:mUBi, mLBj:mUBj, griddata(ng)%global_z))
      allocate (v3dglobal_ids(mLBi:mUBi, mLBj:mUBj, griddata(ng)%global_z))

      i1 = 1
      do k = 1, griddata(ng)%global_z
         do j = mLBj, mUBj
            do i = mLBi, mUBi
               v3dlocal_ids(i, j, k) = i1
               v3dglobal_ids(i, j, k) = i + (j - 1)*griddata(ng)%global_x &
                              + (k - 1)*griddata(ng)%global_x*griddata(ng)%global_y
#ifdef ANA_GRID
               moab_vert_coords(3*i1 - 2) = GRID(ng)%xp(i, j)/1E3
               moab_vert_coords(3*i1 - 1) = GRID(ng)%yp(i, j)/1E3
               moab_vert_coords(3*i1)     = GRID(ng)%z_r(i, j, k)
#else
               moab_vert_coords(3*i1 - 2) = GRID(ng)%lonp(i, j) !xp(i,j)
               moab_vert_coords(3*i1 - 1) = GRID(ng)%latp(i, j) !yp(i,j)
               moab_vert_coords(3*i1)     = GRID(ng)%z_r(i, j, griddata(ng)%global_z-k+1)
#endif
               i1 = i1 + 1
            end do
         end do
      end do

      print *, "Rank: ", mpi_rank, " Total vertices: ", num_verts_in_cells

      ! now let us create the vertices in MOAB
      ierr = iMOAB_CreateVertices(roms3d_pid, ndimcoord, ndimension, moab_vert_coords)
      call errorout(ierr, 'failed to create vertices')

      ! Eg.,
      ! UBi-LBi = 4
      ! UBj-LBj = 4
      ! num_elements = 16
      !
      !     5   10  15   20   25
      !       D   H    L    P
      !     4   9   14   19   24
      !       C   G    K    O
      !     3   8   13   18   23
      !       B   F    J    N
      !     2   7   12   17   22
      !       A   E    I    M
      !     1   6   11   16   21
      !
      ! eloffset = mpi_rank * num_elements
      num_elements = (mUBi - mLBi)*(mUBj - mLBj)*(griddata(ng)%global_z - 1)
      allocate (connectivity(num_elements*8))
      ! allocate (connectivityg(num_elements*8))
      ! allocate element ID arrays
      !allocate (e3dlocal_ids(mLBi:mUBi - 1, mLBj:mUBj - 1, griddata(ng)%global_z - 1))
      allocate (e3dglobal_ids(mLBi:mUBi - 1, mLBj:mUBj - 1, griddata(ng)%global_z - 1))
      i1 = 1
      offset = 1

      do k = 1, griddata(ng)%global_z - 1
         do j = mLBj, mUBj - 1
            do i = mLBi, mUBi - 1
               !elocal_ids(i, j) = i1
               e3dglobal_ids(i, j, k) = i + (j - 1)*griddata(ng)%global_x &
                                     + (k - 1)*griddata(ng)%global_x*griddata(ng)%global_y
               connectivity(offset + 0) = v3dlocal_ids(i, j, k)
               connectivity(offset + 1) = v3dlocal_ids(i + 1, j, k)
               connectivity(offset + 2) = v3dlocal_ids(i + 1, j + 1, k)
               connectivity(offset + 3) = v3dlocal_ids(i, j + 1, k)
               connectivity(offset + 4) = v3dlocal_ids(i, j, k + 1)
               connectivity(offset + 5) = v3dlocal_ids(i + 1, j, k + 1)
               connectivity(offset + 6) = v3dlocal_ids(i + 1, j + 1, k + 1)
               connectivity(offset + 7) = v3dlocal_ids(i, j + 1, k + 1)
               offset = offset + 8
               i1 = i1 + 1
            end do
         end do
      end do

      print *, "Rank: ", mpi_rank, " Total elements: ", num_elements

      block_id = 100*roms3d_pid + mpi_rank ! we should have only one block right now

      ! next let us create the QUAD elements in MOAB with the right connectivity
      ierr = iMOAB_CreateElements(roms3d_pid, num_elements, element_type, nmax_vpere, connectivity, block_id);
      call errorout(ierr, 'failed to create hexes')

      ! set the global id for vertices
      ! first, retrieve the tag
      tagname = 'GLOBAL_ID'//C_NULL_CHAR
      tagtype = 0  ! dense, integer
      numco = 1
      ierr = iMOAB_DefineTagStorage(roms3d_pid, tagname, tagtype, numco, tagindex)
      call errorout(ierr, 'failed to get global id tag')

      ! now set the values
      ent_type = 0 ! vertex type
      ierr = iMOAB_SetIntTagStorage(roms3d_pid, tagname, num_verts_in_cells, ent_type, &
                                    v3dglobal_ids)
      call errorout(ierr, 'failed to set global id tag for vertices')

      ! set global id tag for elements
      ent_type = 1 ! now set the global id tag on elements
      ierr = iMOAB_SetIntTagStorage(roms3d_pid, tagname, num_elements, ent_type, &
                                    e3dglobal_ids)
      call errorout(ierr, 'failed to set global id tag for hexes')

#ifdef MOAB_HAVE_MPI
      if (mpi_size > 1) then
         ! get next block
         ierr = iMOAB_ResolveSharedEntities(roms3d_pid, num_verts_in_cells, &
                                            v3dglobal_ids)
         call errorout(ierr, 'failed to resolve shared entities')

         !   ierr = iMOAB_MergeVertices( roms3d_pid )
         !   call errorout(ierr, 'failed to parallel merge mesh')

         !   ierr = iMOAB_DetermineGhostEntities ( roms3d_pid, 2, 2, 1 )
         !   call errorout(ierr, 'failed to create ghost entities')
      end if
#endif

      allocate (datasets(mLBi:mUBi, mLBj:mUBj, N(ng)))

      ! create a tag for bathymetry data and set in a tag
      tagname = 'Huon'//C_NULL_CHAR
      tagtype = 1  ! dense, double
      numco = 1
      ! first, create the tag
      ierr = iMOAB_DefineTagStorage(roms3d_pid, tagname, tagtype, numco, tagindex)
      call errorout(ierr, 'failed to get bathymetry tag')

      ! now set the values
      ent_type = 0 ! vertex type

      datasets = GRID(ng)%Huon(mLBi:mUBi, mLBj:mUBj, :)
      ierr = iMOAB_SetDoubleTagStorage(roms3d_pid, tagname, num_verts_in_cells, ent_type, &
                                       datasets)
      call errorout(ierr, 'failed to set bathymetry tag on vertices')

      ! create a tag for zlevels data and set in a tag
      tagname = 'z_r'//C_NULL_CHAR
      tagtype = 1  ! dense, double
      numco = 1
      ! first, create the tag
      ierr = iMOAB_DefineTagStorage(roms3d_pid, tagname, tagtype, numco, tagindex)
      call errorout(ierr, 'failed to get coriolis tag')

      ! ! now set the values
      datasets = GRID(ng)%z_r(mLBi:mUBi, mLBj:mUBj, :)
      ent_type = 0 ! vertex type
      ierr = iMOAB_SetDoubleTagStorage(roms3d_pid, tagname, num_verts_in_cells, ent_type, &
                                       datasets)
      call errorout(ierr, 'failed to set coriolis tag on vertices')

      deallocate (v3dglobal_ids)
      deallocate (v3dlocal_ids)
      deallocate (e3dglobal_ids)
      deallocate (datasets)
      deallocate (moab_vert_coords)
      deallocate (connectivity)

      ierr = iMOAB_UpdateMeshInfo(roms3d_pid)
      call errorout(ierr, 'failed to update mesh info')

      filename = 'roms_3d_mesh.h5m'//C_NULL_CHAR
      fileopts = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
      ierr = iMOAB_WriteMesh(roms3d_pid, filename, fileopts)
      call errorout(ierr, 'failed to write mesh to disk')

      ! destroy the application instance
      ierr = iMOAB_DeregisterApplication(roms3d_pid)
      CALL errorout(ierr, 'failed to de-register application')

      RETURN
   END SUBROUTINE MOAB_SetGrid3D
#endif

END MODULE moab_coupler_mod
!
#endif ! HAVE_MOAB
