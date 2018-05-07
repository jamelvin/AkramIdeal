!!****if* source/Simulation/SimulationMain/Sedov/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!! 
!! SYNOPSIS
!!
!!  call Simulation_initBlock(integer(IN) :: blockId)
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Sedov spherical
!!  explosion problem.
!!
!!  References:  Sedov, L. I., 1959, Similarity and Dimensional Methods
!!                 in Mechanics (New York:  Academic)
!!
!!               Landau, L. D., & Lifshitz, E. M., 1987, Fluid Mechanics,
!!                 2d ed. (Oxford:  Pergamon)
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!  
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_tempAmbient    Initial ambient temperature
!!  sim_exl            Electrode lower x-coordinates
!!  sim_exu            Electrode upper x-coordinates
!!  sim_exl            Electrode lower y-coordinates
!!  sim_exu            Electrode upper y-coordinates
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!***

!!REORDER(4): solnData


subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY: sim_xMax, sim_xMin, sim_yMax, sim_yMin, sim_zMax, sim_zMin, &
     &  sim_exl, sim_exu, sim_eyl, sim_eyu, &
     &  sim_pxl, sim_pxu, sim_pyl, sim_pyu, &
     &  sim_tInitial, sim_gamma, sim_pAmbient, sim_rhoAmbient, sim_tempAmbient, &
     &  sim_smallT, sim_smallP, &
     &  sim_vx, sim_vy, sim_vz, &
     &  R_sp

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
    Grid_getCellCoords, Grid_putPointData
  
  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer,intent(IN) ::  blockId
  
  integer  ::  i, j, k
  integer  ::  ii, jj, kk
  real     ::  xx, dxx, yy, dyy, zz, dzz
  real     ::  vx, vy, vz, p, rho, e, ek, eint, ti

  real,allocatable,dimension(:) :: xCoord, yCoord, zCoord
  integer,dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX, sizeY, sizeZ
  integer,dimension(MDIM) :: axis
  real, dimension(:,:,:,:),pointer :: solnData

  logical :: gcell = .true.
 

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  ! get the coordinate information for the current block from the database

  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoord(sizeX)); xCoord = 0.0
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoord(sizeY)); yCoord = 0.0
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoord(sizeZ)); zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER, gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER, gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)

  call Grid_getBlkPtr(blockId,solnData)


  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     ! Find a real difference between z's if problem is >= 3D
     if (NDIM > 2) then
        if (k .eq. 1) then
           dzz = zCoord(2) - zCoord(1) 
        else
           dzz = zCoord(k) - zCoord(k-1) 
        endif
     ! Otherwise this problem is <= 2D, so dzz is meaningless
     else
       dzz = 0.0
     endif
     zz = zCoord(k)

     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        ! Find a real difference between y's if problem is >= 2D
        if (NDIM > 1) then
           if (j .eq. 1) then
              dyy = yCoord(2) - yCoord(1) 
           else
              dyy = yCoord(j) - yCoord(j-1) 
           endif
        ! Otherwise this problem is <= 1D, so dyy is meaningless
        else
          dyy = 0.0
        endif
        yy = yCoord(j)
        
        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
           xx = xCoord(i)
           if (i .eq. 1) then
              dxx = xCoord(2) - xCoord(1) 
           else
              dxx = xCoord(i) - xCoord(i-1) 
           endif

           rho = sim_rhoAmbient
           vx  = sim_vx
           vy  = sim_vy
           vz  = sim_vz
           ek  = 0.5*(vx*vx + vy*vy + vz*vz)
           ti  = sim_tempAmbient  ! ambient temperature

           ! implement the initial temperature in the plasma channel
           if (xx >= sim_pxl .and. xx <= 2.0*sim_pxu .and. yy >= sim_pyl .and. yy <= sim_pyu) then
             call  expFunc (xx, sim_pxl, 18500.0, 2.0*sim_pxu, 300.0, ti)
             !ti = 10000.0 ! set is to 10000 K like in Ekici paper
           endif          
           !write (*,*) ti
 
           p   = sim_rhoAmbient*ti*R_sp !PRAO: P = rho*T*R 
           e   = p/(sim_gamma-1.)
           eint= e/rho
           e   = e/rho + ek
           e   = max (e, sim_smallP)


           ! implement the electrode geometry here
           ! PRAO: is there a need for this?
           if ((xx > sim_exl .and. xx < sim_exu) .and. (yy < sim_eyl .or. yy > sim_eyu) ) then
             ti = sim_tempAmbient ! temp inside electrode is 300 K
           endif

           
           axis(IAXIS)=i
           axis(JAXIS)=j
           axis(KAXIS)=k

           solnData(DENS_VAR,i,j,k)=rho
           solnData(PRES_VAR,i,j,k)=p
           solnData(ENER_VAR,i,j,k)=e
#ifdef EINT_VAR
           solnData(EINT_VAR,i,j,k)=eint
#endif
           solnData(GAME_VAR,i,j,k)=sim_gamma
           solnData(GAMC_VAR,i,j,k)=sim_gamma
           solnData(VELX_VAR,i,j,k)=vx
           solnData(VELY_VAR,i,j,k)=vy
           solnData(VELZ_VAR,i,j,k)=vz
           solnData(TEMP_VAR,i,j,k)=ti
#ifdef BDRY_VAR
           solnData(BDRY_VAR,i,j,k)=    -1.0
#endif

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, p)
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, e)    
#ifdef EINT_VAR
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, axis, eint)
#endif
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, vx)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vy)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, vz)
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, ti)
#ifdef BDRY_VAR
           call Grid_putPointData(blockId, CENTER, BDRY_VAR, EXTERIOR, axis, -1.0)
#endif
        enddo
     enddo
  enddo

  call Grid_releaseBlkPtr(blockID, solnData)
  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

  return
end subroutine Simulation_initBlock


!!!!!!for interpolating to an exponential function!!!!!!

subroutine expFunc (x, x1, y1, x2, y2, y)

 implicit none
 real, intent(IN) :: x  ! interpolate at x
 real, intent(OUT) :: y  ! interpolated value at x 
 real, intent(IN) :: x1, y1, x2, y2 !two given points
 real :: a, b, test, test2 !coefficients of exp function

 test = log(y1)
 test2 = log(y2)
 b = (log(y1) - log(y2)) / (x1-x2)  
 a = y1/(exp(b*x1))

 y = a*exp(b*x)

 return

end subroutine expFunc
