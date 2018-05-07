!!****if* source/Simulation/SimulationMain/Sedov/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for Sedov Spherical Explosion 
!!  problem.
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_tempAmbient    Initial ambient temperature
!!  sim_minRhoInit     Density floor for initial condition
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!  sim_forceCenterDerefine  Try to force low refinement level around explosion center?
!!  sim_centerRefineLevel    Desired refinement level at center (if "forcing")
!!  sim_derefineRadius       Radius of center region to force derefinement
!!  sim_profFileName   Name of file from which to read a 1D Sedov solution for the
!!                     initial condition
!!
!!***

subroutine Simulation_init()

  use Simulation_data 
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface,           ONLY : Logfile_stamp
  use ut_generalInterface,         ONLY : ut_getFreeFileUnit
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none

#include "Flash.h"
#include "constants.h"

  logical :: threadBlockListBuild, threadWithinBlockBuild  
  real    :: vctr
  real    :: R_u, eos_singleSpeciesA 

  call Driver_getMype(MESH_COMM,   sim_meshMe)
  call Driver_getMype(GLOBAL_COMM, sim_globalMe)
  call Driver_getNumProcs(GLOBAL_COMM, sim_globalNumProcs)

  call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
  call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
  call RuntimeParameters_get('sim_tempAmbient', sim_tempAmbient)
  call RuntimeParameters_get('sim_nsubzones',sim_nSubZones)
  call RuntimeParameters_get('sim_exl',sim_exl)
  call RuntimeParameters_get('sim_exu',sim_exu)
  call RuntimeParameters_get('sim_eyl',sim_eyl)
  call RuntimeParameters_get('sim_eyu',sim_eyu)
  call RuntimeParameters_get('sim_pxl',sim_pxl)
  call RuntimeParameters_get('sim_pxu',sim_pxu)
  call RuntimeParameters_get('sim_pyl',sim_pyl)
  call RuntimeParameters_get('sim_pyu',sim_pyu)
  call RuntimeParameters_get('sim_vx',sim_vx)
  call RuntimeParameters_get('sim_vy',sim_vy)
  call RuntimeParameters_get('sim_vz',sim_vz)  
  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('eos_singleSpeciesA', eos_singleSpeciesA)
  call PhysicalConstants_get("ideal gas constant", R_u)
  call RuntimeParameters_get('smallx', sim_smallX)
  call RuntimeParameters_get('smlrho', sim_smallRho)
  call RuntimeParameters_get('sim_minRhoInit', sim_minRhoInit)
  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('smallT', sim_smallT)
  call RuntimeParameters_get('smallu', sim_smallu)
  call RuntimeParameters_get('xmin',sim_xMin)
  call RuntimeParameters_get('ymin',sim_yMin)
  call RuntimeParameters_get('zmin',sim_zMin)
  call RuntimeParameters_get('xmax',sim_xMax)
  call RuntimeParameters_get('ymax',sim_yMax)
  call RuntimeParameters_get('zmax',sim_zMax)
  call RuntimeParameters_get('tinitial',sim_tInitial)
  call RuntimeParameters_get('sim_forceCenterDerefine',sim_forceCenterDerefine)
  call RuntimeParameters_get('sim_centerRefineLevel',sim_centerRefineLevel)
  call RuntimeParameters_get('sim_derefineRadius',sim_derefineRadius)
  call RuntimeParameters_get('sim_profFileName',sim_profFileName)
  call RuntimeParameters_get('sim_bcSetBdryVar',sim_bcSetBdryVar)
  call RuntimeParameters_get('sim_earliestLSTime',sim_earliestLSTime)
  call RuntimeParameters_get('sim_latestLSTime',sim_latestLSTime)
  call RuntimeParameters_get('sim_smallestNormRadius',sim_smallestNormRadius)
  call RuntimeParameters_get('sim_largestNormRadius',sim_largestNormRadius)
  call RuntimeParameters_get('sim_oneLevelIntegralsOnly',sim_oneLevelIntegralsOnly)
  call RuntimeParameters_get('sim_integralsLevel',       sim_integralsLevel)

  sim_useProfileFromFile = .FALSE.
  if (sim_nProfile > 1) then
     sim_useProfileFromFile = (trim(sim_profFileName) .NE. "/dev/null")
  end if

  ! calculate the specific gas constant
  R_sp = R_u/eos_singleSpeciesA


  sim_inSubZones = 1./real(sim_nSubZones)
  sim_inSubzm1   = 1./real(max(sim_nSubZones-1,1))
  sim_inszd      = sim_inSubZones**NDIM
  
  !
  !  Calculate the initial volume and interior pressure.
  !
  !if (NDIM .eq. 1) then
     !vctr = 2. * sim_rInit
  !elseif (NDIM .eq. 2) then
     !vctr = sim_pi * sim_rInit**2
  !else
     !vctr = 4./3.*sim_pi*sim_rInit**3
  !endif
  
  !sim_pExp    = (sim_gamma-1.) * sim_expEnergy / vctr   ! only used if tinitial .LE. 0.

  !sim_pExp = sim_rhoAmbient*15000.0*288.0946 !PRAO


  call RuntimeParameters_get("threadBlockListBuild", threadBlockListBuild)
  call RuntimeParameters_get("threadHydroBlockList", sim_threadBlockList)

  call RuntimeParameters_get("threadWithinBlockBuild", threadWithinBlockBuild)
  call RuntimeParameters_get("threadHydroWithinBlock", sim_threadWithinBlock)

  if (sim_threadBlockList .and. .not. threadBlockListBuild) then
     call Logfile_stamp('WARNING! Turning off block list threading '//&
          'because FLASH is not built appropriately','[Simulation_init]')
     sim_threadBlockList = .false.
  end if
  if (sim_threadWithinBlock .and. .not. threadWithinBlockBuild) then
     call Logfile_stamp('WARNING! Turning off within block threading '//&
          'because FLASH is not built appropriately','[Simulation_init]')
     sim_threadWithinBlock = .false.
  end if

  if (sim_useProfileFromFile) then
     !call sim_readProfile()
     !if (sim_tinitial > 0.0) call sim_scaleProfile(sim_tinitial)
     write (*,*) "GetOut"
  end if
  sim_analyticTime = -HUGE(1.0)
  sim_analyticGen  = -1
  if (sim_meshMe == MASTER_PE) then
     print*,'sim_rhoAmbient is',sim_rhoAmbient
  end if
  
  ! Open file for writing the numerical solution discretized onto the FLASH grid:
  sim_fileUnitOutNum = ut_getFreeFileUnit()
  open(unit=sim_fileUnitOutNum, status='SCRATCH')
!!$
!!$  ! Write the file header:
!!$  if(sim_globalME == MASTER_PE) then
!!$     write(sim_fileUnitOutNum,'(a10,7a15)') &
!!$          '#   cellno', 'x       ', 'y     ', 'dens    ', 'pres    ', 'velx    ', 'eint (spec.)'
!!$  end if
  
  ! Open file for writing the analytical solution discretized onto the FLASH grid:
  sim_fileUnitOutAna = ut_getFreeFileUnit()
!!$  open(unit=sim_fileUnitOutAna, file="sedSol-ana.out", form="formatted", position='append')
!!$
!!$  ! Write the file header:
!!$  if(sim_globalME == MASTER_PE) then
!!$     write(sim_fileUnitOutAna,'(a10,7a15)') &
!!$          '#   cellno', 'x       ', 'y     ', 'dens    ', 'pres    ', 'velx    ', 'eint density'
!!$  end if
  close(sim_fileUnitOutNum)

end subroutine Simulation_init
