!
! Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
!
! *********************************************************************
!
! This is the mock solver program.
!
! It requires the following command line arguments to get
! System Coupling host (string), port (integer),
! and name (string) given to this participant:
! --schost <host>
! --scport <port>
! --scname <name>
!
! This program is going to run as part of the overall coupled
! analysis. System Coupling is going to start this program
! automatically, so it will pass in the appropriate arguments
! as needed.
!
      program mockSolverFortran
!
      implicit none
      include 'syscSystemCouplingF.fi'
!
      integer :: argCount
      character(len=256) :: arg
      integer :: i
!
      character(len=256) :: scHost = ""
      character(len=256) :: scPortString
      integer :: scPort = 0
      character(len=256) :: scName = ""
      character(len=256) :: buildInfo = ""
!
      type(SyscErrorF) :: ret
!
! *********************************************************************
!
! Mesh representation:
!   10 nodes
!   4 quad elements
!   Rectangle, aligned with XY-plane,
!   stretching from the origin to point (1, 0.1, 0)
!
!  (0,0.1,0)                    (1,0.1,0)
!           6----7----8----9----10
!           |    |    |    |    |
!           |    |    |    |    |
!           1----2----3----4----5
!  (0,0,0)                       (1,0,0)
!
! *********************************************************************
!
!     Array of node coordinates in split array format.
      real(kind=8), target :: nodeCoordsX(10) =
     &                        (/0.0D0, 0.25D0, 0.5D0, 0.75D0, 1.0D0,
     &                          0.0D0, 0.25D0, 0.5D0, 0.75D0, 1.0D0 /)
!
      real(kind=8), target :: nodeCoordsY(10) =
     &                        (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,
     &                          0.1D0, 0.1D0, 0.1D0, 0.1D0, 0.1D0 /)
!
      real(kind=8), target :: nodeCoordsZ(10) = 0.0D0
!
!     Array of element node counts.
      integer(kind=4), target :: elemNodeCounts(4) = 4
!
!     Array of node ids for each element.
      integer(kind=4), target :: elemNodeIds(16) = (/ 1, 2, 7, 6,
     &                                                2, 3, 8, 7,
     &                                                3, 4, 9, 8,
     &                                                4, 5, 10, 9 /)
!
!     Array of temperature values (defined at nodes).
      real(kind=8), target :: temperature(10) = 300.0D0
!
!     Array of heat rate values (defined at elements)
!     values are in Watts - total of 20 [W] generated.
      real(kind=8), target :: heatRate(4) = 5.0D0
!
! *********************************************************************
!
!     Parse and print input arguments (host, port, participant name).
      argCount = command_argument_count()
      do i = 1, argCount
        call get_command_argument(i, arg)
        arg = trim(arg)
        if (arg .eq. '--schost' .and. i + 1 .le. argCount) then
          call get_command_argument(i + 1, scHost)
          scHost = trim(scHost)
        else if (arg .eq. '--scport' .and. i + 1 .le. argCount) then
          call get_command_argument(i + 1, scPortString)
          read (scPortString, *) scPort
        else if (arg .eq. '--scname' .and. i + 1 .le. argCount) then
          call get_command_argument(i + 1, scName)
          scName = trim(scName)
        endif
      enddo
!
      print *, 'Host: ', trim(scHost)
      print *, 'Port: ', scPort
      print *, 'Name: ', trim(scName)
      buildInfo = 'Mock Solver, language: FORTRAN'
!
!     Connection to System Coupling.
      ret = syscConnectF(scHost, scPort, scName, buildInfo)
      call processRetCode(ret)
!
!     Heavyweight data access and initialization.
!
!     Register access to mesh and variable data.
!     NOTE: The participant solver retains ownership
!     of mesh and variable field data arrays and is thus
!     responsible for allocating and deallocating memory
!     for those arrays. System Coupling only gets
!     read and/or write access to those arrays.
!
!     Register surface mesh access by providing pointer
!     to a function that provides mesh information for a given region.
      ret = syscRegisterSurfMeshAccessF(getSurfaceMesh)
      call processRetCode(ret)
!
!     Register variable data access by providing pointer
!     to a function that provides access to variable data for
!     a region.
!
!     Input scalar variable access.
      ret = syscRegisterInputScalarDataAccessF(getInputScalarData)
      call processRetCode(ret)
!
!     Output scalar variable access.
      ret = syscRegisterOutputScalarDataAccessF(getOutputScalarData)
      call processRetCode(ret)
!
!     Coupled analysis initialization.
      ret = syscInitializeAnalysisF()
      call processRetCode(ret)
!
!     Coupled analysis loop.
      do while (syscDoIterationF())
!
!       Inputs update.
!       All input variables will be up-to-date after this call.
        ret = syscUpdateInputsF()
        call processRetCode(ret)
!
!       *************************************************************
!       Solver iterations can be peformed here.
!       Input variable values can be consumed.
!       Output variable values must be generated.
!       ...
!       *************************************************************
!
!       Outputs update.
!       Provide convergence state to System Coupling,
!       and notify that outputs can be consumed.
        ret = syscUpdateOutputsF(SyscConverged)
        call processRetCode(ret)
      end do
!
!     Coupled analysis shutdown.
      ret = syscDisconnectF()
      call processRetCode(ret)
!
      contains
!
!     Process error code.
      subroutine processRetCode(ret)
        type(SyscErrorF) :: ret
        if (ret%retcode .NE. SyscStatusOk) then
          print *, "ERROR: ", trim(ret%message)
          stop 1
        endif
      end subroutine processRetCode
!
!       Function that returns surface mesh, given a region.
        function getSurfaceMesh(regName) result(surfaceMesh)
          character(len=SyscStrLen), intent(in) :: regName
          type(SyscSurfaceMeshF) :: surfaceMesh
          type(SyscOutputVectorDataF) :: nodeCoordsData
          type(SyscOutputIntegerDataF) :: elemNodeCountsData
          type(SyscOutputIntegerDataF) :: elemNodeIdsData
          type(SyscNodeDataF) :: nodes
          type(SyscFaceDataF) :: faces
          character(len=SyscStrLen) :: message
          surfaceMesh = syscGetSurfaceMeshF()
          if (regName .eq. "FSI") then
!           node coords
            nodeCoordsData = syscGetOutputVectorDataSplitF(
     &        nodeCoordsX, nodeCoordsY, nodeCoordsZ, 
     &        size(nodeCoordsX, kind=8))
!           elem node counts
            elemNodeCountsData = syscGetOutputIntegerDataF(
     &        elemNodeCounts, size(elemNodeCounts, kind=8))
!           elem node ids
            elemNodeIdsData = syscGetOutputIntegerDataF(
     &        elemNodeIds, size(elemNodeIds, kind=8))
            nodes = syscGetNodeDataF(nodeCoordsData)
            faces = syscGetFaceDataF(
     &        syscGetElementNodeCountDataF(elemNodeCountsData),
     &        syscGetElementNodeConnectivityDataF(elemNodeIdsData))
            surfaceMesh = syscGetSurfaceMeshF(nodes, faces)
          else
            message = "getSurfaceMesh error: unknown region"
            message = trim(message) // " " // trim(regName)
            call syscFatalErrorF(message)
          endif
        end function getSurfaceMesh
!
!       Function that returns scalar input variable values,
!       given a region and a variable.
        function getInputScalarData(regName,
     &                              varName) result(ret)
          character(len=SyscStrLen), intent(in) :: regName
          character(len=SyscStrLen), intent(in) :: varName
          character(len=SyscStrLen) :: message
          type(SyscInputScalarDataF) :: ret
          ret = syscGetInputScalarDataF()
          if (regName .eq. "FSI" .and. varName .eq. "Temperature") then
            ret = syscGetInputScalarDataF(
     &        temperature, size(temperature, kind=8))
          else
            message = "getInputScalarData error: unknown variable"
            message = trim(message) // " " // trim(varName)
            message = trim(message) // " on region "
            message = trim(message) // " " // trim(regName)
            call syscFatalErrorF(message)
          endif
        end function getInputScalarData
!
!       Function that returns scalar output variable values,
!       given a region and a variable.
        function getOutputScalarData(regName, varName) result(ret)
          character(len=SyscStrLen), intent(in) :: regName
          character(len=SyscStrLen), intent(in) :: varName
          character(len=SyscStrLen) :: message
          type(SyscOutputScalarDataF) :: ret
          if (regName .eq. "FSI" .and. varName .eq. "Heat Rate") then
            ret = syscGetOutputScalarDataF(
     &        heatRate, size(heatRate, kind=8))
          else
            message = "getOutputScalarData error: unknown variable"
            message = trim(message) // " " // trim(varName)
            message = trim(message) // " on region "
            message = trim(message) // " " // trim(regName)
            call syscFatalErrorF(message)
          endif
        end function getOutputScalarData
!
      end program mockSolverFortran
!