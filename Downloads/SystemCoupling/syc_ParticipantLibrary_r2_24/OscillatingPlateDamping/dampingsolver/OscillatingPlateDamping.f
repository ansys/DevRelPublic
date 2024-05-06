!
! Copyright ANSYS, Inc. Unauthorized use, distribution, 
! or duplication is prohibited.
!
! *********************************************************************
!
! See Participant Library's Oscillating Plate Damping tutorial
! for a detailed description of the problem.
!
! This is the oscillating plate damping solver program.
!
! It requires the following command line arguments to get
! System Coupling host (string), port (integer),
! name (string), and damping coefficient (double) given to this 
! participant:
! --schost <host>
! --scport <port>
! --scname <name>
! --dampcoeff <damping coefficient>
!
! This program is going to run as part of the overall coupled
! analysis. System Coupling is going to start this program
! automatically, so it will pass in the appropriate arguments as 
! needed.
!
      program oscillatingPlateDampingFortran
!
      implicit none
      include 'syscSystemCouplingF.fi'
!
      integer :: argCount
      character(len=256) :: arg
      integer :: i
!
      character(len=256) :: scHost = ""
      integer :: scPort = 0
      character(len=256) :: scName = ""
      logical :: isSetupMode = .FALSE.
      character(len=SyscFilePathLen) :: restartPoint = ""
      real(kind=8) :: dampCoeff = 0.0
      character(len=256) :: buildInfo
!
      type(SyscErrorF) :: ret
!
      type(SyscRegionF) :: region
      type(SyscVariableF) :: forceVar
      type(SyscVariableF) :: displacementVar
      type(SyscSetupInfoF) :: setupInfo
!
!     Define surface mesh.
!
!     Define node coordinates in a single array (compact form).
!     There are 24 nodes.
      real(kind=8), target :: nodeCoords(72) =
     &          (/1.00600004D1, 1.00000000D0,  4.00000006D-1,
     &            1.00600004D1, 1.00000000D0,  0.00000000D0,
     &            1.00000000D1, 1.00000000D0,  4.00000006D-1,
     &            1.00000000D1, 1.00000000D0,  0.00000000D0,
     &            1.00600004D1, 0.00000000D0,  4.00000006D-1,
     &            1.00600004D1, 2.00000003D-1, 4.00000006D-1,
     &            1.00600004D1, 4.00000006D-1, 4.00000006D-1,
     &            1.00600004D1, 6.00000024D-1, 4.00000006D-1,
     &            1.00600004D1, 8.00000012D-1, 4.00000006D-1,
     &            1.00600004D1, 8.00000012D-1, 0.00000000D0,
     &            1.00600004D1, 6.00000024D-1, 0.00000000D0,
     &            1.00600004D1, 4.00000006D-1, 0.00000000D0,
     &            1.00600004D1, 2.00000003D-1, 0.00000000D0,
     &            1.00600004D1, 0.00000000D0, 0.00000000D0,
     &            1.00000000D1, 8.00000012D-1, 4.00000006D-1,
     &            1.00000000D1, 6.00000024D-1, 4.00000006D-1,
     &            1.00000000D1, 4.00000006D-1, 4.00000006D-1,
     &            1.00000000D1, 2.00000003D-1, 4.00000006D-1,
     &            1.00000000D1, 0.00000000D0, 4.00000006D-1,
     &            1.00000000D1, 0.00000000D0, 0.00000000D0,
     &            1.00000000D1, 8.00000012D-1, 0.00000000D0,
     &            1.00000000D1, 6.00000024D-1, 0.00000000D0,
     &            1.00000000D1, 4.00000006D-1, 0.00000000D0,
     &            1.00000000D1, 2.00000003D-1, 0.00000000D0/)
!
!     Define element node counts (11 quad elements).
      integer(kind=4), target :: elementNodeCounts(11) = 4
!
!     Define element node ids for each element (11 quad elements).
      integer(kind=4), target :: elementNodeIds(44) = 
     &          (/ 13, 4, 5, 12,
     &             12, 5, 6, 11,
     &             11, 6, 7, 10,
     &             10, 7, 8, 9,
     &             9, 8, 0, 1,
     &             23, 17, 18, 19,
     &             22, 16, 17, 23,
     &             21, 15, 16, 22,
     &             20, 14, 15, 21,
     &             3, 2, 14, 20,
     &             2, 3, 1, 0 /)
!
!     Define array of nodal displacements
!     (24 nodes, 3 components per node).
      real(kind=8), target :: displacement(72) = 0.0D0
!
!     Define array of nodal forces
!     (24 nodes, 3 components per node).
      real(kind=8), target :: force(72) = 0.0D0
!
      integer :: currStep = 0
      type(SyscTimeStepF) :: timeStep
      real(kind=8) :: dt
!
! *********************************************************************
!
!     Parse and print input arguments (host, port, participant name).
      argCount = command_argument_count()
      do i = 1, argCount
        call get_command_argument(i, arg)
        arg = trim(arg)
        if (arg .eq. '--schost' .and. i + 1 .le. argCount) then
          call get_command_argument(i + 1, arg)
          scHost = trim(arg)
        else if (arg .eq. '--scport' .and. i + 1 .le. argCount) then
          call get_command_argument(i + 1, arg)
          read (arg, *) scPort
        else if (arg .eq. '--scname' .and. i + 1 .le. argCount) then
          call get_command_argument(i + 1, arg)
          scName = trim(arg)
        else if (arg .eq. '--dampcoeff' .and. i + 1 .le. argCount) then
          call get_command_argument(i + 1, arg)
          read (arg, *) dampCoeff
        else if (arg .eq. '--screstart' .and. i + 1 .le. argCount) then
          call get_command_argument(i + 1, arg)
          restartPoint = trim(arg)
        else if (arg .eq. '--scsetup') then
            isSetupMode = .TRUE.
        endif
      enddo
!
      print *, 'Host: ', trim(scHost)
      print *, 'Port: ', scPort
      print *, 'Name: ', trim(scName)
      if (trim(restartPoint) .ne. '') then
        print *, 'Restart: ', trim(restartPoint)
      endif
      print *, 'Damping coefficient: ', dampCoeff
!
      buildInfo = 'Oscillating Plate Damping Solver, language: Fortran'
!
!     Connection to System Coupling.
      ret = syscConnectF(scHost, scPort, scName, buildInfo)
      call processRetCode(ret)
!
      if (isSetupMode) then

!       Setup mode.
!
!       Create surface region "plate".
        region = syscGetRegionF("plate", SyscSurface)
!
!       Create the force variable.
        forceVar = syscGetVariableF(
     &    "Damping Force", SyscVector, .TRUE., SyscNode)
!
!       Create the displacement variable.
        displacementVar = syscGetVariableF(
     &    "Mesh Displacement", SyscVector, .FALSE., SyscNode)
!
!       Add force as output and displacement as input to plate.
        ret = syscAddOutputVariableF(region, forceVar)
        call processRetCode(ret)
        ret = syscAddInputVariableF(region, displacementVar)
        call processRetCode(ret)
!
!       Add region to the coupling setup.
        ret = syscAddRegionF(region)
        call processRetCode(ret)
!
!       Create setup info structure. Set analysis type to transient.
!       Set restarts supported to true.
        setupInfo = syscGetSetupInfoF(SyscTransient, .TRUE.)
!
!       Complete the setup.
        ret = syscCompleteSetupF(setupInfo)
        call processRetCode(ret)

      else
!       Heavyweight data access and initialization.
!
!       Register access to mesh and variable data.
!       NOTE: The participant solver retains ownership
!       of mesh and variable field data arrays and is thus
!       responsible for allocating and deallocating memory
!       for those arrays. System Coupling only gets
!       read and/or write access to those arrays.
!
!       Register surface mesh access by providing pointer
!       to a function that provides mesh information for 
!       a given region.
        ret = syscRegisterSurfMeshAccessF(getSurfaceMesh)
        call processRetCode(ret)
!
!       Register variable data access by providing pointer
!       to a function that provides access to variable data for
!       a region.
!
!       Input vector variable access.
        ret = syscRegisterInputVectorDataAccessF(getInputVectorData)
        call processRetCode(ret)
!
!       Output vector variable access.
        ret = syscRegisterOutputVectorDataAccessF(getOutputVectorData)
        call processRetCode(ret)
!
!       Restart point creation callback registration.
        ret = syscRegisterRestartPointCreationF(createRestartPoint)
        call processRetCode(ret)
!
        if (trim(restartPoint) .ne. '') then
          call resumeState()
        endif
!
!       Coupled analysis initialization.
        ret = syscInitializeAnalysisF()
        call processRetCode(ret)
!
!       Coupling analysis loop.
!
!       Coupling time step loop.
        do while (syscDoTimeStepF())
!
!         Increment time step.
          currStep = currStep + 1
!
!         Update nodal coordinates.
          nodeCoords = nodeCoords + displacement
!
!         Coupling iteration loop.
          do while (syscDoIterationF())
!
!           Inputs update.
!           All input variables will be up-to-date after this call.
            ret = syscUpdateInputsF()
            call processRetCode(ret)
!
!           Get time step size, which is used to calculate
!           nodal velocities.
            timeStep = syscGetCurrentTimeStepF()
            dt = timeStep%timeStepSize
!
!           Calculate the damping force.
            force = - dampCoeff * displacement / dt
!
!           Outputs update.
!           Note: convergence status is Complete, since this is
!           not an iterative solver.
            ret = syscUpdateOutputsF(SyscComplete)
            call processRetCode(ret)
!
          end do
!
        end do
!
      endif
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
! Define a function that creates a restart point and returns
! a string that identifies that restart point.
! In this implementation, the file with up-to-date
! nodal coordinates and force values
! will be written and the name of that file is returned.
! Only those values are required to perform the restart
! for this problem.
        function createRestartPoint() result(ret)
          character(len=SyscFilePathLen) :: ret
          integer :: n
          write(ret, "(A8,I0,A4)") "restart.",currStep,".txt"
          open(1, file = trim(ret), status = 'replace')
          do n = 0, 23
            write(1, "(ES15.8e2,A1)", advance = 'no')
     &        nodeCoords(n * 3 + 1), ","
            write(1, "(ES15.8e2,A1)", advance = 'no')
     &        nodeCoords(n * 3 + 2), ","
            write(1, "(ES15.8e2,A1)", advance = 'no')
     &        nodeCoords(n * 3 + 3), ","
            write(1, "(ES15.8e2,A1)", advance = 'no')
     &        force(n * 3 + 1), ","
            write(1, "(ES15.8e2,A1)", advance = 'no')
     &        force(n * 3 + 2), ","
            write(1, "(ES15.8e2)", advance = 'yes')
     &        force(n * 3 + 3)
          end do
          close(1)
        end function createRestartPoint
!
! Error checking is skipped here for simplicity.
! Assume the file exists and is in the format
! consistent with how it gets written
! in the createRestartPoint() function.
        subroutine resumeState()
          integer :: split1
          integer :: split2
          integer :: n
!
!         Set current time step.
          split1 = scan(restartPoint, ".")
          split2 = scan(restartPoint(split1+1:), ".")
          read (restartPoint(split1+1:split1+split2-1), *) currStep
!
!         Set node coordinates and forces from the restart file.
          open(1, file = trim(restartPoint), status = 'old')
          do n = 0, 23
            read(1, *) nodeCoords(n * 3 + 1), nodeCoords(n * 3 + 2),
     &                 nodeCoords(n * 3 + 3), force(n * 3 + 1),
     &                 force(n * 3 + 2), force(n * 3 + 3)
          end do
          close(1)
!
        end subroutine resumeState
!
!       Define function that returns surface mesh, given a region.
        function getSurfaceMesh(regName) result(surfaceMesh)
          character(len=256), intent(in) :: regName
          type(SyscSurfaceMeshF) :: surfaceMesh
          type(SyscOutputVectorDataF) :: nodeCoordsData
          type(SyscOutputIntegerDataF) :: elemNodeCountsData
          type(SyscOutputIntegerDataF) :: elemNodeIdsData
          type(SyscNodeDataF) :: nodes
          type(SyscFaceDataF) :: faces
          character(len=SyscStrLen) :: message
          surfaceMesh = syscGetSurfaceMeshF()
          if (regName .eq. "plate") then
!           Node coordinates.
            nodeCoordsData = syscGetOutputVectorDataCompactF(
     &        nodeCoords, size(nodeCoords, kind=8) / 3)
!           Element node counts.
            elemNodeCountsData = syscGetOutputIntegerDataF(
     &        elementNodeCounts, size(elementNodeCounts, kind=8))
!           Element node ids.
            elemNodeIdsData = syscGetOutputIntegerDataF(
     &        elementNodeIds, size(elementNodeIds, kind=8))
!           Nodes.
            nodes = syscGetNodeDataF(nodeCoordsData)
!           Faces.
            faces = syscGetFaceDataF(
     &        syscGetElementNodeCountDataF(elemNodeCountsData),
     &        syscGetElementNodeConnectivityDataF(elemNodeIdsData))
!           Surface mesh.
            surfaceMesh = syscGetSurfaceMeshF(nodes, faces)
          else
            message = "getSurfaceMesh error: unknown region"
            message = trim(message) // " " // trim(regName)
            call syscFatalErrorF(message)
          endif
        end function getSurfaceMesh
!
!       Define a function that returns vector input variable values,
!       given a region and a variable.
        function getInputVectorData(regName,
     &    varName) result(ret)
          character(len=SyscStrLen), intent(in) :: regName
          character(len=SyscStrLen), intent(in) :: varName
          character(len=SyscStrLen) :: message
          type(SyscInputVectorDataF) :: ret
          ret = syscGetInputVectorDataF()
          if (regName .eq. "plate" .and. 
     &      varName .eq. "Mesh Displacement") then
            ret = syscGetInputVectorDataCompactF(
     &        displacement, size(displacement, kind=8) / 3)
          else
            message = "getInputVectorData: incorrect region or "
            message = trim(message) // "variable name"
            call syscFatalErrorF(message)
          endif
        end function getInputVectorData
!
!       Define a function that returns vector output variable values,
!       given a region and a variable.
        function getOutputVectorData(
     &    regName,
     &    varName) result(ret)
          character(len=SyscStrLen), intent(in) :: regName
          character(len=SyscStrLen), intent(in) :: varName
          character(len=SyscStrLen) :: message
          type(SyscOutputVectorDataF) :: ret
          ret = syscGetOutputVectorDataF()
          if (regName .eq. "plate" .and. 
     &        varName .eq. "Damping Force") then
            ret = syscGetOutputVectorDataCompactF(
     &        force, size(force, kind=8) / 3)
          else
            message = "getOutputVectorData error: unknown variable"
            message = trim(message) // " " // trim(varName)
            message = trim(message) // " on region "
            message = trim(message) // " " // trim(regName)
          endif
        end function getOutputVectorData
!
      end program oscillatingPlateDampingFortran
!