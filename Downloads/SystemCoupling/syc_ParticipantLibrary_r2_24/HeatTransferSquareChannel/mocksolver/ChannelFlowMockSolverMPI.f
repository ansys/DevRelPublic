!
! Copyright ANSYS, Inc. Unauthorized use, distribution, 
! or duplication is prohibited.
!
! *********************************************************************
!
! This is the parallel implementation of the mock solver program.
!
! It takes the following command line arguments to get
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
      include 'mpif.h'
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
      integer :: myrank = 0
      integer :: numranks = 1
      integer :: ierror = 0
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
!           | A  | B  | C  | D  |
!           |    |    |    |    |
!           1----2----3----4----5
!  (0,0,0)                       (1,0,0)
!
! In parallel, the above mesh will be partitioned as follows:
!
! numranks : 1,  elements on rank 0 : A, B, C, D
! numranks : 2,  elements on rank 0 : A, B
!                elements on rank 1 : C, D
! numranks : 3,  elements on rank 0 : A
!                elements on rank 1 : B
!                elements on rank 2 : C, D
! numranks : 4+, elements on rank 0 : A
!                elements on rank 1 : B
!                elements on rank 2 : C
!                elements on rank 3 : D
!                elements on rank 4+: no elements
!
! *********************************************************************
!
!     Declare arrays used for this problem.
!     The contents are populated in the initializeArray()
!     subroutine, depending on parallel settings.
!
!     Array of node coordinates in split array format.
      real(kind=8), target, allocatable :: nodeCoordsX(:)
      real(kind=8), target, allocatable :: nodeCoordsY(:)
      real(kind=8), target, allocatable :: nodeCoordsZ(:)
!
!     Array of element node counts.
      integer(kind=4), target, allocatable :: elemNodeCounts(:)
!
!     Array of node ids for each element.
      integer(kind=4), target, allocatable :: elemNodeIds(:)
!
!     Array of temperature values (defined at nodes).
      real(kind=8), target, allocatable :: temperature(:)
!
!     Array of heat rate values (defined at elements)
!     values are in Watts - total of 20 [W] generated.
      real(kind=8), target, allocatable :: heatRate(:)
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
      buildInfo = 'Channel Flow Mock Solver MPI, language: FORTRAN'
!
!     Initialize MPI.
      call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numranks, ierror)
!
!     Connection to System Coupling.
      ret = syscConnectF(scHost, scPort, scName, 
     &  MPI_COMM_WORLD, buildInfo)
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
      ret = syscRegisterOutputScalarDataAccessF(
     &  getOutputScalarData)
      call processRetCode(ret)
!
!     Initialize the arrays.
      call initializeArrays()
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
      deallocate(nodeCoordsX)
      deallocate(nodeCoordsY)
      deallocate(nodeCoordsZ)
      deallocate(elemNodeCounts)
      deallocate(elemNodeIds)
      deallocate(temperature)
      deallocate(heatRate)
!
!     Coupled analysis shutdown.
      ret = syscDisconnectF()
      call processRetCode(ret)
!
      call MPI_FINALIZE(ierror)
!
      contains
!
!       Process error code.
        subroutine processRetCode(ret)
          type(SyscErrorF) :: ret
          if (ret%retcode .NE. SyscStatusOk) then
            print *, "ERROR: ", trim(ret%message)
            call MPI_FINALIZE(ierror)
            stop 1
          endif
        end subroutine processRetCode
!
!       Subroutine that partitions the mesh based
!       on the number of parallel processes, and
!       initializes the arrays.
!
        subroutine initializeArrays()
          if (numranks .eq. 1) then
!           All four faces are present on rank 0.
            allocate(nodeCoordsX(10))
            nodeCoordsX(:) = (/ 0.0D0, 0.25D0, 0.5D0, 0.75D0, 1.0D0,
     &                          0.0D0, 0.25D0, 0.5D0, 0.75D0, 1.0D0 /)
            allocate(nodeCoordsY(10))
            nodeCoordsY(:) = (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,
     &                          0.1D0, 0.1D0, 0.1D0, 0.1D0, 0.1D0 /)
            allocate(nodeCoordsZ(10))
            nodeCoordsZ(:) = 0.0D0
            allocate(elemNodeCounts(4))
            elemNodeCounts(:) = 4
            allocate(elemNodeIds(16))
            elemNodeIds(:) = (/ 1, 2, 7, 6, 2, 3, 8, 7, 
     &                          3, 4, 9, 8, 4, 5, 10, 9/)
          else if (numranks .eq. 2) then
            if (myrank .eq. 0) then
!             Faces A and B are on rank 0.
              allocate(nodeCoordsX(6))
              nodeCoordsX(:) = (/ 0.0D0, 0.25D0, 0.50D0, 
     &                            0.0D0, 0.25D0, 0.50D0 /)
              allocate(nodeCoordsY(6))
              nodeCoordsY(:) = (/ 0.0D0, 0.0D0, 0.0D0,
     &                            0.1D0, 0.1D0, 0.1D0/)
              allocate(nodeCoordsZ(6))
              nodeCoordsZ(:) = 0.0D0
              allocate(elemNodeCounts(2))
              elemNodeCounts(:) = 4
              allocate(elemNodeIds(8))
              elemNodeIds(:) = (/ 1, 2, 7, 6, 2, 3, 8, 7 /)
            else
!             Faces C and D are on rank 1.
              allocate(nodeCoordsX(6))
              nodeCoordsX(:) = (/ 0.5D0, 0.75D0, 1.0D0, 
     &                            0.5D0, 0.75D0, 1.0D0 /)
              allocate(nodeCoordsY(6))
              nodeCoordsY(:) = (/ 0.0D0, 0.0D0, 0.0D0,
     &                            0.1D0, 0.1D0, 0.1D0/)
              allocate(nodeCoordsZ(6))
              nodeCoordsZ(:) = 0.0D0
              allocate(elemNodeCounts(2))
              elemNodeCounts(:) = 4
              allocate(elemNodeIds(8))
              elemNodeIds(:) = (/ 3, 4, 9, 8, 4, 5, 10, 9 /)
            endif
          else if (numranks .eq. 3) then
            if (myrank .eq. 0) then
!             Face A is on rank 0.
              allocate(nodeCoordsX(4))
              nodeCoordsX(:) = (/ 0.0D0, 0.25D0, 0.0D0, 0.25D0 /)
              allocate(nodeCoordsY(4))
              nodeCoordsY(:) = (/ 0.0D0, 0.0D0, 0.1D0, 0.1D0 /)
              allocate(nodeCoordsZ(4))
              nodeCoordsZ(:) = 0.0D0
              allocate(elemNodeCounts(1))
              elemNodeCounts(:) = 4
              allocate(elemNodeIds(4))
              elemNodeIds(:) = (/ 1, 2, 7, 6 /)
            else if (myrank .eq. 1) then
!             Face B is on rank 1.
              allocate(nodeCoordsX(4))
              nodeCoordsX(:) = (/ 0.25D0, 0.5D0, 0.25D0, 0.5D0 /)
              allocate(nodeCoordsY(4))
              nodeCoordsY(:) = (/ 0.0D0, 0.0D0, 0.1D0, 0.1D0 /)
              allocate(nodeCoordsZ(4))
              nodeCoordsZ(:) = 0.0D0
              allocate(elemNodeCounts(1))
              elemNodeCounts(:) = 4
              allocate(elemNodeIds(4))
              elemNodeIds(:) = (/ 2, 3, 8, 7 /)
            else
!             Faces C and D are on rank 2.
              allocate(nodeCoordsX(6))
              nodeCoordsX(:) = (/ 0.5D0, 0.75D0, 1.0D0, 
     &                            0.5D0, 0.75D0, 1.0D0 /)
              allocate(nodeCoordsY(6))
              nodeCoordsY(:) = (/ 0.0D0, 0.0D0, 0.0D0,
     &                            0.1D0, 0.1D0, 0.1D0/)
              allocate(nodeCoordsZ(6))
              nodeCoordsZ(:) = 0.0D0
              allocate(elemNodeCounts(2))
              elemNodeCounts(:) = 4
              allocate(elemNodeIds(8))
              elemNodeIds(:) = (/ 3, 4, 9, 8, 4, 5, 10, 9 /)
            endif
          else if (numranks .ge. 4) then
!           Four or more processes.
!           Ranks 0, 1, 2, 3 get faces A, B, C, D, respectively.
!           Any higher ranks (if exist) have to mesh.
            if (myrank .eq. 0) then
!             Face A is on rank 0.
              allocate(nodeCoordsX(4))
              nodeCoordsX(:) = (/ 0.0D0, 0.25D0, 0.0D0, 0.25D0 /)
              allocate(nodeCoordsY(4))
              nodeCoordsY(:) = (/ 0.0D0, 0.0D0, 0.1D0, 0.1D0 /)
              allocate(nodeCoordsZ(4))
              nodeCoordsZ(:) = 0.0D0
              allocate(elemNodeCounts(1))
              elemNodeCounts(:) = 4
              allocate(elemNodeIds(4))
              elemNodeIds(:) = (/ 1, 2, 7, 6 /)
            else if (myrank .eq. 1) then
!             Face B is on rank 1.
              allocate(nodeCoordsX(4))
              nodeCoordsX(:) = (/ 0.25D0, 0.5D0, 0.25D0, 0.5D0 /)
              allocate(nodeCoordsY(4))
              nodeCoordsY(:) = (/ 0.0D0, 0.0D0, 0.1D0, 0.1D0 /)
              allocate(nodeCoordsZ(4))
              nodeCoordsZ(:) = 0.0D0
              allocate(elemNodeCounts(1))
              elemNodeCounts(:) = 4
              allocate(elemNodeIds(4))
              elemNodeIds(:) = (/ 2, 3, 8, 7 /)
            else if (myrank .eq. 2) then
!             Face C is on rank 2.
              allocate(nodeCoordsX(4))
              nodeCoordsX(:) = (/ 0.5D0, 0.75D0, 0.5D0, 0.75D0 /)
              allocate(nodeCoordsY(4))
              nodeCoordsY(:) = (/ 0.0D0, 0.0D0, 0.1D0, 0.1D0 /)
              allocate(nodeCoordsZ(4))
              nodeCoordsZ(:) = 0.0D0
              allocate(elemNodeCounts(1))
              elemNodeCounts(:) = 4
              allocate(elemNodeIds(4))
              elemNodeIds(:) = (/ 3, 4, 9, 8 /)
            else if (myrank .eq. 3) then
!             Face D is on rank 3.
              allocate(nodeCoordsX(4))
              nodeCoordsX(:) = (/ 0.75D0, 1.0D0, 0.75D0, 1.0D0 /)
              allocate(nodeCoordsY(4))
              nodeCoordsY(:) = (/ 0.0D0, 0.0D0, 0.1D0, 0.1D0 /)
              allocate(nodeCoordsZ(4))
              nodeCoordsZ(:) = 0.0D0
              allocate(elemNodeCounts(1))
              elemNodeCounts(:) = 4
              allocate(elemNodeIds(4))
              elemNodeIds(:) = (/ 4, 5, 10, 9 /)
            else
!             All arrays remain empty for ranks 4 and higher.
              allocate(nodeCoordsX(0))
              allocate(nodeCoordsY(0))
              allocate(nodeCoordsZ(0))
              allocate(elemNodeCounts(0))
              allocate(elemNodeIds(0))
            endif
          endif
!         Resize element heat rate to however many elements there are,
!         initialize all entries to 5.0 [W].
          allocate(heatRate(size(elemNodeCounts)))
          heatRate(:) = 5.0D0
!         Resize node temperature to however many nodes there are,
!         initialize all entries to 300.0 [K].
          allocate(temperature(size(nodeCoordsX)))
          temperature(:) = 300.0D0
        end subroutine initializeArrays
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