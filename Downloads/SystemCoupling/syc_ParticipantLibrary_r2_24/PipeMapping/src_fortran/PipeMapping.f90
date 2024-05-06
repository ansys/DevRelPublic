!
! Copyright ANSYS, Inc. Unauthorized use, distribution, 
! or duplication is prohibited.
!
! *********************************************************************
!
! This is the pipe mapping program.
!
! It generates two different meshes, each for a region representing a 
! pipe geometry. It then performs non-conformal mesh mapping between 
! these two regions.
!
! The mesh refinement level, pipe offset, pipe radius, etc. can be 
! controlled by the command line arguments. The resulting arrays can 
! be written for further post-processing.
!

program pipeMapping

use PipeMeshGenerator, only: meshProperties, initializeMesh, &
  numNodesQuad, numElemsQuad, numNodesTri, numElemsTri, &
  nodeCoordsQuad, elemNodeCountsQuad, elemNodeIdsQuad, &
  nodeCoordsTri, elemNodeCountsTri, elemNodeIdsTri, &
  linear1Quad, linear2Quad, linear1Tri, linear2Tri, &
  clearMesh

implicit none

include 'syscSystemCouplingF.fi'
include 'mpif.h'

integer(kind=4) :: write = 1

type(SyscErrorF) :: ret
type(SyscVariableF) :: variableLinear1
type(SyscVariableF) :: variableLinear2
type(SyscRegionF) :: quadRegion
type(SyscRegionF) :: triRegion
type(SyscCouplingInterfaceF) :: couplingInterface

integer :: ierror = 0

real(kind=8) :: startTime, localTime, totalTime

integer(kind=4) :: argCount
character(len=256) :: arg
character(len=256) :: argVal
integer(kind=4) :: i

! Parse command line arguments.
argCount = command_argument_count()
do i = 1, argCount
  call get_command_argument(i, arg)
  arg = trim(arg)
  if (arg .eq. '--quad' .and. i + 1 .le. argCount) then
    call get_command_argument(i + 1, argVal)
    read (argVal, *) meshProperties%quadRefine
  else if (arg .eq. '--tri' .and. i + 1 .le. argCount) then
    call get_command_argument(i + 1, argVal)
    read (argVal, *) meshProperties%triRefine
  else if (arg .eq. '--offset' .and. i + 1 .le. argCount) then
    call get_command_argument(i + 1, argVal)
    read (argVal, *) meshProperties%pipeOffset
  else if (arg .eq. '--qradius' .and. i + 1 .le. argCount) then
    call get_command_argument(i + 1, argVal)
    read (argVal, *) meshProperties%quadPipeRadius
  else if (arg .eq. '--tradius' .and. i + 1 .le. argCount) then
    call get_command_argument(i + 1, argVal)
    read (argVal, *) meshProperties%triPipeRadius
  else if (arg .eq. '--write' .and. i + 1 .le. argCount) then
    call get_command_argument(i + 1, argVal)
    read (argVal, *) write
  else if (arg .eq. '--overlap' .and. i + 1 .le. argCount) then
    call get_command_argument(i + 1, argVal)
    read (argVal, *) meshProperties%overlapLayers
  endif
enddo

! Initialize MPI.
call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, meshProperties%myrank, ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, meshProperties%numranks, ierror)

! Print setup information.
if (meshProperties%myrank .eq. 0) then
  write(*, "(A6)") 'Setup:'
  write(*, "(A17,I0)") 'Quad refinement: ', meshProperties%quadRefine
  write(*, "(A16,I0)") 'Tri refinement: ', meshProperties%triRefine
  write(*, "(A13,F5.2)") 'Pipe offset: ', meshProperties%pipeOffset
  write(*, "(A13,F5.2)") 'Quad radius: ', meshProperties%quadPipeRadius
  write(*, "(A12,F5.2)") 'Tri radius: ', meshProperties%triPipeRadius
  write(*, "(A7,I0)") 'Write: ', write
  write(*, "(A16,I0)") 'Overlap layers: ', meshProperties%overlapLayers
  write(*, "(A11,I0)") 'Processes: ', meshProperties%numranks
  write(*, "(A10)") '=========='
endif

! Add variables, regions, and the interface.
variableLinear1 = syscGetVariableF("linear1", SyscScalar, .FALSE., SyscNode)
variableLinear2 = syscGetVariableF("linear2", SyscScalar, .FALSE., SyscNode)

quadRegion = syscGetRegionF("quad", SyscSurface)

ret = syscAddOutputVariableF(quadRegion, variableLinear1)
call processRetCode(ret)

ret = syscAddInputVariableF(quadRegion, variableLinear2)
call processRetCode(ret)

triRegion = syscGetRegionF("tri", SyscSurface)
ret = syscAddOutputVariableF(triRegion, variableLinear2)
call processRetCode(ret)

ret = syscAddInputVariableF(triRegion, variableLinear1)
call processRetCode(ret)

couplingInterface = syscGetCouplingInterfaceF("interface")

ret = syscAddSideOneRegionF(couplingInterface, quadRegion)
call processRetCode(ret)

ret = syscAddSideTwoRegionF(couplingInterface, triRegion)
call processRetCode(ret)

ret = syscAddCouplingInterfaceF(couplingInterface)
call processRetCode(ret)

! Register data access functions.
ret = syscRegisterSurfMeshAccessF(getSurfaceMesh)
call processRetCode(ret)

ret = syscRegisterInputScalarDataAccessF(getInputScalarData)
call processRetCode(ret)

ret = syscRegisterOutputScalarDataAccessF(getOutputScalarData)
call processRetCode(ret)

! Initialize mesh and solution data on both sides of the interface.
startTime = MPI_WTIME()
call initializeMesh()
localTime = MPI_WTIME() - startTime
call MPI_REDUCE(localTime, totalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, ierror)
if (meshProperties%myrank .EQ. 0) then
  write(*, "(A23,ES11.4e2,A4)") 'Time to generate mesh: ', totalTime, ' [s]'
endif

! Update inputs - this will perform mapping and bring inputs up-to-date.
startTime = MPI_WTIME()
ret = syscUpdateInputsF()
localTime = MPI_WTIME() - startTime
call processRetCode(ret)
call MPI_REDUCE(localTime, totalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, ierror)
if (meshProperties%myrank .EQ. 0) then
  write(*, "(A23,ES11.4e2,A4)") 'Time to update inputs: ', totalTime, ' [s]'
endif

! Write results, if desired.
if (write .NE. 0) then
  startTime = MPI_WTIME()
  call writeArrays(quadRegion%regionName)
  call writeArrays(triRegion%regionName)
  ret = syscWriteResultsF(syscGetResultsInfoF("pipe"))
  localTime = MPI_WTIME() - startTime
  call processRetCode(ret)
  call MPI_REDUCE(localTime, totalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, ierror)
  if (meshProperties%myrank .EQ. 0) then
    write(*, "(A23,ES11.4e2,A4)") 'Time to write results: ', totalTime, ' [s]'
  endif
endif

call printMaxError(quadRegion%regionName)
call printMaxError(triRegion%regionName)

call clearMesh()

call MPI_FINALIZE(ierror)

contains

! Function that returns surface mesh, given a region.
function getSurfaceMesh(regName) result(ret)
  character(len=SyscStrLen), intent(in) :: regName
  type(SyscOutputVectorDataF) :: nodeCoordsData
  type(SyscOutputIntegerDataF) :: elemNodeCountsData
  type(SyscOutputIntegerDataF) :: elemNodeIdsData
  type(SyscNodeDataF) :: nodes
  type(SyscFaceDataF) :: faces
  character(len=SyscStrLen) :: message
  type(SyscSurfaceMeshF) :: ret

  ret = syscGetSurfaceMeshF()
  if (regName .eq. "quad") then
    ! node coords
    nodeCoordsData = syscGetOutputVectorDataCompactF(&
      nodeCoordsQuad, size(nodeCoordsQuad, kind=8) / 3)
    ! elem node counts
    elemNodeCountsData = syscGetOutputIntegerDataF(&
      elemNodeCountsQuad, size(elemNodeCountsQuad, kind=8))
    ! elem node ids
    elemNodeIdsData = syscGetOutputIntegerDataF(&
      elemNodeIdsQuad, size(elemNodeIdsQuad, kind=8))
  else if (regName .eq. "tri") then
    ! node coords
    nodeCoordsData = syscGetOutputVectorDataCompactF(&
      nodeCoordsTri, size(nodeCoordsTri, kind=8) / 3)
    ! elem node counts
    elemNodeCountsData = syscGetOutputIntegerDataF(&
      elemNodeCountsTri, size(elemNodeCountsTri, kind=8))
    ! elem node ids
    elemNodeIdsData = syscGetOutputIntegerDataF(&
      elemNodeIdsTri, size(elemNodeIdsTri, kind=8))
  else
    message = "getSurfaceMesh error: unknown region"
    message = trim(message) // " " // trim(regName)
    call syscFatalErrorF(message)
  endif
  ! nodes and faces
  nodes = syscGetNodeDataF(nodeCoordsData)
  faces = syscGetFaceDataF(&
    syscGetElementNodeCountDataF(elemNodeCountsData),&
    syscGetElementNodeConnectivityDataF(elemNodeIdsData))
  ! surface mesh
  ret = syscGetSurfaceMeshF(nodes, faces)
end function getSurfaceMesh

! Function that returns scalar input variable values,
! given a region and a variable.
function getInputScalarData(regName, varName) result(ret)
  character(len=SyscStrLen), intent(in) :: regName
  character(len=SyscStrLen), intent(in) :: varName
  character(len=SyscStrLen) :: message
  type(SyscInputScalarDataF) :: ret
  ret = syscGetInputScalarDataF()
  if (regName .eq. "quad" .and. varName .eq. "linear2") then
    ret = syscGetInputScalarDataF(linear2Quad, size(linear2Quad, kind=8))
  else if (regName .eq. "tri" .and. varName .eq. "linear1") then
    ret = syscGetInputScalarDataF(linear1Tri, size(linear1Tri, kind=8))
  else
    message = "getInputScalarData: incorrect region or variable name"
    call syscFatalErrorF(message)
  endif
end function getInputScalarData

! Function that returns scalar output variable values,
! given a region and a variable.
function getOutputScalarData(regName, varName) result(ret)
  character(len=SyscStrLen), intent(in) :: regName
  character(len=SyscStrLen), intent(in) :: varName
  character(len=SyscStrLen) :: message
  type(SyscOutputScalarDataF) :: ret
  ret = syscGetOutputScalarDataF()
  if (regName .eq. "quad" .and. varName .eq. "linear1") then
    ret = syscGetOutputScalarDataF(linear1Quad, size(linear1Quad, kind=8))
  else if (regName .eq. "tri" .and. varName .eq. "linear2") then
    ret = syscGetOutputScalarDataF(linear2Tri, size(linear2Tri, kind=8))
  else
    message = "getOutputScalarData: incorrect region or variable name"
    call syscFatalErrorF(message)
  endif
end function getOutputScalarData

subroutine processRetCode(ret)
  type(SyscErrorF), intent(in) :: ret
  if (ret%retcode .NE. SyscStatusOk) then
    write(*, *) "ERROR: ", trim(ret%message)
    call MPI_FINALIZE(ierror)
    stop 1
  endif
end subroutine processRetCode

! Function to write out mesh and solution data arrays to text files.
subroutine writeArrays(regionName)
  character(len=SyscStrLen), intent(in) :: regionName
  character(len=SyscStrLen) :: fileName
  integer(kind=4) :: numScalarNodeArrays = 2, numNodes, numElems, i,&
    elemNodeIndex, e, en
  real(kind=8), pointer :: nodeCoords(:), linear1Data(:), linear2Data(:)
  integer(kind=4), pointer :: elemNodeCounts(:), elemNodeIds(:)
  if (regionName .EQ. 'quad') then
    numNodes = numNodesQuad
    numElems = numElemsQuad
    nodeCoords => nodeCoordsQuad
    elemNodeCounts => elemNodeCountsQuad
    elemNodeIds => elemNodeIdsQuad
    linear1Data => linear1Quad
    linear2Data => linear2Quad
  else
    numNodes = numNodesTri
    numElems = numElemsTri
    nodeCoords => nodeCoordsTri
    elemNodeCounts => elemNodeCountsTri
    elemNodeIds => elemNodeIdsTri
    linear1Data => linear1Tri
    linear2Data => linear2Tri
  endif

  write(fileName, "(A8,A1,I0,A4)") trim(regionName), '.', meshProperties%myrank, '.dat'
  open(1, file = trim(adjustl(fileName)), status = 'replace')

  ! Write number of nodes, elements, and solution data arrays.
  write(1, "(I0,A2,I0,A2,I0)", advance = 'yes') numNodes, ', ', numElems, ', ', numScalarNodeArrays

  ! Write node coordinates.
  do i = 1, numNodes
    write(1, "(ES15.8e2,A2,ES15.8e2,A2,ES15.8e2)") &
      nodeCoords((i - 1) * 3 + 1), ', ',&
      nodeCoords((i - 1) * 3 + 2), ', ',&
      nodeCoords((i - 1) * 3 + 3)
  enddo

  ! Write element-to-node connectivity.
  elemNodeIndex = 1
  do e = 1, numElems
    ! Write the first node id for the current element.
    write(1, "(I0)", advance = 'no') elemNodeIds(elemNodeIndex)
    elemNodeIndex = elemNodeIndex + 1
    ! Write subsequent node ids for the current element.
    do en = 2, elemNodeCounts(e)
      write(1, "(A2,I0)", advance = 'no') ', ', elemNodeIds(elemNodeIndex)
      elemNodeIndex = elemNodeIndex + 1
    enddo
    write(1, "()", advance = 'yes')
  enddo

  ! Write solution data arrays.
  write(1, "(A7)") "linear1"
  do i = 1, numNodes
    write(1, "(ES15.8e2)") linear1Data(i)
  enddo

  write(1, "(A7)") "linear2"
  do i = 1, numNodes
    write(1, "(ES15.8e2)") linear2Data(i)
  enddo

  close(1)

end subroutine writeArrays

! Function to calculate and print max error on a given region.
subroutine printMaxError(regionName)
  character(len=SyscStrLen), intent(in) :: regionName
  real(kind=8) :: maxError, totalMaxError, testValue, expectedValue,&
    currDiff
  real(kind=8), pointer :: testValues(:), expectedValues(:)
  integer(kind=4) :: n, numVals

  if (regionName .EQ. 'quad') then
    testValues => linear2Quad
    expectedValues => linear1Quad
    numVals = numNodesQuad
  else
    testValues => linear1Tri
    expectedValues => linear2Tri
    numVals = numNodesTri
  endif

  do n = 1, numVals
    testValue = testValues(n)
    expectedValue = expectedValues(n)
    currDiff = abs(expectedValue - testValue)
    maxError = max(currDiff, maxError)
  enddo

  call MPI_REDUCE(maxError, totalMaxError, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, ierror)

  if (meshProperties%myrank .EQ. 0) then
    write(*, "(A17,A,A2,ES11.4e2)") 'Maximum error on ', trim(regionName), ': ', totalMaxError
  endif

end subroutine printMaxError

end program pipeMapping
