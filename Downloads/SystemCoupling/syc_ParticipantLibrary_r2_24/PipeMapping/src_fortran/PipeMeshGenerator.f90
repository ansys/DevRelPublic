!
! Copyright ANSYS, Inc. Unauthorized use, distribution, 
! or duplication is prohibited.
!
! *********************************************************************

module PipeMeshGenerator

implicit none
private

type, public :: PipeMeshProperties
  integer(kind=4) :: quadRefine = 20
  integer(kind=4) :: triRefine = 30
  real(kind=8) :: pipeOffset = 0.0D0
  real(kind=8) :: quadPipeRadius = 0.05D0
  real(kind=8) :: triPipeRadius = 0.05D0
  integer(kind=4) :: overlapLayers = 0
  real(kind=8) :: pipeLength = 1.0D0
  integer(kind=4) :: myrank = 0
  integer(kind=4) :: numranks = 1
end type PipeMeshProperties

! Define structure to hold partition bounds.
type :: PartitionBounds
  real(kind=8) :: axialStart
  real(kind=8) :: axialEnd
  real(kind=8) :: circStartAngle
  real(kind=8) :: circEndAngle
end type PartitionBounds

! Define useful constants.
real(kind=8) :: pi = 3.14159265358979323846D0

type(PipeMeshProperties), public, save :: meshProperties

! Define data structures for mesh arrays.
real(kind=8), target, allocatable, public :: nodeCoordsQuad(:)
real(kind=8), target, allocatable, public :: nodeCoordsTri(:)
integer(kind=4), target, allocatable, public :: elemNodeCountsQuad(:)
integer(kind=4), target, allocatable, public :: elemNodeCountsTri(:)
integer(kind=4), target, allocatable, public :: elemNodeIdsQuad(:)
integer(kind=4), target, allocatable, public :: elemNodeIdsTri(:)

integer(kind=4), public :: numNodesQuad = 0
integer(kind=4), public :: numNodesTri = 0
integer(kind=4), public :: numElemsQuad = 0
integer(kind=4), public :: numElemsTri = 0

! Define data structures to hold solution data.
real(kind=8), target, allocatable, public :: linear1Quad(:)
real(kind=8), target, allocatable, public :: linear1Tri(:)
real(kind=8), target, allocatable, public :: linear2Quad(:)
real(kind=8), target, allocatable, public :: linear2Tri(:)

integer(kind=4), save :: nodeIdOffset = 0

public :: initializeMesh
public :: clearMesh

contains
!
subroutine clearMesh()
  if (numElemsQuad .GT. 0) then
    deallocate(nodeCoordsQuad)
    deallocate(elemNodeCountsQuad)
    deallocate(elemNodeIdsQuad)
    deallocate(linear1Quad)
    deallocate(linear2Quad)
  endif
  
  if (numElemsTri .GT. 0) then
    deallocate(nodeCoordsTri)
    deallocate(elemNodeCountsTri)
    deallocate(elemNodeIdsTri)
    deallocate(linear1Tri)
    deallocate(linear2Tri)
  endif
end subroutine clearMesh
!
! Recursive function to calculate partition bounds for a patch of a
! pipe.
! 
! Arguments are as follows:
!
! myrank - current rank 
! 
! partitions - number of partitions to break the current patch into 
!
! axialStart - start Z-coordinate of a patch to partition
!
! axialEnd - end Z-coordinate of a patch to partition 
!
! circStart - start angle (in radians) of a patch to partition 
!
! circEnd - end angle (in radians) of a patch to partition
!
! For example, patch can be defined as [0.5, 1.0, 0, 2.0 * PI]. This
! defines a patch of the pipe from 0.5 to 1.0 in axial direction and full
! circle in circumferential direction.
recursive function getPartitionBounds(&
  rank, partitions, &
  axialStart, axialEnd, &
  circStart, circEnd, &
  pipeRadius) result(bounds)
  integer(kind=4), intent(in) :: rank, partitions
  real(kind=8), intent(in) :: axialStart, axialEnd, circStart, &
    circEnd, pipeRadius
  real(kind=8) :: newAxialStart, newAxialEnd, newCircStart, newCircEnd
  real(kind=8) :: axialLen, circLen
  integer :: newRank, newPartitions, splitBucket, splitDir
  type(PartitionBounds) :: bounds

  if (partitions .eq. 1) then
    ! Base case for the recursive function.
    bounds%axialStart = axialStart
    bounds%axialEnd = axialEnd
    bounds%circStartAngle = circStart
    bounds%circEndAngle = circEnd
  else

    axialLen = axialEnd - axialStart
    circLen = (circEnd - circStart) * pipeRadius

    ! Determine if the current rank will be in upper (1) or lower (0) division.
    if (rank .LT. (partitions / 2)) then
      splitBucket = 0
    else
      splitBucket = 1
    endif

    if (splitBucket .EQ. 0) then
      newRank = rank
      newPartitions = partitions / 2
    else
      newRank = rank - partitions / 2
      newPartitions = partitions / 2 + modulo(partitions, 2)
    endif

    ! Determine the split direction (0 is axial, 1 is circumferential).
    if (axialLen .GE. circLen) then
      splitDir = 0
    else
      splitDir = 1
    endif

    newAxialStart = axialStart
    newAxialEnd = axialEnd
    newCircStart = circStart
    newCircEnd = circEnd
    if (splitDir .EQ. 0) then
      if (splitBucket .EQ. 0) then
        newAxialEnd = axialStart + 0.5 * axialLen
      else
        newAxialStart = axialStart + 0.5 * axialLen
      endif
    else
      if (splitBucket .EQ. 0) then
        newCircEnd = circStart + 0.5 * (circEnd - circStart)
      else
        newCircStart = circStart + 0.5 * (circEnd - circStart)
      endif
    endif

    bounds = getPartitionBounds(&
      newRank, newPartitions, &
      newAxialStart, newAxialEnd, &
      newCircStart, newCircEnd, &
      pipeRadius)
  endif

end function getPartitionBounds
!
subroutine addBeginning(array, val)
  integer(kind=4), dimension(:), allocatable, intent(inout) :: array
  integer(kind=4), intent(in) :: val
  integer(kind=4), dimension(:), allocatable :: newArray
  integer(kind=4) :: i, oldSize

  oldSize = size(array, kind=4)
  allocate(newArray(oldSize + 1))
  do i = 1, oldSize          
    newArray(i + 1) = array(i)
  end do
  newArray(1) = val
  deallocate(array)
  call move_alloc(newArray, array)

end subroutine addBeginning

subroutine addEnd(array, val)
  integer(kind=4), allocatable, intent(inout) :: array(:)
  integer(kind=4), intent(in) :: val
  integer(kind=4), allocatable :: newArray(:)
  integer(kind=4) :: i, oldSize

  oldSize = size(array, kind=4)
  allocate(newArray(oldSize + 1))
  do i = 1, oldSize          
    newArray(i) = array(i)
  end do
  newArray(oldSize + 1) = val
  deallocate(array)
  call move_alloc(newArray, array)

end subroutine addEnd

! Function to add overlap layers when running in parallel.
subroutine addOverlapLayers(&
  circumferenceElems, axialEdges,&
  elemAxialRange, elemCircRange)
  integer(kind=4), intent(in) :: circumferenceElems, axialEdges
  integer(kind=4), allocatable, intent(inout) :: elemAxialRange(:),&
    elemCircRange(:)
  integer(kind=4) :: overlap, firstAxialRow, lastAxialRow, e,&
    firstCircRow, prevCircRow, lastCircRow, nextCircRow
  integer(kind=4), allocatable :: circElemFlags(:)

  ! Add overlap layers along axial direction.
  if (size(elemAxialRange) .GT. 0) then
    do overlap = 1, meshProperties%overlapLayers
      firstAxialRow = elemAxialRange(1)
      lastAxialRow = elemAxialRange(size(elemAxialRange))
      if (firstAxialRow .GT. 0) then
        call addBeginning(elemAxialRange, firstAxialRow - 1)
      endif
      if (lastAxialRow .LT. axialEdges - 1) then
        call addEnd(elemAxialRange, lastAxialRow + 1)
      endif
    enddo
  endif

  ! Add overlap layers along circumferential direction.
  if (size(elemCircRange) .GT. 0) then
    allocate(circElemFlags(circumferenceElems))
    circElemFlags(:) = 0

    ! Mark elements already in the elemCircRange.
    do e = 1, size(elemCircRange)
      circElemFlags(elemCircRange(e) + 1) = 1
    enddo

    ! Mark elements in the overlap layer.
    do overlap = 1, meshProperties%overlapLayers
      ! Add element at the beginning.
      firstCircRow = elemCircRange(1)
      prevCircRow = modulo(firstCircRow + circumferenceElems - overlap, circumferenceElems)
      circElemFlags(prevCircRow + 1) = 1
      ! Add element at the end.
      lastCircRow = elemCircRange(size(elemCircRange))
      nextCircRow = modulo(lastCircRow + overlap, circumferenceElems)
      circElemFlags(nextCircRow + 1) = 1
    enddo

    ! Clear and re-populate elemCircRange
    deallocate(elemCircRange)
    allocate(elemCircRange(0))
    do e = 1, circumferenceElems
      if (circElemFlags(e) .EQ. 1) then
        call addEnd(elemCircRange, e - 1)
      endif
    enddo

    deallocate(circElemFlags)
  endif
end subroutine addOverlapLayers

! Helper function to be used by initializeTri() function -
! get range of nodes along circumference (local to this partition). 
function getNodeCircRange(&
  axialRow, elemCircRange, circumferenceEdges) result(nodeCircRange)
  integer(kind=4), intent(in) :: axialRow, circumferenceEdges
  integer(kind=4), allocatable, intent(in) :: elemCircRange(:)
  integer(kind=4), allocatable :: nodeCircRange(:)
  integer(kind=4), allocatable :: circNodeFlags(:)
  integer(kind=4) :: eci, e, n

  allocate(circNodeFlags(circumferenceEdges))
  circNodeFlags(:) = 0

  do eci = 1, size(elemCircRange)
    e = elemCircRange(eci)
    if (.NOT. (modulo(axialRow, 2) .EQ. 1 .AND. modulo(e, 2) .EQ. 1)) then
      circNodeFlags(modulo(e / 2, circumferenceEdges) + 1) = 1
    endif
    if (.NOT. (modulo(axialRow, 2) .EQ. 0 .AND. modulo(e, 2) .EQ. 0)) then
      circNodeFlags(modulo(e / 2 + 1, circumferenceEdges) + 1) = 1
    endif
  enddo

  allocate(nodeCircRange(0))
  do n = 1, circumferenceEdges
    if (circNodeFlags(n) .EQ. 1) then
      call addEnd(nodeCircRange, n - 1)
    endif
  enddo

  deallocate(circNodeFlags)

end function getNodeCircRange

! Initializes mesh on "quad" region (quad mesh).
!
! Here is the schematic of the mesh with 4 cirumference edges and 2 axial
! edges.
! 
!    circumrefential direction
!              ^
!              |
!              |
!              |
!              x------> axial direction
! 
! 0--------4--------8
! |        |        |
! |        |        |
! |        |        |
! 3--------7-------11
! |        |        |
! |        |        |
! |        |        |
! 2--------6-------10
! |        |        |
! |        |        |
! |        |        |
! 1--------5--------9
! |        |        |
! |        |        |
! |        |        |
! 0--------4--------8
!
subroutine initializeQuad()
  integer(kind=4) :: circumferenceEdges, axialEdges
  integer(kind=4) :: totalNumNodes
  integer(kind=4) :: nodesPerFace = 4
  real(kind=8) :: axialEdgeLen, circEdgeAngle
  type(PartitionBounds) :: partBounds
  integer(kind=4), allocatable :: elemAxialRange(:), elemCircRange(:), &
    nodeAxialRange(:), nodeCircFlags(:), nodeCircRange(:)
  integer(kind=4) :: e, c, nodesAlongCirc, nodesAlongAxial, ec, n
  real(kind=8) :: midpointAxial, midpointCirc
  integer(kind=4) :: nai, na, nci, nc, nodeIndex, &
    elemIndex, eai, ea, baseNode, nextCirc, nextAxialNextCirc, nextAxial, &
    eci, i
  real(kind=8) :: angle, x, y, z

  circumferenceEdges = meshProperties%quadRefine
  axialEdges = meshProperties%quadRefine * 3

  axialEdgeLen = meshProperties%pipeLength / real(axialEdges, kind=8)
  circEdgeAngle = 2.0D0 * pi / real(circumferenceEdges, kind=8)

  ! Total number of nodes in this region (across all partitions).
  totalNumNodes = circumferenceEdges * (axialEdges + 1)

  ! Get the partition bounds.
  partBounds = getPartitionBounds(&
    meshProperties%myrank, meshProperties%numranks, &
    meshProperties%pipeOffset, &
    meshProperties%pipeOffset + meshProperties%pipeLength, &
    0.0D0, 2.0D0 * pi, meshProperties%quadPipeRadius)

  allocate(elemAxialRange(0))
  do e = 1, axialEdges
    midpointAxial = meshProperties%pipeOffset + 0.5 * axialEdgeLen + (e - 1) * axialEdgeLen
    if (midpointAxial .GE. partBounds%axialStart .AND.&
        midpointAxial .LT. partBounds%axialEnd) then
      call addEnd(elemAxialRange, e - 1)
    endif
  enddo

  allocate(elemCircRange(0))
  do c = 1, circumferenceEdges
    midpointCirc = 0.5 * circEdgeAngle + (c - 1) * circEdgeAngle
    if (midpointCirc .GE. partBounds%circStartAngle .AND.&
        midpointCirc .LT. partBounds%circEndAngle) then
      call addEnd(elemCircRange, c - 1)
    endif
  enddo

  ! Add overlap layers.
  call addOverlapLayers(circumferenceEdges, axialEdges, elemAxialRange, elemCircRange)

  ! Calculate the number of elements and nodes in the current partition.
  numElemsQuad = size(elemCircRange) * size(elemAxialRange)
  numNodesQuad = 0
  if (numElemsQuad .GT. 0) then
    if (size(elemCircRange) .LT. circumferenceEdges) then
      nodesAlongCirc = size(elemCircRange) + 1
    else
      nodesAlongCirc = size(elemCircRange)
    endif
    nodesAlongAxial = size(elemAxialRange) + 1
    numNodesQuad = nodesAlongCirc * nodesAlongAxial
  endif

  ! Allocate and initialize mesh and solution data arrays.
  allocate(nodeCoordsQuad(numNodesQuad * 3))
  allocate(elemNodeCountsQuad(numElemsQuad))
  elemNodeCountsQuad(:) = nodesPerFace
  allocate(elemNodeIdsQuad(numElemsQuad * nodesPerFace))
  allocate(linear1Quad(numNodesQuad))
  allocate(linear2Quad(numNodesQuad))

  ! If there are no elements in the current partition, return.
  if (numElemsQuad .GT. 0) then

    ! Calculate axial node range
    ! (same as element axial range plus one row of nodes at the end).
    allocate(nodeAxialRange(size(elemAxialRange)))
    nodeAxialRange(:) = elemAxialRange(:)
    do i = 1, size(elemAxialRange)
      !nodeAxialRange(i) = elemAxialRange(i)
    enddo
    call addEnd(nodeAxialRange, elemAxialRange(size(elemAxialRange)) + 1)

    ! Calculate circumferential node range.
    allocate(nodeCircFlags(circumferenceEdges))
    nodeCircFlags(:) = 0
    do i = 1, size(elemCircRange)
      ec = elemCircRange(i)
      nodeCircFlags(ec + 1) = 1
      if (ec + 2 .GT. circumferenceEdges) then
        nodeCircFlags(1) = 1
      else
        nodeCircFlags(ec + 2) = 1
      endif
    enddo

    allocate(nodeCircRange(0))
    do n = 1, size(nodeCircFlags)
      if (nodeCircFlags(n) .EQ. 1) then
        call addEnd(nodeCircRange, n - 1)
      endif
    enddo

    ! Fill node coordinates.
    nodeIndex = 0
    do nai = 1, size(nodeAxialRange)
      na = nodeAxialRange(nai)
      do nci = 1, size(nodeCircRange)
        nc = nodeCircRange(nci)
        angle = 2.0D0 * pi * nc / real(circumferenceEdges)
        x = cos(angle) * meshProperties%quadPipeRadius
        y = sin(angle) * meshProperties%quadPipeRadius
        z = na * axialEdgeLen + meshProperties%pipeOffset
        nodeCoordsQuad(nodeIndex * 3 + 1) = x
        nodeCoordsQuad(nodeIndex * 3 + 2) = y
        nodeCoordsQuad(nodeIndex * 3 + 3) = z
        nodeIndex = nodeIndex + 1
      enddo
    enddo

    ! Fill element to node connectivity.
    elemIndex = 0
    do eai = 1, size(elemAxialRange)
      ea = elemAxialRange(eai)
      do eci = 1, size(elemCircRange)
        ec = elemCircRange(eci)
        ! Base node.
        baseNode = ea * circumferenceEdges + ec
        ! Get next node along circumference.
        nextCirc = baseNode + 1
        if (ec + 1 .GE. circumferenceEdges) then
          nextCirc = ea * circumferenceEdges
        endif
        ! Get next node along axial direction, relative to nextCirc node.
        nextAxialNextCirc = nextCirc + circumferenceEdges
        ! Get next node along axial direction, relative to the baseNode.
        nextAxial = baseNode + circumferenceEdges
        ! Fill element node ids array, taking node id offset into account.
        elemNodeIdsQuad(elemIndex * nodesPerFace + 1) = baseNode + nodeIdOffset
        elemNodeIdsQuad(elemIndex * nodesPerFace + 2) = nextCirc + nodeIdOffset
        elemNodeIdsQuad(elemIndex * nodesPerFace + 3) = nextAxialNextCirc + nodeIdOffset
        elemNodeIdsQuad(elemIndex * nodesPerFace + 4) = nextAxial + nodeIdOffset
        elemIndex = elemIndex + 1
      enddo
    enddo

    deallocate(nodeAxialRange)
    deallocate(nodeCircFlags)
    deallocate(nodeCircRange)

    ! Fill solution data for linear1 variable with the linear profile.
    do n = 1, numNodesQuad
      x = nodeCoordsQuad((n - 1) * 3 + 1)
      y = nodeCoordsQuad((n - 1) * 3 + 2)
      z = nodeCoordsQuad((n - 1) * 3 + 3)
      linear1Quad(n) = 1.0D0 * x + 2.0D0 * y + 3.0D0 * z + 4.0D0
      linear2Quad(n) = 0.0D0 ! TODO: remove
    enddo

  endif

  deallocate(elemAxialRange)
  deallocate(elemCircRange)

  nodeIdOffset = totalNumNodes

end subroutine initializeQuad

! Initializes mesh on "tri" region (tri mesh).
!
! Here is the schematic of the mesh with 4 cirumference edges and 2 axial
! edges.
!
!    circumrefential direction
!              ^
!              |
!              |
!              |
!              x------> axial direction
!
! 0                 0
! |  \     |     /  |
! |     \  |  /     |
! |        0        |
! |     /  |  \     |
! |  /     |     \  |
! 3        |        3
! |  \     |     /  |
! |     \  |  /     |
! |        3        |
! |     /  |  \     |
! |  /     |     \  |
! 2        |        2
! |  \     |     /  |
! |     \  |  /     |
! |        2        |
! |     /  |  \     |
! |  /     |     \  |
! 1        |        1
! |  \     |     /  |
! |     \  |  /     |
! |        1        |
! |     /  |  \     |
! |  /     |     \  |
! 0        |        0
!    \     |     /  
!       \  |  /     
!          0        
!
subroutine initializeTri()
  integer(kind=4) :: circumferenceEdges, axialEdges
  integer(kind=4) :: totalNumNodes
  integer(kind=4) :: nodesPerFace = 3
  real(kind=8) :: axialEdgeLen, circEdgeAngle, angle, x, y, z,&
    midpointAxial, midpointCirc1, midpointCirc2, offset
  type(PartitionBounds) :: partBounds
  integer(kind=4), allocatable :: elemAxialRange(:), elemCircRange(:), &
    nodeAxialRange(:), nodeCircRange(:)
  integer(kind=4) :: e, c, n, nai, na, nodeIndex, nci, nc, &
    elemIndex, eai, ea, eci, ec, baseNode, nextCirc, nextAxial, &
    nextAxialNextCirc, a, b

  circumferenceEdges = meshProperties%triRefine
  axialEdges = meshProperties%triRefine * 3

  axialEdgeLen = meshProperties%pipeLength / real(axialEdges, kind=8)
  circEdgeAngle = 2.0D0 * pi / real(circumferenceEdges, kind=8)

  ! Total number of nodes in this region (across all partitions).
  totalNumNodes = circumferenceEdges * (axialEdges + 1)

  ! Get the partition bounds.
  partBounds = getPartitionBounds(&
    meshProperties%myrank, meshProperties%numranks, &
    0.0D0, &
	meshProperties%pipeLength, &
    0.0D0, 2.0D0 * pi, meshProperties%triPipeRadius)

  allocate(elemAxialRange(0))
  do e = 1, axialEdges
    midpointAxial = 0.5D0 * axialEdgeLen + (e - 1) * axialEdgeLen
    if (midpointAxial .GE. partBounds%axialStart .AND. &
      midpointAxial .LT. partBounds%axialEnd) then
      call addEnd(elemAxialRange, e - 1)
    endif
  enddo

  allocate(elemCircRange(0))
  do c = 1, circumferenceEdges
    ! There are two tri elements for every circumference edge.
    midpointCirc1 = (c - 1) * circEdgeAngle
    midpointCirc2 = (real(c - 1, kind=8) + 0.5D0) * circEdgeAngle
    if (midpointCirc1 .GE. partBounds%circStartAngle .AND. &
      midpointCirc1 .LT. partBounds%circEndAngle) then
      call addEnd(elemCircRange, (c - 1) * 2)
    endif
    if (midpointCirc2 .GE. partBounds%circStartAngle .AND. &
      midpointCirc2 .LT. partBounds%circEndAngle) then
      call addEnd(elemCircRange, (c - 1) * 2 + 1)
    endif
  enddo

  ! Add overlap layers.
  call addOverlapLayers(circumferenceEdges * 2, axialEdges, elemAxialRange, elemCircRange)

  ! Calculate the number of elements and nodes in the current partition.
  numElemsTri = size(elemCircRange) * size(elemAxialRange)
  numNodesTri = 0

  ! Calculate axial node range.
  ! (same as element axial range plus one row of nodes at the end).
  allocate(nodeAxialRange(0))
  if (numElemsTri .GT. 0) then
    do n = 1, size(elemAxialRange)
      call addEnd(nodeAxialRange, elemAxialRange(n))
    enddo
    call addEnd(nodeAxialRange, elemAxialRange(size(elemAxialRange)) + 1)

    ! Number of nodes along circumference can vary, depending on axial row.
    ! Loop over axial rows, use helper function to get number of nodes.
    do nai = 1, size(nodeAxialRange)
      na = nodeAxialRange(nai)
      nodeCircRange = getNodeCircRange(na, elemCircRange, circumferenceEdges)
      numNodesTri = numNodesTri + size(nodeCircRange)
      deallocate(nodeCircRange)
    enddo
  endif

  ! Allocate and initialize mesh and solution data arrays.
  allocate(nodeCoordsTri(numNodesTri * 3))
  allocate(elemNodeCountsTri(numElemsTri))
  elemNodeCountsTri(:) = nodesPerFace
  allocate(elemNodeIdsTri(numElemsTri * nodesPerFace))
  allocate(linear1Tri(numNodesTri))
  allocate(linear2Tri(numNodesTri))

  ! If there are no elements in this partition, return.
  if (numElemsTri .NE. 0) then
    ! Fill node coordinates.
    nodeIndex = 0
    do nai = 1, size(nodeAxialRange)
      na = nodeAxialRange(nai)
      nodeCircRange = getNodeCircRange(na, elemCircRange, circumferenceEdges)
      do nci = 1, size(nodeCircRange)
        nc = nodeCircRange(nci)
        offset = real(modulo(na, 2), kind=8) * (-0.5D0)
        angle = 2.0D0 * pi * (nc + offset) / circumferenceEdges
        x = cos(angle) * meshProperties%triPipeRadius
        y = sin(angle) * meshProperties%triPipeRadius
        z = na * axialEdgeLen
        nodeCoordsTri(nodeIndex * 3 + 1) = x
        nodeCoordsTri(nodeIndex * 3 + 2) = y
        nodeCoordsTri(nodeIndex * 3 + 3) = z
        nodeIndex = nodeIndex + 1
      enddo
      deallocate(nodeCircRange)
    enddo

    ! Fill element-to-node connectivity.
    elemIndex = 0
    do eai = 1, size(elemAxialRange)
      ea = elemAxialRange(eai)
      do eci = 1, size(elemCircRange)
        ec = elemCircRange(eci)
        ! Base node.
        baseNode = ea * circumferenceEdges + ec / 2
        ! Next node along circumferential direction.
        nextCirc = baseNode + 1
        if (baseNode + 1 .GE. (ea + 1) * circumferenceEdges) then
          nextCirc = ea * circumferenceEdges
        endif
        ! Next node along axial direction, relative to base node.
        nextAxial = baseNode + circumferenceEdges
        ! Next node along circumferential direction, relative to nextCirc node.
        nextAxialNextCirc = nextCirc + circumferenceEdges
        ! Get element node ids. Connectivity is different, depending on whether
        ! axial row and circumferential row are odd or even.
        if (modulo(ea, 2) .EQ. 0) then
          if (modulo(ec, 2) .EQ. 0) then
            a = baseNode
            b = nextAxialNextCirc
            c = nextAxial
          else
            a = baseNode
            b = nextCirc
            c = nextAxialNextCirc
          endif
        else
          if (modulo(ec, 2) .EQ. 0) then
            a = baseNode
            b = nextCirc
            c = nextAxial
          else
            a = nextCirc
            b = nextAxialNextCirc
            c = nextAxial
          endif
        endif
        elemNodeIdsTri(elemIndex * nodesPerFace + 1) = a + nodeIdOffset
        elemNodeIdsTri(elemIndex * nodesPerFace + 2) = b + nodeIdOffset
        elemNodeIdsTri(elemIndex * nodesPerFace + 3) = c + nodeIdOffset
        elemIndex = elemIndex + 1
      enddo
    enddo

    ! Fill solution data for linear2 variable with linear profile.
    do n = 1, numNodesTri
      x = nodeCoordsTri((n - 1) * 3 + 1)
      y = nodeCoordsTri((n - 1) * 3 + 2)
      z = nodeCoordsTri((n - 1) * 3 + 3)
      linear2Tri(n) = 1.0D0 * x + 2.0D0 * y + 3.0D0 * z + 4.0D0
      linear1Tri(n) = 0.0D0
    enddo

  endif

  deallocate(elemAxialRange)
  deallocate(elemCircRange)
  deallocate(nodeAxialRange)

  nodeIdOffset = totalNumNodes

end subroutine initializeTri

subroutine initializeMesh
  call initializeQuad()
  call initializeTri()
end subroutine initializeMesh

end module PipeMeshGenerator
