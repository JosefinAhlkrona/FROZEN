!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - Scientific Computing Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or 
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!/******************************************************************************
! *
! *  SIASolver: Solver to inquire the velocity and isotropic pressure from 
! *                   the SIA solution
! *  Exported Variables SIAFlow, dof=dim
! *
! *
! ******************************************************************************
! *
! *  Authors: Josefin Ahlkrona
! *  Email:   josefin.ahlkrona@it.uu.se
! *  Web:     http://www.csc.fi/elmer
! *
! *  Original Date: 25 June 2013 
! * 
! *****************************************************************************

! *****************************************************************************
SUBROUTINE NeighbourFinder( Model,Solver,dt,TransientSimulation )
  !DEC$ATTRIBUTES DLLEXPORT :: SIASolver
  !------------------------------------------------------------------------------
  !******************************************************************************
  !
  !  Find neighbours in horizontal direction and write info to file neighbours.txt
  !
  !  ARGUMENTS:
  !
  !  TYPE(Model_t) :: Model,  
  !     INPUT: All model information (mesh, materials, BCs, etc...)
  !
  !  TYPE(Solver_t) :: Solver
  !     INPUT: Linear & nonlinear equation solver options
  !
  !  REAL(KIND=dp) :: dt,
  !     INPUT: Timestep size for time dependent simulations
  !
  !  LOGICAL :: TransientSimulation
  !     INPUT: Steady state or transient simulation
  !
  !******************************************************************************
  USE DefUtils
  USE ElementUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------

  

  INTEGER :: i, n, m, t, istat, DIM, p, COMP, bc, bottomindex, layercounter,nodeatlayer
  

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

  SAVE  DIM, SolverName

  TYPE(Nodes_t) :: ElementNodes
  TYPE(Element_t),POINTER :: Element

  !-----------Variables needed for integration --------------------------------
  TYPE(Solver_t), POINTER :: PSolver
  INTEGER :: j,k,l,Dofs,dof,nsize,TopNodes,BotNodes
  INTEGER, POINTER :: TopPointer(:),BotPointer(:),UpPointer(:),DownPointer(:)
  LOGICAL :: Initialized = .FALSE.,GotVar
  REAL(KIND=dp), POINTER :: Coord(:)
  TYPE(Variable_t), POINTER :: Var, OldVar

  INTEGER, ALLOCATABLE :: NeighbourElementList(:,:),NeighbourNodeList(:,:)
  LOGICAL :: AllocationsDone = .FALSE.

  INTEGER :: nofelemsconnected, nofnodesconnected
  !-----Needed for transient stuff etc.

  SolverName = "Neighbour_Finder"


   IF ( .NOT.AllocationsDone .OR. Solver % Mesh % Changed ) THEN

     N = Solver % Mesh % MaxElementDOFs

     IF( AllocationsDone ) THEN
        DEALLOCATE( NeighbourElementList, NeighbourNodeList,STAT=istat )
     END IF

     ALLOCATE(NeighbourElementList( Model % Mesh % NumberOfNodes,30 ), &
 NeighbourNodeList( Model % Mesh % NumberOfNodes,20 ), &
       STAT=istat)

     IF ( istat /= 0 ) THEN
        CALL Fatal( SolverName,'Memory allocation error, Aborting.' )
     END IF

     !------------------------------------------------------------------------------

     AllocationsDone = .TRUE.
  END IF
  !------------------------------------------------------------------------------


  NeighbourNodeList=0
  NeighbourElementList=0


  WRITE( Message, * ) 'Detecting structure'
    CALL Info( SolverName, Message, Level=4 )
  !------------------------------------------------------------------------------
  !   Initialize the pointers to top and bottom nodes, needed for the integration 
  !------------------------------------------------------------------------------
  IF( .NOT. Initialized ) THEN

     ! Choose active direction coordinate and set corresponding unit vector
     !---------------------------------------------------------------------
     PSolver => Solver
     CALL DetectExtrudedStructure( Solver % Mesh, PSolver, Var, &
          TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
          UpNodePointer = UpPointer, DownNodePointer = DownPointer )

     Coord => Var % Values
     nsize = SIZE( Coord )
     Initialized = .TRUE.

  END IF

 !------------------------------------------------------------------------------
 !    Find connected elements
 !------------------------------------------------------------------------------

  DO i = 1, GetNOFActive()
     Element => GetActiveElement(i)
     n = GetElementNOFNodes()
     DO j=1,n
        k = Element % NodeIndexes(j)
        !find layer index
        bottomindex=BotPointer(k)
	nodeatlayer=bottomindex
	    layercounter=1
            DO WHILE(nodeatlayer .NE. k)  !step up if we're not there yet
               nodeatlayer = UpPointer(nodeatlayer) 
		layercounter=layercounter+1
	    END DO
        NeighbourElementList(k,1)=NeighbourElementList(k,1)+1 !one more element is connected
        NeighbourElementList(k,2)=layercounter !stores what layer this node is at

        nofelemsconnected=NeighbourElementList(k,1)
        NeighbourElementList(k,nofelemsconnected+2)= i 
     END DO
 END DO

 !------------------------------------------------------------------------------
 !    Based on connected Elements, find connected nodes
 !------------------------------------------------------------------------------

  DO i = 1, Model % Mesh % NumberOfNodes !find neighbours for all nodes
  nofelemsconnected=  NeighbourElementList(i,1)
   NeighbourNodeList(i,1)=0
    DO j=1,nofelemsconnected  !check all elements that node is connected to
      Element = GetActiveElement(NeighbourElementList(i,2+j))
      k = GetElementNOFNodes()
      nofnodesconnected=NeighbourNodeList(i,1)
      DO l=1,k   !check all nodes that is in a connecting element 
         m = Element % NodeIndexes(l)
         IF (.NOT. ANY(m .EQ. NeighbourNodeList(i,2:nofnodesconnected+1)) &
	.AND. (m .NE. i) .AND. NeighbourElementList(i,2).EQ.NeighbourElementList(m,2)) THEN !if this node was already loaded, it is not itself and it belongs to this layer
         !  WRITE(*,*) 'neighbour node number', m, ' is not yet written down'
          NeighbourNodeList(i,1)=NeighbourNodeList(i,1)+1 !one more node is connected
          nofnodesconnected=NeighbourNodeList(i,1)
          NeighbourNodeList(i,nofnodesconnected+1)=m
         END IF
      END DO
    END DO
 END DO



  !------------------------------------------------------------------------------
  !    Print list to file
  !------------------------------------------------------------------------------
  DIM = CoordinateSystemDimension()
 
  open (unit=412, file="neighbours.txt")
  WRITE(412,*) MAXVAL(NeighbourNodeList(:,1),1)
  DO i = 1, Model % Mesh % NumberOfNodes 
     WRITE(412,"(I7.2)",advance='no'), i
     DO j=1,NeighbourNodeList(i,1)+1
	WRITE(412,"(I7.2)",advance='no') NeighbourNodeList(i,j)
     END DO
     WRITE(412,*) ' '
  END DO

  close(412)

  !------------------------------------------------------------------------------
END SUBROUTINE NeighbourFinder
!------------------------------------------------------------------------------

