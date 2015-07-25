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
! *
! *  Original Date: 24 July 2015 
! * 
! *****************************************************************************
SUBROUTINE WeightVariable( Model,Solver,dt,TransientSimulation )
  !DEC$ATTRIBUTES DLLEXPORT :: SIAVariable
  !------------------------------------------------------------------------------
  !******************************************************************************
  !
  ! Allow to have the weights as primary variables (and not exported one)
  ! Allow then to have access to the Previous values to construct a stable 
  ! time discretization sheme. 
  !
  !******************************************************************************
  USE DefUtils
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
END SUBROUTINE WeightVariable
!------------------------------------------------------------------------------


! *****************************************************************************
SUBROUTINE WeightFinder( Model,Solver,dt,TransientSimulation )
  !DEC$ATTRIBUTES DLLEXPORT :: SIASolver
  !------------------------------------------------------------------------------
  !******************************************************************************
  !
  !  Solve the SIA Flow solution !
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

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------

  TYPE(ValueList_t), POINTER :: SolverParams, BodyForce, Material
  TYPE(Variable_t), POINTER :: WeightSol, AppSol, PointerToVariable

  LOGICAL :: AllocationsDone = .FALSE., Found, GotIt

  INTEGER :: i, n, m, t, k, l, istat, DIM, p, nsize, bc, bottomindex
  INTEGER, POINTER :: WeightPerm(:), AppPerm(:),Permutation(:) 
  INTEGER, SAVE :: Timestep

  TYPE(Variable_t), POINTER :: TimeVar

  REAL(KIND=dp) :: PhiWeights

  REAL(KIND=dp), POINTER :: Weight(:), App(:), VariableValues(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

  SAVE  AllocationsDone, DIM, SolverName

  TYPE(Element_t),POINTER :: Element


  !-----------Variables needed for integration --------------------------------
  TYPE(Solver_t), POINTER :: PSolver
  INTEGER, POINTER :: TopPointer(:),BotPointer(:),UpPointer(:),DownPointer(:)
  LOGICAL :: Initialized = .FALSE.,GotVar
  REAL(KIND=dp), POINTER :: Coord(:)
  TYPE(Variable_t), POINTER :: Var, OldVar
  SAVE BotPointer, TopPointer, UpPointer, DownPointer

  SolverName = "Weight_Solver"


  !------------------------------------------------------------------------------
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  N = Solver % Mesh % MaxElementDOFs !needed right now only for beta
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
  WRITE(*,*) 'Extruded detected'
  !------------------------------------------------------------------------------
  !    Get variables needed for solution
  !------------------------------------------------------------------------------
  DIM = CoordinateSystemDimension()
  WeightSol => VariableGet( Solver % Mesh % Variables, 'NodalWeights' )
  IF (ASSOCIATED(WeightSol)) THEN
     Weight => WeightSol % Values
     WeightPerm => WeightSol % Perm
  ELSE
     CALL FATAL(SolverName,'Could not find variable >NodalWeights<')
  END IF
  WRITE(*,*) 'Found Nodal Weights Variable'
  AppSol => VariableGet( Solver % Mesh % Variables, 'ApproximationLevel' )
  IF (ASSOCIATED(AppSol)) THEN
     App => AppSol % Values
     AppPerm => AppSol % Perm
  ELSE
     CALL FATAL(SolverName,'Could not find variable >ApproximationLevel<')
  END IF
  WRITE(*,*) 'Found Approximation Level'

  PhiWeights=GetConstReal( Solver % Values, 'Phi', Found)
  !------------------------------------------------------------------------------
  !   Giving FS nodes and SIA nodes different weights. Phi can be anything 
  !   between 0 and 1. phi = 0 gives FS nodes twice as important as SIA nodes. 
  !   phi = 1 gives SIA nodes no importance. 
  !------------------------------------------------------------------------------

  App=(App-1.0)*10.0+PhiWeights
 
  !------------------------------------------------------------------------------
  !   Summation of weight over columns
  !------------------------------------------------------------------------------
  !summation from bottom to top
  DO i=1,Model % Mesh % NumberOfNodes
     Weight(WeightPerm(i)) = 0.0_dp !Integral is zero at the bottom
     IF( i == BotPointer(i) ) THEN !so this loop will be over the extruded lines
        l = i !start at the bottom of the line
        DO k=1,nsize
        Weight(WeightPerm(i))=Weight(WeightPerm(i))+App(AppPerm(l))
        IF( l == TopPointer(l)) EXIT !if reached end of extruded line, go to next one
        l = UpPointer(l)   !step it up  
    	END DO  
     END IF
  END DO

  !------------------------------------------------------------------------------
  !   WRite to file
  !------------------------------------------------------------------------------
 
  open (unit=77, file="mesh.weights")
 
  DO i = 1, Model % Mesh % NumberOfNodes
  IF( i == BotPointer(i) ) THEN
    WRITE(77,*)  i, Weight(WeightPerm(i))
  END IF
  END DO
  close(77)

  !------------------------------------------------------------------------------
END SUBROUTINE WeightFinder
!------------------------------------------------------------------------------

