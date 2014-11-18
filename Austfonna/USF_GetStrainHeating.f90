!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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
!
!/******************************************************************************
! *
! *  Module containing a solver for computing the strain heating and function for friction heat
! *
! ******************************************************************************
! *
! *  Authors: Thomas Zwinger, Juha Ruokolainen, Hakime Seddik
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *  Current date: 20 July 2012 (Martina)
! *
! *****************************************************************************/

! contains subroutine for exporting strain heating (if strain heat is taken into account as
! heat source is defined by "Friction Heat = Logical True" in the
! Body Force Section and independent from running this Solver.
! ATTENTION: consistent units for capacity and conductivy are needed.
! usage:
! Solver 1
!   Equation = "StrainHeating"
!   Variable = String "StrainHeat"
!   Variable DOFs = 1
!   Procedure = File "USF_GetStrainHeating" "getStrainHeating"
!   Nonlinear System Max Iterations = 1
! End


! contains function for adding friction heat to geothermal heat flux at the basis and exporting it
! ATTENTION: consistent units for capacity and conductivy are needed.
! ATTENTION: The keyword "Friction Heat" does make reference to Strainheating, not Friction heat!
! usage:
! Solver 1
!   Exec Solver = "Before Simulation"
!   Equation = "Normal vector"
!   Variable = "Normal Vector"   
!   ! in 3dimensional simulations we have 3 entries
!   Variable DOFs = 3 
!   !NB: does not need to actually solve a matrix
!   !    hence no BW optimization needed
!   Optimize Bandwidth = Logical False 
!   Procedure = "ComputeNormal" "ComputeNormalSolver"
!   ! if set to True, all boundary normals would be computed by default
!   ComputeAll = Logical False
!End
!Solver 2
!  Equation = String "StressSolver"
!  Procedure =  File "ComputeDevStressNS" "ComputeDevStress"
!  ! this is just a dummy, hence no output is needed
!  !-----------------------------------------------------------------------
!  Variable = -nooutput "Sij"
!  Variable DOFs = 1
!  ! the name of the variable containing the flow solution (U,V,W,Pressure)
!  !-----------------------------------------------------------------------
!  Flow Solver Name = String "Flow Solution"
!  Exported Variable 1 = "Stress" ! [Sxx, Syy, Szz, Sxy] in 2D
!                                 ! [Sxx, Syy, Szz, Sxy, Syz, Szx] in 3D
!  Exported Variable 1 DOFs = 6   ! 4 in 2D, 6 in 3D
!  Linear System Solver = "Iterative"
!  Linear System Iterative Method = "BiCGStab"
!  Linear System Max Iterations = 300
!  Linear System Convergence Tolerance = 1.0E-09
!  Linear System Abort Not Converged = True
!  Linear System Preconditioning = "ILU0"
!  Linear System Residual Output = 1
!End
!Solver 3
!  Equation = "Dummy"
!  Procedure = File "DummySolver" "DummySolver"
!  Variable = String "friction"
!  Variable DOFs = 1
!End
!Material 1
! ...
!   Cauchy = Logical True !needed for ComputeDevStressNS
!
!! bedrock:
!Boundary Condition 1
!  Name = "bed"
!  ComputeNormal = Logical True !needed for ComputeNormal
! ....
!  Temp Flux BC = Logical True
!  friction = Variable Coordinate 1
!       Real Procedure "./USF_GetStrainHeating" "getFrictionHeat"
!  Temp Heat Flux = Variable friction
!      Real MATC "0.063*(31556926.0)*1.0E-06 +tx" !assuming 63mW/m2 as geothermal heatflux
!End


!------------------------------------------------------------------------------
!> Subroutine for solving exporting local strain heating
!> \ingroup Solvers

SUBROUTINE getStrainHeating( Model,Solver,dt,TransientSimulation )

  USE Types
  USE Lists
  USE Integration
  USE ElementDescription
  USE MaterialModels
  USE SolverUtils
  USE DefUtils

  IMPLICIT NONE

  !------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !----------------------------------------------------------------------------

  INTEGER :: t,i,j,k,n,m,NSDOFs,istat, DIM,jj
  REAL(KIND=dp),ALLOCATABLE :: Viscosity(:), Density(:), &
       UX(:), UY(:), UZ(:), Pressure(:),dBasisdx(:,:),Basis(:)
  REAL(KIND=dp), POINTER :: Var(:),FlowSolution(:),Variable(:)
  REAL(KIND=dp) :: U,V,W,detJ,LocalViscosity, LocalDensity,&
       LocalU(3), GradU(3,3)
  INTEGER, ALLOCATABLE :: NumberOfVisits(:)
  INTEGER, POINTER :: FlowPerm(:),VarPerm(:)
  LOGICAL :: AllocationsDone=.FALSE., FirstTime=.TRUE., Found, lstat
  TYPE(Variable_t), POINTER :: FlowSol,VarSol
  TYPE(Element_t),POINTER :: Element
  TYPE(Nodes_t)   :: ElementNodes
  TYPE(ValueList_t), POINTER :: Material, SolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

  SAVE AllocationsDone, FirstTime, SolverName, Viscosity, Density,&
       UX, UY, UZ, Pressure, DIM,NumberOfVisits, M, Basis, dBasisdx


  VarSol     => Solver % Variable
  IF ( ASSOCIATED( VarSol ) ) THEN
     VarPerm    => VarSol % Perm
     Variable => VarSol % Values
     Variable = 0.0_dp
  ELSE
     CALL FATAL(SolverName,'Solver variable not associated')
  END IF

  IF (FirstTime) THEN
     WRITE(SolverName, '(A)') 'getStrainHeating'
     DIM = CoordinateSystemDimension()
     FirstTime = .FALSE.
  END IF

  IF ( (.NOT.AllocationsDone) .OR. (Solver % Mesh % Changed) ) THEN
     N = Solver % Mesh % MaxElementNodes
     M = Solver % Mesh % NumberOfNodes
     IF(AllocationsDone) THEN
        DEALLOCATE(Viscosity,Density,UX,UY,UZ,Pressure,&
             Basis,dBasisdx,NumberOfVisits,STAT=istat)
        IF ( istat /= 0 ) THEN
           CALL Fatal( SolverName, 'Memory deallocation error' )
        END IF
     END IF
     ALLOCATE(Viscosity( N ),&
          Density( N ),&
          UX( N ),&
          UY( N ),&
          UZ( N ),&
          Pressure ( N ),&
          Basis( 2*N ),&
          dBasisdx( 2*N, DIM ),&
          NumberOfVisits(M),&
          STAT=istat)
     IF ( istat /= 0 ) THEN
        CALL Fatal( SolverName, 'Memory allocation error' )
     END IF
     AllocationsDone = .TRUE.
  END IF
  
  FlowSol => VariableGet( Solver % Mesh % Variables, 'Flow Solution' )
  IF ( ASSOCIATED( FlowSol ) ) THEN
     FlowPerm     => FlowSol % Perm
     NSDOFs       =  FlowSol % DOFs
     FlowSolution => FlowSol % Values
  ELSE
     CALL Fatal( SolverName, 'Variable <Flow Solution> not associated')
  END IF
  
  NumberOfVisits=0
  DO t=1,Solver % NumberOfActiveElements !loop over elements
     Element => GetActiveElement(t)
     IF (.NOT.ASSOCIATED(Element)) &
          CALL FATAL(SolverName,'Element not found')
     N = GetElementNOFNodes(Element)
     CALL GetElementNodes( ElementNodes, Element )
     Material => GetMaterial(Element)
     IF (.NOT.ASSOCIATED(Material)) &
          CALL FATAL(SolverName,'Material not found')
     Viscosity(1:N) = GetReal( Material,'Viscosity', Found )
     IF(.NOT.FOUND) &
          CALL FATAL(SolverName,'Viscosity not found')
     Density(1:N) = GetReal( Material,'Density', Found )
     IF(.NOT.FOUND) &
          CALL FATAL(SolverName,'Density not found')

     DO i=1,N ! loop over nodes in element to set velocities
        k = FlowPerm(Element % NodeIndexes(i))
        IF ( k > 0 ) THEN
           SELECT CASE( NSDOFs )
           CASE(3)
              UX(i) = FlowSolution( NSDOFs*k-2 )
              UY(i) = FlowSolution( NSDOFs*k-1 )
              UZ(i) = 0.0D0
           CASE(4)
              UX(i) = FlowSolution( NSDOFs*k-3 )
              UY(i) = FlowSolution( NSDOFs*k-2 )
              UZ(i) = FlowSolution( NSDOFs*k-1 )
           END SELECT
        ELSE
           UX(i) = 0.0d0
           UY(i) = 0.0d0
           UZ(i) = 0.0d0
        END IF
     END DO !loop over nodes in elment

     DO i=1,N     ! loop nodes in element
        jj=Element % NodeIndexes(i)
        NumberOfVisits(jj) = NumberOfVisits(jj) + 1
        ! get local coordinates of the point i inside the element
        U = Element % Type % NodeU(i)
        V = Element % Type % NodeV(i)
        W = Element % Type % NodeW(i)
        ! get local information on test-functions and derivatives of the point i
        lstat = ElementInfo( Element,ElementNodes,U,V,W,detJ, &
             Basis,dBasisdx) 
        LocalViscosity = SUM( Viscosity(1:n) * Basis(1:n) )
        LocalDensity = SUM( Density(1:n) * Basis(1:n) )
        LocalViscosity = EffectiveViscosity( LocalViscosity, LocalDensity,&
             Ux, Uy, Uz, Element, ElementNodes, N, N, u,v,w )
        DO j=1,DIM
           GradU(1,j) = SUM( UX(1:n)*dBasisdx(1:n,j) )
           GradU(2,j) = SUM( UY(1:n)*dBasisdx(1:n,j) )
           IF (DIM == 3) GradU(3,j) = SUM( UZ(1:n)*dBasisdx(1:n,j) )
        END DO
        LocalU = 0.0_dp
        LocalU(1) = SUM( UX(1:n)*Basis(1:n) )
        LocalU(2) = SUM( UY(1:n)*Basis(1:n) )
        IF (DIM==3) LocalU(3) = SUM( UZ(1:n)*Basis(1:n) )
        Variable(VarPerm(Element % NodeIndexes(i))) = &
             Variable(VarPerm(Element % NodeIndexes(i)))&
             + 0.5d0 * LocalViscosity *SecondInvariant(LocalU,GradU)
     END DO ! loop nodes in element
  END DO ! loop elements
  
  
  DO i=1,M !loop over all nodes
     IF (NumberOfVisits(i) == 0) CALL FATAL(SolverName,'Division by zero')
     Variable(VarPerm(i)) = Variable(VarPerm(i))/NumberOfVisits(i)
  END DO ! loop over all nodes
  
  
END SUBROUTINE getStrainHeating


!****************************


!-----------
!function to export and define FrictionHeat at bedrock

FUNCTION getFrictionHeat(  Model, Node )RESULT(frictionheat)
  
  USE Types
  USE Lists
  USE Integration
  USE ElementDescription
  USE MaterialModels
  USE SolverUtils
  USE DefUtils

  IMPLICIT NONE
  
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: DummyInput, frictionHeat
  !----------------------------------------------------------------------------
  
  INTEGER :: DIM, i, j, Ind(3,3)
  REAL(KIND=dp), POINTER :: FlowValues(:),NormalValues(:),StressValues(:)
  REAL(KIND=dp) :: normal(3), velo(3), un, ut(3), Sig(3,3), Sn(3)
  INTEGER, POINTER :: FlowPerm(:),StressPerm(:), NormalPerm(:)
  LOGICAL :: FirstTime=.TRUE.
  TYPE(Variable_t), POINTER :: FlowSol,StressVariable, NormalVar
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
  
  SAVE FirstTime, DIM, SolverName
  
  IF (FirstTime) THEN
     WRITE(SolverName, '(A)') 'getFrictionHeat'
     DIM = CoordinateSystemDimension()
     FirstTime = .FALSE.
     DO i=1, 3
        Ind(i,i) = i
     END DO
     Ind(1,2) = 4
     Ind(2,1) = 4
     Ind(2,3) = 5
     Ind(3,2) = 5
     Ind(3,1) = 6
     Ind(1,3) = 6
  END IF
  
  ! Get the variable velocity
  !---------------------------
  FlowSol => VariableGet( Model % Variables, 'Flow Solution' )
  IF ( ASSOCIATED( FlowSol ) ) THEN
     FlowPerm    => FlowSol % Perm
     FlowValues  => FlowSol % Values
  ELSE
     CALL FATAL(SolverName, 'Need NS Solver, Flow Solution not found')
  END IF
  
  ! Get the stress variable
  !------------------------
  StressVariable => VariableGet( Model % Variables, 'Stress' )
  IF ( ASSOCIATED( StressVariable ) ) THEN
     StressPerm    => StressVariable % Perm
     StressValues  => StressVariable % Values
  ELSE
     CALL FATAL(SolverName,'No variable Stress found')   
  END IF
  
  ! Get the variable for normal vector
  NormalVar =>  VariableGet(Model % Variables,'Normal Vector')
  IF ( ASSOCIATED( NormalVar ) ) THEN
     NormalPerm => NormalVar % Perm
     NormalValues => NormalVar % Values
  ELSE
     CALL FATAL(SolverName, 'Normal Vector variable not found')
  END IF
  
  DO i=1, DIM
     normal(i) = NormalValues(DIM*(NormalPerm(Node)-1) + i)      
     velo(i) = FlowValues( (DIM+1)*(FlowPerm(Node)-1) + i )
  END DO
  
  !Tangential velocity
  un = SUM(velo(1:DIM)*(normal(1:DIM))) 
  ut = velo(1:DIM)-un*(normal(1:DIM))
  
  DO i=1, DIM
     DO j= 1, DIM
        Sig(i,j) =  &
             StressValues( 2*DIM *(StressPerm(Node)-1) + Ind(i,j) )
     END DO
  END DO
  
  DO i=1, DIM
     Sn(i) = SUM(Sig(i,1:DIM)*normal(1:DIM)) 
  END DO
  frictionHeat =(SUM(abs(ut(1:DIM)*Sn(1:DIM))))
END FUNCTION getFrictionHeat
