!----------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE InitVeloTemp( Model,Solver,dt,TransientSimulation )
!---------------------------------------------------------------------------------- 

  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve stress equations for one timestep
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations (NOTE: Not used
!            currently)
!
!******************************************************************************

   TYPE(Model_t)  :: Model
   TYPE(Solver_t), TARGET :: Solver

   LOGICAL ::  TransientSimulation
   REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
 
 INTEGER :: nvalue, i, j, k, l, m, n, o, FlowDOFs, SicoVxDOFS, SicoVyDOFS, SicoVzDOFS, DepthDOFs, &
         TempDOFs, HomoTempDOFs, SicoTempDOFs, DIM, istat
 TYPE(Element_t), POINTER :: CurrentElement
 TYPE(ValueList_t), POINTER :: Material
 REAL(KIND=dp), POINTER :: FlowValues(:), SicoVxVal(:),  SicoVyVal(:),  SicoVzVal(:), DepthValues(:), &
                        TempVal(:), HomoTempVal(:), SicoTempVal(:)
 REAL(KIND=dp), ALLOCATABLE :: UpperLimit(:)
 REAL(KIND=dp) :: temperature
 TYPE(Variable_t), POINTER :: FlowSol, SicoVxSol, SicoVySol, SicoVzSol, DepthSol, TempSol, HomoTempSol, &
                       SicoTempSol
 INTEGER, POINTER :: FlowPerm(:), SicoVxPerm(:), SicoVyPerm(:), SicoVzPerm(:), DepthPerm(:), TempPerm(:), &
                 HomoTempPerm(:), SicoTempPerm(:), NodeIndexes(:)
 REAL(KIND=dp) :: rho, g
 CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, TempName
 LOGICAL :: AllocationsDone = .FALSE., Found

 SAVE AllocationsDone, UpperLimit, DIM

 SolverName = 'InitVeloTemp'

 IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
     M = Model % MaxElementNodes
     dim = CoordinateSystemDimension()  

   IF ( AllocationsDone ) THEN
      DEALLOCATE( UpperLimit )
   END IF

   ALLOCATE( UpperLimit(M), STAT=istat )

   IF ( istat /= 0 ) THEN
        CALL FATAL( SolverName, 'Memory allocation error' )
   ELSE
           CALL INFO(SolverName, 'Memory allocation done', level=1 )
   END IF
        
   AllocationsDone = .TRUE.

 END IF

 rho = 910.0_dp*1.0E-06_dp*(31556926.0_dp)**(-2.0_dp)
 g = 9.81_dp*(31556926.0_dp)**(2.0_dp)

 FlowSol => VariableGet( Solver % Mesh % Variables, 'flow solution' )
 IF ( ASSOCIATED( FlowSol ) ) THEN
       FlowPerm    => FlowSol % Perm
       FlowValues  => FlowSol % Values
       FlowDOFs = FlowSol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find velocity field variable')
 END IF

 SicoVxSol => VariableGet( Solver % Mesh % Variables, 'sicovx' )
 IF ( ASSOCIATED( SicoVxSol ) ) THEN
       SicoVxPerm    => SicoVxSol % Perm
       SicoVxVal  => SicoVxSol % Values
       SicoVxDOFs = SicoVxSol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find variable for Sicopolis vx')
 END IF

 SicoVySol => VariableGet( Solver % Mesh % Variables, 'sicovy' )
 IF ( ASSOCIATED( SicoVySol ) ) THEN
       SicoVyPerm    => SicoVySol % Perm
       SicoVyVal  => SicoVySol % Values
       SicoVyDOFs = SicoVySol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find variable for Sicopolis vy')
 END IF

 SicoVzSol => VariableGet( Solver % Mesh % Variables, 'sicovz' )
 IF ( ASSOCIATED( SicoVzSol ) ) THEN
       SicoVzPerm    => SicoVzSol % Perm
       SicoVzVal  => SicoVzSol % Values
       SicoVzDOFs = SicoVzSol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find variable for Sicopolis vz')
 END IF
 
 DepthSol => VariableGet( Solver % Mesh % Variables, 'Depth' )
 IF ( ASSOCIATED( DepthSol ) ) THEN
       DepthPerm => DepthSol % Perm
       DepthValues => DepthSol % Values
       DepthDOFs = DepthSol % DOFs
 ELSE
       CALL FATAL(SolverName,'Could not find Depth field variable')
 END IF

 TempSol => VariableGet( Solver % Mesh % Variables, 'Temp' )
 IF ( ASSOCIATED( TempSol ) ) THEN
       TempPerm => TempSol % Perm
       TempVal => TempSol % Values
       TempDOFs = TempSol % DOFs
 ELSE
       CALL FATAL(SolverName,'Could not find Temperature field variable')
 END IF

 HomoTempSol => VariableGet( Solver % Mesh % Variables, 'Temp Homologous' )
 IF ( ASSOCIATED( HomoTempSol ) ) THEN
       HomoTempPerm => HomoTempSol % Perm
       HomoTempVal => HomoTempSol % Values
       HomoTempDOFs = HomoTempSol % DOFs
 ELSE
       CALL FATAL(SolverName,'Could not find Temperature Homologous field variable')
 END IF

 SicoTempSol => VariableGet( Solver % Mesh % Variables, 'sicotemp' )
 IF ( ASSOCIATED( SicoTempSol ) ) THEN
       SicoTempPerm => SicoTempSol % Perm
       SicoTempVal => SicoTempSol % Values
       SicoTempDOFs = SicoTempSol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find variable for Sicopolis temperature')
 END IF


 DO i=1,Solver % NumberOFActiveElements

   CurrentElement => GetActiveElement(i)
   NodeIndexes => CurrentElement % NodeIndexes

    DO k=1, GetElementNOFNodes(CurrentElement)
   
     j = FlowPerm(NodeIndexes(k))
     l = SicoVxPerm(NodeIndexes(k))
     m = SicoVyPerm(NodeIndexes(k))
     n = SicoVzPerm(NodeIndexes(k))
     o = DepthPerm(NodeIndexes(k))

     FlowValues(FlowDOFs*(j-1)+1) = SicoVxVal(SicoVxDOFs*(l-1)+1)
     FlowValues(FlowDOFs*(j-1)+2) = SicoVyVal(SicoVyDOFs*(m-1)+1)
     FlowValues(FlowDOFs*(j-1)+3) = SicoVzVal(SicoVzDOFs*(n-1)+1)
     FlowValues(FlowDOFs*(j-1)+4) = rho*g*DepthValues(DepthDOFs*(o-1)+1)

    END DO
 END DO

 DO i=1,Solver % NumberOFActiveElements

   CurrentElement => GetActiveElement(i)
   n = GetElementNOFNodes(CurrentElement)
   NodeIndexes => CurrentElement % NodeIndexes

   Material => GetMaterial()

   TempName =  GetString(Material  ,'Temperature Name', Found)
   IF (.NOT. Found) CALL FATAL(SolverName,'No Temperature Name found')

   UpperLimit(1:n) = ListGetReal( Material, TRIM(TempName) // ' Upper Limit',&
                        n, NodeIndexes, Found)
   IF (.NOT. Found) CALL FATAL(SolverName,'No Upper Limit found')

    DO k=1, GetElementNOFNodes(CurrentElement)

     j = HomoTempPerm(NodeIndexes(k))
     l = sicoTempPerm(NodeIndexes(k))
     m = TempPerm(NodeIndexes(k))

     temperature = SicoTempVal(SicoTempDOFs*(l-1)+1)+273.15_dp
     TempVal(TempDOFs*(m-1)+1) = temperature

     HomoTempVal(HomoTempDOFs*(j-1)+1) = temperature-UpperLimit(k)
    END DO

 END DO


WRITE(Message,'(a)') '---------------------------------------------------------------'
CALL Info(SolverName,Message, Level=3)
CALL Info(SolverName,'Initialize values for Velocities and Temperature:..........done', Level=3)
WRITE(Message,'(a)') '---------------------------------------------------------------'
CALL Info(SolverName,Message, Level=3)

!------------------------------------------------------------------------------
END SUBROUTINE InitVeloTemp
!------------------------------------------------------------------------------
