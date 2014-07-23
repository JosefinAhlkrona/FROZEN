MODULE ErrorEstimationSubs

CONTAINS

  SUBROUTINE SolutionErrorEstimate( Model,Solver,dt,TransientSimulation, &
       NodeType2, SIAVelPermuted, NumberOfSIANodes, NumberOfFSNodes)
    !******************************************************************************
    !
    !  Estimate Approximation Error
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
    !     INPUT: Timestep size for time dependent simulations
    !
    !******************************************************************************
    !------------------------------------------------------------------------------
    USE DefUtils

    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Solver_t), TARGET :: Solver

    REAL(KIND=dp) :: dt
    LOGICAL :: TransientSimulation

    INTEGER, ALLOCATABLE, intent(inout) :: NodeType2(:)
    REAL(KIND=dp),ALLOCATABLE, intent(in) :: SIAVelPermuted(:)
    INTEGER, intent(out) :: NumberOfSIANodes, NumberOfFSNodes    

    !------------------------------------------------------------------------------
    !    Local variables
    !------------------------------------------------------------------------------
    TYPE(Element_t),POINTER :: Element
    INTEGER :: ACounter, i, j, k, n, NSDOFs, istat 
    LOGICAL :: OnlyHorizontalError
    LOGICAL :: AllocationsDone = .FALSE.
    LOGICAL :: gotIt

    REAL(KIND=dp),ALLOCATABLE :: CoupledSolution(:), ErrorInSIA(:), &
         NodeWiseError(:)
    REAL(KIND=dp) :: UNorm

    TYPE(Variable_t), POINTER :: FlowSol, NodeType2Variable
    TYPE(Variable_t), POINTER ::  NodeWiseErrorVariable
    INTEGER, POINTER :: FlowPerm(:), NodeType2Perm(:), NodeWiseErrorPerm(:)
    REAL(KIND=dp), POINTER :: FlowSolution(:), ForceVector(:), &
         NodeType2Values(:), NodeWiseErrorValues(:)
    INTEGER, SAVE :: Timestep 
    TYPE(Variable_t), POINTER :: TimeVar

    TYPE(Matrix_t),POINTER :: A

    REAL(KIND=dp), SAVE :: ErrorBound,  LowerHorVelLimit, ErrorBoundAbs, ErrorBoundRel  

    SAVE OnlyHorizontalError, CoupledSolution, &
         ErrorInSIA, NodeWiseError

    !-----------------------------------------------------------------

    WRITE( Message, * ) 'Computing the Error'
    CALL Info( 'FlowSolve', Message, Level=4 )


    FlowSol => Solver % Variable
    NSDOFs         =  FlowSol % DOFs
    FlowPerm       => FlowSol % Perm
    FlowSolution   => FlowSol % Values

    IF ( .NOT.AllocationsDone .OR. Solver % Mesh % Changed ) THEN
       IF( AllocationsDone ) THEN
          DEALLOCATE(                               &
               ErrorInSIA, &
               NodeWiseError, &
               STAT=istat )
       END IF
       ALLOCATE( &
            ErrorInSIA(  SIZE( FlowSolution )), &   
            NodeWiseError(  Model % Mesh % NumberOfNodes ), &         
            STAT=istat)
       AllocationsDone = .TRUE.
    END IF

    ErrorBoundRel = 0.01*GetConstReal(  Solver % Values,  &
         'Relative Error Allowed In Percent', gotIt )

    IF (.NOT. gotIt) THEN
       WRITE( Message, * ) 'Relative Error Tolerance not found, setting to 15 %'
       CALL Info( 'FlowSolve', Message, Level=4 )
       ErrorBoundRel = 0.15
    END IF

    ErrorBoundAbs = GetConstReal(  Solver % Values,  &
         'Absolute Error Allowed', gotIt )
    IF (.NOT. gotIt) THEN
       WRITE( Message, * ) 'Absolute Error Tolerance not found, setting to 1 m/a'
       CALL Info( 'FlowSolve', Message, Level=4 )
       ErrorBoundAbs = 1.0
    END IF

    LowerHorVelLimit  = GetConstReal(  Solver % Values,  &
         'Horizontal Velocity Regarded as Zero', gotIt )
    IF (.NOT. gotIt) THEN
       WRITE( Message, * ) 'Min hor. vel. in error not found, setting to 0.1 m/a'
       CALL Info( 'FlowSolve', Message, Level=4 )
       LowerHorVelLimit = 0.1
    END IF

    OnlyHorizontalError = GetLogical( Solver % Values, &
         'Only Compute Horizontal Error', gotIt ) 
    IF (.NOT. gotIt) THEN
       WRITE( Message, * ) 'Computing error using only horizontal velocity'
       CALL Info( 'FlowSolve', Message, Level=4 )
       OnlyHorizontalError = .TRUE.
    END IF

    NodeType2Variable => VariableGet( Solver % Mesh % Variables, 'ApproximationLevel' )
    IF ( ASSOCIATED( NodeType2Variable ) ) THEN
       NodeType2Perm    => NodeType2Variable % Perm
       NodeType2Values  => NodeType2Variable % Values
    ELSE
       CALL Fatal( 'FlowSolveSIAFS','Cannot find variable <ApproximationLevel>' )
    END IF

    NodeWiseErrorVariable => VariableGet( Solver % Mesh % Variables, 'SIAError' )
    IF ( ASSOCIATED(  NodeWiseErrorVariable ) ) THEN
       NodeWiseErrorPerm     => NodeWiseErrorVariable  % Perm
       NodeWiseErrorValues   => NodeWiseErrorVariable % Values
    ELSE
       CALL Fatal( 'FlowSolveSIAFS','Cannot find variable <SIAError>' )
    END IF

    TimeVar => VariableGet( Solver % Mesh % Variables, 'Timestep')
    Timestep = NINT(Timevar % Values(1))

    ErrorInSIA=FlowSolution-SIAVelPermuted

    WRITE(*,*) 'a'
    NodeWiseError=0.0_dp
    WRITE(*,*) '1'

    !if error is to big then do FS
    NodeType2=0
    NumberOfFSNodes=0
    NumberOfSIANodes=0

    WRITE(*,*) 'b'

    DO i = 1, GetNOFActive()
       Element => GetActiveElement(i)
       n = GetElementNOFNodes()
       !     
       DO j=1,GetElementNOFNOdes()
          k = Element % NodeIndexes(j) 

          IF (NodeType2(k)/=0) CYCLE

          WRITE(*,*) 'hej'
          SELECT CASE( NSDOFs )
          CASE(3) !2D simulation
             IF (ABS(FlowSolution(NSDOFs*(FlowPerm(k)-1)+1))>LowerHorVelLimit) THEN

                NodeWiseError(k)=ABS(ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+1)&
                     /FlowSolution(NSDOFs*(FlowPerm(k)-1)+1))

                !Weighing absolute and relative error
                ErrorBound = MAXVAL((/ErrorBoundAbs &
                     /ABS(FlowSolution(NSDOFs*(FlowPerm(k)-1)+1)), ErrorBoundRel/))


             ELSE
                NodeWiseError(k)=0.0_dp
             END IF
             WRITE(*,*) 'd'

             WRITE(*,*) 'hej 1'

          CASE(4) !3D simulation

             IF (ABS(FlowSolution(NSDOFs*(FlowPerm(k)-1)+1))>LowerHorVelLimit .AND. &
                  ABS(FlowSolution(NSDOFs*(FlowPerm(k)-1)+2))>LowerHorVelLimit) THEN
                NodeWiseError(k)=0.5*SQRT( (ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+1)&
                     /FlowSolution(NSDOFs*(FlowPerm(k)-1)+1))**2.0+ &
                     (ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+2)&
                     /FlowSolution(NSDOFs*(FlowPerm(k)-1)+2))**2.0 )  
             ELSE IF (ABS(FlowSolution(NSDOFs*(FlowPerm(k)-1)+1))>LowerHorVelLimit) THEN
                NodeWiseError(k)=ABS(ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+1)&
                     /FlowSolution(NSDOFs*(FlowPerm(k)-1)+1))
             ELSE IF (ABS(FlowSolution(NSDOFs*(FlowPerm(k)-1)+2))>LowerHorVelLimit) THEN
                NodeWiseError(k)=ABS(ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+2)&
                     /FlowSolution(NSDOFs*(FlowPerm(k)-1)+2))
             ELSE
                NodeWiseError(k)=0.0_dp
             END IF

             !Weighing absolute and relative error
             ErrorBound = MAXVAL((/ ErrorBoundAbs/SQRT(FlowSolution(NSDOFs*(FlowPerm(k)-1)+1)**2.0 &
                  +FlowSolution(NSDOFs*(FlowPerm(k)-1)+2)**2.0), ErrorBoundRel/))

          END SELECT !select dimension

          WRITE(*,*) 'c'


!!!!! SORT NODES

          IF (NodeWiseError(k)> ErrorBound) THEN !FS-Node 
             NumberOfFSNodes=NumberOfFSNodes+1
             NodeType2(k) = 2
             NodeType2Values(NodeType2Perm(k))=REAL(NodeType2(k))
          ELSE    !SIA-Node
             NumberOfSIANodes=NumberOfSIANodes+1
             NodeType2(k) = 1
             NodeType2Values(NodeType2Perm(k))=REAL(NodeType2(k))
          END IF

          NodeWiseErrorValues(NodeWiseErrorPerm(k))=NodeWiseError(k)

       END DO
    END DO

    !FlowSolution=CoupledSolution !necessary??


  END SUBROUTINE SolutionErrorEstimate

  !-------------------------------------------------------------------------------------

  SUBROUTINE ResidualEstimate( Model,Solver,dt,TransientSimulation,SIAVelPermuted, &
       HastighetsError,NodeType2,NumberOfSIANodes,NumberOfFSNodes,ReorderTimeInterval)
    !******************************************************************************
    !
    !  Estimate Approximation Error Based on Residual
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
    !     INPUT: Timestep size for time dependent simulations
    !
    !******************************************************************************
    !------------------------------------------------------------------------------
    USE DefUtils

    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Solver_t), TARGET :: Solver

    REAL(KIND=dp) :: dt
    LOGICAL :: TransientSimulation
    INTEGER, ALLOCATABLE, intent(inout) :: NodeType2(:)
    REAL(KIND=dp),ALLOCATABLE, intent(in) :: SIAVelPermuted(:), HastighetsError(:)
    INTEGER, intent(out) :: NumberOfSIANodes, NumberOfFSNodes
    INTEGER, intent(in) :: ReorderTimeInterval
    !------------------------------------------------------------------------------
    !    Local variables
    !------------------------------------------------------------------------------

    TYPE(Element_t),POINTER :: Element

    REAL(KIND=dp),ALLOCATABLE :: Ax(:), NodeWiseResidual(:),residual(:)

    REAL(KIND=dp) :: ResidualBound, RelationForError, ErrorBound

    INTEGER :: ACounter, i, j, k, n, NSDOFs, istat

    LOGICAL :: UserResidualBound, RelativeResidual, OnlyHorizontalError
    LOGICAL :: AllocationsDone = .FALSE.
    LOGICAL :: gotIt


    TYPE(Variable_t), POINTER :: FlowSol 
    INTEGER, POINTER :: FlowPerm(:) 
    REAL(KIND=dp), POINTER :: FlowSolution(:) 
    INTEGER, SAVE :: Timestep 
    TYPE(Variable_t), POINTER :: TimeVar

    SAVE Ax,residual,NodeWiseResidual, ResidualBound, OnlyHorizontalError, &
         ErrorBound, UserResidualBound, RelationForError, RelativeResidual 


    !------------------------------------------------------------------------------

    TimeVar => VariableGet( Solver % Mesh % Variables, 'Timestep')
    Timestep = NINT(Timevar % Values(1))

    FlowSol => Solver % Variable
    NSDOFs         =  FlowSol % DOFs
    FlowPerm       => FlowSol % Perm
    FlowSolution   => FlowSol % Values


    !------------------------------------------------------------------------------
    !     Allocate some permanent storage, this is done first time only
    !------------------------------------------------------------------------------

    IF ( .NOT.AllocationsDone .OR. Solver % Mesh % Changed ) THEN

       N = Solver % Mesh % MaxElementDOFs

       IF( AllocationsDone ) THEN
          DEALLOCATE( &
               residual, &
               NodeWiseResidual, &
               STAT=istat )
       END IF
       ALLOCATE(  Ax(  SIZE( FlowSolution )), &
            residual(  SIZE( FlowSolution )), &
            NodeWiseResidual(  Model % Mesh % NumberOfNodes ), &
            STAT=istat )    
       IF ( istat /= 0 ) THEN
          CALL Fatal( 'FlowSolve','Memory allocation error, Aborting.' )
       END IF

       AllocationsDone = .TRUE.
    END IF  !------------------------------------------------------------------------------


    ResidualBound = GetConstReal(  Solver % Values,  &
         'Maximum Allowed Residual', UserResidualBound )
    RelativeResidual = GetLogical(  Solver % Values,  &
         'Relative Residual instead of Absolute', gotIt )
    ErrorBound = 0.01*GetConstReal(  Solver % Values,  &
         'Relative Error Allowed In Percent', gotIt )
    OnlyHorizontalError = GetLogical( Solver % Values, &
         'Only Compute Horizontal Error', gotIt ) 

!!!---------------------------------------------------

    WRITE( Message, * ) 'Computing the residual'
    CALL Info( 'FlowSolve', Message, Level=4 )


    IF ( ParEnv % PEs > 1 ) THEN !!!!!!!!!!!!!!!!!!!!!! we have a parallel run

       WRITE(*,*) 'PARALLELL IMPLEMENTATION MISSING FOR RESIDUAL COMPUTATION...'

    ELSE !!!!!!!!!!!!!!!!!!!!!! serial run

       WRITE(*,*) 'Multiplying Ax'

       Ax=0.0

       CALL CRS_MatrixVectorMultiply(Solver % Matrix,SIAVelPermuted,Ax)

    END IF

    WRITE(*,*) 'computing Ax-b'
    residual = Ax - Solver % Matrix % RHS


    WRITE(*,*) 'computing nodewise residual'

    NodeWiseResidual=1.0E25 !Something unreasonable
    RelationForError=0
    ACounter=0


    DO i = 1, GetNOFActive()
       Element => GetActiveElement(i)
       n = GetElementNOFNodes()
       !     

       DO j=1,GetElementNOFNOdes()
          k = Element % NodeIndexes(j)

          IF(NodeWiseResidual(k)/=1.0E25) CYCLE !already computed
          !IF (NodeType2(k)/=0) CYCLE
          !
          SELECT CASE( NSDOFs )
          CASE(3) !2D simulation
             IF (OnlyHorizontalError) THEN !only x-velocity
                IF (RelativeResidual) THEN
                   NodeWiseResidual(k)=SQRT((residual(NSDOFs*FlowPerm(k)-2)/FlowSolution(NSDOFs*FlowPerm(k)-2))**2.0)
                ELSE
                   NodeWiseResidual(k)=SQRT(residual(NSDOFs*FlowPerm(k)-2)**2.0)
                END IF
             ELSE !All velocity componentts
                IF (k==1) THEN
                   WRITE(*,*) 'RelativeResidual', RelativeResidual
                END IF
                IF (RelativeResidual) THEN
                   !WRITE(*,*) 'NSDOFs=',NSDOFs
		   !WRITE(*,*) 'FlowSolution(NSDOFs*FlowPerm(k)-1)', FlowSolution(NSDOFs*FlowPerm(k)-1)
                   NodeWiseResidual(k)=SQRT((residual(NSDOFs*FlowPerm(k)-2)/FlowSolution(NSDOFs*FlowPerm(k)-2))**2.0+&
                        (residual(NSDOFs*FlowPerm(k)-1)/FlowSolution(NSDOFs*FlowPerm(k)-1))**2.0)
                   IF ( ISNAN(NodeWiseResidual(k)) ) THEN
                      NodeWiseResidual(k)=0.0
                   END IF
                   IF ( NodeWiseResidual(k)>=huge(NodeWiseResidual(k)) ) THEN
                      NodeWiseResidual(k)=100.0*ResidualBound
                   END IF

                ELSE  
                   NodeWiseResidual(k)=SQRT(residual(NSDOFs*FlowPerm(k)-1)**2.0+ &
                        residual(NSDOFs*FlowPerm(k)-2)**2.0)
                END IF

             END IF
          CASE(4) !3D simulation              
             IF (OnlyHorizontalError) THEN !only x-velocity
                IF (RelativeResidual) THEN
                   NodeWiseResidual(k)=SQRT( (residual(NSDOFs*FlowPerm(k)-3)/FlowSolution(NSDOFs*FlowPerm(k)-3))**2.0 +&
                        (residual(NSDOFs*FlowPerm(k)-2)/FlowSolution(NSDOFs*FlowPerm(k)-2))**2.0)
                ELSE
                   NodeWiseResidual(k)=SQRT(residual(NSDOFs*FlowPerm(k)-3)**2.0+residual(NSDOFs*FlowPerm(k)-2)**2.0)
                END IF
             ELSE !All velocity componentts
                IF (k==1) THEN
                   WRITE(*,*) 'RelativeResidual', RelativeResidual
                END IF
                IF (RelativeResidual) THEN
                   !WRITE(*,*) 'NSDOFs=',NSDOFs
		   !WRITE(*,*) 'FlowSolution(NSDOFs*FlowPerm(k)-1)', FlowSolution(NSDOFs*FlowPerm(k)-1)
                   NodeWiseResidual(k)=SQRT((residual(NSDOFs*FlowPerm(k)-3)/FlowSolution(NSDOFs*FlowPerm(k)-3))**2.0+&
                        (residual(NSDOFs*FlowPerm(k)-2)/FlowSolution(NSDOFs*FlowPerm(k)-2))**2.0+&
			(residual(NSDOFs*FlowPerm(k)-1)/FlowSolution(NSDOFs*FlowPerm(k)-1))**2.0)
                   IF ( ISNAN(NodeWiseResidual(k)) ) THEN
                      NodeWiseResidual(k)=0.0
                   END IF
                   IF ( NodeWiseResidual(k)>=huge(NodeWiseResidual(k)) ) THEN
                      NodeWiseResidual(k)=100.0*ResidualBound
                   END IF

                ELSE  
                   NodeWiseResidual(k)=SQRT(residual(NSDOFs*FlowPerm(k)-3)**2.0+ &
                        residual(NSDOFs*FlowPerm(k)-2)**2.0+residual(NSDOFs*FlowPerm(k)-1)**2.0)
                END IF

             END IF
          END SELECT



          IF (.NOT. UserResidualBound) THEN
             !Compute relation between residual and error in velocity in FS-field
             IF (NodeType2(k)==2) THEN !FS node
                SELECT CASE( NSDOFs )
		CASE(3)
                   IF (HastighetsError(k)>=ErrorBound) THEN
                      IF (.NOT. ISNAN(NodeWiseResidual(k)) .AND. .NOT. ISNAN(HastighetsError(k)) &
                           .AND. HastighetsError(k)/=0.0 .AND. FlowSolution(NSDOFs*(FlowPerm(k)-1)+1)/=0 &
                           .AND. FlowSolution(NSDOFs*(FlowPerm(k)-1)+2)/=0 ) THEN
                         IF (ABS(NodeWiseResidual(k)/HastighetsError(k))<=huge(NodeWiseResidual(k)/HastighetsError(k))) THEN
                            RelationForError=RelationForError+NodeWiseResidual(k)/HastighetsError(k)
                         END IF
                         ACounter=ACounter+1
                      END IF
                   END IF
		CASE(4)   
                   IF (HastighetsError(k)>=ErrorBound) THEN
                      IF (.NOT. ISNAN(NodeWiseResidual(k)) .AND. .NOT. ISNAN(HastighetsError(k)) &
                           .AND. HastighetsError(k)/=0.0 .AND. FlowSolution(NSDOFs*(FlowPerm(k)-1)+1)/=0 &
                           .AND. FlowSolution(NSDOFs*(FlowPerm(k)-1)+2)/=0 .AND. FlowSolution(NSDOFs*(FlowPerm(k)-1)+3)/=0 ) THEN
                         IF (ABS(NodeWiseResidual(k)/HastighetsError(k))<=huge(NodeWiseResidual(k)/HastighetsError(k))) THEN
                            RelationForError=RelationForError+NodeWiseResidual(k)/HastighetsError(k)
                         END IF
                         ACounter=ACounter+1
                      END IF
                   END IF
                END SELECT
             END IF
          END IF

       END DO

    END DO




    IF (.NOT. UserResidualBound) THEN
       RelationForError=RelationForError/ACounter
       ResidualBound=RelationForError*ErrorBound
    END IF

    WRITE(*,*) ' '
    WRITE(*,*) 'Relation for error: ',RelationForError, 'giving a residual bound of ' ,ResidualBound
    WRITE(*,*) ' '

    NodeType2=0
    NumberOfFSNodes=0
    NumberOfSIANodes=0


    DO i = 1, GetNOFActive()
       Element => GetActiveElement(i)
       n = GetElementNOFNodes()
       !     
       DO j=1,GetElementNOFNOdes()
          k = Element % NodeIndexes(j)

          IF (NodeType2(k)/=0) CYCLE
          !
          IF (NodeWiseResidual(k)> ResidualBound) THEN !FS-Node
             NumberOfFSNodes=NumberOfFSNodes+1
             NodeType2(k) = 2
             ! WRITE(*,*) 'Node ', k, 'is FS and NodeWiseResidual is', NodeWiseResidual(k)
          ELSE    !SIA-Node
             NumberOfSIANodes=NumberOfSIANodes+1
             NodeType2(k) = 1
             ! WRITE(*,*) 'Node ', k, 'is SIA and NodeWiseResidual is', NodeWiseResidual(k)
          END IF

          !Compute relation between residual and error in velocity
       END DO
    END DO

    WRITE(*,*) 'Printing residual to Residual.dat'

    IF ( Timestep==ReorderTimeInterval) THEN !First time residual is computed
       open (unit=91, file="Residual.dat",STATUS='NEW',POSITION='APPEND')
       !Compute the difference between SIA solution and FS solution in FS-nodes
    ELSE
       open (unit=91, file="Residual.dat",STATUS='OLD',POSITION='APPEND')
    END IF

    WRITE(91,*) 'Timestep=',Timestep
    DO i = 1, SIZE(FlowSolution)/NSDOFs
       WRITE(91,*) i,NodeWiseResidual(i)
    END DO

    close (91)

  END SUBROUTINE ResidualEstimate


END MODULE ErrorEstimationSubs
