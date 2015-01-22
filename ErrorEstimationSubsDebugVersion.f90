MODULE ErrorEstimationSubs
     USE SparIterSolve
     USE ParallelUtils
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
    CALL Info( 'Error Estimation', Message, Level=4 )

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

    NodeWiseError=0.0_dp

    !if error is to big then do FS
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

          CASE(4) !3D simulation

             NodeWiseError(k)=0.5*SQRT( (ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+1)**2.0+ ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+2)**2.0 ) &
                  /(FlowSolution(NSDOFs*(FlowPerm(k)-1)+1)**2.0+FlowSolution(NSDOFs*(FlowPerm(k)-1)+2)**2.0))

      !       IF (ABS(FlowSolution(NSDOFs*(FlowPerm(k)-1)+1))>LowerHorVelLimit .AND. &
      !            ABS(FlowSolution(NSDOFs*(FlowPerm(k)-1)+2))>LowerHorVelLimit) THEN
      !          NodeWiseError(k)=0.5*SQRT( (ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+1)&
      !               /FlowSolution(NSDOFs*(FlowPerm(k)-1)+1))**2.0+ &
      !               (ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+2)&
      !               /FlowSolution(NSDOFs*(FlowPerm(k)-1)+2))**2.0 )  
      !       ELSE IF (ABS(FlowSolution(NSDOFs*(FlowPerm(k)-1)+1))>LowerHorVelLimit) THEN
      !          NodeWiseError(k)=ABS(ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+1)&
      !               /FlowSolution(NSDOFs*(FlowPerm(k)-1)+1))
      !       ELSE IF (ABS(FlowSolution(NSDOFs*(FlowPerm(k)-1)+2))>LowerHorVelLimit) THEN
      !          NodeWiseError(k)=ABS(ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+2)&
      !               /FlowSolution(NSDOFs*(FlowPerm(k)-1)+2))
      !       ELSE
      !          NodeWiseError(k)=0.0_dp
      !       END IF

             !Weighing absolute and relative error
             ErrorBound = MAXVAL((/ ErrorBoundAbs/SQRT(FlowSolution(NSDOFs*(FlowPerm(k)-1)+1)**2.0 &
                  +FlowSolution(NSDOFs*(FlowPerm(k)-1)+2)**2.0), ErrorBoundRel/))

          END SELECT !select dimension

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


  SUBROUTINE FunctionalErrorEstimate( Model,Solver,dt,TransientSimulation, &
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
    USE Functionals

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
    TYPE(Matrix_t),POINTER :: A, AT

    INTEGER :: ACounter, i, j, k, n, NSDOFs, istat 
    LOGICAL :: OnlyHorizontalError
    LOGICAL :: AllocationsDone = .FALSE.
    LOGICAL :: gotIt

    REAL(KIND=dp),ALLOCATABLE :: x(:), &
         NodeWiseError(:), &
         Ax(:),residual(:), functional(:)

    REAL(KIND=dp) :: UNorm

    TYPE(Variable_t), POINTER :: FlowSol, NodeType2Variable
    TYPE(Variable_t), POINTER ::  NodeWiseErrorVariable
    INTEGER, POINTER :: FlowPerm(:), NodeType2Perm(:), NodeWiseErrorPerm(:)
    REAL(KIND=dp), POINTER :: FlowSolution(:), ForceVector(:), &
         NodeType2Values(:), NodeWiseErrorValues(:)
    INTEGER, SAVE :: Timestep 
    TYPE(Variable_t), POINTER :: TimeVar

    REAL(KIND=dp), SAVE :: ErrorBound, adv,xdf,av,xf

    REAL(KIND=dp), POINTER :: xx(:)
    TYPE(Matrix_t), POINTER :: ss
    INTEGER, POINTER :: pp(:)

    REAL(KIND=dp), POINTER :: functionalpointer(:)

    CHARACTER(LEN=MAX_NAME_LEN) :: FunctionalName

    TARGET :: x, functional

    SAVE OnlyHorizontalError, &
         NodeWiseError, AT, functional, functionalpointer, &
         pp, ss, xx, &
         Ax, residual, x, FunctionalName

    !-----------------------------------------------------------------
    ! INITIALIZATIONS AND ALLOCATIONS
    !----------------------------------------------------------------

    WRITE( Message, * ) 'Computing the Error'
    CALL Info( 'FlowSolve', Message, Level=4 )

    FlowSol => Solver % Variable
    NSDOFs         =  FlowSol % DOFs
    FlowPerm       => FlowSol % Perm
    FlowSolution   => FlowSol % Values

    A => Solver % Matrix

    AT => AllocateMatrix()  !Will be transponate of A
    AT % Format = MATRIX_LIST

    IF ( .NOT.AllocationsDone .OR. Solver % Mesh % Changed ) THEN
       IF( AllocationsDone ) THEN
          DEALLOCATE(                               &
               NodeWiseError, &
               AT % Rows, &
               AT % Values, &
               AT % Cols, &
               x, &
               functional, &
               Ax, &
               STAT=istat )
       END IF
       ALLOCATE( &
            NodeWiseError(  Model % Mesh % NumberOfNodes ), &  
            AT % Rows(SIZE(A % Rows)), &
            AT % Values(SIZE(A % Values)), &
            AT % Cols(SIZE(A % Cols)), &
            x(SIZE(A % RHS)), &
            functional(SIZE(A % RHS)), &
            Ax(SIZE(A % RHS)), &
            STAT=istat)
       AllocationsDone = .TRUE.
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



    FunctionalName  = GetString(  Solver % Values,  &
         'Functional', gotIt )
    IF (.NOT. gotIt) THEN
       CALL Fatal( 'Error Estimation: ', 'Functional not chosen, aborting.')
    END IF
    !-----------------------------------------------------------------------------
    !   GET TRANSPONATE OF SYSTEM MATRIX A^T, and the functional
    !-----------------------------------------------------------------------------

    !First copy A to AT
    AT = A

    !Transpose AT
    AT = CRS_Transpose(AT)

    !Allocate the right hand side
    ALLOCATE(AT % RHS(SIZE(A % RHS)))

    functional = 0.0
    functionalpointer => functional

    !Get the functional  
    SELECT CASE(FunctionalName)
    CASE('flux across point') 
       CALL FluxThroughLine( Model,Solver,dt,TransientSimulation, &
            functionalpointer)!  1.0_dp
    CASE DEFAULT
       Call FATAL('Error Estimation', 'No valid functional chosen')
    END SELECT

    AT % RHS = functional

    !-----------------------------------------------------------------------------
    !   SOLVE A^T x_imp = a 
    !-----------------------------------------------------------------------------

    !Save old system matrix
    pp => Solver % Variable % Perm
    ss => Solver % Matrix
    xx => Solver % Variable % Values

    !Set system matrix to AT, and values to x
    Solver % Matrix => AT 
    A => Solver % Matrix
    Solver % Variable % Values => x


    !Prepare some parallell stuff
    IF (ParEnv % PEs>1) THEN
       A % Comm = MPI_COMM_WORLD
       IF(.NOT.ASSOCIATED(A % ParMatrix)) THEN              
          CALL ParallelInitMatrix(Solver,A,FlowPerm)
       END IF
    END IF

    !Solve AT x = a
    Solver % Variable  % Values=0._dp
    UNorm = DefaultSolve() 

    !Reset system matrix to the old one ... did I really get something in x now?
    Solver % Matrix => ss
    A => Solver % Matrix
    Solver % Variable % Perm => pp
    Solver % Variable % Values => xx 

    !-----------------------------------------------------------------------------
    !   COMPUTE residual
    !-----------------------------------------------------------------------------

    IF ( ParEnv % PEs > 1 ) THEN ! we have a parallel run
       ss => A
       A => Solver % Matrix
       Solver % Matrix % Comm = MPI_COMM_WORLD

       IF(.NOT.ASSOCIATED(A % ParMatrix)) CALL ParallelInitMatrix(Solver,A,FlowPerm)

       CALL ParallelInitSolve(A,SIAVelPermuted,Ax,Ax)
       CALL ParallelMatrixVector( A, SIAVelPermuted, Ax, .TRUE. )
       CALL ParallelSumVector(A, Ax)

       Solver % Matrix => ss
       A => Solver % Matrix
    ELSE ! serial run

       CALL CRS_MatrixVectorMultiply(Solver % Matrix,SIAVelPermuted,Ax)
    END IF

    residual = Ax - Solver % Matrix % RHS

    !-----------------------------------------------------------------------------
    !   Compute x^T*residual (if functional error is too big)
    !-----------------------------------------------------------------------------

    adv=0
    xdf=0
    av=0
    xf=0
    !compute elements of x^T*residual
    DO i = 1, Model % Mesh % NumberOfNodes
       NodeWiseErrorValues(NodeWiseErrorPerm(i))=0.0
       DO k=1, NSDOFs
          NodeWiseErrorValues(NodeWiseErrorPerm(i))=NodeWiseErrorValues(NodeWiseErrorPerm(i))+ &
               ABS(x(NSDOFs*(FlowPerm(i)-1)+k)*residual(NSDOFs*(FlowPerm(i)-1)+k))

          adv= adv+functional(NSDOFs*(FlowPerm(i)-1)+k)*(SIAVelPermuted(NSDOFs*(FlowPerm(i)-1)+k)- &
               FlowSolution(NSDOFs*(FlowPerm(i)-1)+k))
          xdf=xdf+ x(NSDOFs*(FlowPerm(i)-1)+k)*residual(NSDOFs*(FlowPerm(i)-1)+k)

          av= av+functional(NSDOFs*(FlowPerm(i)-1)+k)*FlowSolution(NSDOFs*(FlowPerm(i)-1)+k)
          xf=xf+ x(NSDOFs*(FlowPerm(i)-1)+k)*Solver % Matrix % RHS(NSDOFs*(FlowPerm(i)-1)+k)

       END DO
    END DO

    !-----------------------------------------------------------------------------
    !   Sort Nodes
    !-----------------------------------------------------------------------------

    open (unit=135, file="FunctionalStuff.dat",POSITION='APPEND')
    WRITE(135,*)  Timestep

    WRITE(135,*)  adv
    WRITE(135,*)  av
    WRITE(135,*)  xdf
    WRITE(135,*)  xf

    ErrorBound = GetConstReal( Solver % Values, 'Nodewise limit for dual problem', gotIt )    
    IF (.NOT. gotIt) THEN
       WRITE( Message, * ) 'Nodewise limit for dual problem not found, setting to 5.0'
       CALL Info( 'Error Estimation: ',Message, Level=4 )
       ErrorBound = 5.0
    END IF

    NodeType2=0
    NumberOfFSNodes=0
    NumberOfSIANodes=0
    DO i = 1, Model % Mesh % NumberOfNodes
       !if error is to big then do FS
       IF (NodeWiseErrorValues(NodeWiseErrorPerm(i))> ErrorBound) THEN !FS-Node 
          NumberOfFSNodes=NumberOfFSNodes+1
          NodeType2(i) = 2
          NodeType2Values(NodeType2Perm(i))=REAL(NodeType2(i))
       ELSE    !SIA-Node
          NumberOfSIANodes=NumberOfSIANodes+1
          NodeType2(i) = 1
          NodeType2Values(NodeType2Perm(i))=REAL(NodeType2(i))
       END IF

       WRITE(135,*) i,NodeWiseErrorValues(NodeWiseErrorPerm(i))

    END DO

    WRITE(135,*) '***************************************************************'
    WRITE(135,*) '                                                               '

    close(135)  

    DEALLOCATE(AT % RHS)

  END SUBROUTINE FunctionalErrorEstimate


  !-------------------------------------------------------------------

  SUBROUTINE ResidualEstimate( Model,Solver,dt,TransientSimulation, &
      NodeType2,SIAVelPermuted,NumberOfSIANodes,NumberOfFSNodes)
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
    REAL(KIND=dp),ALLOCATABLE, intent(in) :: SIAVelPermuted(:)
    INTEGER, intent(out) :: NumberOfSIANodes, NumberOfFSNodes
    !------------------------------------------------------------------------------
    !    Local variables
    !------------------------------------------------------------------------------

    TYPE(Element_t),POINTER :: Element

    TYPE(Matrix_t),POINTER :: A

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

    TYPE(Variable_t), POINTER :: NodeType2Variable
    TYPE(Variable_t), POINTER ::  NodeWiseErrorVariable
    INTEGER, POINTER :: NodeType2Perm(:), NodeWiseErrorPerm(:)
    REAL(KIND=dp), POINTER :: NodeType2Values(:), NodeWiseErrorValues(:)

    REAL(KIND=dp), POINTER :: xx(:)
    TYPE(Matrix_t), POINTER :: ss
    INTEGER, POINTER :: pp(:)

    SAVE Ax,residual,NodeWiseResidual, ResidualBound, OnlyHorizontalError, &
         ErrorBound, UserResidualBound, RelationForError, RelativeResidual 

    !------------------------------------------------------------------------------

    TimeVar => VariableGet( Solver % Mesh % Variables, 'Timestep')
    Timestep = NINT(Timevar % Values(1))

    FlowSol => Solver % Variable
    NSDOFs         =  FlowSol % DOFs
    FlowPerm       => FlowSol % Perm
    FlowSolution   => FlowSol % Values
    
    A => Solver % Matrix
    !-----------------------------------------------------------------
    ! INITIALIZATIONS AND ALLOCATIONS
    !-----------------------------------------------------------------

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
    END IF  

    ResidualBound = GetConstReal(  Solver % Values,  &
         'Maximum Allowed Residual', UserResidualBound )
   
    NodeType2Variable => VariableGet( Solver % Mesh % Variables, 'ApproximationLevel' )
    IF ( ASSOCIATED( NodeType2Variable ) ) THEN
       NodeType2Perm    => NodeType2Variable % Perm
       NodeType2Values  => NodeType2Variable % Values
    ELSE
       CALL Fatal( 'Error Estimate','Cannot find variable <ApproximationLevel>' )
    END IF

    NodeWiseErrorVariable => VariableGet( Solver % Mesh % Variables, 'SIAError' )
    IF ( ASSOCIATED(  NodeWiseErrorVariable ) ) THEN
       NodeWiseErrorPerm     => NodeWiseErrorVariable  % Perm
       NodeWiseErrorValues   => NodeWiseErrorVariable % Values
    ELSE
       CALL Fatal( 'Error Estimate','Cannot find variable <SIAError>' )
    END IF

    !-----------------------------------------------------------------------------
    !   COMPUTE residual
    !-----------------------------------------------------------------------------

    IF ( ParEnv % PEs > 1 ) THEN ! we have a parallel run
       ss => A
       A => Solver % Matrix
       Solver % Matrix % Comm = MPI_COMM_WORLD

       IF(.NOT.ASSOCIATED(A % ParMatrix)) CALL ParallelInitMatrix(Solver,A,FlowPerm)

       CALL ParallelInitSolve(A,SIAVelPermuted,Ax,Ax)
       CALL ParallelMatrixVector( A, SIAVelPermuted, Ax, .TRUE. )
       CALL ParallelSumVector(A, Ax)

       Solver % Matrix => ss
       A => Solver % Matrix
    ELSE ! serial run
       CALL CRS_MatrixVectorMultiply(Solver % Matrix,SIAVelPermuted,Ax)
    END IF

    residual = Ax - Solver % Matrix % RHS

    WRITE( Message, * ) 'Computing the residual'
    CALL Info( 'FlowSolve', Message, Level=4 )

    !-----------------------------------------------------------------------------
    !   COMPUTE nodewise residual
    !-----------------------------------------------------------------------------

    NodeWiseResidual=1.0E25 !Something unreasonable
    RelationForError=0
    ACounter=0

    DO i = 1, GetNOFActive()
       Element => GetActiveElement(i)
       n = GetElementNOFNodes()

       DO j=1,n
          k = Element % NodeIndexes(j)

          IF(NodeWiseResidual(k)/=1.0E25) CYCLE !already computed
          !IF (NodeType2(k)/=0) CYCLE
          !
          SELECT CASE( NSDOFs )

          CASE(3) !2D simulation
             NodeWiseResidual(k)=SQRT(residual(NSDOFs*(FlowPerm(k)-1)+1)**2.0 + &
                  residual(NSDOFs*(FlowPerm(k)-1)+2)**2.0 + &
                  residual(NSDOFs*(FlowPerm(k)-1)+3)**2.0)

          CASE(4) !3D simulation              
                NodeWiseResidual(k)=SQRT(residual(NSDOFs*(FlowPerm(k)-1)+1)**2.0 + &
                  residual(NSDOFs*(FlowPerm(k)-1)+2)**2.0 + &
                  residual(NSDOFs*(FlowPerm(k)-1)+3)**2.0 + &
                  residual(NSDOFs*(FlowPerm(k)-1)+4)**2.0)

          END SELECT
       END DO
    END DO

    !-----------------------------------------------------------------------------
    !   Sort Nodes
    !-----------------------------------------------------------------------------

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
             NodeType2Values(NodeType2Perm(k))=REAL(NodeType2(k))

          ELSE    !SIA-Node
             NumberOfSIANodes=NumberOfSIANodes+1
             NodeType2(k) = 1
             NodeType2Values(NodeType2Perm(k))=REAL(NodeType2(k))

          END IF
          NodeWiseErrorValues(NodeWiseErrorPerm(k))=NodeWiseResidual(k)
       END DO
    END DO

  END SUBROUTINE ResidualEstimate




  SUBROUTINE SaveErrorMeasures( Model,Solver,dt,TransientSimulation, &
      SIAVelPermuted)
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
    USE Functionals

    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Solver_t), TARGET :: Solver

    REAL(KIND=dp) :: dt
    LOGICAL :: TransientSimulation

    REAL(KIND=dp),ALLOCATABLE, intent(in) :: SIAVelPermuted(:)

    !------------------------------------------------------------------------------
    !    Local variables
    !------------------------------------------------------------------------------
    TYPE(Element_t),POINTER :: Element
    TYPE(Matrix_t),POINTER :: A

    INTEGER :: ACounter, i, j, k, n, NSDOFs, istat 
    LOGICAL :: AllocationsDone = .FALSE.
    LOGICAL :: gotIt

    REAL(KIND=dp),ALLOCATABLE :: Ax(:),residual(:), &
         functional(:), nodewiseresidual(:)

    REAL(KIND=dp) :: UNorm

    TYPE(Variable_t), POINTER :: FlowSol
    INTEGER, POINTER :: FlowPerm(:)
    REAL(KIND=dp), POINTER :: FlowSolution(:)
    INTEGER, SAVE :: Timestep 
    TYPE(Variable_t), POINTER :: TimeVar

    REAL(KIND=dp), SAVE :: av

    REAL(KIND=dp), POINTER :: xx(:)
    TYPE(Matrix_t), POINTER :: ss
    INTEGER, POINTER :: pp(:)

    REAL(KIND=dp), POINTER :: functionalpointer(:)

    CHARACTER(LEN=MAX_NAME_LEN) :: FunctionalName

    TARGET :: functional

    CHARACTER(LEN=MAX_NAME_LEN) :: TimeFileName

    SAVE functional, functionalpointer, &
         pp, ss, xx, &
         Ax, residual, FunctionalName, nodewiseresidual

    !-----------------------------------------------------------------
    ! INITIALIZATIONS AND ALLOCATIONS
    !----------------------------------------------------------------

    WRITE( Message, * ) 'Computing the Error'
    CALL Info( 'FlowSolve', Message, Level=4 )

    FlowSol => Solver % Variable
    NSDOFs         =  FlowSol % DOFs
    FlowPerm       => FlowSol % Perm
    FlowSolution   => FlowSol % Values

    A => Solver % Matrix

    IF ( .NOT.AllocationsDone .OR. Solver % Mesh % Changed ) THEN
       IF( AllocationsDone ) THEN
          DEALLOCATE(                               &
               nodewiseresidual, &
               functional, &
               Ax, &
               STAT=istat )
       END IF
       ALLOCATE( &
            nodewiseresidual(  Model % Mesh % NumberOfNodes ), &  
            functional(SIZE(A % RHS)), &
            Ax(SIZE(A % RHS)), &
            STAT=istat)
       AllocationsDone = .TRUE.
    END IF


    TimeVar => VariableGet( Solver % Mesh % Variables, 'Timestep')
    Timestep = NINT(Timevar % Values(1))

    FunctionalName  = GetString(  Solver % Values,  &
         'Functional', gotIt )
    IF (.NOT. gotIt) THEN
       CALL Fatal( 'Error Estimation: ', 'Functional not chosen, aborting.')
    END IF
    !-----------------------------------------------------------------------------
    !   Get the functional
    !-----------------------------------------------------------------------------

    functional = 0.0
    functionalpointer => functional

    !Get the functional  
    SELECT CASE(FunctionalName)
    CASE('flux across point') 
       CALL FluxThroughLine( Model,Solver,dt,TransientSimulation, &
            functionalpointer)!  1.0_dp
    CASE DEFAULT
       Call FATAL('Error Estimation', 'No valid functional chosen')
    END SELECT

  
    !-----------------------------------------------------------------------------
    !   COMPUTE residual
    !-----------------------------------------------------------------------------

    IF ( ParEnv % PEs > 1 ) THEN ! we have a parallel run
       ss => A
       A => Solver % Matrix
       Solver % Matrix % Comm = MPI_COMM_WORLD

       IF(.NOT.ASSOCIATED(A % ParMatrix)) CALL ParallelInitMatrix(Solver,A,FlowPerm)

       CALL ParallelInitSolve(A,SIAVelPermuted,Ax,Ax)
       CALL ParallelMatrixVector( A, SIAVelPermuted, Ax, .TRUE. )
       CALL ParallelSumVector(A, Ax)

       Solver % Matrix => ss
       A => Solver % Matrix
    ELSE ! serial run

       CALL CRS_MatrixVectorMultiply(Solver % Matrix,SIAVelPermuted,Ax)
    END IF

    residual = Ax - Solver % Matrix % RHS

    !-----------------------------------------------------------------------------
    !   Compute flux and nodewise residual
    !-----------------------------------------------------------------------------

    av=0
    nodewiseresidual=0.0
    DO i = 1, Model % Mesh % NumberOfNodes
       DO k=1, NSDOFs
          av= av+functional(NSDOFs*(FlowPerm(i)-1)+k)*FlowSolution(NSDOFs*(FlowPerm(i)-1)+k)
       END DO
       nodewiseresidual(i)= SQRT(residual(NSDOFs*(FlowPerm(i)-1)+1)**2.0 + &
                  residual(NSDOFs*(FlowPerm(i)-1)+2)**2.0 + &
                  residual(NSDOFs*(FlowPerm(i)-1)+3)**2.0)
    END DO

    !-----------------------------------------------------------------------------
    !   WRITE TO FILE
    !-----------------------------------------------------------------------------
    
    TimeFileName=GetString( Solver % Values, 'Error File Name', gotIt )

    open (unit=201, file=TimeFileName,POSITION='APPEND')

    WRITE(201,*) '***************************************************************'

    WRITE(201,*)  Timestep
    WRITE(201,*)  av

    DO i = 1, Model % Mesh % NumberOfNodes
       WRITE(201,*) i,nodewiseresidual(i)
    END DO

    WRITE(201,*) '***************************************************************'
    WRITE(201,*) '                                                               '

    close(201)  

  END SUBROUTINE SaveErrorMeasures




 SUBROUTINE ParallelInitMatrix2( Solver, Matrix, inPerm )
!-------------------------------------------------------------------------------
       INTEGER, POINTER, OPTIONAL :: inPerm(:)
       TYPE(Solver_t) :: Solver
       TYPE(Matrix_t), POINTER :: Matrix
!-------------------------------------------------------------------------------
       TYPE(ParallelInfo_t), POINTER :: MatrixPI, MeshPI
       INTEGER :: i, j, k, l, m, n, DOFs
       LOGICAL :: DGSolver, Found, GB, Global_dof
       TYPE(Mesh_t), POINTER :: Mesh
       TYPE(Element_t), POINTER :: Element
       TYPE(NeighbourList_t), POINTER :: MtrxN, MeshN

       INTEGER :: maxnode, maxedge, maxface, fdofs, maxfdofs, &
                  edofs, maxedofs, maxbdofs, l_beg, g_beg, proc_ind, sz

       INTEGER, POINTER :: Perm(:)

       INTEGER :: pardofs_bef, pardofs_all
       INTEGER, ALLOCATABLE :: pardofs_in(:), pardofs_out(:), ind(:)
       REAL(KIND=dp) :: realtime,tt
!-------------------------------------------------------------------------------

       IF ( ParEnv % PEs <= 1 .OR. .NOT. ASSOCIATED(Matrix) ) RETURN
       Mesh => Solver % Mesh
       DOFs = Solver % Variable % DOFs

       Perm => Solver % Variable % Perm
       IF(PRESENT(inPerm)) THEN
          WRITE(*,*) 'got a perm!'
          Perm=>InPerm
       END IF
       n = SIZE(Perm)
       k = n*DOFs + Matrix % ExtraDOFs

       WRITE(*,*) 'n=', n
       WRITE(*,*) 'k=', k
       ALLOCATE( Matrix % Perm(k), Matrix % InvPerm(k))

              WRITE(*,*) 'Processor no ', ParEnv % MyPE,'MAXVAL(Perm)=',MAXVAL(Perm)


       Matrix % Perm = 0
       DO i=1,n
         IF ( Perm(i) /= 0 )  THEN
            DO j=1,DOFs
               Matrix % Perm((i-1)*DOFs+j) = DOFs * (Perm(i)-1) + j
            END DO
         END IF
       END DO

        DO i=n*DOFs+1,SIZE(Matrix % Perm)
          Matrix % Perm(i) = i
        END DO
       Matrix % INVPerm = 0

              WRITE(*,*) 'Processor no ', ParEnv % MyPE,'SIZE(Matrix % Perm)=', SIZE(Matrix % Perm)

          WRITE(*,*) 'Processor no ', ParEnv % MyPE,'SIZE(Matrix % INVPerm)=', SIZE(Matrix % INVPerm)
              WRITE(*,*) 'Processor no ', ParEnv % MyPE,'finding inverese perm'

              WRITE(*,*) 'Processor no ', ParEnv % MyPE,'MAXVAL(Matrix % Perm)=',MAXVAL(Matrix % Perm)


       DO i=1,SIZE(Matrix % Perm)
          IF ( Matrix % Perm(i) /= 0 ) THEN
             Matrix % INVPerm(Matrix % Perm(i)) = i
          END IF
       END DO
              WRITE(*,*) 'Processor no ', ParEnv % MyPE,'bb'
       IF ( .NOT. Matrix % DGMatrix ) THEN
         n = Matrix % NumberOfRows
         ALLOCATE( Matrix % ParallelInfo )
         ALLOCATE( Matrix % ParallelInfo % NeighbourList(n) )
         CALL AllocateVector( Matrix % ParallelInfo % Interface, n)
         CALL AllocateVector( Matrix % ParallelInfo % GlobalDOFs, n)
              WRITE(*,*) 'Processor no ', ParEnv % MyPE,'cc'
         DO i=1,n
           Matrix % ParallelInfo % Interface(i) = .FALSE.
           Matrix % ParallelInfo % GlobalDOFs(i) = 0
           Matrix % ParallelInfo % NeighbourList(i) % Neighbours => NULL()
         END DO

         DO i=1,Solver % Mesh % NumberOfNodes
           DO j=1,DOFs
              k = Matrix % Perm((i-1)*DOFs+j)
              IF(k<=0) CYCLE
              Matrix % ParallelInfo % GlobalDOFs(k) = &
                DOFs*(Solver % Mesh % ParallelInfo % GlobalDOFs(i)-1)+j
              Matrix % ParallelInfo % Interface(k) = &
                Solver % Mesh % ParallelInfo % Interface(i)
              ALLOCATE( Matrix % ParallelInfo % NeighbourList(k) % Neighbours(SIZE( &
                   Solver % Mesh % ParallelInfo % NeighbourList(i) % Neighbours)) )
              Matrix % ParallelInfo % NeighbourList(k) % Neighbours = &
                Solver % Mesh % ParallelInfo % NeighbourList(i) % Neighbours
           END DO
         END DO
              WRITE(*,*) 'Processor no ', ParEnv % MyPE,'c'

         GB = ListGetLogical( Solver % Values, 'Bubbles in Global System', Found )
         IF (.NOT.Found) GB = .TRUE.

         maxnode = MAXVAL(Solver % Mesh % ParallelInfo % GlobalDOFs)
         maxnode = NINT(ParallelReduction(1._dp*maxnode,2))

         edofs = 0; fdofs = 0; maxedofs = 0; maxfdofs = 0
         maxedge = 0; maxface = 0

         IF ( ASSOCIATED(Solver % Mesh % Edges) ) THEN
           g_beg = maxnode
           l_beg = Solver % Mesh % NumberOfNodes

           n = Solver % Mesh % NumberOfEdges

           edofs = Solver % Mesh % MaxEdgeDOFS
           maxedofs = NINT(ParallelReduction(edofs*1._dp,2))

           maxedge = 0
           DO i=1,n
             maxedge = MAX(maxedge, Solver % Mesh % Edges(i) % GElementindex)
           END DO
           maxedge = NINT(ParallelReduction(1._dp*maxedge,2))

           DO i=1,n
             Element => Solver % Mesh % Edges(i)
             DO j=1,Element % BDOFs
               DO m=1,DOFs
                 l = DOFs*(l_beg + edofs*(i-1)+j-1)+m
                 l=Matrix % Perm(l)
                 IF(l==0) CYCLE
                 Matrix % ParallelInfo % GlobalDOFs(l) = &
                     DOFs*(g_beg+maxedofs*(Element % GelementIndex-1)+j-1)+m
                 Matrix % ParallelInfo % Interface(l) = Mesh % ParallelInfo % EdgeInterface(i)
                
                 ALLOCATE(Matrix % Parallelinfo % NeighbourList(l) % Neighbours(SIZE( &
                      Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours)))
                 Matrix % Parallelinfo % NeighbourList(l) % Neighbours = &
                      Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours
               END DO
             END DO
           END DO
         END IF

         IF ( ASSOCIATED(Solver % Mesh % Faces) ) THEN
           g_beg = maxnode + maxedofs*maxedge
           l_beg = Solver % Mesh % NumberOfNodes + &
                   Solver % Mesh % NumberOfEdges*Solver % Mesh % MaxEdgeDOFs

           n = Solver % Mesh % NumberOfFaces

           fdofs = Solver % Mesh % MaxFaceDOFS
           maxfdofs = NINT(ParallelReduction(fdofs*1._dp,2))

           maxface = 0
           DO i=1,n
             maxface = MAX(maxface, Solver % Mesh % Faces(i) % GElementindex)
           END DO
           maxface = NINT(ParallelReduction(1._dp*maxface,2))

           DO i=1,n
             Element => Solver % Mesh % Faces(i)
             DO j=1,Element % BDOFs
               DO m=1,DOFs
                 l = Matrix % Perm(DOFs*(l_beg + fdofs*(i-1)+j-1)+m)
                 IF(l==0) CYCLE
                 Matrix % ParallelInfo % GlobalDOFs(l) = &
                     DOFs*(g_beg+maxfdofs*(Element % GelementIndex-1)+j-1)+m
                 Matrix % ParallelInfo % Interface(l) = Mesh % ParallelInfo % FaceInterface(i)

                 ALLOCATE(Matrix % Parallelinfo % NeighbourList(l) % Neighbours(SIZE( &
                      Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours)))

                 Matrix % Parallelinfo % NeighbourList(l) % Neighbours = &
                      Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours
               END DO
             END DO
           END DO
         END IF

         IF ( GB ) THEN
           l_beg = Solver % Mesh % NumberOfNodes + &
                   Solver % Mesh % NumberOfEdges*Solver % Mesh % MaxEdgeDOFs + &
                   Solver % Mesh % NumberOfFaces*Solver % Mesh % MaxFaceDOFs

           g_beg = Maxnode +  maxedge*maxedofs + maxface*maxfdofs
           maxbdofs = NINT(ParallelReduction(1._dp*Solver % Mesh % MaxBDOFs,2))

           DO i=1,Solver % Mesh % NumberOfBulkElements
             Element=>Solver % Mesh % Elements(i)
             DO l=1,Element % BDOFs
               DO j=1,DOFs 
                 k = Matrix % Perm(DOFs*(l_beg+Element % BubbleIndexes(l)-1)+j)
                 IF(k==0) CYCLE
                 Matrix % ParallelInfo % GlobalDOFs(k) = &
                   DOFs*(g_beg+maxbdofs*(Element % GElementIndex-1)+l-1)+j
                 Matrix % ParallelInfo % Interface(k) = .FALSE.
                 ALLOCATE(Matrix % ParallelInfo % NeighbourList(k) % Neighbours(1))
                 Matrix % ParallelInfo % NeighbourList(k) % Neighbours=ParEnv % MyPE
               END DO
             END DO
           END DO
         END IF
              WRITE(*,*) 'Processor no ', ParEnv % MyPE,'d'

         ! Add extra degrees of freedom to parallel structures. The additional
         ! variables are assingned to task zero, and are assumed to be shared by
         ! all tasks (TODO: to be optimized if need be...)
         ! --------------------------------------------------------------------
         g_beg = NINT(ParallelReduction(1._dp*MAXVAL(Matrix % ParallelInfo % GlobalDOFs),2))
         pardofs_all = NINT(ParallelReduction(Matrix % ParallelDOFs*1._dp,2))

         Found = ASSOCIATED(Solver % Matrix % ConstraintMatrix)
         IF(Found) Found = ASSOCIATED(Solver % Matrix % ConstraintMatrix % InvPerm)

         j=0
         proc_ind = ParEnv % PEs
         DO i=Matrix % NumberOFRows-Matrix % ExtraDOFs+1,Matrix % NumberOfRows
           j=j+1
           Global_dof = j <= Matrix % ParallelDOFs
           IF (Global_dof) THEN
             Matrix % Parallelinfo % Interface(i) = .TRUE.
             ALLOCATE(Matrix % ParallelInfo % NeighbourList(i) % Neighbours(ParEnv % PEs))

             IF(proc_ind<=0) proc_ind=ParEnv % PEs
             proc_ind=proc_ind-1

!            Matrix % ParallelInfo % NeighbourList(i) % Neighbours(1)=proc_ind
!            l = -1
!            DO k=2,ParEnv % PEs
!              l = l + 1
!              IF(l==proc_ind) l = l + 1
!              Matrix % ParallelInfo % NeighbourList(i) % Neighbours(k)=l
!            END DO
             DO k=1,ParEnv % PEs
               Matrix % ParallelInfo % NeighbourList(i) % Neighbours(k)=ParEnv % PEs-k
             END DO
             k = j + g_beg
           ELSE IF(Found) THEN
             k = Solver % Matrix % ConstraintMatrix % InvPerm(j-Matrix % ParallelDOFs)
             IF(k<=0) CYCLE
             Matrix % Parallelinfo % Interface(i) = Matrix % ParallelInfo % Interface(k)

             l = SIZE(Matrix % ParallelInfo % NeighbourList(k) % Neighbours)
             ALLOCATE(Matrix % ParallelInfo % NeighbourList(i) % Neighbours(l))

             IF (Matrix % Parallelinfo % Interface(i)) &
               CALL Sort(l,Matrix % ParallelInfo % NeighbourList(k) % Neighbours)

             Matrix % Parallelinfo % NeighbourList(i) % Neighbours = &
                Matrix % ParallelInfo % NeighbourList(k) % Neighbours

             k = Matrix % ParallelInfo % GlobalDOFs(k) + g_beg + pardofs_all
           ELSE
             ALLOCATE(Matrix % ParallelInfo % NeighbourList(i) % Neighbours(1))
             Matrix % Parallelinfo % NeighbourList(i) % Neighbours(1)=ParEnv % myPE
             k = j + g_beg + pardofs_all
           END IF
           Matrix % ParallelInfo % GlobalDofs(i) = k
         END DO

         l=Matrix % ExtraDOFs
         n=SIZE(Perm)
         IF(l>0) THEN
            ALLOCATE(Ind(l))
            k = 0
            DO i=Matrix % NumberOfRows-l+1,Matrix % NumberOfRows
              k = k + 1
              Matrix % Perm(n*DOFs+k) = i
              Ind(k) = Matrix % Parallelinfo % Globaldofs(i)
            END DO
            CALL SortI(l,Ind,Matrix % Perm(n*DOFs+1:))
         END IF

       ELSE
         MeshPI => Solver % Mesh % ParallelInfo

         ALLOCATE( Matrix % ParallelInfo )
         MatrixPI => Matrix % ParallelInfo

         n = 0
         DO i=1,Mesh % NumberOfBulkElements
           Element => Mesh % Elements(i)
           IF ( .NOT. ASSOCIATED(Element % DGIndexes) ) CYCLE
           n = MAX(n,MAXVAL(Element % DGIndexes))
         END DO

         ALLOCATE( MatrixPI % GlobalDOFs(n) ); MatrixPI % GlobalDOFs=0

         DO i=1,Mesh % NumberOfBulkElements
           Element => Mesh % Elements(i)
           IF ( .NOT. ASSOCIATED(Element % DGIndexes) ) CYCLE
           DO j=1,SIZE(Element % DGIndexes)
             k = Matrix % Perm(Element % DGIndexes(j))
             IF(K==0) CYCLE
             MatrixPI % GlobalDOFs(k) = 8*(Element % GElementIndex-1)+j
           END DO
         END DO
         ALLOCATE( MatrixPI % Interface(n), MatrixPI % NeighbourList(n) )

         MatrixPI % Interface = .FALSE.
         DO i=1,n
           MtrxN => MatrixPI % NeighbourList(i)
           MtrxN % Neighbours => NULL()
         END DO

         DO i=1,Mesh % NumberOfBulkElements
           Element => Mesh % Elements(i)
           IF ( .NOT. ASSOCIATED(Element % DGIndexes) ) CYCLE
           IF ( ALL(MeshPI % Interface(Element % NodeIndexes)) ) THEN
             DO j=1,Element % Type % NumberOfNodes
                K = Matrix % Perm(Element % DGIndexes(j))
                IF(K==0) CYCLE
                L = Element % Nodeindexes(j)

                MeshN => MeshPI % NeighbourList(L)
                MtrxN => MatrixPI % NeighbourList(K)

                MatrixPI % Interface(k) = .TRUE.

                CALL AllocateVector( MtrxN % Neighbours, &
                      SIZE(MeshN % Neighbours) )
                MtrxN % Neighbours = MeshN % Neighbours
                DO m=1,SIZE(MeshN % Neighbours)
                 IF ( MeshN % Neighbours(m) == Element % PartIndex ) THEN
                   MtrxN % Neighbours(1) = MeshN % Neighbours(m)
                   MtrxN % Neighbours(m) = MeshN % Neighbours(1)
                   EXIT
                 END IF
                END DO
             END DO
           END IF
         END DO
         DO i=1,n
           MtrxN => MatrixPI % NeighbourList(i)
           IF ( .NOT.ASSOCIATED( MtrxN % Neighbours) ) THEN
             CALL AllocateVector(MtrxN % Neighbours,1)
             MtrxN % Neighbours(1) = ParEnv % myPE
           END IF
         END DO
       END IF
              WRITE(*,*) 'Processor no ', ParEnv % MyPE,'e'

       Matrix % ParMatrix => &
          ParInitMatrix( Matrix, Matrix % ParallelInfo )
!if(parenv%mype==0) print*,'MATRIX INIT TIME: ', realtime()-tt
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelInitMatrix2
!-------------------------------------------------------------------------------



END MODULE ErrorEstimationSubs
