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

    REAL(KIND=dp),ALLOCATABLE :: ErrorInSIA(:)
    REAL(KIND=dp) :: UNorm

    TYPE(Variable_t), POINTER :: FlowSol, NodeType2Variable
    TYPE(Variable_t), POINTER ::  NodeWiseErrorVariable
    INTEGER, POINTER :: FlowPerm(:), NodeType2Perm(:), NodeWiseErrorPerm(:)
    REAL(KIND=dp), POINTER :: FlowSolution(:), ForceVector(:), &
         NodeType2Values(:), NodeWiseErrorValues(:)
    INTEGER, SAVE :: Timestep 
    TYPE(Variable_t), POINTER :: TimeVar

    TYPE(Matrix_t),POINTER :: A

    REAL(KIND=dp), SAVE :: ErrorBound, ErrorBoundAbs, &
         ErrorBoundRel 
    REAL(KIND=dp),ALLOCATABLE :: Bound(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: DivisionMethod

    real(KIND=dp) :: FlowSolution_velocity, ErrorInSIA_magnitude, Base_Velocity_Cutoff
    integer :: ExtrudeLevels, node_layer, nodes_per_layer, file_unit

    logical :: LOGICAL, Include_z_velocity

    SAVE OnlyHorizontalError, ErrorInSIA, Bound
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
               Bound, &
               STAT=istat )
       END IF
       ALLOCATE( &
            ErrorInSIA(  SIZE( FlowSolution )), &   
            Bound(  Model % Mesh % NumberOfNodes ),&
       STAT=istat)
       AllocationsDone = .TRUE.
    END IF

    NodeWiseErrorVariable => VariableGet( Solver % Mesh % Variables, 'SIAError' )
    IF ( ASSOCIATED(  NodeWiseErrorVariable ) ) THEN
       NodeWiseErrorPerm     => NodeWiseErrorVariable  % Perm
       NodeWiseErrorValues   => NodeWiseErrorVariable % Values
    ELSE
       CALL Fatal( 'FlowSolveSIAFS','Cannot find variable <SIAError>' )
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


    TimeVar => VariableGet( Solver % Mesh % Variables, 'Timestep')
    Timestep = NINT(Timevar % Values(1))

    ErrorInSIA=FlowSolution-SIAVelPermuted

    NodeWiseErrorValues=0.0

    !if error is to big then do FS
    NodeType2=0
    NumberOfFSNodes=0
    NumberOfSIANodes=0





  ! see whether or not the z velocity should be included in the error estimation

   Include_z_velocity = GetLogical( Solver % Values, 'Z velocity Included ', gotIt )
    IF (.NOT. gotIt) THEN
      WRITE( Message, * ) 'Z velocity Included  not found, setting to false by default'
      CALL Info( 'FlowSolve', Message, Level=4 )
      Include_z_velocity = .false.
    END IF

   if(NSDOFs < 4 .and. Include_z_velocity) then
      WRITE( Message, * ) 'Z velocity Included  was set to true, but there are fewer than 4 DOF'
      CALL Info( 'FlowSolve', Message, Level=4 )
      WRITE( Message, * ) 'Setting Z velocity Included  to false'
      CALL Info( 'FlowSolve', Message, Level=4 )
      Include_z_velocity = .false.
   endif

! determine the error for each node
    DO k = 1, Model % Mesh % NumberOfNodes

	  if (NSDOFs == 3 .or. NSDOFs == 4) THEN

		if (Include_z_velocity) THEN
			FlowSolution_velocity = sqrt(FlowSolution(NSDOFs*(FlowPerm(k)-1)+1)**2.0 + &
                                               FlowSolution(NSDOFs*(FlowPerm(k)-1)+2)**2.0 + &
                                               FlowSolution(NSDOFs*(FlowPerm(k)-1)+3)**2.0)

			ErrorInSIA_magnitude = SQRT(ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+1)**2.0 + &
                                              ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+2)**2.0 +  &
                                              ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+3)**2.0  )

		else
			FlowSolution_velocity = sqrt(FlowSolution(NSDOFs*(FlowPerm(k)-1)+1)**2.0 + &
                                               FlowSolution(NSDOFs*(FlowPerm(k)-1)+2)**2.0)



			ErrorInSIA_magnitude = SQRT(ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+1)**2.0 + &
                                              ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+2)**2.0  )
		endif

       	NodeWiseErrorValues(NodeWiseErrorPerm(k)) = ErrorInSIA_magnitude / FlowSolution_velocity

	  else

	    CALL Fatal( 'FlowSolveSIAFS','NODOFs is not 3 or 4' )
	  endif	
		



    END DO


!!!!! SORT NODES
    !! Finding bounds and tolerances

    DivisionMethod  = GetString( Solver % Values, &
         'Approximation Level Determination', gotIt ) 
    IF (.NOT. gotIt) THEN
       WRITE( Message, * ) 'Dividing approximation levels according to tolerance method'
       CALL Info( 'FlowSolve', Message, Level=4 )
       DivisionMethod = 'tolerance'
    END IF

    ErrorBoundAbs = GetConstReal(  Solver % Values,  &
         'Absolute Error Allowed', gotIt )
    IF (.NOT. gotIt) THEN
       WRITE( Message, * ) 'Absolute Error Tolerance not found, setting to 1 m/a'
       CALL Info( 'FlowSolve', Message, Level=4 )
       ErrorBoundAbs = 1.0
    END IF

    SELECT CASE(DivisionMethod)

    CASE('tolerance')

       ErrorBoundRel = 0.01*GetConstReal(  Solver % Values,  &
            'Relative Error Allowed In Percent', gotIt )

       IF (.NOT. gotIt) THEN
          WRITE( Message, * ) 'Relative Error Tolerance not found, setting to 15 %'
          CALL Info( 'FlowSolve', Message, Level=4 )
          ErrorBoundRel = 0.15
       END IF

       DO k = 1, Model % Mesh % NumberOfNodes

          !Weighing absolute and relative error

	   if (NSDOFs == 3 .or. NSDOFs == 4) THEN

		if (Include_z_velocity) THEN

			FlowSolution_velocity = sqrt(FlowSolution(NSDOFs*(FlowPerm(k)-1)+1)**2.0 + &
                                               FlowSolution(NSDOFs*(FlowPerm(k)-1)+2)**2.0 + &
                                               FlowSolution(NSDOFs*(FlowPerm(k)-1)+3)**2.0)


	      else
			FlowSolution_velocity = sqrt(FlowSolution(NSDOFs*(FlowPerm(k)-1)+1)**2.0 + &
                                               FlowSolution(NSDOFs*(FlowPerm(k)-1)+2)**2.0)

		endif

    		Bound(k) = MAXVAL((/ErrorBoundAbs / FlowSolution_velocity, ErrorBoundRel/))


	   else

	     CALL Fatal( 'FlowSolveSIAFS','NODOFs is not 3 or 4' )
	   endif

       END DO
    CASE('sorting')
       Bound = ErrorBoundAbs
    END SELECT

    CALL SortNodes(Model,Solver,NodeType2, NumberOfFSNodes, NumberOfSIANodes, Bound, ErrorInSIA)

    NodeType2Values(NodeType2Perm)=REAL(NodeType2)


	! not sure if this should be done
       DEALLOCATE( &
            ErrorInSIA, &   
            Bound,&
       STAT=istat)
	AllocationsDone = .FALSE.
	write(6,*) "******************************************"
	write(6,*) "SolutionErrorEstimate deallocated:", istat
	write(6,*) "******************************************"	


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

    CHARACTER(LEN=MAX_NAME_LEN) :: FunctionalName,linsysmateth

    CHARACTER(LEN=MAX_NAME_LEN) :: DivisionMethod
    REAL(KIND=dp),ALLOCATABLE :: Bound(:)

!---
 INTEGER, ALLOCATABLE :: Row(:)
       INTEGER :: NVals
       
!---
    TARGET :: x, functional

    SAVE OnlyHorizontalError, &
         NodeWiseError, AT, functional, functionalpointer, &
         pp, ss, xx, &
         Ax, residual, x, FunctionalName, Bound

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
            !   AT % Rows, &
            !   AT % Values, &
            !   AT % Cols, &
               x, &
               functional, &
               Ax, &
               Bound, &
               STAT=istat )
       END IF
       ALLOCATE( &
            NodeWiseError(  Model % Mesh % NumberOfNodes ), &  
           ! AT % Rows(SIZE(A % Rows)), &
           ! AT % Values(SIZE(A % Values)), &
           ! AT % Cols(SIZE(A % Cols)), &
            x(SIZE(A % RHS)), &
            functional(SIZE(A % RHS)), &
            Ax(SIZE(A % RHS)), &            
            Bound(  Model % Mesh % NumberOfNodes ),&
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

    CALL List_toCRSMatrix(AT)


    AT = CRS_Transpose(AT)
    CALL CRS_SortMatrix(AT)

    ALLOCATE(AT % RHS(SIZE(A % RHS)))


    !Allocate the right hand side
    !END SELECT

    functional = 0.0

    functionalpointer => functional

    !Get the functional  
    SELECT CASE(FunctionalName)
    CASE('flux across point') 
       CALL FluxAcrossPoint( Model,Solver,dt,TransientSimulation, &
            functionalpointer)!  1.0_dp 
    CASE('flux across line') 
       CALL FluxAcrossLine( Model,Solver,dt,TransientSimulation, &
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
!stoppa in ISCAL här
       CALL CRS_MatrixVectorMultiply(Solver % Matrix,SIAVelPermuted,Ax)
    END IF

    residual = Ax - Solver % Matrix % RHS

    !-----------------------------------------------------------------------------
    !   Compute x^T*residual (if functional error is too big)
    !-----------------------------------------------------------------------------

    !compute elements of x^T*residual
    DO i = 1, Model % Mesh % NumberOfNodes
       NodeWiseErrorValues(NodeWiseErrorPerm(i))=0.0
       DO k=1, NSDOFs
          NodeWiseErrorValues(NodeWiseErrorPerm(i))=NodeWiseErrorValues(NodeWiseErrorPerm(i))+ &
               ABS(x(NSDOFs*(FlowPerm(i)-1)+k)*residual(NSDOFs*(FlowPerm(i)-1)+k))
       END DO
    END DO

    !-----------------------------------------------------------------------------
    !   Sort Nodes
    !-----------------------------------------------------------------------------
    
    DivisionMethod  = GetString( Solver % Values, &
         'Approximation Level Determination', gotIt ) 
    IF (.NOT. gotIt) THEN
       WRITE( Message, * ) 'Dividing approximation levels according to tolerance method'
       CALL Info( 'FlowSolve', Message, Level=4 )
       DivisionMethod = 'tolerance'
    END IF

    SELECT CASE(DivisionMethod)
    CASE('sorting')
    Bound = 0
    CASE('tolerance')
    ErrorBound = GetConstReal( Solver % Values, 'Nodewise limit for dual problem', gotIt )    
    IF (.NOT. gotIt) THEN
       WRITE( Message, * ) 'Nodewise limit for dual problem not found, setting to 5.0'
       CALL Info( 'Error Estimation: ',Message, Level=4 )
       ErrorBound = 5.0
    END IF
    Bound = ErrorBound
    END SELECT
    NumberOfFSNodes=0
    NumberOfSIANodes=0

  CALL SortNodes(Model,Solver,NodeType2, NumberOfFSNodes, NumberOfSIANodes, Bound)

    NodeType2Values(NodeType2Perm)=REAL(NodeType2)

    DEALLOCATE(AT % RHS)
    DEALLOCATE(AT) 

    ! not sure if this should be done
    AllocationsDone = .FALSE.
    DEALLOCATE(                               &
         NodeWiseError, x, &
         functional, &
         Ax, &
         Bound, &
         STAT=istat )

	write(6,*) "******************************************"
	write(6,*) "FunctionalErrorEstimate deallocated:", istat
	write(6,*) "******************************************"	


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

    CHARACTER(LEN=MAX_NAME_LEN) :: DivisionMethod
    REAL(KIND=dp),ALLOCATABLE :: Bound(:)


    SAVE Ax,residual,NodeWiseResidual, ResidualBound, OnlyHorizontalError, &
         ErrorBound, UserResidualBound, RelationForError, RelativeResidual, &
         Bound

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
               Bound, &
               STAT=istat )
       END IF
       ALLOCATE(  Ax(  SIZE( FlowSolution )), &
            residual(  SIZE( FlowSolution )), &
            NodeWiseResidual(  Model % Mesh % NumberOfNodes ), &
            Bound(Model % Mesh % NumberOfNodes), &
            STAT=istat )    
       IF ( istat /= 0 ) THEN
          CALL Fatal( 'FlowSolve','Memory allocation error, Aborting.' )
       END IF

       AllocationsDone = .TRUE.
    END IF

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
    NodeWiseErrorValues=1.0E25 !Something unreasonable
    RelationForError=0
    ACounter=0

    DO i = 1, GetNOFActive()
       Element => GetActiveElement(i)
       n = GetElementNOFNodes()
  
       DO j=1,n
          k = Element % NodeIndexes(j)

          IF(NodeWiseErrorValues(NodeWiseErrorPerm(k))/=1.0E25) CYCLE !already computed
          !IF (NodeType2(k)/=0) CYCLE
          !
          SELECT CASE( NSDOFs )

          CASE(3) !2D simulation
            NodeWiseErrorValues(NodeWiseErrorPerm(k)) =SQRT(residual(NSDOFs*(FlowPerm(k)-1)+1)**2.0 + &
                  residual(NSDOFs*(FlowPerm(k)-1)+2)**2.0 + &
                  residual(NSDOFs*(FlowPerm(k)-1)+3)**2.0)

          CASE(4) !3D simulation              
            NodeWiseErrorValues(NodeWiseErrorPerm(k)) =SQRT(residual(NSDOFs*(FlowPerm(k)-1)+1)**2.0 + &
                  residual(NSDOFs*(FlowPerm(k)-1)+2)**2.0 + &
                  residual(NSDOFs*(FlowPerm(k)-1)+3)**2.0 + &
                  residual(NSDOFs*(FlowPerm(k)-1)+4)**2.0)
          END SELECT
       END DO
    END DO

 !   NodeType2=0
   NumberOfFSNodes=0
   NumberOfSIANodes=0

 !   DO i = 1, GetNOFActive()
  !     Element => GetActiveElement(i)
   !    n = GetElementNOFNodes()
   !    !     
   !    DO j=1,GetElementNOFNOdes()
    !      k = Element % NodeIndexes(j)
!
 !         IF (NodeType2(k)/=0) CYCLE
  !        !
   !       IF (NodeWiseResidual(k)> ResidualBound) THEN !FS-Node
    !         NumberOfFSNodes=NumberOfFSNodes+1
     !        NodeType2(k) = 2
      !       NodeType2Values(NodeType2Perm(k))=REAL(NodeType2(k))
!
 !         ELSE    !SIA-Node
  !           NumberOfSIANodes=NumberOfSIANodes+1
   !          NodeType2(k) = 1
    !         NodeType2Values(NodeType2Perm(k))=REAL(NodeType2(k))
!
 !         END IF
  !        NodeWiseErrorValues(NodeWiseErrorPerm(k))=NodeWiseResidual(k)
   !    END DO
   ! END DO
   !-----------------------------------------------------------------------------
    !   Sort Nodes
    !-----------------------------------------------------------------------------
    
    DivisionMethod  = GetString( Solver % Values, &
         'Approximation Level Determination', gotIt ) 
    IF (.NOT. gotIt) THEN
       WRITE( Message, * ) 'Dividing approximation levels according to tolerance method'
       CALL Info( 'FlowSolve', Message, Level=4 )
       DivisionMethod = 'tolerance'
    END IF

    SELECT CASE(DivisionMethod)
    CASE('sorting')
       Bound = 0.0
    CASE('tolerance')
       ResidualBound = GetConstReal(  Solver % Values,  &
            'Maximum Allowed Residual', gotIt)
       IF (.NOT. gotIt) THEN
          WRITE( Message, * ) 'Bound for residual not found, setting to 10.0'
          CALL Info( 'Error Estimation: ',Message, Level=4 )
          ResidualBound = 10.0
       END IF
       Bound = ResidualBound
    END SELECT

    CALL SortNodes(Model,Solver,NodeType2, NumberOfFSNodes, NumberOfSIANodes, Bound)

    NodeType2Values(NodeType2Perm)=REAL(NodeType2)

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

    functional = 0.0 !a vector describing the functional  functional(x)=functional^T*x
    functionalpointer => functional

    !Get the functional  
    SELECT CASE(FunctionalName)
    CASE('flux across point') 
       CALL FluxAcrossPoint( Model,Solver,dt,TransientSimulation, &
            functionalpointer)!  1.0_dp
    CASE('flux across line') 
       CALL FluxAcrossLine( Model,Solver,dt,TransientSimulation, &
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

  SUBROUTINE SortNodes(Model,Solver, NodeType2, NumberOfFSNodes, &
       NumberOfSIANodes, Bound, ErrorInSIA) 
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
    !
    !  Sorting nodes into FS and SIA nodes according to ErrorInSIA and Bound,
    !  outputting NumberOfSIAnodes, NumberOfFSNodes, and NodeType2
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
    INTEGER, intent(inout) :: NumberOfSIANodes, NumberOfFSNodes    
    REAL(KIND=dp),ALLOCATABLE, OPTIONAL, intent(in) :: ErrorInSIA(:)
    REAL(KIND=dp),ALLOCATABLE, intent(in) :: Bound(:)


    !------------------------------------------------------------------------------
    !    Local variables
    !------------------------------------------------------------------------------
    TYPE(Element_t),POINTER :: Element
    INTEGER :: ACounter, i, j, k, n, NSDOFs, istat 
    LOGICAL :: OnlyHorizontalError,CantDealWithLoneliness=.FALSE.
    LOGICAL :: AllocationsDone = .FALSE.
    LOGICAL :: gotIt

    TYPE(Variable_t), POINTER :: FlowSol
    TYPE(Variable_t), POINTER ::  NodeWiseErrorVariable
    INTEGER, POINTER :: FlowPerm(:), NodeWiseErrorPerm(:)
    REAL(KIND=dp), POINTER :: FlowSolution(:), NodeWiseErrorValues(:)


    INTEGER, SAVE :: AntalFSnodes,mink,si,counter
    REAL(KIND=dp), SAVE :: FractionFSnodes, MinErrorOfFSNodes, AbsError
    INTEGER, ALLOCATABLE :: FSNodes(:)
    LOGICAL :: BeenSet = .FALSE.
    CHARACTER(LEN=MAX_NAME_LEN) :: DivisionMethod


    integer :: ExtrudeLevels, nodes_per_layer, node_layer
    real(KIND=dp) :: FlowSolution_velocity, Base_Velocity_Cutoff
    logical :: Include_z_velocity

    integer :: minimum_FS_nodes, max_velocity_index
    logical, dimension(:), allocatable :: velocity_mask
    REAL(KIND=dp), dimension(:), allocatable :: velocity_array

    !-----------------------------------------------------------------
    WRITE( Message, * ) '** Sorting Nodes **'
       CALL Info( 'FlowSolve', Message, Level=4 )
    FlowSol => Solver % Variable
    NSDOFs         =  FlowSol % DOFs
    FlowPerm       => FlowSol % Perm
    FlowSolution   => FlowSol % Values

    DivisionMethod  = GetString( Solver % Values, &
         'Approximation Level Determination', gotIt ) 
    IF (.NOT. gotIt) THEN
       DivisionMethod = 'tolerance'
    END IF

    NodeWiseErrorVariable => VariableGet( Solver % Mesh % Variables, 'SIAError' )
    IF ( ASSOCIATED(  NodeWiseErrorVariable ) ) THEN
       NodeWiseErrorPerm     => NodeWiseErrorVariable  % Perm
       NodeWiseErrorValues   => NodeWiseErrorVariable % Values
    ELSE
       CALL Fatal( 'FlowSolveSIAFS','Cannot find variable <SIAError>' )
    END IF



   ! get the number of extruded levels in the mesh

    ExtrudeLevels = GetInteger(Model % Simulation,'Extruded Mesh Levels',gotIt)
    nodes_per_layer = Model % Mesh % NumberOfNodes / ExtrudeLevels


   ! get the basal velocity cutoff. Any node at the base of the ice sheet with a velocity that is 
   ! higher than this will automatically be set to Full Stokes

    Base_Velocity_Cutoff = GetConstReal( Solver % Values, 'Base Velocity Cutoff', gotIt )
    IF (.NOT. gotIt) THEN
      WRITE( Message, * ) 'Base Velocity Cutoff not found, setting to 100 m/yr'
      CALL Info( 'FlowSolve', Message, Level=4 )
      Base_Velocity_Cutoff = 100.
    END IF

    CantDealWithLoneliness = GetLogical( Solver % Values, 'Remove Lonely Pillars', gotIt )
      
  ! see whether or not the z velocity should be included in the error estimation

   Include_z_velocity = GetLogical( Solver % Values, 'Z velocity Included ', gotIt )
    IF (.NOT. gotIt) THEN
      WRITE( Message, * ) 'Z velocity Included  not found, setting to false by default'
      CALL Info( 'FlowSolve', Message, Level=4 )
      Include_z_velocity = .false.
    END IF

   if(NSDOFs < 4 .and. Include_z_velocity) then
      WRITE( Message, * ) 'Z velocity Included  was set to true, but there are fewer than 4 DOF'
      CALL Info( 'FlowSolve', Message, Level=4 )
      WRITE( Message, * ) 'Setting Z velocity Included  to false'
      CALL Info( 'FlowSolve', Message, Level=4 )
      Include_z_velocity = .false.
   endif

! Read in the minimum amount of Full Stokes nodes included in Approximation Levels
   minimum_FS_nodes = GetInteger( Solver % Values, 'Minimum FS Nodes ', gotIt )
   IF (.NOT. gotIt) THEN
     WRITE( Message, * ) 'Minimum FS Nodes  not found, setting to 0 by default'
     CALL Info( 'FlowSolve', Message, Level=4 )
     minimum_FS_nodes = 0
   END IF


    SELECT CASE(DivisionMethod)

    CASE('tolerance')
    WRITE( Message, * ) 'Dividing approximation levels according to tolerance method'
       CALL Info( 'FlowSolve', Message, Level=4 )

       allocate(velocity_array(Model % Mesh % NumberOfNodes), velocity_mask(Model % Mesh % NumberOfNodes))

       velocity_mask = .true.

       DO k = 1, Model % Mesh % NumberOfNodes

		if (Include_z_velocity) THEN
			FlowSolution_velocity = sqrt(FlowSolution(NSDOFs*(FlowPerm(k)-1)+1)**2.0 + &
                                               FlowSolution(NSDOFs*(FlowPerm(k)-1)+2)**2.0 + &
                                               FlowSolution(NSDOFs*(FlowPerm(k)-1)+3)**2.0)


		else
			FlowSolution_velocity = sqrt(FlowSolution(NSDOFs*(FlowPerm(k)-1)+1)**2.0 + &
                                               FlowSolution(NSDOFs*(FlowPerm(k)-1)+2)**2.0)

		endif

	    velocity_array(k) = FlowSolution_velocity

           node_layer = (k-1) / nodes_per_layer + 1


          !Weighing absolute and relative error
	    ! also FS if the base velocity is above the threshold. Needs to be included otherwise things get unstable.
          IF (NodeWiseErrorValues(NodeWiseErrorPerm(k))> Bound(k) .or. &
	      (FlowSolution_velocity > Base_Velocity_Cutoff .and. node_layer == 1)) THEN !FS-Node 
             NumberOfFSNodes=NumberOfFSNodes+1
             NodeType2(k) = 2
		 velocity_mask(k) = .false.
          ELSE    !SIA-Node
             NumberOfSIANodes=NumberOfSIANodes+1
             NodeType2(k) = 1
          END IF

       END DO

	IF (CantDealWithLoneliness) THEN
           CALL RemoveLonelyPillars(Model,Solver, NodeType2, NumberOfFSNodes)
        END IF

	 if(Model % Mesh % NumberOfNodes < minimum_FS_nodes) THEN ! all nodes should be Full Stokes

	   NodeType2 = 2
	   NumberOfFSNodes = Model % Mesh % NumberOfNodes

       else

	   do while (NumberOfFSNodes < minimum_FS_nodes)  ! must include more nodes as full stokes nodes

		max_velocity_index = maxloc(velocity_array,1,velocity_mask)

		NodeType2(max_velocity_index) = 2
		velocity_mask(max_velocity_index) = .false.
            NumberOfFSNodes = NumberOfFSNodes+1

    	   end do

	 endif

       deallocate(velocity_array, velocity_mask)

       NumberOfSIANodes=Model % Mesh % NumberOfNodes-NumberOfFSNodes

    CASE('sorting')
      WRITE( Message, * ) 'Dividing approximation levels according to tolerance method'
       CALL Info( 'FlowSolve', Message, Level=4 )
       FractionFSnodes = GetConstReal(  Solver % Values,  &
            'Number of FS nodes', gotIt )
       IF (.NOT. gotIt) THEN
          WRITE( Message, * ) 'Number of FS nodes not found, setting to 10 %'
          CALL Info( 'FlowSolve', Message, Level=4 )
          FractionFSnodes = 0.1
       END IF
       FractionFSnodes=FractionFSnodes/100.0 !converting from percent
       AntalFSnodes=FLOOR(FractionFSnodes*Model % Mesh % NumberOfNodes)    

       ALLOCATE(FSNodes(AntalFSNodes), STAT=istat)
       FSNodes=0

       DO k = 1, Model % Mesh % NumberOfNodes
          IF (PRESENT(ErrorInSIA)) THEN

		if( NSDOFs == 3) THEN

             	AbsError = SQRT(ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+1)**2.0+ ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+2)**2.0)
		elseif ( NSDOFs == 4) THEN
             	AbsError = SQRT(ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+1)**2.0+ ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+2)**2.0 + &
		                 ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+3)**2.0)
		else
			CALL Fatal( 'FlowSolveSIAFS','NSDOFs not 3 or 4' )
		endif
          ELSE
             AbsError=1.0 !just something that is bigger than zero
          END IF

          IF (NumberOfFSNodes .GE. SIZE(FSNodes)) THEN
             IF (NodeWiseErrorValues(NodeWiseErrorPerm(k))>MinErrorOfFSNodes .AND. AbsError>Bound(k)) THEN     !FS     

                NodeType2(FSNodes(mink))=1 !throwing out the one of the Fs nodew with the smalles error
            
                FSNodes(mink)=k !replacing the one with smallest error with this one
                NodeType2(k)=2  !So this one is a FS node

                !find out who is now closest to be kicked out next time
                MinErrorOfFSNodes=MINVAL(NodeWiseErrorValues(NodeWiseErrorPerm(FSNodes)))  
                DO si=1,size(FSNodes)
                   IF ( NodeWiseErrorValues(NodeWiseErrorPerm(FSNodes(si))).EQ.MinErrorOfFSNodes) THEN
                      mink=si
                   END IF
                END DO
             ELSE
                NodeType2(k)=1
             END IF
          ELSE !didn't fill up all of FSNodes yet

             IF( AbsError > Bound(k) ) THEN !Filling up FS area
                NumberOfFSNodes=NumberOfFSNodes+1
                FSNodes(NumberOfFSNodes)=k !filling up, but only if absolute value is high enough
                NodeType2(k)=2

                IF (NodeWiseErrorValues(NodeWiseErrorPerm(k))<MinErrorOFFSNodes .OR. &
.NOT.BeenSet) THEN
                   MinErrorOfFSNodes=NodeWiseErrorValues(NodeWiseErrorPerm(k))
                   mink=NumberOfFSNodes
                   BeenSet=.TRUE.
                END IF

             ELSE !SIA
                NodeType2(k)=1
             END IF

             !only if we're filling up the last spot
             IF (NumberOfFSNodes .EQ. SIZE(FSNodes)) THEN !This is the last time we fill up
                IF (.NOT.BeenSet) THEN
                   mink=1
                   MinErrorOfFSNodes=0.0
                END IF
             END IF
          END IF

       END DO
       
       IF (CantDealWithLoneliness) THEN
       	   CALL RemoveLonelyPillars(Model,Solver, NodeType2, NumberOfFSNodes)
       END IF

       NumberOfSIANodes=Model % Mesh % NumberOfNodes-NumberOfFSNodes
       DEALLOCATE(FSNodes)
    END SELECT


  END SUBROUTINE SortNodes

SUBROUTINE RemoveLonelyPillars(Model,Solver, NodeType2, NumberOfFSNodes)
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
    !
    !  Removing some lonely SIA or FS nodes to make areas less dotty
    !  
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
    INTEGER, intent(inout) :: NumberOfFSNodes
    INTEGER :: node, numberofneighbours, neighbours(20),i,maxnofneighbour

  INTEGER :: processor_number, file_unit, file_number_length
  character(len=80) :: file_in, format1
  character(len=20) :: file_number
    !-----------------------------------------------------------------
    !-----------------------------------------------------------------

    !Retrieve information about which nodes are neighbours
    !if alone, set nodetype to same as others	

  processor_number = ParEnv % MyPE +1

  file_unit = 1000 + processor_number

  write(file_number,*) processor_number
  file_number_length = log10(real(processor_number)) + 1
  write(format1,*)  '(A10,I', file_number_length, ",A4)"

  write(file_in,format1)  "neighbours", processor_number, ".txt"


  open (unit=file_unit, file=file_in, access="sequential", form="formatted", status="old")
    READ(file_unit,*) maxnofneighbour
    DO i = 1, Model % Mesh % NumberOfNodes !find neighbours for all nodes
	    READ(file_unit,"(I7.2,I7.2)",advance='no')  node,numberofneighbours
	    READ(file_unit,*) neighbours(1:numberofneighbours)
	  !  WRITE(*,*) i,node,numberofneighbours,neighbours(1:numberofneighbours)
	 !   WRITE(*,*) NodeType2(neighbours(1:numberofneighbours))
	    IF ( ALL( NodeType2(neighbours(1:numberofneighbours)).NE. NodeType2(i)) ) THEN
!		write(558,*) i, node, numberofneighbours,neighbours(1:numberofneighbours), "| ", NodeType2(i), "| ",&
!		 NodeType2(neighbours(1:numberofneighbours))

		NodeType2(i)=NodeType2(neighbours(1))
	    END IF
    END DO

    close(file_unit)

 
END SUBROUTINE RemoveLonelyPillars 





END MODULE ErrorEstimationSubs
