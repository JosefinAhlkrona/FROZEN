MODULE ErrorEstimationSubs

CONTAINS

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



  SUBROUTINE SolutionErrorEstimate( Model,Solver,dt,TransientSimulation, NodeType2, &
       SIAVelPermuted, NumberOfSIANodes, NumberOfFSNodes, ReorderTimeInterval)
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
    INTEGER, intent(in) :: ReorderTimeInterval

    !------------------------------------------------------------------------------
    !    Local variables
    !------------------------------------------------------------------------------
    TYPE(Element_t),POINTER :: Element
    REAL(KIND=dp) :: ErrorBound   
    INTEGER :: ACounter, i, j, k, n, NSDOFs, istat 
    LOGICAL :: OnlyHorizontalError
    LOGICAL :: AllocationsDone = .FALSE.
    LOGICAL :: gotIt

    REAL(KIND=dp),ALLOCATABLE :: CoupledSolution(:), ErrorInSIA(:), &
         NodeWiseError(:)
    REAL(KIND=dp) :: UNorm

    TYPE(Variable_t), POINTER :: FlowSol 
    INTEGER, POINTER :: FlowPerm(:) 
    REAL(KIND=dp), POINTER :: FlowSolution(:), ForceVector(:) 
    INTEGER, SAVE :: Timestep 
    TYPE(Variable_t), POINTER :: TimeVar

    TYPE(Matrix_t),POINTER :: A


    SAVE OnlyHorizontalError, ErrorBound, CoupledSolution, &
         ErrorInSIA, NodeWiseError



!!!---------------------------------------------------

    WRITE( Message, * ) 'Computing the Error'
    CALL Info( 'FlowSolve', Message, Level=4 )

    ErrorBound = 0.01*GetConstReal(  Solver % Values,  &
         'Relative Error Allowed In Percent', gotIt )
    OnlyHorizontalError = GetLogical( Solver % Values, &
         'Only Compute Horizontal Error', gotIt ) 

    TimeVar => VariableGet( Solver % Mesh % Variables, 'Timestep')
    Timestep = NINT(Timevar % Values(1))

    FlowSol => Solver % Variable
    NSDOFs         =  FlowSol % DOFs
    FlowPerm       => FlowSol % Perm
    FlowSolution   => FlowSol % Values
    A => Solver % Matrix
    ForceVector => A % RHS
    UNorm = Solver % Variable % Norm

    !------------------------------------------------------------------------------
    !     Allocate some permanent storage, this is done first time only
    !------------------------------------------------------------------------------

    IF ( .NOT.AllocationsDone .OR. Solver % Mesh % Changed ) THEN

       N = Solver % Mesh % MaxElementDOFs

       IF( AllocationsDone ) THEN
          DEALLOCATE( &
               CoupledSolution, &
               ErrorInSIA, &
               NodeWiseError, &
               STAT=istat )
       END IF
       ALLOCATE(&
            CoupledSolution(  SIZE( FlowSolution )), &
            ErrorInSIA(  SIZE( FlowSolution )), &            
            NodeWiseError(  Model % Mesh % NumberOfNodes ), &
            STAT=istat )    
       IF ( istat /= 0 ) THEN
          CALL Fatal( 'FlowSolve','Memory allocation error, Aborting.' )
       END IF

       AllocationsDone = .TRUE.
    END IF  !------------------------------------------------------------------------------


    CoupledSolution=FlowSolution

    !Solve System

    Solver % Variable  % Values=0._dp

    !A => Solver % Matrix
    !ForceVector => A % RHS
    !FlowSol => Solver % Variable
    !FlowSolution => FlowSol % Values

    UNorm = DefaultSolve()


    ErrorInSIA=FlowSolution-SIAVelPermuted

    NodeWiseError=0.0_dp
WRITE(*,*) '************IN ERRORSUBS*********'

    DO i = 1, GetNOFActive()
       Element => GetActiveElement(i)
       n = GetElementNOFNodes()
       !     
       DO j=1,GetElementNOFNOdes()
          k = Element % NodeIndexes(j)
          SELECT CASE( NSDOFs )
          CASE(3)
             IF (FlowSolution(NSDOFs*(FlowPerm(k)-1)+1)/=0) THEN

                NodeWiseError(k)=ABS(ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+1)&
                     /FlowSolution(NSDOFs*(FlowPerm(k)-1)+1))


             ELSE
                ! WRITE(*,*) 'Error would be NaN'
                NodeWiseError(k)=0.0_dp
             END IF
          CASE(4)

             WRITE(*,*) 'vx=', FlowSolution(NSDOFs*(FlowPerm(k)-1)+1), 'for node ', k 
             WRITE(*,*) 'vy=', FlowSolution(NSDOFs*(FlowPerm(k)-1)+2), 'for node ', k
             WRITE(*,*) 'siavx=', SIAVelPermuted(NSDOFs*(FlowPerm(k)-1)+1), 'for node ', k
             WRITE(*,*) 'siavy=', SIAVelPermuted(NSDOFs*(FlowPerm(k)-1)+2), 'for node ', k 
             WRITE(*,*) 'coupledvx=',  CoupledSolution(NSDOFs*(FlowPerm(k)-1)+1), 'for node ', k
             WRITE(*,*) 'coupledvy=', CoupledSolution(NSDOFs*(FlowPerm(k)-1)+2), 'for node ', k 

             IF (FlowSolution(NSDOFs*(FlowPerm(k)-1)+1)/=0 .AND. FlowSolution(NSDOFs*(FlowPerm(k)-1)+2)/=0) THEN
                NodeWiseError(k)=ABS(ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+1)&
                     /FlowSolution(NSDOFs*(FlowPerm(k)-1)+1)+ &
                     ErrorInSIA(NSDOFs*(FlowPerm(k)-1)+2)&
                     /FlowSolution(NSDOFs*(FlowPerm(k)-1)+2))  
             ELSE
                ! WRITE(*,*) 'Error would be NaN'
                NodeWiseError(k)=0.0_dp
             END IF
             WRITE(*,*) 'NodeWiseError(k)=',NodeWiseError(k)
          END SELECT
       END DO
    END DO

WRITE(*,*) '***************************************'

    FlowSolution=CoupledSolution !necessary???


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
          !
          IF (NodeWiseError(k)> ErrorBound) THEN !FS-Node
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

    WRITE(*,*) 'Error estimation, number of FS nodes: ', NumberOfFSNodes, ' and number of SIA nodes: ', NumberOfSIANodes 

    WRITE(*,*) 'Printing Solution n Error to Solution.dat and Error.dat'

    IF ( Timestep==ReorderTimeInterval) THEN !First time residual is computed
       open (unit=91, file="Solution.dat",STATUS='NEW',POSITION='APPEND')
       open (unit=97, file="Error.dat",STATUS='NEW',POSITION='APPEND')
       !Compute the difference between SIA solution and FS solution in FS-nodes
    ELSE
       open (unit=91, file="Solution.dat",STATUS='OLD',POSITION='APPEND')
       open (unit=97, file="Error.dat",STATUS='OLD',POSITION='APPEND')
    END IF

    WRITE(91,*) 'Timestep=',Timestep
    WRITE(97,*) 'Timestep=',Timestep
    DO i = 1, SIZE(FlowSolution)/NSDOFs
       WRITE(91,*) i,FlowSolution(NSDOFs*(FlowPerm(i)-1)+1)
       WRITE(97,*) i,NodeWiseError(i)
    END DO

    close (91)
    close (97)
  END SUBROUTINE SolutionErrorEstimate

END MODULE ErrorEstimationSubs
