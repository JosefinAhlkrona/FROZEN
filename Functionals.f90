MODULE Functionals

CONTAINS

  SUBROUTINE FluxAcrossPoint ( Model,Solver,dt,TransientSimulation, &
       a  )
    !******************************************************************************
    !
    !  Write a function for computing the horizontal flux over a vertical line,
    !  in the form of a vector, a. a = integral of vx over the line. The 
    !  trapezoidal method is used to compute the integral. 
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
    !******************************************************************************
    !------------------------------------------------------------------------------


    USE DefUtils
    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Solver_t), TARGET :: Solver

    REAL(KIND=dp) :: dt
    LOGICAL :: TransientSimulation

    TYPE(Variable_t), POINTER :: FlowSol 
    INTEGER, POINTER :: FlowPerm(:) 
    REAL(KIND=dp), POINTER :: FlowSolution(:) 

    REAL(KIND=dp) :: xcoord
    REAL(KIND=dp), POINTER :: a(:)
    INTEGER :: i, mini
    REAL(KIND=dp) :: dist, olddist

    !-----------Variables needed for integration --------------------------------
    TYPE(Solver_t), POINTER :: PSolver
    INTEGER :: j,k,l,Dofs,dof,nsize,TopNodes,BotNodes
    INTEGER, POINTER :: TopPointer(:),BotPointer(:),UpPointer(:),DownPointer(:)
    LOGICAL :: Initialized = .FALSE.,GotVar
    REAL(KIND=dp) :: dx,Level,q
    REAL(KIND=dp), POINTER :: Coord(:)
    TYPE(Variable_t), POINTER :: Var, OldVar

    SAVE :: BotPointer, TopPointer, UpPointer, DownPointer, Coord, TopNodes, &
         BotNodes

    !------------------------------------------------------------------------------
    !    Get the solution  vx,vy,p
    !------------------------------------------------------------------------------
    FlowSol => Solver % Variable
    FlowPerm       => FlowSol % Perm
    FlowSolution   => FlowSol % Values


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
    !    Figure out what line you wanna integrate
    !------------------------------------------------------------------------------
    !get the coordinate  
    xcoord = GetConstReal( Solver % Values, 'Point x-coord', GotVar)    
    IF (.NOT. GotVar) THEN
       CALL FATAL( 'Error Estimation: ','Point x-coord not set')
    END IF


    olddist = 10e10 !big number
    !figure out which line is closest
    DO i=1,Model % Mesh % NumberOfNodes
       dist = ABS(Model % Nodes % x(i) - xcoord)
       IF (dist < olddist) THEN
          mini=i
          olddist=dist
       END IF
    END DO

    WRITE( Message, * ) 'Computing flux over x = ', Model % Nodes % x(mini) 
    CALL Info( 'Error Estimation: ',Message, Level=4 )
    !------------------------------------------------------------------------------
    !    Compute a
    !------------------------------------------------------------------------------

    a=0.0      

    !Get the bottom pointer for the node you found above
    i = BotPointer(mini)

    dx = (Coord(UpPointer(i)) - Coord(i))

    a(3*(FlowPerm(i)-1)+1)= 0.5*dx

    DO WHILE (i /= TopPointer(i))
       i = UpPointer(i)
       dx = (Coord(UpPointer(i)) - Coord(DownPointer(i)))
       a(3*(FlowPerm(i)-1)+1)=0.5*dx
    END DO

    dx = (Coord(i) - Coord(DownPointer(i)))
    a(3*(FlowPerm(i)-1)+1)= 0.5*dx

  END SUBROUTINE FluxAcrossPoint


END MODULE Functionals
