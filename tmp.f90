
  SUBROUTINE FluxAcrossLine ( Model,Solver,dt,TransientSimulation, &
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

    CHARACTER(LEN=MAX_NAME_LEN) ::
    INTEGER :: linelength
    INTEGER, ALLOCATABLE :: linenode
    REAL(KIND=dp), POINTER :: a(:)

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
    !    Figure out what line you wanna integrate over
    !------------------------------------------------------------------------------
    !get the description of the line
    filename = GetString( Solver % Values, 'Line Description', GotVar)    
    IF (.NOT. GotVar) THEN
       CALL FATAL( 'Error Estimation: ','File describing line not found')
    END IF

    OPEN(unit = 102, file = FileName, status = 'old', action = 'read')
    read(102,*) linelength

    ALLOCATE( linenode(linelength), STAT=istat )

    IF ( istat /= 0 ) THEN
       CALL Fatal( 'Error Estimation','Memory allocation error, Aborting.' )
    END IF

    DO i=1,linelength
       READ(102,*) linenode(i)
    END DO

    CLOSE(102)
    
    !------------------------------------------------------------------------------
    !    Compute a
    !------------------------------------------------------------------------------

    a=0.0      

    DO k=1,linelength !looping over the nodes in the line, computing an integral from bottom to top
    
    i = BotPointer(linenode(k)) !bottom node
!-----------
    dx = (Coord(UpPointer(i)) - Coord(i)) 

    a(3*(FlowPerm(i)-1)+1)= 0.5*dx

    DO WHILE (i /= TopPointer(i))
       i = UpPointer(i)
       dx = (Coord(UpPointer(i)) - Coord(DownPointer(i)))
       a(3*(FlowPerm(i)-1)+1)=0.5*dx
    END DO

    dx = (Coord(i) - Coord(DownPointer(i)))
    a(3*(FlowPerm(i)-1)+1)= 0.5*dx
!----------
END DO

    DEALLOCATE(linenode)

  END SUBROUTINE FluxAcrossLine

