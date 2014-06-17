!------------------------------------------------------------------------------
  RECURSIVE FUNCTION ListGetRealAtNode( List, Name, Node, Found ) RESULT(s)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*)  :: Name
     INTEGER :: Node
     LOGICAL, OPTIONAL :: Found
     REAL(KIND=dp) :: s
!-----------------------------------------------------------------------------
     INTEGER, TARGET, SAVE :: Dnodes(1)
     INTEGER, POINTER :: NodeIndexes(:) 
     REAL(KIND=dp) :: x(1)
     INTEGER, PARAMETER :: one = 1

     IF ( PRESENT( Found ) ) Found = .FALSE.

     IF ( ASSOCIATED(List) ) THEN
       NodeIndexes => Dnodes
       NodeIndexes(one) = Node
       
       x(1:one) = ListGetReal( List, Name, one, NodeIndexes, Found )
       s = x(one)
     ELSE
       s = 0.0_dp
     END IF

!------------------------------------------------------------------------------
