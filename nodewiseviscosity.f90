



 DO i = 1, GetNOFActive()
           Element => GetActiveElement(i)
           n = GetElementNOFNodes()

           SiaNodes(1:n) = GetReal(GetMaterial(), 'SIA node',GotIt)


           DO j=1,GetElementNOFNOdes()
              k = Element % NodeIndexes(j)
              IF(NodeType2(k)/=0) CYCLE !already sorted

              IF(SiaNodes(j)>0) THEN
                 NumberOfSIANodes=NumberOfSIANodes+1
                 NodeType2(k) = 1
                 NodeType2Values(NodeType2Perm(k))=REAL(NodeType2(k))
              ELSE
                 NumberOfFSNodes=NumberOfFSNodes+1
                 NodeType2(k) = 2
                 NodeType2Values(NodeType2Perm(k))=REAL(NodeType2(k))
              END IF
           END DO
        END DO !nofactiveelemes
