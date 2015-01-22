
  !Integration from bottom to z-level, trapezoidal method
  DO i=1,nsize
     IF( i == BotPointer(i) ) THEN !so this loop will be over the extruded lines
        l = i !start at the bottom of the line
        A3hminusz3(l) = 0.0_dp !Integral is zero at the bottom
        DO k=1,nsize

           !              Surf = FreeSurf(FreeSurfPerm(TopPointer(l)))
           IF (dim == 2) THEN
              Surf= Model % Nodes % y(TopPointer(l))
              Position =  Model % Nodes % y(l)
              PositionUnder=Model % Nodes % y(DownPointer(l))
           ELSEIF (dim ==3) THEN
              Surf= Model % Nodes % z(TopPointer(l))
              Position = Model % Nodes % z(l)
              PositionUnder= Model % Nodes % z(DownPointer(l))
           END IF

           IF( k > 1 ) THEN !above bottom

              dx = (Coord(l) - Coord(DownPointer(l)))

              A3hminusz3(l)=A3hminusz3(DownPointer(l))+0.5_dp*dx*( &
                   2*ArrheniusFactor(l)*(rho*g)**nGlen*& 
                   (Surf - Position)**nGlen + &
                   2*ArrheniusFactor(DownPointer(l))*(rho*g)**nGlen*&
                   (Surf - PositionUnder)**nGlen)
           END IF

           IF( l == TopPointer(l)) EXIT !if reached end of extruded line, go to next one
           l = UpPointer(l)   !step up    

        END DO
     END IF
  END DO
