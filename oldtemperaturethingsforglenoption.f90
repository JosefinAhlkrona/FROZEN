ConstTemp



SetArrheniusFactor = GetLogical(Material, 'Set Arrhenius Factor', GotIt)

     IF (SetArrheniusFactor) THEN
        ConstTemp= .TRUE.
     END IF

     IF ( (.NOT.GotIt) .OR. .NOT.(SetArrheniusFactor) ) THEN

        Temp = GetConstReal(Model % Materials(1) % Values, "Constant Temperature", ConstTemp)

        IF(.NOT.ConstTemp) THEN !we have to find a temperature field

           TemperatureName = GetString(Material, 'Temperature Field Variable', GotIt)
           IF (.NOT.GotIt) WRITE(TemperatureName,'(A)') 'Temperature'

           TempSol => VariableGet( CurrentModel % Variables,TRIM(TemperatureName))
 	   IF (ASSOCIATED(TempSol)) THEN
		   TempPerm    => TempSol % Perm
		   Temperature => TempSol % Values  
           ELSE
     	     CALL FATAL(SolverName,'Could not find temperature variable')
           END IF

        END IF

        R = GetConstReal( CurrentModel % Constants,'Gas Constant',GotIt)
        IF (.NOT.GotIt) R = 8.314_dp
        ! lets for the time being have this hardcoded
        Tlimit = GetConstReal(Material, 'Limit Temperature', GotIt)
        IF (.NOT.GotIt) THEN
           Tlimit = -10.0_dp
           CALL INFO('SIASolver','Limit Temperature not found. Setting to -10', Level=5)
        END IF
        A1 = GetConstReal(Material, 'Rate Factor 1', GotIt)
        IF (.NOT.GotIt) THEN
           A1 = 3.985d-13
           CALL INFO('SIASolver','Rate Factor 1 not found. Setting to 3.985e-13', Level=5)
        END IF
        A2 = GetConstReal(Material, 'Rate Factor 2', GotIt)
        IF (.NOT.GotIt) THEN
           A2 = 1.916d03
           CALL INFO('SIASolver','Rate Factor 2 not found. Setting to 1.916E03', Level=5)
        END IF
        Q1 = GetConstReal(Material, 'Activation Energy 1', GotIt)
        IF (.NOT.GotIt) THEN
           Q1 = 60.0d03
           CALL INFO('SIASolver','Activation Energy 1 not found. Setting to 60.0E03', Level=5)
        END IF
        Q2 = GetConstReal(Material, 'Activation Energy 2', GotIt)
        IF (.NOT.GotIt) THEN
           Q2 = 139.0d03
           CALL INFO('SIASolver','Activation Energy 2 not found. Setting to 139.0d03', Level=5)
        END IF

        IF (ConstTemp) THEN
           IF (Temp.LE. Tlimit) THEN
              ArrheniusFactor = A1 * EXP( -Q1/(R * (273.15 + Temp)))
              dArrheniusFactordT = ArrheniusFactor*R/Q1 !eh how did i think here?
           ELSE IF((Tlimit<Temp) .AND. (Temp .LE. 0.0_dp)) THEN
              ArrheniusFactor = A2 * EXP( -Q2/(R * (273.15 + Temp)))
              dArrheniusFactordT = ArrheniusFactor*R/Q2

           ELSE
              ArrheniusFactor = A2 * EXP( -Q2/(R * (273.15)))
              dArrheniusFactordT = ArrheniusFactor*R/Q2
              CALL INFO('SIASolver','Positive Temperature detected in Glen - limiting to zero!', Level = 5)
           END IF
        ELSE

           DO i = 1, Model % Mesh % NumberOfNodes
              Temp=Temperature(TempPerm(i))
              IF (Temp.LE. Tlimit) THEN
                 ArrheniusFactor(i) = A1 * EXP( -Q1/(R * (273.15 + Temp)))
                 dArrheniusFactordT(i) = ArrheniusFactor(i)*R/Q1
              ELSE IF((Tlimit<Temp) .AND. (Temp .LE. 0.0_dp)) THEN
                 ArrheniusFactor(i) = A2 * EXP( -Q2/(R * (273.15 + Temp)))
                 dArrheniusFactordT(i) = ArrheniusFactor(i)*R/Q2
              ELSE
                 ArrheniusFactor(i) = A2 * EXP( -Q2/(R * (273.15)))
                 dArrheniusFactordT(i) = ArrheniusFactor(i)*R/Q2
                 CALL INFO('SIASolver','Positive Temperature detected in Glen - limiting to zero!', Level = 5)
              END IF
           END DO

        END IF





----




 TempGradSol => VariableGet( Solver % Mesh % Variables, 'Temperature grad')
        IF (ASSOCIATED(TempGradSol)) THEN
           TempGrad => TempGradSol % Values
           TempGradPerm => TempGradSol % Perm
        ELSE
           CALL FATAL(SolverName,'Could not find variable >Temperature grad<')
        END IF

