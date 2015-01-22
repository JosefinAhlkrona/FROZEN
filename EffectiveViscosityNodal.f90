 FUNCTION EffectiveViscosityNodal( Viscosity,Density,Ux,Uy,Uz,Element, &
        Nodes,n,nd,u,v,w, muder ) RESULT(mu)

      USE ModelDescription
      USE MaterialModels

   REAL(KIND=dp)  :: Viscosity,Density,u,v,w,mu,Ux(:),Uy(:),Uz(:)
     REAL(KIND=dp), OPTIONAL :: muder
     TYPE(Nodes_t)  :: Nodes
     INTEGER :: n,nd
     TYPE(Element_t),POINTER :: Element

     !------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3)
     REAL(KIND=dp) :: ss,s,SqrtMetric,SqrtElementMetric,Velo(3)
     REAL(KIND=dp) :: Metric(3,3), dVelodx(3,3), CtrMetric(3,3), &
          Symb(3,3,3), dSymb(3,3,3,3)

     INTEGER :: i,j,k
     LOGICAL :: stat,GotIt

     CHARACTER(LEN=MAX_NAME_LEN) :: ViscosityFlag, TemperatureName
     TYPE(ValueList_t), POINTER :: Material
     REAL(KIND=dp) :: x, y, z, c1n(n), c2n(n), c3n(n), c4n(n), &
          c1, c2, c3, c4, c5, c6, c7, Temp, NodalTemperature(n), Tlimit, TempCoeff, &
          h, A1, A2, Q1, Q2, R, NodalEhF(n), EhF, ArrheniusFactor

     ! Temperature is needed for thermal models
     TYPE(Variable_t), POINTER :: TempSol 
     REAL(KIND=dp), POINTER :: Temperature(:)
     INTEGER, POINTER :: TempPerm(:)

     INTEGER(KIND=AddrInt) :: Fnc

     TYPE(Variable_t), POINTER :: Var

     REAL(KIND=dp) :: dist,F2,F3
     REAL(KIND=dp) :: KE_K, KE_E, KE_Z, CT, TimeScale,Clip, Cmu, Vals(n)

     CHARACTER(LEN=MAX_NAME_LEN) :: str

     LOGICAL :: SetArrheniusFactor=.FALSE.
     REAL(KIND=dp),ALLOCATABLE :: viskositet(:)
     mu = Viscosity

    ALLOCATE(viskositet( CurrentModel % Mesh % NumberOfNodes), STAT=istat )

  IF ( PRESENT(muder) ) muder=0

     k = ListGetInteger( CurrentModel % Bodies(Element % BodyId) % Values, 'Material', &
          minv=1, maxv=CurrentModel % NumberOFMaterials )

     Material => CurrentModel % Materials(k) % Values

     ViscosityFlag = ListGetString( Material,'Viscosity Model', GotIt)

     IF(.NOT. gotIt) RETURN

 !------------------------------------------------------------------------------
     !    Basis function values & derivatives at the calculation point
     !------------------------------------------------------------------------------
     stat = ElementInfo( Element,Nodes,u,v,w, &
          SqrtElementMetric, Basis,dBasisdx )


      DO i=1,n
        IF viskositet(Element % NodeIndexes(i)).NE.0.0 CYCLE !Already set

     !------------------------------------------------------------------------------
     !   Coordinate system dependent information
     !------------------------------------------------------------------------------
     x =  Nodes % x(i)
     y =  Nodes % y(i)
     z =  Nodes % z(i)

     CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
     !------------------------------------------------------------------------------
     DO j=1,3
        dVelodx(1,j) =Ux(i)*dBasisdx(i,j) 
        dVelodx(2,j) =Uy(i)*dBasisdx(i,j) 
        dVelodx(3,j) =Uz(i)*dBasisdx(i,j) 
     END DO

     Velo(1) =  Basis(i) * Ux(i) 
     Velo(2) =  Basis(i) * Uy(i) 
     Velo(3) =  Basis(i) * Uz(i) 

     ! This is the square of shearrate which results to 1/2 in exponent 
     ! Also the derivative is taken with respect to the square
     !-------------------------------------------------------------------
     ss = SecondInvariant(Velo,dVelodx,Metric,Symb)/2

     SELECT CASE( ViscosityFlag )

CASE('power law')

        c2n = ListGetReal( Material, 'Viscosity Exponent', 1, Element % NodeIndexes(i) )

  
        s = ss
        IF (PRESENT(muder)) muder = Viscosity * (c2-1)/2 * s**((c2-1)/2-1)

        c3n = ListGetReal( Material, 'Critical Shear Rate',1, Element % NodeIndexes(i),gotIt )
        IF (GotIt) THEN
           c3 = c3n
           IF(s < c3**2) THEN
              s = c3**2
              IF (PRESENT(muder)) muder = 0._dp
           END IF
        END IF
        mu = Viscosity * s**((c2-1)/2)

END SELECT

        viskositet(Element % NodeIndexes(i))=mu

END DO
WRITE(*,*) '---'
    WRITE(*,*) 'viskositet=',viskositet
WRITE(*,*) '---'
    DEALLOCATE(viskositet)

END FUNCTION EffectiveViscosityNodal
