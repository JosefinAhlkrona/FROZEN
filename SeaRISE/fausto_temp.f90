FUNCTION fausto_temp( Model, node, zs) RESULT(temp)
    USE Types
    USE CoordinateSystems
    USE SolverUtils
    USe ElementDescription

    !-------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------external variables-----------------------
    TYPE(Model_t) :: Model 
    INTEGER :: node
    TYPE(Variable_t), POINTER :: phiSol, lambdaSol 
    REAL(KIND=dp), POINTER :: phii(:), lambda(:)
    INTEGER, POINTER :: phiPerm(:), lambdaPerm(:)
    REAL(KIND=dp) :: temp, zs
    REAL(KIND=dp) :: pi_180_inv = 180.0_dp/pi
    INTEGER ::  phiDOFs, lambdaDOFs

    REAL(KIND=dp) :: theta_ma, c_ma, kappa_ma, gamma_ma, &
                     theta_mj, c_mj, kappa_mj, gamma_mj

    theta_ma = 41.83_dp
    gamma_ma = -6.309e-03_dp
    c_ma     = -0.7189_dp
    kappa_ma = -0.0672_dp

    theta_mj = 14.70_dp
    gamma_mj = -5.426e-03_dp
    c_mj     = -0.1585_dp
    kappa_mj = -0.0518_dp

    phiSol => VariableGet( Model % Variables, 'Latitude' )
    IF ( ASSOCIATED( phiSol ) ) THEN
       phiPerm    => phiSol % Perm
       phii       => phiSol % Values
       phiDOFs = phiSol % DOFs
    ELSE
       CALL FATAL('fausto_temp', ' Could not find Latitude field variable')
    END IF

    lambdaSol => VariableGet( Model % Variables, 'Longitude' )
    IF ( ASSOCIATED( lambdaSol ) ) THEN
       lambdaPerm    => lambdaSol % Perm
       lambda        => lambdaSol % Values
       lambdaDOFs = lambdaSol % DOFs
    ELSE
       CALL FATAL('fausto_temp', ' Could not find Longitude field variable')
    END IF

    temp = theta_ma &
           + gamma_ma*zs &
           + c_ma*phii(phiDOFs*(phiPerm(node)-1)+1)*pi_180_inv &
           + kappa_ma*(modulo(lambda(lambdaDOFS*(lambdaPerm(node)-1)+1)+pi, 2.0_dp*pi)-pi) &
           *pi_180_inv
           ! west longitudes counted negatively
    temp = temp + 273.15

END FUNCTION fausto_temp

