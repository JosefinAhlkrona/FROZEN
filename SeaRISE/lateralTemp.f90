FUNCTION lateralTemp(Model, Node, z) RESULT(temp)

USE Types  
USE CoordinateSystems  
USE SolverUtils  
USE ElementDescription

IMPLICIT NONE

TYPE(Model_t) :: Model
INTEGER :: Node, zsDOFs, phiDOFs, lambdaDOFs

REAL(KIND=dp) :: z
REAL(KIND=dp) :: temp
REAL(KIND=dp) :: theta_ma, c_ma, kappa_ma, gamma_ma, &
                  theta_mj, c_mj, kappa_mj, gamma_mj

INTEGER, POINTER :: phiPerm(:), lambdaPerm(:)
REAL(KIND=dp), POINTER :: phii(:), lambda(:) 
TYPE(Variable_t), POINTER :: phiSol, lambdaSol

REAL(KIND=dp), PARAMETER :: pi_180 = pi/180.0_dp
REAL(KIND=dp), PARAMETER :: pi_180_inv = 180.0_dp/pi

CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

SolverName = 'lateralTemp' 

phiSol => VariableGet( Model % Variables, 'Latitude' )
IF ( ASSOCIATED( phiSol ) ) THEN
      phiPerm    => phiSol % Perm
      phii       => phiSol % Values
      phiDOFs = phiSol % DOFs
ELSE
      CALL FATAL(SolverName, ' Could not find Latitude field variable')
END IF

lambdaSol => VariableGet( Model % Variables, 'Longitude' )
IF ( ASSOCIATED( lambdaSol ) ) THEN
     lambdaPerm    => lambdaSol % Perm
     lambda        => lambdaSol % Values
     lambdaDOFs = lambdaSol % DOFs
ELSE
     CALL FATAL(SolverName, ' Could not find Longitude field variable')
END IF

! Hard coded Parameterization by Fausto et al. (2009). Should be consistent with what is used in pdd

theta_ma = 41.83_dp
gamma_ma = -6.309e-03_dp
c_ma     = -0.7189_dp
kappa_ma = -0.0672_dp

theta_mj = 14.70_dp
gamma_mj = -5.426e-03_dp
c_mj     = -0.1585_dp
kappa_mj = -0.0518_dp

temp = theta_ma &
       + gamma_ma*z &
       + c_ma*phii(phiDOFs*(phiPerm(Node)-1)+1)*pi_180_inv &
       + kappa_ma*(modulo(lambda(lambdaDOFS*(lambdaPerm(Node)-1)+1)+pi, 2.0_dp*pi)-pi) &
                 *pi_180_inv

temp = temp+273.15_dp

END FUNCTION lateralTemp
