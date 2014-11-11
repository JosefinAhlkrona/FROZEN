FUNCTION linkZsGreen_newData (Model, nodenumber, Array) RESULT(TopSurface)
  USE DefUtils
  IMPLICIT NONE
  ! in
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL (KIND=dp) :: TopSurface, Array(4), SurfaceElevation, BedrockElevation
  ! internal
  REAL (KIND=dp) :: Time, FS, X, Y, Z, R, L, H0, N, EPS, a, pr, hpr, lim
  INTEGER :: DIM
  LOGICAL :: FirstTime = .TRUE.

  SAVE FirstTime, DIM


  DIM = CoordinateSystemDimension()
  EPS = 10.0

  !PRINT *, "Node=", nodenumber, "A=", Array(1:4)
  Time = Array(1)
  IF (Time == 1) THEN 
     FirstTime = .TRUE.
  ELSE
     FirstTime = .FALSE.
  END IF
  SurfaceElevation = Array(2)
  BedrockElevation = Array(3)
  IF (SurfaceElevation-BedrockElevation < EPS) THEN
        SurfaceElevation = BedrockElevation+ EPS
  END IF
  FS=Array(4)



  IF (FirstTime) THEN
        TopSurface= SurfaceElevation
        FirstTime = .FALSE.
      !WRITE(*,*) 'FS=', FS
  ELSE
     TopSurface = FS
   IF (FS-BedrockElevation < EPS) THEN
        TopSurface = BedrockElevation+ EPS
  END IF
  END IF
  RETURN
END FUNCTION linkZsGreen_newData
