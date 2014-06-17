FUNCTION linkZsGreen_newData (Model, nodenumber, Array) RESULT(Slide)
  USE DefUtils
  IMPLICIT NONE
  ! in
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL (KIND=dp) :: Array(4)
  ! internal
  REAL (KIND=dp) :: Time, x, y,
  INTEGER :: DIM

  SAVE FirstTime, DIM


  DIM = CoordinateSystemDimension()

  x=Array(1)
  x=Array(2)
  Time = Array(3)

  Slide=1.0

  RETURN
END FUNCTION linkZsGreen_newData
