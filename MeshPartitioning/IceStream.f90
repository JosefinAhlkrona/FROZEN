FUNCTION IceStream (Model, nodenumber, Array) RESULT(Slide)
  USE DefUtils
  IMPLICIT NONE
  ! in
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL (KIND=dp) :: Array(3), a, b, c, d, theta
  ! internal
  REAL (KIND=dp) :: Time, x, y, Slide
  INTEGER :: DIM

  SAVE DIM


  DIM = CoordinateSystemDimension()

  x=Array(1)
  y=Array(2)
  Time = Array(3)
    a=-3.0!1e-5-0.25
    b=0.5*Time
    c=0.3!0.07
    d=-0.7

IF (x==0) THEN
   IF (y>0) THEN 
	theta = 3.14159265/2.0
   ELSE 
	theta = -3.14159265/2.0
   ENDIF
ELSEIF (x>0) THEN
  IF (y==0) THEN 
    theta=0
  ELSE 
   theta=atan(y/x)
  ENDIF
ELSE !x<0
   IF (y>=0) THEN
      theta = atan(y/x)+3.14159265359  
   ELSE  
      theta=atan(y/x)-3.14159265359       
   ENDIF
ENDIF


Slide = a*exp(-(theta-b)**2.0/(2.0*c**2.0))+d

Slide=10**Slide


IF (Slide < 1.0E-04) THEN
   Slide = 1.0E-04
ENDIF

  RETURN
END FUNCTION IceStream
