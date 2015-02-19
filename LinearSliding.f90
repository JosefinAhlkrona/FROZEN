SUBROUTINE LinearSliding( Beta1, Beta2, n1,n2,n3, txz,tyz,p, uB1,uB2,uB3)

 
  !------------------------------------------------------------------------------
  ! This routine is called for each bedrock node
  !------------------------------------------------------------------------------


  USE ElementUtils

  IMPLICIT NONE

  REAL(KIND=dp) ::  Beta1,Beta2 
  REAL(KIND=dp) :: n1,n2,n3
  REAL(KIND=dp) :: txz,tyz,p
  REAL(KIND=dp) :: Tangent(3),Tangent2(3)
  REAL(KIND=dp) :: uB_t1,uB_t2,TauB_t1,TauB_t2
  REAL(KIND=dp), INTENT(OUT) :: uB1,uB2,uB3


  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------

 

 !------------------------------------------------------------------------------
 !    set the bottom velocity for current node, 
 !       input: beta, normal vec, nodenumber. return: uB
 !------------------------------------------------------------------------------

 

!     Tangent(1) =  Normal(2)   !2D case
!     Tangent(2) = -Normal(1)   !2D case
!     Tangent(3) =  0.0_dp      !2D case
!     Tangent2   =  0.0_dp      !2D case

  CALL TangentDirections((/ n1, n2, n3/), Tangent, Tangent2 ) 


!WRITE(*,*) 'n1=', n1
!WRITE(*,*) 'n2=', n2
!WRITE(*,*) 'n3=', n3
!WRITE(*,*) 'Tangent(1)=', Tangent(1)
!WRITE(*,*) 'Tangent(2)=', Tangent(2)
!WRITE(*,*) 'Tangent(3)=', Tangent(3)
!WRITE(*,*) 'Tangent2(1)=', Tangent2(1)
!WRITE(*,*) 'Tangent2(2)=', Tangent2(2)
!WRITE(*,*) 'Tangent2(3)=', Tangent2(3)
!WRITE(*,*) 'txz=', txz
!WRITE(*,*) 'tyz=', tyz
!WRITE(*,*) 'Beta1=', Beta1
!WRITE(*,*) 'Beta2=', Beta2



  TauB_t1 = Tangent(1)*(-n1*p+n3*txz)+Tangent(2)*(-n2*p+n3*tyz)+ &
       Tangent(3)*(n1*txz+n2*tyz-n3*p)

  TauB_t2 = Tangent2(1)*(-n1*p+n3*txz)+Tangent2(2)*(-n2*p+n3*tyz)+ &
       Tangent2(3)*(n1*txz+n2*tyz-n3*p)


!WRITE(*,*) 'TauB_t1=', TauB_t1
!WRITE(*,*) 'TauB_t2=', TauB_t2

  uB_t1=-TauB_t1/Beta1
  uB_t2=-TauB_t2/Beta2

  uB1=uB_t1*Tangent(1)+uB_t2*Tangent2(1)
  uB2=uB_t1*Tangent(2)+uB_t2*Tangent2(2)
  uB3=uB_t1*Tangent(3)+uB_t2*Tangent2(3)

!WRITE(*,*) 'ub1=',uB1
!WRITE(*,*) 'ub2=',uB2
!WRITE(*,*) 'ub3=',uB3

  !------------------------------------------------------------------------------
END SUBROUTINE LinearSliding
!------------------------------------------------------------------------------

