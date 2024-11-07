MODULE  FFTW3
  use, intrinsic :: iso_c_binding
  include 'fftw3.f03'
END MODULE

MODULE COMM_DATA

	integer, parameter 			:: NX = 2048,NX2=NX/2
	double precision, parameter :: LX=20.48D0
	double precision, parameter :: PI= 4.0D0*atan(1.0d0),TWO_PI=2.0D0*PI
	
	double complex	, parameter :: CI=(0.0D0,1.0D0)

END MODULE COMM_DATA

MODULE GPE_DATA
  
  USE COMM_DATA, ONLY : TWO_PI, LX, NX, PI
  
  double precision, parameter :: dx=LX/real(NX)
  double precision, parameter :: dt=dx*dx
  double precision, parameter :: dkx = TWO_PI/LX
  			
  double precision, parameter :: NU = 1.0D0, NU2 = NU * NU
  double precision, parameter :: G=0.0D0

  double precision ::x(NX),x2(NX)
  double precision ::kx(NX),kx2(NX)
  double precision ::V(NX),R2(NX)
  double complex   ::CP(NX),KIN(NX),HX(NX),S(NX),STMP(NX)
  
END MODULE GPE_DATA

program spectral_gpe

!------------------------------------------------------------------------------------------------------------------------
	use FFTW3
	use COMM_DATA 
	use GPE_DATA 
	implicit none
!----------------------------------------------------------------------------------------------------------------------------	
	integer :: NRUN=2*10000,LSTORE=1000
	integer :: i,L,step_counter=0.0
	double precision :: t,t_max,tsteps,T1,T2,ZNORM, MU, EN,RMS
	double precision ::CP2(NX)
	
	TYPE(C_PTR) :: PLAN1, PLAN2
	INTEGER :: CLCK_COUNTS_BEG, CLCK_COUNTS_END, CLCK_RATE
	
!-----------------------------------------------------------------------------------------------------------------------------
	
	INTERFACE 
    
		SUBROUTINE INITIALIZE()
		END SUBROUTINE INITIALIZE
 
		SUBROUTINE CALCULATE_TRAP()
		END SUBROUTINE CALCULATE_TRAP
  
		SUBROUTINE EODE(DT, CP, HX)
			COMPLEX (8), DIMENSION(:), INTENT(INOUT) :: CP
			COMPLEX (8), DIMENSION(:), INTENT(OUT) :: HX
			REAL (8), INTENT(IN) :: DT
		END SUBROUTINE EODE
		
		SUBROUTINE SIMP(F, DX,VALUE1)
			REAL (8), DIMENSION(0:), INTENT(IN) :: F
			REAL (8), INTENT(IN) :: DX
			REAL (8) :: VALUE1
		END SUBROUTINE SIMP
		
		SUBROUTINE DIFF(P,DX,DP) 
			REAL (8), DIMENSION(0:), INTENT(IN) :: P
			REAL (8), INTENT(IN) :: DX
			REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
		END SUBROUTINE DIFF
		
		SUBROUTINE NORM(CP, ZNORM)
			COMPLEX (8), DIMENSION(:), INTENT(INOUT) :: CP
			REAL (8), INTENT(OUT) :: ZNORM
		END SUBROUTINE NORM
		
		SUBROUTINE CHEM(CP, MU, EN)
			  COMPLEX (8), DIMENSION(:), INTENT(IN) :: CP
			  REAL (8), INTENT(OUT) :: MU, EN
		END SUBROUTINE CHEM
		
		SUBROUTINE RAD(CP2, R)
			REAL (8), DIMENSION(:), INTENT(IN) :: CP2
			REAL (8), INTENT(OUT) :: R
		END SUBROUTINE RAD
		
		SUBROUTINE EFFTW(PLAN1, PLAN2, KIN, U,HX)
			USE FFTW3, ONLY : C_PTR, C_INT, FFTW_EXECUTE_DFT
			TYPE(C_PTR) :: PLAN1, PLAN2
			COMPLEX (8), DIMENSION(:), INTENT(IN) :: KIN,HX
			COMPLEX (8), DIMENSION(:), INTENT(INOUT) :: U
		END SUBROUTINE EFFTW

	END INTERFACE
	
!-----------------------------------------------------------------------------------------------------------------------------

	CALL SYSTEM_CLOCK (CLCK_COUNTS_BEG, CLCK_RATE)
	CALL CPU_TIME (T1)
	
	OPEN(7, file='real1d-out.txt')
	!------------------------------------------------------------------------
	WRITE(7,*) ' ****************** realtime propagation ******************'
	
	WRITE(7,*) '  Nonlinearity G =',G
	WRITE(7,*) '  Parameters of trap: NU = ',NU
	WRITE(7,*)
	WRITE(7,*)' # Space Stp: NX = ', NX
	WRITE(7,*) '  Space Step: DX = ',dx, 'Time Step:   DT = ', dt
	
	WRITE(7,*)
	!---------------------------------------------------------------------------
	
	CALL INITIALIZE()
	CALL CALCULATE_TRAP()
	CALL NORM(CP, ZNORM)
	CALL CHEM(CP, MU, EN)
	CP2 = REAL(CP) * REAL(CP) + AIMAG(CP) * AIMAG(CP)
	CALL RAD(CP2, RMS)
	! 

	t=0.0D0
  	WRITE (7,*)'--------------------------------------------------------------------------------------------------'
	WRITE (7,*)"t","						       ",'Norm',"				       ",'MU'," 			      ",' EN',"			       ",'RMS',"			       ",'|Psi(0)|^2'
	WRITE (7,*)'---------------------------------------------------------------------------------------------------'
	WRITE (7,*) t,ZNORM, MU, EN, RMS, CP2(NX2) 
	
	!INITIAL PROFILE
	OPEN(50, file='real1d-initial.txt')
	DO i = 1, NX
    WRITE(50,*) x(i),  ABS(CP(i))
	END DO
	CLOSE(50)
	
	PLAN1 = fftw_plan_dft_1d(NX, S,  STMP, FFTW_FORWARD, FFTW_ESTIMATE)
	PLAN2 = fftw_plan_dft_1d(NX, S,  STMP, FFTW_BACKWARD, FFTW_ESTIMATE)
	
	 DO i = 1, NX
		KIN(i) = EXP(-CI * dt * kx2(i) / 2.0D0)      
	END DO
	
    OPEN(15, file='real1d-dyn.txt', status='unknown')
    OPEN(16, file='real1d-g.txt', status='unknown')
      
	  t=0.0D0
	  
	  
	  DO i = 1, NX, 5
        WRITE(15, *) x(i),t, ABS(CP(i))
      END DO
	

!	do while (t < t_max)
	do L=1,NRUN
	
        CALL EODE(dt, CP,HX)  
        CALL EFFTW(PLAN1, PLAN2, KIN, CP,HX)
		CALL NORM(CP, ZNORM)
       
        t = t + dt
		!step_counter=step_counter+1
		!write(100,*)L,"		",t
		!WRITE (7,*) t,ZNORM,CP2(NX2) 
		
        IF (MOD(L, LSTORE) == 0) THEN 
          
		  CP2 = REAL(CP) * REAL(CP) + AIMAG(CP) * AIMAG(CP)
          CALL CHEM(CP, MU, EN)
          CALL RAD(CP2, RMS)
		  WRITE (7,*) t,ZNORM, MU, EN, RMS, CP2(NX2)
		  WRITE (16,*) t,RMS,CP2(NX2)
		  
          DO i = 1, NX, 5
            WRITE(15, *) x(i),t,ABS(CP(i))
          END DO
			WRITE(15, *)
        END IF
		
	end do
	
	CLOSE(16)
	CLOSE(15)
	
	
	!FINAL PROFILE
	OPEN(60, file='real1d-final.txt')
	DO i = 1, NX
    WRITE(60,*) x(i),  ABS(CP(i))
	END DO
	CLOSE(60)
	
	
  CALL FFTW_DESTROY_PLAN(PLAN1)
  CALL FFTW_DESTROY_PLAN(PLAN2)
!
  CALL SYSTEM_CLOCK (CLCK_COUNTS_END, CLCK_RATE)
  CALL CPU_TIME(T2)
  WRITE (7,*)
  WRITE (7,'(A,I7,A)') ' Clock Time: ', (CLCK_COUNTS_END - CLCK_COUNTS_BEG)/INT (CLCK_RATE,8), ' seconds'
  WRITE (7,'(A,I7,A)') '   CPU Time: ', INT(T2-T1), ' seconds' 
  CLOSE(7)

end program
!-----------------------------------------------------------------------------------------------------------------------------

SUBROUTINE INITIALIZE()
 
	USE COMM_DATA, ONLY : LX, NX, NX2, PI, TWO_PI
	USE GPE_DATA, ONLY : x,x2,kx,kx2,R2,CP,NU,NU2,dx
  
	IMPLICIT NONE
	REAL (8) :: PI2, PI4
	INTEGER :: i

	!PI2 = SQRT(PI/NU)
	!PI4 = SQRT(PI2)
	PI4=SQRT(SQRT(4*PI))
	do i=1,NX
		x(i) = - LX / 2.0D0 + real(i-1) * dx
		x2(i) = x(i) * x(i)
	end do	
  
	do i = 1, NX2
		kx(i) = real(i - 1) * TWO_PI / LX 
	end do
	do i=NX2+1,NX
		kx(i)= real((i - 1)-NX)*TWO_PI/LX
	end do

	kx2=kx*kx

	do i = 1, NX
		R2(i) = x2(i)
	end do
  
	do i = 1, NX
		CP(i) = EXP(-NU2 * x2(i) / 2.0D0) / PI4
		CP(i)=(1.0D0+sqrt(2.0D0)*x(i))*CP(i)
	end do

END SUBROUTINE INITIALIZE


SUBROUTINE CALCULATE_TRAP()
  USE COMM_DATA, ONLY : NX
  USE GPE_DATA, ONLY :  NU2, V, x2
  IMPLICIT NONE
  INTEGER :: i
  
  DO i = 1, NX
    V(i) = NU2 * x2(i) / 2.0D0
  END DO
  
END SUBROUTINE CALCULATE_TRAP
!-------------------------------------------------------------------------------------------------------
SUBROUTINE SIMP(F, DX,VALUE1)

  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: F
  REAL (8), INTENT(IN) :: DX
  REAL (8) :: VALUE1

  REAL (8) :: F1, F2
  INTEGER :: I, N

  N = SIZE(F) - 1
  
  F1 = F(1) + F(N-1) 
  F2 = F(2) 
  DO I = 3, N-3, 2
     F1 = F1 + F(I)
     F2 = F2 + F(I+1)
  END DO
  VALUE1 = DX*(F(0) + 4.0D0*F1 + 2.0D0*F2 + F(N))/3.0D0

END SUBROUTINE SIMP

SUBROUTINE DIFF(P,DX,DP) 
  
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: P
  REAL (8), INTENT(IN) :: DX
  REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
  INTEGER :: I, N
  N = SIZE(P) - 1
  DP(0) = 0.0D0
  DP(1) = (P(2) - P(0))/(2.0D0*DX)
  DO I=2, N-2
    DP(I) = (P(I-2)-8.0D0*P(I-1)+8.0D0*P(I+1)-P(I+2))/(12.0D0*DX)
  END DO
  DP(N-1) = (P(N) - P(N-2))/(2.0D0*DX)
  
  DP(N) = 0.0D0
END SUBROUTINE DIFF
!-----------------------------------------------------------------------------------------------------------

SUBROUTINE NORM(CP, ZNORM)
!  Calculates the normalization of the wave function and sets it to unity.
  USE COMM_DATA, ONLY : NX
  USE GPE_DATA, ONLY : dx
 
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: ZNORM
  REAL (8), DIMENSION(SIZE(CP)) :: P2
  
  interface
	
	SUBROUTINE SIMP(F, DX,VALUE1)
		REAL (8), DIMENSION(0:), INTENT(IN) :: F
		REAL (8), INTENT(IN) :: DX
		REAL (8) :: VALUE1
	END SUBROUTINE SIMP
  
  end interface
  
!
  P2 = REAL(CP) * REAL(CP) + AIMAG(CP) * AIMAG(CP)
  CALL SIMP(P2,dx,ZNORM)
 
 END SUBROUTINE NORM

SUBROUTINE CHEM(CP, MU, EN)
  USE COMM_DATA, ONLY : CI
  USE GPE_DATA, ONLY : V, G, dx
  !USE UTIL, ONLY : DIFF, SIMP
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(:), INTENT(IN) :: CP
  REAL (8), INTENT(OUT) :: MU, EN
!----------------------------------------------------------------------
  INTERFACE 
    SUBROUTINE DIFF(P, DX,DP)
      REAL (8), DIMENSION(0:), INTENT(IN) :: P
      REAL (8), INTENT(IN) :: DX
      REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
    END SUBROUTINE DIFF
  END INTERFACE
  
 interface
	
	SUBROUTINE SIMP(F, DX,VALUE1)
		REAL (8), DIMENSION(0:), INTENT(IN) :: F
		REAL (8), INTENT(IN) :: DX
		REAL (8) :: VALUE1
	END SUBROUTINE SIMP
  
  end interface
!----------------------------------------------------------------------
  REAL (8), DIMENSION(SIZE(CP)) :: DCPR,DCPI
  REAL (8), DIMENSION(SIZE(CP)) :: P2,DP2, DPR2,DPI2, GP2, G1P2
  REAL (8), DIMENSION(SIZE(CP)) :: TMP1D, EMP1D
!
  CALL DIFF(REAL(CP),dx,DCPR)
  DPR2 = DCPR * DCPR
  CALL DIFF(AIMAG(CP),dx,DCPI)
  DPI2 = DCPI * DCPI
  
  DP2 = (DPR2  + DPI2) / 2.0D0
  P2 = REAL(CP) * REAL(CP) + AIMAG(CP) * AIMAG(CP)
  GP2 = G * P2  
  TMP1D = (V + GP2 )*P2 + DP2 
  EMP1D = (V + (GP2/2.0D0))*P2 + DP2
!    
  CALL SIMP(TMP1D, dx,MU)
  CALL SIMP(EMP1D, dx,EN)
   
END SUBROUTINE CHEM


SUBROUTINE RAD(CP2, R)
! Calculates the root mean square size RMS
  USE COMM_DATA, ONLY : NX
  USE GPE_DATA, ONLY : dx, x2
  !USE UTIL, ONLY : SIMP
  IMPLICIT NONE
  REAL (8), DIMENSION(:), INTENT(IN) :: CP2
  REAL (8)							 :: R1
  REAL (8), INTENT(OUT) :: R

  INTEGER :: I
  REAL (8), DIMENSION(SIZE(CP2)) :: TMP1D
  
    interface
	
	SUBROUTINE SIMP(F, DX,VALUE1)
		REAL (8), DIMENSION(0:), INTENT(IN) :: F
		REAL (8), INTENT(IN) :: DX
		REAL (8) :: VALUE1
	END SUBROUTINE SIMP
  
    end interface
 
  DO I = 1, NX
    TMP1D(I) = x2(I) * CP2(I)
  END DO
  
   CALL SIMP(TMP1D, dx,R1)
   R=sqrt(R1)
  
  
END SUBROUTINE RAD

SUBROUTINE EODE(DT, CP,HX) ! Exact solution
 
  USE COMM_DATA, ONLY : CI, NX
  USE GPE_DATA, ONLY : V, G
  IMPLICIT NONE
!-------------------------------------------------
  COMPLEX (8), DIMENSION(:), INTENT(INOUT) :: CP
  REAL (8), INTENT(IN) :: DT
  COMPLEX (8), DIMENSION(:), INTENT(OUT) :: HX
  REAL (8), DIMENSION(SIZE(CP)) :: P2
  COMPLEX (8), DIMENSION(SIZE(CP)) :: TMP1D
  COMPLEX (8) :: CDT
  INTEGER :: i
!--------------------------------------------------   
  P2 = REAL(CP) * REAL(CP) + AIMAG(CP) * AIMAG(CP)
  CDT = (CI * DT)/2.0D0

  DO i = 1, NX
    TMP1D(i) = V(i) + G * P2(i)
    HX(i) = EXP(-CDT * TMP1D(i))
    CP(i) = HX(i)*CP(i)
  END DO
 
END SUBROUTINE EODE

SUBROUTINE EFFTW(PLAN1, PLAN2, KIN, U,HX) 
  USE FFTW3, ONLY : C_PTR, C_INT, FFTW_EXECUTE_DFT
  USE COMM_DATA, ONLY : NX
  USE GPE_DATA, ONLY : dt
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(:), INTENT(INOUT) :: U
  COMPLEX (8), DIMENSION(:), INTENT(IN) :: KIN,HX
  COMPLEX (8), DIMENSION(NX) :: S, STMP
  TYPE(C_PTR) :: PLAN1, PLAN2
 
  S = U
  CALL FFTW_EXECUTE_DFT(PLAN1, S, STMP)
  S = KIN * STMP
  CALL FFTW_EXECUTE_DFT(PLAN2, S, STMP)
  U = STMP / NX
  
  U=HX*U
  
END SUBROUTINE EFFTW
