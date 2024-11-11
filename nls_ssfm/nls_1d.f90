
module FFTW3
  use, intrinsic :: iso_c_binding
  include 'fftw3.f03'
end module

MODULE COMM_DATA
  INTEGER, PARAMETER :: NRUN = 2*10000, LSTORE = 50, IMSTORE = NRUN / 5
  INTEGER, PARAMETER :: NX = 2048, NXX = NX-1, NX2 = NX/2
  REAL (8), PARAMETER :: LX = 20.48D0
  REAL (8), PARAMETER :: PI = 3.14159265358979D0, TWO_PI = 2.0D0 * PI 
  COMPLEX (8), PARAMETER :: CI_ = (0.0D0, 1.0D0)
  COMPLEX :: CI
END MODULE COMM_DATA

MODULE GPE_DATA
  USE COMM_DATA, ONLY : TWO_PI, LX, NX, PI
  REAL (8), PARAMETER :: DX = LX / NX
  REAL (8), PARAMETER :: DKX = TWO_PI/LX
!  			
  REAL (8), PARAMETER :: NU = 1.0D0, NU2 = NU * NU
  REAL (8), PARAMETER :: G0 = 1.0D0,TOL = 1E-5,MUL=PI/8.0D0

!
  REAL (8)    ::X(NX), X2(NX)
  REAL (8)    :: KX(NX), KX2(NX)
  REAL (8)    :: V(NX), R2(NX)
  COMPLEX (8) :: CP(NX), KIN(NX)
  
  REAL (8) :: DT, G
END MODULE GPE_DATA

PROGRAM ODGPE_FFT1D
  USE FFTW3 
  USE COMM_DATA, ONLY : NRUN, LSTORE, IMSTORE, CI, NX, NX2
  USE GPE_DATA
  IMPLICIT NONE
!------------------------ interface blocks -----------------------       
  INTERFACE 
   
    SUBROUTINE INITIALIZE()
    END SUBROUTINE INITIALIZE
 
    SUBROUTINE CALCULATE_TRAP()
    END SUBROUTINE CALCULATE_TRAP
  
    SUBROUTINE EODE(DT, CP)
      COMPLEX (8), DIMENSION(:), INTENT(INOUT) :: CP
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
      REAL (8)								   :: ZNORM
    END SUBROUTINE NORM
 
    SUBROUTINE EFFTW(PLAN1, PLAN2, KIN, U) 
      USE FFTW3, ONLY : C_PTR
      TYPE(C_PTR) :: PLAN1, PLAN2
      COMPLEX (8), DIMENSION(:), INTENT(IN) :: KIN
      COMPLEX (8), DIMENSION(:), INTENT(INOUT) :: U
    END SUBROUTINE EFFTW


  END INTERFACE
!-------------------- END INTERFACE BLOCKS -----------------------
  TYPE(C_PTR) :: PLAN1, PLAN2
  INTEGER :: I, L
  REAL (8):: CP2(NX)
  COMPLEX (8) :: S(NX), STMP(NX)
  REAL (8) :: RMS
  REAL (8) :: ZNORM, T, T1, DT2, T2 
  REAL (8) :: TI, TMPR, TMPC 
  INTEGER :: CLCK_COUNTS_BEG, CLCK_COUNTS_END, CLCK_RATE
!  
  CALL SYSTEM_CLOCK (CLCK_COUNTS_BEG, CLCK_RATE)
  CALL CPU_TIME (T1)
!
  DT = DX * DX
!  
  OPEN(7, file='nls-out.txt', status='unknown')
  WRITE(7, *) ' ****************** realtime propagation ******************'
  WRITE(7,*) 

  WRITE(7,*)'  Nonlinearity G =', G0
  WRITE(7,*)
  WRITE(7,*) ' # Space Stp: NX = ',NX
  WRITE(7,*) '  Space Step: DX = ',DX,'Time Step:   DT = ', DT
  WRITE(7,*)'--------------------------------------------------------------'

!
  CALL INITIALIZE()
!  CALL CALCULATE_TRAP()
!
  CALL NORM(CP, ZNORM)
!  
  G = G0
!  
  WRITE (7, *)
  WRITE (7, *) 'ZNORM', 'ABS(CP(0))' 
  WRITE (7, *)  ZNORM, ABS(CP(NX2)) 


! Initial Profile
  OPEN(50, file='nls-initial-profile.txt', status='unknown')
  DO I = 1, NX
    WRITE(50, 2001) X(I), ABS(CP(I)) 
  END DO
  CLOSE(50)
! 

   CI = (0.0D0, 1.0D0)
  
! 
  PLAN1 = FFTW_PLAN_DFT_1D(NX, S, STMP, FFTW_FORWARD, FFTW_ESTIMATE)
  PLAN2 = FFTW_PLAN_DFT_1D(NX, S, STMP, FFTW_BACKWARD, FFTW_ESTIMATE)
!  
  DO I = 1, NX
    KIN(I) = EXP(CI * DT * KX2(I) / 2.0D0)      
  END DO
!
  TI = 0.0D0
  DT2 = DT / 2.0D0 

!
  
	OPEN(15, file='real1d-dyn.txt', status='unknown')
      
	DO I = 1, NX, 1
        WRITE(15,*) X(I), TI, ABS(CP(I))
        WRITE(200,*) X(I), TI, ABS(CP(I))
    END DO
    WRITE (15,*)
    WRITE (200,*)
 
!TIME EVOLUTION

    DO L = 1, NRUN
        
		CALL EODE(DT2, CP)  
        CALL EFFTW(PLAN1, PLAN2, KIN, CP)
		CALL EODE(DT2, CP)  
		CALL NORM(CP, ZNORM)		
       
	    TI = TI + DT
		
		IF (ABS(MOD(TI, MUL)) < TOL) THEN
       
        DO I = 1, NX
            WRITE(200,*) X(I), TI, ABS(CP(I))
        END DO
        WRITE (200, *)

		END IF

		
		
		
		
        
		IF (MOD(L, LSTORE) == 0) THEN 
                   
		  DO I = 1, NX, 1
            WRITE(15,*) X(I), TI, ABS(CP(I))
          END DO
          WRITE (15, *)
           
          WRITE (7,*) ZNORM, ABS(CP(NX2)) 
        
		END IF
      
	END DO
      
	CLOSE(15)
    
!
  2001 FORMAT(F12.6, 4G17.8E3)

  WRITE(*,*) 'Strength of GF = ', G
    
      OPEN(45, file='nls-final-profile.txt', status='unknown')
      DO I = 1, NX
        WRITE(45,*) X(I), ABS(CP(I))
      END DO

	  CLOSE(45)

!
  CALL FFTW_DESTROY_PLAN(PLAN1)
  CALL FFTW_DESTROY_PLAN(PLAN2)
!
  CALL SYSTEM_CLOCK (CLCK_COUNTS_END, CLCK_RATE)
  CALL CPU_TIME(T2)
  WRITE (7,*)
  WRITE (7,'(A,I7,A)') ' Clock Time: ', (CLCK_COUNTS_END - CLCK_COUNTS_BEG)/INT (CLCK_RATE,8), ' seconds'
  WRITE (7,'(A,I7,A)') '   CPU Time: ', INT(T2-T1), ' seconds' 
  CLOSE(7)
  
END PROGRAM ODGPE_FFT1D
!

! 
SUBROUTINE INITIALIZE()
  USE COMM_DATA, ONLY : LX, NX, NXX, NX2, PI, TWO_PI
  USE GPE_DATA, ONLY : DX, NU, X, X2, R2, CP, DKX, KX, KX2
  IMPLICIT NONE
  REAL (8) :: PI2, PI4
  INTEGER :: I,A
!
  FORALL (I=1:NX) X(I) = - LX / 2.0D0 + (I-1) * DX
  X2 = X * X  
  DO I = 1, NX2
    KX(I) = (I - 1) * TWO_PI / LX
    KX(NX2+I) = -(NX2 - I + 1) * TWO_PI / LX 
  END DO
!
  KX2 = KX * KX
  DO I = 1, NX
    R2(I) = X2(I)
  END DO
  
  A=1
  
  DO I = 1, NX
	CP(I)=2.0D0/(EXP(X(I))+EXP(-X(I)))
    CP(I) = A*CP(I)
  END DO
END SUBROUTINE INITIALIZE

SUBROUTINE CALCULATE_TRAP()
  USE COMM_DATA, ONLY : NX
  USE GPE_DATA, ONLY :  NU2, V, X2
  IMPLICIT NONE
  INTEGER :: I
  
  DO I = 1, NX
    V(I) = NU2 * X2(I) / 2.0D0
  END DO
  
END SUBROUTINE CALCULATE_TRAP

SUBROUTINE EODE(DT, CP) ! Exact solution
  USE COMM_DATA, ONLY : CI_, CI, NX, NXX
  USE GPE_DATA, ONLY : V, G
  IMPLICIT NONE
!-------------------------------------------------
  COMPLEX (8), DIMENSION(:), INTENT(INOUT) :: CP
  REAL (8), INTENT(IN) :: DT
  REAL (8), DIMENSION(SIZE(CP)) :: P2
  COMPLEX (8), DIMENSION(SIZE(CP)) :: TMP1D
  COMPLEX (8) :: GMCIG1, CDT
  INTEGER :: I
!--------------------------------------------------------
   
  P2 = REAL(CP) * REAL(CP) + AIMAG(CP) * AIMAG(CP)
  GMCIG1 = G
  CDT = CI * DT
  
  DO I = 1, NX
  
    TMP1D(I) = GMCIG1 * P2(I)
    CP(I) = CP(I) * EXP(-CDT * TMP1D(I))
	
  END DO
  
END SUBROUTINE EODE

SUBROUTINE EFFTW(PLAN1, PLAN2, KIN, U) 
  USE COMM_DATA, ONLY : NX
  USE FFTW3, ONLY : C_PTR, C_INT, FFTW_EXECUTE_DFT
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(:), INTENT(INOUT) :: U
  COMPLEX (8), DIMENSION(:), INTENT(IN) :: KIN
  COMPLEX (8), DIMENSION(NX) :: S, STMP
  TYPE(C_PTR) :: PLAN1, PLAN2
  S = U
  CALL FFTW_EXECUTE_DFT(PLAN1, S, STMP)
  S = KIN * STMP
  CALL FFTW_EXECUTE_DFT(PLAN2, S, STMP)
  U = STMP / NX
END SUBROUTINE EFFTW
!-------------------------------------------------------------------------------------------------------------
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

!----------------------------------------------------------------------------------------------
SUBROUTINE NORM(CP, ZNORM)
!  Calculates the normalization of the wave function and sets it to unity.
  USE COMM_DATA, ONLY :  NX
  USE GPE_DATA, ONLY : DX
  !USE UTIL, ONLY : SIMP
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(:), INTENT(INOUT) :: CP
  REAL (8) :: ZNORM

  REAL (8), DIMENSION(SIZE(CP)) :: P2
!
    interface
	
	SUBROUTINE SIMP(F, DX,VALUE1)
		REAL (8), DIMENSION(0:), INTENT(IN) :: F
		REAL (8), INTENT(IN) :: DX
		REAL (8) :: VALUE1
	END SUBROUTINE SIMP
  
   end interface



  P2 = REAL(CP) * REAL(CP) + AIMAG(CP) * AIMAG(CP)
  !ZNORM = SIMP(P2, DX)
  CALL SIMP(P2,DX,ZNORM)
  
  
 END SUBROUTINE NORM




