
C     The values of these AMSR-E parameters come from E:\AMSR\Resampling_weights\PM1_AMSR_MOD.for
C     The values of th Midori2   parameters come from E:\AMSR\Resampling_weights\ADEOSII_AMSR_MOD.for
C     For these antenna parameters, Peter found one set for 89 GHZ, which was the averaged of 89A and 89B
C     This value was put in the IFREQ=8 slot
C     The IFREQ=6 and 7 slots were set to 0
C     For Midori2, the IFREQ=6 and 7 slots contain the 50.3 and  52.8 patterns respectively.

C     Here, I use slot 7 for 89A and slot 8 for 89B, so they have the same value
C     When I do Midori-2, I will need to decide what to do about the 50 GHz channels (maybe put in average value into slot 6)
  
 	  REAL(8), PARAMETER :: SIN_ANGLIM(8)=(/  0.105D0,  0.105D0,  0.0524D0, 0.0524D0, 0.0262D0, 0.000D0,  0.0262D0 ,  0.0262D0 /)
      REAL(8), PARAMETER :: ant_approx_a(8)=(/4.343d-6, 2.096d-6,	1.890d-6, 1.623d-6, 7.248d-7, 0.000d0,  2.070d-6,   2.070d-6 /)
      REAL(8), PARAMETER :: ant_approx_b(8)=(/6.892d-4, 4.059d-4,	3.727d-4, 7.251d-4, 3.051d-4, 0.000d0,  2.381d-4,   2.381d-4 /)
      REAL(8), PARAMETER :: ant_approx_c(8)=(/0.503d0,  0.662d0,  1.391d0,  1.804d0,  1.964d0,  0.000d0,  4.593d0,    4.593d0  /)
      REAL(8), PARAMETER :: ant_approx_d(8)=(/0.651d0,  1.345d0, 	4.844d0,  4.721d0,  15.18d0,  0.000d0, 79.785d0,   79.785d0  /)

	real(8) coeff_a,coeff_b,coeff_c,coeff_d,anglim,qsum,gain

	CALL OPENBIG(4,'show_antenna_patterns.dat','replace')

	DO JFREQ=1,8
	ANGLIM=asind(SIN_ANGLIM(JFREQ))
	coeff_a =  ant_approx_a(JFREQ)
	coeff_b =  ant_approx_b(JFREQ)
	coeff_c =  ant_approx_c(JFREQ)
	coeff_d =  ant_approx_d(JFREQ)

	QSUM=0; DELTASV=-1
	DO I=0,70000
	DELTA=0.0001D0*I

      IF (DELTA.LE.ANGLIM) THEN
	gain = coeff_a + coeff_b*exp(-coeff_c*delta) + 	exp(-coeff_d*delta*delta)
	IF(I.EQ.0) GAIN0=GAIN
	else
	gain=1.e-6
	endif

	WRITE(4) GAIN

	QSUM=QSUM + GAIN*SIND(DELTA)
	IF(ABS(GAIN/GAIN0 -0.5).LT.0.001) DELTASV=DELTA
	ENDDO !I

	SUMGAIN=QSUM*6.283185D0*1.745329D-6	  !TWO PI TIMES INTEGRATION STEP IN RADIANS (0.001 DEG)

	WRITE(6,1001) JFREQ, SUMGAIN, 2*DELTASV,
     & 10*DLOG10((coeff_a + coeff_b*exp(-coeff_c*ANGLIM) + exp(-coeff_d*ANGLIM*ANGLIM))/GAIN0) 
 1001 FORMAT(I3,E15.5,F8.3,F8.2)
	ENDDO !JFREQ

      STOP 'NORM END'
	END

	 					 		 		     
	 
	 													  
	INCLUDE 'X:\SYSPGM\OPENBIG.FOR'
