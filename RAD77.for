C=======================================================================

      PROGRAM RAD TWO

C This program is an interactive numerical analysis code for the study
C of linear electroacoustic systems using the T-matrix technique.
C
C Program RAD was written by Michael Lampton as part of a study of
C acoustic resonances in transmission-line and bass-reflex loudspeaker
C enclosures. It provides examples of techniques that are useful in
C interactively modifying parameter arrays, building a system T-matrix,
C and displaying the results.
C
C A screen plot of each configuration is standard. A printer plot is
C optional.
C
C Assign printer output to logical unit = 3.
C Assign user keyboard  to logical unit = 5.
C Assign user screen    to logical unit = 6.
C
C NOTES:
C
C This version of the program is a modified version of the original
C code. The main changes relate to RAD's conversion to use Fortran 77.
C The use of Holleriths has been eliminated and hopefully improves
C clarity for the user.
C
C REFERENCES:
C
C Lampton, M. An interactive numerical analysis code for linear
C electroacoustic systems. Acustica, Volume 42, Number 2, 1979,
C pages 115-120.
C
C Lampton, M. Transmission matrices in electroacoustics. Acta
C Acustica united with Acustica, Volume 39, Number 4, 1978,
C pp. 239-251.

      IMPLICIT NONE

      COMPLEX   UNITY(4), MINUS(4)
      COMPLEX   ZRAD, T(4), U(4), V(4)
      COMPLEX   RESP(128), ZIN(128), DET(128)
      REAL      ENCL(4), VENT(4), WOOF(7), ENTRY(8)
      REAL      A, F, FMAX, FMIN, RRAD, PI, TWOPI, FOURPI, W, RHOAIR
      CHARACTER HC*1, HD*1, HF*1, HL*1, HP*1, HR*1, HZ*1
      CHARACTER RZ*1, BC*1, PR*1, WHAT*1
      INTEGER   I, IBAD, J, NUM

      DATA UNITY /( 1.0,0.0), (0.0,0.0), (0.0,0.0), ( 1.0,0.0)/
      DATA MINUS /(-1.0,0.0), (0.0,0.0), (0.0,0.0), (-1.0,0.0)/
      DATA WOOF  /8.0, 0.001, 8.8, 0.0225, 1.4, 800.0, 0.0320/
      DATA ENCL  /1.00, 0.1270, 100.0, 0.0/
      DATA VENT  /0.26, 0.0100, 100.0, 0.0/

      DATA FMIN /10.0/ FMAX/1000.0/
      DATA NUM  /101/

C.... CHARACTERS BELOW ARE TO BE RECOGNIZED FROM KEYBOARD.

      DATA HC/'C'/,HD/'D'/,HF/'F'/,HL/'L'/,HP/'P'/,HR/'R'/,HZ/'Z'/

C.... SOME USEFUL CONSTANTS.

      PI     = 4.0*ATAN(1.0)
      TWOPI  = 2.0*PI
      FOURPI = 4.0*PI
      RHOAIR = 1.2041

C.... FLASH TITLE AND INSTRUCTIONS ONTO THE SCREEN.

      WRITE(6, 65)
   65 FORMAT(//' RAD TWO INITIALIZATION'/
     & ' ENTER CONSTANT ARRAYS USING 8F8 FORMAT.'/
     & ' ZEROS OR BLANKS LEAVE DEFAULT OR PREVIOUS VALUES.'/
     & ' TO ZERO A PARAMETER, ENTER ANY NEGATIVE NUMBER.'/)

C.... COLLECT ALL OF THE NEEDED ARRAYS.

  110 WRITE(6,'(1X,A)') 'WOOFER:'
      WRITE(6,112) WOOF
  112 FORMAT(/,F8.1, F8.3, F8.1, F8.4, F8.1, F8.1, F8.4)
      WRITE(6,115)
  115 FORMAT('      RE       L      BL    MASS      RM   STIFF    AREA')
      READ(5,100) ENTRY
  100 FORMAT(8F8.4)
      DO J=1,7
        IF(ENTRY(J).LT.0.0) WOOF(J) = 0.0
        IF(ENTRY(J).GT.0.0) WOOF(J) = ENTRY(J)
      END DO

  120 WRITE(6,'(1X,A)') 'ENCLOSURE:'
      WRITE(6,122) ENCL
  122 FORMAT(/, F8.2, F8.4, F8.1, F8.2)
      WRITE(6,125)
  125 FORMAT('       L      AM     RES       M')
      READ(5,100) ENTRY
      DO J=1,4
        IF(ENTRY(J).LT.0.0) ENCL(J) = 0.0
        IF(ENTRY(J).GT.0.0) ENCL(J) = ENTRY(J)
      END DO

      WRITE(6,'(1X,A)') 'VENT:'
      WRITE(6,122) VENT
      WRITE(6,125)
      READ(5,100) ENTRY
      DO J=1,4
        IF(ENTRY(J).LT.0.0) VENT(J) = 0.0
        IF(ENTRY(J).GT.0.0) VENT(J) = ENTRY(J)
      END DO

C.... COMPLETES FILLING THE PARAMETER ARRAYS.

  170 WRITE(6,'(1X,A)') 'FREQUENCY RANGE:'
      WRITE(6,180) FMIN, FMAX
  180 FORMAT(/,2F8.1,/,'    FMIN    FMAX')
      READ(5,100) ENTRY
      IF (ENTRY(1).GT.0.0) FMIN = ENTRY(1)
      IF (ENTRY(2).GT.0.0) FMAX = ENTRY(2)

C.... COMMENCE THE OVERALL FREQUENCY STEP LOOP.

      DO I=1,NUM
        F = FMIN * 10.0**(ALOG10(FMAX/FMIN)*FLOAT(I-1)/FLOAT(NUM))
        W = TWOPI * F
C...... T-MATRIX BUILDUP BEGINS HERE.
        CALL SPKR(WOOF, W, T, IBAD)
        IF(IBAD.GE.1) THEN
          WRITE(6,'(A)') ' ERROR IN SPKR SUBROUTINE.'
          GO TO 110
        END IF
        CALL EXPO(ENCL, W, U)
        CALL EXPO(VENT, W, V)
        CALL MULT(U, V, U)
        CALL MULT(MINUS, U, U)
        CALL SIPO(UNITY, U, V, IBAD)
        IF(IBAD.GE.1) THEN
          WRITE(6,'(A)') ' ERROR IN SIPO SUBROUTINE.'
          GO TO 120
        END IF
        CALL MULT(T, V, T)
C...... CONSTRUCT RADIATION LOAD.
        A = WOOF(7)
        RRAD = RHOAIR*W*W/TWOPI/342.0/(1.0 + 0.5*W*W*A*A/(342.0*342.0))
        ZRAD = CMPLX(RRAD, 0.0)
        CALL ZMAT(ZRAD, U)
        CALL MULT(T, U, T)
C...... RESPONSE IS IN UNITS OF PASCALS PER VOLT AT ONE METER.
        RESP(I) = CMPLX(0.0, RHOAIR/FOURPI) * W / T(2)
        ZIN(I) = T(2) / T(4)
        DET(I) = (T(1)*T(4)) - (T(2)*T(3))
      END DO

C.... FINISHED WITH THE CALCULATION. FIND OUT WHAT PLOTS ARE WANTED.

  300 WRITE(6,'(A)') ' RESP OR ZIN? TYPE R OR Z.'
      READ(5,'(A1)') RZ
      WRITE(6,'(A)') ' BODE OR COMPLEX-PLANE? TYPE B OR C.'
      READ(5,'(A1)') BC
      IF(RZ.EQ.HZ.AND.BC.EQ.HC) CALL CPLN(6,51,21, ZIN,NUM,FMIN,FMAX)
      IF(RZ.NE.HZ.AND.BC.EQ.HC) CALL CPLN(6,51,21,RESP,NUM,FMIN,FMAX)
      IF(RZ.EQ.HZ.AND.BC.NE.HC) CALL BODE(6,51,21, ZIN,NUM,FMIN,FMAX)
      IF(RZ.NE.HZ.AND.BC.NE.HC) CALL BODE(6,51,21,RESP,NUM,FMIN,FMAX)
      IF(RZ.EQ.HD)              CALL CPLN(6,51,21, DET,NUM,FMIN,FMAX)

C.... CRT PLOT DONE. WANT IT PRINTED?

      WRITE(6,*)
      WRITE(6,'(A)') ' WANT IT PRINTED TO LOGICAL UNIT 3? TYPE P.'
      READ(5,'(A1)') PR
      IF(PR.EQ.HP) THEN
        WRITE(3,'(1X,A)') 'RAD VERSION 2.4, REVISED 10 OCT 1977, MLL.'
        WRITE(3,112) WOOF
        WRITE(3,115)
        WRITE(3,122) ENCL
        WRITE(3,125)
        WRITE(3,122) VENT
        WRITE(3,125)
        WRITE(3,*)
        IF(RZ.EQ.HZ.AND.BC.EQ.HC) CALL CPLN(3, 72,51, ZIN,NUM,FMIN,FMAX)
        IF(RZ.NE.HZ.AND.BC.EQ.HC) CALL CPLN(3, 72,51,RESP,NUM,FMIN,FMAX)
        IF(RZ.EQ.HZ.AND.BC.NE.HC) CALL BODE(3,101,51, ZIN,NUM,FMIN,FMAX)
        IF(RZ.NE.HZ.AND.BC.NE.HC) CALL BODE(3,101,51,RESP,NUM,FMIN,FMAX)
      END IF

C.... DONE WITH PLOTTING. WHAT NEXT?

  350 WRITE(6,'(A)')
     &      ' MORE PLOTS? P. NEW FREQS? F. NEW CONSTANTS? C. DONE? D.'
      READ(5,'(A1)') WHAT
      IF(WHAT.EQ.HP) GO TO 300
      IF(WHAT.EQ.HF) GO TO 170
      IF(WHAT.EQ.HC) GO TO 110
      IF(WHAT.NE.HD) GO TO 350

      STOP
      END

C=======================================================================

      SUBROUTINE MULT(X,Y, XY)

      IMPLICIT NONE

      COMPLEX X(4), Y(4), XY(4)

      COMPLEX Z(4)
      INTEGER I

      Z(1) = X(1)*Y(1) + X(2)*Y(3)
      Z(2) = X(1)*Y(2) + X(2)*Y(4)
      Z(3) = X(3)*Y(1) + X(4)*Y(3)
      Z(4) = X(3)*Y(2) + X(4)*Y(4)

C.... CLEVER USE OF AN INTERMEDIATE MATRIX, Z, IN ORDER TO
C.... PERMIT CALLING SEQUENCES SUCH AS: CALL MULT(X,Y,X).

      DO I=1,4
        XY(I) = Z(I)
      END DO

      RETURN
      END

C=======================================================================

      SUBROUTINE ZMAT(Z, T)

      IMPLICIT NONE

      COMPLEX Z, T(4)

      T(1) = CMPLX(1.0,0.0)
      T(2) = Z
      T(3) = CMPLX(0.0,0.0)
      T(4) = CMPLX(1.0,0.0)

      RETURN
      END

C=======================================================================

      SUBROUTINE YMAT(Y, T)

      IMPLICIT NONE

      COMPLEX Y, T(4)

      T(1) = CMPLX(1.0,0.0)
      T(2) = CMPLX(0.0,0.0)
      T(3) = Y
      T(4) = CMPLX(1.0,0.0)

      RETURN
      END

C=======================================================================

      SUBROUTINE EXPO(HORN,W, T)

C CONSTANTS IN HORN(4) ARE: LENGTH, MOUTH AREA, RESISTIVITY, FLARE RATE.

      IMPLICIT NONE

      REAL    HORN(4), W
      COMPLEX T(4)

      COMPLEX K, KSQ, COSKL, SINKL
      REAL    M, L, MOUTH, R, GPO, RHO, EML

      L     = AMAX1(0.001, HORN(1))
      MOUTH = AMAX1(1.0E-9, HORN(2))
      R     = HORN(3)
      M     = HORN(4)

      GPO = 1.40E5
      RHO = 1.2041
      EML = EXP(M*L/2.0)
      KSQ = CMPLX(RHO*W*W/GPO-M*M/4.0, -W*R/GPO)
      K = CSQRT(KSQ)
      IF(AIMAG(K).GT.0.0) K = CMPLX(-1.0,0.0)*K
      COSKL = CCOS(K*L)
      SINKL = CSIN(K*L)

      T(1) = (COSKL - CMPLX(0.5,0.0)*M/K*SINKL) * EML
      T(2) = CMPLX(1.0,0.0) + CMPLX(0.25,0.0)*M*M/KSQ
      T(2) = CMPLX(0.0,1.0)*SINKL*GPO*K/W/MOUTH * EML * T(2)
      T(3) = CMPLX(0.0,1.0)*SINKL/GPO/K*W*MOUTH / EML
      T(4) = (COSKL + CMPLX(0.5,0.0)*M/K*SINKL) / EML

      RETURN
      END

C=======================================================================

      SUBROUTINE SPKR(P,W, T,IERR)

C CALLING QUANTITIES MUST BE:
C
C  R, OHMS L,HY     BL, NIA M, KG    R,NS/M S,N/M AREA,SOM
C    P(1)   P(2)     P(3)    P(4)     P(5)    P(6)     P(7)

      IMPLICIT NONE

      REAL    P(7), W
      COMPLEX T(4)
      INTEGER IERR

      COMPLEX ZE, ZM
      REAL    BL, AREA

C.... CATCHES ERRONEOUS CALLING VALUES, FOR TRACING ERRORS.

      IF(P(3).EQ.0.0 .OR. P(7).LE.0.0) THEN
        IERR = 1
        RETURN
      END IF

      IERR = 0
      ZE = CMPLX(P(1), P(2)*W)
      ZM = CMPLX(P(5), P(4)*W - P(6)/W)
      BL = P(3)
      AREA = P(7)
      T(1) = ZE * AREA / BL
      T(2) = ZE*ZM/BL/AREA + BL/AREA
      T(3) = CMPLX(AREA/BL,0.0)
      T(4) = ZM / BL/ AREA

      RETURN
      END

C=======================================================================

      SUBROUTINE FREESP(D,W, T,IERR)

C COMPUTES FREE SPACE COUPLING IN FOUR PI STERADIANS.
C D = DISTANCE BETWEEN PORTS, IN METERS.

      IMPLICIT NONE

      REAL    D, W
      COMPLEX T(4)
      INTEGER IERR

      COMPLEX ZR
      REAL    C, RHOAIR, FOURPI

      IF(W.EQ.0.0 .OR. D.EQ.0.0) THEN
        IERR = 1
        RETURN
      END IF

      C      = 343.21
      RHOAIR = 1.2041
      FOURPI = 4.0*(4.0*ATAN(1.0))

      IERR = 0
      ZR = CMPLX(0.0,RHOAIR/FOURPI)*W/D * CEXP(CMPLX(0.0,-1.0/C)*W*D)
      T(1) = CMPLX(0.0,0.0)
      T(2) = CMPLX(-1.0,0.0)*ZR
      T(3) = CMPLX(1.0,0.0)/ZR
      T(4) = CMPLX(0.0,0.0)

      RETURN
      END

C=======================================================================

      SUBROUTINE SIPO(X,Y, Z,IERR)

C COMBINES T-MATRICES X AND Y WITH SERIES INPUTS AND PARALLEL OUTPUTS.

      IMPLICIT NONE

      COMPLEX X(4), Y(4), Z(4)
      INTEGER IERR

      COMPLEX DEN

      DEN = X(4) + Y(4)
      IF(CABS(DEN).EQ.0.0) THEN
        IERR = 1
        RETURN
      END IF

      IERR = 0
      Z(1) = (X(2)*Y(3)+X(3)*Y(2)-X(2)*X(3)-Y(2)*Y(3))/DEN+X(1)+Y(1)
      Z(2) = (X(4)*Y(2)+X(2)*Y(4))/DEN
      Z(3) = (X(3)*Y(4)+X(4)*Y(3))/DEN
      Z(4) = X(4)*Y(4)/DEN

      RETURN
      END

C=======================================================================

      SUBROUTINE SIZO(U,V,Z, T,IERR)

C                  * U * Z1
C  INPUT ** SERIES         Z3 ** OUTPUT
C                  * V * Z2

      IMPLICIT NONE

      COMPLEX U(4), V(4), Z(3,3), T(4)
      INTEGER IERR

      COMPLEX W(4)
      COMPLEX M11, M12, M13, M21, M22, M23, M33
      COMPLEX DETU, DETV, DETZ, DENOM
      INTEGER J

C.... COMPUTE THE MINORS OF Z.

      M11 = Z(2,2)*Z(3,3) - Z(2,3)*Z(3,2)
      M12 = Z(2,1)*Z(3,3) - Z(2,3)*Z(3,1)
      M13 = Z(2,1)*Z(3,2) - Z(2,2)*Z(3,1)
      M21 = Z(1,2)*Z(3,3) - Z(1,3)*Z(3,2)
      M22 = Z(1,1)*Z(3,3) - Z(1,3)*Z(3,1)
      M23 = Z(1,1)*Z(3,2) - Z(1,2)*Z(3,1)
      M33 = Z(1,1)*Z(2,2) - Z(1,2)*Z(2,1)
      DETZ = Z(1,1)*M11 - Z(1,2)*M12 + Z(1,3)*M13

C.... CHECK FOR ZERO DENOMINATOR.

      DENOM = U(4)*Z(3,2) + V(4)*Z(3,1) + U(3)*M23 - V(3)*M13
      IF(CABS(DENOM).EQ.0.0) THEN
        IERR = 1
        RETURN
      END IF

      IERR = 0
      DETU = U(1)*U(4) - U(2)*U(3)
      DETV = V(1)*V(4) - V(2)*V(3)

C.... CALCULATE THE NUMERATOR OF THE T MATRIX.

      W(1) = (U(1)*V(4)+U(3)*V(2))*Z(1,1) + (U(4)*V(1)+U(2)*V(3))*Z(2,2)
     &     + (U(1)*V(3) + U(3)*V(1))*M33 + DETU*Z(1,2) + DETV*Z(2,1)
     &     + U(2)*V(4) + U(4)*V(2)
      W(2) = (U(1)*V(3)+U(3)*V(1))*DETZ + DETU*M21 + DETV*M12
     &     + (U(2)*V(4)+U(4)*V(2))*Z(3,3)
     &     + (U(2)*V(3)+U(4)*V(1))*M11 + (U(1)*V(4)+U(3)*V(2))*M22
      W(3) = U(4)*V(4) + U(3)*V(4)*Z(1,1) + U(4)*V(3)*Z(2,2)
     &     + U(3)*V(3)*M33
      W(4) = U(4)*V(4)*Z(3,3) + U(3)*V(3)*DETZ
     &     + U(3)*V(4)*M22 + U(4)*V(3)*M11

C.... CALCULATE THE T MATRIX.

      DO J=1,4
        T(J) = W(J)/DENOM
      END DO

      RETURN
      END

C=======================================================================

      SUBROUTINE PIZO(U,V,Z, T,IERR)

C                    * U * Z1
C  INPUT ** PARALLEL          Z3 ** OUTPUT
C                    * V * Z2

      IMPLICIT NONE

      COMPLEX U(4), V(4), Z(3,3), T(4)
      INTEGER IERR

      COMPLEX W(4)
      COMPLEX M11, M12, M13, M21, M22, M23, M33
      COMPLEX DETU, DETV, DETZ, DENOM
      INTEGER J

C.... COMPUTE THE MINORS OF Z.

      M11 = Z(2,2)*Z(3,3) - Z(2,3)*Z(3,2)
      M12 = Z(2,1)*Z(3,3) - Z(2,3)*Z(3,1)
      M13 = Z(2,1)*Z(3,2) - Z(2,2)*Z(3,1)
      M21 = Z(1,2)*Z(3,3) - Z(1,3)*Z(3,2)
      M22 = Z(1,1)*Z(3,3) - Z(1,3)*Z(3,1)
      M23 = Z(1,1)*Z(3,2) - Z(1,2)*Z(3,1)
      M33 = Z(1,1)*Z(2,2) - Z(1,2)*Z(2,1)
      DETZ = Z(1,1)*M11 - Z(1,2)*M12 + Z(1,3)*M13

C.... CHECK FOR ZERO DENOMINATOR.

      DENOM = U(2)*Z(3,2) + V(2)*Z(3,1) + U(1)*M23 - V(1)*M13
      IF(CABS(DENOM).EQ.0.0) THEN
        IERR = 1
        RETURN
      END IF

      IERR = 0
      DETU = U(1)*U(4) - U(2)*U(3)
      DETV = V(1)*V(4) - V(2)*V(3)

C.... CALCULATE THE NUMERATOR OF THE T MATRIX.

      W(1) = U(1)*V(1)*M33+U(1)*V(2)*Z(1,1)+U(2)*V(1)*Z(2,2)+U(2)*V(2)
      W(2) = U(1)*V(1)*DETZ+U(1)*V(2)*M22+U(2)*V(1)*M11+U(2)*V(2)*Z(3,3)
      W(3) = (U(1)*V(3)+U(3)*V(1))*M33 + U(4)*V(2) + U(2)*V(4)
     &     + (U(1)*V(4)+U(3)*V(2))*Z(1,1) + (U(4)*V(1)+U(2)*V(3))*Z(2,2)
     &     - DETU*Z(1,2) - DETV*Z(2,1)
      W(4) = (U(1)*V(3)+U(3)*V(1))*DETZ + (U(4)*V(2)+U(2)*V(4))*Z(3,3)
     &     + (U(1)*V(4)+U(3)*V(2))*M22 + (U(2)*V(3)+U(4)*V(1))*M11
     &     - DETU*M21 - DETV*M12

C.... CALCULATE THE T MATRIX.

      DO J=1,4
        T(J) = W(4)/DENOM
      END DO

      RETURN
      END

C=======================================================================

      SUBROUTINE BODE(LU,NX,NY,DAT,NDAT,X1,X2)

C CONSTRUCTS A PLOT OF LOG MAGNITUDE AND PHASE OF A COMPLEX ARRAY.
C
C CALLING PARAMETERS ARE:
C
C LU = LOGIC UNIT ON WHICH OUTPUT WILL BE PLOTTED.
C NX, NY = WIDTH AND HEIGHT OF DESIRED PLOT.
C THIS VERSION IS DIMENSIONED TO ALLOW NX AS BIG AS 101.
C DAT, NDAT = COMPLEX ARRAY TO BE PLOTTED AND ITS SIZE.
C X1, X2 = REAL NUMBER LABELS PLACED BENEATH THE PLOT.
C AMPLITUDE RANGE IS FROM MAXIMUM DOWN TO 25 DB BELOW MAXIMUM.
C PHASE RANGE IS +180 DEGREES AT TOP TO -180 DEGREES AT BOTTOM.

      IMPLICIT NONE

      INTEGER LU, NX, NY, NDAT
      COMPLEX DAT(NDAT)
      REAL    X1, X2

      CHARACTER MARK(3)*1, LINE(101)*1
      INTEGER   I, IDAT, J
      REAL      DBDAT(128), DGDAT(128)
      REAL      DBBOT, DBLABL, DBMAX, DGBOT, DGLABL, DGTOP, PI, R2D
      REAL      YRANGE

      DATA MARK   /'.','X','*'/
      DATA YRANGE /25.0/

      DBMAX = -400.0
      PI = 4.0*ATAN(1.0)
      R2D = 180.0/PI
      DO I = 1,NDAT
        DBDAT(I) = 20.0*ALOG10(CABS(DAT(I)))
        DGDAT(I) = R2D*ATAN2(AIMAG(DAT(I)),REAL(DAT(I)))
        DBMAX = AMAX1(DBMAX,DBDAT(I))
      ENDDO

C.... COMMENCE PRINTING, LINE BY LINE.

      WRITE(LU,'(1X,A6,A7)') 'DEG', 'DB'
      DO J = 1,NY
        DBLABL = DBMAX - FLOAT(J-1)/FLOAT(NY-1)*YRANGE
        DGLABL = 180.0 - FLOAT(J-1)/FLOAT(NY-1)*360.0
        DBBOT = DBLABL - YRANGE/2.0/FLOAT(NY-1)
        DGTOP = DGLABL + 360.0/2.0/FLOAT(NY-1)
        DGBOT = DGLABL - 360.0/2.0/FLOAT(NY-1)
        DO I = 1,NX
          IDAT = 1 + ((NDAT-1)*(I-1))/(NX-1)
          LINE(I) = MARK(1)
          IF(DBDAT(IDAT).GT.DBBOT) LINE(I) = MARK(2)
          IF(DGDAT(IDAT).GT.DGBOT .AND. DGDAT(IDAT).LE.DGTOP)
     &      LINE(I) = MARK(3)
        ENDDO
        WRITE(LU,'(1X,F6.1,F7.1,2X,101A1)')
     &    DGLABL, DBLABL, (LINE(I),I=1,NX)
      ENDDO

      IF(NX.LT.75) THEN
        WRITE(LU,'(11X,F8.1,40X,F8.1)') X1, X2
      ELSE IF(NX.GE.75) THEN
        WRITE(LU,'(11X,F8.1,90X,F8.1)') X1, X2
      END IF

      END

C=======================================================================

      SUBROUTINE CPLN(LU,NX,NY,DAT,NDAT,X1,X2)

C CREATES A PLOT OF A GIVEN ARRAY IN THE COMPLEX PLANE.
C
C CALLING PARAMETERS ARE:
C
C LU = LOGIC UNIT ON WHICH OUTPUT WILL BE PLOTTED,
C NX, NY = WIDTH AND HEIGHT OF DESIRED PLOT,
C DAT, NDAT = COMPLEX ARRAY TO BE PLOTTED AND ITS SIZE.
C X1, X2 = REAL NUMBER LABELS PLACED BENEATH THE PLOT.

      IMPLICIT NONE

      INTEGER LU, NX, NY, NDAT
      COMPLEX DAT(NDAT)
      REAL X1, X2

      CHARACTER MARK(14)*1, LINE(101)*1
      REAL XMAX, XMIN, YMAX, YMIN, XSPAN, YSPAN, SPAN, XZERO
      REAL X, Y, YLABL, YTOP, YBOT
      INTEGER I, J, K, IZERO, INDX

      DATA MARK/'.','+','X','0','1','2','3','4','5','6','7','8','9','A'/

C.... FIND THE RANGE OF THE DATA.

      XMAX = REAL(DAT(1))
      XMIN = XMAX
      YMAX = AIMAG(DAT(1))
      YMIN = YMAX
      DO I=1,NDAT
        XMAX = AMAX1(XMAX, REAL(DAT(I)))
        XMIN = AMIN1(XMIN, REAL(DAT(I)))
        YMAX = AMAX1(YMAX, AIMAG(DAT(I)))
        YMIN = AMIN1(YMIN, AIMAG(DAT(I)))
      END DO

C.... COMPUTE THE X AND Y SPANS.

      XSPAN = AMAX1(1.0E-9, XMAX-XMIN)
      YSPAN = AMAX1(1.0E-9, YMAX-YMIN)
      SPAN = AMAX1(XSPAN, YSPAN)

C.... COMMENCE WORKING LINE BY LINE.

      DO J=1,NY
        YLABL = YMAX  - SPAN*FLOAT(J-1)/FLOAT(NY-1)
        YTOP  = YLABL + SPAN/2.0/FLOAT(NY-1)
        YBOT  = YLABL - SPAN/2.0/FLOAT(NY-1)
        DO I=1,NX
          LINE(I) = MARK(1)
          IF(YTOP.GE.0.0 .AND. YBOT.LE.0.0) LINE(I) = MARK(2)
        END DO
        XZERO = AMIN1(2.0, AMAX1(-2.0, XMIN/SPAN))
        IZERO = IFIX(1.5-FLOAT(NX-1) * XZERO)
        IF(IZERO.GE.1 .AND. IZERO.LE.NX) LINE(IZERO) = MARK(2)
C...... NOW RUN THROUGH DATA AND FILL IN POINTS.
        DO K=1,NDAT
          Y = AIMAG(DAT(K))
          IF(.NOT. (Y.LE.YBOT .OR. Y.GT.YTOP) ) THEN
            X = REAL(DAT(K))
            I = IFIX(1.5 + FLOAT(NX-1)*(X-XMIN)/SPAN)
            INDX = 4 + MOD(K/10, 11)
            IF(I.GE.1 .AND. I.LE.NX) LINE(I) = MARK(INDX)
          END IF
        END DO
        WRITE(LU,'(3X,F12.4,2X,101A1)') YLABL, (LINE(I),I=1,NX)
      END DO

      XMAX = XMIN + SPAN

      IF(NX.LE.70) THEN
        WRITE(LU,'(10X,F12.4,9X,F7.1,'' TO '',F7.1,'' HZ'',4X,F12.4)')
     &    XMIN, X1, X2, XMAX
      ELSE IF(NX.GT.70) THEN
        WRITE(LU,'(10X,F12.4,19X,F7.1,'' TO '',F7.1,'' HZ'',15X,F12.4)')
     &    XMIN, X1, X2, XMAX
      END IF

      RETURN
      END
