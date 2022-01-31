C     3-DIMENSIONAL WCA FOR A BOUNDARY DRIVEN SHEARING SYSTEM OF
C     PARITICLES OR MOLECULES  
C     
C     Last changes: OCTOBER 2010
C     Based on program kawascii
C
C     NEMD transient trajectory simulation 
C     with a constant colour field applied to the waslls after equilibration
C     Colour field equations - 4th order Runge-Kutta
C     Peculiar internal energy or kinetic energy
c     thermostatted with Gaussian thermostat
C     Periodic boundary conditions in the x & y direction
C     WCA potential
C     FENE potential for intramolecular forces
C     Included NCELL AND ABILITY TO ALTER CUBE LENGTH IN Z DIRECTION
C
C    (i) Equilibration is carried out
C    (ii) A constant colour field to the walls is then applied and 
C         various properties
C        and tcf calculated. Some tau-averaged properties are
C        stored to disk for tau=MAXTAU
C    (ii) The time-correlation function is integrated (if NTYPE.gt.1).
C
C     The accumulated properties from are stored in MEAN 
C   and the number of data accumulated in NMEAN 
C     The accumulated autocorrelation functions are stored in TCF(I,0..NTCF)
C   and the number of data accumulated in KOUNT 
C     The integrated autocorrelation function is stored in INTCF(I,0..NTCF) 
C     The integral of properties along a trajectory is stored in TAUAV(0..NTCF)
C     The average of properties at points on a trajectory is stored in
C   AVER(I,0..NTCF)
C
C-----------------------------------------------------------------------
C
C      MODULE TTCF 

C      CONTAINS 
      
      SUBROUTINE SETUP()

      INCLUDE "SSGK.inc"

C
C***** LOCAL VARIABLES
C
      REAL (KIND=prec) TEST,R6,R12,R4,RSI,TOTPX,TOTPY,TOTPZ
      REAL (KIND=prec) AR1PX,AR1PY,AR1PZ,XX,XX2,RX,RXEQ
      REAL (KIND=prec) AR2PX,AR2PY,AR2PZ,TESTL,TESTH
      REAL (KIND=prec) KENER
      INTEGER I,J,K,M,N,IXX
      INTEGER RANK, NUM_PROCESSORS, IERR
      CHARACTER*9 DATNOW
      REAL (KIND=prec) TIMNOW,SECNDS,RZERO

C
C***** DISPLAY TITLE
C
C      WRITE(6,'(//''NEMD WITH COLOUR: RK4,LJ'')')
C      CALL DATE(DATNOW)
      RZERO=0.0E0
C      TIMNOW=SECNDS(RZERO)
C      WRITE(6,*) DATNOW,TIMNOW
      CONSER=0.0D0
      E00=0.0D0
      KB=0
      TOTPX=0.0D0
      TOTPY=0.0D0
      TOTPZ=0.0D0
      AR1PX=0.0D0
      AR1PY=0.0D0
      AR1PZ=0.0D0
      AR2PX=0.0D0
      AR2PY=0.0D0
      AR2PZ=0.0D0

      FIELD=0.0D0

      DO_PRESSURE = .FALSE.
C
C***** ANALYSE INPUT
C
       IF (MAXTAU.GT.NK) THEN
          WRITE(6,*)'MAXTAU TOO HIGH FOR DIMENSION:',NK
          STOP
       END IF
       IF (NTYPE.LT.1.OR.NTYPE.GT.3) THEN
          WRITE(*,*) NTYPE
          WRITE(6,*)'NO SUCH NTYPE IN THIS PROGRAM, CHOOSE:1,2,3'
          STOP
       END IF
C
       TEST=DBLE(NPART)
       IF(NTYPE.EQ.1) THEN
        IF(LATT.EQ.1)THEN
            IF (MOD(TEST,4.0D0).NE.0.0D0) THEN
               WRITE(6,*)'NPART,MUST BE DIVISIBLE BY 4'
               STOP
            ELSE
               TEST=(NPART/4.0D0*YZDIVX)**(1.0D0/3.0D0)
               TESTL=TEST-0.0000001
               TESTH=TEST+0.0000001
              IF (ANINT(TEST).LT.TESTL.OR.ANINT(TEST).GT.TESTH) THEN
                 WRITE(6,*)'NPART MUST BE 4*N^3/YZDIVX; N INTEGER'
                 WRITE(6,*)TESTL,TESTH
                 STOP
              ENDIF
            ENDIF
        ENDIF 
C C
       TEST=DBLE(NPART*YZDIVX)
         IF(LATT.EQ.2)THEN
           IF (MOD(TEST,2.0D0).NE.0.0D0) THEN
              WRITE(6,*)'NPART,MUST BE AN EVEN NO'
              STOP
           ELSE
              TEST=(NPART*YZDIVX/2.0D0)**(1.0D0/3.0D0)
              TESTL=TEST-0.0000001
              TESTH=TEST+0.0000001

             IF (ANINT(TEST).LT.TESTL.OR.ANINT(TEST).GT.TESTH) THEN
                WRITE(6,*)'NPART MUST BE 2*N^3, WHERE N IS AN INTEGER'
                WRITE(6,*)TESTL,TESTH
                STOP
             END IF
           END IF
          END IF
       END IF

C*****DEFINE NPART2
       NPART2 = NPART/2
C***** ASSIGN MOLECULE NO TO FLUID PARTICLES
    
      DO J=1,NPART2
         IMOL(J)=1
      ENDDO

      DO J=NPART2+1, NPART 
         IMOL(J)=2
      ENDDO
C       
C C***** CHECK ALL FLUID PARTICLES ARE PART OF A MOLECULE
C        IF (MOD(NFLUID,LIMOL).NE.0.0D0) THEN
C           WRITE(6,*)'NFLUID,MUST BE DIVISIBLE BY LIMOL NO REMAINDER'
C           STOP
C        END IF
C
C***** ADD CHARGE TO WALLS PARTICLES ONLY
C
      DO I =1,NPART
         IKOL(I)= (-1)**(I)
      ENDDO
      DO I =1,NPART2
         IMOL(I)= 1
      ENDDO
       DO I =NPART2+1,NPART
         IMOL(I)= -1
      ENDDO

      IF (NGAUS.EQ.6) THEN
         DF=3.0D0*DBLE(NPART)-3.0D0
      ELSE 
         DF=3.0D0*DBLE(NPART)-4.0D0 
      ENDIF
C
C***** ZERO ARRAYS
C
      DENSAV1 = 0.0D0
      DENSAV2 = 0.0D0
      TCFL = 0.0D0
      TCFLL = 0.0D0
      INTCFL = 0.0D0
      INTCFLL = 0.0D0
      KOUNT = 0
      AVER = 0.0D0
      TCF = 0.0D0
      INTCF = 0.0D0
      TCF1 = 0.0D0
      INTCF1 = 0.0D0
      NMEAN = 0
      MEAN = 0.0D0
      M0 = 0.0D0
      M1 = 0.0D0
      M2 = 0.0D0
      M3 = 0.0D0
C
C***** SET INITIAL VALUES
C      RCUT SET TO 2**(1/6) FOR WCA POTENTIAL
C      SEED set as unix time 
        ISEED = 10101
C       ISEED = TIME()
#ifdef DEBUG
       WRITE(6,*) 'SEED', ISEED  
#endif
C
C***** CALCULATE RUN PARAMETERS
C
      IF (NTYPE.EQ.1) THEN
         VOL = NPART/DRF
      ELSE 
         VOL = NPART/DRF
      END IF
C      write(6,*)'VOLF,VOL,DRW,DRF',VOLF,VOL,DRW,DRF
C***** CALCULATE LENGTHS BASED ON YZDIVX RATIO 
      CUBEZ = (YZDIVX*VOL)**(1.0D0/3.0D0)
      CUBEY = CUBEZ
      CUBEX = CUBEZ/YZDIVX
      CUBEX2 = CUBEX/2.0D0
      CUBEY2 = CUBEY/2.0D0                                     
      CUBEZ2 = CUBEZ/2.0D0
      BINS = FLOOR(CUBEX)
      TAUNBINS = FLOOR(CUBEX)
C
C****** IF RCUT TOO LARGE RESET TO HALF BOXLENGTH
C      MIX applies mixing rules  based on EPS1 and EPS2 
      IF (RCUT.GT.CUBEY2) THEN
         RCUT = CUBEX2
         WRITE(6,*)'RCUT RESET TO HALF BOXLENGTH',RCUT
      END IF
      RMAX  = RCUT**2
      RSI    = 1.D0/RMAX
      R4     = RSI*RSI
      R6     = R4*RSI
      R12    = R6*R6
      SHIFT  = R12 - 2.0D0*R6 + EPS_LJG*EXP(((-1.0D0*(SQRT(1.0D0/RSI) - R_LJG)**2.0D0)/(2.0D0*0.02D0)))
C
C***** WRITE SIMULATION PARAMETERS
C
#ifdef DEBUG
      WRITE(6,150) NPART
      WRITE(6,151) TR,E0,DRW,DRF,RCUT,FE0,NTYPE,
     &             DELTA,NON,RCUT,NGAUS,CUBE,MAXTAU,
     &             EQTIM,NCYC,NLAYER,NLP,KH,LATT,NPRINT,
     &             KF,R0,LIMOL,YZDIVX,CUBEX,CUBEY,CUBEZ,
     &             DXXDIV
  150 FORMAT(/,10X,'SIMULATION PARAMETERS',I6,' PARTICLES')
  151 FORMAT(/,5X,'TR      =',F14.7,10X,'E/N     =',F14.7,10X,
     &       /,5X,'WALL DENSITY =',F14.7,10X,'FL DEN=',F14.7,10X,
     &       /,5X,'RCUT    =',F14.7,
     &       /,5X,'FE0   =',F14.7,10X,'NTYPE  =',I6,
     &       /,5X,'DELTA  =',F14.7,10X,'NON    =',I6,
     &       /,5X,'RCUT   =',F14.7,10X,'NGAUS  =',I6,
     &       /,5X,'CUBE   =',F14.7,10X,'MAXTAU =',I6,
     &       /,5X,'EQTIM  =',I6,10X,'NCYC    =',I6, 
     &       /,5X,'NLAYER =',F14.7,10X,'NLP =',I6,
     &       /,5X,'KH     =',F14.7,10X,'LATT =',I6,
     &       /,5X,'NPRINT =',I6
     &       /,5X,'KF     =',F14.7,10X,'R0  =',F14.7,
     &       /,5X,'LIMOL =',I6
     &       /,5X,'YZDIVX =',F14.7,10X,'CUBEX =',F14.7,
     &       /,5X,'CUBEY =',F14.7,10X,'CUBEZ =',F14.7,
     &       /,5X,'DXXDIV =',F14.7)
#endif

C
C****************************************************************
C
C***** INITIAL STARTUP
C
C
C***** NTYPE=1  INITIAL FROM FCC
C
         CALL FCC
C Store initial particle positions for mean-square displacement calculations
         X0 = X
         Y0 = Y
         Z0 = Z
C
C***** CHECK ALL FLUID PARTICLES ARE PART OF A MOLECULE
C        IF (MOD(NFLUID,LIMOL).NE.0.0D0) THEN
C           WRITE(6,*)'NPART DIFFERENT FROM INITIAL INPUT, NFLUID',
C      &               ' MUST BE DIVISIBLE BY LIMOL NO REMAINDER'
C           STOP
C        END IF
C
C***** UPDATE INPUT UNIT
C
C
#ifdef DEBUG
      WRITE(6,*) TR,DRW,DRF,DELTA,LATT,NPART,NLP
      WRITE(6,*) FE0,RCUT,NLAYER,KH,NPRINT
      WRITE(6,*) MIX, EPS1, EPS2, QVOL
      WRITE(6,*) KF,R0,LIMOL,YZDIVX
      WRITE(6,*) DXXDIV
      WRITE(6,*) NTYPE,NON,NGAUS,E0
      WRITE(6,*) NPLOT,MAXTAU,EQTIM,NCYC
      WRITE(6,*) 'TR,DRW,DRF,DELTA,LATT,NPART,NLP'
      WRITE(6,*) 'FE0,RCUT,NLAYER,KH,NPRINT'
      WRITE(6,*) 'KF,R0,LIMOL,YZDIVX'
      WRITE(6,*) 'MIX, EPS1, EPS2, QVOL'
      WRITE(6,*) 'DXXDIV'
      WRITE(6,*) 'NTYPE,NON,NGAUS,E0'
      WRITE(6,*) 'NPLOT,MAXTAU,EQTIM,NCYC'

C***** CHECK IMOL IKOL AND SWITCHES
        WRITE(6,*)'I,IMOL(I),IKOL(I),S1(I),S2(I),S3(I),S4(I)'
        DO I=1,NPART
         WRITE(6,*)I,IMOL(I),IKOL(I),S1(I),S2(I),S3(I),S4(I)
        ENDDO
#endif
C
C CALCULATE PROPERTIES AT TIME ZERO
C
         KE2=0.0D0
         DO 565 I=1, NPART
            KE2 = KE2+(PX(I)**2+PY(I)**2+PZ(I)**2)
 565     CONTINUE
c         write(6,*)'force 1'
         CALL FORCECELL 
         E00=UPOT+KE2/2.0D0
         FIELD=0.0D0

C***** UPDATE INITIAL TEMPERATURE
         TEMP = (KE2/2.0D0)/(DF/2.0D0)
       
 45   FORMAT(/,2X,'EQUILIBRATION - INSTANTANEOUS VALUES',/,8X,
     &  'STEP    TEMP    UPOT/NP    TOTE/NP  ALPHA')
 75   FORMAT(/,2X,'EQUILIBRATION - INSTANTANEOUS VALUES',/,8X,
     &  'STEP    TEMP    UPOT/NP    TOTE/NP  ALPHA')
 55   FORMAT(/,2X,'EQUILIBRATION - INSTANTANEOUS VALUES',/,8X,
     &  'STEP    TEMP    UPOT/NP    TOTE/NP  ALPHA')
 85    FORMAT(/,2X,'EQUILIBRATION - INSTANTANEOUS VALUES',/,8X,
     &  'STEP    TEMP    UPOT/NP    TOTE/NP  ALPHA')
 65   FORMAT(6E14.6)
      END SUBROUTINE


C
C**********************************************************
C
      SUBROUTINE MD(KPROP,IFLAG)
C***** DO COLOUR DIFFUSION MD USING RK4 FOR KPROP STEPS
C***** ACCUMULATE VALUE OF VARIOUS PROPERTIES MEAN AND
C***** ACCUMULATE TCF IN TCF 
C
Cf2py threadsafe
      INCLUDE "SSGK.inc"
C
C***** LOCAL VARIABLES
C
      INTEGER I,J,ISTEP,IX,IY,IZ,IFLAG,KPROP, BIN, NUM,LIM1,LIM2,TAUBIN 
      REAL (KIND=prec) KENER,KTRAN,KENERW,AHEAT
      REAL (KIND=prec) PP,PPF,PXY,PXZ,PYZ,JX,JY,JZ,JXL,JYL,JZL
      REAL (KIND=prec) VX,VY,VZ,VXL,VYL,VZL
      REAL (KIND=prec) VVX,VVY,VVZ,VVXL,VVYL,VVZL,VVXLL,VVYLL,
     &                 VVZLL,XR(NP),YR(NP),ZR(NP)
      REAL (KIND=prec) PX3,PXYF,PXZF,PYZF,P2,NTOTE
      INTEGER TMAX,K,NINBIN,NINBINL,NINBIN0L
C
      MSDX = 0.0D0
      MSDY = 0.0D0
      MSDZ = 0.0D0
C Turn on the field as appropriate
      IF (IFLAG .EQ. 0) THEN
        FIELD = 0.0D0
      ELSE
        FIELD = FE0
      ENDIF

C
C CALCULATE PROPERTIES AT TIME ZERO
C
      KE2=0.0D0
      DO 811 I=1,NPART
         KE2=KE2+(PX(I)**2+PY(I)**2+PZ(I)**2)
 811  CONTINUE
      CALL FORCECELL
      E00=UPOT+KE2/2.0D0
      KENER = KE2/2.0D0

C       DO I=0,TAUNBINS 
C         TAUDENS1(I) = 0.0D0
C         TAUDENS2(I) = 0.0D0
C       ENDDO 

C***** KINETIC PART OF AVERAGES
C     
      KTRAN = 0.0D0
      KENER = 0.0D0
      KENERW = 0.0D0

      DO 81 I = 1, NPART
            S0L(I) = 0.0D0
            SL(I) = 0.0D0
            KTRAN = KTRAN + PY(I)**2
            KENER = KENER + (PX(I)**2+PY(I)**2+PZ(I)**2)/2.D0
            T0V(I,1) = PX(I)   
            T0V(I,2) = PY(I)
            T0V(I,3) = PZ(I)

            IF (FLOOR(X(I)).LT.1.0D0) THEN 
                  S0L(I) = 1.0D0 
                  SL(I) = 1.0D0
            ENDIF  

            PX3 = PX3+PX(I)**3

 81   CONTINUE

         KTRAN = KTRAN/(DF/2.0D0)
         TEMP  = KENER/(DF/2.0D0)
         TOTE  = KENER + UPOT
         TOTE  = TOTE/NPART
C          DO I=0, BINS
C             DISS0L(I) = DISS0L(I) + (JXL(I)*FIELD)*TR
C          ENDDO 


      DO ISTEP=1,KPROP
         KE2=0.0D0
         DO 11 I=1,NPART
            KE2=KE2+(PX(I)**2+PY(I)**2+PZ(I)**2)
 11      CONTINUE
         E00=UPOT+KE2*0.5D0
c         IF (NON.EQ.1.OR.((NGAUS.EQ.2).AND.(IFLAG.EQ.0))) THEN
            CALL FORCECELL
c         ENDIF

         CALL DOMD
C        
C THE FOLLOWING CALL FORCECELL IS REQUIRED DUE TO PBCS
         CALL FORCECELL
C
C***** RESCALE TEMP/ENERGY IF NON.EQ.1 OR DURING EQUILIBRATION
C

         KENER  = 0.0D0
         DO 3 I = 1, NPART
            KENER = KENER + (PX(I)**2 + PY(I)**2
     &               +PZ(I)**2)/2.0D0
 3       CONTINUE

C THERMOSTAT CORRECTIONS
         IF (MOD(NGAUS,2).EQ.1) THEN
             CONSER=MAX(CONSER,ABS(E0*NPART-KENER-UPOT))
             AHEAT   = SQRT((E0*NPART-UPOT)/KENER)
          IF (NON.EQ.1.OR.((NGAUS.EQ.1).AND.(IFLAG.EQ.0)))
     &        THEN
                DO 4 I = 1, NPART
                   PX(I) = PX(I)*AHEAT
                   PY(I) = PY(I)*AHEAT
                   PZ(I) = PZ(I)*AHEAT
 4              CONTINUE
            END IF               ! NON.EQ.1....
C*****************************************
         ELSE IF (MOD(NGAUS,2).EQ.0) THEN
            CONSER=MAX(CONSER,ABS((TR*(DF))-(2*KENER)))
            AHEAT   = SQRT((TR*(DF))/(2*KENER))
         IF (NON.EQ.1.OR.((NGAUS.EQ.2).AND.(IFLAG.EQ.0))) THEN
               DO 5 I = 1, NPART
                  PX(I) = PX(I)*AHEAT
                  PY(I) = PY(I)*AHEAT
                  PZ(I) = PZ(I)*AHEAT
    5          CONTINUE
            ENDIF
         ENDIF

C
C***** CALCULATE AVERAGES 
C
            KB=KB+1
C
C***** KINETIC PART OF AVERAGES
          DO 1 I = 1, NPART
               SL(I) = 0.0D0
               KTRAN = KTRAN + PY(I)**2
               KENER = KENER + (PX(I)**2+PY(I)**2+PZ(I)**2)/2.D0

              IF (FLOOR(X(I)).LT.1.0D0) THEN 
                  SL(I) = 1.0D0
              ENDIF 

             MSDX = MSDX + (X(I) - X0(I))**2
             MSDY = MSDX + (Y(I) - Y0(I))**2
             MSDZ = MSDX + (Z(I) - Z0(I))**2
 1          CONTINUE
            MSDX = MSDX / (REAL(NPART, prec))
            MSDY = MSDY / (REAL(NPART, prec))
            MSDZ = MSDZ / (REAL(NPART, prec))
C*********** Pressure
         PXY   = 0.5D0*(PT(1,2)+PT(2,1))
         PXZ   = 0.5D0*(PT(1,3)+PT(3,1))
         PYZ   = 0.5D0*(PT(2,3)+PT(3,2))
         PXYF   = 0.5D0*(PTF(1,2)+PTF(2,1))
         PXZF   = 0.5D0*(PTF(1,3)+PTF(3,1))
         PYZF   = 0.5D0*(PTF(2,3)+PTF(3,2))
         PP    = 0.5D0*(PT(1,1)+PT(2,2)+PT(3,3))/3
         PPF   = 0.5D0*(PTF(1,1)+PTF(2,2)+PTF(3,3))/3
C********** ISTEP
      END DO
C
 55   FORMAT(I12,7F9.5)
 75   FORMAT(I12,7F9.5)
 65   FORMAT(6E14.6)
      

C
      RETURN
      END
C
C**********************************************************             
C
      SUBROUTINE DOMD
C
C***** DO A SINGLE STEP OF A SIMULATION
C      USE COLOUR DIFFUSION ALGORITHM,FOURTH ORDER RUNGE-KUTTA METHOD,
C      PERIODIC BOUNDARY CONDITIONS AND
C      GAUSSIAN THERMOSTAT OR ERGOSTAT
C*****
C
      INCLUDE "SSGK.inc"
C
C***** LOCAL VARIABLES
C
      INTEGER I
      REAL (KIND=prec) K1X(NP),K1Y(NP),K1Z(NP),K1XEQ(NP),
     &                 K1PX(NP),K1PY(NP),K1PZ(NP),
     &                 K2X(NP),K2Y(NP),K2Z(NP),K2XEQ(NP),
     &                 K2PX(NP),K2PY(NP),K2PZ(NP),            
     &                 K3X(NP),K3Y(NP),K3Z(NP),K3XEQ(NP),
     &                 K3PX(NP),K3PY(NP),K3PZ(NP),
     &                 K4X(NP),K4Y(NP),K4Z(NP),K4XEQ(NP),
     &                 K4PX(NP),K4PY(NP),K4PZ(NP),
     &                 K1NH,K2NH,K3NH,K4NH, KENER          
      REAL (KIND=prec) XOLD(NP),YOLD(NP),ZOLD(NP),XEQOLD(NP),
     &                 PXOLD(NP),PYOLD(NP),PZOLD(NP) 
      REAL (KIND=prec) NHOLD
      REAL (KIND=prec) RX,RY,RZ,RXEQ,RYEQ,CIS,SUMRAN
      REAL (KIND=prec) SUMX,SUMY,SUMZ
C
      SUMRAN=0.0D0
      KENER =0.0D0
      NHOLD = NHALPH 
C
C***** 4TH ORDER RK DE SOLVER
C
      DO 10 I = 1, NPART                                                   
       XOLD(I)  = X(I)                                                  
       YOLD(I)  = Y(I)
       ZOLD(I)  = Z(I)                                                   
       PXOLD(I) = PX(I)                                                 
       PYOLD(I) = PY(I)
       PZOLD(I) = PZ(I)
       XEQOLD(I)= XEQ(I)   
   10 CONTINUE                          
           
C                                                                       
C***** K1                                                                
C                                                                       
      DO 15 I = 1, NPART                                                   
       K1X(I)  = PX(I) 
       K1Y(I)  = PY(I)
       K1Z(I)  = PZ(I)                                                  
       K1PX(I) = FX(I) - ALPH*PX(I) + IKOL(I)*FIELD       
       K1PY(I) = FY(I) - ALPH*PY(I)
       K1PZ(I) = FZ(I) - ALPH*PZ(I)
       KENER = KENER + (PX(I)**2+PY(I)**2+PZ(I)**2)                
   15 CONTINUE
      K1NH = (KENER-((DF)/TR))/QVOL
C                                                                       
C***** K2                                                                
C                                                                       
      DO 20 I = 1,NPART                                                    
       X(I)  = XOLD(I)  + K1X(I) *DELTA/2.D0                            
       Y(I)  = YOLD(I)  + K1Y(I) *DELTA/2.D0
       Z(I)  = ZOLD(I)  + K1Z(I) *DELTA/2.D0                             
       PX(I) = PXOLD(I) + K1PX(I)*DELTA/2.D0                            
       PY(I) = PYOLD(I) + K1PY(I)*DELTA/2.D0
       PZ(I) = PZOLD(I) + K1PZ(I)*DELTA/2.D0
       XEQ(I)= XEQOLD(I)+ K1XEQ(I)*DELTA/2.D0   

   20 CONTINUE

      NHALPH = NHOLD + K1NH*DELTA/2.0D0

C
C***** GET NEW FORCECELLS,ALPHA
C
      CALL FORCECELL
      KENER=0.0D0                                                        
C                                                                       
      DO 25 I = 1,NPART                                                    
       K2X(I)  = PX(I)
       K2Y(I)  = PY(I)
       K2Z(I)  = PZ(I)                                                  
       K2PX(I) = FX(I) - ALPH*PX(I) + IKOL(I)*FIELD       
       K2PY(I) = FY(I) - ALPH*PY(I)
       K2PZ(I) = FZ(I) - ALPH*PZ(I)  
       KENER = KENER + (PX(I)**2+PY(I)**2+PZ(I)**2)                 
   25 CONTINUE   

      K2NH = (KENER-((DF)/TR))/QVOL                                               
C                                                                       
C***** K3                                                                
C                                                                       
      DO 30 I = 1,NPART                                                    
       X(I)  = XOLD(I)  + K2X(I)  *DELTA/2.D0                           
       Y(I)  = YOLD(I)  + K2Y(I)  *DELTA/2.D0
       Z(I)  = ZOLD(I)  + K2Z(I)  *DELTA/2.D0                           
       PX(I) = PXOLD(I) + K2PX(I) *DELTA/2.D0                           
       PY(I) = PYOLD(I) + K2PY(I) *DELTA/2.D0
       PZ(I) = PZOLD(I) + K2PZ(I) *DELTA/2.D0 
       XEQ(I)= XEQOLD(I)+ K2XEQ(I)*DELTA/2.D0                          
   30 CONTINUE

      NHALPH = NHOLD + K2NH*DELTA/2.0D0     

C
C***** GET NEW FORCECELLS, ALPH
C
      CALL FORCECELL 
      KENER=0.0D0                                                       
C                                                                       
      DO 35 I = 1,NPART                                                    
       K3X(I)  = PX(I)
       K3Y(I)  = PY(I)
       K3Z(I)  = PZ(I)                                                  
       K3PX(I) = FX(I) - ALPH*PX(I) + IKOL(I)*FIELD       
       K3PY(I) = FY(I) - ALPH*PY(I)
       K3PZ(I) = FZ(I) - ALPH*PZ(I)
       KENER = KENER + (PX(I)**2+PY(I)**2+PZ(I)**2)                 
   35 CONTINUE
      
      K3NH = (KENER-((DF)/TR))/QVOL                                                               
C***** K4                                                                
C                                                                       
      DO 40 I = 1, NPART                                                   
       X(I)  = XOLD(I)  + K3X(I)*DELTA                                  
       Y(I)  = YOLD(I)  + K3Y(I)*DELTA
       Z(I)  = ZOLD(I)  + K3Z(I)*DELTA                                  
       PX(I) = PXOLD(I) + K3PX(I)*DELTA                                 
       PY(I) = PYOLD(I) + K3PY(I)*DELTA  
       PZ(I) = PZOLD(I) + K3PZ(I)*DELTA
       XEQ(I)= XEQOLD(I)+ K3XEQ(I)*DELTA    

   40 CONTINUE  

      NHALPH = NHOLD + K3NH*DELTA                                                      
C                                                                       
C
C***** GET NEW FORCECELLS,ALPH
C
      CALL FORCECELL  
      KENER=0.0D0                                                      
C                                                                       
      DO 45 I = 1, NPART                                                   
       K4X(I)  = PX(I) 
       K4Y(I)  = PY(I)
       K4Z(I)  = PZ(I)                                                  
       K4PX(I) = FX(I) - ALPH*PX(I) + IKOL(I)*FIELD       
       K4PY(I) = FY(I) - ALPH*PY(I)
       K4PZ(I) = FZ(I) - ALPH*PZ(I) 
       KENER = KENER + (PX(I)**2+PY(I)**2+PZ(I)**2)                     
   45 CONTINUE
      K4NH= (KENER-((DF)/TR))/QVOL
C                                                                        
      DO 50 I = 1, NPART                                                   
       X(I)  = XOLD(I)                                                  
     &      + DELTA/6.D0*(K1X(I) + 2.D0*K2X(I) + 2.D0*K3X(I) + K4X(I))  
       Y(I)  = YOLD(I)                                                  
     &      + DELTA/6.D0*(K1Y(I) + 2.D0*K2Y(I) + 2.D0*K3Y(I) + K4Y(I))
       Z(I)  = ZOLD(I)
     &      + DELTA/6.D0*(K1Z(I) + 2.D0*K2Z(I) + 2.D0*K3Z(I) + K4Z(I))  
       PX(I) = PXOLD(I)                                                 
     &      + DELTA/6.D0*(K1PX(I)+ 2.D0*K2PX(I)+ 2.D0*K3PX(I)+ K4PX(I)) 
       PY(I) = PYOLD(I)                                                 
     &      + DELTA/6.D0*(K1PY(I)+ 2.D0*K2PY(I)+ 2.D0*K3PY(I)+ K4PY(I)) 
       PZ(I) = PZOLD(I)
     &      + DELTA/6.D0*(K1PZ(I)+ 2.D0*K2PZ(I)+ 2.D0*K3PZ(I)+ K4PZ(I)) 
       XEQ(I) = XEQOLD(I)
     &      +DELTA/6.D0*(K1XEQ(I)+2.D0*K2XEQ(I)+2.D0*K3XEQ(I)+K4XEQ(I)) 
   50 CONTINUE

      NHALPH = NHOLD                                                 
     &      + DELTA/6.D0*(K1NH+ 2.D0*K2NH+ 2.D0*K3NH+ K4NH)   
C***** PERIODIC BOUNDARY CONDITIONS                                      
C                                                                       
      DO 60 I = 1, NPART                                                   
         RX =   X(I) - CUBEX2
         RY =   Y(I) - CUBEY2
         RZ =   Z(I) - CUBEZ2
         RX   = RX - CUBEX*ANINT(RX/CUBEX)
         RY   = RY - CUBEY*ANINT(RY/CUBEY)
         RZ   = RZ - CUBEZ*ANINT(RZ/CUBEZ)
         X(I) = RX + CUBEX2
         Y(I) = RY + CUBEY2
         Z(I) = RZ + CUBEZ2

C          IF (X(I).GT.CUBEX) THEN
C             WRITE(6,*) '1',KB,X(I),Y(I),Z(I),PX(I),PY(I),PZ(I)
C             STOP
C          END IF
C          IF (Y(I).GT.CUBEY) THEN
C             WRITE(6,*) '2',KB,X(I),Y(I),Z(I),PX(I),PY(I),PZ(I)
C             STOP
C          END IF
C          IF (Z(I).GT.CUBEZ) THEN
C             WRITE(6,*) '3',KB,I,X(I),Y(I),Z(I),PX(I),PY(I),PZ(I)
C             STOP
C          END IF
C          IF (X(I).LT.0.0D0) THEN
C             WRITE(6,*) '4',KB,X(I),Y(I),Z(I),PX(I),PY(I),PZ(I)
C             STOP
C          END IF
C          IF (Y(I).LT.0.0D0) THEN
C             WRITE(6,*) '5',KB,X(I),Y(I),Z(I),PX(I),PY(I),PZ(I)
C             STOP
C          END IF
C           IF (Z(I).LT.0.0D0) THEN
C             WRITE(6,*) '6',KB,I,X(I),Y(I),Z(I),PX(I),PY(I),PZ(I)
C             STOP
C          END IF
   60 CONTINUE                                                      
C                                                                       
      RETURN                                                            
      END
C
C   ********************************************************
C
      SUBROUTINE FCC
      INCLUDE "SSGK.inc"
C
C***** LOCAL VARIABLES
C
      INTEGER NX,NY,NZ,M,J,K,IJ,I
      REAL (KIND=prec) SUMX,SUMY,SUMZ,SUMX2,SUMY2,SUMZ2
      REAL (KIND=prec) RTR,XX,YY,ZZ,XYZ,PS,PS2,HEAT,HEAT2
      REAL (KIND=prec) DX,DY,DZ,DNY,CUBEXY

C
C***** CALCULATE NUMBER OF TWO PARTICLE CELLS IN X,Y AND Z DIRECTION
C

      IF(LATT.EQ.1)THEN
        DNY=DFLOAT(NPART)
        DNY=DNY/4.0D0
        DNY = DNY*YZDIVX
        DNY=DNY**(1.0D0/3.0D0)
        NX = NINT(DNY/YZDIVX)
        NY = NINT(DNY)
        NZ = NY
C
C***** DX,DY,DZ ARE THE SIZE OF EACH SIDE OF A SINGLE TWO-PARTICLE
C***** CELL ASSUMING BOXLENGTH IS 1
C
        DX = 1.D0/NX
        DY = 1.D0/NY
        DZ = 1.D0/NZ
C
C***** SET POSITION OF THE FOUR PARTICLES FOR FCC(LATT=1) IN FIRST CELL
C
        X(1) = 0.25D0*DX
        Y(1) = 0.25D0*DY
        Z(1) = 0.25D0*DZ
        X(2) = 0.75D0*DX
        Y(2) = 0.75D0*DY
        Z(2) = 0.25D0*DZ
        X(3) = 0.25D0*DX
        Y(3) = 0.75D0*DY
        Z(3) = 0.75D0*DZ
        X(4) = 0.75D0*DX
        Y(4) = 0.25D0*DY
        Z(4) = 0.75D0*DZ
C
C***** SET POSITION OF ALL OTHER PARTICLES
C
         M = 0
          DO 1  K = 1, NX
          DO 1  J = 1, NY
          DO 1  I = 1, NZ
          DO 2 IJ = 1, 4
             X(IJ + M) = X(IJ) + DX*(K-1)
             Y(IJ + M) = Y(IJ) + DY*(J-1)
             Z(IJ + M) = Z(IJ) + DZ*(I-1)
   2      CONTINUE
           M = M + 4
   1     CONTINUE
      ENDIF
C
C***** SET POSITIONS FOR  PARTICLES IN  BCC(LATT=2)
C
      IF(LATT.EQ.2)THEN
         DNY=DFLOAT(NPART)
         DNY=DNY/2.0D0
         DNY = DNY*YZDIVX
         DNY=DNY**(1.0D0/3.0D0)
         NX = NINT(DNY/YZDIVX)
         NY = NINT(DNY)
         NZ = NY
         NL = NX

C*****CHECKING
#ifdef DEBUG
         WRITE(6,*)'DNY,NX,NY,NZ,NL',DNY,NX,NY,NZ,NL
         WRITE(6,*)'CUBEZ/CUBEX',CUBEZ/CUBEX
#endif
C
C***** DX,DY,DZ ARE THE SIZE OF EACH SIDE OF A SINGLE TWO-PARTICLE
C***** CELL ASSUMING BOXLENGTH IS 1
C
      DX = 1.0D0/NX
      DY = 1.0D0/NY
      DZ = 1.0D0/NZ
C
C***** SET POSITION OF THE  PARTICLES FOR BCC(LATT=2) IN FIRST CELL
C
         X(1) = 0.25D0*DX
         Y(1) = 0.25D0*DY
         Z(1) = 0.25D0*DZ
         X(2) = 0.75D0*DX
         Y(2) = 0.75D0*DY
         Z(2) = 0.75D0*DZ
C
C***** SET POSITION OF ALL OTHER PARTICLES
C
         M = 0
         DO 11  I = 1, NZ
         DO 11  J = 1, NY
         DO 11  K = 1, NX
          DO 12 IJ = 1, 2
           X(IJ + M) = X(IJ) + DX*(K-1)
           Y(IJ + M) = Y(IJ) + DY*(J-1)
           Z(IJ + M) = Z(IJ) + DZ*(I-1)
   12     CONTINUE
           M = M + 2
   11    CONTINUE
      ENDIF
C
C*****ETCH OUT FLUID SITES BY REDEFIENING NPART
C*****check
#ifdef DEBUG
       WRITE(6,*)'NPART,NPART,DRW,DRF'
       WRITE(6,*)NPART,NPART,DRW,DRF
#endif
       NPART2  = NPART/2
C        VOLF = NPART*DRF
C        NPART=NPART+INT(DRF/DRW*NFLUID/2.0D0)*2
C C       NFLUID=NPART-NPART
C        DRF=NFLUID/VOLF

C*******CHECK
#ifdef DEBUG
       WRITE(6,*)'NPART,NPART, DRW,DRF'
       WRITE(6,*)NPART,NPART,DRW,DRF
#endif

C
C***** REDUCE TO COORDINATES OF THE CUBE
C
      DO 3 I = 1, NPART
       X(I) = X(I)*CUBEX
       Y(I) = Y(I)*CUBEY 
       Z(I) = Z(I)*CUBEZ
    3 CONTINUE
C
C***** DEFINE DXX
C
      DXX=DX*CUBEX
*****CHECKING
C
C***** DEFINE NPART AND SWITCHES
      
         DO I=1, NPART2
               S1(I)=1.0D0
               S2(I)=0.0D0
               S3(I)=1.0D0
               S4(I)=0.0D0
         ENDDO
         DO I=NPART2+1, NPART
               S1(I)=0.0D0
               S2(I)=1.0D0
               S3(I)=0.0D0
               S4(I)=1.0D0
         ENDDO

C***** SET EQILIBRIUM POSITIONS 
C      
         DO I=1, NPART
            IF(I.LE.NPART) THEN 
               XEQ(I)=X(I)
               YEQ(I)=Y(I)
               ZEQ(I)=Z(I)
            ENDIF
          ENDDO

C
C***** WRITE INITIAL POSITIONS OF ATOMS
C
C      WRITE(9,100)
C  100 FORMAT(///,5X,'INITIAL COORDINATES',//)
C     WRITE(9,101) (X(I),Y(I),Z(I),PX(I),PY(I),P(I),I=1,NPART)
C 101 FORMAT(/,2(2X,F10.4))
C
C***** GENERATE RANDOM VELOCITIES WITH A GAUSSIAN DISTRUBUTION
C***** ABOUT SQRT(TR,E0)
C
      SUMX = 0.0D0
      SUMY = 0.0D0
      SUMZ = 0.0D0

         IF (MOD(NGAUS,2).EQ.1) THEN
         RTR  = SQRT(2.0D0*E0)
      ELSE IF (MOD(NGAUS,2).EQ.0) THEN
         RTR  = SQRT(2.0D0*TR)
      END IF
      DO 4 I = 1,NPART,2
       CALL RAN1(ISEED,XX)
       CALL RAN1(ISEED,YY)
       XYZ   = 1.D0/SQRT(XX*XX+YY*YY)
       PX(I) = XX*XYZ*RTR
       PY(I) = YY*XYZ*RTR

       CALL RAN1(ISEED,XX)
       CALL RAN1(ISEED,YY)
       XYZ   = 1.D0/SQRT(XX*XX+YY*YY)
       PZ(I) = XX*XYZ*RTR
       PX(I+1) = YY*XYZ*RTR

       CALL RAN1(ISEED,XX)
       CALL RAN1(ISEED,YY)
       XYZ   = 1.D0/SQRT(XX*XX+YY*YY)
       PY(I+1) = XX*XYZ*RTR
       PZ(I+1) = YY*XYZ*RTR

    4 CONTINUE

      DO I=1,NPART
        SUMX = SUMX +PX(I)
        SUMY = SUMY +PY(I)
        SUMZ = SUMZ +PZ(I)
      ENDDO

C
C***** WALL - CORRECT SO SUM MOMENTA IS ZERO IN EACH DIRECTION 
C***** AND SUM OF SQUARES OF MOMENTA EQUALS 2.0*E0*NPART 
C
      PS     = 0.0D0
      DO 5 I = 1,NPART
         PX(I) = PX(I) - SUMX/NPART
         PY(I) = PY(I) - SUMY/NPART
         PZ(I) = PZ(I) - SUMZ/NPART
  105    CONTINUE
         PS = PS + PX(I)**2 + PY(I)**2+ PZ(I)**2
    5 CONTINUE

      IF (MOD(NGAUS,2).EQ.1) THEN
C***** THIS BIT IS NOT YET PROPERLY IMPLEMENTED
         CALL FORCECELL
         CALL FORCECELL
         HEAT = SQRT(2.0D0*(E0*NPART-UPOT)/PS)
C***********************************************
      ELSE IF (MOD(NGAUS,2).EQ.0) THEN
         HEAT = SQRT(TR*(DF)/(PS))
      END IF
      DO I = 1, NPART
         PX(I) = PX(I)*HEAT
         PY(I) = PY(I)*HEAT
         PZ(I) = PZ(I)*HEAT
      ENDDO

C***** OUTPUT INITIAL POSTIONS TO FILE
#ifdef DEBUG
      WRITE(6,*)'X(I),Y(I),Z(I),XEQ(I),NPART,NL'
      DO I=1,NPART
         WRITE(6,*) X(I),Y(I),Z(I),XEQ(I),NPART,NL,PX(I),
     &               PY(I),PZ(I)
      ENDDO
#endif
      SUMX=0.0D0
      SUMY=0.0D0
      SUMZ=0.0D0
      DO I=1,NPART
         SUMX = SUMX+PX(I)
         SUMY = SUMY+PY(I)
         SUMZ = SUMZ+PZ(I)
      ENDDO 

      NHALPH=0.0D0

C***** OUTPUT INITIAL VELOCITIES TO SCREEN

#ifdef DEBUG
      WRITE(6,103) (PX(I),PY(I),PZ(I),I=1,NPART)

  103 FORMAT(/,2(2X,F10.4))
#endif
      RETURN
      END

C
C   ***********************************************
C
      SUBROUTINE FORCECELL
C
C***** CALCULATES FORCEL ON EACH PARTICLE, TOTAL POTENTIAL ENERGY
C      THERMOSTAT:ALPH, ADN CONFIGURATIONAL PART OF PRESSURE TENSOR
C***** WITH A WCA POTENTIAL
C
      INCLUDE "SSGK.inc"
C
C***** LOCAL VARIABLES
C
      INTEGER I,J,NA,NC
      PARAMETER (NA=1000,NC=60)
      INTEGER NCELLX,NCELLY,NCELLZ
      INTEGER ICX,ICY,ICZ,NIC,NJC
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: NUMCELL
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:) :: ICELL
      INTEGER JCX,JCY,JCZ,JJCX,JJCY,JJCZ
      INTEGER I1,I2
      INTEGER IBIN
      REAL (KIND=prec) RX,RY,RZ,CIS,RSQ,R02,RSI,R4,R6,R12,RSCALAR,EW
      REAL (KIND=prec) FIJX,FIJY,FIJZ,ANUM,ADEN,ANUMB,ADENB
      REAL (KIND=prec) PXPY,PXPZ,PYPZ
      REAL (KIND=prec) KX(NP),KY(NP),KZ(NP)
      REAL (KIND=prec) CUBINVX,CUBINVY,CUBINVZ
      REAL (KIND=prec) RY2,RYC,RX2,RXC,RZ2,RZC,ROUT2
      REAL (KIND=prec) ICXD,ICYD,ICZD,RXCL,RYCL,RZCL,RXCL2
      REAL (KIND=prec) RYCL2,RZCL2,RSQCL,UCL
C
C***** ZERO ARRAYS
C
      FX = 0.0D0
      FY = 0.0D0
      FZ = 0.0D0
      RIJ_HIST = 0.D0
      PT = 0.0D0

C***** This is the number to sort the particle pairs' radial separations into the correct histogram bins
      RBIN_INV = 1.D0 / (RMAX/REAL(NBINS))

C***** DEFINE USEFUL TERMS
C
      R02 = R0**2
      CUBINVX    = 1.0D0/CUBEX
      CUBINVY    = 1.0D0/CUBEY
      CUBINVZ    = 1.0D0/CUBEZ
C
      IF(R0.GT.RCUT)THEN
        NCELLX     = INT(CUBEX/R0)
        NCELLY     = INT(CUBEY/R0)
        NCELLZ     = INT(CUBEZ/R0)
        ROUT2      = R0**2
      ELSE
        NCELLX     = INT(CUBEX/RCUT)
        NCELLY     = INT(CUBEY/RCUT)
        NCELLZ     = INT(CUBEZ/RCUT)
        ROUT2      = RMAX
      ENDIF 

C
C      IF (NCELLX.LT.3) STOP 'Do not use cell code, NCELLX<3'
C      IF (NCELLY.LT.3) STOP 'Do not use cell code, NCELLY<3'
C      IF (NCELLZ.LT.3) STOP 'Do not use cell code, NCELLZ<3'

C
C***** ZERO NUMBER-IN-CELL COUNTER
C     
       
      ALLOCATE(NUMCELL(NCELLX, NCELLY, NCELLZ))
      ALLOCATE(ICELL(NA, NCELLX, NCELLY, NCELLZ))
      NUMCELL = 0
C
C**** SORT PARTICLES INTO CELLS
C
      DO 11 I=1,NPART
            ICXD  =(X(I)*DFLOAT(NCELLX))*CUBINVX
               IF (ICXD.LT.0) THEN
               ICX =INT(ICXD)+NCELLX
               ELSE
               ICX =INT(ICXD)+1
               ENDIF

            ICYD  =(Y(I)*DFLOAT(NCELLY))*CUBINVY
               IF (ICYD.LT.0) THEN
               ICY =INT(ICYD)+NCELLY
               ELSE
               ICY =INT(ICYD)+1
               ENDIF

            ICZD  =(Z(I)*DFLOAT(NCELLZ))*CUBINVZ
               IF (ICZD.LT.0) THEN
               ICZ =INT(ICZD)+NCELLZ
               ELSE
               ICZ =INT(ICZD)+1
               ENDIF
c         write (6,*)'icx,icy,icz',icx,icy,icz
C
C**** MAKE SURE PARTICLES WITHIN NCELL BOXES
C
            IF (ICX.GT.NCELLX) ICX=ICX-NCELLX
            IF (ICY.GT.NCELLY) ICY=ICY-NCELLY
            IF (ICZ.GT.NCELLZ) ICZ=ICZ-NCELLZ

C
C**** NUMCELL(ICX,ICY,ICZ) NUMBER IN CELL 'IC'
C
            NUMCELL(ICX,ICY,ICZ) = NUMCELL(ICX,ICY,ICZ)+1
            NIC = NUMCELL(ICX,ICY,ICZ)

C**** SO THE iTH PATICLE OF A CELL CAN BE EXPRESSED AS THE NICth
C     PARTICLE IN BOX ICX,ICY,ICZ. EACH PARTICLE CAN ONLY BE IN
C     ONE CELL AT A TIME, SO THIS WRITE SHOULD BE THREAD-SAFE
C
            ICELL(NIC,ICX,ICY,ICZ)=I
C
  11  CONTINUE

C      IF (NIC.GT.NA) STOP

C
C**** START GOING THROUGH CELLS
C
      DO 22 ICX = 1,NCELLX
      DO 22 ICY = 1,NCELLY
      DO 22 ICZ = 1,NCELLZ
         IF (NUMCELL(ICX,ICY,ICZ).EQ.0) GOTO 22

C
C**** SET UP CELL JC BY SCANNING AROUND IC
C
       DO 21 JJCX=-1,1
          JCX = ICX+JJCX
             IF (JCX.LT.1)      JCX=JCX+NCELLX
             IF (JCX.GT.NCELLX) JCX=1
       DO 21 JJCY=-1,1
          JCY = ICY+JJCY
             IF (JCY.LT.1)      JCY=JCY+NCELLY
             IF (JCY.GT.NCELLY) JCY=1
       DO 21 JJCZ=-1,1
          JCZ = ICZ+JJCZ
             IF (JCZ.LT.1)      JCZ=JCZ+NCELLZ
             IF (JCZ.GT.NCELLZ) JCZ=1
C
      IF (NUMCELL(JCX,JCY,JCZ).EQ.0) GOTO 21
C
      I1 = ICX+NCELLX*(ICY-1)+NCELLX*NCELLY*(ICZ-1)
      I2 = JCX+NCELLX*(JCY-1)+NCELLX*NCELLY*(JCZ-1)
      IF (I1.GT.I2) GOTO 21
C
C**** NOW FIND PARTICLES IN SELECTED CELLS
C
      DO 3  NIC = 1,NUMCELL(ICX,ICY,ICZ)
            I = ICELL(NIC,ICX,ICY,ICZ)
      DO 61 NJC = 1,NUMCELL(JCX,JCY,JCZ)
            J = ICELL(NJC,JCX,JCY,JCZ)

C
C**** MAKE SURE SAME CELL PARTICLES ARENT COUNTED TWICE
C
      IF (I1.EQ.I2) THEN
      IF (I.GE.J) GOTO 61
      ENDIF
C
C
C***** CALCULATE DISTANCE BETWEEN PARTICLE I AND J
C
         RX = X(I) - X(J)
         RX  = RX - CUBEX*ANINT(RX/CUBEX)
         RX2 = RX*RX
         IF (RX2.GT.ROUT2) GOTO 61

         RY = Y(I) - Y(J)
         RY  = RY - CUBEY*ANINT(RY/CUBEY)

         RY2 = RY*RY
         IF (RY2.GT.ROUT2) GOTO 61

         RZ = Z(I) - Z(J)
         RZ  = RZ - CUBEZ*ANINT(RZ/CUBEZ)
         RZ2 = RZ*RZ
         IF (RZ2.GT.ROUT2) GOTO 61

         RSQ = RX2 + RY2 + RZ2

C
C***** FIJ = FORCECELL ON I DUE TO J
C***** CALCULATE FORCECELLS.24 AND POTENTIAL ENERGY WITH LJ POTENTIAL
C
C******CALCULATE INTERMOLECULAR FORECES

        IF (RSQ.GT.RMAX) GO TO 61

C***** Add the particle separation to the histogram before we do anything else
C      This saves a lot of time later on
        IBIN = INT(SQRT(RSQ)*RBIN_INV)
        IF (IBIN < NBINS) THEN
          RIJ_HIST(IBIN) = RIJ_HIST(IBIN) + 1
        ENDIF

C***** Force calculation starts here
        IF (IMOL(I).EQ.IMOL(J)) THEN 
          EPS12 = 1  
        ELSE 
          EPS12 = MIX*SQRT(EPS1*EPS2)
        ENDIF

          RSI    = 1.D0/RSQ
          R4     = RSI*RSI
          R6     = R4*RSI
          R12    = R6*R6
          UPOT   = UPOT + ((R12 - 2.0D0*R6) + EPS_LJG*EXP((-1.0D0*(SQRT(1/RSI) - R_LJG)**2.0D0)/(2.0D0*0.02D0)))*EPS12 - EPS12*SHIFT
          RSCALAR= EPS12*(6.0D0*RSI*(2.0D0*R12 - R6) + EPS_LJG*(SQRT(1/RSI) - R_LJG)*EXP((-1.0D0*(SQRT(1/RSI) - R_LJG)**2.0D0)/(2.0D0*0.02D0))/(SQRT(1/RSI)*0.02D0))
          FIJX   = RSCALAR*RX
          FIJY   = RSCALAR*RY
          FIJZ   = RSCALAR*RZ

         FX(I)  = FX(I) + FIJX
         FY(I)  = FY(I) + FIJY
         FZ(I)  = FZ(I) + FIJZ
         FX(J)  = FX(J) - FIJX
         FY(J)  = FY(J) - FIJY
         FZ(J)  = FZ(J) - FIJZ
C******* Update the elements of the pressure tensor
         PT(1,1) = PT(1,1) + RX*FIJX
         PT(1,2) = PT(1,2) + RX*FIJY
         PT(1,3) = PT(1,3) + RX*FIJZ
         PT(2,1) = PT(2,1) + RY*FIJX
         PT(2,2) = PT(2,2) + RY*FIJY
         PT(2,3) = PT(2,3) + RY*FIJZ
         PT(3,1) = PT(3,1) + RZ*FIJX
         PT(3,2) = PT(3,2) + RZ*FIJY
         PT(3,3) = PT(3,3) + RZ*FIJZ
 61   CONTINUE
  3   CONTINUE
 21   CONTINUE
 22   CONTINUE
C
      ANUM = 0.0D0
      ADEN = 0.0D0
      ALPH = 0.0D0
      ANUMB= 0.0D0
      ADENB= 0.0D0
      ALPHB= 0.0D0
      PXPY = 0.0D0
      PXPZ = 0.0D0
      PYPZ = 0.0D0 
C
C     CALCULATE GAUSSIAN ERGOSTAT ALPH
C
C*****ERGOSTAT FOR CONSTANT ENERGY IS NOT IMPLEMENTED FOR A CONFINED SYSTEM
      IF (NGAUS.EQ.1.OR.NGAUS.EQ.3) THEN
         DO 6 I=1,NPART
            ADEN = ADEN + (PX(I)**2 + PY(I)**2+ PZ(I)**2)
            ANUM = ANUM + FIELD*IKOL(I)*PX(I)
 6       CONTINUE
         IF (ADEN.NE.0.0D0) THEN
         ALPH = ANUM/ADEN      
         ENDIF
C
C    ADD FEEDBACK TO ACCOUNT FOR NUMERICAL DRIFT
C
         IF (NGAUS.EQ.3) THEN
            ALPH = ALPH + (E00-E0*NPART)/(E0*NPART-UWALL)
     *                  /10.0D0/DELTA
         END IF
C
C     CALCULATE GAUSSIAN THERMOSTAT ALPH
C     TOP AND BOTTOM WALLS ARE THRMOSTATTED SEPARATELY
C
      ELSE IF (NGAUS.EQ.2.OR.NGAUS.EQ.4) THEN
         DO 7 I=1,NPART
          ADEN = ADEN + (PX(I)**2 + PY(I)**2+
     &                PZ(I)**2)
          ANUM = ANUM + (FX(I)*PX(I)+FY(I)*PY(I)
     &                +FZ(I)*PZ(I)+FIELD*IKOL(I)*PX(I))
C
 7       CONTINUE 
         IF (ADEN.NE.0.0D0) THEN
         ALPH = ANUM/ADEN
         ENDIF
C
C    ADD FEEDBACK TO ACCOUNT FOR NUMERICAL DRIFT
C
         IF (NGAUS.EQ.4) THEN
            ALPH = ALPH+(ADEN-(DF)*TR)/
     *                 (DF)/10.0/DELTA 
         END IF

      ELSE IF (NGAUS.EQ.6) THEN
        ALPH=NHALPH
      ELSE
        ALPH = 0.0D0
      ENDIF
C
C Clean up arrays
C
      DEALLOCATE(NUMCELL)
      DEALLOCATE(ICELL)
      RETURN
      END
C
C----------------------------------------------------------------
C
C***** RANDOM NUMBER GENERATOR
C***** GENERATE RANDOM NUMBERS IN THE RANGE 0,1
C
      SUBROUTINE RAN1(SEED,R)
      INCLUDE "SSGK.inc"
      REAL (KIND=prec) R
      DOUBLE PRECISION D2P31M,D2P31,DSEED
      INTEGER SEED
      DATA D2P31M/2147483647.D0/
      DATA D2P31/2147483648.D0/
      DSEED = SEED
      DSEED = DMOD(16807.D0*DSEED,D2P31M)
      R     = DSEED/D2P31
      SEED = INT(DSEED)
      RETURN
      END
C
C   *****************************************************
C
      SUBROUTINE READ2
C
C***** READ IN COORDINATES AND MOMENTA
C
      INCLUDE "SSGK.inc"
C
C***** LOCAL VARIABLES
C
      INTEGER I,N,NN,NNN
C
      OPEN(UNIT=2,FORM='UNFORMATTED',STATUS='UNKNOWN'
     &,ACCESS='SEQUENTIAL',FILE='F02')
       REWIND 2
       READ(2)N,NN,NNN
       IF (N.NE.NPART) THEN
         WRITE(6,*)'The no of particles in F02 is not equal to NPART'
         REWIND(UNIT=2)
         CLOSE(UNIT=2)
         STOP
      END IF
      IF (NN.NE.NP) THEN
         WRITE(6,*)'The dimension NP is not equal to that used in F02',NN
         REWIND(UNIT=2)
         CLOSE(UNIT=2)
         STOP
      END IF
      IF (NNN.NE.NPART) THEN
         NPART=NNN
      END IF

      READ(2) X,Y,Z,PX,PY,PZ,XEQ,YEQ,ZEQ
      READ(2) S1,S2,S3,S4
      CLOSE(UNIT=2)
      
C
C***** COORDINATES WERE SAVED FOR UNIT BOXLENGTH
C***** THEREFORE SCALE TO SIZE OF BOX
C
      DO 1 I = 1, NPART
         X(I) = X(I)*CUBEX
         Y(I) = Y(I)*CUBEY
         Z(I) = Z(I)*CUBEZ
    1 CONTINUE
C
      RETURN
      END
C
C   ******************************************************
C
      SUBROUTINE READ3
C
C***** OUTPUT ACCUMULATED AVERAGES AND MEANS
C
      INCLUDE "SSGK.inc"
      INTEGER N,NN,NNN
C
      OPEN(UNIT=3,FORM='UNFORMATTED',STATUS='UNKNOWN'
     &,ACCESS='SEQUENTIAL',FILE='F03')
      REWIND 3
      READ(3) N,NN,NNN
      IF (N.NE.NK) THEN
         WRITE(6,*)'Dimension of NK is not equal to that used in F03',N
         REWIND 3
         CLOSE (UNIT=3)
         STOP
      END IF
      IF (NN.NE.NPROP) THEN
         WRITE(6,*)'Dimension NPROP is not equal to that in F03',NN
         REWIND 3
         CLOSE(UNIT=3)
         STOP
      END IF
      IF (NNN.NE.NTCF) THEN
         WRITE(6,*)'Dimension NTCF is not equal to that in F03',NNN
         REWIND 3
         CLOSE(UNIT=3)
         STOP
      END IF
      READ(3) KB,NMEAN,MEAN,KOUNT,TCF,AVER
      CLOSE(UNIT=3)
C
      RETURN
      END
C
C   ****************************************************
C
      SUBROUTINE OUT3
      INCLUDE "SSGK.inc"
C
      OPEN(UNIT=3,FORM='UNFORMATTED',STATUS='UNKNOWN'
     &,ACCESS='SEQUENTIAL',FILE='F03')
      REWIND 3
      WRITE(3) NK,NPROP,NTCF
      WRITE(3) KB,NMEAN,MEAN,KOUNT,TCF,AVER
      CLOSE(UNIT=3)
C
      RETURN
      END
C
C   ***********************************************************
C
      SUBROUTINE OUT4
      INCLUDE "SSGK.inc"
C
      OPEN(UNIT=4,FORM='UNFORMATTED',STATUS='UNKNOWN'
     &,ACCESS='SEQUENTIAL',FILE='F04')
      REWIND 4
      WRITE(4) M0,M1,M2,M3
      CLOSE(UNIT=4)
C
      RETURN
      END
C
C   ***********************************************************
C
      SUBROUTINE READ4
      INCLUDE "SSGK.inc"
C
      OPEN(UNIT=4,FORM='UNFORMATTED',STATUS='UNKNOWN'
     &,ACCESS='SEQUENTIAL',FILE='F04')
      REWIND 4
      READ(4) M0,M1,M2,M3
      CLOSE(UNIT=4)
C
      RETURN
      END
C   ***********************************************************
C
      SUBROUTINE OUT2
C
C***** OUTPUT COORDINATES AND MOMENTA WITH COORDINATES
C***** SCALED TO UNIT BOXLENGTH
C
      INCLUDE "SSGK.inc"
C
C***** LOCAL VARIABLES
C
      INTEGER I
C
C***** RESCALE TO CUBE=1
C
      DO 1 I = 1, NPART
       X(I) = X(I)/CUBEX
       Y(I) = Y(I)/CUBEY
       Z(I) = Z(I)/CUBEZ
    1 CONTINUE
C
      OPEN(UNIT=2,FORM='UNFORMATTED',STATUS='UNKNOWN'
     &,ACCESS='SEQUENTIAL',FILE='F02')
      REWIND 2
      WRITE(2) NPART,NP,NPART
      WRITE(2) X,Y,Z,PX,PY,PZ,XEQ,YEQ,ZEQ
      WRITE(2) S1,S2,S3,S4
      CLOSE(UNIT=2)
C
      DO 2 I = 1, NPART
         X(I) = X(I)*CUBEX
         Y(I) = Y(I)*CUBEY
         Z(I) = Z(I)*CUBEZ
    2 CONTINUE


C
      RETURN
      END
C
C  *********************************************************
C
C----------------------------------------------------------------------
C
C Calculate variance, skewness and kurtosis for the distribution
C of tau-averaged properties
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MOMENTS 
      INCLUDE "SSGK.inc"
C
C***** LOCAL VARIABLES
C
      REAL (KIND=prec) SUMX(NTCF),SUMXX(NTCF),SUMXXX(NTCF),SUMXXXX(NTCF)
      REAL (KIND=prec) AVER1(NTCF),VAR(NTCF),SKEW(NTCF)
      INTEGER I
C
C***** DISPLAY TITLE
C
      DO 30 I=1,NTCF
         VAR(I)=0.0D0 
         SKEW(I)=0.0D0 
         AVER1(I)=0.0D0 
 30   CONTINUE 
      DO 40 I=1,NTCF
         IF (M0(I).GT.0) THEN
            AVER1(I)=M1(I)/M0(I)
            VAR(I) =M2(I)/M0(I)-AVER1(I)**2
            IF (ABS(VAR(I)).GT.1.0D-15) THEN
               SKEW(I)=(M3(I)/M0(I)-3.0D0*AVER1(I)*M2(I)/M0(I)
     *            +2.0D0*AVER1(I)**3)
     *            /VAR(I)**1.5D0
            ELSE 
               SKEW(I)=0.0D0
            END IF
         END IF
 40   CONTINUE
      WRITE(6,*)'No. samples',INT(M0(1))
      IF (M0(1).GT.0) THEN
      WRITE(6,*)'JX/NPART'
      WRITE(6,*)'mean:',AVER1(1)
      WRITE(6,*)'variance:',VAR(1)
      IF (ABS(VAR(1)).GT.1.0D-15) THEN
         WRITE(6,*)'skewness:',SKEW(1)
      ENDIF
      WRITE(6,*)'ALPHA'
      WRITE(6,*)'mean:',AVER1(2)
      WRITE(6,*)'variance:',VAR(2)
      IF (ABS(VAR(2)).GT.1.0D-15) THEN
         WRITE(6,*)'skewness:',SKEW(2)
      ENDIF
      WRITE(6,*)'PRESSURE'
      WRITE(6,*)'mean:',AVER1(3)
      WRITE(6,*)'variance:',VAR(3)
      IF (ABS(VAR(3)).GT.1.0D-15) THEN
         WRITE(6,*)'skewness:',SKEW(3)
      ENDIF
      WRITE(6,*)
      WRITE(6,*)'CHECK DISSIPATION'
      WRITE(6,*)'from alpha:',(DF)*AVER1(2)
      WRITE(6,*)'from Jx:',DBLE(NPART)/TR*FE0*AVER1(1)
     
      END IF
      END
C
      SUBROUTINE VIN
      INCLUDE "SSGK.inc"
      INTEGER I

      DO 10 I=1,NP
         X(I)=XSAV(I)
         Y(I)=YSAV(I)
         Z(I)=ZSAV(I)
         XEQ(I)=XEQSAV(I)
         YEQ(I)=YEQSAV(I)
         ZEQ(I)=ZEQSAV(I)
         PX(I)=PXSAV(I)
         PY(I)=PYSAV(I)
         PZ(I)=PZSAV(I)
         FX(I)=FXSAV(I)
         FY(I)=FYSAV(I)
         FZ(I)=FZSAV(I)
         FXF(I)=FXFSAV(I)
         FYF(I)=FYFSAV(I)
         FZF(I)=FZFSAV(I)
         FXW(I)=FXWSAV(I)
         FYW(I)=FYWSAV(I)
         FZW(I)=FZWSAV(I)     
 10   CONTINUE
      UPOT=UPOTSAV
      UWALL=UWALLSAV
      ALPH=ALSAV
      ALPHB=ALBSAV
      NHALPH=NHSAV

      RETURN
      END
      
      SUBROUTINE VOUT
      INCLUDE "SSGK.inc"
      INTEGER I


      DO 10 I=1,NP
         XSAV(I)=X(I)
         YSAV(I)=Y(I)
         ZSAV(I)=Z(I) 
         XEQSAV(I)=XEQ(I)
         YEQSAV(I)=YEQ(I)
         ZEQSAV(I)=ZEQ(I)
         PXSAV(I)=PX(I)
         PYSAV(I)=PY(I)
         PZSAV(I)=PZ(I) 
         FXSAV(I)=FX(I)
         FYSAV(I)=FY(I)
         FZSAV(I)=FZ(I)
         FXFSAV(I)=FXF(I)
         FYFSAV(I)=FYF(I)
         FZFSAV(I)=FZF(I)
         FXWSAV(I)=FXW(I)
         FYWSAV(I)=FYW(I)
         FZWSAV(I)=FZW(I) 
 10   CONTINUE
      UPOTSAV=UPOT
      UWALLSAV=UWALL
      ALSAV=ALPH
      ALBSAV=ALPHB
      NHSAV=NHALPH

      RETURN
      END


C
      SUBROUTINE FORCENOCELL
C
C***** CALCULATES FORCE ON EACH PARTICLE, TOTAL POTENTIAL ENERGY
C      THERMOSTAT:ALPH, ADN CONFIGURATIONAL PART OF PRESSURE TENSOR
C***** WITH A WCA POTENTIAL AND SIGMA 1
C
      INCLUDE "SSGK.inc"
C
C***** LOCAL VARIABLES
C
      INTEGER I,J
      REAL (KIND=prec) RX,RY,RZ,CIS,RSQ,RSI,R4,R6,R12,RSCALAR
      REAL (KIND=prec) FIJX,FIJY,FIJZ,ANUM,ADEN,PXPY,PXPZ,PYPZ
C
C***** ZERO ARRAYS
C
      FX = 0.0D0
      FY = 0.0D0
      FZ = 0.0D0
C
      PT       = 0.0D0
      UPOT     = 0.0D0
C
C***** CALCULATE DISTANCE BETWEEN PARTICLE I AND J
C
      DO 2 I = 1, NPART-1
      DO 2 J = I+1, NPART
         RX = X(I) - X(J)
         RY = Y(I) - Y(J)
         RZ = Z(I) - Z(J)
C
         RX  = RX - CUBEX*ANINT(RX/CUBEX)
         RY  = RY - CUBEX*ANINT(RY/CUBEX)
         RZ  = RZ - CUBEX*ANINT(RZ/CUBEX)
C
         RSQ =  RX*RX + RY*RY + RZ*RZ
         IF (RSQ.GT.RMAX) GO TO 2
C
C***** FIJ = FORCECELL ON I DUE TO J
C***** CALCULATE FORCECELLS.24 AND POTENTIAL ENERGY WITH WCA POTENTIAL
C
C
         RSI    = 1.0D0/RSQ
         R4     = RSI*RSI
         R6     = R4*RSI
         R12    = R6*R6
         UPOT      = UPOT + 4.0D0*(R12-R6) - SHIFT
         RSCALAR= RSI*(2.D0*R12-R6)
         FIJX   = RSCALAR*RX
         FIJY   = RSCALAR*RY
         FIJZ   = RSCALAR*RZ
         FX(I)  = FX(I) + FIJX
         FY(I)  = FY(I) + FIJY
         FZ(I)  = FZ(I) + FIJZ
         FX(J)  = FX(J) - FIJX
         FY(J)  = FY(J) - FIJY
         FZ(J)  = FZ(J) - FIJZ
C
C***** COMPUTE CONFIGURATIONAL PART OF PRESSURE TENSOR/12
C
         PT(1,1) = PT(1,1) + RX*FIJX
         PT(1,2) = PT(1,2) + RX*FIJY
         PT(1,3) = PT(1,3) + RX*FIJZ
         PT(2,1) = PT(2,1) + RY*FIJX
         PT(2,2) = PT(2,2) + RY*FIJY
         PT(2,3) = PT(2,3) + RY*FIJZ
         PT(3,1) = PT(3,1) + RZ*FIJX
         PT(3,2) = PT(3,2) + RZ*FIJY
         PT(3,3) = PT(3,3) + RZ*FIJZ
C
    2 CONTINUE
C
      DO 5 I=1,NPART
         FX(I) = 24.D0*FX(I)
         FY(I) = 24.D0*FY(I)
         FZ(I) = 24.D0*FZ(I)
c         write(6,*)i,fx(i),fy(i),fz(i)
    5 CONTINUE
C
C     CALCULATE TEMP/ENERGY CONSTRAINT
C
      ANUM = 0.0D0
      ADEN = 0.0D0
      ALPH = 0.0D0
      PXPY = 0.0D0
      PXPZ = 0.0D0
      PYPZ = 0.0D0
C
C     CALCULATE GAUSSIAN ERGOSTAT ALPH
C
      IF (NGAUS.EQ.1.OR.NGAUS.EQ.3) THEN
         DO 6 I=1,NPART
            ADEN = ADEN + PX(I)**2 + PY(I)**2 + PZ(I)**2
            ANUM = ANUM + FIELD*DBLE(IKOL(I))*PX(I)
 6       CONTINUE
         ALPH = ANUM/ADEN      
C
C    ADD FEEDBACK TO ACCOUNT FOR NUMERICAL DRIFT
C
         IF (NGAUS.EQ.3) THEN
            ALPH = ALPH + (E00-E0*NPART)/(E0*NPART-UPOT)
     *                  /10.0D0/DELTA
         END IF
C
C     CALCULATE GAUSSIAN THERMOSTAT ALPH
C
C
      ELSE IF (NGAUS.EQ.2.OR.NGAUS.EQ.4) THEN
         DO 7 I=1,NPART
            ADEN = ADEN + PX(I)**2 + PY(I)**2 + PZ(I)**2
            ANUM = ANUM + FX(I)*PX(I)+FY(I)*PY(I)+FZ(I)*PZ(I)+
     *             DBLE(IKOL(I))*FIELD*PX(I)
 7       CONTINUE
         ALPH = ANUM/ADEN      
C
C    ADD FEEDBACK TO ACCOUNT FOR NUMERICAL DRIFT
C
         IF (NGAUS.EQ.4) THEN
            ALPH = ALPH+(KE2-(3.0D0*NPART-4.0D0)*TR)/
     *                 (3.0D0*NPART-4.0D0)/10.0D0/DELTA
         END IF
      ELSE 
         ALPH = 0.0D0
      ENDIF
C
C**** SCALE PRESSURE TENSOR
C
      PT(1,1) = 24.D0*PT(1,1)
      PT(1,2) = 24.D0*PT(1,2)
      PT(1,3) = 24.D0*PT(1,3)
      PT(2,1) = 24.D0*PT(2,1)
      PT(2,2) = 24.D0*PT(2,2)
      PT(2,3) = 24.D0*PT(2,3)
      PT(3,1) = 24.D0*PT(3,1)
      PT(3,2) = 24.D0*PT(3,2)
      PT(3,3) = 24.D0*PT(3,3)
C
      RETURN
      END
C
C

C      END MODULE