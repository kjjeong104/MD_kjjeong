C       LIQUID STRUCTURE CALCULATION
C       WRITTEN BY CHANDRA N. PATRA
C       MODIFIED BY KYEONG-JUN JEONG
C       2CM DIFFERENT SIZE IN PLANAR GEOMETRY
C       WRITTEN BY CHANDRA N. PATRA
C       ALL RIGHTS RESERVED
C       v03: revised image charge int. from v02
C       self-part:MEP, pairwise:corr fxn
C       ljdft v03 : from noljdft v03, recalled lj adsorption formalism
C       currently has only LJ att potential by electrode surface(PSILA)
C       note:EPP,EPM,EPS:LJ epsilon for cation,anion,surface
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER(PI = 3.141592653589793D0, FPI = 4.0D0*PI)
        PARAMETER(ALH = 30.0D0, DZ = 0.02D0, NMAX=1+(ALH/DZ))
        DIMENSION R(NMAX),DENP(NMAX),DENM(NMAX),DENS(NMAX)
        DIMENSION DENBARP(NMAX),DENBARM(NMAX),DENBARS(NMAX)
        DIMENSION DELDENP(NMAX),DELDENM(NMAX),DELDENS(NMAX)
        DIMENSION FPP(NMAX),FMM(NMAX),FSS(NMAX)
        DIMENSION PSI(NMAX),PSIS(NMAX),PSILAP(NMAX),PSILAM(NMAX)
        DIMENSION C1EL11(NMAX),C1EL12(NMAX),C1EL21(NMAX),C1EL22(NMAX)
        DIMENSION C1IM11(NMAX),C1IM12(NMAX),C1IM21(NMAX),C1IM22(NMAX)
        DIMENSION C1EL1P(NMAX),C1EL1M(NMAX),PIMP(NMAX),PIMM(NMAX)
C
        COMMON/DIST/R
        COMMON/DRATIO/DBD2,DBD3
        COMMON/DENRATIO/XP,XM
        COMMON/DENSITY/DENAVT
        COMMON/MBREAK/MESHBREAK
        COMMON/CHARGE/QP,QM,SIGMA,ESTAR
        COMMON/LJ/EPP,EPM,EPS,DENWALL
C
        CHARACTER(len=64) :: BUFFER
        REAL :: tstart,telapsed
        INTEGER :: RECFREQ
        DO I=1,iargc()
         CALL GETARG(I,BUFFER)
        ENDDO

        OPEN(4,FILE=trim(BUFFER) // '.inp',STATUS='UNKNOWN')
        READ(4,*)CONCAVP,CONCAVM
        READ(4,*)QP,QM,SIGMA
        READ(4,*)EPP,EPM,EPS,DENWALL
        READ(4,*)DVAL
        READ(4,*)DBD2,DBD3
        READ(4,*)TEMP,DIEL
        READ(4,*)CMIX,DEVTOL
        READ(4,*)NOPT,RECFREQ
        READ(4,*)OLDDEN
C
        CLOSE(4,STATUS='KEEP')
C
        BOLTZ = 1.38066D-16
        AVOGAD = 6.023D23
        ELECHARGE = 4.8D-10
        ELEC = 1.602D-19
C
        SUMT_PREV = 1.0D23
C
        BETA = 1.0D0/BOLTZ/TEMP
        ESTAR=ELECHARGE*DSQRT(BETA/(DIEL*DVAL))
        DVALM = DVAL*1.D-2
        SIGMAOLD = SIGMA
        SIGMA = SIGMA*DVALM*DVALM/ELEC
        FACPSI = 1.D3/BETA/1.602D-12
        write(*,*) 'sigma,facpsi,estar=',sigma,facpsi,estar
        CALL cpu_time(tstart)
C
        DENAVP = CONCAVP*AVOGAD*DVAL*DVAL*DVAL/1000.0D0
        DENAVM = CONCAVM*AVOGAD*DVAL*DVAL*DVAL/1000.0D0
C
        CDENTOTAL = QP*DENAVP + QM*DENAVM
        DENAVT = DENAVP + DENAVM
        XP=DENAVP/DENAVT
        XM=DENAVM/DENAVT
C
        MESH = NMAX
        MESHBREAK = MESH - (5.0D0/DZ)
C     
        OPEN(9,FILE=trim(BUFFER) // '-ds.out',STATUS='UNKNOWN')
C
        WRITE(9,111)ALH, MESH
        WRITE(9,*)
        WRITE(9,112)QP,QM,SIGMA
        WRITE(9,*)
        WRITE(9,1121)EPP,EPM,EPS,DENWALL
        WRITE(9,*)
        WRITE(9,1199)DVAL
        WRITE(9,*)
        WRITE(9,1191)DBD2,DBD3
        WRITE(9,*)
        WRITE(9,113)DENAVP,DENAVM,CDENTOTAL,DENAVT,XP,XM
        WRITE(9,*)
        WRITE(9,114)TEMP,DIEL
        WRITE(9,*)
        WRITE(9,1141)CMIX,DEVTOL
        WRITE(9,*)
        WRITE(9,115)NOPT
        WRITE(9,*)
        WRITE(9,116)OLDDEN
C
111     FORMAT(2X,'WALL SEPARATION',4X,F18.8 /
     C         2X,'MESH POINTS',2X,I10)
1199    FORMAT(2X,'HARD SPHERE DIAMETER',4X,D18.10)
1191    FORMAT(2X,'DIAMETER RATIO(BIGION)',4X,E18.8/
     C         2X,'DIAMETER RATIO(SURFACE)',4X,E18.8)
112     FORMAT(2X,'POSITIVE CHARGE ',4X,D18.10/
     C         2X,'NEGATIVE CHARGE ',4X,D18.10/
     C         2X,'SURF CHARGE DENSITY',2X,D18.10)
1121    FORMAT(2X,'CATION LJ WELL DEPTH',4X,D18.10/
     C         2X,'ANION LJ WELL DEPTH',4X,D18.10/
     C         2X,'SURFACE LJ WELL DEPTH',4X,D18.10/
     C         2X,'SURFACE ATOM DENSITY(red)',4X,D18.10)
113     FORMAT(2X,'AVERAGE DENSITY',4X,(7E18.8,2X))
114     FORMAT(2X,'TEMPERATURE        ',4X,D18.10/
     C         2X,'DIELECTRIC CONSTANT',2X,F8.4)
1141     FORMAT(2X,'MIXING COEFFICIENT',2X,E18.8/
     C         2X,'TOLERANCE LIMIT',5X,E18.8)
115     FORMAT(2X,'WRITING OPTION',2X,I10)
116     FORMAT(2X,'OLD DENSITY USED',2X,F18.8)
C
C       INITIALIZE THE MESH POINTS
        DO I=1,MESH
         R(I)=DFLOAT(I-1)*DZ
        ENDDO
C
        AKAPP=DSQRT(FPI*ESTAR*ESTAR*(QP*QP*DENAVP+QM*QM*DENAVM))
        ATERM=0.5D0*FPI*ESTAR*ESTAR*SIGMA*QP/AKAPP
        BTERM=2.0D0*DLOG(ATERM+DSQRT(ATERM*ATERM+1.0D0))
        UTERM=(1.0D0-DEXP(-2.0D0*BTERM/4.0D0))
     C       /(1.0D0+DEXP(-2.0D0*BTERM/4.0D0))
C
        DO I=1,MESH
         UZ=UTERM*DEXP(-AKAPP*(R(I)-R(1)))
         BQPSIZ=2.0*DLOG((1.0+UZ)/(1.0-UZ))
         DENP(I)=DENAVP*DEXP(-BQPSIZ)
         DENM(I)=DENAVM*DEXP(BQPSIZ)
        ENDDO
C
        DO I=1,MESH
         IF(R(I).LT.(0.5D0)) DENP(I) = 0.0D0
         IF(R(I).LT.(0.5D0*DBD2)) DENM(I) = 0.0D0
         DENBARP(I) = DENAVP+DENAVM
         DENBARM(I) = DENAVP+DENAVM
        ENDDO 
C
        IF(OLDDEN .EQ. 1.0D0) THEN
         OPEN(1,FILE=trim(BUFFER) // '-ds.ic',STATUS='UNKNOWN')
         DO I = 1, MESH
          READ(1,117)R(I),DENP(I),DENM(I),DENBARP(I),DENBARM(I)
         ENDDO
         CLOSE(1,STATUS='KEEP')
        ENDIF

C       BEFORE STARTING ITERATION, WRITE INITIAL DENS
        NITER=0
        OPEN(8, FILE =trim(BUFFER) // '-ds.fc',STATUS='UNKNOWN')
        WRITE(8,199)NITER
        DO I = 1, MESH
         WRITE(8,117)R(I),DENP(I),DENM(I),DENS(I),
     C               DENBARP(I),DENBARM(I),DENBARS(I)
        ENDDO
 
C       START ITERATION 
C
105     CONTINUE
        NITER=NITER+1
C        CALL cpu_time(telapsed)
C        WRITE(*,*) 'Step started ',NITER,
C     C ' elapsed time(sec) ',telapsed-tstart
C
        IF(NOPT .EQ. 1)THEN
         IF(MOD(NITER,5) .EQ. 0)THEN
          WRITE(6,*)'GIVE CMIX,DEVTOL'
          READ(5,*)CMIX,DEVTOL
         ENDIF
        ENDIF
C
C       CALCULATE POTENTIAL (\beta*e*\psi(z))
C       Calculate LJ attraction pot(\beta*\psi_ljatt(z))

        CALL COUL(DENP,DENM,PSI,PSIS)
        CALL VDWPOT(EPP,EPS,1.0D0,DBD3,DENWALL,PSILAP)
        CALL VDWPOT(EPM,EPS,DBD2,DBD3,DENWALL,PSILAM)
C
C        CALCULATE C2ELINT
         CALL C2ELINT(DENP,DENM,DENAVP,DENAVM,DENBARP,DENBARM,
     C   C1EL11,C1EL12,C1EL21,C1EL22,C1IM11,C1IM12,C1IM21,C1IM22)
C       DEFINE \rho_0

        DENP1=DENAVP+DENAVM
        DENM2=DENAVM+DENAVP
      
C       CALCULATE C_1(\rho_0)

        CALL CFUNAB(DENP1,DENM2,C11P0,C12P0,C21M0,C22M0)
C
        C1PA=C11P0+C12P0
        C1MA=C21M0+C22M0

C        write(*,*)C11P0,C12P0,C13P0,C21M0,C22M0,C23M0,C31S0,C32S0,C33S0
C
C       CALCULATE DENBAR (\bar{\rho})

        CALL DENBAR(DENP,DENM,DENAVP,DENAVM,DENBARP,DENBARM)
C
         DO 50 I=1,MESH
          DENP1=DENBARP(I)
          DENM2=DENBARM(I)
C
          CALL CFUNAB(DENP1,DENM2,C11P,C12P,C21M,C22M)
C
          C1PB=C11P+C12P
          C1MB=C21M+C22M
C         write(*,*)C11P,C12P,C13P,C21M,C22M,C23M,C31S,C32S,C33S
C        CALCULATE DENSITY (MAIN EQUATION)

C          FP=-QP*PSI(I)+C1PB-C1PA+C1EL11(I)+C1EL12(I)
C          FM=-QM*PSI(I)+C1MB-C1MA+C1EL21(I)+C1EL22(I)
          FP=-QP*PSI(I)-PSILAP(I)+C1PB-C1PA+C1EL11(I)+C1EL12(I)
C     C       +C1LA11(I)+C1LA12(I)
          FM=-QM*PSI(I)-PSILAM(I)+C1MB-C1MA+C1EL21(I)+C1EL22(I)
C     C       +C1LA21(I)+C1LA22(I)
C
          FPP(I) = FP
          FMM(I) = FM
C
50       CONTINUE     
C
        FPPMESH=FPP(MESH)
        FMMMESH=FMM(MESH)
        PSIMESH=PSI(MESH)
C
        DO I=1,MESH
         FP=FPP(I) - FPPMESH
         FM=FMM(I) - FMMMESH
         PSI(I) = PSI(I) - PSIMESH
C
         DELDENP(I) = 0.0D0
         IF(DENP(I) .GT. 0.0D0) THEN
          DELDENP(I)=DENP(I)*(FP-DLOG(DABS(DENP(I)/DENAVP)))
         ENDIF
         DELDENM(I) = 0.0D0
         IF(DENM(I).GT.0.0)THEN
          DELDENM(I)=DENM(I)*(FM-DLOG(DABS(DENM(I)/DENAVM)))
         ENDIF
C
         DENP(I)=(1.0-CMIX)*DENP(I)+CMIX*(DENP(I)+DELDENP(I))
         DENM(I)=(1.0-CMIX)*DENM(I)+CMIX*(DENM(I)+DELDENM(I))
        ENDDO
C
        SUMP=0.0D0
        SUMM=0.0D0
        SUMT=0.0D0
C
        DO I=1,MESH
         IF(DENAVP .GT. 0.0D0) THEN
          SUMP = SUMP + (DELDENP(I)/DENAVP)*(DELDENP(I)/DENAVP)
         ENDIF
         IF(DENAVM .GT. 0.0D0) THEN
          SUMM = SUMM + (DELDENM(I)/DENAVM)*(DELDENM(I)/DENAVM)
         ENDIF
         SUMT = SUMT + SUMP + SUMM 
        ENDDO

C        write(*,*)SUMP,SUMM,SUMS,SUMT
C        
        SUMP=DSQRT(0.5*SUMP/FLOAT(MESH))
        SUMM=DSQRT(0.5*SUMM/FLOAT(MESH))
        SUMT=DSQRT(0.5*SUMT/FLOAT(MESH))
C       Dynamic tuning of mixing coef
        IF((NITER.GT.100).AND.(SUMT_PREV.LT.SUMT))THEN
          CMIX=0.90D0*CMIX
        ENDIF
        SUMT_PREV = SUMT
C
        IF(NOPT.EQ.1)THEN
         WRITE(6,121)NITER,SUMP,SUMM,SUMT
        ENDIF
        IF(MOD(NITER,100).EQ.0)THEN
         CALL cpu_time(telapsed)
         WRITE(9,121)NITER,SUMP,SUMM,SUMT,CMIX,telapsed-tstart
         WRITE(6,121)NITER,SUMP,SUMM,SUMT,CMIX,telapsed-tstart
        ENDIF
121     FORMAT(2X,'NITER',I10,4X,'SUMP',4X,E18.8,4X,
     C        'SUMM',4X,E18.8,4X,'SUMT',4X,E18.8,2X,
     C        'CMIX',2X,E11.4,2X,'time(sec)',4X,F10.4)
C
       IF(MOD(NITER,RECFREQ) .EQ. 0)THEN
C        OPEN(1, FILE = 'pdl2cm-ds.fc',STATUS='UNKNOWN')
        WRITE(8,199)NITER
        DO I = 1, MESH
         WRITE(8,117)R(I),DENP(I),DENM(I),DENBARP(I),DENBARM(I)
        ENDDO
C        CLOSE(1,STATUS = 'KEEP')
       ENDIF
C
       IF((SUMT.GT.DEVTOL).AND.(CMIX.GT.1.0D-6)) GO TO 105
C
C       OPEN(1, FILE = 'pdl2cm-ds.fc',STATUS='UNKNOWN')
       WRITE(8,199)NITER
       DO I = 1, MESH
        WRITE(8,117)R(I),DENP(I),DENM(I),DENBARP(I),DENBARM(I)
       ENDDO
       CLOSE(8,STATUS = 'KEEP')
C
C       WRITE(9,*)'R,RHO(I),GRI,RHOBAR(I),FS(I)'
C
C      for Helmholtz layer analysis, define counterion first.
       ONECI=DBD2
C
       OPEN(1, FILE =trim(BUFFER) // '.out', STATUS = 'UNKNOWN')
       OPEN(2, FILE =trim(BUFFER) // '.pot', STATUS = 'UNKNOWN') 
       CUMULCD=0.0D0
       PREVCCD=0.0D0
       DO I=1,MESH
        C1EL1P(I) = C1EL11(I) + C1EL12(I)
        C1EL1M(I) = C1EL21(I) + C1EL22(I)
        PIMP(I) = -QP*PSIS(I) + C1IM11(I) + C1IM12(I)
        PIMM(I) = -QM*PSIS(I) + C1IM21(I) + C1IM22(I)
        IF(DENP(I).GT.0.0D0) THEN
          GPI = DENP(I)/(DENP(I)+DENM(I))
          GMI = DENM(I)/(DENP(I)+DENM(I))
        ENDIF
        RI = 2.D0*R(I)/(0.5+DBD2)
        DENPI = DENP(I)/DENAVP
        DENMI = DENM(I)/DENAVM
        CDENI = (QP*DENPI*CONCAVP+QM*DENMI*CONCAVM)*AVOGAD*ELEC*1000.0D0
        CUMULCD=CUMULCD+CDENI*DZ*DVALM
        IF(DABS(R(I)-ONECI).LT.DZ) THEN
          CD1CI=CUMULCD
        ENDIF
C        WRITE(1,122)R(I),DENPI,DENMI,PSI,DENP(I),DENM(I),DENBARP(I),
C     C              DENBARM(I),C1EL1P(I),C1EL1M(I)
C       output modified. x axis is R/sigma_1. Don't use variable RI.
        WRITE(1,122)R(I),DENPI,DENMI,PSI(I),PSILAP(I),PSILAM(I),CUMULCD
        WRITE(2,122)R(I),PSI(I),PSIS(I),C1EL1P(I),C1EL1M(I),
     C              PIMP(I),PIMM(I),PSILAP(I),PSILAM(I)
C        WRITE(9,122)R(I),DENPI,DENMI,PSI,DENP(I),DENM(I),DENBARP(I),
C     C              DENBARM(I),C1EL1P(I),C1EL1M(I)
       ENDDO
       WRITE(9,123)ONECI,CD1CI
       WRITE(9,124)PSI(1)
       IF(CMIX.LE.1.0D-6) THEN
         WRITE(9,*)'Converging Failure'
       ENDIF
       CLOSE(1,STATUS = 'KEEP')
       CLOSE(2,STATUS = 'KEEP')
122    FORMAT(1X,14(E18.8,2X))
123    FORMAT(1X,'POS_and_fluidICdist_at_R=1d_ci',2X,F18.8,2X,F18.8)
124    FORMAT(1X,'Elst.potential at z=0',2X,F18.8)
117    FORMAT(14(E18.8,2X))
199    FORMAT(2X,'R,DEN(P,M),DENBAR(P,M) of ITER ',I10)
C
        STOP
        END
C
        SUBROUTINE SIMPBAR(F,VALINT,MESH1,MESH2)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER(ALH = 30.0D0, DZ = 0.02D0, NMAX=1+(ALH/DZ))
        DIMENSION F(NMAX)
C
        A1=0.0
        DO I=MESH1+1,MESH2
         A1=A1+(F(I)+F(I-1))
        ENDDO
        VALINT=A1*DZ/2.0D0
        RETURN
        END
C
C       v02: image charge interaction (self + mutual?)
C       v03: revised scheme. include here only self-part.
C JPC Lett 2016,7,2753
C beta*U_self(d) = - l_B/2d = - beta e^2/(2d epsilon) (for charge e)
C VALC is the storage for image charge integrand
        SUBROUTINE COUL(DENP,DENM,PSI,PSIS)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER(PI = 3.141592653589793D0, FPI = 4.0D0*PI)
        PARAMETER(ALH = 30.0D0, DZ = 0.02D0, NMAX=1+(ALH/DZ))
        DIMENSION R(NMAX),DENP(NMAX),DENM(NMAX)
        DIMENSION VALA(NMAX),VALB(NMAX),VALC(NMAX)
        DIMENSION PSI(NMAX),PSIS(NMAX)
        COMMON/DIST/R
        COMMON/CHARGE/QP,QM,SIGMA,ESTAR
C
        MESH = NMAX
        AP1=0.0D0
        BP1=0.0D0
        AM1=0.0D0
        BM1=0.0D0
        C1=0.0D0

        DO I=2,MESH
         AP1=AP1+QP*DENP(I)+QP*DENP(I-1)
         AM1=AM1+QM*DENM(I)+QM*DENM(I-1)
         BP1=BP1+QP*(DENP(I)*R(I)+DENP(I-1)*R(I-1))
         BM1=BM1+QM*(DENM(I)*R(I)+DENM(I-1)*R(I-1))
         IF(R(I-1).NE.0.0D0)THEN
         C1=(QP*DENP(I)+QM*DENM(I))/(2.0D0*R(I)) +
     C      (QP*DENP(I-1)+QM*DENM(I-1))/(2.0D0*R(I-1))
         ENDIF
         VALA(I)=(AP1+AM1)*DZ/2.0D0
         VALB(I)=(BP1+BM1)*DZ/2.0D0
         VALC(I)=-C1*DZ/2.0D0
        ENDDO
        VALA(1)=0.0D0
        VALB(1)=0.0D0
        VALC(1)=0.0D0
C
        DO I=1,MESH
         PSI(I) = 0.0D0
          PSI(I)=FPI*ESTAR*ESTAR*
     C        ((-VALA(I)-SIGMA)*R(I)-(VALB(MESH)-VALB(I)))
         PSIS(I) =ESTAR*ESTAR*VALC(I)
         PSI(I) = PSI(I) + PSIS(I) 
        ENDDO
C
        RETURN
        END
C
C       VDWPOT version 01. LJ field by electrode surface
C       without consideration of liquid layers
        SUBROUTINE VDWPOT(EP,EPS,DBD,DBD3,DENWALL,PSILA)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER(PI = 3.141592653589793D0, FPI = 4.0D0*PI)
        PARAMETER(ALH = 30.0D0, DZ = 0.02D0, NMAX=1+(ALH/DZ))
        DIMENSION R(NMAX)
        DIMENSION VALA(NMAX),VALB(NMAX)
        DIMENSION PSILA(NMAX)
        COMMON/DIST/R
C        COMMON/DRATIO/DBD3
C        COMMON/LJ/EPS,DENWALL
C
        MESH = NMAX
        SNUM = 1.0D0
        IF(EP.LT.0.0D0)THEN
         SNUM = -1.0D0
        ENDIF
        EPMIX = DSQRT(DABS(EP)*EPS)*SNUM
        SMIX = (DBD+DBD3)/2.0D0
        HSLIM = DBD/2.0D0
        AP1=0.0D0
        BP1=0.0D0
        AM1=0.0D0
        BM1=0.0D0
        SM3=SMIX**3
        SM9=SM3*SM3*SM3
C       This version:single-body potential energy.(wall-pot)
C       terms here are not integrands.
C       Also, regard surface atom radius to calc LJ int.
C       (coord z=0 is "surface", not "center" of surface atom layer)
        DO I=2,MESH
         R3=(R(I)+DBD3/2.0D0)**3
         R9=R3**R3**R3
         REP=SM3*(SM9/R9)/15.0D0
         ATT=SM3*(SM3/R3)/2.0D0
         VALA(I)=REP
         VALB(I)=ATT
        ENDDO
        VALA(1)=0.0D0
        VALB(1)=0.0D0
C       Keep HS constraint to prevent glitch. Set PSILA=0 for inside
        DO I=1,MESH
         PSILA(I) = 0.0D0
         IF(R(I).GE.HSLIM)THEN
          PSILA(I)=FPI*DENWALL*EPMIX*(VALA(I)-VALB(I))/3.0D0
         ENDIF
        ENDDO
C
        RETURN
        END
C
        SUBROUTINE DENBAR(DENP,DENM,DENAVP,DENAVM,DENBARP,DENBARM)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER(PI = 3.141592653589793D0, FPI = 4.0D0*PI)
        PARAMETER(ALH=30.0D0, DZ = 0.02D0, NMAX=1+(ALH/DZ)) 
        DIMENSION R(NMAX),DENP(NMAX),DENM(NMAX)
        DIMENSION DENBARP(NMAX),DENBARM(NMAX)
        DIMENSION WP11(NMAX),WP12(NMAX)
        DIMENSION WM21(NMAX),WM22(NMAX)
        DIMENSION FUNCP(NMAX),FUNCM(NMAX)
C
        COMMON/DIST/R
        COMMON/DRATIO/DBD2
        COMMON/MBREAK/MESHBREAK
C        
        MESH = NMAX
C
        AM12 = (DBD2+1.0D0)/2.0D0
        AN12 = DABS(DBD2-1.0D0)/2.0D0
C
        DBD22 = DBD2*DBD2
        DBD23 = DBD22*DBD2
        DBD24 = DBD23*DBD2
        DBD25 = DBD23*DBD22
C
        DENNDT=DENAVP + DENAVM 
        CALL DENBARAB(DENNDT,PC11T,PC12T,PC21T,PC22T,
     C           PA1T,PA2T,PB1T,PB2T,PB0T,PD0T)
C
C        write(*,*)pa1t,pa2t,pb1t,pb2t,pb12t,pb23t,pd0t,
C     C            pc11t,pc12t,pc22t
C
        DO 100 I=1,MESH
C
         DO 200 J=1,MESH
C
          WP11(J)=0.0
          WP12(J)=0.0
          WM21(J)=0.0
          WM22(J)=0.0
          FUNCP(J)=0.0
          FUNCM(J)=0.0
          C11TRYP = 0.0D0
          C12TRYP = 0.0D0
          C21TRYM = 0.0D0
          C22TRYM = 0.0D0
C
          DSR=DABS(R(I)-R(J))
          DSR2=DSR*DSR
          DSR3=DSR2*DSR
          DSR5=DSR3*DSR2
          DSR12=DSR-AN12
          DSR122=DSR12*DSR12
          DSR123=DSR122*DSR12
          DSR124=DSR123*DSR12
          DSR125=DSR123*DSR122
C
          IF(DSR.LE.1.0D0) THEN
           C11TRYP=-2.0D0*PI*(PA1T*(1.0D0-DSR2)/2.0D0
     C          +PB1T*(1.0D0-DSR3)/3.0D0+PD0T*(1.0D0-DSR5)/5.0D0)
          ENDIF
          IF(DSR.LE.DBD2)THEN
           C22TRYM=-2.0*PI*(PA2T*(DBD22-DSR2)/2.0D0
     C          +PB2T*(DBD23-DSR3)/3.0D0+PD0T*(DBD25-DSR5)/5.0D0)
          ENDIF
          IF(DSR.LE.AN12)THEN
           C12TRYP=-2.0*PI*(PA1T*(AM12*AM12-DSR2)/2.0D0
     C             +PB0T/3.0+PD0T*(AN12+1.0D0/5.0D0))
           C21TRYM=C12TRYP
          ENDIF
          IF(DSR.GT.AN12 .AND. DSR.LE.AM12)THEN
           C12TRYP=-2.0*PI*(PA1T*(AM12*AM12-DSR2)/2.0D0
     C             +PB0T*(1.0D0-DSR123)/3.0D0+AN12*PD0T*(1.0D0-DSR124)
     C             +PD0T*(1.0D0-DSR125)/5.0D0)
           C21TRYM=C12TRYP
          ENDIF
C
          WP11(J)=C11TRYP/PC11T
          WP12(J)=C12TRYP/PC12T
          WM21(J)=C21TRYM/PC21T
          WM22(J)=C22TRYM/PC22T
C
          FUNCP(J)=WP11(J)*DENP(J)+WP12(J)*DENM(J)
          FUNCM(J)=WM21(J)*DENP(J)+WM22(J)*DENM(J)
C
200      CONTINUE
C
         CALL SIMPBAR(FUNCP,VALP,1,MESH)
         CALL SIMPBAR(FUNCM,VALM,1,MESH)
C
         CALL SIMPBAR(WP11,WNORMP11,1,MESH)
         CALL SIMPBAR(WP12,WNORMP12,1,MESH)
         CALL SIMPBAR(WM21,WNORMM21,1,MESH)
         CALL SIMPBAR(WM22,WNORMM22,1,MESH)
C
         DENBARP(I)=VALP
         DENBARM(I)=VALM
         IF(I .GT. MESHBREAK) THEN
          DENBARP(I) = DENAVP + DENAVM
          DENBARM(I) = DENAVP + DENAVM
         ENDIF
C         IF(I .EQ. 400) write(*,*)WNORMP11,WNORMP12,WNORMP13,
C     C   WNORMM21,WNORMM22,WNORMM23,WNORMS31,WNORMS32,WNORMS33

100     CONTINUE
C
        RETURN
        END
C
C       v02 C2ELINT : introduced additional integral for image charge
       SUBROUTINE C2ELINT(DENP,DENM,DENAVP,DENAVM,DENBARP,DENBARM,
     C  C1EL11,C1EL12,C1EL21,C1EL22,C1IM11,C1IM12,C1IM21,C1IM22)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       PARAMETER(TOL=1.D-5)
       PARAMETER(PI = 3.141592653589793D0, FPI = 4.0D0*PI)
       PARAMETER(ALH = 30.0D0, DZ = 0.02D0, NMAX=1+(ALH/DZ))
       DIMENSION R(NMAX),DENP(NMAX),DENM(NMAX)
       DIMENSION DENBARP(NMAX),DENBARM(NMAX)
       DIMENSION C2EL11(NMAX),C2EL12(NMAX),C2EL21(NMAX),C2EL22(NMAX)
       DIMENSION C1El11(NMAX),C1EL12(NMAX),C1EL21(NMAX),C1EL22(NMAX)
       DIMENSION C2IM11(NMAX),C2IM12(NMAX),C2IM21(NMAX),C2IM22(NMAX)
       DIMENSION C1IM11(NMAX),C1IM12(NMAX),C1IM21(NMAX),C1IM22(NMAX)
C
       COMMON/DIST/R
       COMMON/CHARGE/QP,QM,SIGMA,ESTAR
       COMMON/DRATIO/DBD2
       COMMON/DENRATIO/XP,XM
       COMMON/MBREAK/MESHBREAK
       EXTERNAL GAMMA
C
       MESH = NMAX
       X1I = -4.0D0
       X2I = 4.0D0
       XCAL= BRENT(GAMMA,X1I,X2I,TOL)
C       WRITE (6,*) XCAL
       XG = XCAL
C
       AM1 = (DBD2+1.0D0)
       AN1 = DABS(DBD2-1.0D0)
       AM12 = (DBD2+1.0D0)/2.0D0
       AN12 = DABS(DBD2-1.0D0)/2.0D0
       DBD22 = DBD2**2
       DBD23 = DBD2**3
       DBD25 = DBD2**5
       ETA = (PI/6.0)*(DENAVP+DBD23*DENAVM)
       ETA1M = 1.D0-ETA
       CVAL = PI/ETA1M/2.0D0
       CHID = 1.D0 + CVAL*(DENAVP/(1.D0+XG) 
     C        + DENAVM*DBD23/(1.D0+XG*DBD2))
       CHIN = -CVAL*(DENAVP*QP/(1.D0+XG) 
     C        + DENAVM*QM*DBD2/(1.D0+XG*DBD2))
       CHI = CHIN/CHID
       XI1 = (QP + CHI)/(1+XG)
       XI2 = (QM + CHI*DBD22)/(1+XG*DBD2)
       DGAMMA = DENAVP*XI1*XI1 + DENAVM*XI2*XI2
       Y1 = (1.D0+XG)
       Y2 = (1.D0+XG*DBD2)
       YVAL = Y1*Y2
       ADA=ESTAR*ESTAR
C
C      COEFFICIENTS********
       AB0 = 2.D0*ADA*(CHI*CHI*(XG-2.D0)/3.D0/Y1 + CHI*QP*AN1/YVAL
     C      - QP*QM*XG/Y2)
       AL0 = (ADA*AN1*AN1/16.D0/YVAL)*
     C (CHI*CHI*(4.D0*(1.D0+DBD22) - 4.D0*XG*XG*DBD22 - AN1*AN1*YVAL) 
     C       + 4.D0*(QP+QM)*CHI + 4.D0*QP*QM*XG*XG)  
       AL1 = (ADA/YVAL)*((QP*AN1-QM*AN1)*CHI - 2.D0*XG*QP*QM 
     C       - AM1*QP*QM*XG*XG 
     C       +(Y2*(XG-2.D0)/3.D0+DBD23*Y1*(XG*DBD2-2.D0)/3.D0)*CHI*CHI) 
       AL2 = (ADA/YVAL)*((QP+QM)*CHI + XG*XG*QP*QM
     C       + (1.D0+DBD22-XG*XG*DBD22-AN1*AN1*YVAL/2.D0)*CHI*CHI)
       AL3 = ADA*CHI*CHI/3.D0
       P12 = ADA*QP*QM
C
       AI0 = (ADA/Y1/Y1)*(-2.D0*XG*QP*QP - 2.D0*QP*QP*XG*XG
     C       +(Y1*(XG-2.D0)*2.D0/3.D0)*CHI*CHI)  
       AI1 = (ADA/Y1/Y1)*
     C      ((2.D0-XG*XG)*CHI*CHI + 2.D0*QP*CHI + XG*XG*QP*QP)
       AI2 = ADA*CHI*CHI/3.D0
       P11 = ADA*QP*QP
C
       AJ0 = (ADA/Y2/Y2)*(-2.D0*XG*QM*QM - 2.D0*DBD2*QM*QM*XG*XG
     C       +(Y2*(XG*DBD2-2.D0)*DBD23*2.D0/3.D0)*CHI*CHI)  
       AJ1 = (ADA/Y2/Y2)*
     C       (2.D0*QM*CHI + XG*XG*QM*QM 
     C      + (2.D0*DBD22-XG*XG*DBD22*DBD22)*CHI*CHI)
       AJ2 = ADA*CHI*CHI/3.D0
       P22 = ADA*QM*QM
C
       DO 9100 I=1,MESH
C        DENXP = DABS(QM)*DENBARP(I)/(DABS(QP)+DABS(QM))
C        DENXM = DABS(QP)*DENBARP(I)/(DABS(QP)+DABS(QM))
C        XVAL=DSQRT(FPI*ESTAR*ESTAR*
C     C            (QP*QP*DENXP+QM*QM*DENXM))
C        BVAL=(XVAL+1.0D0-DSQRT(1.0D0+2.0D0*XVAL))/XVAL
C        ADA=ESTAR*ESTAR
C        AEA=ADA*BVAL
C        AFA=AEA*BVAL
        C1EL11(I) = 0.0D0
        C1EL12(I) = 0.0D0
        C1EL21(I) = 0.0D0
        C1EL22(I) = 0.0D0
        DO 9200 J=1,MESH
         C2EL11(J)=0.0
         C2EL12(J)=0.0
         C2EL21(J)=0.0
         C2EL22(J)=0.0
         C2IM11(J)=0.0
         C2IM12(J)=0.0
         C2IM21(J)=0.0
         C2IM22(J)=0.0
         TERM11=0.0D0
         TERM12=0.0D0
         TERM22=0.0D0
         TIM11=0.0D0
         TIM12=0.0D0
         TIM22=0.0D0
C
         DSR=DABS(R(I)-R(J))
         DSR2=DSR*DSR
         DSR3=DSR2*DSR
         DSR5=DSR3*DSR2
         DSR12=DSR-AN12
         DSR122=DSR12*DSR12
         DSR123=DSR122*DSR12
         DSR124=DSR123*DSR12
         DSR125=DSR123*DSR122
         DAR=DABS(R(I)+R(J))
         DAR2=DAR*DAR
         DAR3=DAR2*DAR
         DAR5=DAR3*DAR2
         DAR12=DAR-AN12
         DAR122=DAR12*DAR12
         DAR123=DAR122*DAR12
         DAR124=DAR123*DAR12
         DAR125=DAR123*DAR122
C       
         IF(DSR.LE.1.0D0) THEN
          TERM11 = 2.0D0*PI*(P11*(1.0-DSR)
     C          + AI0*(1.D0-DSR2)/2.D0
     C          + AI1*(1.D0-DSR3)/3.0D0
     C          + AI2*(1.D0-DSR5)/5.D0)
         ENDIF
C
         IF(DSR.LE.DBD2)THEN
          TERM22 = 2.0D0*PI*(P22*(DBD2-DSR)
     C          + AJ0*(DBD22-DSR2)/2.D0
     C          + AJ1*(DBD23-DSR3)/3.0D0
     C          + AJ2*(DBD25-DSR5)/5.D0)
         ENDIF
C
         IF(DSR.LE.AN12)THEN
          TERM12 = 2.0D0*PI*(P12*(AM12-DSR)
     C           + AB0*(AN12*AN12-DSR2)/2.D0
     C           + AL0*(AM12-AN12)
     C           + AL1*(AM12*AM12-AN12*AN12)/2.D0
     C           + Al2*(AM12**3-AN12**3)/3.D0
     C           + AL3*(AM12**5-AN12**5)/5.D0)
         ENDIF
C
         IF(DSR.GT.AN12 .AND. DSR.LE.AM12)THEN
          TERM12 = 2.0D0*PI*(P12*(AM12-DSR)
     C           + AL0*(AM12-DSR)
     C           + AL1*(AM12*AM12-DSR2)/2.0D0
     C           + AL2*(AM12**3-DSR3)/3.D0
     C           + AL3*(AM12**5-DSR5)/5.D0)

         ENDIF
CCCCCCCCCCCCCCCCCCCCCC
         IF(DAR.LE.1.0D0) THEN
          TIM11 = 2.0D0*PI*(P11*(1.0-DAR)
     C          + AI0*(1.D0-DAR2)/2.D0
     C          + AI1*(1.D0-DAR3)/3.0D0
     C          + AI2*(1.D0-DAR5)/5.D0)
         ENDIF
C
         IF(DAR.LE.DBD2)THEN
          TIM22 = 2.0D0*PI*(P22*(DBD2-DAR)
     C          + AJ0*(DBD22-DAR2)/2.D0
     C          + AJ1*(DBD23-DAR3)/3.0D0
     C          + AJ2*(DBD25-DAR5)/5.D0)
         ENDIF
C
         IF(DAR.LE.AN12)THEN
          TIM12 = 2.0D0*PI*(P12*(AM12-DAR)
     C           + AB0*(AN12*AN12-DAR2)/2.D0
     C           + AL0*(AM12-AN12)
     C           + AL1*(AM12*AM12-AN12*AN12)/2.D0
     C           + Al2*(AM12**3-AN12**3)/3.D0
     C           + AL3*(AM12**5-AN12**5)/5.D0)
         ENDIF
C
         IF(DAR.GT.AN12 .AND. DAR.LE.AM12)THEN
          TIM12 = 2.0D0*PI*(P12*(AM12-DAR)
     C           + AL0*(AM12-DAR)
     C           + AL1*(AM12*AM12-DAR2)/2.0D0
     C           + AL2*(AM12**3-DAR3)/3.D0
     C           + AL3*(AM12**5-DAR5)/5.D0)

         ENDIF
         C2EL11(J)=TERM11*(DENP(J)-DENAVP)
         C2EL22(J)=TERM22*(DENM(J)-DENAVM)
         C2EL12(J)=TERM12*(DENM(J)-DENAVM)
         C2EL21(J)=TERM12*(DENP(J)-DENAVP)
C
         C2IM11(J)=(-2.0D0)*TIM11*(DENP(J)-DENAVP)
         C2IM22(J)=(-2.0D0)*TIM22*(DENM(J)-DENAVM)
         C2IM12(J)=(-2.0D0)*TIM12*(DENM(J)-DENAVM)
         C2IM21(J)=(-2.0D0)*TIM12*(DENP(J)-DENAVP)
         IF(I.EQ.J)THEN
          C2IM11(J)=0.0D0
          C2IM22(J)=0.0D0
          C2IM12(J)=0.0D0
          C2IM21(J)=0.0D0
         ENDIF
9200    CONTINUE
C
         CALL SIMPBAR(C2EL11,VALP11,1,MESH)
         CALL SIMPBAR(C2EL12,VALM12,1,MESH)
         CALL SIMPBAR(C2EL21,VALP21,1,MESH)
         CALL SIMPBAR(C2EL22,VALM22,1,MESH)
         CALL SIMPBAR(C2IM11,VALIMP11,1,MESH)
         CALL SIMPBAR(C2IM12,VALIMM12,1,MESH)
         CALL SIMPBAR(C2IM21,VALIMP21,1,MESH)
         CALL SIMPBAR(C2IM22,VALIMM22,1,MESH)
         C1EL11(I)=VALP11+VALIMP11
         C1EL12(I)=VALM12+VALIMM12
         C1EL21(I)=VALP21+VALIMP21
         C1EL22(I)=VALM22+VALIMM22
         C1IM11(I)=VALIMP11
         C1IM12(I)=VALIMM12
         C1IM21(I)=VALIMP21
         C1IM22(I)=VALIMM22
         IF(I.GT.MESHBREAK)THEN
          C1EL11(I) = 0.0D0
          C1EL12(I) = 0.0D0
          C1EL21(I) = 0.0D0
          C1EL22(I) = 0.0D0
          C1IM11(I) = 0.0D0
          C1IM12(I) = 0.0D0
          C1IM21(I) = 0.0D0
          C1IM22(I) = 0.0D0
         ENDIF
9100    CONTINUE
C
        RETURN
        END
C
        SUBROUTINE DENBARAB(DEN,PC11,PC12,PC21,PC22,
     C             A1,A2,B1,B2,B0,D0)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
        PARAMETER(PI = 3.141592653589793D0, FPI = 4.0D0*PI)
        COMMON/DRATIO/DBD2
        COMMON/DENRATIO/XP,XM
C
        DBDN2 = 1.0D0/DBD2
        DBD23 = DBD2**3
        DBDN23 = 1.0D0/DBD23
        DBD24 = DBD23*DBD2
        DBD26 = DBD23*DBD23
C
        AM12 = (DBD2+1.0D0)/2.0D0
        AN12 = DABS(DBD2-1.0D0)/2.0D0
C
        AMN12 = AM12*2.0D0/DBD2
        ANN12 = AN12*2.0D0/DBD2
C
        ETA1F=XP/(XP+DBD23*XM)
        ETA2F=DBD23*XM/(XP+DBD23*XM)
C
        ETA=(PI/6.0)*DEN*(XP+DBD23*XM)
        ETA1=ETA*ETA1F
        ETA2=ETA*ETA2F
        ETA1M=1.0-ETA
        ETA12=ETA1M*ETA1M
        ETA13=ETA12*ETA1M
        ETA14=ETA12*ETA12
C
        BETAPC = ((1.0D0+ETA+ETA**2)*(ETA1+DBDN23*ETA2)
     C         - 3.0D0*(ETA1*ETA2*ANN12*ANN12*
     C          (AMN12+ETA1+DBDN2*ETA2)))/ETA13
C
        A1=((1.0+ETA+ETA*ETA)+(ETA1+DBDN23*ETA2)*
     C      (1.0+2.0*ETA)-3.0D0*(ETA2*ANN12*ANN12*
     C      (AMN12+ETA1+DBDN2*ETA2)
     C      +ETA1*ETA2*ANN12*ANN12))/ETA13
     C      +3.0D0*BETAPC/ETA1M
      A2=DBD23*(DBDN23*(1.0+ETA+ETA*ETA)+(ETA1+DBDN23*ETA2)*
     C     (1.0+2.0*ETA)-3.0D0*(ETA1*ANN12*ANN12*
     C      (AMN12+ETA1+DBDN2*ETA2)
     C      +ETA1*ETA2*DBDN2*ANN12*ANN12))/ETA13
     C      +DBD23*3.0D0*BETAPC/ETA1M
        G11=(1.0-ETA+1.5D0*(ETA1+DBDN2*ETA2))/ETA12
        G22=(1.0-ETA+1.5*DBD2*(ETA1+DBDN2*ETA2))/ETA12
        G12=(G11+DBDN2*G22)/AMN12
C
        B1=-6.0*(ETA1*G11*G11+ETA2*DBDN2*AMN12*AMN12*G12*G12/4.0D0)
        B2=-6.0*DBDN2*(ETA2*G22*G22+ETA1*DBD23*AMN12**2*G12*G12/4.0D0)
        B0=-6.0*(ETA1*AM12*G11*G12+ETA2*DBDN2*AMN12*G12*G22/2.0D0)
        D0=(ETA1*A1+DBDN23*ETA2*A2)/2.0D0
        PC11=-FPI*(A1/3.0+B1/4.0+D0/6.0)
        PC22=-FPI*(DBD23*A2/3.0+DBD24*B2/4.0+DBD26*D0/6.0)
C
        AM122=AM12*AM12
        AM123=AM122*AM12
        AM124=AM123*AM12
        AM125=AM123*AM122
        AM126=AM124*AM122
        AN122=AN12*AN12
        AN123=AN122*AN12
        AN124=AN123*AN12
        AN125=AN123*AN122
        AN126=AN124*AN122
C
        X11=AM124/4.0D0-2.0*AM123*AN12/3.0D0+AM122*AN122/2.0D0
     C    -AN124/12.0D0
        X12=4.0*AN12*(AM125/5.0D0-3.0D0*AM124*AN12/4.0D0+AM123*AN122
     C    -AM122*AN123/2.0D0+AN125/20.0D0)
        X13=AM126/6.0D0-4.0D0*AM125*AN12/5.0D0+3.0D0*AM124*AN122/2.0D0
     C    -4.0D0*AM123*AN123/3.0D0+AM122*AN124/2.0D0-AN126/30.0D0
        PC12=-FPI*(A1*AM123/3.0+B0*X11+D0*(X12+X13))
        PC21=PC12
C
        RETURN
        END
C
        SUBROUTINE CFUNAB(DENP1,DENM2,C11P0,C12P0,C21M0,C22M0)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER(PI = 3.141592653589793D0, FPI = 4.0D0*PI)
        COMMON/DENRATIO/XP,XM
C
        DENPL = 0.0D0
        DENPU = DENP1
        DDENP = 0.001D0
        NDP = (DENPU-DENPL)/DDENP + 1
        DDDENP = (DENPU - DENPL)/(NDP)
        SUM11P = 0.0D0
        SUM12P = 0.0D0
C 
        DO I = 1, NDP
         RHOP = DENPL + (I-1)*DDDENP
         CALL DENBARAB(RHOP,PC11P,PC12P,PC21P,PC22P,
     C    A1P,A2P,B1P,B2P,B0P,D0P)
         SUM11P = SUM11P + PC11P
         SUM12P = SUM12P + PC12P
         IF(I.NE. 1 .OR. I .NE. NDP) THEN
          SUM11P = SUM11P + PC11P
          SUM12P = SUM12P + PC12P
         ENDIF
        ENDDO
        C11P0 = XP*SUM11P*0.5D0*DDDENP
        C12P0 = XM*SUM12P*0.5D0*DDDENP
C
        DENML = 0.0D0
        DENMU = DENM2
        DDENM = 0.001D0
        NDM = (DENMU-DENML)/DDENM + 1
        DDDENM = (DENMU - DENML)/(NDM)
        SUM21M = 0.0D0
        SUM22M = 0.0D0
        DO I = 1, NDM
         RHOM = DENML + (I-1)*DDDENM 
         CALL DENBARAB(RHOM,PC11M,PC12M,PC21M,PC22M,
     C                 A1M,A2M,B1M,B2M,B0M,D0M)
         SUM21M = SUM21M + PC21M
         SUM22M = SUM22M + PC22M
         IF(I .NE. 1 .OR. I .NE. NDM) THEN
          SUM21M = SUM21M + PC21M
          SUM22M = SUM22M + PC22M
         ENDIF
        ENDDO
        C21M0 = XP*SUM21M*0.5D0*DDDENM
        C22M0 = XM*SUM22M*0.5D0*DDDENM
C
        RETURN
        END
C
      DOUBLE PRECISION FUNCTION GAMMA(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(PI = 3.141592653589793D0)
       COMMON/CHARGE/QP,QM,SIGMA,ESTAR
       COMMON/DRATIO/DBD2
       COMMON/DENRATIO/XP,XM
       COMMON/DENSITY/DENAVT
C
      DBD22 = DBD2**2
      DBD23 = DBD2**3
       DENAVP = XP*DENAVT
       DENAVM = XM*DENAVT
       ETA = (PI/6.0)*(DENAVP+DBD23*DENAVM)
       ETA1M = 1.D0-ETA
       CVAL = PI/ETA1M/2.0D0
       CHID = 1.D0 + CVAL*(DENAVP/(1.D0+X) + DENAVM*DBD23/(1.D0+X*DBD2))
       CHIN = -CVAL*(DENAVP*QP/(1.D0+X) + DENAVM*QM*DBD2/(1.D0+X*DBD2))
       CHI = CHIN/CHID
       XI1 = (QP + CHI)/(1+X)
       XI2 = (QM + CHI*DBD22)/(1+X*DBD2)
       DGAMMA = DENAVP*XI1*XI1 + DENAVM*XI2*XI2
       GAMMA = X*X - PI*ESTAR*ESTAR*DGAMMA
C
      RETURN
      END
C
        DOUBLE PRECISION FUNCTION BRENT(FUNC,X1,X2,TOL)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EXTERNAL FUNC
C       PROGRAM TO FIND THE ROOT OF FUNC KNOWN TO LIE BETWEEN
C       X1 AND X2
C       PARAPHRASED FROM `NUMERICAL RECIPES', W. H. PRESS ET AL.
        ITMAX = 100     
        EPS = 3.D-8
        A = X1
        B = X2
        FA = FUNC(A)
        FB = FUNC(B)
        FC = FB
        DO 11 ITER = 1, ITMAX
         IF(FB*FC.GT.0.0)THEN
          C = A
          FC = FA
          D = B-A
          E = D
         ENDIF
         IF(DABS(FC).LT.DABS(FB))THEN
          A = B
          B = C
          C = A
          FA = FB
          FB = FC
          FC = FA
         ENDIF
         TOL1 = 2.*EPS*DABS(B)+0.5*TOL
         XM  = 0.5*(C-B)
         IF(DABS(XM).LE.TOL1.OR.FB.EQ.0.0)THEN
          BRENT = B
          RETURN
         ENDIF
         IF(DABS(E).GE.TOL1.AND.DABS(FA).GT.DABS(FB))THEN
          S = FB/FA
          IF(A.EQ.C)THEN
           P = 2.*XM*S
           Q = 1. -S
          ELSE
           Q = FA/FC
           R = FB/FC
           P = S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
           Q = (Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q = -Q
          P = DABS(P)
          IF(2.*P.LT.MIN(3.*XM*Q-DABS(TOL1*Q),DABS(E*Q)))THEN
           E = D
           D = P/Q
          ELSE
           D = XM
           E = D
          ENDIF
         ELSE
          D = XM
          E = D
         ENDIF
         A = B
         FA = FB
         IF(DABS(D).GT.TOL1)THEN
          B = B + D
         ELSE
          B = B + SIGN(TOL1,XM)
         ENDIF
         FB = FUNC(B)
11      CONTINUE
        WRITE(6,*)'BRENT: EXCEEDED MAXIMUM ITERATIONS'
        BRENT = B
        RETURN
        END             
