*** SUPCRTBL - Calculates the standard molal thermodynamic properties 
***            of reactions among minerals, gases, and aqueous species
***            using equations and data given by Helgeson et al. (1978),
***            Tanger and Helgeson (1988), Shock and Helgeson 
***            (1988, 1990), Shock et al. (1989, 1991), Johnson and
***            Norton (1991), Johnson et al. (1991), Sverjensky
***            et al. (1991), and Holland and Powell (2011). 
*** 
************************************************************************
***
*** Author:     James W. Johnson
***             Earth Sciences Department, L-219
***             Lawrence Livermore National Laboratory
***             Livermore, CA 94550
***             johnson@s05.es.llnl.gov
***
*** Revised 2015:   
***             Kurt Zimmer(1), Yilun Zhang(2), Chen Zhu(2), etc.
***             (1) Department of Computer Science 
***             (2) Department of Geology
***             Indiana University
***             Contact: chenzhu@indiana.edu
***
************************************************************************

      PROGRAM supcrt
 
      PARAMETER (NPLOTF = 8)

      LOGICAL wetrun, unirun
      INTEGER reac, rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF)

      COMMON /io/ rterm, wterm, iconf, reacf, pronf, tabf, plotf

      SAVE

      CALL banner

      CALL readin(nreac,wetrun,unirun) 

      WRITE(wterm,10)
 10   FORMAT(/,' execution in progress ... ',/)

      IF (wetrun) THEN 
           WRITE(wterm,20)
 20        FORMAT(' calculating H2O properties ...',/)
           CALL getH2O(unirun)
      END IF

      CALL tabtop

      DO 30  reac = 1,nreac
           WRITE(wterm,40) reac, nreac
 40        FORMAT(' calculating properties for reaction ',i2,
     1            ' of ',i2,' ...')

           CALL getmgi(reac)

           CALL wrtrxn(reac)

 30        CALL runrxn(reac,wetrun)

      WRITE(wterm,50)
 50   FORMAT(/,' ... execution completed.',/)

      END

********************************************************************

*** banner - Write program banner to the terminal screen.

      SUBROUTINE banner

      PARAMETER (NPLOTF = 8)

      INTEGER rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF)

      COMMON /io/ rterm, wterm, iconf, reacf, pronf, tabf, plotf

      SAVE


      WRITE(wterm,10)
 10   FORMAT(/,5x,' Welcome to SUPCRTBL',
     1       /,5x,' Originally Developed By:',
     1       /,5x,' James W. Johnson,',
     1       ' Eric Oelkers,',
     1       ' Harold Helgeson',
     1       /,5x,' On 13 November 1991',
     1       /,5x,' Updated by Kurt Zimmer, Yilun Zhang and Chen Zhu',
     2       /,5x,' in 2016 in Bloomington ("BL"), Indiana, USA',/)

      END 

********************************************************************

*** consts - Constants

      BLOCK DATA consts

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXISO = 11, MAXINC = 21, NPLOTF = 8)

      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     1         tempf, mapiso(2,3), mapinc(2,3), mapv3(2,3),
     2         rec1m1, rec1m2, rec1m3, rec1m4, rec1aa, rec1gg 

      DOUBLE PRECISION  mwH2O, satmin(2)
      DOUBLE PRECISION  dsvar(MAXINC,MAXISO), Vw(MAXINC,MAXISO),
     1                  bew(MAXINC,MAXISO), alw(MAXINC,MAXISO),
     2                  dalw(MAXINC,MAXISO), Sw(MAXINC,MAXISO),
     3                  Cpw(MAXINC,MAXISO), Hw(MAXINC,MAXISO),
     4                  Gw(MAXINC,MAXISO), Zw(MAXINC,MAXISO),
     5                  Qw(MAXINC,MAXISO), Yw(MAXINC,MAXISO),
     6                  Xw(MAXINC,MAXISO) 

      LOGICAL lvdome(MAXINC,MAXISO), H2Oerr(MAXINC,MAXISO),
     1        EQ3run, lv1bar

      CHARACTER*4  incvar(2,3)
      CHARACTER*10 isov(2,3), incv(2,3), var3(2,3), isosat(2)
      CHARACTER*12 isovar(2,3)
      CHARACTER*20 namecf, namerf, nametf, namepf(NPLOTF)

      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /io2/    tempf
      COMMON /stvars/ isosat, isovar, incvar
      COMMON /TPDmap/ mapiso, mapinc, mapv3
      COMMON /headmp/ isov, incv, var3
      COMMON /refval/ mwH2O, R, Pref, Tref, ZPrTr, YPrTr
      COMMON /rlimit/ nmin1,  nmin2,  nmin3,  nmin4,  ngas,   naqs,
     1                rec1m1, rec1m2, rec1m3, rec1m4, rec1gg, rec1aa
      COMMON /tranm2/ ntrm2
      COMMON /aqscon/ eta, theta, psi, anion, cation, gref
      COMMON /qtzcon/ aa, ba, ca, VPtTta, VPrTtb, Stran
      COMMON /satend/ satmin
      COMMON /defval/ DPMIN,  DPMAX,  DPINC, DTMIN, DTMAX, DTINC,
     1                DTSMIN, DTSMAX, DTSINC 
      COMMON /null/   XNULLM, XNULLA
      COMMON /badtd/  lvdome, H2Oerr
      COMMON /fnames/ namecf, namerf, nametf, namepf
      COMMON /EQ36/   EQ3run
      COMMON /lv1b/   lv1bar
      COMMON /H2Ogrd/ dsvar, Vw, bew, alw, dalw, Sw, Cpw, Hw, Gw, 
     1                Zw, Qw, Yw, Xw
      COMMON /H2Oss/  Dwss, Vwss, bewss, alwss, dalwss, Swss,
     1                Cpwss, Hwss, Gwss, Zwss, Qwss, Ywss, Xwss


      SAVE

      DATA EQ3run, lv1bar / 2*.FALSE. /

***   8 = NPLOTF
      DATA namepf / 8*'                    ' /

***   13*MAXISO*MAXINC = 3003
      DATA dsvar, Vw, bew, alw, dalw, Sw, Cpw, Hw, Gw, Zw, Qw, Yw, Xw 
     1     / 3003*0.0d0 /

***   2*MAXISO*MAXINC = 462 
      DATA lvdome, H2Oerr / 462*.FALSE. /

      DATA Dwss, Vwss, bewss, alwss, dalwss, Swss,
     1     Cpwss, Hwss, Gwss, Zwss, Qwss, Ywss, Xwss / 13*0.0d0 /

      DATA XNULLM, XNULLA / 999999.0d0, 999.0d0 /
      DATA DPMIN,  DPMAX,  DPINC  / 500.0d0, 5000.0d0, 500.0d0 /
      DATA DTMIN,  DTMAX,  DTINC  /   0.0d0, 1000.0d0, 100.0d0 /
      DATA DTSMIN, DTSMAX, DTSINC /   0.0d0,  350.0d0,  25.0d0 /

      DATA satmin / 0.01d0, 0.006117316772d0 /

      DATA aa, ba, ca / 0.549824d3,  0.65995d0, -0.4973d-4 /
      DATA VPtTta, VPrTtb, Stran / 0.23348d2, 0.2372d2, 0.342d0 /
      DATA eta, theta, psi / 0.166027d6, 0.228d3, 0.26d4 /
      DATA anion, cation, gref / 0.0d0,      0.94d0,  0.0d0 /

      DATA mwH2O, R     / 18.0152d0, 1.9872d0 /
      DATA Pref, Tref   /  0.1d1, 0.29815d3 /
***   ZPrTr, YPrTr calculated in SUBR getH2O as f(epseqn)

      DATA rterm, wterm, iconf, reacf, pronf, tabf, tempf
     1     / 5,     6,     40,    41,    42,   43,    44 /

***   8 = NPLOTF
      DATA plotf / 61, 62, 63, 64, 65, 66, 67, 68 /

      DATA isovar / 'CHORES(g/cc)', 'BARS(bars)  ', 3*'THERMS(degC)',
     1              'BARS(bars)  ' /
      DATA incvar / 2*'TEMP', 'DENS'  , 2*'PRES',   'TEMP' /
      DATA isosat / 'TEMP(degC)',   'PRES(bars)' /
   
      DATA isov   / 'DH2O(g/cc)', 'PRES(bars)', 3*'TEMP(degC)',
     1              'PRES(bars)' /

      DATA incv   / 2*'TEMP(degC)', 'DH2O(g/cc)', 'PRES(bars)',
     1                'PRES(bars)', 'TEMP(degC)' /

      DATA var3   / 'PRES(bars)', 'DH2O(g/cc)', 'PRES(bars)',
     1              3*'DH2O(g/cc)' /

      DATA mapiso / 3, 2, 1, 1, 1, 2 /
      DATA mapinc / 1, 1, 3, 2, 2, 1 /
      DATA mapv3  / 2, 3, 2, 3, 4, 4 /

      END

************************************************************************

*** readin - Open user-specified, direct-access data file (STOP if none
***          can be located); open/read or create[/store] an input file 
***          containing i/o specifications and state conditions; 
***          open/read line 1 of an existing file containing reaction 
***          titles and stoichiometries or create[/store] such a file
***          in its entirety.

      SUBROUTINE readin(nreac,wetrun,unirun)

      LOGICAL wetrun, wetcon, wetrxn, unirun, getdf

      SAVE


      IF (.NOT. getdf()) STOP
      CALL getcon(wetcon,unirun)
      CALL getrxn(nreac,wetrxn)
      wetrun = (wetcon) .OR. (wetrxn)
      CALL getout

      END

************************************************************************

*** getdf - Returns .TRUE. if an appropriate direct-access
***         data file can be opened; otherwise returns .FALSE.

      LOGICAL FUNCTION getdf()

      PARAMETER (MAXTRY = 5, NPLOTF = 8)

      CHARACTER*1   ans
      CHARACTER*20  pfname
      LOGICAL  openf
      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF), try,
     1         rec1m1, rec1m2, rec1m3, rec1m4, rec1gg, rec1aa

      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /rlimit/ nmin1,  nmin2,  nmin3,  nmin4,  ngas,   naqs,
     1                rec1m1, rec1m2, rec1m3, rec1m4, rec1gg, rec1aa
      COMMON /dapron/ pfname

      SAVE


   1  WRITE(wterm,10) 
  10  FORMAT(/,' would you like to use the default thermodynamic'
     1        ,' database? (y/n)',/)
      READ(rterm,20) ans
  20  FORMAT(a1)
      IF ((ans .NE. 'y') .AND. (ans .NE. 'Y') .AND.
     1    (ans .NE. 'n') .AND. (ans .NE. 'N')) THEN
           GO TO 1
      END IF
      IF ((ans .EQ. 'y') .OR. (ans .EQ. 'Y')) THEN
           pfname = 'dpronsbl.dat'
      ELSE
           WRITE(wterm,30) 
  30       FORMAT(/,' specify filename for thermodynamic database: ',/)
           READ(rterm,40) pfname
  40       FORMAT(a20)
      END IF

      try = 0

*JP*  2   IF (openf(wterm,pronf,pfname,1,2,1,90)) THEN
  2   IF (openf(wterm,pronf,pfname,1,2,1,114)) THEN
           READ(pronf,50,REC=1) nmin1, nmin2, nmin3, nmin4, 
     1                          ngas, naqs
  50       FORMAT(6(1x,i4))
           READ(pronf,50,REC=2) rec1m1, rec1m2, rec1m3, rec1m4, 
     1                          rec1gg, rec1aa
           getdf = .TRUE. 
           RETURN
      ELSE
           try = try + 1
           IF (try .LT. MAXTRY) THEN
***             prompt for alternative file ***
                WRITE(wterm,60) pfname
  60            FORMAT(/,' Cannot find ',a20,
     1                 /,' enter correct filename: ',/)
                READ(rterm,40) pfname
                GO TO 2
           ELSE
***             give up ***
                WRITE(wterm,70)
  70            FORMAT(/,' I am tired of looking for this file;',
     1                 /,' please do the legwork yourself!',
     2                //,' Bye for now ...',/) 
                getdf = .FALSE.
                RETURN
           END IF
      END IF

      END

************************************************************************

*** getcon - Open/read or create[/store] an input (CON) file that 
***          contains i/o specifications and state conditions.  

      SUBROUTINE getcon(wetcon,unirun)

      PARAMETER (NPLOTF = 8)

      LOGICAL wetcon, unirun
      INTEGER rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF)

      COMMON /io/ rterm, wterm, iconf, reacf, pronf, tabf, plotf

      SAVE


   1  WRITE(wterm,10)
  10  FORMAT(/,' choose file option for specifying',
     1         ' reaction-independent parameters: ',
     1       /,'      1 = select one of three default files',
     2       /,'      2 = select an existing non-default file',
     3       /,'      3 = build a new file:',/)
      READ(rterm,*) ifopt
      IF ((ifopt .LT. 1) .OR. (ifopt .GT. 3)) GO TO 1

      IF (ifopt .EQ. 1) THEN
	   unirun = .FALSE.
	   CALL defaul(wetcon)
	   RETURN
      END IF

      IF (ifopt .EQ. 2) THEN
           CALL readcf(wetcon,unirun)
      ELSE
           CALL makecf(wetcon,unirun)
      END IF

      END

************************************************************************

*** defaul - Set default options / state conditions.

      SUBROUTINE defaul(wetcon)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXODD = 21, NPLOTF = 8)

      DOUBLE PRECISION  isomin, isomax, isoinc, Kmin, Kmax, Kinc,
     1                  oddv1(MAXODD), oddv2(MAXODD) 
      INTEGER           rterm, wterm, reacf, pronf, tabf, 
     1                  plotf(NPLOTF), univar, useLVS, epseqn, geqn

      LOGICAL wetcon, EQ3run, savecf, saverf

      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /grid/   isomin, isomax, isoinc, v2min, v2max, v2inc,
     1                oddv1, oddv2, Kmin, Kmax, Kinc, niso, nv2, nlogK
      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /defval/ DPMIN,  DPMAX,  DPINC, DTMIN, DTMAX, DTINC,
     1                DTSMIN, DTSMAX, DTSINC 
      COMMON /EQ36/   EQ3run
      COMMON /saveif/ savecf, saverf

      SAVE


      univar = 0
      noninc = 0
      useLVS = 1
      epseqn = 4
      geqn   = 3
      EQ3run = .FALSE.
      savecf = .FALSE.

***** prompt for / read isat *****

  1   WRITE(wterm,10) 
 10   FORMAT(/,' input solvent phase region',
     1       /,'      1 = one-phase region ',
     2       /,'      2 = liq-vap saturation curve:',
     3       /,'      3 = EQ3/6 one-phase/sat grid:',/)
      READ(rterm,*) isat
      IF ((isat .LT. 1) .OR. (isat .GT. 3)) THEN
           GO TO 1
      ELSE
           isat = isat - 1 
           wetcon = (isat .EQ. 1) 
      END IF

      IF (isat .EQ. 0) THEN
           iopt   = 2
           iplot  = 1 
           isomin = DPMIN 
           isomax = DPMAX
           isoinc = DPINC
           niso   = 1 + NINT((isomax - isomin)/isoinc)
           v2min  = DTMIN
           v2max  = DTMAX
	   v2inc  = DTINC
           nv2    = 1 + NINT((v2max - v2min)/v2inc)
           RETURN
      END IF

      IF (isat .EQ. 1) THEN
           iopt   = 1
           iplot  = 3
           v2min  = DTSMIN 
           v2max  = DTSMAX
           v2inc  = DTSINC
           nv2    = 1 + NINT((v2max - v2min)/v2inc)
           isomin = 0.0d0
           isomax = 0.0d0
	   isoinc = 0.0d0
           niso   = 1
           RETURN
      END IF

      IF (isat .EQ. 2) THEN
           isat   = 0
           iopt   = 2
           iplot  = 2
           niso   = 0
           nv2    = 0
           noninc = 8
           EQ3run = .TRUE.

           oddv1(1) =   0.01d0
           oddv1(2) =  25.00d0
           oddv1(3) =  60.00d0
           oddv1(4) = 100.00d0
           oddv1(5) = 150.00d0
           oddv1(6) = 200.00d0
           oddv1(7) = 250.00d0
           oddv1(8) = 300.00d0

           oddv2(1) =  1.01322d0
           oddv2(2) =  1.01322d0
           oddv2(3) =  1.01322d0
           oddv2(4) =  1.01322d0
           oddv2(5) =  4.75717d0
           oddv2(6) = 15.53650d0
           oddv2(7) = 39.73649d0
           oddv2(8) = 85.83784d0

           RETURN
      END IF

      END

************************************************************************

*** readcf - Read options / state conditions (CON) file.

      SUBROUTINE readcf(wetcon,unirun)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXISO = 11, MAXINC = 21, MAXODD = 21, 
     1           NPLOTF = 8, TOL = 1.0d-6)

      CHARACTER*1   TP(2)
      CHARACTER*20  namecf, namerf, nametf, namepf(NPLOTF)
      LOGICAL  openf, wetcon, unirun, savecf, saverf
      DOUBLE PRECISION  isomin, isomax, isoinc, Kmin, Kmax, Kinc,
     1                  oddv1(MAXODD), oddv2(MAXODD) 
      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     1         univar, useLVS, epseqn, geqn


      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /grid/   isomin, isomax, isoinc, v2min, v2max, v2inc,
     1                oddv1, oddv2, Kmin, Kmax, Kinc, niso, nv2, nlogK
      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /fnames/ namecf, namerf, nametf, namepf
      COMMON /saveif/ savecf, saverf

      SAVE

      DATA TP / 'T', 'P' /


  1   WRITE(wterm,10)
 10   FORMAT(/,' specify file name:',/)
      READ(rterm,20) namecf
 20   FORMAT(a20)
      IF (.NOT. openf(wterm,iconf,namecf,1,1,1,132)) GO TO 1

      savecf = .TRUE.

***** skip first 4 comment lines
      READ(iconf,21)
 21   FORMAT(///) 

**********************************************************
*** READ and hardwire statements for distribution version

      READ(iconf,*) isat, iopt, iplot, univar, noninc
      useLVS = 1
      epseqn = 4
      geqn   = 3

*** READ statement for development version

*     READ(iconf,*) isat, iopt, iplot, univar, noninc,
*    1              useLVS, epseqn, geqn
**********************************************************

*** insert validity checker for /icon/
*** variables here if desired later

      wetcon = (isat .EQ. 1) .OR. (iopt .EQ. 1) 
      unirun = (univar .EQ. 1)

      IF (noninc .NE. 0) THEN 
***        univar = 0
***        read non-incremental state conditions 
           IF (noninc .GT. MAXODD) THEN
                WRITE(wterm,131) noninc, MAXODD
                WRITE(tabf,131) noninc, MAXODD
 131            FORMAT(/,' Number of specified odd-increment pairs '
     1                  ,'(',i3,') exceeds MAXODD (',i3,').',
     2                 /,' Revise specifications.')
                STOP
           END IF
           DO 30 i = 1,noninc
                IF (isat .EQ. 1) THEN
                     READ(iconf,*) oddv1(i)
                ELSE
                     READ(iconf,*) oddv1(i), oddv2(i)
                END IF
 30             CONTINUE
           RETURN
      END IF

      IF (isat .EQ. 0) THEN
           READ(iconf,*) isomin, isomax, isoinc
           IF (isomin .EQ. isomax) THEN
                niso = 1
           ELSE
                IF (isoinc .EQ. 0.0d0) THEN
                     WRITE(wterm,935)
                     WRITE(tabf,935)
 935                 FORMAT(/,' Ill-defined ',
     1               ' min,max,increment  trio',/,
     2               ' Revise specifications.')
                     STOP
                END IF
                fpniso = 1.0d0 + ((isomax - isomin)/isoinc)
                niso   = NINT(fpniso)
                IF (DABS(fpniso-DBLE(niso)) .GT. TOL) THEN
                     WRITE(wterm,935)
                     WRITE(tabf,935)
                     STOP
                END IF
           END IF
           IF (niso .GT. MAXISO) THEN
                WRITE(wterm,31) niso, MAXISO
                WRITE(tabf,31) niso, MAXISO
 31             FORMAT(/,' Number of specified isopleths (',i4,')'
     1                  ,' exceeds MAXISO (',i3,').',
     2                 /,' Revise specifications.')
                STOP
           END IF
      ELSE
           READ(iconf,*) v2min, v2max, v2inc
           IF (v2min .EQ. v2max) THEN
                nv2 = 1
           ELSE
                IF (v2inc .EQ. 0.0d0) THEN
                     WRITE(wterm,935)
                     WRITE(tabf,935)
                     STOP
                END IF
                fpnv2 = 1.0d0 + ((v2max - v2min)/v2inc)
                nv2   = NINT(fpnv2)
                IF (DABS(fpnv2-DBLE(nv2)) .GT. TOL) THEN
                     WRITE(wterm,935)
                     WRITE(tabf,935)
                     STOP
                END IF
           END IF
           IF (nv2 .GT. MAXINC) THEN
                WRITE(wterm,32) nv2, MAXINC
                WRITE(tabf,32) nv2, MAXINC
 32             FORMAT(/,' Number of specified increments '
     1                  ,'(',i3,') exceeds MAXINC (',i3,').',
     2                 /,' Revise specifications.')
                STOP
           END IF
           niso = 1
           isomin = 0.0d0
           isomax = 0.0d0
           isoinc = 0.0d0
           RETURN
      END IF

      IF (univar .EQ. 1) THEN
***        univariant curve option enabled
           READ(iconf,*) Kmin, Kmax, Kinc
           IF (Kmin .EQ. Kmax) THEN
                nlogK = 1
           ELSE
                IF (Kinc .EQ. 0.0d0) THEN
                     WRITE(wterm,935)
                     WRITE(tabf,935)
                     STOP
                END IF
                fplK  = 1.0d0 + ((Kmax - Kmin)/Kinc)
                nlogK = NINT(fplK)
                IF (DABS(fplK-DBLE(nlogK)) .GT. TOL) THEN
                     WRITE(wterm,935)
                     WRITE(tabf,935)
                     STOP
                END IF
           END IF
           READ(iconf,*) v2min, v2max
           IF (v2min .LT. v2max) THEN
                v2inc = 0.0d0
           ELSE
                WRITE(wterm,152) TP(iplot), TP(iplot)
 152            FORMAT(/,1x,a1,'min >= ',a1,'max ',
     1                 /,1x,' revise specifications')
                STOP
           END IF
      ELSE
***        univariant curve option disabled
           READ(iconf,*) v2min, v2max, v2inc
           IF (v2min .EQ. v2max) THEN
                nv2 = 1
           ELSE
                IF (v2inc .EQ. 0.0d0) THEN
                     WRITE(wterm,935)
                     WRITE(tabf,935)
                     STOP
                END IF
                fpnv2 = 1.0d0 + ((v2max - v2min)/v2inc)
                nv2   = NINT(fpnv2)
                IF (DABS(fpnv2-DBLE(nv2)) .GT. TOL) THEN
                     WRITE(wterm,935)
                     WRITE(tabf,935)
                     STOP
                END IF
           END IF
           IF (nv2 .GT. MAXINC) THEN
                WRITE(wterm,32) nv2, MAXINC
                WRITE(tabf,32) nv2, MAXINC
                STOP
           END IF
      END IF

      END

************************************************************************

*** makecf - Prompt for and create options / state conditions 
***          (CON) file.

      SUBROUTINE makecf(wetcon,unirun)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXISO = 11, MAXINC = 21, MAXODD = 21, 
     1           NPLOTF = 8, TOL = 1.0d-6)

      CHARACTER*1  ptype2(2), TP(2), ans
      CHARACTER*4  incvar(2,3)
      CHARACTER*6  ptype1(2)
      CHARACTER*10 isov(2,3), incv(2,3), var3(2,3), isosat(2)
      CHARACTER*12 isovar(2,3) 
      CHARACTER*20 namecf, namerf, nametf, namepf(NPLOTF)
      LOGICAL  openf, wetcon, unirun, savecf, saverf
      DOUBLE PRECISION isomin, isomax, isoinc, Kmin, Kmax, Kinc,
     1                 oddv1(MAXODD), oddv2(MAXODD)
      INTEGER rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     1        univar, useLVS, epseqn, geqn

      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /grid/   isomin, isomax, isoinc, v2min, v2max, v2inc,
     1                oddv1, oddv2, Kmin, Kmax, Kinc, niso, nv2, nlogK
      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /stvars/ isosat, isovar, incvar
      COMMON /headmp/ isov, incv, var3
      COMMON /fnames/ namecf, namerf, nametf, namepf
      COMMON /saveif/ savecf, saverf

      SAVE

      DATA TP     / 'T',      'P'      /
      DATA ptype2 / 'D',      'P'      /
      DATA ptype1 / 'CHORIC', 'BARIC ' /


***** prompt for / read isat *****

  1   WRITE(wterm,10) 
 10   FORMAT(/,' specify solvent phase region ',
     1       /,'      1 = one-phase region ',
     2       /,'      2 = liq-vap saturation curve:',/)
      READ(rterm,*) isat
      IF ((isat .NE. 1) .AND. (isat .NE. 2)) THEN
           GO TO 1
      ELSE
           isat = isat - 1 
      END IF

***** prompt for / read iopt *****

  2   IF (isat .EQ. 0) THEN
           WRITE(wterm,20) 
 20        FORMAT(/,' specify independent state variables: ',
     1            /,'      1 = temperature (degC), density[H2O] (g/cc) ',
     2            /,'      2 = temperature (degC), pressure (bars)',/)
      ELSE
           WRITE(wterm,30) 
 30        FORMAT(/,' specify independent liq-vap saturation variable:',
     1            /,'      1 = temperature (degC)',
     2            /,'      2 = pressure (bars)',/)
      END IF 
      READ(rterm,*) iopt
      IF ((iopt .NE. 1) .AND. (iopt .NE. 2)) GO TO 2

      wetcon = (isat .EQ. 1) .OR. (iopt .EQ. 1)

      IF (isat .EQ. 1) THEN
***** saturation curve option enabled *****
***** set univar and iplot *****
           univar = 0
           iplot  = 3
***** prompt for / read noninc *****
  3        WRITE(wterm,40) 
 40        FORMAT(/,' specify table-increment option: ',
     1     /,'      1 = calculate tables having uniform increments',
     2     /,'      2 = calculate tables having unequal increments',/)
           READ(rterm,*) noninc
           IF ((noninc .NE. 1) .AND. (noninc .NE. 2)) THEN
                GO TO 3
           ELSE
                noninc = noninc - 1
                IF (noninc .EQ. 0) THEN
***** prompt for / read state condition range along 
***** the saturation curve curve isopleth 
 444                 WRITE(wterm,50) isosat(iopt)
 50                  FORMAT(/,' specify ',a10, ' range:',/,
     1                        ' min, max, increment:'
     1                      ,/)
                     READ(rterm,*) v2min, v2max, v2inc
                     IF (v2max .GT. 373.917d0) THEN
                          WRITE(wterm,899) v2max
 899                      FORMAT(/,' Maximum saturation temperature ',
     1                           '(',f4.0,') > critical temperature',
     2                           ' (373.917).',
     3                           /,' Revise specifications.')  
                          GO TO 444
                     END IF
                     IF (v2min .EQ. v2max) THEN
                          nv2 = 1
                     ELSE
                          IF (v2inc .EQ. 0.0d0) THEN
                               WRITE(wterm,935)
 935                           FORMAT(/,' Ill-defined ',
     1                                  ' min,max,increment  trio',/,
     2                                  ' Revise specifications.')
                               GO TO 444
                          END IF
                          fpnv2 = 1.0d0 + ((v2max - v2min)/v2inc)
                          nv2   = NINT(fpnv2)
                          IF (DABS(fpnv2-DBLE(nv2)) .GT. TOL) THEN
                               WRITE(wterm,935)
                               GO TO 444
                          END IF
                          IF (nv2 .GT. MAXINC) THEN
                               WRITE(wterm,31) nv2, MAXINC
 31                            FORMAT(/,' Number of specified isopleths'
     1                             ,' (',i4,') exceeds MAXINC (',i3,').'
     2                             ,/,' Revise specifications.')
                               GO TO 444
                          END IF
                     END IF
                     niso = 1
                     isomin = 0.0d0
                     isomax = 0.0d0
                     isoinc = 0.0d0
                ELSE
***** prompt for / read [noninc] non-incremental state ***** 
***** condition points along saturation curve          *****
                     WRITE(wterm,60) isosat(iopt)
 60                  FORMAT(/,' specify liq-vap saturation ',a10, 
     1                        ' values', 
     2               /,' one per line, concluding with a zero:',/)
  4                  READ(rterm,*) oddv1(noninc)
                     IF ((oddv1(noninc) .NE. 0.0d0) .AND.
     1                   (noninc .LT. MAXODD)) THEN
                          noninc = noninc + 1
                          GO TO 4
                     END IF
                     IF (oddv1(noninc) .EQ. 0.0d0) THEN
                          noninc = noninc - 1
                     ELSE 
                          WRITE(wterm,241) MAXODD
 241                      FORMAT(/,' Only ',i3,' coordinates separated',
     1                    ' by unequal increments',/,' can be',
     2                    ' processed during one SUPCRTBL execution',/)
                     END IF
                END IF
           END IF
      ELSE
***** saturation curve option curve disabled *****
           IF (iopt .EQ. 1) THEN
                univar = 0
           ELSE
***** prompt for / read univar *****
  5             WRITE(wterm,70) 
 70             FORMAT(/,' would you like to use the univariant curve',
     1                   ' option;',
     2                 /,' i.e., calculate T(logK,P) or P(logK,T) ',
     3                   ' (y/n)',/)
                READ(rterm,75) ans
 75             FORMAT(a1)
                IF ((ans .NE. 'y') .AND. (ans .NE. 'Y') .AND.
     1              (ans .NE. 'n') .AND. (ans .NE. 'N')) THEN
                     GO TO 5
                END IF
                IF ((ans .EQ. 'y') .OR. (ans .EQ. 'Y')) THEN
                     univar = 1
                ELSE
                     univar = 0
                END IF
           END IF
           
           IF (univar .EQ. 0) THEN
***** univariant curve option disabled *****
***** prompt for / read iplot *****
  6             WRITE(wterm,80) ptype1(iopt), ptype2(iopt)
 80             FORMAT(/,' specify tablulation option:',
     1                 /,'      1 = calculate ISO',a6,'(T) tables, ',
     2                 /,'      2 = calculate ISOTHERMAL(',a1,') ',
     3                              'tables ',/)
                READ(rterm,*) iplot
                IF ((iplot .NE. 1) .AND. (iplot .NE. 2)) THEN
                     GO TO 6
                END IF

***** prompt for / read noninc *****
  7             WRITE(wterm,40) 
                READ(rterm,*) noninc
                IF ((noninc .NE. 1) .AND. (noninc .NE. 2)) THEN
                     GO TO 7
                ELSE
                     noninc = noninc - 1
                END IF

                IF (noninc .EQ. 0) THEN
***** prompt for / read state condition ranges in one-phase region *****
 445                 WRITE(wterm,100) isovar(iopt,iplot)
 100                 FORMAT(/,' specify ISO',a12,
     1                      /,' min, max, increment',/)
                     READ(rterm,*) isomin, isomax, isoinc
                     IF (isomin .EQ. isomax) THEN
                          niso = 1
                     ELSE
                          IF (isoinc .EQ. 0.0d0) THEN
                               WRITE(wterm,935)
                               GO TO 445
                          END IF
                          fpniso = 1.0d0 + ((isomax - isomin)/isoinc)
                          niso   = NINT(fpniso)
                          IF (DABS(fpniso-DBLE(niso)) .GT. TOL) THEN
                               WRITE(wterm,935)
                               GO TO 445
                          END IF
                     END IF
                     IF (niso .GT. MAXISO) THEN
                          WRITE(wterm,31) niso, MAXISO
                          GO TO 445
                     END IF

 446                 WRITE(wterm,110) incv(iopt,iplot)
 110                 FORMAT(/,' specify ',a10,' range',
     1                      /,' min, max, increment',/)
                     READ(rterm,*) v2min, v2max, v2inc
                     IF (v2min .EQ. v2max) THEN
                          nv2 = 1 
                     ELSE
                          IF (v2inc .EQ. 0.0d0) THEN
                               WRITE(wterm,935)
                               GO TO 446
                          END IF
                          fpnv2 = 1.0d0 + ((v2max - v2min)/v2inc)
                          nv2    = NINT(fpnv2)
                          IF (DABS(fpnv2-DBLE(nv2)) .GT. TOL) THEN
                               WRITE(wterm,935)
                               GO TO 446
                          END IF
                     END IF
                     IF (nv2 .GT. MAXINC) THEN
                          WRITE(wterm,32) nv2, MAXINC
 32                       FORMAT(/,' Number of specified increments'
     1                           ,' (',i4,') exceeds MAXINC (',i3,').',
     2                           /,' Revise specifications.')
                          GO TO 446
                     END IF
                ELSE
***** prompt for / read [noninc] non-incremental state ***** 
***** condition points in the one-phase region        *****
                     WRITE(wterm,120) isov(iopt,iplot), incv(iopt,iplot)
 120                 FORMAT(/,' specify ',a10,', ',a10,' values; ', 
     1               /,' one pair per line, concluding with 0,0:',/)
  8                  READ(rterm,*) oddv1(noninc), oddv2(noninc)
                     IF ((oddv1(noninc) .NE. 0.0d0) .AND.
     1                   (noninc .LT. MAXODD)) THEN
                          noninc = noninc + 1
                          GO TO 8
                     END IF
                     IF (oddv1(noninc) .EQ. 0.0d0) THEN
                          noninc = noninc - 1
                     ELSE
                          WRITE(wterm,241)
                     END IF
                END IF
           ELSE
***** univariant curve option enabled *****
***** set noninc *****
                noninc = 0                                 
***** prompt for / read iplot *****
  9             WRITE(wterm,130) 
 130            FORMAT(/,' specify univariant calculation option:',
     1                 /,'      1 = calculate T(logK,isobars), ',
     2                 /,'      2 = calculate P(logK,isotherms): ',/)
                READ(rterm,*) iplot
                IF ((iplot .NE. 1) .AND. (iplot .NE. 2)) THEN
                     GO TO 9
                END IF
***** prompt for / read state condition ranges in one-phase region *****
 447            WRITE(wterm,140) isovar(iopt,iplot)
 140            FORMAT(/,' specify ISO',a12,
     1                 /,' min, max, increment ',/)
                READ(rterm,*) isomin, isomax, isoinc
                IF (isomin .EQ. isomax) THEN
                     niso = 1
                ELSE
                     IF (isoinc .EQ. 0.0d0) THEN
                          WRITE(wterm,935)
                          GO TO 447
                     END IF
                     fpniso = 1.0d0 + ((isomax - isomin)/isoinc)
                     niso   = NINT(fpniso)
                     IF (DABS(fpniso-DBLE(niso)) .GT. TOL) THEN
                          WRITE(wterm,935)
                          GO TO 447
                     END IF
                END IF
                IF (niso .GT. MAXISO) THEN
                     WRITE(wterm,31) niso, MAXISO
                     GO TO 447
                END IF
 448            WRITE(wterm,150)
 150            FORMAT(/,' specify logK range: ',
     1                 /,' Kmin, Kmax, Kincrement: ',/)
                READ(rterm,*) Kmin, Kmax, Kinc
                IF (Kmin .EQ. Kmax) THEN
                     nlogK = 1
                ELSE
                     IF (Kinc .EQ. 0.0d0) THEN
                          WRITE(wterm,935)
                          GO TO 448
                     END IF
                     fpnK  = 1.0d0 + ((Kmax - Kmin)/Kinc)
                     nlogK = NINT(fpnK)
                     IF (DABS(fpnK-DBLE(nlogK)) .GT. TOL) THEN
                          WRITE(wterm,935)
                          GO TO 448
                     END IF
                END IF
 449            WRITE(wterm,151) incv(iopt,iplot),
     1                           TP(iplot), TP(iplot)
 151            FORMAT(/,' specify bounding ',a10,' range:',
     1                 /,1x,a1,'min, ',a1,'max: ',/)
                READ(rterm,*) v2min, v2max
                IF (v2min .LT. v2max) THEN
                     v2inc = 0.0d0
                ELSE
                     WRITE(wterm,152) TP(iplot), TP(iplot)
 152                 FORMAT(/,1x,a1,'min >= ',a1,'max ',
     1                      /,1x,' revise specifications')
                     GO TO 449
                END IF
           END IF
      END IF

***************************************************************
*** variable assignments for distribution version
      useLVS = 1
      epseqn = 4
      geqn   = 3 

*** select equation options for development version
*     CALL geteqn(useLVS,epseqn,geqn)
***************************************************************

***** set unirun ******

      unirun = (univar .EQ. 1)

***** write input parameters to new file if desired *****

 16   WRITE(wterm,210)
 210  FORMAT(/,' would you like to save these reaction-independent',
     1       /,' parameters to a file (y/n):',/)
      READ(rterm,75) ans
      IF ((ans .NE. 'y') .AND. (ans .NE. 'Y') .AND.
     1    (ans .NE. 'n') .AND. (ans .NE. 'N')) GO TO 16

      savecf = ((ans .EQ. 'y') .OR. (ans .EQ. 'Y')) 

      IF (savecf) THEN
 17        WRITE(wterm,230)
 230       FORMAT(/,' specify file name:',/)
           READ(rterm,240) namecf
 240       FORMAT(a20)
           IF (.NOT. openf(wterm,iconf,namecf,2,1,1,132)) THEN
                GO TO 17
           END IF

******************************************************************
           WRITE(iconf,250) 

*** statement 250 for distribution versions

 250       FORMAT(' Line 1 (free format):',
     1            ' isat, iopt, iplot, univar, noninc')

*** statement 250 for development versions

*250       FORMAT(' Line 1 (free format): isat, iopt, iplot,',
*    1            ' univar, noninc, useLVS, epseqn, geqn') 
******************************************************************

           IF (noninc .EQ. 0) THEN
                IF (isat .EQ. 1) THEN
                     WRITE(iconf,251)
 251            FORMAT(' Line 2 (free format): v2min, v2max, v2inc')
                     WRITE(iconf,256) 
 256                 FORMAT(66('*'))
                ELSE
                     WRITE(iconf,249)
 249            FORMAT(' Line 2 (free format): isomin, isomax, isoinc')
                END IF
           ELSE
                IF (isat .EQ. 1) THEN
                     WRITE(iconf,252) noninc + 1
 252                 FORMAT(' Lines i=2..',i2,
     1                      ' (free format): oddv1(i)')
                ELSE
                     WRITE(iconf,253) noninc + 1
 253                 FORMAT(' Lines i=2..',i2,
     1                      ' (free format): oddv1(i), oddv2(i)')
                END IF
                WRITE(iconf,256)
           END IF
          
           IF ((isat .EQ. 0) .AND. (noninc .EQ. 0)) THEN
                IF (univar .EQ. 0) THEN
                     WRITE(iconf,254) 
 254                 FORMAT(' Line 3 (free format):',
     1                      ' v2min, v2max, v2inc')
                ELSE
                     WRITE(iconf,255) 
 255                 FORMAT(' Line 3 (free format):',
     1                      ' Kmin, Kmax, Kinc')
                     WRITE(iconf,259) 
 259                 FORMAT(' Line 4 (free format):',
     1                      ' v2min, v2max')
                END IF
           END IF

           IF (univar .EQ. 0) THEN
                WRITE(iconf,256)
           END IF

*************************************************************
*** WRITE statement for distribution version

           WRITE(iconf,350) isat, iopt, iplot, univar, noninc
 350       FORMAT(5(1x,i3))

*** WRITE statement for development version

*          WRITE(iconf,350) isat, iopt, iplot, univar, noninc, 
*    1                      useLVS, epseqn, geqn 
*350       FORMAT(8(1x,i3))
*************************************************************

           IF (noninc .EQ. 0) THEN
                IF (isat .EQ. 1) THEN
                     WRITE(iconf,*) v2min, v2max, v2inc
                ELSE
                     WRITE(iconf,*) isomin, isomax, isoinc
                END IF
           ELSE
                DO 360 i = 1,noninc
                     IF (isat .EQ. 1) THEN
                          WRITE(iconf,*) oddv1(i)
                     ELSE
                          WRITE(iconf,*) oddv1(i), oddv2(i)
                     END IF
 360                 CONTINUE
           END IF
          
           IF ((isat .EQ. 0) .AND. (noninc .EQ. 0)) THEN
                IF (univar .EQ. 0) THEN
                     WRITE(iconf,*) v2min, v2max, v2inc
                ELSE
                     WRITE(iconf,*) Kmin, Kmax, Kinc
                     WRITE(iconf,*) v2min, v2max    
                END IF
           END IF
      END IF

      END

************************************************************************

*** geteqn - prompt for / read useLVS, epseqn, geqn. 

      SUBROUTINE geteqn(useLVS,epseqn,geqn)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NPLOTF = 8)

      CHARACTER*1 ans 

      INTEGER useLVS, epseqn, geqn, rterm, wterm, reacf, 
     1        pronf, tabf, plotf(NPLOTF)

      COMMON /io/ rterm, wterm, iconf, reacf, pronf, tabf, plotf

      SAVE


 11   WRITE(wterm,160) 
 160  FORMAT(/,' would you like to use the Levelt Sengers et al. (1983)'
     1      ,/,' equation of state for H2O in the critical region (y/n)'
     2      ,/)
      READ(rterm,165) ans
 165  FORMAT(a1)
      IF ((ans .NE. 'y') .AND. (ans .NE. 'Y') .AND.
     1    (ans .NE. 'n') .AND. (ans .NE. 'N')) GO TO 11
     
      IF ((ans .EQ. 'y') .OR. (ans .EQ. 'Y')) THEN
           useLVS = 1
      ELSE
           useLVS = 0
      END IF

 12   WRITE(wterm,170)
 170  FORMAT(/,' specify dielectric option: ',
     1       /,'      1 = use Helgeson-Kirkham (1974) equation',
     2       /,'      2 = use Pitzer (1983) equation',
     3       /,'      3 = use Uematsu-Franck (1980) equation',
     4       /,'      4 = use Johnson-Norton (1991) equation', 
     5       /,'      5 = use Archer-Wang (1990) equation',/) 
      READ(rterm,*) epseqn
      IF ((epseqn .LT. 1) .OR. (epseqn .GT. 5)) GO TO 12

 13   WRITE(wterm,180) 
 180  FORMAT(/,' specify g-function option',
     1       /,'      1 = use Tanger-Helgeson (1988) equation',
     2       /,'      2 = use Shock et al. (in prep.) equation',
     3       /,'      3 = use modified Shock et al. equation',/) 
      READ(rterm,*) geqn
      IF ((geqn .LT. 1) .OR. (geqn .GT. 3)) GO TO 13

      END

************************************************************************

*** getrxn - Open and read an existing reaction (RXN) file or 
***          prompt for, create, [and save] a new reaction file.

      SUBROUTINE getrxn(nreac,wetrxn)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NPLOTF = 8)

      LOGICAL  wetrxn
      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF)

      COMMON /io/ rterm, wterm, iconf, reacf, pronf, tabf, plotf

      SAVE


   1  WRITE(wterm,10)
  10  FORMAT(/,' choose file option for specifying reactions ',
     1       /,'      1 = use an existing reaction file',
     2       /,'      2 = build a new reaction file:',/)
      READ(rterm,*) ifopt
      IF ((ifopt .NE. 1) .AND. (ifopt .NE. 2)) GO TO 1

      IF (ifopt .EQ. 1) THEN
           CALL readrf(nreac,wetrxn)
      ELSE
           CALL makerf(nreac,wetrxn)
      END IF

      END

********************************************************************

*** parse - If the first non-blank substring of the input character 
***         string (chrstr) represents a valid integer or 
***         non-exponential floating-point number, parse returns 
***         .TRUE. and converts this first substring into the 
***         corresponding real number (r8num), then transfers the 
***         second such subset into a CHAR*20 variable (name); 
***         otherwise, parse returns .FALSE. 
 
      LOGICAL FUNCTION parse(chrstr,r8num,name)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXLEN = 20) 
 
      CHARACTER*(*) chrstr
      CHARACTER*20  numstr, name
      LOGICAL sign, deci
      INTEGER chrlen, tempf

      COMMON /io2/ tempf

      SAVE


*** calculate length of chrstr ***

      chrlen = LEN(chrstr)

*** read through leading blanks ***

      nblank = 0
      DO 10 i = 1,chrlen
	  IF (chrstr(i:i) .EQ. ' ') THEN
	       nblank = nblank + 1
          ELSE
	       GO TO 2
	  END IF
 10       CONTINUE

*** initialize local variables ***

  2   sign = .FALSE.
      deci = .FALSE.

*** extract numerical string (integer or 
*** non-exponentiated floating-point numbers only) 

      numlen = 0

      DO 20 i = nblank+1,chrlen

	   IF (chrstr(i:i) .EQ. ' ') THEN
	        IF (((numlen .EQ. 1) .AND. (sign .OR. deci)) .OR. 
     1              ((numlen .EQ. 2) .AND. (sign .AND. deci))) THEN
		     parse = .FALSE.
		     RETURN
                ELSE
***             valid integer or non-exponentiated floating-point 
***             number has been read; pad numerical string with blanks; 
***             read numerical numerical character string numstr into 
***             real*8 variable r8num; jump below to read in name.
                     parse = .TRUE.
                     DO 30 j = numlen+1,MAXLEN
                          numstr(j:j) = ' '
 30                       CONTINUE
*** the following CHARACTER-to-DOUBLE PRECISION conversion is acceptable
*** to most compilers ... but not all
*                    READ(numstr,*) r8num
*** hence, portability considerations require use of the following
*** procedure, which is equivalent and universally acceptable
*** ... albeit ugly
                     OPEN(UNIT=tempf,FILE='zero.dat')
                     WRITE(tempf,*) numstr
                     BACKSPACE(tempf)
                     READ(tempf,*) r8num
                     CLOSE(UNIT=tempf)
		     GO TO 3
		END IF
           END IF

           IF ((chrstr(i:i) .EQ. '-') .OR. (chrstr(i:i) .EQ. '+')) THEN
                IF ((.NOT. sign) .AND. (numlen .EQ. 0)) THEN
		     sign = .TRUE.
		     numlen = numlen + 1
		     numstr(numlen:numlen) = chrstr(i:i)
                ELSE
		     parse = .FALSE.
		     RETURN
                END IF
           ELSE IF (chrstr(i:i) .EQ. '.') THEN
                IF (.NOT. deci) THEN
		     deci = .TRUE.
		     numlen = numlen + 1
		     numstr(numlen:numlen) = chrstr(i:i)
                ELSE
		     parse = .FALSE.
		     RETURN
                END IF
           ELSE IF ((chrstr(i:i) .GE. '0') .AND. 
     1               (chrstr(i:i) .LE. '9')) THEN
		     numlen = numlen + 1
		     numstr(numlen:numlen) = chrstr(i:i)
           ELSE
		     parse = .FALSE.
		     RETURN
           END IF

 20        CONTINUE

*** read through blanks that separate the 
*** number string from the name string 

  3   DO 40  name1 = nblank+numlen+1,chrlen
	  IF (chrstr(name1:name1) .NE. ' ') GO TO 4
 40       CONTINUE

*** transfer non-blank substring beginning 
*** at chrstr(name1:name1) into name

  4   j = 0
      DO 50 i = name1,chrlen

	   IF (chrstr(i:i) .NE. ' ') THEN
		j = j + 1
		name(j:j) = chrstr(i:i)
           ELSE
                IF (j .NE. 0) THEN
***             valid non-blank substring has been read into
***             CHAR*20 variable name; pad name with blanks; 
***             return
                     GO TO 5
                ELSE
		     parse = .FALSE.
		     RETURN
                END IF
           END IF
 50        CONTINUE

  5   DO 60 i = j+1,MAXLEN
           name(i:i) = ' ' 
 60        CONTINUE

      RETURN
       
      END  

************************************************************************

*** getout - Prompt for and read names for output files.

      SUBROUTINE getout

      PARAMETER (NPLOTF = 8)

      LOGICAL  openf, EQ3run
      CHARACTER*4   suffx(NPLOTF)
      CHARACTER*13  prefx2
      CHARACTER*16  prefx1
      CHARACTER*20  namecf, namerf, nametf, namepf(NPLOTF)
      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     1         univar, useLVS, epseqn, geqn, xyplot, end

      COMMON /io/   rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /icon/ isat, iopt, iplot, univar, noninc,
     1              useLVS, epseqn, geqn
      COMMON /fnames/ namecf, namerf, nametf, namepf
      COMMON /plottr/ xyplot, end, nplots
      COMMON /EQ36/   EQ3run

      SAVE
 
      DATA suffx / '.kxy', '.gxy', '.hxy', '.sxy', 
     1             '.cxy', '.vxy', '.dxy', '.2xy' /


  1   WRITE(wterm,10)
 10   FORMAT(/,' specify name for tabulated output file:',/)
      READ(rterm,20) nametf
 20   FORMAT(a20)
      IF (.NOT. openf(wterm,tabf,nametf,2,1,1,132)) THEN
           GO TO 1
      END IF

      IF ((noninc .GT. 0) .AND. (.NOT. EQ3run)) THEN
	   xyplot = 0
           RETURN
      ELSE
  2        WRITE(wterm,30)
  30       FORMAT(/,' specify option for x-y plot files:',
     1     /,' logK, G, H, S, Cp, and V of reaction: ',
     1     /,'      1 = do not generate plot files ',
     2     /,'      2 = generate plot files in generic format',
     3     /,'      3 = generate plot files in KaleidaGraph format',/)
           READ(rterm,*) xyplot
           IF ((xyplot .LT. 1) .OR. (xyplot .GT. 3)) THEN 
		GO TO 2
	   ELSE
		xyplot = xyplot - 1
           END IF
      END IF

      IF (xyplot .EQ. 0) RETURN

      IF (xyplot .EQ. 1) THEN
           IF (EQ3run) THEN
                nplots = NPLOTF
           ELSE
                IF (univar .EQ. 1) THEN
                     nplots = 1 
                ELSE
                     nplots = (NPLOTF-1)+isat 
                END IF
           END IF
           IF (univar .EQ. 1) THEN
                WRITE(wterm,35)
 35             FORMAT(/,' specify prefix for name of x-y plot file;',
     1                 /,' suffix will be ".uxy"',/)
           ELSE
                WRITE(wterm,40)
 40             FORMAT(/,' specify prefix for names of x-y plot files;',
     1                 /,' suffix will be ".[d,[2],k,g,h,s,c,v]xy"',/)
           END IF
           READ(rterm,50) prefx1 
 50        FORMAT(a16)
           DO 60 i = 1,LEN(prefx1)
                IF (prefx1(i:i) .EQ. ' ') THEN
                     end = i-1
                     GO TO 65
                END IF
 60             CONTINUE
 65        IF (univar .EQ. 1) THEN
                namepf(1)(1:end) = prefx1(1:end)
                namepf(1)(end+1:end+4) = '.uxy'
           ELSE
                DO 70 i = 1,nplots
                     namepf(i)(1:end) = prefx1(1:end)
                     namepf(i)(end+1:end+4) = suffx(i)
 70                  CONTINUE
           END IF
	   RETURN
      END IF

*** xyplot = 2

      IF ((isat .EQ. 1) .OR. EQ3run) THEN
	   nplots = 1
           WRITE(wterm,80)
 80        FORMAT(/,' specify prefix for names of x-y plot files;',
     1            /,' suffix will be "R#.axy"',/)
           READ(rterm,90) prefx2 
 90        FORMAT(a13)
      ELSE 
	   IF (univar .EQ. 1) THEN
	        nplots = 1
                WRITE(wterm,100)
 100            FORMAT(/,' specify prefix for names of x-y plot files;',
     1                 /,' suffix will be "R#.uxy"',/)
                READ(rterm,90) prefx2
           ELSE 
	        nplots = NPLOTF-1
                WRITE(wterm,110)
 110            FORMAT(/,' specify prefix for names of x-y plot files;',
     1                 /,' suffix will be "R#.[d,[2],k,g,h,s,c,v]xy"',/)
                READ(rterm,90) prefx2 
	   END IF
      END IF

      DO 120 i = 1,LEN(prefx2)
           IF (prefx2(i:i) .EQ. ' ') THEN
                end = i-1
                GO TO 125
           END IF
 120       CONTINUE

 125  IF ((isat .EQ. 1) .OR. EQ3run) THEN
           namepf(1)(1:end) = prefx2(1:end)
           namepf(1)(end+1:end+3) = 'R01'
           namepf(1)(end+4:end+7) = '.axy'
           RETURN
      END IF

      IF (univar .EQ. 1) THEN
           namepf(1)(1:end) = prefx2(1:end)
           namepf(1)(end+1:end+3) = 'R01'
           namepf(1)(end+4:end+7) = '.uxy'
           RETURN
      END IF 

      DO 130 i = 1,nplots
          namepf(i)(1:end) = prefx2(1:end)
          namepf(i)(end+1:end+3) = 'R01'
          namepf(i)(end+4:end+7) = suffx(i)
 130      CONTINUE

      END

************************************************************************

*** getH2O - Calculate/store requisite H2O properties over the
***          user-specified state condition grid.

      SUBROUTINE getH2O(unirun)

      IMPLICIT DOUBLE PRECISION  (a-h,o-z)

      PARAMETER (NPROP2 = 46)

      LOGICAL unirun, error

      INTEGER univar, useLVS, epseqn, geqn, specs(10)

      DOUBLE PRECISION  states(4), props(NPROP2), mwH2O

      COMMON /refval/ mwH2O, R, Pref, Tref, ZPrTr, YPrTr
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn

      SAVE

      DATA specs  / 2,2,2,5,1,0,0,0,0,0 /
      DATA states / 4*0.0d0 /


      specs(8) = useLVS
      specs(9) = epseqn

*****************************************************************
*** assignment of [Z,Y]PrTr to Johnson-Norton (1991) 
*** values for distribution version

      ZPrTr = -0.1278034682d-1
      YPrTr = -0.5798650444d-4 

*** set ZPrTr and YPrTR per espeqn value for development version

*     CALL seteps(Tref-273.15d0,Pref,epseqn,ZPrTr,YPrTr)
*****************************************************************

***** calculate H2O properties at standard state of 25 degC, 1 bar

      states(1) = Tref-273.15d0
      states(2) = Pref
      specs(6) = 0
      specs(7) = 2
      CALL H2O92(specs,states,props,error)
      CALL H2Ostd(states,props)

      IF (unirun) RETURN
 
      IF (noninc .GT. 0) THEN
           CALL oddH2O
           RETURN
      END IF

      IF (isat .EQ. 0) THEN
           CALL oneH2O
      ELSE
           CALL twoH2O
      END IF

      END 

************************************************************************

*** oddH2O - Calculate/store requisite H2O properties over the
***          user-specified set of state conditions.

      SUBROUTINE oddH2O

      IMPLICIT DOUBLE PRECISION  (a-h,o-z)

      PARAMETER (MAXODD = 21, MAXINC = 21, MAXISO = 11, NPROP2 = 46)

      LOGICAL error, lvdome(MAXINC,MAXISO),
     1        H2Oerr(MAXINC,MAXISO), EQ3run
      INTEGER mapiso(2,3), mapinc(2,3), mapv3(2,3),
     1        univar, useLVS, epseqn, geqn, specs(10)

      DOUBLE PRECISION  states(4), props(NPROP2), 
     1                  isomin, isomax, isoinc, Kmin, Kmax, Kinc,
     2                  oddv1(MAXODD), oddv2(MAXODD)

      COMMON /TPDmap/ mapiso, mapinc, mapv3
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /grid/   isomin, isomax, isoinc, v2min, v2max, v2inc,
     1                oddv1, oddv2, Kmin, Kmax, Kinc, niso, nv2, nlogK
      COMMON /badtd/  lvdome, H2Oerr
      COMMON /EQ36/   EQ3run

      SAVE

      DATA specs  / 2,2,2,5,1,0,0,0,0,0 /
      DATA states / 4*0.0d0 /


      specs(6) = isat
      specs(7) = iopt
      specs(8) = useLVS
      specs(9) = epseqn

      DO 30 iodd = 1,noninc
           states(mapiso(iopt,iplot)) = oddv1(iodd)
           IF (isat .EQ. 0) THEN
                states(mapinc(iopt,iplot)) = oddv2(iodd)
           END IF
           CALL H2O92(specs,states,props,error)
           H2Oerr(iodd,1) = error
           IF (.NOT. error) THEN
	        lvdome(iodd,1) = ((iplot .NE. 3) .AND. 
     1                            (specs(6) .EQ. 1) .AND.
     2                            (.NOT. EQ3run))
		IF (lvdome(iodd,1)) THEN
		     specs(6) = 0
		ELSE
                     IF (EQ3run .AND. (specs(6) .EQ. 1)) THEN
                          isat = 1
                          specs(7) = 1
                     END IF
                     CALL H2Osav(iodd,1,states,props)
		END IF 
           END IF
 30        CONTINUE

      IF (EQ3run) isat = 0 
 
      END

************************************************************************

*** oneH2O - Calculate/store requisite H2O properties over the
***          user-specified state condition grid in the 
***          one-phase region.

      SUBROUTINE oneH2O

      IMPLICIT DOUBLE PRECISION  (a-h,o-z)

      PARAMETER (MAXODD = 21, MAXINC = 21, MAXISO = 11, NPROP2 = 46)

      LOGICAL error, lvdome(MAXINC,MAXISO), H2Oerr(MAXINC,MAXISO) 

      INTEGER mapiso(2,3), mapinc(2,3), mapv3(2,3),
     1        univar, useLVS, epseqn, geqn, specs(10)

      DOUBLE PRECISION  states(4), props(NPROP2), 
     1                  isomin, isomax, isoinc, Kmin, Kmax, Kinc,
     2                  oddv1(MAXODD), oddv2(MAXODD) 

      COMMON /TPDmap/ mapiso, mapinc, mapv3
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /grid/   isomin, isomax, isoinc, v2min, v2max, v2inc,
     1                oddv1, oddv2, Kmin, Kmax, Kinc, niso, nv2, nlogK
      COMMON /badtd/  lvdome, H2Oerr

      SAVE

      DATA specs  / 2,2,2,5,1,0,0,0,0,0 /
      DATA states / 4*0.0d0 /


      specs(6) = isat
      specs(7) = iopt
      specs(8) = useLVS
      specs(9) = epseqn

      DO 10 iso = 1,niso
           states(mapiso(iopt,iplot)) = isomin + (iso-1)*isoinc
           DO 10 inc = 1,nv2
                specs(6) = isat
                specs(7) = iopt
                states(mapinc(iopt,iplot)) = v2min + (inc-1)*v2inc
                CALL H2O92(specs,states,props,error)
                H2Oerr(inc,iso) = error
                IF (error) THEN
                     states(mapiso(iopt,iplot)) = isomin + 
     1                                            (iso-1)*isoinc
                ELSE
                     lvdome(inc,iso) = (specs(6) .EQ. 1)
                     IF (lvdome(inc,iso)) THEN
                          specs(6) = 0
                          states(mapiso(iopt,iplot)) = 
     1                           isomin + (iso-1)*isoinc
		     ELSE
                          CALL H2Osav(inc,iso,states,props)
                     END IF 
                END IF
 10             CONTINUE

      END

************************************************************************

*** twoH2O - Calculate/store requisite H2O properties over the
***          user-specified state condition grid along the 
***          vaporization boundary. 

      SUBROUTINE twoH2O

      IMPLICIT DOUBLE PRECISION  (a-h,o-z)

      PARAMETER (MAXODD = 21, MAXINC = 21, MAXISO = 11, NPROP2 = 46)
      PARAMETER (TS1BAR = 99.6324d0)

      LOGICAL error, lvdome(MAXINC,MAXISO),
     1        H2Oerr(MAXINC,MAXISO), lv1bar
      INTEGER mapiso(2,3), mapinc(2,3), mapv3(2,3),
     1        univar, useLVS, epseqn, geqn, specs(10)

      DOUBLE PRECISION  states(4), props(NPROP2), mwH2O, satmin(2),
     1                  isomin, isomax, isoinc, Kmin, Kmax, Kinc,
     2                  oddv1(MAXODD), oddv2(MAXODD)

      COMMON /refval/ mwH2O, R, Pref, Tref, ZPrTr, YPrTr
      COMMON /TPDmap/ mapiso, mapinc, mapv3
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /grid/   isomin, isomax, isoinc, v2min, v2max, v2inc,
     1                oddv1, oddv2, Kmin, Kmax, Kinc, niso, nv2, nlogK
      COMMON /satend/ satmin
      COMMON /badtd/  lvdome, H2Oerr
      COMMON /lv1b/   lv1bar

      SAVE

      DATA specs  / 2,2,2,5,1,0,0,0,0,0 /
      DATA states / 4*0.0d0 /


      specs(6) = isat
      specs(7) = iopt
      specs(8) = useLVS
      specs(9) = epseqn

      lv1bar   = (iopt .EQ. 1) .AND. (v2min .LE. TS1BAR)

      DO 10 inc = 1,nv2
           IF ((inc .EQ. 1) .AND. (v2min. EQ. 0.0d0)) THEN
                states(mapiso(iopt,iplot)) = satmin(iopt)
           ELSE
                states(mapiso(iopt,iplot)) = v2min+(inc-1)*v2inc
           END IF
           IF (lv1bar .AND. (states(mapiso(iopt,iplot)) 
     1         .LE. TS1BAR)) THEN
                isat      = 0
                specs(6)  = 0
                specs(7)  = 2
                states(2) = Pref
           ELSE
                IF (lv1bar) THEN
                     isat = 1
                END IF
                specs(6) = isat
                specs(7) = iopt
           END IF
           CALL H2O92(specs,states,props,error)
           H2Oerr(inc,1) = error
           IF (.NOT. error) CALL H2Osav(inc,1,states,props)
 10        CONTINUE

      END

************************************************************************

*** seteps - Set ZPrTr and YPrTR per espeqn value. 

      SUBROUTINE seteps(TCref,Pref,epseqn,ZPrTr,YPrTr)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      LOGICAL           error
      INTEGER           epseqn, specs(10)
      DOUBLE PRECISION  states(4), props(46)

      SAVE

      DATA specs / 2,2,2,5,1,0,2,0,0,0 /


      specs(9)  = epseqn
      states(1) = TCref
      states(2) = Pref
      states(3) = 0.0d0

      CALL H2O92(specs,states,props,error)

      ZPrTr = props(37)
      YPrTr = props(39)

      END

************************************************************************

*** H2Ostd - Archive requisite H2O properties for the 
***          standard state of 25 degC, 1 bar.

      SUBROUTINE H2Ostd(states,props)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NPROP2 = 46)

      INTEGER A, G, S, U, H, Cv, Cp, vs, al, be,
     1        di, vi, tc, st, td, Pr, vik, albe,
     2        Z, Y, Q, daldT, X

      DOUBLE PRECISION  states(4), props(NPROP2), mwH2O

      DOUBLE PRECISION  Dwss, Vwss, bewss, alwss, dalwss, Swss,
     1                  Cpwss, Hwss, Gwss, Zwss, Qwss, Ywss, Xwss 

      COMMON /refval/ mwH2O, R, Pref, Tref, ZPrTr, YPrTr

      COMMON /H2Oss/ Dwss, Vwss, bewss, alwss, dalwss, Swss,
     1               Cpwss, Hwss, Gwss, Zwss, Qwss, Ywss, Xwss

      SAVE

      DATA A, G, S, U, H, Cv, Cp, vs, al, be, di, vi,
     1     tc, st, td, Pr, vik, albe, Z, Y, Q, daldT, X
     2   /  1,  3,  5,  7,  9, 11, 13, 15, 17, 19, 21, 23, 
     3     25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45 /


*** archive requisite properties ***

      Dwss   = states(3)
      Vwss   = mwH2O/states(3)
      bewss  = props(be)
      alwss  = props(al)
      dalwss = props(daldT)

      Swss  = props(S)
      Cpwss = props(Cp)
      Hwss  = props(H)
      Gwss  = props(G)

      Zwss = props(Z)
      Qwss = props(Q)
      Ywss = props(Y)
      Xwss = props(X)

      END

************************************************************************

*** H2Osav - Archive requisite H2O properties over the
***          user-specified state condition grid.

      SUBROUTINE H2Osav(row,col,states,props)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NPROP2 = 46, MAXINC = 21, MAXISO = 11)

      LOGICAL EQ3run, lv1bar
      INTEGER row, col
      INTEGER mapiso(2,3), mapinc(2,3), mapv3(2,3),
     1        univar, useLVS, epseqn, geqn 
      INTEGER A, G, S, U, H, Cv, Cp, vs, al, be,
     1        di, vi, tc, st, td, Pr, vik, albe,
     2        Z, Y, Q, daldT, X

      DOUBLE PRECISION  states(4), props(NPROP2), mwH2O

      DOUBLE PRECISION  dsvar(MAXINC,MAXISO), Vw(MAXINC,MAXISO),
     1                  bew(MAXINC,MAXISO), alw(MAXINC,MAXISO),
     2                  dalw(MAXINC,MAXISO), Sw(MAXINC,MAXISO),
     3                  Cpw(MAXINC,MAXISO), Hw(MAXINC,MAXISO),
     4                  Gw(MAXINC,MAXISO), Zw(MAXINC,MAXISO),
     5                  Qw(MAXINC,MAXISO), Yw(MAXINC,MAXISO),
     6                  Xw(MAXINC,MAXISO) 

      COMMON /EQ36/   EQ3run
      COMMON /lv1b/   lv1bar
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /refval/ mwH2O, R, Pref, Tref, ZPrTr, YPrTr
      COMMON /TPDmap/ mapiso, mapinc, mapv3
      COMMON /H2Ogrd/ dsvar, Vw, bew, alw, dalw, Sw, Cpw, Hw, Gw, 
     1                Zw, Qw, Yw, Xw

      SAVE

      DATA A, G, S, U, H, Cv, Cp, vs, al, be, di, vi,
     1     tc, st, td, Pr, vik, albe, Z, Y, Q, daldT, X
     2   /  1,  3,  5,  7,  9, 11, 13, 15, 17, 19, 21, 23, 
     3     25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45 /


*** archive dependent state variables ***

      IF (isat .EQ. 1) THEN
           IF (EQ3run) THEN
                dsvar(row,col) = states(4)
           ELSE
                dsvar(row,col) = states(2/iopt)
           END IF
      ELSE
           IF (lv1bar) THEN
                dsvar(row,col) = states(2)
           ELSE
                dsvar(row,col) = states(mapv3(iopt,iplot))
           END IF
      END IF

*** archive requisite properties ***
      
      Vw(row,col)   = mwH2O/states(3+isat)
      bew(row,col)  = props(be+isat)
      alw(row,col)  = props(al+isat)
      dalw(row,col) = props(daldT+isat)

      Sw(row,col)   = props(S+isat)
      Cpw(row,col)  = props(Cp+isat)
      Hw(row,col)   = props(H+isat)
      Gw(row,col)   = props(G+isat)

      Zw(row,col) = props(Z+isat)
      Qw(row,col) = props(Q+isat)
      Yw(row,col) = props(Y+isat)
      Xw(row,col) = props(X+isat)
 
      END

************************************************************************

*** getmgi - Read standard state properties, equation of state 
***          parameters, and heat capacity coefficients for all
***          mineral, gas, and aqueous species in the current
***          reaction.

      SUBROUTINE getmgi(ireac)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMIN = 10, MAXGAS = 10, MAXAQS = 10, MAXRXN = 3000)

      CHARACTER*80  rtitle(MAXRXN)

      LOGICAL  m2reac(MAXRXN)

      INTEGER  nm(MAXRXN), na(MAXRXN), ng(MAXRXN), nw(MAXRXN),
     1         rec1m(MAXRXN,MAXMIN), rec1a(MAXRXN,MAXAQS), 
     2         rec1g(MAXRXN,MAXGAS)

      DOUBLE PRECISION  coefm(MAXRXN,MAXMIN), coefa(MAXRXN,MAXAQS),
     1                  coefg(MAXRXN,MAXGAS), coefw(MAXRXN)

      COMMON /reac1/ rtitle
      COMMON /reac2/ coefm, coefa, coefg, coefw, nm, na, ng, nw,
     1               rec1m, rec1a, rec1g, m2reac

      SAVE


*** retrieve thermodynamic data for minerals

      DO 10 i = 1,nm(ireac)

           CALL getmin(i,rec1m(ireac,i))

 10        CONTINUE


*** retrieve thermodynamic data for gases

      DO 20  i = 1,ng(ireac)

           CALL getgas(i,rec1g(ireac,i))

 20        CONTINUE

*** retrieve thermodynamic data for aqueous species

      DO 30  i = 1,na(ireac)
           CALL getaqs(i,rec1a(ireac,i))

 30        CONTINUE
      END

************************************************************************

*** getmin - Read, from dprons.dat or an analogous database (starting
***          at record rec1), standard state parameters for the i[th]
***          one-phase mineral species in the current reaction;
***          set ntran(i) to zero.

      SUBROUTINE getmin(i,rec1)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

**JP*      PARAMETER (MAXMIN = 10, MXTRAN = 3, IABC = 3, MAXMK = 4,
      PARAMETER (MAXMIN = 10, MXTRAN = 3, IABC = 4, MAXMK = 4,
     1           NPLOTF = 8)

      INTEGER rec1, ntran(MAXMIN)
      INTEGER rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     1        rec1m1, rec1m2, rec1m3, rec1m4, rec1gg, rec1aa

      CHARACTER*20  mname(MAXMIN)
      CHARACTER*30  mform(MAXMIN)

      DOUBLE PRECISION  Gfmin(MAXMIN), Hfmin(MAXMIN), VPrTrm(MAXMIN),
     1                  SPrTrm(MAXMIN), Alpha(MAXMIN), Kappa(MAXMIN),
     2                  KappaP(MAXMIN), KappaDP(MAXMIN),
     3                  MK1(IABC,MAXMIN), MK2(IABC,MAXMIN),
     4                  MK3(IABC,MAXMIN), MK4(IABC,MAXMIN),
     5                  Ttran(MXTRAN,MAXMIN), Htran(MXTRAN,MAXMIN),
     6                  Vtran(MXTRAN,MAXMIN), dPdTtr(MXTRAN,MAXMIN),
     7                  Tmaxm(MAXMIN), NumAtoms(MAXMIN), Tc(MAXMIN),
     8                  Smax(MAXMIN), Vmax(MAXMIN)

      COMMON /io/ rterm, wterm, iconf, reacf, pronf, tabf, plotf

      COMMON /rlimit/ nmin1,  nmin2,  nmin3,  nmin4,  ngas,   naqs,
     1                rec1m1, rec1m2, rec1m3, rec1m4, rec1gg, rec1aa

      COMMON /mnames/ mname, mform
      COMMON /minref/ Gfmin, Hfmin, SPrTrm, VPrTrm, MK1, MK2, MK3, MK4, 
     1                Ttran, Htran, Vtran, dPdTtr, Tmaxm, ntran, Alpha,
     2                Kappa, KappaP, KappaDP, NumAtoms, Tc, Smax, Vmax

      SAVE

**JP* following statements added to allow additional heat capacity coef d.
**JP* Need to reset to zero when used by minerals that don't use d.

        MK1(4,i)=0.
        MK2(4,i)=0.
        MK3(4,i)=0.
        MK4(4,i)=0.

      IF (rec1 .LT. rec1m2) THEN
           ntran(i) = 0
           GO TO 1
      END IF

      IF (rec1 .LT. rec1m3) THEN 
           ntran(i) = 1
           GO TO 1
      END IF

      IF (rec1 .LT. rec1m4) THEN 
           ntran(i) = 2
           GO TO 1
      END IF

      ntran(i) = 3

  1   READ(pronf,10,REC=rec1)   mname(i), mform(i)
 10   FORMAT(1x,a20,a30)

      READ(pronf,20,REC=rec1+3) Gfmin(i), Hfmin(i), SPrTrm(i), VPrTrm(i)
 20   FORMAT(4x,2(2x,f12.2),2x,f8.2,2x,f8.3)
      Gfmin(i) = Gfmin(i)*1d3/4.184
      Hfmin(i) = Hfmin(i)*1d3/4.184
      SPrTrm(i) = SPrTrm(i)*1/4.184
      VPrTrm(i) = VPrTrm(i)*1d1


**JP*           READ(pronf,30,REC=rec1+4) (MK1(j,i), j=1,3)
**JP* 30        FORMAT(4x,3(2x,f12.6))

      READ(pronf,30,REC=rec1+4) (MK1(j,i), j=1,4)
 30   FORMAT(4x,2(2x,f10.4),2x,f10.2,2x,f10.4)

CMR*** adjust magnitude for coefficient d MK1(4,i), also
**JP* adjustment of MK1(4,i) not neccessary since already adjusted in
**JP* HPRONS.FOR---however MK1(2,i) and MK1(3,i) need to 
**JP* be "deadjusted" in SUPREFIX.FOR (CPB and CPC in that program)
**JP* so as to avoid modifying the code here. Line for MK1(4,i) with 
**JP* coef of 1.0 added anyway to facilitate later modification

****      adjust magnitude for Cp coeffs
      MK1(1,i) = MK1(1,i)*1d3/4.184
      MK1(2,i) = MK1(2,i)*1.0d-3*1d1/4.184
      MK1(3,i) = MK1(3,i)*1.0d5*1d-2/4.184
      MK1(4,i) = MK1(4,i)*1d3/4.184
**JP*      MK1(4,i) = MK1(4,i)*1.0d0
      READ(pronf,35,REC=rec1+5) Alpha(i), Kappa(i), KappaP(i),
     1     KappaDP(i), NumAtoms(i)
 35   FORMAT(8x, f6.2, 2x, f8.1, 2x, f6.2, 2x, f8.4, 2x, f7.3)
      Alpha(i) = Alpha(i)*.00001
 40   FORMAT(8x,f7.2)
      IF (ntran(i) .EQ. 1) THEN
           READ(pronf,45,REC=rec1+6) Tc(i), Smax(i), Vmax(i)
 45        FORMAT(6x, f8.2, 2x, f6.2, 2x, f7.4)
           READ(pronf,40,REC=rec1+7) Tmaxm(i)
      ELSE
           READ(pronf,40,REC=rec1+6) Tmaxm(i)
      END IF


      END

************************************************************************

*** getgas - Read, from dprons.dat or an analogous database (starting
***          at record rec1), standard state parameters for the i[th]
***          gas species in the current reaction.

      SUBROUTINE getgas(i,rec1)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

**JP*      PARAMETER (MAXGAS = 10, IABC = 3, NPLOTF = 8)
      PARAMETER (MAXGAS = 10, IABC = 4, NPLOTF = 8)


      INTEGER  rec1
      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF)

      CHARACTER*20  gname(MAXGAS)
      CHARACTER*30  gform(MAXGAS)

      DOUBLE PRECISION Gfgas(MAXGAS), Hfgas(MAXGAS), VPrTrg(MAXGAS), 
     1                 SPrTrg(MAXGAS), MKg(IABC,MAXGAS), Tmaxg(MAXGAS)


      COMMON /io/ rterm, wterm, iconf, reacf, pronf, tabf, plotf

      COMMON /gnames/ gname, gform 
      COMMON /gasref/ Gfgas, Hfgas, SPrTrg, VPrTrg, MKg, Tmaxg

      SAVE

**JP* following statement added to allow additional heat capacity coef d.
**JP* Need to reset to zero when used by minerals that don't use d.

        MKg(4,i)=0.

      READ(pronf,20,REC=rec1)   gname(i), gform(i)
 20   FORMAT(1x,a20,a30)

      READ(pronf,30,REC=rec1+3) Gfgas(i), Hfgas(i), SPrTrg(i), VPrTrg(i)
 30   FORMAT(4x,2(2x,f12.2),2(2x,f8.2))

      Gfgas(i) = Gfgas(i)*1d3/4.184
      Hfgas(i) = Hfgas(i)*1d3/4.184
      SPrTrg(i) = SPrTrg(i)*1/4.184
      VPrTrg(i) = VPrTrg(i)*1d1

**JP*      READ(pronf,40,REC=rec1+4) MKg(1,i), MKg(2,i), MKg(3,i)
**JP* 40   FORMAT(4x,3(2x,f12.6))

      READ(pronf,40,REC=rec1+4) MKg(1,i), MKg(2,i), MKg(3,i), MKg(4,i)
 40   FORMAT(4x,2(2x,f10.4),2x,f10.2,2x,f10.4)

      READ(pronf,50,REC=rec1+5) Tmaxg(i)
 50   FORMAT(8x,f7.2)

***** adjust magnitude for Cp coeffs

      MKg(1,i) = MKg(1,i)*1d3/4.184
      MKg(2,i) = MKg(2,i)*1.0d-3*1d1/4.184
      MKg(3,i) = MKg(3,i)*1.0d5*1d-2/4.184
      MKg(4,i) = MKg(4,i)*1d3/4.184

      END

************************************************************************

*** getaqs - Read, from dprons.dat or an analogous database (starting 
***          at record rec1), standard state parameters for the i[th] 
***          aqueous species in the current reaction.

      SUBROUTINE getaqs(i,rec1)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXAQS = 10, NPLOTF = 8)

      INTEGER  rec1
      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF)

      CHARACTER*20  aname(MAXAQS)
      CHARACTER*30  aform(MAXAQS)

      DOUBLE PRECISION  Gfaqs(MAXAQS), Hfaqs(MAXAQS), SPrTra(MAXAQS), 
     1                  a(4,MAXAQS), c(2,MAXAQS), 
     2                  wref(MAXAQS), chg(MAXAQS)
 
      COMMON /io/ rterm, wterm, iconf, reacf, pronf, tabf, plotf

      COMMON /anames/ aname, aform
      COMMON /aqsref/ Gfaqs, Hfaqs, SPrTra, c, a, wref, chg

      SAVE


      READ(pronf,20,REC=rec1)   aname(i), aform(i)
 20   FORMAT(1x,a20,a30)

      READ(pronf,30,REC=rec1+3) Gfaqs(i), Hfaqs(i), SPrTra(i)
 30   FORMAT(4x,2(2x,f10.3),2x,f8.3)

      Gfaqs(i) = Gfaqs(i)*1d3/4.184
      Hfaqs(i) = Hfaqs(i)*1d3/4.184
      SPrTra(i) = SPrTra(i)*1/4.184

      READ(pronf,40,REC=rec1+4) a(1,i), a(2,i), a(3,i), a(4,i)
 40   FORMAT(4x,4(2x,f9.4))

      READ(pronf,50,REC=rec1+5) c(1,i), c(2,i), wref(i), chg(i)
 50   FORMAT(4x,3(2x,f9.4),2x,f3.0)
      
***** adjust magnitude for e-o-s coefficients and omega

      a(1,i)  = a(1,i)*1.0d-1*1/4.184
      a(2,i)  = a(2,i)*1.0d2*1/4.184
      a(3,i)  = a(3,i)*1.0d0*1/4.184
      a(4,i)  = a(4,i)*1.0d4*1/4.184
      c(1,i)  = c(1,i)*1.0d0*1/4.184
      c(2,i)  = c(2,i)*1.0d4*1/4.184
      wref(i) = wref(i)*1.0d5*1/4.184

      END

************************************************************************

*** runrxn - Calculate the standard molal thermodynamic properties of
***          the i[th] reaction over the range of user-specified state
***          conditions. 

      SUBROUTINE runrxn(i,wetrun)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INTEGER univar, useLVS, epseqn, geqn
      LOGICAL wetrun

      COMMON /icon/ isat, iopt, iplot, univar, noninc,
     1              useLVS, epseqn, geqn

      SAVE


      IF (univar .EQ. 1) THEN
***** univariant curve option enabled *****
           CALL rununi(i)
           RETURN
      END IF

      IF (noninc .EQ. 0) THEN
***** run orthogonal T-d or T-P grid *****
           CALL rungrd(i,wetrun)
      ELSE
***** run "oddball" T,P or T,d pairs *****
           CALL runodd(i)
      END IF

      END

************************************************************************

*** rungrd - Calculate the standard molal thermodynamic properties of
***          the i[th] reaction over the user-specified
***          state-condition grid. 

      SUBROUTINE rungrd(i,wetrun)

      IMPLICIT DOUBLE PRECISION  (a-h,o-z)

      PARAMETER (MAXMIN = 10, MAXGAS = 10, MAXAQS = 10, MXTRAN = 3, 
     1           MAXINC = 21, MAXISO = 11, MAXODD = 21, MAXRXN = 3000,
     2           NPLOTF = 8)

      LOGICAL  m2reac(MAXRXN), rptran, newiso, wetrun,
     1         lvdome(MAXINC,MAXISO), H2Oerr(MAXINC,MAXISO)

      CHARACTER*80  rtitle(MAXRXN)

      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     1         mapiso(2,3), mapinc(2,3), mapv3(2,3),
     2         univar, useLVS, epseqn, geqn, ptrans(MAXMIN)

      INTEGER  nm(MAXRXN), na(MAXRXN), ng(MAXRXN), nw(MAXRXN),
     1         rec1m(MAXRXN,MAXMIN), rec1a(MAXRXN,MAXAQS), 
     2         rec1g(MAXRXN,MAXGAS)

      DOUBLE PRECISION  coefm(MAXRXN,MAXMIN), coefa(MAXRXN,MAXAQS),
     1                  coefg(MAXRXN,MAXGAS), coefw(MAXRXN)

      DOUBLE PRECISION  isomin, isomax, isoinc, Kmin, Kmax, Kinc,
     1                  logKr, oddv1(MAXODD), oddv2(MAXODD),
     2                  TPD(3), TPDtrn(MAXMIN,MXTRAN,3), 
     3                  mwH2O, satmin(2)

      DOUBLE PRECISION  dsvar(MAXINC,MAXISO), Vw(MAXINC,MAXISO),
     1                  bew(MAXINC,MAXISO), alw(MAXINC,MAXISO),
     2                  dalw(MAXINC,MAXISO), Sw(MAXINC,MAXISO),
     3                  Cpw(MAXINC,MAXISO), Hw(MAXINC,MAXISO),
     4                  Gw(MAXINC,MAXISO), Zw(MAXINC,MAXISO),
     5                  Qw(MAXINC,MAXISO), Yw(MAXINC,MAXISO),
     6                  Xw(MAXINC,MAXISO) 

      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /TPDmap/ mapiso, mapinc, mapv3
      COMMON /refval/ mwH2O, R, Pref, Tref, ZPrTr, YPrTr
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /grid/   isomin, isomax, isoinc, v2min, v2max, v2inc,
     1                oddv1, oddv2, Kmin, Kmax, Kinc, niso, nv2, nlogK
      COMMON /H2Ogrd/ dsvar, Vw, bew, alw, dalw, Sw, Cpw, Hw, Gw, 
     1                Zw, Qw, Yw, Xw
      COMMON /reac1/  rtitle
      COMMON /reac2/  coefm, coefa, coefg, coefw, nm, na, ng, nw,
     1                rec1m, rec1a, rec1g, m2reac
      COMMON /fmeq/   dVr, dSr, dCpr, dHr, dGr, logKr, dlogKT, dlogKP
      COMMON /satend/ satmin
      COMMON /badtd/  lvdome, H2Oerr

      SAVE


      DO 10 iso = 1,niso
           IF (isat .EQ. 0)
     1          TPD(mapiso(iopt,iplot)) = isomin + (iso-1)*isoinc
           DO 20 inc = 1,nv2
                IF (isat .EQ. 0) THEN
                     TPD(mapinc(iopt,iplot)) = v2min + (inc-1)*v2inc
                     TPD(mapv3(iopt,iplot)) = dsvar(inc,iso)
                ELSE
                     IF ((inc .EQ. 1) .AND. (v2min .EQ. 0.0d0)) THEN
                         TPD(mapiso(iopt,iplot)) = satmin(iopt)
                     ELSE
                         TPD(mapiso(iopt,iplot)) = v2min + (inc-1)*v2inc
                     END IF
                     TPD(mapinc(iopt,iplot)) = dsvar(inc,iso)
                     TPD(mapv3(iopt,iplot)-isat) = mwH2O/Vw(inc,iso)
                END IF
                IF (.NOT. (lvdome(inc,iso) .OR. H2Oerr(inc,iso))) THEN
                     CALL reac92(i,TPD(2),TPD(1),TPD(3),Vw(inc,iso),
     1                      bew(inc,iso), alw(inc,iso), dalw(inc,iso),
     2                      Sw(inc,iso), Cpw(inc,iso), Hw(inc,iso), 
     3                      Gw(inc,iso), Zw(inc,iso), Qw(inc,iso),
     4                      Yw(inc,iso), Xw(inc,iso), geqn)
                END IF
                IF (.NOT. m2reac(i)) THEN
                     rptran = .FALSE.
                ELSE
                     newiso = ((inc .EQ. 1) .OR. lvdome(inc-1,iso) .OR.
     1                         H2Oerr(inc-1,iso))
                     CALL m2tran(inc,iso,newiso,nm(i),
     1                           rptran,ptrans,TPD,TPDtrn,wetrun)
                END IF
                CALL report(i,iso,inc,TPD,TPDtrn,rptran,ptrans, 
     1                      dVr,dSr,dCpr,dHr,dGr,logKr, 
     2                      lvdome(inc,iso),H2Oerr(inc,iso),
     3                      .FALSE.)
 20             CONTINUE
 10        CONTINUE    

       END

*******************************************************************

*** m2tran - Returns rptran = .TRUE. if a phase transition occurs 
***          for one or more minerals in the current reaction between 
***          the immediately previous and current state conditions. 

      SUBROUTINE m2tran(inc,iso,newiso,nmreac,rptran,ptrans,TPD,TPDtrn,
     1                  wetrun)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MXTRAN = 3, MAXMIN = 10)

      LOGICAL  rptran, newiso, wetrun

      INTEGER  phaser(MAXMIN), prprev(MAXMIN), ptrans(MAXMIN)

      DOUBLE PRECISION  TPD(3), TPDtrn(MAXMIN,MXTRAN,3)
      DOUBLE PRECISION  Vmin(MAXMIN), Smin(MAXMIN), Cpmin(MAXMIN),
     2                  Hmin(MAXMIN), Gmin(MAXMIN)

      COMMON /minsp/  Vmin, Smin, Cpmin, Hmin, Gmin, phaser

      SAVE
       

      rptran = .FALSE.
      IF (newiso) THEN
           DO 10 imin = 1,nmreac 
                prprev(imin) = phaser(imin)      
                ptrans(imin) = 0
 10             CONTINUE
      ELSE
           DO 20 imin = 1,nmreac
                IF (prprev(imin) .EQ. phaser(imin)) THEN
                     ptrans(imin) = 0
                ELSE
                     rptran = .TRUE.
                     ptrans(imin) = IABS(phaser(imin) - prprev(imin))
                     prprev(imin) = phaser(imin)      
                     CALL getsct(inc,iso,imin,phaser(imin),
     1                           ptrans(imin),TPD,TPDtrn,wetrun)
                END IF
 20             CONTINUE                     
      END IF

      END

*********************************************************************

*** getsct - Get s[tate] c[onditions of phase] t[ransition] 
***          iphase for mineral imin.

      SUBROUTINE getsct(inc,iso,imin,iphase,ntrans,TPD,TPDtrn,wetrun)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

**JP*      PARAMETER (MXTRAN =  3, MAXMIN = 10, IABC = 3, NPROP2 = 46,
      PARAMETER (MXTRAN =  3, MAXMIN = 10, IABC = 4, NPROP2 = 46,
     1           MAXINC = 21, MAXISO = 11, MAXODD = 21, NPLOTF = 8)

      LOGICAL  error, wetrun

      CHARACTER*20  mname(MAXMIN)
      CHARACTER*30  mform(MAXMIN)

      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF)
      INTEGER  univar, useLVS, epseqn, geqn, ntran(MAXMIN),
     1         mapiso(2,3), mapinc(2,3), mapv3(2,3), specs(10)

      DOUBLE PRECISION TPDtrn(MAXMIN,MXTRAN,3), TtranP(MXTRAN,MAXMIN), 
     2                 PtranT(MXTRAN,MAXMIN), states(4), props(NPROP2), 
     3                 isomin, isomax, isoinc, Kmin, Kmax, Kinc,
     4                 oddv1(MAXODD), oddv2(MAXODD), TPD(3)

      DOUBLE PRECISION  Gfmin(MAXMIN), Hfmin(MAXMIN), VPrTrm(MAXMIN),
     1                  SPrTrm(MAXMIN), Alpha(MAXMIN), Kappa(MAXMIN),
     2                  KappaP(MAXMIN), KappaDP(MAXMIN),
     3                  MK1(IABC,MAXMIN), MK2(IABC,MAXMIN),
     4                  MK3(IABC,MAXMIN), MK4(IABC,MAXMIN),
     5                  Ttran(MXTRAN,MAXMIN), Htran(MXTRAN,MAXMIN),
     6                  Vtran(MXTRAN,MAXMIN), dPdTtr(MXTRAN,MAXMIN),
     7                  Tmaxm(MAXMIN), NumAtoms(MAXMIN), Tc(MAXMIN),
     8                  Smax(MAXMIN), Vmax(MAXMIN)

      DOUBLE PRECISION  dsvar(MAXINC,MAXISO), Vw(MAXINC,MAXISO),
     1                  bew(MAXINC,MAXISO), alw(MAXINC,MAXISO),
     2                  dalw(MAXINC,MAXISO), Sw(MAXINC,MAXISO),
     3                  Cpw(MAXINC,MAXISO), Hw(MAXINC,MAXISO),
     4                  Gw(MAXINC,MAXISO), Zw(MAXINC,MAXISO),
     5                  Qw(MAXINC,MAXISO), Yw(MAXINC,MAXISO),
     6                  Xw(MAXINC,MAXISO) 


      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /TPDmap/ mapiso, mapinc, mapv3
      COMMON /PTtran/ TtranP, PtranT
      COMMON /H2Ogrd/ dsvar, Vw, bew, alw, dalw, Sw, Cpw, Hw, Gw, 
     1                Zw, Qw, Yw, Xw
      COMMON /mnames/ mname, mform
      COMMON /minref/ Gfmin, Hfmin, SPrTrm, VPrTrm, MK1, MK2, MK3, MK4, 
     1                Ttran, Htran, Vtran, dPdTtr, Tmaxm, ntran, Alpha,
     2                Kappa, KappaP, KappaDP, NumAtoms, Tc, Smax, Vmax
      COMMON /grid/   isomin, isomax, isoinc, v2min, v2max, v2inc,
     1                oddv1, oddv2, Kmin, Kmax, Kinc, niso, nv2, nlogK

      SAVE

      DATA specs  / 2,2,2,5,1,0,0,0,0,0 /
      DATA states / 4*0.0d0 /
      DATA Tfssat / 139.8888149d0 /

      specs(6) = isat
      specs(7) = iopt
      specs(8) = useLVS
      specs(9) = epseqn      

*** ntrans = # phase transitions for mineral imin between
***          current and last isopleth locations.

      DO 10 itran = ntrans,1,-1
           IF (isat .EQ. 1) THEN
***             vaporization boundary
                IF (mname(imin)(1:11) .EQ. 'FERROSILITE') THEN
                     states(1) = Tfssat
                ELSE
                     states(1) = TtranP(iphase-itran,imin) - 273.15d0
                END IF
                IF (specs(7) .EQ. 2) THEN
                     specs(7) = 1
                END IF
           ELSE
                IF (iplot .EQ. 2) THEN
***                  isotherms(pres or dens)
                     states(1) = TPD(mapiso(iopt,iplot))
                     states(2) = PtranT(iphase-itran+1,imin)
                     IF (specs(7) .EQ. 1) THEN
                          specs(7) = 2
                     END IF
                ELSE
                     IF (iopt .EQ. 2) THEN
***                       isobars(temp)
                          states(1) = TtranP(iphase-itran,imin) 
     1                                - 273.15d0
                          states(2) = TPD(mapiso(iopt,iplot))
                     ELSE
***                       isochores(temp)
                          states(3) = TPD(mapiso(iopt,iplot))
                          IF (dPdTtr(iphase-1,imin) .EQ. 0.0d0) THEN
                               states(1) = TtranP(iphase-itran,imin)  
     1                                     - 273.15d0
                          ELSE
***                            special case, make 
***                            appropriate approximation
                               P1 = dsvar(inc-1,iso)
                               P2 = dsvar(inc,iso)
                               T1 = v2min + (inc-2)*v2inc
                               T2 = v2min + (inc-1)*v2inc
                               states(1) = Tint(P1,P2,T1,T2,
     1                         TtranP(iphase-itran,imin)-273.15d0,
     2                         dPdTtr(iphase-itran,imin))
                          END IF
                     END IF
                END IF
           END IF

           IF (wetrun) THEN
                CALL H2O92(specs,states,props,error)
           ELSE
                error = .FALSE.
                states(isat+3) = 0.0d0
           END IF

           IF (error) THEN
                WRITE(wterm,20) (states(jjj), jjj=1,3)
                WRITE(tabf,20) (states(jjj), jjj=1,3)
 20        format(/,' State conditions fall beyond validity limits of',
     1            /,' the Haar et al. (1984) H2O equation of state:',
     2            /,' T < Tfusion@P; T > 2250 degC; or P > 30kb.',
     3            /,' SUPCRTBL stopped in SUBROUTINE getsct:',
     4           //,' T = ',e12.5,
     5            /,' P = ',e12.5,
     6            /,' D = ',e12.5,/)
                STOP
           ELSE
                TPDtrn(imin,itran,1) = states(1)
                TPDtrn(imin,itran,2) = states(2)
                TPDtrn(imin,itran,3) = states(isat+3)
           END IF

 10        CONTINUE

      END

*********************************************************************

*** Tint - Returns the temperature intersection of isochore(T) 
***        with a mineral phase transition boundary where
***        (dP/dT)tr .NE. 0.  Approximation involves assumption 
***        that (dP/dT)isochore is linear between P1,T1,D
***        and P2,T2,D (consecutive locations on isochore D(T))
***        that bridge the phase transition. 

      DOUBLE PRECISION FUNCTION Tint(P1,P2,T1,T2,TtrnP2,dPdTtr)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      SAVE

      
      bmin  = P2 - dPdTtr*TtrnP2
      dPdTi = (P2 - P1)/(T2 - T1)
      biso  = P2 - dPdTi*T2 


      Tint = (bmin - biso)/(dPdTi - dPdTtr)

      END

************************************************************************

*** runodd - Calculate the standard molal thermodynamic properties of
***          the i[th] reaction over the user-specified set of
***          nonincremental state condition pairs. 

      SUBROUTINE runodd(i)

      IMPLICIT DOUBLE PRECISION  (a-h,o-z)

      PARAMETER (MAXMIN = 10, MAXGAS = 10, MAXAQS =  10, MXTRAN = 3,
     1           MAXINC = 21, MAXISO = 11, MAXODD = 21, NPLOTF = 8)

      LOGICAL  rptdum, lvdome(MAXINC,MAXISO), H2Oerr(MAXINC,MAXISO)

      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     1         mapiso(2,3), mapinc(2,3), mapv3(2,3),
     2         univar, useLVS, epseqn, geqn, ptdumb(MAXMIN)

      DOUBLE PRECISION  isomin, isomax, isoinc, Kmin, Kmax, Kinc,
     1                  logKr, oddv1(MAXODD), oddv2(MAXODD),
     2                  TPD(3), mwH2O

      DOUBLE PRECISION  dsvar(MAXINC,MAXISO), Vw(MAXINC,MAXISO),
     1                  bew(MAXINC,MAXISO), alw(MAXINC,MAXISO),
     2                  dalw(MAXINC,MAXISO), Sw(MAXINC,MAXISO),
     3                  Cpw(MAXINC,MAXISO), Hw(MAXINC,MAXISO),
     4                  Gw(MAXINC,MAXISO), Zw(MAXINC,MAXISO),
     5                  Qw(MAXINC,MAXISO), Yw(MAXINC,MAXISO),
     6                  Xw(MAXINC,MAXISO),
     7                  TPDdum(MAXMIN,MXTRAN,3)

      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /TPDmap/ mapiso, mapinc, mapv3
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /grid/   isomin, isomax, isoinc, v2min, v2max, v2inc,
     1                oddv1, oddv2, Kmin, Kmax, Kinc, niso, nv2, nlogK
      COMMON /refval/ mwH2O, R, Pref, Tref, ZPrTr, YPrTr
      COMMON /H2Ogrd/ dsvar, Vw, bew, alw, dalw, Sw, Cpw, Hw, Gw, 
     1                Zw, Qw, Yw, Xw
      COMMON /fmeq/   dVr, dSr, dCpr, dHr, dGr, logKr, dlogKT, dlogKP
      COMMON /badtd/  lvdome, H2Oerr

      SAVE


*** MAXMIN*MXTRAN*3 = 90 ***
      DATA  TPDdum / 90*0.0d0 /

      DATA  rptdum / .FALSE. /
*** MAXMIN*0
      DATA  ptdumb / 10*0 /


      DO 10 iodd = 1,noninc
           TPD(mapiso(iopt,iplot)) = oddv1(iodd)
           IF (isat .EQ. 0) THEN
                TPD(mapinc(iopt,iplot)) = oddv2(iodd)
                TPD(mapv3(iopt,iplot)) = dsvar(iodd,1)
           ELSE
                TPD(mapinc(iopt,iplot)) = dsvar(iodd,1)
                TPD(mapv3(iopt,iplot)-isat) = mwH2O/Vw(iodd,1)
           END IF
           IF (.NOT. (lvdome(iodd,1) .OR. H2Oerr(iodd,1))) THEN
                CALL reac92(i,TPD(2),TPD(1),TPD(3),Vw(iodd,1),
     1                 bew(iodd,1), alw(iodd,1), dalw(iodd,1),
     2                 Sw(iodd,1), Cpw(iodd,1), Hw(iodd,1), 
     3                 Gw(iodd,1), Zw(iodd,1), Qw(iodd,1),
     4                 Yw(iodd,1), Xw(iodd,1), geqn)
           END IF
           CALL report(i, 1, iodd, TPD, TPDdum, rptdum, ptdumb,
     1                 dVr, dSr, dCpr, dHr, dGr, logKr, 
     2                 lvdome(iodd,1), H2Oerr(iodd,1),.FALSE.)
 10        CONTINUE

      END

************************************************************************

*** rununi - Calculate the standard molal thermodynamic properties of
***          the i[th] reaction over the user-specified set of T,logK
***          or P,logK pairs.

      SUBROUTINE rununi(i)

      IMPLICIT DOUBLE PRECISION  (a-h,o-z)

      PARAMETER (MAXMIN = 10, MAXAQS = 10, MAXGAS = 10, MAXRXN = 3000,
     1           MAXODD = 21, MXTRAN =  3, NPLOTF = 8)

      LOGICAL foundK, Kfound, wetrxn, m2reac(MAXRXN), rptdum

      CHARACTER*80  rtitle(MAXRXN)

      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     1         mapiso(2,3), mapinc(2,3), mapv3(2,3),
     2         univar, useLVS, epseqn, geqn, ptdumb(MAXMIN)
      INTEGER  nm(MAXRXN), na(MAXRXN), ng(MAXRXN), nw(MAXRXN),
     1         rec1m(MAXRXN,MAXMIN), rec1a(MAXRXN,MAXAQS), 
     2         rec1g(MAXRXN,MAXGAS)

      DOUBLE PRECISION  isomin, isomax, isoinc, Kmin, Kmax, Kinc, 
     1                  Kfind, logKr, isoval, oddv1(MAXODD), 
     2                  oddv2(MAXODD), TPD(3), TPDdum(MAXMIN,MXTRAN,3)
      DOUBLE PRECISION  coefm(MAXRXN,MAXMIN), coefa(MAXRXN,MAXAQS),
     1                  coefg(MAXRXN,MAXGAS), coefw(MAXRXN)

      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf 
      COMMON /TPDmap/ mapiso, mapinc, mapv3
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /grid/   isomin, isomax, isoinc, v2min, v2max, v2inc,
     1                oddv1, oddv2, Kmin, Kmax, Kinc, niso, nv2, nlogK
      COMMON /reac1/  rtitle
      COMMON /reac2/  coefm, coefa, coefg, coefw, nm, na, ng, nw,
     1                rec1m, rec1a, rec1g, m2reac
      COMMON /fmeq/   dVr, dSr, dCpr, dHr, dGr, logKr, dlogKT, dlogKP

      SAVE


*** MAXMIN*MXTRAN*3 = 90 ***
      DATA TPDdum / 90*0.0d0 /

      DATA rptdum / .FALSE. /
*** MAXMIN*0
      DATA ptdumb / 10*0 /


      nv2 = nlogK
      wetrxn = ((nw(i) .GT. 0) .OR. (na(i) .GT. 0))
      DO 10 iso = 1,niso
           isoval = isomin + (iso-1)*isoinc
           DO 10 inc = 1,nlogK 
                Kfind = Kmin + (inc-1)*Kinc
                Kfound = foundK(i,wetrxn,Kfind,isoval,
     1                          v2min,v2max,v2val,dH2O)
                IF (.NOT. Kfound) logKr = Kfind  
                TPD(mapiso(iopt,iplot)) = isoval
                TPD(mapinc(iopt,iplot)) = v2val
                TPD(3) = dH2O 
                CALL report(i, iso, inc, TPD, TPDdum, rptdum, 
     1                      ptdumb, dVr, dSr, dCpr, dHr, dGr, logKr, 
     2                      .FALSE.,.FALSE.,Kfound)
 10             CONTINUE
      END

********************************************************************

*** SUBRs report, wrtrxn, wrtssp, report
*** SUBR  blanks

********************************************************************

*** foundK - Returns '.TRUE.' and v2Kfnd[T|P](isoval[P|T],Kfind) if 
***          (1) logK(isoval,var2=v2min..v2max) for the i[th] reaction
***          is unimodal, and (2) logK value Kfind at isoval occurs
***          within v2min..v2max; otherwise returns '.FALSE.'.   
***          v2Kfnd(usival,Kfind) is isolated using a straightforward
***          implementation of the golden section search algorithm 
***          (e.g., Miller (1984), pp. 130-133.)


      LOGICAL FUNCTION foundK(i,wetrxn,Kfind,isoval,v2min,v2max,
     1                        v2val,dH2O)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NPROP2 = 46, MXEVAL = 50, TOL = 1.0d6, NPLOTF = 8)

      LOGICAL wetrxn, error
      INTEGER rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF)
      INTEGER mapiso(2,3), mapinc(2,3), mapv3(2,3),
     1        univar, useLVS, epseqn, geqn, specs(10)
      INTEGER AA, G, S, U, H, Cv, Cp, vs, al, be,
     1        di, vi, tc, st, td, Pr, vik, albe,
     2        Z, Y, Q, daldT, X
      DOUBLE PRECISION isoval, Kfind, major, logKr, mwH2O,
     1                 states(4), props(NPROP2)

      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /refval/ mwH2O, RR, Pref, Tref, ZPrTr, YPrTr
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /TPDmap/ mapiso, mapinc, mapv3
      COMMON /fmeq/   dVr, dSr, dCpr, dHr, dGr, 
     1                logKr, dlogKT, dlogKP

      SAVE

      DATA AA, G, S, U, H, Cv, Cp, vs, al, be, di, vi,
     1     tc, st, td, Pr, vik, albe, Z, Y, Q, daldT, X
     2   /  1,  3,  5,  7,  9, 11, 13, 15, 17, 19, 21, 23, 
     3     25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45 /

      DATA specs / 2,2,2,5,1,0,2,0,0,0 /
      DATA props / 46*0.0d0 /
      DATA error / .FALSE. /


      foundK = .TRUE.
      j = 0
      a = v2min
      b = v2max
      r = (3.0d0 - DSQRT(5.0d0)) / 2.0d0
      major = r * (b-a)
      c = a + major
      d = b - major

*** set acceptance tolerance per TOL
      accept = (1.0d0 + DABS(Kfind))/TOL

      states(mapiso(iopt,iplot)) = isoval
      states(mapinc(iopt,iplot)) = c
      states(3) = 1.0d0
      IF (wetrxn) THEN 
           specs(8) = useLVS
           specs(9) = epseqn
           CALL H2O92(specs,states,props,error)
           IF (error) THEN
                WRITE(wterm,10) (states(jjj), jjj=1,3)
                WRITE(tabf,10) (states(jjj), jjj=1,3)
 10        format(/,' State conditions fall beyond validity limits of',
     1            /,' the Haar et al. (1984) H2O equation of state:',
     2            /,' T < Tfusion@P; T > 2250 degC; or P > 30kb.',
     3            /,' SUPCRTBL stopped in LOGICAL FUNCTION foundK:',
     4           //,' T = ',e12.5,
     5            /,' P = ',e12.5,
     6            /,' D = ',e12.5,/)
                STOP
           END IF
      END IF

      CALL reac92(i,states(2),states(1),states(3),
     1     mwH2O/states(3),props(be),props(al),props(daldT),
     2     props(S),props(Cp),props(H),props(G),props(Z),props(Q),
     3     props(Y),props(X),geqn)
      fc = DABS(logKr - Kfind)

      states(mapiso(iopt,iplot)) = isoval
      states(mapinc(iopt,iplot)) = d
      IF (wetrxn) THEN
           CALL H2O92(specs,states,props,error)
           IF (error) THEN
                WRITE(wterm,10)
                STOP
           END IF
      END IF

      CALL reac92(i,states(2),states(1),states(3),
     1     mwH2O/states(3),props(be),props(al),props(daldT),
     2     props(S),props(Cp),props(H),props(G),props(Z),props(Q),
     3     props(Y),props(X),geqn)
      fd = DABS(logKr - Kfind)

 1    IF (fc .LE. accept) THEN
           states(mapiso(iopt,iplot)) = isoval
           states(mapinc(iopt,iplot)) = c
           IF (wetrxn) THEN
                CALL H2O92(specs,states,props,error)
                IF (error) THEN
                     WRITE(wterm,10)
                     STOP
                END IF
           END IF
           CALL reac92(i,states(2),states(1),states(3),
     1          mwH2O/states(3),props(be),props(al),props(daldT),
     2          props(S),props(Cp),props(H),props(G),props(Z),props(Q),
     3          props(Y),props(X),geqn) 
           v2val = c
           IF (wetrxn) THEN
                dH2O = states(3)
           ELSE
                dH2O = 0.0d0
           END IF
           RETURN
      END IF

      IF (fd .LE. accept) THEN
           states(mapiso(iopt,iplot)) = isoval
           states(mapinc(iopt,iplot)) = d
           IF (wetrxn) THEN
                CALL H2O92(specs,states,props,error)
                IF (error) THEN
                     WRITE(wterm,10)
                     STOP
                END IF
           END IF
           CALL reac92(i,states(2),states(1),states(3),
     1          mwH2O/states(3),props(be),props(al),props(daldT),
     2          props(S),props(Cp),props(H),props(G),props(Z),props(Q),
     3          props(Y),props(X),geqn) 
           v2val = d
           IF (wetrxn) THEN
                dH2O = states(3)
           ELSE
                dH2O = 0.0d0
           END IF
           RETURN
      END IF

      IF (j .GT. MXEVAL) THEN
           foundK = .FALSE.
           IF (wetrxn) THEN
                dH2O = states(3)
           ELSE
                dH2O = 0.0d0
           END IF
           RETURN
      ELSE 
           j = j + 1
      END IF

      IF (fc .LT. fd) THEN
           b = d
           d = c
           fd = fc
           c = a + r*(b-a)
           states(mapiso(iopt,iplot)) = isoval
           states(mapinc(iopt,iplot)) = c
           IF (wetrxn) THEN
                CALL H2O92(specs,states,props,error)
                IF (error) THEN
                     WRITE(wterm,10)
                     STOP
                END IF
           END IF
           CALL reac92(i,states(2),states(1),states(3),
     1          mwH2O/states(3),props(be),props(al),props(daldT),
     2          props(S),props(Cp),props(H),props(G),props(Z),
     3          props(Q),props(Y),props(X),geqn)
           fc = DABS(logKr - Kfind)
      ELSE
           a = c
           c = d
           fc = fd
           d = b - r*(b-a)
           states(mapiso(iopt,iplot)) = isoval
           states(mapinc(iopt,iplot)) = d
           IF (wetrxn) THEN
                CALL H2O92(specs,states,props,error)
                IF (error) THEN
                     WRITE(wterm,10)
                     STOP
                END IF
           END IF 
           CALL reac92(i,states(2),states(1),states(3),
     1          mwH2O/states(3),props(be),props(al),props(daldT),
     2          props(S),props(Cp),props(H),props(G),props(Z),
     3          props(Q),props(Y),props(X),geqn)
           fd = DABS(logKr - Kfind)
      END IF

      GO TO 1

      END

*********************************************************************

*** makerf - Prompt for and create a reaction (RXN) file.

      SUBROUTINE makerf(nreac,wetrxn)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMIN = 10, MAXAQS =  10, MAXGAS =  10, MAXRXN = 3000,
     1           MAXBAD = 10, NPLOTF = 8)

      CHARACTER*1  ans
      CHARACTER*20 namecf, namerf, nametf, namepf(NPLOTF)
      CHARACTER*20 specie, namem(MAXRXN,MAXMIN),
     1             namea(MAXRXN,MAXAQS), nameg(MAXRXN,MAXGAS),
     2             sbad(MAXBAD)
      CHARACTER*30 form, formm(MAXRXN,MAXMIN), formg(MAXRXN,MAXGAS),
     1             forma(MAXRXN,MAXAQS)
      CHARACTER*80 rtitle(MAXRXN), string
      LOGICAL  openf, wetrxn, m2reac(MAXRXN), parse,
     1         savecf, saverf, rxnok, match
      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF), rec1,
     1         nm(MAXRXN), na(MAXRXN), ng(MAXRXN), nw(MAXRXN),
     2         rec1m(MAXRXN,MAXMIN), rec1a(MAXRXN,MAXAQS), 
     3         rec1g(MAXRXN,MAXGAS), univar, useLVS, epseqn, geqn,
     4         rec1m1, rec1m2, rec1m3, rec1m4, rec1gg, rec1aa
      DOUBLE PRECISION  coefm(MAXRXN,MAXMIN), coefa(MAXRXN,MAXAQS),
     1                  coefg(MAXRXN,MAXGAS), coefw(MAXRXN)

      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /rlimit/ nmin1,  nmin2,  nmin3,  nmin4,  ngas,   naqs,
     1                rec1m1, rec1m2, rec1m3, rec1m4, rec1gg, rec1aa
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /reac1/  rtitle
      COMMON /reac2/  coefm, coefa, coefg, coefw, nm, na, ng, nw,
     1                rec1m, rec1a, rec1g, m2reac
      COMMON /fnames/ namecf, namerf, nametf, namepf
      COMMON /saveif/ savecf, saverf

      SAVE


      nm1234 = nmin1 + nmin2 + nmin3 + nmin4 
      nmga   = nm1234 + ngas + naqs

***** prompt for / read nreac *****

   1  WRITE(wterm,5) 
   5  FORMAT(/,' specify number of reactions to be processed: ',/)
      READ(rterm,*) nreac
      IF (nreac .LE. 0) GO TO 1

      DO 10 ireac = 1,nreac

***** prompt for / read specifications for next reaction *****

           WRITE(wterm,15) ireac, nreac
  15       FORMAT(/,' input title for reaction ',i2,' of ',i2,':',/)
           READ(rterm,25) rtitle(ireac)
  25       FORMAT(a80)

 333       WRITE(wterm,35) ireac
  35       FORMAT(/,' enter [coeff  species] pairs, separated by'
     1            /,' blanks, one pair per line, for reaction ',i2,
     2            /,' (conclude with [0 done]) ',/)

           ibad = 0
           m2reac(ireac) = .FALSE.
           nm(ireac) = 0
           ng(ireac) = 0
           na(ireac) = 0
           nw(ireac) = 0

 111       READ(rterm,112) string
 112       FORMAT(a80)

           IF (.NOT. parse(string,coeff,specie)) THEN
                WRITE(wterm,113)
 113            FORMAT(/,' ill-defined [coeff species] pair; ',
     1                   'try again',/) 
                GO TO 111
           END IF

           IF (coeff .EQ. 0.0d0) THEN
******          reaction stoichiometry complete ******
                IF (ibad .NE. 0) THEN
                     CALL wrtbad(ibad,sbad)       
                     GO TO 111 
                ELSE
******               ensure that stoichiometry is correct
                     CALL chkrxn(ireac,namem,namea,nameg,
     1                           formm,forma,formg,rxnok)
                     IF (.NOT. rxnok) THEN
                          GO TO 333
                     END IF
		END IF
           ELSE
******          determine disposition of current specie: either H2O, 
******          found or not found within the current database
		IF (specie .EQ. 'H2O') THEN
		     nw(ireac) = 1
                     coefw(ireac) = coeff
                ELSE
                     IF (match(specie,form,rec1,rec1m1,1,nmga,nm1234))
     1               THEN
******                    update [n|coef|rec1][m|g|a]; continue
                          CALL umaker(ireac,coeff,specie,form,rec1,
     1                    namem,namea,nameg,formm,forma,formg)
                     ELSE
                          ibad = ibad + 1
		          sbad(ibad) = specie
                     END IF
                END IF
                GO TO 111
           END IF 
  10       CONTINUE

****** set wetrxn variable ******

      iwet = 0
      wetrxn = .FALSE.
      IF ((isat .EQ. 1) .OR. (iopt .EQ. 1)) THEN
	   wetrxn = .TRUE.
	   iwet = 1
      ELSE
           DO 70 ireac = 1,nreac
                IF ((nw(ireac) .EQ. 1) .OR. (na(ireac) .GT. 0)) THEN	
	             wetrxn = .TRUE.
		     iwet = 1
		     GO TO 444
                END IF
  70            CONTINUE
      END IF

****** save reaction file if desired ******

 444  WRITE(wterm,125) 
 125  FORMAT(/,' would you like to save these reactions to a file ',
     1         '(y/n)',/)
      READ(rterm,135) ans
 135  FORMAT(a1)
      IF ((ans .NE. 'y') .AND. (ans .NE. 'Y') .AND.
     1    (ans .NE. 'n') .AND. (ans .NE. 'N')) GO TO 444

      saverf = ((ans .EQ. 'y') .OR.  (ans .EQ. 'Y'))

      IF (saverf) THEN
 555       WRITE(wterm,145)
 145       FORMAT(/,' specify file name:',/)
           READ(rterm,155) namerf
 155       FORMAT(a20)
           IF (.NOT. openf(wterm,reacf,namerf,2,1,1,132)) THEN
                GO TO 555
           ELSE
***             write generic header
                WRITE(reacf,205) 
 205            FORMAT(' Line 1:  nreac, iwet',12x,
     1                 '(free format)')
                WRITE(reacf,210) 
 210            FORMAT(' Line 2:  [blank]',16x,
     1                 '(free format)')
                WRITE(reacf,215) 
 215            FORMAT(' Line 3:  descriptive title',6x,
     1                 '(a80)')
                WRITE(reacf,220) 
 220            FORMAT(' Line 4:  nm, na, ng, nw',9x,
     1                 '(free format)')
                WRITE(reacf,225) 
 225            FORMAT(' nm Lines:  coeff  mname  mform',2x,
     1                 '(1x,f9.3,2x,a20,2x,a30)')
                WRITE(reacf,230) 
 230            FORMAT(' ng Lines:  coeff  aname  aform',2x,
     1                 '(1x,f9.3,2x,a20,2x,a30)')
                WRITE(reacf,235) 
 235            FORMAT(' na Lines:  coeff  gname  gform',2x,
     1                 '(1x,f9.3,2x,a20,2x,a30)')
                WRITE(reacf,240) 
 240            FORMAT(' [1 Line:   coeff  H2O    H2O] ',2x,
     1                 '(1x,f9.3,2x,a20,2x,a30)',/)
                WRITE(reacf,245) 
 245            FORMAT('*** each of the nreac reaction blocks',/,
     1                 '*** contains 3+nm+ng+na+nw lines',/)
                WRITE(reacf,250)
 250            FORMAT(55('*'),/)

***             write reaction information

                WRITE(reacf,165) nreac, iwet
 165            FORMAT(2(1x,i2))
           END IF

	   DO 80 ireac = 1,nreac
                WRITE(reacf,175) rtitle(ireac), nm(ireac),
     1                           na(ireac), ng(ireac), nw(ireac)
 175            FORMAT(/,1x,a80,/,4(1x,i3))

                IF (nm(ireac) .GT. 0) WRITE(reacf,185) 
     1          (coefm(ireac,imin), namem(ireac,imin), 
     2          formm(ireac,imin), imin = 1,nm(ireac))
 185            FORMAT(1x,f9.3,2x,a20,2x,a30)

                IF (na(ireac) .GT. 0) WRITE(reacf,185) 
     1          (coefa(ireac,iaqs), namea(ireac,iaqs), 
     2          forma(ireac,iaqs), iaqs = 1,na(ireac))

                IF (ng(ireac) .GT. 0) WRITE(reacf,185) 
     1          (coefg(ireac,igas), nameg(ireac,igas), 
     2          formg(ireac,igas), igas = 1,ng(ireac))

                IF (nw(ireac) .EQ. 1) WRITE(reacf,195) 
     1          coefw(ireac), 'H2O                 ', 
     2                        'H2O                           ' 
 195            FORMAT(1x,f9.3,2x,a20,2x,a30)

  80            CONTINUE

      END IF

      END

***************************************************************

*** nxtrec - Get rec1 for next database species.

      INTEGER FUNCTION nxtrec(irec,mga,nm1234)

      INTEGER rec1m1, rec1m2, rec1m3, rec1m4, rec1gg, rec1aa

      COMMON /rlimit/ nmin1,  nmin2,  nmin3,  nmin4,  ngas,   naqs,
     1                rec1m1, rec1m2, rec1m3, rec1m4, rec1gg, rec1aa

      SAVE


      IF (mga .LE. nmin1) THEN
***        one-phase mineral
           nxtrec = irec + 7
           RETURN
      END IF

      IF (mga .GT. nm1234) THEN
***        gas or aqueous species
           nxtrec = irec + 6
           RETURN
      END IF

      IF (mga .LE. (nmin1 + nmin2)) THEN
***        two-phase mineral
           nxtrec = irec + 8
           RETURN
      END IF

      IF (mga .LE. (nmin1 + nmin2 + nmin3)) THEN
***        three-phase mineral
           nxtrec = irec + 9
      ELSE
***        four-phase mineral
           nxtrec = irec + 10
      END IF

      RETURN

      END

*******************************************************************

*** readrf - Open/read user-specified reaction (RXN) file.

      SUBROUTINE readrf(nreac,wetrxn)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMIN = 10, MAXAQS = 10, MAXGAS = 10, MAXRXN = 3000,
     1           NPLOTF = 8)

      CHARACTER*20  namecf, namerf, nametf, spname, namepf(NPLOTF),
     1              pfname
      CHARACTER*30  form
      CHARACTER*80  rtitle(MAXRXN)
      LOGICAL  openf, wetrxn, m2reac(MAXRXN), savecf, saverf, match
      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     1         rec1m1, rec1m2, rec1m3, rec1m4, rec1gg, rec1aa
      INTEGER  nm(MAXRXN), na(MAXRXN), ng(MAXRXN), nw(MAXRXN),
     1         rec1m(MAXRXN,MAXMIN), rec1a(MAXRXN,MAXAQS), 
     2         rec1g(MAXRXN,MAXGAS), rec1
      DOUBLE PRECISION  coefm(MAXRXN,MAXMIN), coefa(MAXRXN,MAXAQS),
     1                  coefg(MAXRXN,MAXGAS), coefw(MAXRXN)

      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /reac1/  rtitle
      COMMON /reac2/  coefm, coefa, coefg, coefw, nm, na, ng, nw,
     1                rec1m, rec1a, rec1g, m2reac
      COMMON /fnames/ namecf, namerf, nametf, namepf
      COMMON /saveif/ savecf, saverf
      COMMON /rlimit/ nmin1,  nmin2,  nmin3,  nmin4,  ngas,   naqs,
     1                rec1m1, rec1m2, rec1m3, rec1m4, rec1gg, rec1aa
      COMMON /dapron/ pfname

      SAVE


      nm1234 = nmin1 + nmin2 + nmin3 + nmin4 
      nmga   = nm1234 + ngas + naqs

  1   WRITE(wterm,10)
 10   FORMAT(/,' specify name of reaction file:',/)
      READ(rterm,20) namerf
 20   FORMAT(a20)
      IF (.NOT. openf(wterm,reacf,namerf,1,1,1,132)) GO TO 1

      saverf = .TRUE.

***** read number of reactions and their wet/dry character ******

***** skip first 13 comment lines
      READ(reacf,25)
 25   FORMAT(////////////)

      READ(reacf,*) nreac, iwet
      wetrxn = (iwet .EQ. 1)

      DO 30 ireac = 1,nreac

***** read title, nm, na, ng, nw for next reaction *****

           READ(reacf,40) rtitle(ireac)
 40        FORMAT(/,1x,a80)
           READ(reacf,*) nm(ireac), na(ireac), ng(ireac), nw(ireac)
        
***** read mineral, aqueous species, gas, H2O stoichiometry *****

           m2reac(ireac) = .FALSE.
           IF (nm(ireac) .GT. 0) THEN
                DO 50  imin = 1,nm(ireac)
                     READ(reacf,51) coefm(ireac,imin), spname, form 
 51                  FORMAT(1x,f9.3,2x,a20,2x,a30)
                     IF (.NOT. match(spname,form,rec1,rec1m1,1,
     1                               nm1234,nm1234)) THEN
                          GO TO 999
                     ELSE
                          rec1m(ireac,imin) = rec1
                          IF (rec1m(ireac,imin) .GE. rec1m2) THEN
                               m2reac(ireac) = .TRUE.
                          END IF
                     END IF
 50                  CONTINUE
           END IF

           IF (na(ireac) .GT. 0) THEN
                istart = nm1234 + ngas + 1
                DO 60  iaqs = 1,na(ireac)
                     READ(reacf,51) coefa(ireac,iaqs), spname, form
                     IF (.NOT. match(spname,form,rec1,rec1aa,istart,
     1                               nmga,nm1234)) THEN
                          GO TO 999
                     ELSE
                          rec1a(ireac,iaqs) = rec1
                     END IF
 60                  CONTINUE
           END IF

           IF (ng(ireac) .GT. 0) THEN
                istart = nm1234 + 1
                iend = nm1234 + ngas
                DO 70  igas = 1,ng(ireac)
                     READ(reacf,51) coefg(ireac,igas), spname, form
                     IF (.NOT. match(spname,form,rec1,rec1gg,istart,
     1                               iend,nm1234)) THEN
                          GO TO 999
                     ELSE
                          rec1g(ireac,igas) = rec1
                     END IF
 70                  CONTINUE
           END IF

           IF (nw(ireac) .EQ. 0) THEN
                coefw(ireac) = 0
           ELSE
                READ(reacf,*) coefw(ireac)
           END IF
  
 30        CONTINUE

      RETURN

 999  WRITE(wterm,1000) ireac, spname, pfname 
 1000 FORMAT(//,' Reaction ',i2,' species ',a20,
     1        /,' not found in database ',a20,/
     1        /,' re-run with correct database or re-create',
     2        /,' reaction file from this database.',/)
      STOP

      END

************************************************************************

*** wrtbad - Write the list of species not found in database pfname;
***          prompt for repeats.

      SUBROUTINE wrtbad(ibad,sbad)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NPLOTF = 8, MAXBAD = 10)

      INTEGER rterm, wterm, iconf, reacf, pronf, tabf, 
     1        plotf(NPLOTF)
      CHARACTER*20 sbad(MAXBAD), pfname

      COMMON /io/    rterm, wterm, iconf, reacf, pronf, tabf, plotf 
      COMMON /dapron/ pfname

      SAVE


      WRITE(wterm,45) pfname 
  45  FORMAT(/,' the following species were not', 
     1       /,' found in database ',a20,/)

      DO 20 i = 1,ibad
           WRITE(wterm,55) sbad(i)
  55       FORMAT(5x,a20)
  20       CONTINUE

      WRITE(wterm,65) 
  65  FORMAT(/,' input new [coeff  species] pairs',
     1       /,' to replace these incorrect entries',
     2       /,' (conclude with [0 done]) ',/)

      ibad = 0

      END

************************************************************************

*** chkrxn - Give the user a chance to look over rxn stoichiometry;
***          if it's ok, then rxnok returns .TRUE.; otherwise,
***          returns .FALSE.

      SUBROUTINE chkrxn(ireac,namem,namea,nameg,formm,forma,formg,rxnok)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMIN = 10, MAXAQS =  10, MAXGAS =  10, MAXRXN = 3000,
     1           NPLOTF = 8)

      CHARACTER*1  ans
      CHARACTER*20 namem(MAXRXN,MAXMIN), namea(MAXRXN,MAXAQS),
     1             nameg(MAXRXN,MAXGAS), namew
      CHARACTER*30 formm(MAXRXN,MAXMIN), formg(MAXRXN,MAXGAS),
     1             forma(MAXRXN,MAXAQS), formw
      CHARACTER*80 rtitle(MAXRXN)

      LOGICAL m2reac(MAXRXN), rxnok

      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     1         nm(MAXRXN), na(MAXRXN), ng(MAXRXN), nw(MAXRXN),
     2         rec1m(MAXRXN,MAXMIN), rec1a(MAXRXN,MAXAQS), 
     3         rec1g(MAXRXN,MAXGAS)

      DOUBLE PRECISION  coefm(MAXRXN,MAXMIN), coefa(MAXRXN,MAXAQS),
     1                  coefg(MAXRXN,MAXGAS), coefw(MAXRXN)

      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /reac1/  rtitle
      COMMON /reac2/  coefm, coefa, coefg, coefw, nm, na, ng, nw,
     1                rec1m, rec1a, rec1g, m2reac
      
      SAVE

      DATA namew, formw /
     1    'H2O                 ',
     2    'H2O                           '/


      WRITE(wterm,75) ireac
  75  FORMAT(/,' reaction ',i2,' stoichiometry:',/) 

***** write reactants

      DO 30 imin = 1,nm(ireac)
           IF (coefm(ireac,imin) .LT. 0.0d0) THEN
                WRITE(wterm,85) coefm(ireac,imin),
     1          namem(ireac,imin), formm(ireac,imin)
  85            FORMAT(6x,f7.3,3x,a20,3x,a30)
           END IF
  30       CONTINUE

      DO 40 igas = 1,ng(ireac)
           IF (coefg(ireac,igas) .LT. 0.0d0) THEN
                WRITE(wterm,86) coefg(ireac,igas),
     1          formg(ireac,igas)(1:20), nameg(ireac,igas)
  86            FORMAT(6x,f7.3,3x,a20,3x,a20)
           END IF
  40       CONTINUE

      DO 50 iaqs = 1,na(ireac)
           IF (coefa(ireac,iaqs) .LT. 0.0d0) THEN
                WRITE(wterm,85) coefa(ireac,iaqs),
     1          namea(ireac,iaqs), forma(ireac,iaqs)
           END IF
  50       CONTINUE

      IF ((nw(ireac) .EQ. 1) .AND. (coefw(ireac) .LT. 0.0d0)) THEN 
           WRITE(wterm,85) coefw(ireac), namew, formw 
      END IF

***** write products

      DO 31 imin = 1,nm(ireac)
           IF (coefm(ireac,imin) .GT. 0.0d0) THEN
                WRITE(wterm,85) coefm(ireac,imin),
     1          namem(ireac,imin), formm(ireac,imin)
           END IF
  31       CONTINUE

      DO 41 igas = 1,ng(ireac)
           IF (coefg(ireac,igas) .GT. 0.0d0) THEN
                WRITE(wterm,86) coefg(ireac,igas),
     1          formg(ireac,igas)(1:20), nameg(ireac,igas)
           END IF
  41       CONTINUE

      DO 51 iaqs = 1,na(ireac)
           IF (coefa(ireac,iaqs) .GT. 0.0d0) THEN
                WRITE(wterm,85) coefa(ireac,iaqs),
     1          namea(ireac,iaqs), forma(ireac,iaqs)
           END IF
  51       CONTINUE

      IF ((nw(ireac) .EQ. 1) .AND. (coefw(ireac) .GT. 0.0d0)) THEN 
           WRITE(wterm,85) coefw(ireac), namew, formw
      END IF
      
 222  WRITE(wterm,95)
  95  FORMAT(/,' is this correct? (y/n)',/) 
      READ(rterm,105) ans
 105  FORMAT(a1)
	
      IF ((ans .NE. 'Y') .AND. (ans .NE. 'y') .AND.
     1    (ans .NE. 'N') .AND. (ans .NE. 'n')) THEN
          GO TO 222
      ELSE 
          rxnok = ((ans .EQ. 'Y') .OR. (ans .EQ. 'y')) 
      END IF

      END

*******************************************************************

*** match - Returns .TRUE. (and rec1sp) if specie is found in
***         database pfname; otherwise returns .FALSE.    


      LOGICAL FUNCTION match(specie,form,rec1sp,rec1ty,first,last,
     1                       nm1234)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NPLOTF = 8)

      CHARACTER*20 specie, name
      CHARACTER*30 form
      INTEGER      rterm, wterm, iconf, reacf, pronf, tabf, 
     1             plotf(NPLOTF), rec1sp, rec1ty, first, last

      COMMON /io/  rterm, wterm, iconf, reacf, pronf, tabf, plotf

      SAVE


      irec = rec1ty 
      DO 60 mga = first,last
           READ(pronf,115,REC=irec) name, form
 115       FORMAT(1x,a20,a30)
           IF (specie .EQ. name) THEN 
                match = .TRUE.
                rec1sp = irec
                RETURN
           ELSE
                irec = nxtrec(irec,mga,nm1234)
           END IF
  60       CONTINUE

       match = .FALSE.

       RETURN
       END 

*******************************************************************

*** umaker - Update /reac/ arrays to include current species.

      SUBROUTINE umaker(ireac,coeff,specie,form,rec1,
     1                  namem,namea,nameg,formm,forma,formg)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMIN = 10, MAXAQS =  10, MAXGAS =  10, MAXRXN = 3000)

      CHARACTER*20 specie, namem(MAXRXN,MAXMIN), namea(MAXRXN,MAXAQS),
     1             nameg(MAXRXN,MAXGAS)
      CHARACTER*30 form, formm(MAXRXN,MAXMIN), formg(MAXRXN,MAXGAS),
     1             forma(MAXRXN,MAXAQS)
      CHARACTER*80 rtitle(MAXRXN)

      LOGICAL  m2reac(MAXRXN)

      INTEGER  nm(MAXRXN), na(MAXRXN), ng(MAXRXN), nw(MAXRXN),
     2         rec1m(MAXRXN,MAXMIN), rec1a(MAXRXN,MAXAQS), 
     3         rec1g(MAXRXN,MAXGAS), 
     4         rec1m1, rec1m2, rec1m3, rec1m4, rec1gg, rec1aa, rec1

      DOUBLE PRECISION coefm(MAXRXN,MAXMIN), coefa(MAXRXN,MAXAQS),
     1                 coefg(MAXRXN,MAXGAS), coefw(MAXRXN)

      COMMON /rlimit/ nmin1,  nmin2,  nmin3,  nmin4,  ngas,   naqs,
     1                rec1m1, rec1m2, rec1m3, rec1m4, rec1gg, rec1aa
      COMMON /reac1/  rtitle
      COMMON /reac2/  coefm, coefa, coefg, coefw, nm, na, ng, nw,
     1                rec1m, rec1a, rec1g, m2reac

      SAVE


      IF (rec1 .GE. rec1aa) THEN
           na(ireac) = na(ireac) + 1
           coefa(ireac,na(ireac)) = coeff
           rec1a(ireac,na(ireac)) = rec1
           namea(ireac,na(ireac)) = specie
           forma(ireac,na(ireac)) = form
           RETURN
      END IF

      IF (rec1 .GE. rec1gg) THEN
           ng(ireac) = ng(ireac) + 1
           coefg(ireac,ng(ireac)) = coeff
           rec1g(ireac,ng(ireac)) = rec1
           nameg(ireac,ng(ireac)) = specie
           formg(ireac,ng(ireac)) = form
      ELSE
           nm(ireac) = nm(ireac) + 1
           coefm(ireac,nm(ireac)) = coeff
           rec1m(ireac,nm(ireac)) = rec1
           namem(ireac,nm(ireac)) = specie
           formm(ireac,nm(ireac)) = form
           IF (rec1 .GE. rec1m2) THEN
                m2reac(ireac) = .TRUE.
           END IF
      END IF
       
      RETURN
      END

************************************************************************

*** openf -  Returns .TRUE. and opens the file specified by fname, 
***          fstat, facces, fform, and frecl if this file exists and is
***          accessible; otherwise, returns .FALSE. and prints an
***          appropriate error message to the device specified by iterm.
***

      LOGICAL FUNCTION openf(iterm,iunit,fname,istat,iacces,iform,irecl)      

      CHARACTER*11  fform(2)
      CHARACTER*10  facces(2)
      CHARACTER*20  fname
      CHARACTER*3   fstat(2)

      SAVE

      DATA fform  / 'FORMATTED  ',  'UNFORMATTED' /
      DATA facces / 'SEQUENTIAL',   'DIRECT    '  /
      DATA fstat  / 'OLD',          'NEW'         /


      openf = .FALSE.
      
      IF ((iacces .LT. 1) .OR. (iacces .GT. 2) .OR.
     1    (iform  .LT. 1) .OR. (iform  .GT. 2) .OR.
     2    (istat  .LT. 1) .OR. (istat  .GT. 2)) GO TO 10

      IF (iacces .EQ. 1) THEN
           OPEN(UNIT=iunit,FILE=fname,ACCESS=facces(iacces),
     1          FORM=fform(iform),STATUS=fstat(istat),ERR=10)
           openf = .TRUE.
           RETURN
      ELSE
           OPEN(UNIT=iunit,FILE=fname,ACCESS=facces(iacces),
     1          FORM=fform(iform),STATUS=fstat(istat),RECL=irecl,
     2          ERR=10)
           openf = .TRUE.
           RETURN
      END IF

 10   WRITE(iterm,20)
 20   FORMAT(/,' nonexistant file or invalid specifications',
     1         ' ... try again',/)
      RETURN

      END
