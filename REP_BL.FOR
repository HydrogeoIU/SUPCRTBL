*** REPBL - Collection of routines that write the calculated standard 
***         molal thermodynamic properties of reactions to the TAB
***         file and, optionally, the PLOT files.
***
*********************************************************************
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
*********************************************************************

*** tabtop - Write global header for output file.

      SUBROUTINE tabtop

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXODD = 21, NPLOTF = 8)

      LOGICAL      savecf, saverf
      INTEGER      rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     1             mapiso(2,3), mapinc(2,3), mapv3(2,3),
     2             univar, useLVS, epseqn, geqn, xyplot, end 
      DOUBLE PRECISION  isomin, isomax, isoinc, Kmin, Kmax, Kinc,
     1                  oddv1(MAXODD), oddv2(MAXODD)
      CHARACTER*4  incvar(2,3)
      CHARACTER*10 isov(2,3), incv(2,3), var3(2,3), isosat(2)
      CHARACTER*12 isovar(2,3)
      CHARACTER*20 pfname, namecf, namerf, nametf, namepf(NPLOTF),
     1             nosave

      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /dapron/ pfname
      COMMON /stvars/ isosat, isovar, incvar
      COMMON /TPDmap/ mapiso, mapinc, mapv3
      COMMON /headmp/ isov, incv, var3
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /grid/   isomin, isomax, isoinc, v2min, v2max, v2inc,
     2                oddv1, oddv2, Kmin, Kmax, Kinc, niso, nv2, nlogK
      COMMON /fnames/ namecf, namerf, nametf, namepf
      COMMON /saveif/ savecf, saverf
      COMMON /plottr/ xyplot, end, nplots

      SAVE

      DATA nosave / 'file not saved      ' /


      WRITE(tabf,5) 
  5   FORMAT(/,'********** SUPCRTBL: input/output specifications for',
     1         ' this run',/)

      IF (savecf) THEN
           WRITE(tabf,15) namecf
      ELSE
           WRITE(tabf,15) nosave
      END IF
 15   FORMAT(  '            USER-SPECIFIED  CON FILE  containing ',
     1       /,'            T-P-D grid & option switches: ',a20,/)

      IF (saverf) THEN
           WRITE(tabf,25) namerf
      ELSE
           WRITE(tabf,25) nosave
      END IF
 25   FORMAT(  '            USER-SPECIFIED  RXN FILE  containing ',
     1       /,'            chemical reactions: ',a20,/)

      WRITE(tabf,35) pfname
 35   FORMAT(  '            THERMODYNAMIC DATABASE: ',a20,/)

      WRITE(tabf,45) nametf
 45   FORMAT(  '            SUPCRT-GENERATED  TAB FILE  containing ',
     1       /,'            tabulated reaction properties ',
     2                      '(this file): ',a20)

      IF (xyplot .GT. 0) THEN 
        WRITE(tabf,47) namepf(1)
 47     FORMAT(/,'            SUPCRT-GENERATED  PLT FILES  containing ',
     1         /,'            reaction properties for x-y plots: '
     2          ,a20,' etc.')
      END IF

      CALL wrtopt

      END

*******************************************************************

*** wrtopt - Write various switch options to tabular output file.

      SUBROUTINE wrtopt

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXODD = 21, NPLOTF = 8)

      INTEGER      rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     2             univar, useLVS, epseqn, geqn 
      DOUBLE PRECISION  isomin, isomax, isoinc, Kmin, Kmax, Kinc,
     1                  oddv1(MAXODD), oddv2(MAXODD)
      CHARACTER*4  incvar(2,3)
      CHARACTER*10 isov(2,3), incv(2,3), var3(2,3), isosat(2)
      CHARACTER*12 isovar(2,3)
      CHARACTER*20 namecf, namerf, nametf, namepf(NPLOTF)

      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /stvars/ isosat, isovar, incvar
      COMMON /headmp/ isov, incv, var3
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /grid/   isomin, isomax, isoinc, v2min, v2max, v2inc, 
     2                oddv1, oddv2, Kmin, Kmax, Kinc, niso, nv2, nlogK
      COMMON /fnames/ namecf, namerf, nametf, namepf

      SAVE


      WRITE(tabf,75) 
 75   FORMAT(/,'********** summary of option switches ',/)

      WRITE(tabf,85) isat,iopt,iplot,univar,noninc
 85   FORMAT(  '            isat, iopt, iplot, univar, noninc: ',5i3)

*** useLVS, epseqn, geqn not written to TAB for distribution copies
***
*      WRITE(tabf,105) useLVS,epseqn,geqn
* 105  FORMAT(  '            useLVS, epseqn, geqn:              ',3i3)

      WRITE(tabf,115) 
 115  FORMAT(/,'********** summary of state conditions ',/)

      IF (noninc .EQ. 0) THEN
           IF (isat .EQ. 0) THEN
                WRITE(tabf,125) isovar(iopt,iplot),isomin,isomax,isoinc
 125            FORMAT(12x,'ISO',a12,':  min, max, increment:',
     1                 3(2x,f10.4))
                WRITE(tabf,135) incv(iopt,iplot),v2min, v2max, v2inc
 135            FORMAT(12x,a10,' range: min, max, increment:',
     1                 3(2x,f10.4))
           ELSE
                WRITE(tabf,145) isosat(iopt),v2min, v2max, v2inc
 145            FORMAT(12x,'saturation ',a10,' range: min, max,',
     1                     ' increment:',3(2x,f10.4))
           END IF
      ELSE
           IF (isat .EQ. 0) THEN
                WRITE(tabf,155) isov(iopt,iplot), incv(iopt,iplot), 
     1                          noninc
 155            FORMAT(12x,'nonincremental ',a10,', ',a10,
     1                 ' coordinates: ',i2,' pair')
           ELSE
                WRITE(tabf,165) isosat(iopt), noninc
 165            FORMAT(12x,'nonincremental saturation ',a10,': ',i2,
     1                 ' points')
           END IF
      END IF

      END

*********************************************************************

*** wrtop2 - Write various switch options to plot file k.

      SUBROUTINE wrtop2(k)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXODD = 21, NPLOTF = 8)

      INTEGER      rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     2             univar, useLVS, epseqn, geqn 
      DOUBLE PRECISION  isomin, isomax, isoinc, Kmin, Kmax, Kinc,
     1                  oddv1(MAXODD), oddv2(MAXODD)
      CHARACTER*4  incvar(2,3)
      CHARACTER*10 isov(2,3), incv(2,3), var3(2,3), isosat(2)
      CHARACTER*12 isovar(2,3)
      CHARACTER*20 namecf, namerf, nametf, namepf(NPLOTF)

      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /stvars/ isosat, isovar, incvar
      COMMON /headmp/ isov, incv, var3
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /grid/   isomin, isomax, isoinc, v2min, v2max, v2inc, 
     2                oddv1, oddv2, Kmin, Kmax, Kinc, niso, nv2, nlogK
      COMMON /fnames/ namecf, namerf, nametf, namepf

      SAVE


      WRITE(plotf(k),75) 
 75   FORMAT(/,'*********** summary of option switches ',/)

      WRITE(plotf(k),85) isat,iopt,iplot,univar,noninc
 85   FORMAT(  '            isat, iopt, iplot, univar, noninc: ',5i3)

*** useLVS, epseqn, geqn not written to PLOT for distribution copies
***
*      WRITE(plotf(k),105) useLVS,epseqn,geqn
* 105  FORMAT(  '            useLVS, epseqn, geqn:              ',3i3)

      WRITE(plotf(k),115) 
 115  FORMAT(/,'*********** summary of state conditions ',/)

      IF (noninc .EQ. 0) THEN
           IF (isat .EQ. 0) THEN
                WRITE(plotf(k),125) isovar(iopt,iplot),
     1                              isomin,isomax,isoinc
 125            FORMAT(12x,'ISO',a12,':  min, max, increment:',
     1                 3(2x,f10.4))
                WRITE(plotf(k),135) incv(iopt,iplot),v2min, v2max, v2inc
 135            FORMAT(12x,a10,' range: min, max, increment:',
     1                 3(2x,f10.4))
           ELSE
                WRITE(plotf(k),145) isosat(iopt),v2min, v2max, v2inc
 145            FORMAT(12x,'saturation ',a10,' range: min, max,',
     1                     ' increment:',3(2x,f10.4))
           END IF
      ELSE
           IF (isat .EQ. 0) THEN
                WRITE(plotf(k),155) isov(iopt,iplot), incv(iopt,iplot), 
     1                          noninc
 155            FORMAT(12x,'nonincremental ',a10,', ',a10,
     1                 ' coordinates: ',i2,' pair')
           ELSE
                WRITE(plotf(k),165) isosat(iopt), noninc
 165            FORMAT(12x,'nonincremental saturation ',a10,': ',i2,
     1                 ' points')
           END IF
      END IF
           
      WRITE(plotf(k),175)
 175  FORMAT(/,86('*'))

      END

***************************************************************

*** wrtrxn - Write header information for the i[th] reaction to 
***          tabulated output file tabf and (if appropriate) to 
***          the plot files.

      SUBROUTINE wrtrxn(i)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

**JP**      PARAMETER (MAXMIN = 10, MAXGAS = 10, MAXAQS = 10,
**JP**     1           MXTRAN =  3, IABC   =  3, MAXRXN = 3000, NPLOTF = 8)
      PARAMETER (MAXMIN = 10, MAXGAS = 10, MAXAQS = 10,
     1           MXTRAN =  3, IABC   =  4, MAXRXN = 3000, NPLOTF = 8)
      CHARACTER*4   incvar(2,3)
      CHARACTER*10  isov(2,3), incv(2,3), var3(2,3), isosat(2)
      CHARACTER*12  isovar(2,3)
      CHARACTER*20  mname(MAXMIN), gname(MAXGAS), aname(MAXAQS)
      CHARACTER*30  mform(MAXMIN), gform(MAXGAS), aform(MAXAQS)
      CHARACTER*80  rtitle(MAXRXN)

      LOGICAL  m2reac(MAXRXN), nullrx

      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     1         univar, useLVS, epseqn, geqn, xyplot, end

      INTEGER  nm(MAXRXN), na(MAXRXN), ng(MAXRXN), nw(MAXRXN),
     1         rec1m(MAXRXN,MAXMIN), rec1a(MAXRXN,MAXAQS), 
     2         rec1g(MAXRXN,MAXGAS)

      INTEGER ntran(MAXMIN)

      DOUBLE PRECISION coefm(MAXRXN,MAXMIN), coefa(MAXRXN,MAXAQS),
     1                 coefg(MAXRXN,MAXGAS), coefw(MAXRXN)


      DOUBLE PRECISION  Gfmin(MAXMIN), Hfmin(MAXMIN), VPrTrm(MAXMIN),
     1                  SPrTrm(MAXMIN), Alpha(MAXMIN), Kappa(MAXMIN),
     2                  KappaP(MAXMIN), KappaDP(MAXMIN),
     3                  MK1(IABC,MAXMIN), MK2(IABC,MAXMIN),
     4                  MK3(IABC,MAXMIN), MK4(IABC,MAXMIN),
     5                  Ttran(MXTRAN,MAXMIN), Htran(MXTRAN,MAXMIN),
     6                  Vtran(MXTRAN,MAXMIN), dPdTtr(MXTRAN,MAXMIN),
     7                  Tmaxm(MAXMIN), NumAtoms(MAXMIN), Tc(MAXMIN),
     8                  Smax(MAXMIN), Vmax(MAXMIN)

      DOUBLE PRECISION Gfgas(MAXGAS), Hfgas(MAXGAS), VPrTrg(MAXGAS), 
     1                 SPrTrg(MAXGAS), MKg(IABC,MAXGAS), Tmaxg(MAXGAS)

      DOUBLE PRECISION Gfaqs(MAXAQS), Hfaqs(MAXAQS), SPrTra(MAXAQS), 
     1                 a(4,MAXAQS), c(2,MAXAQS), 
     2                 wref(MAXAQS), chg(MAXAQS)

      COMMON /stvars/ isosat, isovar, incvar
      COMMON /headmp/ isov, incv, var3
      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn

      COMMON /reac1/  rtitle
      COMMON /reac2/  coefm, coefa, coefg, coefw, nm, na, ng, nw,
     1                rec1m, rec1a, rec1g, m2reac

      COMMON /mnames/ mname, mform 

      COMMON /minref/ Gfmin, Hfmin, SPrTrm, VPrTrm, MK1, MK2, MK3, MK4, 
     1                Ttran, Htran, Vtran, dPdTtr, Tmaxm, ntran, Alpha,
     2                Kappa, KappaP, KappaDP, NumAtoms, Tc, Smax, Vmax

      COMMON /gnames/ gname, gform
      COMMON /gasref/ Gfgas, Hfgas, SPrTrg, VPrTrg, MKg, Tmaxg
 
      COMMON /anames/ aname, aform
      COMMON /aqsref/ Gfaqs, Hfaqs, SPrTra, c, a, wref, chg

      COMMON /plottr/ xyplot, end, nplots

      SAVE

***** write header to TAB file *****

      WRITE(tabf,5) i
  5   FORMAT(//,36('*'),' REACTION ',i3,1x,36('*'),/)

***** write reaction title and stoichiometry to TAB file *****

      WRITE(tabf,15) rtitle(i)
 15   FORMAT(//,a80,//)

      WRITE(tabf,25) 
 25   FORMAT(' REACTION STOICHIOMETRY: ',/)

      WRITE(tabf,26)
 26   FORMAT(8x,' COEFF.',3x,'NAME',16x,3x,'FORMULA',22x,/,
     3       8x,7('-'),3x,20('-'),3x,20('-'))

***** write reactants *****

      IF (nm(i) .GT. 0) THEN
           DO 10 j = 1,nm(i)
                IF (coefm(i,j) .LT. 0.0d0) THEN
                     WRITE(tabf,35) coefm(i,j), mname(j), mform(j)
 35                  FORMAT(6x,f9.3,3x,a20,3x,a30)
                END IF
 10             CONTINUE
      END IF

      IF (ng(i) .GT. 0) THEN
           DO 20 j = 1,ng(i)
                IF (coefg(i,j) .LT. 0.0d0) THEN
                     WRITE(tabf,36) coefg(i,j), gform(j)(1:20), gname(j)
 36                  FORMAT(6x,f9.3,3x,a20,3x,a20)
                END IF
 20             CONTINUE
      END IF

      IF (na(i) .GT. 0) THEN
           DO 30 j = 1,na(i)
                IF (coefa(i,j) .LT. 0.0d0) THEN
                     WRITE(tabf,35) coefa(i,j), aname(j), aform(j)
                END IF
 30             CONTINUE
      END IF

      IF ((nw(i) .GT. 0) .AND. (coefw(i) .LT. 0.0d0)) THEN
           WRITE(tabf,55) coefw(i)
 55        FORMAT(6x,f9.3,3x,'H2O',17x,3x,'H2O')
      END IF

***** write products *****

      IF (nm(i) .GT. 0) THEN
           DO 11 j = 1,nm(i)
                IF (coefm(i,j) .GT. 0.0d0) THEN
                     WRITE(tabf,35) coefm(i,j), mname(j), mform(j)
                END IF
 11             CONTINUE
      END IF

      IF (ng(i) .GT. 0) THEN
           DO 21 j = 1,ng(i)
                IF (coefg(i,j) .GT. 0.0d0) THEN
                     WRITE(tabf,36) coefg(i,j), gform(j)(1:20), gname(j)
                END IF
 21             CONTINUE
      END IF

      IF (na(i) .GT. 0) THEN
           DO 31 j = 1,na(i)
                IF (coefa(i,j) .GT. 0.0d0) THEN
                     WRITE(tabf,35) coefa(i,j), aname(j), aform(j)
                END IF
 31             CONTINUE
      END IF

      IF ((nw(i) .GT. 0) .AND. (coefw(i) .GT. 0.0d0)) THEN
           WRITE(tabf,55) coefw(i)
      END IF

***** write standard state properties, equation-of-state
***** parameters, and heat capacity coefficients  

      CALL wrtssp(i)

***** write header for property tabulation *****

      WRITE(tabf,65)
 65   FORMAT(//,' STANDARD STATE PROPERTIES OF THE REACTION', 
     1          ' AT ELEVATED TEMPERATURES AND PRESSURES ',//)

      CALL zero(i,nullrx)
  
      IF (nullrx) THEN
          WRITE(tabf,888)
 888      FORMAT(' CAUTION: INCOMPLETE DATA FOR ONE OR MORE SPECIES',//)
      END IF


      WRITE(tabf,75) isov(iopt,iplot), incv(iopt,iplot), 
     1               var3(iopt,iplot)
 75   FORMAT(50x,' DELTA G  ',1x,
     1        1x,' DELTA H  ',1x,
     1        1x,' DELTA S  ',1x,
     1        1x,' DELTA V  ',1x,
     1           ' DELTA Cp ',1x,/,
     1        2x,a10,2x,a10,2x,a10,
     1        2x,'  LOG K   ',1x,
     2        1x,' (kJ/mol) ',1x,
     3        1x,' (kJ/mol) ',1x,
     4        1x,' (J/mol/K)',1x,
     5        1x,' (J/bar)  ',1x,
     6        1x,' (J/mol/K)',1x,/,
     7        3(2x,10('-')),1x,6(1x,10('-'),1x))

      IF (xyplot .GT. 0) CALL pltrxn(i)

      END 

*****************************************************************

***** zero - Zero-out NULL values for reaction i to eliminate their
*****        contribution to standard molal properties at elevated 
*****        temperatures and pressures; set nullrx to .TRUE. if Gf 
*****        missing for mineral species or a1..4 for aqueous species.

      SUBROUTINE zero(i,nullrx)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

**JP**      PARAMETER (MAXMIN = 10, MAXGAS = 10, MAXAQS = 10, MXTRAN =  3,
**JP**     1           IABC   =  3, MAXRXN = 3000)
      PARAMETER (MAXMIN = 10, MAXGAS = 10, MAXAQS = 10, MXTRAN =  3,
     1           IABC   =  4, MAXRXN = 3000)
      CHARACTER*20  mname(MAXMIN), aname(MAXAQS)
      CHARACTER*30  mform(MAXMIN), aform(MAXAQS)
      CHARACTER*80  rtitle(MAXRXN)

      LOGICAL  m2reac(MAXRXN), nullrx

      INTEGER  nm(MAXRXN), na(MAXRXN), ng(MAXRXN), nw(MAXRXN),
     1         rec1m(MAXRXN,MAXMIN), rec1a(MAXRXN,MAXAQS), 
     2         rec1g(MAXRXN,MAXGAS)

      INTEGER ntran(MAXMIN)

      DOUBLE PRECISION coefm(MAXRXN,MAXMIN), coefa(MAXRXN,MAXAQS),
     1                 coefg(MAXRXN,MAXGAS), coefw(MAXRXN)

      DOUBLE PRECISION  Gfmin(MAXMIN), Hfmin(MAXMIN), VPrTrm(MAXMIN),
     1                  SPrTrm(MAXMIN), Alpha(MAXMIN), Kappa(MAXMIN),
     2                  KappaP(MAXMIN), KappaDP(MAXMIN),
     3                  MK1(IABC,MAXMIN), MK2(IABC,MAXMIN),
     4                  MK3(IABC,MAXMIN), MK4(IABC,MAXMIN),
     5                  Ttran(MXTRAN,MAXMIN), Htran(MXTRAN,MAXMIN),
     6                  Vtran(MXTRAN,MAXMIN), dPdTtr(MXTRAN,MAXMIN),
     7                  Tmaxm(MAXMIN), NumAtoms(MAXMIN), Tc(MAXMIN),
     8                  Smax(MAXMIN), Vmax(MAXMIN)

      DOUBLE PRECISION Gfaqs(MAXAQS), Hfaqs(MAXAQS), SPrTra(MAXAQS), 
     1                 a(4,MAXAQS), c(2,MAXAQS), 
     2                 wref(MAXAQS), chg(MAXAQS)

      COMMON /reac1/  rtitle
      COMMON /reac2/  coefm, coefa, coefg, coefw, nm, na, ng, nw,
     1                rec1m, rec1a, rec1g, m2reac

      COMMON /mnames/ mname, mform 
      COMMON /minref/ Gfmin, Hfmin, SPrTrm, VPrTrm, MK1, MK2, MK3, MK4, 
     1                Ttran, Htran, Vtran, dPdTtr, Tmaxm, ntran, Alpha,
     2                Kappa, KappaP, KappaDP, NumAtoms, Tc, Smax, Vmax

      COMMON /anames/ aname, aform
      COMMON /aqsref/ Gfaqs, Hfaqs, SPrTra, c, a, wref, chg
 
      COMMON /null/   XNULLM, XNULLA

      SAVE
                                                       
      nullrx = .FALSE.

      DO 10  j = 1,nm(i)
           IF (Gfmin(j) .EQ. XNULLM) THEN
                nullrx = .TRUE.
                Gfmin(j) = 0.0d0
                Hfmin(j) = 0.0d0
           END IF
           IF (ntran(j) .GT. 0) THEN
                DO 20 k = 1,ntran(j)
                     IF (Htran(k,j) .EQ. XNULLM) THEN
                          Htran(k,j) = 0.0d0
                     END IF
                     IF (Vtran(k,j) .EQ. XNULLM) THEN
                          Vtran(k,j) = 0.0d0
                     END IF
                     IF (dPdTtr(k,j) .EQ. XNULLM) THEN
                          dPdTtr(k,j) = 0.0d0
                     END IF
 20                  CONTINUE 
           END IF
 10        CONTINUE 

      DO 30 j = 1,na(i)
           IF (a(3,j) .EQ. XNULLA) THEN
                nullrx = .TRUE.
                DO 40 k = 1,4
                     a(k,i) = 0.0d0
 40                  CONTINUE  
           END IF
 30        CONTINUE  

      END

*****************************************************************

*** wrtssp - Write, to tabf, standard molal thermodynamic 
***          properties at 25 C and 1 bar, equation-of-state
***          parameters, and heat capacity coefficients for all
***          species in the i[th] reaction.

      SUBROUTINE wrtssp(i)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMIN = 10, MAXGAS = 10, MAXAQS = 10, MAXRXN = 3000,
**JP**     1           MXTRAN =  3, IABC   =  3, MAXMK  =  4, NPLOTF = 8,
     1           MXTRAN =  3, IABC   =  4, MAXMK  =  4, NPLOTF = 8,              
     2           ACON  = 4.184d-3,  BCON  = 4.184d2, CCON = 4.184d-3,
     3           DCON = 4.184d-3,
     3           A1CON = 1.0d1,  A2CON = 1.0d-2,
     4           A3CON = 1.0d0,  A4CON = 1.0d-4,
     5           C1CON = 1.0d0,  C2CON = 1.0d-4,
     6           WCON  = 1.0d-5)

      CHARACTER*4   incvar(2,3)
      CHARACTER*10  isov(2,3), incv(2,3), var3(2,3), isosat(2)
      CHARACTER*12  isovar(2,3)
      CHARACTER*20  mname(MAXMIN), gname(MAXGAS), aname(MAXAQS), 
     1              wname
      CHARACTER*30  mform(MAXMIN), gform(MAXGAS), aform(MAXAQS)
      CHARACTER*80  rtitle(MAXRXN)

      LOGICAL  m2reac(MAXRXN), nullrx

      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF), 
     1         univar, useLVS, epseqn, geqn,
     2         mapiso(2,3), mapinc(2,3), mapv3(2,3)

      INTEGER  nm(MAXRXN), na(MAXRXN), ng(MAXRXN), nw(MAXRXN),
     1         rec1m(MAXRXN,MAXMIN), rec1a(MAXRXN,MAXAQS), 
     2         rec1g(MAXRXN,MAXGAS), phaser(MAXMIN)

      INTEGER ntran(MAXMIN)

      DOUBLE PRECISION logKr, TPDref(4)

      DOUBLE PRECISION coefm(MAXRXN,MAXMIN), coefa(MAXRXN,MAXAQS),
     1                 coefg(MAXRXN,MAXGAS), coefw(MAXRXN)

      DOUBLE PRECISION  Gfmin(MAXMIN), Hfmin(MAXMIN), VPrTrm(MAXMIN),
     1                  SPrTrm(MAXMIN), Alpha(MAXMIN), Kappa(MAXMIN),
     2                  KappaP(MAXMIN), KappaDP(MAXMIN),
     3                  MK1(IABC,MAXMIN), MK2(IABC,MAXMIN),
     4                  MK3(IABC,MAXMIN), MK4(IABC,MAXMIN),
     5                  Ttran(MXTRAN,MAXMIN), Htran(MXTRAN,MAXMIN),
     6                  Vtran(MXTRAN,MAXMIN), dPdTtr(MXTRAN,MAXMIN),
     7                  Tmaxm(MAXMIN), NumAtoms(MAXMIN), Tc(MAXMIN),
     8                  Smax(MAXMIN), Vmax(MAXMIN)

      DOUBLE PRECISION Gfgas(MAXGAS), Hfgas(MAXGAS), VPrTrg(MAXGAS), 
     1                 SPrTrg(MAXGAS), MKg(IABC,MAXGAS), Tmaxg(MAXGAS)

      DOUBLE PRECISION Gfaqs(MAXAQS), Hfaqs(MAXAQS), SPrTra(MAXAQS), 
     1                 a(4,MAXAQS), c(2,MAXAQS), 
     2                 wref(MAXAQS), chg(MAXAQS)

      DOUBLE PRECISION Vmin(MAXMIN), Smin(MAXMIN), Cpmin(MAXMIN),
     2                 Hmin(MAXMIN), Gmin(MAXMIN)

      DOUBLE PRECISION Vgas(MAXGAS), Sgas(MAXGAS), Cpgas(MAXGAS), 
     2                 Hgas(MAXGAS), Ggas(MAXGAS)

      DOUBLE PRECISION Vaqs(MAXAQS), Saqs(MAXAQS), Cpaqs(MAXAQS),
     2                 Haqs(MAXAQS), Gaqs(MAXAQS),
     3                 VQterm(MAXAQS), SYterm(MAXAQS), CpXtrm(MAXAQS),
     4                 HYterm(MAXAQS), GZterm(MAXAQS)

      DOUBLE PRECISION mwH2O, Gftemp(MAXMIN), Hftemp(MAXMIN),
     1                 a3temp(MAXAQS)

      COMMON /stvars/ isosat, isovar, incvar
      COMMON /headmp/ isov, incv, var3
      COMMON /TPDmap/ mapiso, mapinc, mapv3
      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /refval/ mwH2O, R, Pref, Tref, ZPrTr, YPrTr
      COMMON /reac1/  rtitle
      COMMON /reac2/  coefm, coefa, coefg, coefw, nm, na, ng, nw,
     1                rec1m, rec1a, rec1g, m2reac
      COMMON /mnames/ mname, mform 
      COMMON /minref/ Gfmin, Hfmin, SPrTrm, VPrTrm, MK1, MK2, MK3, MK4, 
     1                Ttran, Htran, Vtran, dPdTtr, Tmaxm, ntran, Alpha,
     2                Kappa, KappaP, KappaDP, NumAtoms, Tc, Smax, Vmax
      COMMON /gnames/ gname, gform
      COMMON /gasref/ Gfgas, Hfgas, SPrTrg, VPrTrg, MKg, Tmaxg
      COMMON /anames/ aname, aform
      COMMON /aqsref/ Gfaqs, Hfaqs, SPrTra, c, a, wref, chg
      COMMON /H2Oss/  Dwss, Vwss, bewss, alwss, dalwss, Swss,
     1                Cpwss, Hwss, Gwss, Zwss, Qwss, Ywss, Xwss

      COMMON /minsp/ Vmin, Smin, Cpmin, Hmin, Gmin, phaser
      COMMON /gassp/ Vgas, Sgas, Cpgas, Hgas, Ggas 
      COMMON /aqsp/  Vaqs, Saqs, Cpaqs, Haqs, Gaqs
      COMMON /solvn/ VQterm, SYterm, CpXtrm, HYterm, GZterm
      COMMON /fmeq/  dVr, dSr, dCpr, dHr, dGr, logKr, dlogKT, dlogKP
      COMMON /null/  XNULLM, XNULLA

      SAVE

      DATA wname / 'H2O                 ' /



***** remove NULL contributions to standard state calculations

            nullrx = .FALSE.
      DO 456 iii = 1, nm(i) 
           Gftemp(iii) = 0.0d0
           IF (Gfmin(iii) .EQ. XNULLM) THEN
                nullrx = .TRUE.
                Gftemp(iii) = Gfmin(iii)
                Hftemp(iii) = Hfmin(iii)
                Gfmin(iii) = 0.0d0
                Hfmin(iii) = 0.0d0
           END IF
 456       CONTINUE

      DO 556 iii = 1, na(i) 
           a3temp(iii) = 0.0d0
           IF (a(3,iii) .EQ. XNULLA) THEN
                 nullrx = .TRUE.
                 a3temp(iii) = a(3,iii)
                 DO 557 jjj = 1,4
                      a(jjj,iii) = 0.0d0
 557                  CONTINUE
           END IF
 556       CONTINUE

***** calculate all reaction species heat capacities and reactant 
***** aqueous species standard partial molal volumes at 
***** 25 degC, 1 bar

      CALL reac92(i,Pref,Tref-273.15d0,Dwss,Vwss,bewss,alwss,dalwss,
     1            Swss,Cpwss,Hwss,Gwss,Zwss,Qwss,Ywss,Xwss,geqn)

***** return NULL contributions to faciltate blanking

      DO 457 iii = 1, nm(i) 
           IF (Gftemp(iii) .EQ. XNULLM) THEN
                Gfmin(iii) = Gftemp(iii)
                Hfmin(iii) = Hftemp(iii)
           END IF
 457       CONTINUE

      DO 656 iii = 1, na(i) 
           IF (a3temp(iii) .EQ. XNULLA) a(3,iii) = a3temp(iii) 
 656       CONTINUE

      WRITE(tabf,5)
  5   FORMAT(//,' STANDARD STATE PROPERTIES OF THE SPECIES AT',
     1          ' 25 DEG C AND 1 BAR')

      IF (nm(i) .GT. 0) THEN
           WRITE(tabf,6)
  6        FORMAT(//,42x,' ...... MINERALS ...... ',///,
     1              24x,'   DELTA G   ',1x,
     1               1x,'   DELTA H   ',1x,
     1               1x,'      S      ',1x,
     1               1x,'      V      ',1x,
     1               1x,'      Cp     ',1x,/,
     1               7x,'NAME',12x,
     1               1x,'   (kJ/mol)  ',1x,
     2               1x,'   (kJ/mol)  ',1x,
     3               1x,'  (J/mol/K)  ',1x,
     4               1x,'   (J/bar)   ',1x,
     5               1x,'  (J/mol/K)  ',/,
     6               2x,20('-'),2x,13('-'),2x,13('-'),2x,13('-'),
     7               2x,13('-'),2x,13('-'))

***** write mineral G, H, S, V, Cp at 25 C, 1 bar *****

           DO 10  j = 1,nm(i)
                IF (Gfmin(j) .EQ. XNULLM) THEN
                     WRITE(tabf,14) mname(j), SPrTrm(j)*4.184, 
     1                              VPrTrm(j)*1d-1, Cpmin(j)*4.184
 14                  FORMAT(2x,a20,30x,5x,f8.3,3x,4x,f8.3,4x,4x,f8.3,4x)
                ELSE
                     WRITE(tabf,15) mname(j), Gfmin(j)*4.184d-3,
     1                              Hfmin(j)*4.184d-3,
     2                              SPrTrm(j)*4.184, VPrTrm(j)*1d-1,
     3                              Cpmin(j)*4.184
 15                  FORMAT(2x,a20,3x,f10.3,2x,3x,f10.3,2x,5x,f8.3,3x,
     1                      4x,f8.3,4x,4x,f8.3,4x)
                END IF
 10             CONTINUE
                                                   
           WRITE(tabf,7)
  7        FORMAT(//,24x,'HOLLAND & POWELL COEFFICIENTS',
     1               32x,'LANDAU PARAMETERS',/,
     1                7x,'NAME',11x,
     1                1x,'  a(10**0)  ',
     2                   '  b(10**5)  ',
     3                   '  c(10**0) ',1x
     3                   '  d(10**0) ',2x
     4                   ' T limit (C)',
     5                1x,'     Tc (C)   ',1x,
     6                1x,'Smax (J/mol/K)',1x,
     7                1x,' Vmax (J/bar) ',1x,/,
     8                2x,20('-'),
     9                2x,10('-'),1x,
     1                1x,10('-'),1x,
     1                1x,10('-'),1x,
     1                1x,10('-'),1x,
     2                1x,11('-'),1x,
     3                1x,13('-'),1x,
     4                1x,14('-'),1x,
     5                1x,14('-'),1x)

           DO 21  j = 1,nm(i)

***** write mineral heat capacity coefficients
***** a, b, c and phase transition T, H, V, dPdT

                IF (ntran(j) .EQ. 0) THEN

                     WRITE(tabf,16) mname(j), MK1(1,j)*ACON, 
     1                              MK1(2,j)*BCON, MK1(3,j)*CCON, 
     3                              MK1(4,j)*DCON, Tmaxm(j)-273.15d0
 16                  FORMAT(2x,a20,1x,f9.4,2x,1x,f9.4,2x,1x,f11.3,3x,
     1                      f9.4,1x,5x,f6.0)

                ELSE
                            WRITE(tabf,219) mname(j), MK1(1,j)*ACON, 
     1                           MK1(2,j)*BCON, MK1(3,j)*CCON, 
     2                           MK1(4,j)*DCON,
     3                           Tmaxm(j)-273.15d0, Tc(j)-273.15d0,
     4                           Smax(j), Vmax(j)
 219                       FORMAT(2x,a20,1x,f9.4,2x,1x,f9.4,1x,2x,f11.3,
     1                             3x,f9.4,1x,5x,f6.0,5x,4x,f8.2,4x,
     2                             4x,f7.3,5x,f9.4)
                END IF
 21             CONTINUE
           WRITE(tabf,123)
  123        FORMAT(//,24x,'PRESSURE TERMS',/,
     1                7x,'NAME',11x,
     1                1x,'alpha(10**5)',
     2                   '   kappa    ',
     3                   "   kappa'   ",
     3                   "   kappa''  ",
     4                   ' # of Atoms ',/
     8                2x,20('-'),
     9                2x,10('-'),1x,
     1                1x,10('-'),1x,
     1                1x,10('-'),1x,
     1                1x,10('-'),1x,
     2                1x,10('-'),1x)
           DO 19  j = 1,nm(i)
                     WRITE(tabf,18) mname(j), Alpha(j)*1d5, 
     1                              Kappa(j), KappaP(j), 
     3                              KappaDP(j), NumAtoms(j)
 18                  FORMAT(2x,a20,3x,f5.2,6x,1x,f7.1,4x,1x,f5.2,5x,3x,
     1                      f7.4,4x,f7.3)
 19             CONTINUE


      END IF

      IF (ng(i) .GT. 0) THEN
           WRITE(tabf,8)
  8        FORMAT(///,42x,' ...... GASES ...... ',///,
     1              24x,'   DELTA G   ',1x,
     1               1x,'   DELTA H   ',1x,
     1               1x,'      S      ',1x,
     1               1x,'      V      ',1x,
     1               1x,'      Cp     ',1x,/,
     1               7x,'NAME',12x,
     1               1x,'  (kJ/mol)   ',1x,
     2               1x,'  (kJ/mol)   ',1x,
     3               1x,'  (J/mol/K)  ',1x,
     4               1x,'    (J/bar)  ',1x,
     5               1x,'  (J/mol/K)  ',/,
     6               2x,20('-'),2x,13('-'),2x,13('-'),2x,13('-'),
     7               2x,13('-'),2x,13('-'))

***** write gas G, H, S, V, Cp at 25 C, 1 bar and
***** heat capacity coefficients a, b, c

           DO 30  j  = 1,ng(i)
                WRITE(tabf,17) gname(j), Gfgas(j)*4.184d-3, 
     1                         Hfgas(j)*4.184d-3, 
     2                         SPrTrg(j)*4.184, VPrTrg(j)*1d-1, 
     3                         Cpgas(j)*4.184
 17             FORMAT(2x,a20,
     1                 3x,f10.3,2x,
     2                 3x,f10.3,2x,
     3                 5x,f8.3,3x,
     4                 7x,f3.1,6x,
     5                 4x,f6.3,6x)
 30             CONTINUE


           WRITE(tabf,9)
  9        FORMAT(//,29x,'HOLLAND-POWELL COEFFICIENTS',/,
     1                7x,'NAME',11x,
     1                1x,'  a(10**0)  ',
     2                   '  b(10**5)  ',
     3                   '  c(10**0) ',
     3                   '  d(10**0) ',
     4                   ' T limit (C)',/,
     4                2x,20('-'),
     5                2x,10('-'),1x,
     5                1x,10('-'),1x,
     6                1x,10('-'),1x,
     6                1x,10('-'),1x,
     7                1x,10('-'))

           DO 31  j = 1,ng(i)
                WRITE(tabf,16) gname(j), MKg(1,j)*ACON, 
     1                         MKg(2,j)*BCON, MKg(3,j)*CCON,
     2                         MKg(4,j)*DCON,
     3                         Tmaxg(j)-273.15d0
 31             CONTINUE
      END IF

      IF ((na(i) .GT. 0) .OR. (nw(i) .GT. 0)) THEN
           WRITE(tabf,46)
 46        FORMAT(///,36x,' ...... AQUEOUS SPECIES ...... ',///,
     1              24x,'   DELTA G   ',1x,
     1               1x,'   DELTA H   ',1x,
     1               1x,'      S      ',1x,
     1               1x,'      V      ',1x,
     1               1x,'      Cp     ',1x,/,
     1               7x,'NAME',12x,
     1               1x,'   (kJ/mol)  ',1x,
     2               1x,'   (kJ/mol)  ',1x,
     3               1x,'   (J/mol/K) ',1x,
     4               1x,'    (J/bar)  ',1x,
     5               1x,'   (J/mol/K) ',/,
     6               2x,20('-'),2x,13('-'),2x,13('-'),2x,13('-'),
     7               2x,13('-'),2x,13('-'))

***** write aqueous species G, H, S, V, Cp at 25 C, 1 bar

           DO 40  j = 1,na(i)
                WRITE(tabf,79) aname(j), Gfaqs(j)*4.184d-3, 
     1                         Hfaqs(j)*4.184d-3, 
     2                         SPrTra(j)*4.184, Vaqs(j)*1d-1, 
     3                         Cpaqs(j)*4.184
 79             FORMAT(2x,a20,
     1                 3x,f10.3,2x,
     2                 3x,f10.3,2x,
     3                 5x,f8.3,3x,
     4                 5x,f6.3,5x,
     5                 4x,f9.3,3x)
 40             CONTINUE
           IF (nw(i) .GT. 0) THEN

***** write H2O G, H, S, V, Cp at 25 C, 1 bar ******

                WRITE(tabf,79) wname, Gwss*4.184d-3, Hwss*4.184d-3, 
     1                         Swss*4.184, Vwss*1d-1, Cpwss*4.184
           END IF


           IF (na(i) .GT. 0) THEN
                WRITE(tabf,56)
 56             FORMAT(//,50x,'EQUATION-OF-STATE COEFFICIENTS',/,
     1                     7x,'NAME',11x,
     1                     2x,' a1(10**1)  ',1x,
     1                     1x,' a2(10**-2) ',1x,
     1                     1x,' a3(10**0)  ',1x,
     1                     1x,' a4(10**-4) ',1x,
     2                     1x,' c1(10**0)  ',1x,
     3                     1x,' c2(10**-4) ',
     1                     2x,'omega(10**-5)',/,
     4                     2x,20('-'),
     5                     2x,12('-'),1x,
     6                     1x,12('-'),1x,
     6                     1x,12('-'),1x,
     6                     1x,12('-'),1x,
     6                     1x,12('-'),1x,
     6                     1x,12('-'),1x,
     6                     1x,13('-'))

***** write aqueous species equation-of-state parameters
***** wref, c[1..4], a[1..2]

                DO 50  j = 1,na(i)
                     IF (a(3,j) .EQ. XNULLA) THEN
                          WRITE(tabf,64) aname(j), 
     2                    c(1,j)*C1CON, c(2,j)*C2CON, wref(j)*WCON 
 64                       FORMAT(2x,a20,1x,56x,3(3x,f8.4,3x))
                     ELSE
                          WRITE(tabf,65) aname(j), 
     1                    a(1,j)*A1CON, a(2,j)*A2CON, a(3,j)*A3CON,
     2                    a(4,j)*A4CON, c(1,j)*C1CON, c(2,j)*C2CON,
     3                    wref(j)*WCON 
 65                       FORMAT(2x,a20,1x,7(3x,f8.4,3x))
                     END IF

 50                  CONTINUE
           END IF
      END IF

***** write reaction properties at 25 C, 1 bar

      WRITE(tabf,74) 
 74   FORMAT(///,' STANDARD STATE PROPERTIES OF THE REACTION AT',
     1           ' 25 DEG C AND 1 BAR',//)

      IF (nullrx) THEN
          WRITE(tabf,888)
 888      FORMAT(' CAUTION: INCOMPLETE DATA FOR ONE OR MORE SPECIES',//)
      END IF

      WRITE(tabf,75) isov(iopt,iplot), incv(iopt,iplot), 
     1               var3(iopt,iplot)
 75   FORMAT(50x,' DELTA G  ',1x,
     1        1x,' DELTA H  ',1x,
     1        1x,' DELTA S  ',1x,
     1        1x,' DELTA V  ',1x,
     1           ' DELTA Cp ',1x,/,
     1        2x,a10,2x,a10,2x,a10,
     1        2x,'  LOG K   ',1x,
     2        1x,' (kJ/mol) ',1x,
     3        1x,' (kJ/mol) ',1x,
     4        1x,' (J/mol/K)',1x,
     5        1x,'  (J/bar) ',1x,
     6        1x,' (J/mol/K)',1x,/,
     7        3(2x,10('-')),1x,6(1x,10('-'),1x))

      TPDref(1) = Tref - 273.15d0 
      TPDref(2) = Pref 
      TPDref(3) = Dwss 
      TPDref(4) = Dwss 

      WRITE(tabf,85) TPDref(mapiso(iopt,iplot)),
     1               TPDref(mapinc(iopt,iplot)),
     2               TPDref(mapv3(iopt,iplot)), 
     3               logKr, dGr*4.184d-3, dHr*4.184d-3
     4               , dSr*4.184, dVr*1d-1, dCpr*4.184

  85  FORMAT(3(2x,f9.3,1x),
     1         1x,f10.3,1x,
     2         1x,f10.3,1x,
     3         1x,f10.3,1x,
     4         1x,f9.3,2x,
     5         1x,f9.3,2x,
     6         1x,f11.3,/)

      END

******************************************************************

*** pltrxn - Write header information and reaction titles 
***          to plot files.

      SUBROUTINE pltrxn(ireac)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NPLOTF = 8, MAXRXN = 3000, MXRPLT = 10)

      CHARACTER*80 rtitle(MAXRXN)
      CHARACTER*20 namecf, namerf, nametf, namepf(NPLOTF)
      CHARACTER*3  rstamp(MXRPLT)
      INTEGER      rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     1             xyplot, end, rlen(2) 
      LOGICAL      openf

      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /fnames/ namecf, namerf, nametf, namepf
      COMMON /plottr/ xyplot, end, nplots
      COMMON /reac1/  rtitle

      SAVE

      DATA rlen / 80, 249 /

      DATA rstamp / 'R01', 'R02', 'R03', 'R04', 'R05',
     1              'R06', 'R07', 'R08', 'R09', 'R10' /


*** close files if necessary
      IF ((xyplot .EQ. 2) .AND. (ireac .GT. 1)) THEN
           DO 10 i = 1,nplots
		CLOSE(UNIT=plotf(i))
                namepf(i)(end+1:end+3) = rstamp(ireac)
 10             CONTINUE
      END IF

*** if necessary, open files and write header  
      IF ((xyplot .EQ. 2) .OR. (ireac .EQ. 1)) THEN
	   DO 20 i = 1,nplots
                IF (openf(wterm,plotf(i),namepf(i),2,1,1,
     1               rlen(xyplot))) THEN
                     CALL plttop(i)  
                ELSE
                     WRITE(wterm,*) ' cannot open plot file ',i
                END IF 
 20             CONTINUE
      END IF

*** write reaction number and title
      DO 30 i = 1,nplots 
           WRITE(plotf(i),40) ireac, rtitle(ireac)
 40        FORMAT(//,' REACTION ',i2,
     1             /,' TITLE: ',a80)
 30        CONTINUE

      END

**********************************************************************

*** plttop - Write global banner to plot file i.

      SUBROUTINE plttop(i)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NPLOTF = 8)

      CHARACTER*20 pfname, namecf, namerf, nametf, namepf(NPLOTF),
     1             nosave
      LOGICAL      savecf, saverf
      INTEGER      rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     1             xyplot, end 

      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /fnames/ namecf, namerf, nametf, namepf
      COMMON /saveif/ savecf, saverf
      COMMON /plottr/ xyplot, end, nplots
      COMMON /dapron/ pfname

      SAVE

      DATA nosave / 'file not saved      ' /


      WRITE(plotf(i),20) 
 20   FORMAT(/,' ***** SUPCRTBL: input/output specifications for',
     1         ' this run',/)
      IF (savecf) THEN
           WRITE(plotf(i),30) namecf
 30        FORMAT(  '            USER-SPECIFIED  CON FILE  containing ',
     1            /,'            T-P-D grid & option switches: ',a20,/)
      ELSE
           WRITE(plotf(i),30) nosave
      END IF
      IF (saverf) THEN
           WRITE(plotf(i),40) namerf
 40   FORMAT(  '            USER-SPECIFIED  RXN FILE  containing ',
     1       /,'            chemical reactions: ',a20,/)
      ELSE
           WRITE(plotf(i),40) nosave
      END IF

      WRITE(plotf(i),50) pfname
 50   FORMAT(  '            THERMODYNAMIC DATABASE: ',a20,/)
      WRITE(plotf(i),60) nametf
 60   FORMAT('            SUPCRT-GENERATED  TAB FILE  containing ',
     1     /,'            tabulated reaction properties: ',a20)

      CALL wrtop2(i)

      END

*********************************************************************

*** report - Report computed reaction properties.

      SUBROUTINE report(ireac, iso, inc, TPD, TPDtrn, rptran, ptrans,
     1                  dVr, dSr, dCpr, dHr, dGr, logKr, 
     2                  lvdome, H2Oerr, Kfound)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMIN = 10, MAXAQS = 10, MAXGAS = 10, MAXRXN = 3000,
**JP**     1           MXTRAN =  3, IABC   =  3, NPLOTF  = 8, MAXODD = 21)
     1           MXTRAN =  3, IABC   =  4, NPLOTF  = 8, MAXODD = 21)
      CHARACTER*1   PT(2)
      CHARACTER*5   pprop(NPLOTF)
      CHARACTER*10  isov(2,3), incv(2,3), var3(2,3)
      CHARACTER*20  mname(MAXMIN)
      CHARACTER*30  mform(MAXMIN)
      CHARACTER*80  rtitle(MAXRXN)

      LOGICAL  m2reac(MAXRXN), xall, xHSVCp, xCp,
     1         MKwarn, MKdone, Pwarn, Pdone, Kfound, Klost
      LOGICAL  rptran, lvdome, H2Oerr, EQ3run

      INTEGER  phaser(MAXMIN), ptrans(MAXMIN), ntran(MAXMIN),
     1         xyplot, end

      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     1         mapiso(2,3), mapinc(2,3), mapv3(2,3),
     2         univar, useLVS, epseqn, geqn

      INTEGER  nm(MAXRXN), na(MAXRXN), ng(MAXRXN), nw(MAXRXN),
     1         rec1m(MAXRXN,MAXMIN), rec1a(MAXRXN,MAXAQS), 
     2         rec1g(MAXRXN,MAXGAS)

      DOUBLE PRECISION  coefm(MAXRXN,MAXMIN), coefa(MAXRXN,MAXAQS),
     1                  coefg(MAXRXN,MAXGAS), coefw(MAXRXN)

      DOUBLE PRECISION TPD(3), TPDtrn(MAXMIN,MXTRAN,3), logKr

      DOUBLE PRECISION Vmin(MAXMIN), Smin(MAXMIN), Cpmin(MAXMIN),
     2                 Hmin(MAXMIN), Gmin(MAXMIN)

      DOUBLE PRECISION  Gfmin(MAXMIN), Hfmin(MAXMIN), VPrTrm(MAXMIN),
     1                  SPrTrm(MAXMIN), Alpha(MAXMIN), Kappa(MAXMIN),
     2                  KappaP(MAXMIN), KappaDP(MAXMIN),
     3                  MK1(IABC,MAXMIN), MK2(IABC,MAXMIN),
     4                  MK3(IABC,MAXMIN), MK4(IABC,MAXMIN),
     5                  Ttran(MXTRAN,MAXMIN), Htran(MXTRAN,MAXMIN),
     6                  Vtran(MXTRAN,MAXMIN), dPdTtr(MXTRAN,MAXMIN),
     7                  Tmaxm(MAXMIN), NumAtoms(MAXMIN), Tc(MAXMIN),
     8                  Smax(MAXMIN), Vmax(MAXMIN)

      DOUBLE PRECISION isomin, isomax, isoinc, Kmin, Kmax, Kinc, 
     1                 oddv1(MAXODD), oddv2(MAXODD)


      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /grid/   isomin, isomax, isoinc, v2min, v2max, v2inc,
     1                oddv1, oddv2, Kmin, Kmax, Kinc, niso, nv2, nlogK
      COMMON /TPDmap/ mapiso, mapinc, mapv3
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /reac1/  rtitle
      COMMON /reac2/  coefm, coefa, coefg, coefw, nm, na, ng, nw,
     1                rec1m, rec1a, rec1g, m2reac
      COMMON /mnames/ mname, mform 
      COMMON /minref/ Gfmin, Hfmin, SPrTrm, VPrTrm, MK1, MK2, MK3, MK4, 
     1                Ttran, Htran, Vtran, dPdTtr, Tmaxm, ntran, Alpha,
     2                Kappa, KappaP, KappaDP, NumAtoms, Tc, Smax, Vmax
      COMMON /minsp/  Vmin, Smin, Cpmin, Hmin, Gmin, phaser
      COMMON /headmp/ isov, incv, var3
      COMMON /EQ36/   EQ3run
      COMMON /plottr/ xyplot, end, nplots 

      SAVE 

      DATA PT            / 'P', 'T' /
      DATA MKdone, Pdone / 2*.FALSE. /
      DATA pprop         / 'logK ', 'delG ', 'delH ', 'delS ',
     1                     'delCp', 'delV ', 'dvar1', 'dvar2' / 

      IF ((inc .EQ. 1)) THEN 
           MKdone = .FALSE.
           Pdone  = .FALSE.
           WRITE(tabf,10)
 10        FORMAT()
           IF (xyplot .EQ. 1) THEN 
                CALL pltln1(TPD)
           END IF
      END IF

      IF ((xyplot .EQ. 2) .AND. (iso .EQ. 1) .AND. (inc .EQ. 1)) THEN
           CALL pltln1(TPD)
      END IF

      IF (lvdome) THEN
*****      T,DH2O location falls within liquid-vapor dome 
           WRITE(tabf,15) TPD(mapiso(iopt,iplot)),
     1                    TPD(mapinc(iopt,iplot))
 15        FORMAT(2(2x,f9.3,1x),23x,
     1     ' T-DH2O FALLS WITHIN LIQ-VAP T-DH2O DOME ',/)
      END IF

      IF (H2Oerr) THEN
*****      T-P-d beyond fluid H2O equation of state 
           WRITE(tabf,16) TPD(mapiso(iopt,iplot)),
     1                    TPD(mapinc(iopt,iplot))
 16        FORMAT(2(2x,f9.3,1x),14x,' *** BEYOND RANGE OF',
     1     ' APPLICABILITY OF H2O EQUATION OF STATE **',/)
      END IF

      Klost = (univar .EQ. 1) .AND. (.NOT. Kfound)
      IF (Klost) THEN
           WRITE(tabf,17) TPD(mapiso(iopt,iplot)), logKr, 
     1                    v2min, PT(2/iplot), v2max
 17        FORMAT(2x,f9.3,26x,f10.3,4x,
     1            ' LOG K NOT FOUND: ',f10.3,' <= ',a1,' <= ',f10.3,/)
      END IF

      IF (rptran) THEN
           DO 20 imin = 1,nm(ireac) 
                DO 25 itran = ptrans(imin),1,-1
                     IF ((iplot .EQ. 2) .AND. (isat .EQ. 0)) THEN
***                       isotherms(P,D): phase c -> b -> a
                          iptnum = phaser(imin)
                     ELSE
***                       iso[bars,chores](T): phase a -> b -> c
                          iptnum = phaser(imin)-itran
                     END IF
                     WRITE(tabf,30) 
     1               TPDtrn(imin,itran,mapiso(iopt,iplot)),
     2               TPDtrn(imin,itran,mapinc(iopt,iplot)),
     3               TPDtrn(imin,itran,mapv3(iopt,iplot)-isat),
     4               iptnum, mname(imin)
  30                  FORMAT(3(2x,f9.3,1x),11x,
     1               ' PHASE TRANSITION #',i1,' for ',
     2               a20,/)
  25                 CONTINUE
  20            CONTINUE
       END IF

*******************************************************************
*** Assignment of "blank" variables below facilitates full
*** reporting of reaction properties beyond certain limits
*** of certain equations when the CALL to SUBR blank is
*** commented-out for the development version.  

      xall   = .FALSE.
      xHSVCp = .FALSE.
      xCp    = .FALSE.
      MKwarn = .FALSE.
      Pwarn =  .FALSE.

*** SUBR blanks to be called for distribution version;
*** not called for development version

      IF (.NOT. (lvdome .OR. H2Oerr)) 
     1     CALL blanks(TPD(1),TPD(2),TPD(3),m2reac(ireac),nm(ireac),
     2            ng(ireac),na(ireac),xall,xHSVCp,xCp,MKwarn,Pwarn)

***** write reaction properties to plot files

      IF (xyplot .GT. 0) 
     1     CALL pltrep(TPD,iso,inc,logKr,dGr,dHr,dSr,dCpr,dVr,
     2            lvdome,H2Oerr,Klost,xall,xHSVCp,xCp)

      IF (lvdome .OR. H2Oerr .OR. Klost) RETURN

*******************************************************************

***** write reaction properties to tab file

      IF (xall) THEN
***** P-T beyond validity limits of 
***** aqueous species equations
           WRITE(tabf,40) TPD(mapiso(iopt,iplot)), 
     1                    TPD(mapinc(iopt,iplot)), 
     2                    TPD(mapv3(iopt,iplot)-isat)
  40       FORMAT(3(2x,f9.3,1x),2x,' *** BEYOND RANGE OF',
     1     ' APPLICABILITY OF AQUEOUS SPECIES EQNS ***',/)
           RETURN
      END IF

      IF (xHSVCp) THEN
***** P-T within region where only Gibbs free energy can 
***** be computed for charged aqueous species
           WRITE(tabf,45) TPD(mapiso(iopt,iplot)), 
     1                    TPD(mapinc(iopt,iplot)), 
     2                    TPD(mapv3(iopt,iplot)-isat), logKr
     3                    , dGr*4.184d-3
  45       FORMAT(3(2x,f9.3,1x),
     1              1x,f10.3,1x,
     2              1x,f10.3,3x,
     1    '*** DELTA G ONLY (CHARGED AQUEOUS SPECIES) ***',/)
           RETURN
      END IF

      IF (MKwarn .AND. (.NOT. MKdone)) THEN
***** beyond temperature limit of Maier-Kelly Cp coefficients; 
***** issue warning to this effect. 
           MKdone = .TRUE.
           WRITE(tabf,50) TPD(mapiso(iopt,iplot)), 
     1                    TPD(mapinc(iopt,iplot)), 
     2                    TPD(mapv3(iopt,iplot)-isat) 
  50       FORMAT(3(2x,f9.3,1x),2x,' *** CAUTION: BEYOND T LIMIT', 
     1     ' OF CP COEFFS FOR A MINERAL OR GAS ***',/)
      END IF

      IF (Pwarn .AND. (.NOT. Pdone)) THEN
***** beyond qualitative pressure limit of mineral/gas calculations; 
***** issue warning to this effect. 
           Pdone = .TRUE.
           WRITE(tabf,55) TPD(mapiso(iopt,iplot)), 
     1                    TPD(mapinc(iopt,iplot)), 
     2                    TPD(mapv3(iopt,iplot)-isat) 
  55       FORMAT(3(2x,f9.3,1x),2x,' *** CAUTION: BEYOND P LIMIT', 
     1     ' OF APPROXIMATIONS IN MINERAL/GAS EQNS ***',/)
      END IF

      IF (xCp) THEN
***** within +/- 25 degC of mineral phase transition; 
***** report all property values except dCpr 
           WRITE(tabf,60) TPD(mapiso(iopt,iplot)), 
     1                    TPD(mapinc(iopt,iplot)), 
     2                    TPD(mapv3(iopt,iplot)-isat), 
     3                    logKr, dGr*4.184d-3, dHr*4.184d-3,
     4                    dSr*4.184, dVr*1d-1
  60       FORMAT(3(2x,f9.3,1x),
     1         1x,f10.3,1x,
     2         1x,f10.3,1x,
     3         1x,f10.3,1x,
     4         1x,f9.1,2x,
     5         1x,f9.3,2x,'  TRANSITION',/)
           RETURN
      END IF

***** .NOT. (xall .OR. xHSVCp .OR. xCp); 
***** report all property values
      WRITE(tabf,70) TPD(mapiso(iopt,iplot)), 
     1               TPD(mapinc(iopt,iplot)), 
     2               TPD(mapv3(iopt,iplot)-isat), 
     3               logKr, dGr*4.184d-3, dHr*4.184d-3, 
     4               dSr*4.184, dVr*1d-1, dCpr*4.184
  70  FORMAT(3(2x,f9.3,1x),
     1         1x,f10.3,1x,
     2         1x,f10.3,1x,
     3         1x,f10.3,1x,
     4         1x,f9.3,2x,
     5         1x,f9.3,2x,
     6         1x,f11.3,/)

      RETURN
      END

*******************************************************************

*** pltln1 - Write top line of next property block.

      SUBROUTINE pltln1(TPD)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NPLOTF = 8, MXIPLT = 11, MAXODD = 21, MAXINC = 21)

      DOUBLE PRECISION TPD(3), isomin, isomax, isoinc, Kmin, Kmax, Kinc,
     1                 oddv1(MAXODD), oddv2(MAXODD),
     2                 pbuff(MAXINC,MXIPLT,NPLOTF)
      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     1         mapiso(2,3), mapinc(2,3), mapv3(2,3),
     2         univar, useLVS, epseqn, geqn, xyplot, end
      LOGICAL EQ3run

      CHARACTER*7  isov2(2,3), incv2(2,3), var32(2,3)
      CHARACTER*10 isov(2,3), incv(2,3), var3(2,3)
      CHARACTER*20 pprop(NPLOTF)

      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /grid/   isomin, isomax, isoinc, v2min, v2max, v2inc,
     2                oddv1, oddv2, Kmin, Kmax, Kinc, niso, nv2, nlogK
      COMMON /headmp/ isov, incv, var3
      COMMON /TPDmap/ mapiso, mapinc, mapv3
      COMMON /EQ36/   EQ3run
      COMMON /plottr/ xyplot, end, nplots
      COMMON /pbuffr/ pbuff

      SAVE

      DATA pprop / 'logK         ', 'delG (kJ/mol)', 'delH (kJ/mol)', 
     1          'delS (J/mol/K)', 'delCp (J/mol/K)', 'delV (J/bar)', 
     2             'dvar1        ', 'dvar2        ' /

      DATA isov2  / 'D(g/cc)', 'P(bars)', 3*'T(degC)', 'P(bars)' /

      DATA incv2  / 2*'T(degC)', 'D(g/cc)', 'P(bars)', 'P(bars)', 
     1                'T(degC)' /

      DATA var32  / 'P(bars)', 'D(g/cc)', 'P(bars)', 3*'D(g/cc)' /


      IF (xyplot .EQ. 1) THEN
           IF (univar .EQ. 1) THEN
                RETURN
           END IF
           DO 10 i = 1,nplots   
                IF ((isat .EQ. 0) .AND. (.NOT. EQ3run)) THEN
                     IF (i .EQ. 7) THEN
                          WRITE(plotf(i),20) isov(iopt,iplot),
     1                    TPD(mapiso(iopt,iplot)),
     2                    incv(iopt,iplot), var3(iopt,iplot)
 20                       FORMAT(/,1x,a10,' = ',f10.4,' ; ',
     1                           a10,', ',a10,/) 
                     ELSE
                          WRITE(plotf(i),25) isov(iopt,iplot),
     1                    TPD(mapiso(iopt,iplot)),
     2                    incv(iopt,iplot), pprop(i) 
 25                       FORMAT(/,1x,a10,' = ',f10.4,' ; ',
     1                           a10,', ',a13,/) 
                     END IF
                ELSE
                     IF (i .LT. 7) THEN
                          WRITE(plotf(i),30) isov(iopt,iplot), pprop(i) 
 30                       FORMAT(/,1x,a10,', ',a13,/) 
                     END IF
                     IF (i .EQ. 7) THEN
                          WRITE(plotf(i),35) isov(iopt,iplot),
     1                                       incv(iopt,iplot) 
 35                       FORMAT(/,1x,a10,', ',a10,/) 
                     END IF
                     IF (i .EQ. 8) THEN
                          WRITE(plotf(i),35) isov(iopt,iplot),
     1                                       var3(iopt,iplot) 
                     END IF
                END IF
 10             CONTINUE
      ELSE
           DO 40 i = 1,nplots
                IF ((isat .EQ. 1) .OR. EQ3run) THEN
                     WRITE(plotf(i),50) isov2(iopt,iplot), 
     1               incv2(iopt,iplot), var32(iopt,iplot), 
     1               (pprop(j), j = 1,6) 
 50                  FORMAT(/,3(a7,2x),6(a13,2x)) 
                     GO TO 40
                END IF
                IF (univar .EQ. 0) THEN
                     IF (i .LT. 7) THEN
                          WRITE(plotf(i),60) isov2(iopt,iplot)(1:1),
     1                    incv2(iopt,iplot)(1:1), pprop(i)
 60                       FORMAT(/,1x,a1,'-',a1,' grid: ',a13)
                     ELSE
                          WRITE(plotf(i),65) isov2(iopt,iplot)(1:1),
     1                    incv2(iopt,iplot)(1:1), var32(iopt,iplot) 
 65                       FORMAT(/,1x,a1,'-',a1,' grid: ',a7)
                     END IF
                     WRITE(plotf(i),70) incv2(iopt,iplot),
     1               (isov2(iopt,iplot), isomin+(j-1)*isoinc,
     2               j=1,MIN(niso,MXIPLT))
***                  MXIPLT = 11 
 70                  FORMAT(/,a7,11(2x,a7,' =',e11.4))
                ELSE
                     WRITE(plotf(i),80) isov2(iopt,iplot)(1:1),
     1                                  incv2(iopt,iplot)
 80                  FORMAT(/,1x,a1,'-logK grid: ',a7)
                     WRITE(plotf(i),90) isov2(iopt,iplot),
     1                    (incv2(iopt,iplot)(1:1), 
     2                     Kmin+(k-1)*Kinc, k=1,MIN(nlogK,MXIPLT))
***                  MXIPLT = 11 
 90                  FORMAT(/,a7,11(2x,a1,'(logK =',e11.4,')'))
                END IF
 40             CONTINUE
	 
*** initialize property buffer
           DO 100 k = 1,nplots
                DO 100 j = 1,MIN(niso,MXIPLT)
                     DO 100 i = 1,nv2
                          pbuff(i,j,k) = 0.0d0
 100                      CONTINUE

      END IF

      END

***********************************************************************
 
*** pltrep - Report calculated property values to plot files or buffers

      SUBROUTINE pltrep(TPD,iso,inc,logKr,dGr,dHr,dSr,dCpr,dVr,
     1                  lvdome,H2Oerr,Klost,xall,xHSVCp,xCp)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXINC = 21, MAXODD = 21, MXIPLT = 11, NPLOTF = 8)

      LOGICAL lvdome, H2Oerr, Klost, xall, xHSVCp, xCp, EQ3run

      CHARACTER*10 isov(2,3), incv(2,3), var3(2,3)

      INTEGER  rterm, wterm, reacf, pronf, tabf, plotf(NPLOTF),
     1         mapiso(2,3), mapinc(2,3), mapv3(2,3), univar, 
     2         useLVS, epseqn, geqn, xyplot, end

      DOUBLE PRECISION TPD(3), pbuff(MAXINC,MXIPLT,NPLOTF),
     1                 isomin, isomax, isoinc, Kmin, Kmax, Kinc, logKr,
     2                 oddv1(MAXODD), oddv2(MAXODD), satmin(2)

      COMMON /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
      COMMON /TPDmap/ mapiso, mapinc, mapv3
      COMMON /headmp/ isov, incv, var3
      COMMON /icon/   isat, iopt, iplot, univar, noninc,
     1                useLVS, epseqn, geqn
      COMMON /grid/   isomin, isomax, isoinc, v2min, v2max, v2inc,
     2                oddv1, oddv2, Kmin, Kmax, Kinc, niso, nv2, nlogK
      COMMON /satend/ satmin
      COMMON /EQ36/   EQ3run
      COMMON /plottr/ xyplot, end, nplots 
      COMMON /pbuffr/ pbuff

      SAVE

***************************  xyplot = 1  *******************************
			WRITE(*,31) xyplot, univar
  31  FORMAT(2(1x,i14.7))
      IF ((xyplot .EQ. 1) .AND. (lvdome .OR. H2Oerr .OR. xall)) RETURN

      IF ((xyplot .EQ. 1) .AND. (univar .EQ. 1)) THEN
           IF (Klost) THEN
                pbuff(inc,iso,1) = 0.0d0
           ELSE
                pbuff(inc,iso,1) = TPD(mapinc(iopt,iplot))
           END IF
           GO TO 20
      END IF

      IF ((xyplot .EQ. 1) .AND. xHSVCp) THEN
           IF ((isat .EQ. 1) .OR. EQ3run) THEN
                WRITE(plotf(1),10) TPD(mapiso(iopt,iplot)), logKr
                WRITE(plotf(2),10) TPD(mapiso(iopt,iplot)), dGr*4.184d-3
                WRITE(plotf(7),10) TPD(mapiso(iopt,iplot)), 
     1                             TPD(mapinc(iopt,iplot)) 
                WRITE(plotf(8),10) TPD(mapiso(iopt,iplot)), 
     1                             TPD(mapv3(iopt,iplot)-isat) 
  10            FORMAT(2(1x,e14.7))
           END IF
           IF (univar .EQ. 0) THEN
                WRITE(plotf(1),10) TPD(mapinc(iopt,iplot)), logKr
                WRITE(plotf(2),10) TPD(mapinc(iopt,iplot)), dGr*4.184d-3
                WRITE(plotf(7),10) TPD(mapinc(iopt,iplot)), 
     1                             TPD(mapv3(iopt,iplot)-isat) 
           END IF
           RETURN
      END IF

      IF ((xyplot .EQ. 1) .AND. xCp) THEN
           IF ((isat .EQ. 1) .OR. EQ3run) THEN
                WRITE(plotf(1),10) TPD(mapiso(iopt,iplot)), logKr
                WRITE(plotf(2),10) TPD(mapiso(iopt,iplot)), dGr*4.184d-3
                WRITE(plotf(3),10) TPD(mapiso(iopt,iplot)), dHr*4.184d-3
                WRITE(plotf(4),10) TPD(mapiso(iopt,iplot)), dSr*4.184
                WRITE(plotf(6),10) TPD(mapiso(iopt,iplot)), dVr*1d-1
                WRITE(plotf(7),10) TPD(mapiso(iopt,iplot)), 
     1                             TPD(mapinc(iopt,iplot)) 
                WRITE(plotf(8),10) TPD(mapiso(iopt,iplot)), 
     1                             TPD(mapv3(iopt,iplot)-isat) 
                RETURN
           END IF
           IF (univar .EQ. 0) THEN
                WRITE(plotf(1),10) TPD(mapinc(iopt,iplot)), logKr
                WRITE(plotf(2),10) TPD(mapinc(iopt,iplot)), dGr*4.184d-3
                WRITE(plotf(3),10) TPD(mapinc(iopt,iplot)), dHr*4.184d-3
                WRITE(plotf(4),10) TPD(mapinc(iopt,iplot)), dSr*4.184
                WRITE(plotf(6),10) TPD(mapinc(iopt,iplot)), dVr*1d-1
                WRITE(plotf(7),10) TPD(mapinc(iopt,iplot)), 
     1                             TPD(mapv3(iopt,iplot)-isat) 
           END IF
           RETURN
      END IF

      IF (xyplot .EQ. 1) THEN
***        .NOT. (xall .OR. xHSVCp .OR. xCp) 
           IF ((isat  .EQ. 1) .OR. EQ3run) THEN
                WRITE(plotf(1),10) TPD(mapiso(iopt,iplot)), logKr
                WRITE(plotf(2),10) TPD(mapiso(iopt,iplot)), dGr*4.184d-3
                WRITE(plotf(3),10) TPD(mapiso(iopt,iplot)), dHr*4.184d-3
                WRITE(plotf(4),10) TPD(mapiso(iopt,iplot)), dSr*4.184
                WRITE(plotf(5),10) TPD(mapiso(iopt,iplot)), dCpr*4.184
                WRITE(plotf(6),10) TPD(mapiso(iopt,iplot)), dVr*1d-1
                WRITE(plotf(7),10) TPD(mapiso(iopt,iplot)), 
     1                             TPD(mapinc(iopt,iplot)) 
                WRITE(plotf(8),10) TPD(mapiso(iopt,iplot)), 
     1                             TPD(mapv3(iopt,iplot)-isat) 
           END IF
           IF ((isat .EQ. 0) .AND. (.NOT. EQ3run)) THEN
                WRITE(plotf(1),10) TPD(mapinc(iopt,iplot)), logKr
                WRITE(plotf(2),10) TPD(mapinc(iopt,iplot)), dGr*4.184d-3
                WRITE(plotf(3),10) TPD(mapinc(iopt,iplot)), dHr*4.184d-3
                WRITE(plotf(4),10) TPD(mapinc(iopt,iplot)), dSr*4.184
                WRITE(plotf(5),10) TPD(mapinc(iopt,iplot)), dCpr*4.184
                WRITE(plotf(6),10) TPD(mapinc(iopt,iplot)), dVr*1d-1
                WRITE(plotf(7),10) TPD(mapinc(iopt,iplot)), 
     1                             TPD(mapv3(iopt,iplot)-isat) 
           END IF
	   RETURN
      END IF

***************************  xyplot = 2  *******************************

      IF (lvdome .OR. H2Oerr .OR. Klost .OR. xall) GO TO 20
     
      IF (univar .EQ. 1) THEN
           pbuff(inc,iso,1) = TPD(mapinc(iopt,iplot)) 
           GO TO 20 
      END IF

      IF (xHSVCp) THEN
           pbuff(inc,iso,1) = logKr
	   pbuff(inc,iso,2) = dGr*4.184d-3
	   pbuff(inc,iso,7) = TPD(mapv3(iopt,iplot)-isat) 
	   pbuff(inc,iso,8) = TPD(mapinc(iopt,iplot)) 
	   GO TO 20
      END IF

      IF (xCp) THEN
           pbuff(inc,iso,1) = logKr
	   pbuff(inc,iso,2) = dGr*4.184d-3
	   pbuff(inc,iso,3) = dHr*4.184d-3
	   pbuff(inc,iso,4) = dSr*4.184
	   pbuff(inc,iso,6) = dVr*1d-1
	   pbuff(inc,iso,7) = TPD(mapv3(iopt,iplot)-isat) 
	   pbuff(inc,iso,8) = TPD(mapinc(iopt,iplot)) 
      ELSE
           pbuff(inc,iso,1) = logKr
	   pbuff(inc,iso,2) = dGr*4.184d-3
	   pbuff(inc,iso,3) = dHr*4.184d-3
	   pbuff(inc,iso,4) = dSr*4.184
	   pbuff(inc,iso,5) = dCpr*4.184
	   pbuff(inc,iso,6) = dVr*1d-1
	   pbuff(inc,iso,7) = TPD(mapv3(iopt,iplot)-isat) 
	   pbuff(inc,iso,8) = TPD(mapinc(iopt,iplot)) 
      END IF

 20   IF ((iso .LT. niso) .OR. (inc .LT. nv2) .OR.
     1    (EQ3run .AND. (inc .LT. noninc))) RETURN

*** flush buffers
      IF (xyplot .EQ. 1) THEN
           DO 25 i = 1,nlogk
                WRITE(plotf(1),22) Kmin + (i-1)*Kinc, 
     1          isov(iopt,iplot), incv(iopt,iplot)
 22             FORMAT(/,' LOG K = ',e12.5,' ; ',a10,', ',a10,' = ',/)
                DO 25 j = 1,niso
                     IF (pbuff(i,j,1) .NE. 0.0d0) THEN
                          WRITE(plotf(1),10) isomin + (j-1)*isoinc, 
     1                    pbuff(i,j,1)
                     END IF
 25                  CONTINUE 
           RETURN
      END IF

      IF (EQ3run) THEN
	   DO 30 i = 1,noninc
                WRITE(plotf(1),40) oddv1(i), oddv2(i),
     1          pbuff(i,1,7), (pbuff(i,1,k), k = 1,6) 
 40             FORMAT(9(e14.7,2x))
 30             CONTINUE
           RETURN
      END IF

      IF (isat .EQ. 1) THEN
	   DO 50 i = 1,nv2
                IF ((i .EQ. 1) .AND. (v2min .EQ. 0.0d0)) THEN
                     WRITE(plotf(1),40) satmin(iopt),  
     1               pbuff(i,1,8), pbuff(i,1,7), 
     1               (pbuff(i,1,k), k = 1,6) 
                ELSE
                     WRITE(plotf(1),40) v2min+(i-1)*v2inc, 
     1               pbuff(i,1,8), pbuff(i,1,7), 
     1               (pbuff(i,1,k), k = 1,6) 
                END IF
 50             CONTINUE
           RETURN
      END IF

      IF (univar .EQ. 0) THEN
           DO 60 k = 1,nplots
                DO 60 i = 1,nv2  
                     WRITE(plotf(k),70) v2min+(i-1)*v2inc, 
     1               (pbuff(i,j,k), j = 1,MIN(niso,MXIPLT)) 
 70                  FORMAT(12(e14.7,2x))
 60                  CONTINUE
      ELSE
           DO 80 j = 1,niso
                WRITE(plotf(1),90) isomin+(j-1)*isoinc, 
     1          (pbuff(i,j,1), i = 1,MIN(nlogK,MXIPLT)) 
 90             FORMAT(12(e14.7,2x))
 80             CONTINUE
      END IF

      END

***********************************************************************

*** blanks - Set xall, xCp, and MKwarn per current state conditions.

      SUBROUTINE blanks(T,P,D,m2reac,nm,ng,na,xall,xHSVCp,xCp,
     1                  MKwarn,Pwarn)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

**JP**      PARAMETER (MAXMIN = 10, MXTRAN = 3, IABC = 3,
      PARAMETER (MAXMIN = 10, MXTRAN = 3, IABC = 4,
     1           MAXGAS = 10, MAXAQS = 10,
     1           TCPM2  =   25.0d0,  PMAXMG = 50000.0d0,
     2           TMAXA  = 1000.0d0,  PMAXA  =  5000.0d0,  
     3           TMAXX  =  400.0d0,  PMAXX1 =   500.0d0,  
     4           TMINX  =  350.0d0,  PMAXX2 =  1000.0d0,  
     5           DMINCA =  0.35d0, DMINNA = 0.05d0, TOL = 1.0d-5)

      CHARACTER*20  mname(MAXMIN), aname(MAXAQS), gname(MAXGAS)
      CHARACTER*30  mform(MAXMIN), aform(MAXAQS), gform(MAXGAS)

      LOGICAL  m2reac, xall, xHSVCp, xCp, MKwarn, Pwarn, aqschg

      INTEGER ntran(MAXMIN), phaser(MAXMIN)

      DOUBLE PRECISION  TtranP(MXTRAN,MAXMIN), 
     1                  PtranT(MXTRAN,MAXMIN)

      DOUBLE PRECISION  Vmin(MAXMIN), Smin(MAXMIN), Cpmin(MAXMIN),
     2                  Hmin(MAXMIN), Gmin(MAXMIN)

      DOUBLE PRECISION Gfaqs(MAXAQS), Hfaqs(MAXAQS), SPrTra(MAXAQS), 
     1                 a(4,MAXAQS), c(2,MAXAQS), 
     2                 wref(MAXAQS), chg(MAXAQS)

      DOUBLE PRECISION  Gfmin(MAXMIN), Hfmin(MAXMIN), VPrTrm(MAXMIN),
     1                  SPrTrm(MAXMIN), Alpha(MAXMIN), Kappa(MAXMIN),
     2                  KappaP(MAXMIN), KappaDP(MAXMIN),
     3                  MK1(IABC,MAXMIN), MK2(IABC,MAXMIN),
     4                  MK3(IABC,MAXMIN), MK4(IABC,MAXMIN),
     5                  Ttran(MXTRAN,MAXMIN), Htran(MXTRAN,MAXMIN),
     6                  Vtran(MXTRAN,MAXMIN), dPdTtr(MXTRAN,MAXMIN),
     7                  Tmaxm(MAXMIN), NumAtoms(MAXMIN), Tc(MAXMIN),
     8                  Smax(MAXMIN), Vmax(MAXMIN)

      DOUBLE PRECISION Gfgas(MAXGAS), Hfgas(MAXGAS), VPrTrg(MAXGAS), 
     1                 SPrTrg(MAXGAS), MKg(IABC,MAXGAS), Tmaxg(MAXGAS)

      COMMON /PTtran/ TtranP, PtranT
      COMMON /minsp/  Vmin, Smin, Cpmin, Hmin, Gmin, phaser

      COMMON /anames/ aname, aform
      COMMON /aqsref/ Gfaqs, Hfaqs, SPrTra, c, a, wref, chg

      COMMON /mnames/ mname, mform 
      COMMON /minref/ Gfmin, Hfmin, SPrTrm, VPrTrm, MK1, MK2, MK3, MK4, 
     1                Ttran, Htran, Vtran, dPdTtr, Tmaxm, ntran, Alpha,
     2                Kappa, KappaP, KappaDP, NumAtoms, Tc, Smax, Vmax

      COMMON /gnames/ gname, gform
      COMMON /gasref/ Gfgas, Hfgas, SPrTrg, VPrTrg, MKg, Tmaxg

      SAVE

      xall   = .FALSE.
      xHSVCp = .FALSE.
      xCp    = .FALSE.
      MKwarn = .FALSE.
      Pwarn  = .FALSE.
      aqschg = .FALSE.

      IF ((nm .GT. 0) .OR. (ng .GT. 0)) THEN
***** consider validity limitations of mineral/gas equations

           IF (P .GT. PMAXMG+TOL) THEN
*****      P exceeds pressure limit; issue warning
                Pwarn = .TRUE.
           END IF

           TK = T + 273.15d0 
           DO 10 i = 1,ng
                IF (TK .GT. Tmaxg(i)) THEN
*****           T beyond limit of Maier-Kelly coefficients;
*****           issue appropriate warning in calling routine
                     MKwarn = .TRUE.
                END IF
 10             CONTINUE

           DO 11 i = 1,nm
                IF (TK .GT. Tmaxm(i)) THEN
*****           T beyond limit of Maier-Kelly coefficients;
*****           issue appropriate warning in calling routine
                     MKwarn = .TRUE.
                END IF
 11             CONTINUE

      END IF

 25   IF (na .GT. 0) THEN
***** consider validity limitations of aqueous species equations
 
           IF ((P .GT. PMAXA+TOL) .OR. (T .GT. TMAXA+TOL)) THEN
*****      P or T exceeds validity limits of equations
                xall = .TRUE.
                RETURN
           END IF

           DO 30 j = 1,na
                IF (chg(j) .NE. 0.0d0) THEN
                     aqschg = .TRUE.
                END IF
 30             CONTINUE

           IF (aqschg) THEN
*****      check limits for charged aqueous species     
                xall = ((D .LT. DMINCA) .OR. ((P. LT. PMAXX1) .AND.
     1                  (T .GT. TMINX) .AND. (T .LT. TMAXX)))
                xHSVCp = xall .OR.
     1                   ((P .LT. PMAXX2) .AND. (T .GT. TMINX))
           ELSE
*****      check limits for neutral aqueous species     
                xall = (D .LT. DMINNA)
           END IF

      END IF

      END
