      SUBROUTINE FORCEHEFLASH(TARGETCOREMASS, MODELNMESH)
C This subroutine should force the code around the helium flash by
C taking a 3 solar mass star that has just started non-degenerate 
C core He, and then adjusting it so it "matches" the "expected" 
C post-Flash core He burning model
      IMPLICIT REAL*8 (A-H,N-Z)
      PARAMETER (MAXMSH = 2000)
      REAL*8 RSTORE(100), HSTORE(60, MAXMSH), HNUCSTORE(100,MAXMSH)
      INTEGER ISTORE(100)
      CHARACTER*200 FILENAMESTORE(20)
      CHARACTER*43 M3MODELDIR
      
      SAVE
C Required Common Blocks***********************************************
      
      COMMON /NUCMAT/ HNUC(100,MAXMSH), DHNUC(100,MAXMSH)
      COMMON /PREVM / HPR(60,MAXMSH),DHPR(60,MAXMSH),PR(14),MS(9999),
     &                ST(9999),NPR,KPR, HPPR(60,MAXMSH),
     &                HNUCPR(100,MAXMSH), DHNUCPR(100, MAXMSH)
      COMMON H(60,MAXMSH), DH(60,MAXMSH), EP(3), NH, JIN, ID(200)
      COMMON /AUXIN / ICL, ION, IAM, IOP, INUC, IBC, ICN, IML(2), ISGTH,
     :     IMO, IDIFF
      COMMON /SODDS / ALPHA, RML, CMG, CSI, CFE, CT(10), AGE, DT, M1, 
     :  EC, BM, ANG, CM, MTA, MTB, TM(2), T0, M0, TC(2), OS, AC, RCD,
     :  RMG, RHL, XF, DR, AK1 ,RMT, AK2, IZ(4), IB, ISX(45),
     :  TRB
      COMMON /STAT1 / CSX(10), CS(90, 127, 10), HAT(23320), NCSX
      COMMON /OP    / ZS, LEDD, VM, GR, GRAD, ETH, RLF, EGR, R, QQ
      COMMON /ATDATA/ DH2(4), CHI(26,9), OMG(27), AM(10), BN(10), JZ(10)
      COMMON /NDATA / RATEN(9000)
      COMMON /CNSTS / CPI, PI4, CLN10, CDUM(11), CSECYR, LSUN, MSUN,
     &                RSUN, TSUNYR
      COMMON /YUK1  / PX(34), WMH, WMHE, VMH, VME, VMC, VMG, BE, VLH, 
     :                VLE, VLC, VLN, VLT, MCB(12),WWW(100)
      COMMON /OPDAT / cbase,obase,opT(141),opR(31),fZ
      COMMON /XOPDAT/ opac(4,4,141,31,5)
      COMMON /COPDAT/ opacCO(4,4,141,31,305)
      COMMON /EVMODE/ IMODE
      COMMON /MIXFUD/ SGTHFAC, FACSGMIN, FACSG, ISGFAC
* extra common for mesh-spacing
      COMMON /PMESH / PMH(2), PME(2), IAGB
* first guess of pressure at H, He-burning shell (should be in input file!)
*      data pmh, pme /1.0e17, 7.5e19/
*
* Common for artifical energy injection
*      REAL TENG, SMASS, FMASS
      COMMON /ARTENG/ TENG, SMASS, FMASS, SHIENG, INJMD
      COMMON /INJTIMECTRL/ STARTTIMEINJ, ENDTIMEINJ
      COMMON /AGECTRL/ ENDAGE

* Extra COMMON for main-sequence evolution.
*     
      
      COMMON /ZAMS  / TKH(2)
      COMMON /MESH  / TRC1,TRC2,DD,DT1,DT2,MWT,MWTS,IVMC,IVMS
      COMMON /DHBLOC/ IDREDGE
      COMMON /DTCONT/ VLHP(2), VLEP(2), VLCP(2), RLFP(2), TOTMP(2), VLHC(2),
     :     VLEC(2), VLCC(2), RLFC(2), TOTMC(2)
      COMMON /ANGMOM/ VROT1, VROT2, FMAC, FAM, IRAM, IRS1, IRS2
      COMMON /DIFCOE/ DC(50,4,3), DCD(50,4)
      COMMON /FLASHSTORE/ SM, NM
C******COMMON BLOCKS FINISHED ****************************************
C***************FORMAT STUFF******************************************
99003 FORMAT (12I4,/,12I4,/,7I4,/,1P,5E8.1,0P,/,2(10I3,/,3(30I3,/)),3(15I
     :3,/), 9F5.2, 1P, 3E8.1,
     :/, E9.2, 0P, 9F6.3, /, 1P, 2(7E9.2, /), 0P, I2, 2(I2,1X,E8.2),2(1X,F4.2)
     : ,/, I2,F6.1,I2,F6.1, 1X, F4.2, I2, I2, 2(1X, E8.2)
     : ,/, 3E14.6, I2, 2E14.6, /, E14.6)
C*************FORMAT STUFF DONE***************************************
      
C Setup extra vars for the he flash
      M3MODELMASS = 3.0
      TARGETDTY = 1.5E+05
      MASSGAIN = -1E-6
      COREGROWDTY = 1E+06
      
C Parameters of the target model
      TARGETMASS = DEXP(H(4, 1))
      TIMENEEDED = (TARGETMASS - M3MODELMASS*MSUN) / MASSGAIN
      NTTIMESTEPS = 1 + INT(TIMENEEDED/TARGETDTY)
      IF (NTTIMESTEPS.LT.10) NTIMESTEPS = 10
      MDTY = TIMENEEDED / DBLE(NTIMESTEPS)
      
C Store a whole load of data variables
C BY FAR the easiest way to get this to work is just to dump
C the data into a temp-ish file, it'll be fine sure
C Create the temp file
      OPEN(UNIT=111,FILE='data_store_for_flash',ACCESS='sequential')
      WRITE (111,99003) NH2,IT1,IT2,JIN,JOUT,NCH,JP,IZ,IMODE,
     :ICL,ION,IAM,IOP,IBC,INUC,ICN,IML(1),IML(2),ISGTH,IMO,IDIFF,
     :NT1,NT2,NT3,NT4,NT5,NSV,NMONT,
     :EP,DT3,DD,ID,ISX,DT1,DT2,CT,ZS,ALPHA,CH,CC,CN,CO, 
     :CNE,CMG,CSI,CFE,RCD,OS,RML,RMG,ECA,XF,DR,RMT,RHL,AC,AK1,AK2,ECT,
     :TRB,
     :IRAM, IRS1, VROT1, IRS2, VROT2, FMAC, FAM,
     :IVMC, TRC1, IVMS, TRC2, MWTS, IAGB, ISGFAC, FACSGMIN, SGTHFAC,
     :TENG, SMASS, FMASS, INJMD, STARTTIMEINJ, ENDTIMEINJ,
     :ENDAGE
C Store the model
      DO j = 1,MAXMSH
        DO i = 1, 60
            HSTORE(i,j) = H(i,j)
        ENDDO
      ENDDO
      
C Store the nuke data
      DO j = 1, MAXMSH
         DO i = 1, 100
            HNUCSTORE(i,j) = HNUC(i,j)
         ENDDO
      ENDDO
      DT = CSECYR * MDTY
      
C Now, need to read in the 3 solar mass for masslos, IEND = -2
      CALL PRINTA ( -2, NSTEP, ITER1, ITER2, NWRT5 )
*      SM, DTY, AGE, PER, BMS, EC,NH,NP,NMOD,IB,PMH(1),PME(1)
      
C Now, we have a standard evolution loop
      ITER = 100
      NM = 0
      NWRT5 = 0
      
      DO WHILE (NTTIMESTEPS.GT.0)
         CALL SOLVER ( 1, ITER, KTER, ERR, ID, NWRT5 )
         CALL PRINTA ( 0, NSTEP, ITER1, ITER2, NWRT5 )
         NTIMESTEPS = NTIMESTEPS - NM
         IF (ERR.GT.EPS) WRITE (*,*) '!!Error Exceeded in massloss!! FIX THE DAMN CODE'
C TODO Now, need to retune the timestep, in case we failed to converge TODO
        
            
      ENDDO
      
C Now, load in the core burning stuff
      DT = CSECYR * COREGROWDTY
      
C Blank out the DH and DHNUC stuff to cancel out mass loss
      DO j=1, MAXMSH
        DO i=1,60
            DH(i,j) = 0.D0
        ENDDO
      ENDDO
      DO j=1,MAXMSH
        DO i=1,100
            DHNUC(i,j) = 0.D0
        ENDDO
      ENDDO
      
C Load in the core burning model 
      CALL PRINTA ( -3, NSTEP, ITER1, ITER2, NWRT5 )
      COREMASS = VMH
      DO WHILE ((COREMASS.LT.TARGETCOREMASS)
     &         .AND.(((TARGETCOREMASS - COREMASS)/TARGETCOREMASS).GT.1.D-3))
         
         PREVCOREMASS = COREMASS
         CALL SOLVER ( 1, ITER, KTER, ERR, ID, NWRT5 )
         CALL PRINTA ( 0, NSTEP, ITER1, ITER2, NWRT5 )
         IF (ERR.GT.EPS) WRITE (*,*) '!!Error Exceeded in core burn!! FIX THE DAMN CODE'
         COREMASS = VMH
         TIMENEEDED = DT * (TARGETCOREMASS - COREMASS) / 
     &      (COREMASS - PREVCOREMASS)
         DT = TIMENEEDED / 
     &        DMAX1(DNINT(TIMENEEDED / (CSECYR*COREGROWDTY)),1.D0)
         
      ENDDO
      
C Okay, with all this done, lets restore the envelope compositions
      i=1
      j=1
      DO WHILE ((i.LT.MAXSH).AND.((DEXP(H(4,i)) / MSUN).GT.(COREMASS + 1.D-2)))
C As Ross did, lets perform a linear interpolation in log mass space
         DO WHILE (H(4,i).LT.HSTORE(4,j+1))
            j = j+1
         ENDDO
C This bit copied wholesale from Ross, hope it works!
         frac = (H(4,i)-Hstore(4,j))/(Hstore(4,j+1)-Hstore(4,j))
         H(3,i) = Hstore(3,j)*(1-frac) + Hstore(3,j+1)*frac
         H(5,i) = Hstore(5,j)*(1-frac) + Hstore(5,j+1)*frac
         H(9,i) = Hstore(9,j)*(1-frac) + Hstore(9,j+1)*frac
         H(10,i) = Hstore(10,j)*(1-frac) + Hstore(10,j+1)*frac
         H(11,i) = Hstore(11,j)*(1-frac) + Hstore(11,j+1)*frac
         i=i+1
      ENDDO
      
C Okay okay, restore everything from that datafile 
      READ (111,99003) NH2,IT1,IT2,JIN,JOUT,NCH,JP,IZ,IMODE,
     :ICL,ION,IAM,IOP,IBC,INUC,ICN,IML(1),IML(2),ISGTH,IMO,IDIFF,
     :NT1,NT2,NT3,NT4,NT5,NSV,NMONT,
     :EP,DT3,DD,ID,ISX,DT1,DT2,CT,ZS,ALPHA,CH,CC,CN,CO, 
     :CNE,CMG,CSI,CFE,RCD,OS,RML,RMG,ECA,XF,DR,RMT,RHL,AC,AK1,AK2,ECT,
     :TRB,
     :IRAM, IRS1, VROT1, IRS2, VROT2, FMAC, FAM,
     :IVMC, TRC1, IVMS, TRC2, MWTS, IAGB, ISGFAC, FACSGMIN, SGTHFAC,
     :TENG, SMASS, FMASS, INJMD, STARTTIMEINJ, ENDTIMEINJ,
     :ENDAGE
C and nuke the temp file
      CLOSE(111,STATUS='DELETE')
      
            
      
      
      RETURN
      END
      
      

      
      
      
            
      
      
      
      
      
      
      
      
      
       