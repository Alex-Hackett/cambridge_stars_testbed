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
      REAL*8 IZ_STORE(4), EP_STORE(3), ID_STORE(200), ISX_STORE(45),
     :       CT_STORE(10) 
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
*      OPEN(UNIT=111,FILE='data_store_for_flash',ACCESS='sequential')
*      WRITE (111,*) NH2,IT1,IT2,JIN,JOUT,NCH,JP,IZ,IMODE,
*     :ICL,ION,IAM,IOP,IBC,INUC,ICN,IML(1),IML(2),ISGTH,IMO,IDIFF,
*     :NT1,NT2,NT3,NT4,NT5,NSV,NMONT,
*     :EP,DT3,DD,ID,ISX,DT1,DT2,CT,ZS,ALPHA,CH,CC,CN,CO, 
*     :CNE,CMG,CSI,CFE,RCD,OS,RML,RMG,ECA,XF,DR,RMT,RHL,AC,AK1,AK2,ECT,
*     :TRB,
*     :IRAM, IRS1, VROT1, IRS2, VROT2, FMAC, FAM,
*     :IVMC, TRC1, IVMS, TRC2, MWTS, IAGB, ISGFAC, FACSGMIN, SGTHFAC,
*     :TENG, SMASS, FMASS, INJMD, STARTTIMEINJ, ENDTIMEINJ,
*     :ENDAGE
      NH2_STORE = NH2
      IT1_STORE = IT1
      IT2_STORE = IT2
      JIN_STORE = JIN
      JOUT_STORE = JOUT
      NCH_STORE = NCH
      JP_STORE = JP
      IZ_STORE = IZ
      IMODE_STORE = IMODE
      ICL_STORE = ICL
      ION_STORE = ION
      IAM_STORE = IAM
      IOP_STORE = IOP
      IBC_STORE = IBC
      INUC_STORE = INUC
      ICN_STORE = ICN
      IML_STORE1 = IML(1)
      IML_STORE2 = IML(2)
      ISGTH_STORE = ISGTH
      IMO_STORE = IMO
      IDIFF_STORE = IDIFF
      NT1_STORE = NT1
      NT2_STORE = NT2
      NT3_STORE = NT3
      NT4_STORE = NT4
      NT5_STORE = NT5
      NSV_STORE = NSV
      NMONT_STORE = NMONT
      EP_STORE = EP
      DT3_STORE = DT3
      DD_STORE = DD
      ID_STORE = ID
      ISX_STORE = ISX
      DT1_STORE = DT1
      DT2_STORE = DT2
      CT_STORE = CT
      ZS_STORE = ZS
      ALPHA_STORE = ALPHA
      CH_STORE = CH
      CC_STORE = CC
      CN_STORE = CN
      CO_STORE = CO
      CNE_STORE = CNE
      CMG_STORE = CMG
      CSI_STORE = CSI
      CFE_STORE = CFE
      RCD_STORE = RCD
      OS _STORE = OS
      RML_STORE = RML
      RMG_STORE = RMG
      ECA_STORE = ECA
      XF_STORE = XF
      DR_STORE = DR
      RMT_STORE = RMT
      RHL_STORE = RHL
      AC_STORE = AC
      AK1_STORE = AK1
      AK2_STORE = AK2
      ECT_STORE = ECT
      TRB_STORE = TRB
      IRAM_STORE = IRAM
      IRS1_STORE = IRS1
      VROT1_STORE= VROT1
      IRS2_STORE = IRS2
      VROT2_STORE = VROT2
      FMAC_STORE = FMAC
      FAM_STORE = FAM
      IVMC_STORE = IVMC
      TRC1_STORE = TRC1
      IVMS_STORE = IVMS
      TRC2_STORE = TRC2
      MWTS_STORE = MWTS
      IAGB_STORE = IAGB
      ISGFAC_STORE = ISGFAC
      FACSGMIN_STORE = FACSGIM
      SGTHFAC_STORE = SGTHFAC
      TENG_STORE = TENG
      SMASS_STORE = SMASS
      FMASS_STORE = FMASS
      INJMD_STORE = INJMD
      STARTTIMEINJ_STORE = STARTTIMEINJ
      ENDTIMEINJ_STORE = ENDTIMEINJ
      ENDAGE_STORE = ENDAGE


C Store the model
      DO j = 1,MAXMSH
        DO i = 1, 60
            HSTORE(i,j) = H(i,j)
        ENDDO
      ENDDO
      
C Store the nuke data
      DO j = 1, MAXMSH
         DO i = 1, 50
            HNUCSTORE(i,j) = HNUC(i,j)
         ENDDO
      ENDDO
      DT = CSECYR * MDTY

C Now, need to read in the 3 solar mass for masslos, IEND = -2
      CALL PRINTA ( -2, NSTEP, ITER1, ITER2, NWRT5 )

*      SM, DTY, AGE, PER, BMS, EC,NH,NP,NMOD,IB,PMH(1),PME(1)
      
C Now, we have a standard evolution loop
      ITER = 100
      NMNEW = 0
*      NWRT5 = 0
      WRITE (*,*) NTIMESTEPS
      DO WHILE (NTTIMESTEPS.GT.0)
         WRITE (*,*) 'WE"VE REACHED THE MASS LOSS'
         CALL SOLVER ( 1, ITER, KTER, ERR, ID, NWRT5 )
         CALL PRINTA ( 0, NSTEP, ITER1, ITER2, NWRT5 )
         NMNEW = NMNEW + 1
         ITER = ITER + 1
         NTIMESTEPS = NTIMESTEPS - NMNEW
         IF (ERR.NE.0) WRITE (*,*) '!!Error Exceeded in massloss!! FIX THE DAMN CODE'
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
*      DO j=1,MAXMSH
*        DO i=1,100
*            DHNUC(i,j) = 0.D0
*        ENDDO
*      ENDDO
      
C Load in the core burning model 
      CALL PRINTA ( -3, NSTEP, ITER1, ITER2, NWRT5 )
      COREMASS = VMH
      DO WHILE ((COREMASS.LT.TARGETCOREMASS)
     &         .AND.(((TARGETCOREMASS - COREMASS)/TARGETCOREMASS).GT.1.D-3))
         WRITE (*,*) 'WE"VE REACHED THE CORE BURN'
         PREVCOREMASS = COREMASS
         CALL SOLVER ( 1, ITER, KTER, ERR, ID, NWRT5 )
         CALL PRINTA ( 0, NSTEP, ITER1, ITER2, NWRT5 )
         IF (ERR.NE.0) WRITE (*,*) '!!Error Exceeded in core burn!! FIX THE DAMN CODE'
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
*      READ (111,*) NH2,IT1,IT2,JIN,JOUT,NCH,JP,IZ,IMODE,
*     :ICL,ION,IAM,IOP,IBC,INUC,ICN,IML(1),IML(2),ISGTH,IMO,IDIFF,
*     :NT1,NT2,NT3,NT4,NT5,NSV,NMONT,
*     :EP,DT3,DD,ID,ISX,DT1,DT2,CT,ZS,ALPHA,CH,CC,CN,CO, 
*     :CNE,CMG,CSI,CFE,RCD,OS,RML,RMG,ECA,XF,DR,RMT,RHL,AC,AK1,AK2,ECT,
*     :TRB,
*     :IRAM, IRS1, VROT1, IRS2, VROT2, FMAC, FAM,
*     :IVMC, TRC1, IVMS, TRC2, MWTS, IAGB, ISGFAC, FACSGMIN, SGTHFAC,
*     :TENG, SMASS, FMASS, INJMD, STARTTIMEINJ, ENDTIMEINJ,
*     :ENDAGE
C and nuke the temp file
*      CLOSE(111,STATUS='DELETE')

     

C Idiot proofing -- otherwise the logic in solver will fail
*      FACSGMIN = DMIN1(1d0, FACSGMIN)
      
      NH2 = NH2_STORE
      IT1 = IT1_STORE
      IT2 = IT2_STORE
      JIN = JIN_STORE
      JOUT = JOUT_STORE
      NCH = NCH_STORE
      JP = JP_STORE
      IZ = IZ_STORE
      IMODE = IMODE_STORE
      ICL = ICL_STORE
      ION = ION_STORE
      IAM = IAM_STORE
      IOP = IOP_STORE
      IBC = IBC_STORE
      INUC = INUC_STORE
      ICN = ICN_STORE
      IML(1) = IML_STORE1
      IML(2) = IML_STORE2
      ISGTH = ISGTH_STORE
      IMO = IMO_STORE
      IDIFF = IDIFF_STORE
      NT1 = NT1_STORE
      NT2 = NT2_STORE
      NT3 = NT3_STORE
      NT4 = NT4_STORE
      NT5 = NT5_STORE
      NSV = NSV_STORE
      NMONT = NMONT_STORE
      EP = EP_STORE
      DT3 = DT3_STORE
      DD = DD_STORE
      ID = ID_STORE
      ISX = ISX_STORE
      DT1 = DT1_STORE
      DT2 = DT2_STORE
      CT = CT_STORE
      ZS = ZS_STORE
      ALPHA = ALPHA_STORE
      CH = CH_STORE
      CC = CC_STORE
      CN = CN_STORE
      CO = CO_STORE
      CNE = CNE_STORE
      CMG = CMG_STORE
      CSI = CSI_STORE
      CFE = CFE_STORE
      RCD = RCD_STORE
      OS = OS_STORE
      RML = RML_STORE
      RMG = RMG_STORE
      ECA = ECA_STORE
      XF = XF_STORE
      DR = DR_STORE
      RMT = RMT_STORE
      RHL = RHL_STORE
      AC = AC_STORE
      AK1 = AK1_STORE
      AK2 = AK2_STORE
      ECT = ECT_STORE
      TRB = TRB_STORE
      IRAM = IRAM_STORE
      IRS1 = IRS1_STORE
      VROT1= VROT1_STORE
      IRS2 = IRS2_STORE
      VROT2 = VROT2_STORE
      FMAC = FMAC_STORE
      FAM = FAM_STORE
      IVMC = IVMC_STORE
      TRC1 = TRC1_STORE
      IVMS = IVMS_STORE
      TRC2 = TRC2_STORE
      MWTS = MWTS_STORE
      IAGB = IAGB_STORE
      ISGFAC = ISGFAC_STORE
      FACSGMIN = FACSGIM_STORE
      SGTHFAC = SGTHFAC_STORE
      TENG = TENG_STORE
      SMASS = SMASS_STORE
      FMASS = FMASS_STORE
      INJMD = INJMD_STORE
      STARTTIMEINJ = STARTTIMEINJ_STORE
      ENDTIMEINJ = ENDTIMEINJ_STORE
      ENDAGE = ENDAGE_STORE
      
C Restore the nuke data
      DO j = 1, MAXMSH
         DO i = 1, 50
            HNUC(i,j) = HNUCSTORE(i,j)
         ENDDO
      ENDDO      
            
      
      
      RETURN
      END
      
      

      
      
      
            
      
      
      
      
      
      
      
      
      
       