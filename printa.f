**==PRINTA.FOR
      SUBROUTINE PRINTA(IEND, NP, IT1, IT2, NT5)
      IMPLICIT REAL*8(A-H, L, M, O-Z)
      real*8 MAT(4,141),Xcompos(3,305),COcompos(8)
      SAVE
      INTEGER MAXMSH
          
      PARAMETER (MAXMSH = 2000)
      CHARACTER*128 FILENAMESTORE(20)
      CHARACTER*128 COtableNAME
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
      COMMON /STAT1 / CSX(10), CS(121,191,10), HAT(23320), NCSX
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
      COMMON /MASSCTRL/ MINMASS, MAXMASS
      SAVE
* Extra COMMON for main-sequence evolution.
*     
      COMMON /FLASHSTORE/ SM, NM
      COMMON /ZAMS  / TKH(2)
      COMMON /MESH  / TRC1,TRC2,DD,DT1,DT2,MWT,MWTS,IVMC,IVMS
      COMMON /DHBLOC/ IDREDGE
      COMMON /DTCONT/ VLHP(2), VLEP(2), VLCP(2), RLFP(2), TOTMP(2), VLHC(2),
     :     VLEC(2), VLCC(2), RLFC(2), TOTMC(2)
      COMMON /ANGMOM/ VROT1, VROT2, FMAC, FAM, IRAM, IRS1, IRS2
      COMMON /DIFCOE/ DC(50,4,3), DCD(50,4)
      CBRT(VX) = DEXP(DLOG(VX)/3.0D0)
      RLOBE(VX) = 0.49D0*VX*VX/(0.6D0*VX*VX+DLOG(1.0D0+VX)) 
      DIMENSION WW(16),WX(52),DELDAT(22)
      data COcompos /0.0d0,0.01d0,0.03d0,0.1d0,0.2d0,0.4d0,0.6d0,1.0d0/
      
C Common blocks for TZO stuff
      COMMON /ITZO/ itzo_yn,itzo_cmass_pre, itzo_stripcorehe, itzo_stophighburn,
     :          itzo_noneutburn, itzo_zerocore,
     :          itzo_ct_1, itzo_ct_2, itzo_ct_3
      COMMON /RTZO/ rtzo_mod_emass, rtzo_dcmassdt, rtzo_maxdt,
     :          rtzo_cut_non_degen_hburn, rtzo_EC, rtzo_nucap,
     :          rtzo_nucap_per_yr, rtzo_nucap_min, rtzo_degen_cutoff,
     :          rtzo_degen_cutoff_per_yr, rtzo_degen_cutoff_max,
     :          rtzo_RCD_per_yr, rtzo_RCD_max,
     :          rtzo_meshfluid, rtzo_meshfluid_per_yr, rtzo_meshfluid_min,
     :          rtzo_alpha, rtzo_alpha_per_yr, rtzo_alpha_max,
     :          rtzo_ct_1, rtzo_ct_1_per_yr, rtzo_ct_1_max,
     :          rtzo_ct_2, rtzo_ct_2_per_yr, rtzo_ct_2_max,
     :          rtzo_ct_3, rtzo_ct_3_per_yr, rtzo_ct_3_max
      COMMON /TZOSTUFF/ cmass
      COMMON /PRETZO/ S_DTZO(100)
      
*99002 FORMAT (1X, 1P, 12E13.5, 0P) 33

99002 FORMAT (1P, 50E15.8, 0P)
99003 FORMAT (12I4,/,12I4,/,7I4,/,1P,5E8.1,0P,/,2(10I3,/,3(30I3,/)),34I
     :3,/,11I3,/, 9F5.2, 1P, 3E8.1,
     :/, E9.2, 0P, 9F6.3, /, 1P, 2(7E9.2, /), 0P, I2, 2(I2,1X,E8.2),2(1X,F4.2)
     : ,/, I2,F6.1,I2,F6.1, 1X, F4.2, I2, I2, 2(1X, E8.2)
     : ,/, 3E14.6, I2, 2E14.6, /, E14.6, /, E14.6, E14.6)
99004 FORMAT (1X, 10F7.3)
99005 FORMAT (1X, 1P, 2E14.6, E17.9, 3E14.6, 0P, 4I6, 1P, 2E11.3)
C Are we loading in a Helium flash model?
!      IF ((IEND.EQ.-2).OR.(IEND.EQ.-3)) GOTO 42
C Or just loading in our initial model?
      IF ( IEND.NE.-1 ) GOTO 30
* Initialize physical constants
      CALL CONSTS
C Read opacity, nuclear reaction and neutrino loss rate data
      READ (11,'(I4)') NCSX
      READ (11,99004) CSX
      DO N = 1, NCSX
         READ (11,99004) ((CS(JR,JT,N),JR=1,90),JT=1,127)
      ENDDO
      READ (11,99004) HAT
C Old table read, file pointer in right place, now overwrite with
C Ross's table
      OPEN (UNIT=98, FILE='dat/phys02-v4.dat', ACCESS='SEQUENTIAL', STATUS='OLD')
      READ (98, '(I4)') NCSX
      READ (98, 99004) CSX
      DO N = 1, NCSX
        READ (98, 99004) ((CS(JR,JT,N),JR=1,121),JT=1,191)
      ENDDO
      CLOSE(98)
      READ (13,99004) RATEN
C RJS 18/4/08 - read in spline coefficients for diffusion
99006 FORMAT (E12.5,3(1X,E12.5))
      DO K = 1,3
         DO I = 1,50
            READ (14,99006) (DC(I,J,K),J=1,4)
         END DO
      END DO
C d-coefficients
      DO I = 1,50
         READ (14,99006) (DCD(I,J),J=1,4)
      END DO
      DO 2 J = 1, 60
C      CT(J) = 0.0D0
         DO 2 K = 1,MAXMSH
            H(J, K) = 0.0D0
    2       DH(J, K) = 0.0D0
      DO 21 J=1,9999
         MS(J) = 0.0D0
   21    ST(J) = J
C Read miscellaneous data, usually unchanged during one evol run         
      READ  (1,99003) NH2,IT1,IT2,JIN,JOUT,NCH,JP,IZ,IMODE,
     :ICL,ION,IAM,IOP,IBC,INUC,ICN,IML(1),IML(2),ISGTH,IMO,IDIFF,
     :NT1,NT2,NT3,NT4,NT5,NSV,NMONT,
     :EP,DT3,DD,ID,ISX,DT1,DT2,CT,ZS,ALPHA,CH,CC,CN,CO, 
     :CNE,CMG,CSI,CFE,RCD,OS,RML,RMG,ECA,XF,DR,RMT,RHL,AC,AK1,AK2,ECT,
     :TRB,
     :IRAM, IRS1, VROT1, IRS2, VROT2, FMAC, FAM,
     :IVMC, TRC1, IVMS, TRC2, MWTS, IAGB, ISGFAC, FACSGMIN, SGTHFAC,
     :TENG, SMASS, FMASS, INJMD, STARTTIMEINJ, ENDTIMEINJ,
     :ENDAGE,
     :MINMASS, MAXMASS
     
C Read the TZO Control data
      !OPEN (UNIT=70, FILE='tzo_data', ACCESS='SEQUENTIAL', STATUS='OLD')
      READ (70, 99101) itzo_yn,itzo_cmass_pre, itzo_stripcorehe, itzo_stophighburn,
     :          itzo_noneutburn, itzo_zerocore,
     :          itzo_ct_1, itzo_ct_2, itzo_ct_3
      READ (70, 99102) rtzo_mod_emass, rtzo_dcmassdt,
     :          rtzo_maxdt,
     :          rtzo_cut_non_degen_hburn,
     :          rtzo_EC,
     :          rtzo_nucap, rtzo_nucap_per_yr, rtzo_nucap_min,
     :          rtzo_degen_cutoff,rtzo_degen_cutoff_per_yr, rtzo_degen_cutoff_max,
     :          rtzo_RCD_per_yr, rtzo_RCD_max,
     :          rtzo_meshfluid, rtzo_meshfluid_per_yr, rtzo_meshfluid_min,
     :          rtzo_alpha, rtzo_alpha_per_yr, rtzo_alpha_max,
     :          rtzo_ct_1, rtzo_ct_1_per_yr, rtzo_ct_1_max,
     :          rtzo_ct_2, rtzo_ct_2_per_yr, rtzo_ct_2_max,
     :          rtzo_ct_3, rtzo_ct_3_per_yr, rtzo_ct_3_max
      !CLOSE (70)
99101 FORMAT (I4,/,I4,/,I4,/,I4,/,I4,/,I4,/,3I4)
99102 FORMAT (2E14.6,/,E14.6,/,E14.6,/,E14.6,/,3E14.6,/,
     :        3E14.6,/,2E14.6,/,3E14.6,/,3E14.6,/,
     :        3E14.6,/,3E14.6,/,3E14.6)  
C Idiot proofing -- otherwise the logic in solver will fail
      FACSGMIN = DMIN1(1d0, FACSGMIN)
C Read data for initial model (often last model of previous run)
C e.g. SM = stellar mass, solar units; DTY = next timestep, years	 
      READ  (30, 99005) SM, DTY, AGE, PER, BMS, EC,NH,NP,NMOD,IB,PMH(1),PME(1)
      WRITE (333,99003) NH2,IT1,IT2,JIN,JOUT,NCH,JP,IZ,IMODE,
     :ICL,ION,IAM,IOP,IBC,INUC,ICN,IML(1),IML(2), ISGTH, IMO, IDIFF,
     :NT1,NT2,NT3,NT4,NT5,NSV,NMONT,
     :EP,DT3,DD,ID,ISX,DT1,DT2,CT,ZS,ALPHA,CH,CC,CN,CO, 
     :CNE,CMG,CSI,CFE,RCD,OS,RML,RMG,ECA,XF,DR,RMT,RHL,AC,AK1,AK2,ECT,
     :TRB,
     :IRAM, IRS1, VROT1, IRS2, VROT2, FMAC, FAM,
     :IVMC, TRC1, IVMS, TRC2, MWTS, IAGB, ISGFAC, FACSGMIN, SGTHFAC,
     :TENG, SMASS, FMASS, INJMD, STARTTIMEINJ, ENDTIMEINJ,
     :ENDAGE
      WRITE (333, 99005)
      WRITE (333,*) '==========================TZO CONTROL FILE START=========================='
      WRITE (333, 99101) itzo_yn,itzo_cmass_pre, itzo_stripcorehe, itzo_stophighburn,
     :          itzo_noneutburn, itzo_zerocore,
     :          itzo_ct_1, itzo_ct_2, itzo_ct_3
      WRITE (333,*)
      WRITE (333,99102) rtzo_mod_emass, rtzo_dcmassdt,
     :          rtzo_maxdt,
     :          rtzo_cut_non_degen_hburn,
     :          rtzo_EC,
     :          rtzo_nucap, rtzo_nucap_per_yr, rtzo_nucap_min,
     :          rtzo_degen_cutoff,rtzo_degen_cutoff_per_yr, rtzo_degen_cutoff_max,
     :          rtzo_RCD_per_yr, rtzo_RCD_max,
     :          rtzo_meshfluid, rtzo_meshfluid_per_yr, rtzo_meshfluid_min,
     :          rtzo_alpha, rtzo_alpha_per_yr, rtzo_alpha_max,
     :          rtzo_ct_1, rtzo_ct_1_per_yr, rtzo_ct_1_max,
     :          rtzo_ct_2, rtzo_ct_2_per_yr, rtzo_ct_2_max,
     :          rtzo_ct_3, rtzo_ct_3_per_yr, rtzo_ct_3_max
      WRITE (333,*) '==========================TZO CONTROL FILE END=========================='
      WRITE (333, 99005) SM, DTY, AGE, PER, BMS, EC,NH,NP,NMOD,IB,PMH(1),PME(1)
      
C Print out some basic details about the run
      WRITE (*,*) '================================================================================================================================'
      WRITE (*,*) '_______ _______ _______ ______   ______ _____ ______   ______ _______      _______ _______ _______  ______ _______'
      WRITE (*,*) '|       |_____| |  |  | |_____] |_____/   |   |     \ |  ____ |______      |______    |    |_____| |_____/ |______'
      WRITE (*,*) '|_____  |     | |  |  | |_____] |    \_ __|__ |_____/ |_____| |______      ______|    |    |     | |    \_ ______|'
                                                                                                                   
      WRITE (*,*) '================================================================================================================================'
      WRITE (*,*) '================================================================================================================================'

      IF (SM.EQ.1) WRITE (*,*) 'Initial Mass of Star: ', SM, 'Solar Mass'
      IF (SM.NE.1) WRITE (*,*) 'Initial Mass of Star: ', SM, 'Solar Masses'
      
      WRITE (*,*) 'Minimum Allowed Star Mass: ', MINMASS, 'Msun'
      WRITE (*,*) 'Maximun Allowed Star Mass: ', MAXMASS, 'Msun'
      
      IF (ENDAGE.LE.10**6) WRITE (*,*) 'Desired Final Age of Star: ', ENDAGE, 'Yrs'
      IF ((ENDAGE.GT.10**6).AND.(ENDAGE.LT.10**9)) WRITE (*,*) 'Desired Final Age of Star: ', ENDAGE / 10**6, 'Myrs'
      IF (ENDAGE.GE.10**9) WRITE (*,*) 'Desired Final Age of Star: ', ENDAGE / 10**9, 'Gyrs'
      WRITE (*,*) 'Maximum Number of Allowed Models: ', NP
      WRITE (*,*) '===================================================================================================='
      
           
C Print out the artifical energy stuff, see if it works
      IF (TENG.GT.0) THEN
      WRITE (*,*) '===================================================================================================='
      WRITE (*,*) 'Desired Total Artificial Energy Injection: ', TENG
      WRITE (*,*) 'Desired Mass Below Base of Artificial Energy Region: ', SMASS
      WRITE (*,*) 'Desired Mass Below Top of Artificial Energy Region: ', FMASS
      
      IF (STARTTIMEINJ.LT.10**6) WRITE (*,*) 'Artificial Energy Injection to Start at: ', STARTTIMEINJ, 'Yrs'
      IF ((STARTTIMEINJ.GT.10**6).AND.(STARTTIMEINJ.LT.10**9)) WRITE (*,*) 'Artificial Energy Injection to Start at: ', STARTTIMEINJ / 10**6, 'Myrs'
      IF (STARTTIMEINJ.GE.10**9) WRITE (*,*) 'Artificial Energy Injection to Start at: ', STARTTIMEINJ / 10**9, 'Gyrs'
      
      IF (ENDTIMEINJ.LT.10**6) WRITE (*,*) 'Artificial Energy Injection to End at: ', ENDTIMEINJ, 'Yrs'
      IF ((ENDTIMEINJ.GT.10**6).AND.(ENDTIMEINJ.LT.10**9)) WRITE (*,*) 'Artificial Energy Injection to End at: ', ENDTIMEINJ / 10**6, 'Myrs'
      IF (ENDTIMEINJ.GE.10**9) WRITE (*,*) 'Artificial Energy Injection to End at: ', ENDTIMEINJ / 10**9, 'Gyrs'
      
      IF (INJMD.EQ.0) WRITE (*,*) 'Artificial Energy Injection Profile in Mass: ', 'TOP-HAT'
      IF (INJMD.EQ.1) WRITE (*,*) 'Artificial Energy Injection Profile in Mass: ', 'TRIANGULAR'
      IF (INJMD.EQ.2) WRITE (*,*) 'Artificial Energy Injection Profile in Mass: ', 'SINE'
      IF (INJMD.EQ.3) WRITE (*,*) 'Artificial Energy Injection Profile in Mass: ', 'EXP'
      WRITE (*,*) '===================================================================================================='
      ENDIF
C Convert RML from eta to coefficient required
      RML = 4d-13*RML
*     
* Create the spline interpolation data.
* 
      IF (IOP .EQ. 1) CALL OPSPLN
!extra lines for COopac bit
      WRITE (*,*) '===================================================================================================='
      write(*,*) 'Selection for opacity is:',IOP
      IF (IOP.EQ.1) WRITE (*,*) 'Using V4 Style Extended Opacity Tables'
      cbase=ZS*CC
      obase=ZS*CO
      fZ=ZS
C READ IN NEW OPACITY DATA and SETUP STUFF - JJ 4/11/02
c     Read in Opal Data
c     Setup format statements
99042 FORMAT (F5.2, 31F7.3)
99043 FORMAT (5F7.3)
99045 FORMAT (3F7.3)
      if(IOP.ge.2) then
         write(*,*) 'Reading in base tables and setting up splines'
         do I=1,141
            opT(I)=3d0+0.05d0*(I-1)
         enddo
         do J=1,31
            opR(J)=-8d0+0.5d0*(J-1)
         enddo
c     Load in CO tables and setup splines
         write(*,*) 'Reading in Variable tables and setting up splines'
*         OPEN(10,FILE='COtables',STATUS='unknown',ACCESS='SEQUENTIAL')
C Convert the metallicity to a string, so that we can read the
C correct opacity table automatically, it's mad that this wasn't a thing before 5 and 6
         WRITE (COtableNAME, '(F6.5)') ZS
         COtableNAME = COtableNAME(2:6)
         WRITE (*,*) 'Reading in CO Table: ', COtableNAME
         COtableNAME = 'dat/COtables_z'//COtableNAME
         OPEN (UNIT=10, FILE = COtableNAME, ACCESS = 'SEQUENTIAL'
     :   ,STATUS = 'OLD', ERR=1001)
         do K=1,305
            read(10,99045) b3,b1,b2
C            write (*,*) K,b3
            do I=1,141
               read (10,99042) temp,(opacCO(1,1,I,J,K),J=1,31)
            enddo
         enddo
         CLOSE(10)
         WRITE (*,*) '===================================================================================================='
c     Bit to add in variable molecular bits from old paper in Marigo
c     Setup composition matrix
         if(IOP.eq.4.or.IOP.eq.6) then
            write(*,*) "Not Available"
         endif
c     Setup CO spline tables
         do K=1,305
C            write(*,*) K
c     Construct splines in T direciton
            do J=1,31
               do I=1,141
                  MAT(1,I)=opacCO(1,1,I,J,K)
               enddo
               call spline(141,opT,MAT)
               do I=1,140
                  opacCO(2,1,I,J,K)=MAT(2,I)
                  opacCO(3,1,I,J,K)=MAT(3,I)
                  opacCO(4,1,I,J,K)=MAT(4,I)
               enddo
            enddo
c     Construct splines in R direction
            do I=1,140
               do IC=1,4
                  do J=1,31
                     MAT(1,J)=opacCO(IC,1,I,J,K)
                  enddo
                  MAT(2,1)=0d0
                  MAT(3,1)=0d0
                  MAT(2,31)=0d0
                  MAT(3,31)=0d0
                  call spline(31,opR,MAT)
                  do J=1,30
                     opacCO(IC,2,I,J,K)=MAT(2,J)
                     opacCO(IC,3,I,J,K)=MAT(3,J)
                     opacCO(IC,4,I,J,K)=MAT(4,J)
                  enddo
               enddo
            enddo
         enddo
C         write (*,*) 'DONE!!'
      endif
                                !! END of new opac tables bit - JJ - 4/11/02
*
* If IAM=0, use integer atomic weights
*
      IF(IAM.EQ.0)THEN
         DO J = 1, 9
            AM(J) = BN(J)
         END DO
      END IF
C Read the initial model         
      DO  K = 1, NH
         READ (30, 99002) (H(J,K), J=1, JIN)
      END DO
* If available, read initial (last converged) changes
      DO K = 1, NH
         READ (30, 99002, END = 61, ERR = 61) (DH(J,K), J=1, JIN)
         DO 15 J = 1,JIN
            DHPR(J,K) = DH(J,K)
   15    CONTINUE
      END DO
 61   CONTINUE
C Read in first line of star 2 - most of this gets ignored
      IF (IMODE.EQ.2) THEN
         READ (50, 99005) SM2, DTY2, AGE2, PER2, BMS2, EC2,NNH2,NP2,NMOD2,IB2,PMH(2),PME(2)
         DO K = 1, NH
            READ (50, 99002) (H(J,K), J=16, JIN+15)
         END DO
      END IF
C Attempt to read in nucleosynthesis input - but don't worry if it doesn't exist.
      DO J = 1, 100
         DO K = 1, NH
            HNUC(J,K) = 0d0
            DHNUC(J,K) = 0d0
         END DO
      END DO
C Star 1 nucleosynthesis data
      READ (31, 99005, ERR = 12, END = 12)
      DO K = 1, NH
         READ(31, 99002, ERR = 12, END = 12) (HNUC(J,K), J=1, 50)
      END DO
      IF (IMODE.EQ.2) THEN
C Star 2 nucleosynthesis data
         READ (51, 99005, ERR = 12, END = 12)
         DO K = 1, NH
            READ(51, 99002, ERR = 12, END = 12) (HNUC(J,K), J=51, 100)
         END DO
      END IF
***********************************************************************
***********HELIUM FLASH MODEL READ*************************************      
!   42 CONTINUE
!      FILENAMESTORE(1) = 'src/flashdata/modin_for_3_msun' ! Unit 100
!      FILENAMESTORE(2) = 'src/flashdata/data_for_mass_strip'! Unit 101
!      FILENAMESTORE(3) = 'src/flashdata/data_for_core_growth'! Unit 102
!      FILENAMESTORE(4) = 'src/flashdata/nucmodin_for_3_msun'! Unit 103
!      IF (IEND.EQ.-3) GOTO 43
!C Read the mass loss bit****
!C Read the data file
!      OPEN(UNIT = 101,FILE=FILENAMESTORE(2),ACCESS='sequential',STATUS='old')
!      READ  (101,99003) NH2,IT1,IT2,JIN,JOUT,NCH,JP,IZ,IMODE,
!     :ICL,ION,IAM,IOP,IBC,INUC,ICN,IML(1),IML(2),ISGTH,IMO,IDIFF,
!     :NT1,NT2,NT3,NT4,NT5,NSV,NMONT,
!     :EP,DT3,DD,ID,ISX,DT1,DT2,CT,ZS,ALPHA,CH,CC,CN,CO, 
!     :CNE,CMG,CSI,CFE,RCD,OS,RML,RMG,ECA,XF,DR,RMT,RHL,AC,AK1,AK2,ECT,
!     :TRB,
!     :IRAM, IRS1, VROT1, IRS2, VROT2, FMAC, FAM,
!     :IVMC, TRC1, IVMS, TRC2, MWTS, IAGB, ISGFAC, FACSGMIN, SGTHFAC,
!     :TENG, SMASS, FMASS, INJMD, STARTTIMEINJ, ENDTIMEINJ,
!     :ENDAGE
!      CLOSE(101)
!C Read the initial 3 solar mass model
!      OPEN(UNIT=100,FILE=FILENAMESTORE(1),ACCESS='sequential',STATUS='old')
!      DO  K = 1, NH
!         READ (100, 99002) (H(J,K), J=1, JIN)
!      END DO
!      CLOSE(100)  
!      GOTO 44  
!   43 CONTINUE
!C Read the core grow bit****** 
!      OPEN(UNIT=102,FILE=FILENAMESTORE(3),ACCESS='sequential',STATUS='old')
!      READ  (102,99003) NH2,IT1,IT2,JIN,JOUT,NCH,JP,IZ,IMODE,
!     :ICL,ION,IAM,IOP,IBC,INUC,ICN,IML(1),IML(2),ISGTH,IMO,IDIFF,
!     :NT1,NT2,NT3,NT4,NT5,NSV,NMONT,
!     :EP,DT3,DD,ID,ISX,DT1,DT2,CT,ZS,ALPHA,CH,CC,CN,CO, 
!     :CNE,CMG,CSI,CFE,RCD,OS,RML,RMG,ECA,XF,DR,RMT,RHL,AC,AK1,AK2,ECT,
!     :TRB,
!     :IRAM, IRS1, VROT1, IRS2, VROT2, FMAC, FAM,
!     :IVMC, TRC1, IVMS, TRC2, MWTS, IAGB, ISGFAC, FACSGMIN, SGTHFAC,
!     :TENG, SMASS, FMASS, INJMD, STARTTIMEINJ, ENDTIMEINJ,
!     :ENDAGE
!      CLOSE(102)
!C Read the 3 solar mass nucleosynthsis bit TODO
!   44 CONTINUE
!      OPEN(UNIT=103,FILE=FILENAMESTORE(4),ACCESS='sequential',STATUS='old')
!      READ (103, 99005)
!      DO K = 1, NH
!         READ(103, 99002) (HNUC(J,K), J=1, 50)
!      END DO
!      CLOSE(103)   
!      GOTO 33
*************FINISHED THE HELIUM FLASH BIT*****************************
    
C Convert some things to `cgs' units: 10**11 cm, 10**33 gm, 10**33 erg/s 
   12 DT = CSECYR*DTY
      TM(1) = MSUN*SM
      TM(2) = MSUN*SM2
      RMG = RMG/CSECYR
      RMT = MSUN*RMT/CSECYR
      RML = RML*MSUN**2/LSUN/RSUN/CSECYR
C Optionally, re-initialise mass
! Old Version with mass resolution bug!!
      IF ( NCH.GE.1 ) THEN
         H(4,1) = DLOG(TM(1)) 
         HPR(4,1) = DLOG(TM(1)) 
         IF (IMODE.EQ.2) THEN
            H(19,1) = DLOG(TM(2))
            HPR(19,1) = DLOG(TM(2))
         END IF
      END IF
      IF (IMODE.EQ.2) THEN
         BM = MSUN*(SM + SM2)   !MSUN*BMS
      ELSE
         BM = MSUN*BMS
      END IF
      M0=BM-TM(1)
      ANG = TM(1)*(BM-TM(1))*(3.55223D0*PER/BM)**(1D0/3D0)
*
* CAT 30th August 2009
* This is the problematic bit for remeshing because TM(1/2) just doesn't
* have the required accuracy!
*
*      IF ( NCH.GE.1 ) THEN
*         IF ( NCH .LE. 3 ) THEN
*            NEWMS = DLOG(TM(1))
*            OLDMS = H(4,1)
*            NEWMSR = NEWMS
*            OLDMSR = OLDMS
*
* CAT 7th September 2009
* Scale all mass shells when changing mass.  This replicates the old
* situation where relative mass was stored rather than actual, or so we think
*
*            DO SHELL = 1,NH
*               H(4,SHELL) = H(4,SHELL)*NEWMS/OLDMS
*               HPR(4,SHELL) = HPR(4,SHELL)*NEWMSR/OLMSR
*            ENDDO
*            IF (IMODE.EQ.2) THEN
*               H(19,1) = DLOG(TM(2))
*               HPR(19,1) = DLOG(TM(2))
*            END IF
*         ELSE
*            NCH = NCH - 2
*         END IF
*      END IF
*      IF (IMODE.EQ.2) THEN
*         BM = MSUN*(SM + SM2)   !MSUN*BMS
*      ELSE
*         BM = MSUN*BMS
*      END IF
*      M0=BM-TM(1)
*      ANG = TM(1)*(BM-TM(1))*(3.55223D0*PER/BM)**(1D0/3D0)


 13   CONTINUE
  102 FORMAT(2(F8.4,1PE16.9,4E10.3,0P7F8.5,F8.3,2F8.4,/),3F8.4,F10.4, 
     :1P2E10.3,0PF10.4,7F8.5,F8.3,2F8.4,/,F8.4,12F8.3,5F8.4,/,3F8.4,I6)
C REMESH optionally rezones the model, e.g. for different no. of meshpoints
   14 CALL REMESH ( NH2, NCH, CH, CO, CC, CNE )
      DO K=1,NH
         DO J = 1,JOUT
            HPR(J,K) = H(J,K)
         END DO
C Store nucleosynthesis
         DO J=1,100
            HNUCPR(J,K) = HNUC(J,K)
         END DO
      END DO
C COMPOS puts composition variables to zero if they are very small
      CALL COMPOS
      CALL PRINTB ( DTY, PER, NT1, NT2, NT3, NT4, NMONT, NMOD, IEND )
*
* If initial timestep is negative calculate DT as a fraction of the
* Kelvin-Helmholtz timescale and scale mass loss to evolve up the main
* sequence.
*
C Does this still work??? I never use it...
      IF (ICN .EQ. 1) DT = CSECYR*5D0*TKH(ISTAR)/(SM*SM)
      IF (ICN .EQ. 1) RMG = CLN10/(DT*2D2)
      IHOLD = 4
      CLOSE (17)
      CLOSE (20)
      CLOSE (25)
c Store certain previous values, for possible emergency restart
      PR(1) = AGE
      PR(2) = DT
      PR(3) = M1
      PR(4) = EC
      PR(5) = BM
      PR(6) = ANG      
      PR(7) = CM
      PR(8) = MTA
      PR(9) = MTB
      PR(10) = TM(1)
      PR(11) = T0
      PR(12) = M0
      PR(13) =  TC(1)
      PR(14) = TC(2)
C      DO 16 J = 1,13
C   16    PR(J) = CT(J+10)
      ! Store the "real" tzo values
      S_DTZO(1) = rtzo_mod_emass
      S_DTZO(2) = rtzo_dcmassdt
      S_DTZO(3) = rtzo_maxdt
      S_DTZO(4) = rtzo_cut_non_degen_hburn
      S_DTZO(5) = rtzo_EC
      S_DTZO(6) = rtzo_nucap
      S_DTZO(7) = rtzo_nucap_per_yr
      S_DTZO(8) = rtzo_nucap_min
      S_DTZO(9) = rtzo_degen_cutoff
      S_DTZO(10) = rtzo_degen_cutoff_per_yr
      S_DTZO(11) = rtzo_degen_cutoff_max
      S_DTZO(12) = rtzo_RCD_per_yr
      S_DTZO(13) = rtzo_RCD_max
      S_DTZO(14) = rtzo_meshfluid
      S_DTZO(15) = rtzo_meshfluid_per_yr
      S_DTZO(16) = rtzo_meshfluid_min
      S_DTZO(17) = rtzo_alpha
      S_DTZO(18) = rtzo_alpha_per_yr
      S_DTZO(19) = rtzo_alpha_max
      
      S_DTZO(20) = rtzo_ct_1
      S_DTZO(21) = rtzo_ct_1_per_yr
      S_DTZO(22) = rtzo_ct_1_max
      
      S_DTZO(23) = rtzo_ct_2
      S_DTZO(24) = rtzo_ct_2_per_yr
      S_DTZO(25) = rtzo_ct_2_max
      
      S_DTZO(26) = rtzo_ct_3
      S_DTZO(27) = rtzo_ct_3_per_yr
      S_DTZO(28) = rtzo_ct_3_max
      
      NPR = NMOD
      KPR = KS
      GOTO 40
C Almost end of initial input section. Start of regular update section
   30 IF ( IEND.NE.0 ) GO TO 31
      CALL COMPOS
c Store certain previous values, for possible emergency restart
      PR(1) = AGE
      PR(2) = DT
      PR(3) = M1
      PR(4) = EC
      PR(5) = BM
      PR(6) = ANG      
      PR(7) = CM
      PR(8) = MTA
      PR(9) = MTB
      PR(10) = TM(1)
      PR(11) = T0
      PR(12) = M0
      PR(13) =  TC(1)
      PR(14) = TC(2)
C      DO 10 J = 1,13
C   10    PR(J) = CT(J+10)
      S_DTZO(1) = rtzo_mod_emass
      S_DTZO(2) = rtzo_dcmassdt
      S_DTZO(3) = rtzo_maxdt
      S_DTZO(4) = rtzo_cut_non_degen_hburn
      S_DTZO(5) = rtzo_EC
      S_DTZO(6) = rtzo_nucap
      S_DTZO(7) = rtzo_nucap_per_yr
      S_DTZO(8) = rtzo_nucap_min
      S_DTZO(9) = rtzo_degen_cutoff
      S_DTZO(10) = rtzo_degen_cutoff_per_yr
      S_DTZO(11) = rtzo_degen_cutoff_max
      S_DTZO(12) = rtzo_RCD_per_yr
      S_DTZO(13) = rtzo_RCD_max
      S_DTZO(14) = rtzo_meshfluid
      S_DTZO(15) = rtzo_meshfluid_per_yr
      S_DTZO(16) = rtzo_meshfluid_min
      S_DTZO(17) = rtzo_alpha
      S_DTZO(18) = rtzo_alpha_per_yr
      S_DTZO(19) = rtzo_alpha_max
      
      S_DTZO(20) = rtzo_ct_1
      S_DTZO(21) = rtzo_ct_1_per_yr
      S_DTZO(22) = rtzo_ct_1_max
      
      S_DTZO(23) = rtzo_ct_2
      S_DTZO(24) = rtzo_ct_2_per_yr
      S_DTZO(25) = rtzo_ct_2_max
      
      S_DTZO(26) = rtzo_ct_3
      S_DTZO(27) = rtzo_ct_3_per_yr
      S_DTZO(28) = rtzo_ct_3_max
      
      NPR = NMOD
      KPR = KS
C PRINTB prints out every NT2'th meshpoint of every NT1'th model; NT3
C `pages' per printed model; also 4-line summary for every NT4'th model
      DTY = DT/CSECYR
      AGE = AGE + DTY
      NMOD = NMOD + 1
      PPER = PER
      CALL PRINTB ( DTY, PER, NT1, NT2, NT3, NT4, NMONT, NMOD, IEND )
      ANG = ANG/(1.0D0 + RHL*DTY)
      EC = EC*(1.0D0 + DTY*ECT)/(1.0D0 - DTY*ECA*EC) 
      
C TZO Evolve TZO control parameters with time, as needed...      
      IF (itzo_yn.EQ.1) THEN
        ! Do the neutron luminosity cap rtzo_nucap, rtzo_nucap_per_yr, rtzo_nucap_min
        rtzo_nucap = rtzo_nucap * (1.D0 + DTY*rtzo_nucap_per_yr)
        IF (rtzo_nucap.LT.rtzo_nucap_min) rtzo_nucap = rtzo_nucap_min
        
        ! Do the RCD (Convective turnover timescale)
        RCD = RCD * (1.D0 + DTY*rtzo_RCD_per_yr)
        IF (RCD.GT.rtzo_RCD_max) RCD = rtzo_RCD_max
        
        ! Do the Mesh Fluidity stuff
        rtzo_meshfluid = rtzo_meshfluid * (1.D0 + DTY*rtzo_meshfluid_per_yr)
        rtzo_meshfluid = DMAX1(0.D0, DMIN1(1.D0, rtzo_meshfluid))
      
      
      ENDIF
      
C TZO, evolve the MSF coefficients, as needed 
      IF (itzo_yn.EQ.1) THEN
        
        IF(itzo_ct_1.GT.0) THEN
            CT(itzo_ct_1) = MAX(rtzo_ct_1, 
     :      MIN(rtzo_ct_1_max, CT(itzo_ct_1) + DTY*rtzo_ct_1_per_yr))
        ENDIF
        
        IF(itzo_ct_2.GT.0) THEN
            CT(itzo_ct_2) = MAX(rtzo_ct_2, 
     :      MIN(rtzo_ct_2_max, CT(itzo_ct_2) + DTY*rtzo_ct_2_per_yr))
        ENDIF
            
        IF(itzo_ct_3.GT.0) THEN
            CT(itzo_ct_3) = MAX(rtzo_ct_3, 
     :      MIN(rtzo_ct_3_max, CT(itzo_ct_3) + DTY*rtzo_ct_3_per_yr))
        ENDIF
      
      ENDIF
*     TRB = TRB*(1.0D0 + DTY*ECT)
C FUDGE TO DEAL WITH KS outside range. THIS ***WILL*** SCREW UP BINARIES!
C This is no longer used, so I don't care whether it works or not. RJS
      KS = 1
      IF ( AGE.GT.ST(KS+1) ) KS = KS+1
      IF ( IB.EQ.2 .AND. (ST(KS+2).EQ.0.0D0.OR.RLF.GT.0.0D0) ) STOP
      DELTA = 0.0D0
      DO 5 K = 1, NH
         DO 4 J = 1, 30 !60
C Don't use L, HORB in delta
            IF (J.NE.8.AND.J.NE.23.AND.J.NE.13.AND.J.NE.28.AND.J.NE.14.AND.J.NE.29) then
               DELTA = DELTA + DABS(DH(J,K))
            end if
            HPPR(J,K) = HPR(J,K)
            HPR(J, K) = H(J, K)
            DHPR(J,K) = DH(J,K)
    4       H(J, K) = H(J, K) + DH(J, K)
    5       IF (DTY.GT.3D-4) DELTA = DELTA + AC*DABS(DH(8,K)/HPR(8,1))
C Update nucleosynthesis matrix
      DO K = 1, NH
         DO J=1,100
            HNUCPR(J, K) = HNUC(J, K)
            DHNUCPR(J, K) = DHNUC(J, K)
            HNUC(J, K) = HNUC(J, K) + DHNUC(J, K)
C Blank DHNUC each time
C            DHNUC(J,K) = 0d0
         END DO
      END DO
      write (333,*) "DELTA =",DELTA, " DD = ", DD
      DTF = DMIN1 (DT2, DD/DELTA)      
C     IF ( IHOLD .LE. 3 ) DTF = 1.0D0
      IF ( IHOLD .LE. 2 ) DTF = 1.0D0
      DTY = DMAX1(DT1,DTF)*DTY
      DO ISTAR = 1,IMODE
         IF (IAGB.EQ.1) THEN
C Reduce timestep if He luminosity is increasing too fast -- useful on AGB
            IF ((VLEC(ISTAR) - VLEP(ISTAR))/VLEP(ISTAR).GT.0.05.AND.
     :           VLEC(ISTAR).GT.1d3) DTY = 0.8*DTY
         END IF
         
C Trying this again, we need to reduce the timestep while we're injecting energy
         IF ((AGE.GE.STARTTIMEINJ).AND.(AGE.LE.ENDTIMEINJ)) THEN
            IF (DTY.GT.((ENDTIMEINJ - STARTTIMEINJ)/100)) THEN
                DTY = 0.8 * DTY
*                WRITE (*,*) 'Timestep Cut'
            ENDIF
         ENDIF
         
C TZO, max timestep
         IF (itzo_yn.EQ.1) THEN
            IF (DTF*DTY.GT.rtzo_maxdt) THEN
                DTF = rtzo_maxdt / DTY
                DTY = DTF*DTY
            ENDIF
         ENDIF
         
         
         
         
         
C Control mechanism for dealing with artificial energy generation
C Works to ensure that the timestep is of a reasonable size by the time
C We reach the injection time
C Need to see if we're within 50 current timesteps of the injection starttime
C Make sure we're not already done with the injection though!!

*         IF (((AGE + (DTY * 100)).GT.(STARTTIMEINJ)).AND.(AGE.LT.ENDTIMEINJ).AND.(AGE.LT.STARTTIMEINJ)) THEN
*            WRITE (*,*) 'APPROACHING INJECTION TIMESTEP'
C           Work out how many "fine" we need the timestep
C           At the time of injection, we want the timestep to be 
C           at max, 1/100 of the injection time
*            IF ((DTY*100).LE.(ENDTIMEINJ - STARTTIMEINJ)) GOTO 96
*            DTY = DTY * 80 ! Cut DTY to 80% to ensure we get there!   
*            WRITE (*,*) 'TIMESTEP CUT TO: ', DTY
*            
*         ENDIF
C        Okay, have we reached the injection timesteps?         
*         IF ((AGE.GT.STARTTIMEINJ).AND.(AGE.LT.ENDTIMEINJ)) THEN  
*            WRITE (*,*) 'INJECTING ENERGY'
C           Okay, is our timestep small enough? 
*            IF (DTY.GT.((ENDTIMEINJ - STARTTIMEINJ) / 100)) THEN
*                DTY = (ENDTIMEINJ - STARTTIMEINJ) / 100
*                WRITE (*,*) '!!TIMESTEP TOO LARGE AT INJECTION, REDUCING NOW!!'
*            ENDIF         
*         ENDIF
*   96    CONTINUE
*         WRITE (*,*) 'Current Timestep: ', DTY, 'Yrs'
C Extra control mechanisms that can be uncommented as necessary - RJS
C         IF ((VLHC(ISTAR) - VLHP(ISTAR))/VLHP(ISTAR).GT.0.10) DTY = 0.8*DTY
C         IF ((VLCC(ISTAR) - VLCP(ISTAR))/VLCP(ISTAR).GT.0.10) DTY = 0.8*DTY
C         IF (((RLFC(ISTAR)-RLFP(ISTAR))/RLFP(ISTAR)).GT.0.1
CC     :        .AND.RLFP(ISTAR).GT.-1d-1.AND.RLFP(ISTAR).LT.-1d-3) THEN
C     :        .AND.RLFP(ISTAR).GT.-1d-1) THEN
C            DTY = 0.1*DTY
C            write (333,*) "RLF issues - reducing timestep"
C         END IF
      END DO
C      IF (DABS((PER - PPER)/PPER).GT.0.01) DTY = 0.5*DTY
      IF ( IB.EQ.2 ) DTY = DMIN1(DTY,ST(KS+2)-AGE)
      DT=CSECYR*DTY
C      IF (DT1.EQ.1d0) GO TO 6
C      IF (IDREDGE.EQ.3) GO TO 6      
      IF ( (JP.EQ.1 .AND. DTF.GE.DT1).OR.DTY.LT.6d-5 ) GO TO 6
c clear DH in some circumstances
      WRITE (333,*) "Clearing DH..."
      DO 7 K = 1, NH
         DO 7 J = 1, 60
    7       DH(J, K) = 0.0D0
    6 IHOLD = IHOLD + 1 
c CNO equilibrium on the main sequence.
C Does this still work???
      IF (ICN .EQ. 1) DT = MSUN*MSUN*CSECYR*5D0*TKH(ISTAR)/(VM*VM)
      IF (ICN .EQ. 1) RMG = CLN10/(DT*2D2)
      CALL COMPOS
      
C TZO Stuff
!      OPEN (UNIT=70, FILE='tzo_data', ACCESS='SEQUENTIAL', STATUS='OLD')
!      READ (70, 99101) itzo_yn,itzo_cmass_pre, itzo_stripcorehe, itzo_stophighburn,
!     :          itzo_noneutburn, itzo_zerocore
!      READ (70, 99102) rtzo_mod_emass, rtzo_dcmassdt, rtzo_maxdt,
!     :          rtzo_cut_non_degen_hburn, rtzo_EC, rtzo_nucap,
!     :          rtzo_nucap_per_yr, rtzo_nucap_min, rtzo_degen_cutoff,
!     :          rtzo_degen_cutoff_per_yr, rtzo_degen_cutoff_max
!      CLOSE (70)
* Modify TZO modified electron mass bit
      IF (itzo_yn.EQ.1) THEN
        dlcm = rtzo_dcmassdt * DTY
        cm = MIN(rtzo_mod_emass * (1.D0 + dlcm), 1839.D0)
        rtzo_mod_emass = cm
        
* modify degeneracy param that we switch to neutronic material at
        rtzo_degen_cutoff = MIN(rtzo_degen_cutoff + 
     &      DTY*rtzo_degen_cutoff_per_yr, rtzo_degen_cutoff_max)
      ENDIF
      
C TZO stuff, need to zero out compositions in the "core"
      IF (itzo_yn.EQ.1) THEN
        IF (itzo_zerocore.EQ.1) THEN
            MCORE = VMH + VME
            OLDMCORE = MCORE
            DO JJ = 1, NH
                J = NH + 1 - JJ
                XM = DEXP(H(4,J))
                
                IF (XM.LT.MCORE) THEN
                    H(5,J) = 0.D0 ! Zero H
                    H(9,J) = 0.D0 ! Zero He
                    H(10,J) = 0.D0 ! Zero C12
                ELSE IF ((H(5,J) + H(9,J) + H(10,J)).EQ.0.D0) THEN
                    MCORE = XM
                ENDIF
            ENDDO

        ENDIF
      ENDIF
      
      
c For *2, some factors relating to accretion from *1. Ignored if this is *1
 40   CONTINUE
C   40 T0 = CSECYR*ST(KS+1)
C      M0 = MSUN*MS(KS+1)
C RJS added to allow -C compile and run
C      IF (KS.NE.0) THEN
C         MTA = MSUN/CSECYR*(MS(KS+1)-MS(KS))/(ST(KS+1)-ST(KS))
C         MTB = MSUN/CSECYR*(MS(KS+2)-MS(KS+1))/(ST(KS+2)-ST(KS+1))
C      END IF
      IF ( MOD(NMOD,NSV).NE.0 .OR. IEND.EQ.-1 ) RETURN
C End of regular update section. Intermediate or final output section
   31 IF (IEND.EQ.2) GO TO 32
      SM = DEXP(H(4,1))/MSUN 
      SM2 = DEXP(H(19,1))/MSUN
      BMS = (SM + SM2) !BM/MSUN
      IF (IMODE.EQ.1) BMS = BM/MSUN
      DTY = DT/CSECYR
      
C Restart bit!!      

C=====================================================================
C     This write the final files to the last modin and modout files
      OPEN (UNIT=134,FILE='temp_last_modin',ACCESS='SEQUENTIAL',STATUS='REPLACE')
      WRITE (134, 99005) SM, DTY, AGE, PER, BMS, EC,NH,NP,NMOD,IB,PMH(1),PME(1)
      DO K = 1, NH
        WRITE (134, 99002) (H(J,K), J = 1, 15)
      ENDDO  
      DO K = 1,NH
        WRITE (134, 99002) (DH(J,K), J=1, 15)
      ENDDO
      FLUSH (134)
      CLOSE (134)
      
      OPEN (UNIT=135,FILE='temp_last_nucmodin',ACCESS='SEQUENTIAL',STATUS='REPLACE')
      WRITE (135, 99005) SM, DTY, AGE, PER, BMS, EC,NH,NP,NMOD,IB,PMH(1),PME(1)
      DO K = 1, NH
         WRITE (135, 99002) (HNUC(J,K),J=1,50)
      END DO
      FLUSH (135)
      CLOSE (135)
C=====================================================================       
      
      WRITE (34, 99005) SM, DTY, AGE, PER, BMS, EC,NH,NP,NMOD,IB,PMH(1),PME(1)
      DO K = 1, NH
C Need to do this better - at present I'm writing out blanks
        WRITE (34, 99002) (H(J,K), J=1, 15)
      END DO
      DO K = 1, NH
         WRITE (34, 99002) (DH(J,K), J=1, 15)
      END DO
      CALL FLUSH(34)
      IF (IMODE.EQ.2) THEN
C Write out star 2
         WRITE (54, 99005) SM2, DTY, AGE, PER, BMS, EC,NH,NP,NMOD,IB,PMH(2),PME(2)
         DO K = 1, NH
C Copy HORB from star 1
            H(28,K) = H(13,K)
            DH(28,K) = H(28,K)
C Need to do this better - at present I'm writing out blanks
            WRITE (54, 99002) (H(J,K), J=16, 30)
         END DO
         DO K = 1, NH
            WRITE (54, 99002) (DH(J,K), J=16, 30)
         END DO
         CALL FLUSH(54)
      END IF
C write out nucleosynthesis files - star 1 first
      WRITE (35, 99005) SM, DTY, AGE, PER, BMS, EC,NH,NP,NMOD,IB,PMH(1),PME(1)
      DO K = 1, NH
         WRITE (35, 99002) (HNUC(J,K),J=1,50)
      END DO
      CALL FLUSH(35)
      IF (IMODE.EQ.2) THEN
         WRITE (55, 99005) SM2, DTY, AGE, PER, BMS, EC,NH,NP,NMOD,IB,PMH(1),PME(2)
         DO K = 1, NH
            WRITE (55, 99002) (HNUC(J,K),J=51,100)
         END DO
         CALL FLUSH(55)
      END IF
      RETURN
C End of final output section. Start of emergency restart section 
   32 IF ( IHOLD .LE. 0 ) GO TO 34
      AGE = PR(1)
      DT = PR(2)
      M1 = PR(3)
      EC = PR(4)
      BM = PR(5)
      ANG = PR(6)      
      CM = PR(7)
      MTA = PR(8)
      MTB = PR(9)
      TM(1) = PR(10)
      T0 = PR(11)
      M0 = PR(12)
      TC(1) = PR(13)
      TC(2) = PR(14)
C      DO 11 J = 1,13
C   11    CT(J+10) = PR(J)
      rtzo_mod_emass = S_DTZO(1)
      rtzo_dcmassdt = S_DTZO(2)
      rtzo_maxdt = S_DTZO(3)
      rtzo_cut_non_degen_hburn = S_DTZO(4) 
      rtzo_EC = S_DTZO(5)
      rtzo_nucap = S_DTZO(6)
      rtzo_nucap_per_yr = S_DTZO(7)
      rtzo_nucap_min = S_DTZO(8)
      rtzo_degen_cutoff = S_DTZO(9)
      rtzo_degen_cutoff_per_yr = S_DTZO(10)
      rtzo_degen_cutoff_max = S_DTZO(11)
      rtzo_RCD_per_yr = S_DTZO(12)
      rtzo_RCD_max = S_DTZO(13)
      rtzo_meshfluid = S_DTZO(14)
      rtzo_meshfluid_per_yr = S_DTZO(15)
      rtzo_meshfluid_min = S_DTZO(16)
      rtzo_alpha = S_DTZO(17)
      rtzo_alpha_per_yr = S_DTZO(18)
      rtzo_alpha_max = S_DTZO(19)
      
      rtzo_ct_1 = S_DTZO(20)
      rtzo_ct_1_per_yr = S_DTZO(21)
      rtzo_ct_1_max = S_DTZO(22)
      
      rtzo_ct_2 = S_DTZO(23)
      rtzo_ct_2_per_yr = S_DTZO(24)
      rtzo_ct_2_max = S_DTZO(25)
      
      rtzo_ct_3 = S_DTZO(26)
      rtzo_ct_3_per_yr = S_DTZO(27)
      rtzo_ct_3_max = S_DTZO(28)
      
      NMOD = NPR
   34 DT=0.8*DT
      IF ( DT .LT. 0.01*PR(2) ) STOP
   33 DO 9 K = 1, NH
         DO 9 J = 1, JOUT
            DH(J,K) = DHPR(J,K)
    9       H(J,K) = HPR(J,K) 
            DTOLD = DTOLDP
C Sort out nucleosynthesis for restart
            DO K = 1, NH
               DO J = 1,100
                  DHNUC(J, K) = DHNUCPR(J, K)
                  HNUC(J, K) = HNUCPR(J, K)
               END DO
            END DO
      IHOLD = 0
      RETURN
C I/O error handling
 1001 WRITE (*,*) '!!Specified Metallicity Has No CO Opacity Table!!'
      STOP
      
      
      
      
      END
