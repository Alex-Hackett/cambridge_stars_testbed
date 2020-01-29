C Between the comments added to the code and the file 'bswriteup.tex'
C you can hopefully understand how this all works... RJS 9/6/06
C PPEs original comment:
C The working of this programme is described in a file `writeup.tex';
C Consequently I do not include many comments in the body of the programme.
C Following is the main routine
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MINMASS, MAXMASS
      SAVE
      PARAMETER (MAXMSH = 2000)
*      REAL MSUM
      COMMON H(60,MAXMSH),DH(60,MAXMSH),EPS,V(2),NMESH,JIN,ID(100),IE(100)
*      COMMON /NUCMAT/ HNUC(100,MAXMSH), DHNUC(100,MAXMSH)
      COMMON /EVMODE/ IMODE
      COMMON /SODDS / ALPHA, RML, CMG, CSI, CFE, CT(10), AGE, DT, M1, 
     :  EC, BM, ANG, CM, MTA, MTB, TM(2), T0, M0, TC(2), OS, AC, RCD,
     :  RMG, RHL, XF, DR, AK1 ,RMT, AK2, IZ(4), IB, ISX(45),
     :  TRB
      COMMON /MASTRK/ AMASSSH(MAXMSH), YSHELLS(MAXMSH), SCALEMASS(MAXMSH)
      COMMON /FLASHSTORE/ SM, NM
      COMMON /ARTENG/ TENG, SMASS, FMASS, SHIENG, INJMD
      COMMON /INJTIMECTRL/ STARTTIMEINJ, ENDTIMEINJ
      COMMON /AGECTRL/ ENDAGE
      COMMON /MASSCTRL/ MINMASS, MAXMASS
      SAVE
      COMMON /OP    / ZS, LEDD, VM, GR, GRAD, ETH, RLF, EGR, R, QQ
      COMMON /YUK1  / PX(34), WMH, WMHE, VMH, VME, VMC, VMG, BE, VLH,
     :                VLE, VLC, VLN, VLT, MCB(12),WWW(100)
     
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
      
      real dtime,cpu(2),dt,tcpu
C Read physical data and an initial model
      CALL PRINTA ( -1, NSTEP, ITER1, ITER2, NWRT5 )
      !IF (NSTEP.EQ.0) GO TO 3
      ITER = ITER1
      nter = 0
      nm = 0
      npr = nm
      tcpu = dtime(cpu)
      tcpu = 0.0
      IHAVEFLASH = 0
C Begin evolutionary loop of NSTEP time steps
      WRITE (*,*) 'BEGINNING COMPUTATION'
      CALL CPU_TIME(START)
    2 IF ((NM.EQ.NSTEP).AND.(NSTEP.NE.0)) GOTO 3

C Have we reached the set age?
*      IF (AGE.GT.(1.27038904*10)**10) THEN
*        WRITE (*,*) 'Set Age Reached'
*        GOTO 3
*      ENDIF
         IPRINTAGE = 1
         IF (IPRINTAGE.EQ.1) THEN
         IF (AGE.LT.10**6) WRITE (*,*) 'Solving Model: ', NM, '     Model Age: ', AGE, 'Yrs'
         IF ((AGE.LT.10**9).AND.(AGE.GT.10**6)) WRITE (*,*) 'Solving Model: ', NM, '     Model Age: ', AGE/(10**6), 'Myr'
         IF (AGE.GT.10**9) WRITE (*,*) 'Solving Model: ', NM, '     Model Age: ', AGE/(10**9), 'Gyr'
         ENDIF
         IF ((AGE.GE.ENDAGE).AND.(AGE.NE.0)) THEN
            WRITE (*,*) '!!AGE CUTOFF REACHED: TERMINATING!!'
            GOTO 3
         ENDIF
* Have we reached the max or the min mass?
         IF ((PX(9).GE.MAXMASS).OR.(PX(9).LE.MINMASS)) THEN
            IF (PX(9).GE.MAXMASS) WRITE (*,*) 'MASS TOO LARGE'
            IF (PX(9).LE.MINMASS) WRITE (*,*) 'MASS TOO SMALL'
            WRITE (*,*) '!!STAR MASS REACHED: ', PX(9), 'TERMINATING!!'
            GOTO 3
         ENDIF  
      
         
C Solve for structure, mesh, and major composition variables
         CALL SOLVER ( 1, ITER, KTER, ERR, ID, NWRT5 )
         VLEOLD = VLE
C Grab the old VLE so we can check for the flash
            
C We need to do timestep "messing" in printa.f, otherwise, this will
C be awfully messy, so hope over there now for timestep control regarding
CArtificial energy injection    
C Check if we are at the right age to perform the injection, if not,
C Just set the injection energy to zero and skip all of this 
         IF ((AGE.LT.STARTTIMEINJ).OR.(AGE.GT.ENDTIMEINJ)) THEN
            SHIENG = 0
*            WRITE (*,*) 'Skipping, Age is Wrong'
            GOTO 10
         ENDIF            
C Search the shell mass array AMASSSH, to see where the valid mass coordinates for energy injection are
         DO i = 1, MAXMSH
            IF ((AMASSSH(i).GT.SMASS).AND.(AMASSSH(i).LT.FMASS)) THEN
                YSHELLS(i) = 1
            ELSE
                YSHELLS(i) = 0
            ENDIF
         ENDDO
*         WRITE (*,*) SMASS, FMASS
*         WRITE (*,*) YSHELLS
C     Sum the number of valid shells, to determine how much energy should
C     be injected into each one!
*         WRITE (*,*) SUM(YSHELLS)
         MSUM = 0
         DO i = 1, MAXMSH
            IF (YSHELLS(i).EQ.1) MSUM = MSUM + (AMASSSH(i) - AMASSSH(i + 1))
         ENDDO
         
         SHIENG = TENG / REAL(MSUM)
C     That'll do for the tophat injection, but for the sin
C     function, we'll need to precalc this, store it in the scalemass 
C     array so that we can multiply this in later...   
********SIN INJECTION****************      
         IF (INJMD.EQ.2) THEN
             DO i = 1, MAXMSH
                IF (YSHELLS(i).EQ.1) THEN                 
                    SCALEMASS(i) = SIN(((AMASSSH(i)) - (AMASSSH(i+1)))
     &                * (3.14159 / MAXVAL(AMASSSH)))
                ENDIF
            ENDDO
*            WRITE (*,*) 'SIN Scaled, but not yet scaled by energy', SCALEMASS
*       Sum up these masses, and scale to the energy needed
            EUNSCALE = 0
            DO i = 1, MAXMSH
                IF (YSHELLS(i).EQ.1) THEN
                    EUNSCALE = EUNSCALE + (SCALEMASS(i) * SHIENG)
                ENDIF
            ENDDO
*        EUNSCALE needs to be equal to TENG, but it'll obviously not be!
            SCALEFAC = 0
            SCALEFAC = REAL(TENG) / REAL(EUNSCALE)
*       TEST, print the difference between TENG and EUNSCALE * SCALEFAC
*            WRITE (*,*) 'DEBUG: DIFFERENCE BETWEEN TENG AND EUNSCALE * SCALEFAC SHOULD BE ZERO, is:  ', TENG - (EUNSCALE * SCALEFAC)            
*       SCALEFAC * each element in SCALEMASS will produce a scaled list of masses
*       SCALEMASS will be used instead of AMASSH in funcs1.f to inject energy
            DO i = 1, MAXMSH
                SCALEMASS(i) = SCALEMASS(i) * SCALEFAC
            ENDDO                    
         ENDIF
*         WRITE (*,*) 'SCALED MASSES', SCALEMASS
*         WRITE (*,*) SHIENG
         
            
         
*         WRITE (*,*) 'Masses of Shells', AMASSSH
********** END OF ARTIFICAL ENERGY INJECTION STUFF***********
   10    CONTINUE 
         IF (ERR.LT.EPS) npr = nm
         IF (ERR.LT.EPS) nm = nm + 1
         nter = nter + kter
         dt = dtime(cpu)
         tcpu = tcpu + dt
         IF (ERR.GT.EPS) kter = -kter
C         write(61,99000) nm, kter, dt, dt/kter, nter, tcpu, tcpu/nter,
C     &        tcpu/nm
99000    format(i6,i3,2f8.4,i7,f10.2,2f9.5)
         call flush(31)
C If model didnt converge, restart from 2 steps back with DT halved 
         IF (ERR.GT.EPS) nm = npr ! nm - 2
         IF (ERR.GT.EPS) WRITE (*,*) '!!FAILED TO CONVERGE, HALVING DT!!'
         IF (ERR.GT.EPS) CALL PRINTA ( 2, NSTEP, ITER1, ITER2, NWRT5 ) 
         IF (ERR.GT.EPS) GO TO 1
C If required, solve for minor composition vbles, with mesh, structure fixed.
C First pass for star 1
         CALL SOLVER ( 2, ITER, KTER, ERR, IE, NWRT5 )
C Restart if failed on composition
         IF (ERR.GT.EPS) nm = npr ! nm - 2
         IF (ERR.GT.EPS) WRITE (*,*) '!!COMPOSITION FAILED, RESTARTING!!'
         IF (ERR.GT.EPS) CALL PRINTA ( 2, NSTEP, ITER1, ITER2, NWRT5 ) 
         IF (ERR.GT.EPS) GO TO 1
         IF (IMODE.EQ.2) THEN
C Second nucleosynthesis pass - this time for star 2
            CALL SOLVER ( 3, ITER, KTER, ERR, IE, NWRT5 )
C Restart if failed on composition
            IF (ERR.GT.EPS) nm = npr ! nm - 2
            IF (ERR.GT.EPS) WRITE (*,*) '!!COMPOSITION FAILED, RESTARTING!!'
            IF (ERR.GT.EPS) CALL PRINTA ( 2, NSTEP, ITER1, ITER2, NWRT5 ) 
            IF (ERR.GT.EPS) GO TO 1
         END IF
C If model didn't converge, give up
         IF (ERR.GT.EPS) WRITE (*,*) '!!FAILED TO CONVERGE, TERMINATING!!'
         IF (ERR.GT.EPS) GO TO 3
         IF (ERR.LT.EPS) CALL PRINTA ( 0, NSTEP, ITER1, ITER2, NWRT5 )
         
    1 ITER = ITER2
C Slightly insane bit where we compute dL/dt
      VLENEW = VLE
C Make sure we don't get a weird floating point error
      DLDT = (VLENEW - VLEOLD)/DTY/(VLENEW+(1.D-15))
C Grab the Helium core mass from the printb routines
      COREMASS = VMH      
C Compute the central degen param
      WF = DSQRT(1.0D0 + DEXP(H(1, NMESH)))
      PSICENTRAL = 2.0D0 * (WF - DLOG(WF + 1.0D0)) + H(1, NMESH)     
C Check to see if the model is going through the helium flash TODO
!      IF (((H(9, NMESH) + ZS).GT.(0.96D0)).AND.(PSICENTRAL.GT.6.D0).AND.
!     &   (DLDT.GT.1.D-4).AND.(COREMASS.GT.0).AND.(IHAVEFLASH.EQ.0)) THEN 
!        WRITE (*,*) 'ENTERING HELIUM FLASH'
!        IHAVEFLASH = 1
!        TARGETCOREMASS = COREMASS
!        CALL FORCEHEFLASH(TARGETCOREMASS, NMESH)
!        WRITE (*,*) 'EXITING HELIUM FLASH'
!      ENDIF  
      
      !IF (H(9,NMESH).LE.1.D-5) THEN
      !  WRITE(*,*) '!!!CORE HELIUM DEPLETED: ENTERING E-AGB: TERMINATING!!!'
      !  GOTO 3
      !ENDIF
!      IF (VMH.GE.0.2D0 .AND. iexhaust.EQ.0) THEN
!      iexhaust = 1
!      WRITE (*,*) '!!! 0.2 Solar Mass Exhausted Core Reached !!!'
!      !itzo_yn = 1
!      GOTO 3
!      ENDIF
       

      
      GOTO 2
C Output the last converged model. 
    3 CALL PRINTA ( 1, NSTEP, ITER1, ITER2, NWRT4 )
      CALL CPU_TIME(FINISH)
      WRITE (*,*) 'COMPUTATION TERMINATED'
      WRITE (*,*) FINISH-START, ' CPU SECONDS ELAPSED'
      STOP
      END
      
      
      
