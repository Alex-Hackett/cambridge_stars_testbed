C Between the comments added to the code and the file 'bswriteup.tex'
C you can hopefully understand how this all works... RJS 9/6/06
C PPEs original comment:
C The working of this programme is described in a file `writeup.tex';
C Consequently I do not include many comments in the body of the programme.
C Following is the main routine
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE
      PARAMETER (MAXMSH = 2000)
      REAL MSUM
      COMMON H(60,MAXMSH),DH(60,MAXMSH),EPS,V(2),NMESH,JIN,ID(100),IE(100)
      COMMON /EVMODE/ IMODE
      COMMON /SODDS / ALPHA, RML, CMG, CSI, CFE, CT(10), AGE, DT, M1, 
     :  EC, BM, ANG, CM, MTA, MTB, TM(2), T0, M0, TC(2), OS, AC, RCD,
     :  RMG, RHL, XF, DR, AK1 ,RMT, AK2, IZ(4), IB, ISX(45),
     :  TRB
      COMMON /MASTRK/ AMASSSH(MAXMSH), YSHELLS(MAXMSH), SCALEMASS(MAXMSH)
      COMMON /ARTENG/ TENG, SMASS, FMASS, SHIENG, INJMD
      real dtime,cpu(2),dt,tcpu
C Read physical data and an initial model
      CALL PRINTA ( -1, NSTEP, ITER1, ITER2, NWRT5 )
      IF (NSTEP.EQ.0) GO TO 3
      ITER = ITER1
      nter = 0
      nm = 0
      npr = nm
      tcpu = dtime(cpu)
      tcpu = 0.0
C Begin evolutionary loop of NSTEP time steps
    2 IF (NM .EQ. NSTEP) GOTO 3
C Have we reached the set age?
*      IF (AGE.GT.(1.27038904*10)**10) THEN
*        WRITE (*,*) 'Set Age Reached'
*        GOTO 3
*      ENDIF
         IF (AGE.LT.10**6) WRITE (*,*) 'Solving Model: ', NM, '     Model Age: ', AGE, 'Yrs'
         IF ((AGE.LT.10**9).AND.(AGE.GT.10**6)) WRITE (*,*) 'Solving Model: ', NM, '     Model Age: ', AGE/(10**6), 'Myr'
         IF (AGE.GT.10**9) WRITE (*,*) 'Solving Model: ', NM, '     Model Age: ', AGE/(10**9), 'Gyr'
         
C Solve for structure, mesh, and major composition variables
         CALL SOLVER ( 1, ITER, KTER, ERR, ID, NWRT5 )
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
      GOTO 2
C Output the last converged model. 
    3 CALL PRINTA ( 1, NSTEP, ITER1, ITER2, NWRT4 )
      STOP
      END
