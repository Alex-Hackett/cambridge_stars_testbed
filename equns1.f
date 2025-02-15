**==EQUNS1.FOR
      SUBROUTINE EQUNS1(K, K1, K2, ISTAR, IVAR)
      IMPLICIT REAL*8(A-H, L, M, O-Z)
      PARAMETER (NMAXMSH = 2000)
      PARAMETER (DPI=3.141592653589793238D0)
      COMMON /INE   / BC1(3), BC2(3), BC3(3), BCHORB(3), BCHSPIN(3),
     :            VP(3), VPK(3), R2(3), 
     :            R2K(3), VT(3), VTK(3), L(3), LK(3), LQ(3), GTA(3), 
     :            MT(3), VM(3), VMK(3), QK(3), SG(3), T(3), X1(3),  
     :            X1T(3), X16(3), X16T(3), X4(3), X4T(3), X12(3),  
     :            X12T(3), X20(3), X20T(3), X14(3), X14T(3), X24(3),
     :            SGTH(3), MU(3), X3(3), X3T(3), SGLEV(3), DA4(3),
     :            DA12(3), DA14(3), DA16(3), DA20(3), DA3(3),
     :            D4(3), D12(3), D14(3), D16(3), D20(3), D3(3),
     :            WI(309)
      COMMON /OP    / ZS, LEDD, VVM, GR, GRAD, ETH, RLF, EGR, R, QQ
C VAR(3),(2),(1) are values of VAR at current, previous and anteprevious meshpts
      COMMON /OUTE  / EQU(50)
      COMMON /MESH  / TRC1,TRC2,DD,DT1,DT2,MWT,MWTS, IVMC,IVMS
      COMMON /ACCRET/ ACCOMPOS(7,31, 2)
      COMMON /TRANS / HT(26,NMAXMSH,2)
      COMMON /EVMODE/ IMODE
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
     :          rtzo_ct_3, rtzo_ct_3_per_yr, rtzo_ct_3_max,
     :          rtzo_core_mass
     
C Common block for Quasi-star stuff
      COMMON /QSM/ QSM_COREMASS, QSM_CORERAD, QSM_CORELUM, QSM_DT
      COMMON /IQSM/ IQSM_FLAG
      COMMON /STAT2 / PL, RL, U, P, RHO, FK, TEMPERATURE, SF, ST, ZT, GRADA, CP, 
     :                CH, S, PR, PG, PF, PT, EN, WR(41)
      COMMON /CNSTS / CPI, PI4, CLN10, CDUM(11), CSECYR, LSUN, MSUN,
     &                RSUN, TSUNYR
      COMMON H(60,NMAXMSH),DH(60,NMAXMSH),EPS,V(2),NMESH,JIN,ID(100),IE(100) 
      COMMON /QSM_AGE/ Q_AGE
      PS(VX) = 0.5D0*(VX+DABS(VX))
      DT = QSM_DT
C 30/5/03 RJS Smooth viscous mesh
      WTM = 0.5 + 0.5*tanh((K - TRC1)/1.5)
      WTM = MWT*WTM
C Surface mesh viscosity
      IF (IVMS.EQ.1) THEN
         WTM = WTM + MWTS*(0.5 - 0.5*tanh((K - TRC2)/1.5))
      END IF
      IF (WTM.GT.1.0) WTM = 1.0
      IF (WTM.LT.0.0) WTM = 0.0
C     TZO mesh fluidity stuff!      
      IF (itzo_yn.EQ.1) THEN
        WTM = rtzo_meshfluid
      ENDIF
      
      IF ( K.LE.K1 ) THEN
C surface boundary conditions
         EQU(1) = BC1(3)
         EQU(2) = BC2(3)
         EQU(3) = BC3(3)
C Orbital angular momentum
         EQU(4) = BCHORB(3)
C Spin period of star
         EQU(5) = BCHSPIN(3)
         RETURN
      ELSE IF ( K.LE.K2 ) THEN
C first-order difference equations at interior points
         WT3 = 0.5D0
C weighted alternative to central differencing
C        WT3 = 0.5D0*WT(3)/(1.0D0+WT(3))
         WT2 = 1.0D0 - WT3
         EQU(1) = VP(3) - VP(2) - WT3*VPK(3)-WT2*VPK(2)
         EQU(2) = R2(3) - R2(2) - WT2*R2K(3)-WT3*R2K(2)
         EQU(3) = VT(3) - VT(2) - WT3*VTK(3)-WT2*VTK(2)
         EQU(4) = L(3) - L(2) - WT2*LK(3)-WT3*LK(2) 
     :              - LQ(2)*GTA(2)*PS(MT(2)) + LQ(3)*GTA(3)*PS(-MT(3))
         EQU(5) = VM(3) - VM(2) - 0.5D0*(VMK(3)+VMK(2))
C 22/3/03 RJS Added viscous mesh
         EQU(6)=(1.0-WTM)*(QK(3) - QK(2))+3.0d7*WTM*MT(3)
         EQU(13) = BCHSPIN(3) - BCHSPIN(2)
         IF ( K.EQ.K1+1 ) THEN
C next-to-surface boundary conditions for second-order equations
C Attempt at variable composition accretion - only if in binary mode
            IF (ISTAR.EQ.1) IOTHER = 2
            IF (ISTAR.EQ.2) IOTHER = 1
            IF ((HT(24,1,IOTHER).GT.0d0.OR.HT(23,1,ISTAR).LT.0d0).AND.IMODE.EQ.2) THEN
               EQU(7) = ACCOMPOS(1,IVAR+1,IOTHER) - X1(2)
               EQU(8) = ACCOMPOS(5,IVAR+1,IOTHER) - X16(2)
               EQU(9) = ACCOMPOS(2,IVAR+1,IOTHER) - X4(2)
               EQU(10) = ACCOMPOS(3,IVAR+1,IOTHER) - X12(2)
               EQU(11) = ACCOMPOS(6,IVAR+1,IOTHER) - X20(2)
               EQU(12) = ACCOMPOS(4,IVAR+1,IOTHER) - X14(2)
               EQU(14) = ACCOMPOS(7,IVAR+1,IOTHER) - X3(2)
C The following lines are a more accurate but less stable implementation of the
C next to surface boundary condition -- should be used if you want to do thermohaline mixing
C               SG2 = -(PS(MT(2))+1d-5) !1d-5
C               EQU(7) = SG2*(X1(3)-X1(2)) + PS(MT(2))*(X1(2)-ACCOMPOS(1,IVAR+1,IOTHER))
C     :              - X1T(2)
C               EQU(8) = SG2*(X16(3)-X16(2)) + PS(MT(2))*(X16(2)-ACCOMPOS(5,IVAR+1,IOTHER))
C     :              - X16T(2)
C               EQU(9) = SG2*(X4(3)-X4(2)) + PS(MT(2))*(X4(2)-ACCOMPOS(2,IVAR+1,IOTHER))
C     :              - X4T(2)
C               EQU(10) = SG2*(X12(3)-X12(2)) + PS(MT(2))*(X12(2)-ACCOMPOS(3,IVAR+1,IOTHER))
C     :              - X12T(2)
C               EQU(11) = SG2*(X20(3)-X20(2)) + PS(MT(2))*(X20(2)-ACCOMPOS(6,IVAR+1,IOTHER))
C     :              - X20T(2)
C               EQU(12) = SG2*(X14(3)-X14(2)) + PS(MT(2))*(X14(2)-ACCOMPOS(4,IVAR+1,IOTHER))
C     :              - X14T(2)
C               EQU(14) = SG2*(X3(3)-X3(2)) + PS(MT(2))*(X3(2)-ACCOMPOS(7,IVAR+1,IOTHER))
C     :              - X3T(2)
            ELSE
               EQU(7) = X1(3) - X1(2)
               EQU(8) = X16(3) - X16(2)
               EQU(9) = X4(3) - X4(2)
               EQU(10) = X12(3) - X12(2)
               EQU(11) = X20(3) - X20(2)
               EQU(12) = X14(3) - X14(2)
               EQU(14) = X3(3) - X3(2)
            END IF
            RETURN
         ELSE
C second-order difference equations at interior points
            SG1 = 0.5D0*(SG(1)+SG(2)) - PS(MT(2))
            SG2 = 0.5D0*(SG(2)+SG(3)) - PS(-MT(3))
C Add in thermohaline mixing
            SG1 = SG1 + 0.5*(SGTH(1)+SGTH(2))*PS((MU(1)-MU(2)))
            SG2 = SG2 + 0.5*(SGTH(2)+SGTH(3))*PS((MU(2)-MU(3)))
C Note gravitational settling is hard-wired into having H as the dominant element!
            EQU(7) = (SG2 + 0.5*(DA4(2)+DA4(3)))*(X1(3)-X1(2))
     :           - (SG1+0.5*(DA4(1)+DA4(2)))*(X1(2)-X1(1)) - X1T(2)
     :           + D4(2)*X4(2) - D4(3)*X4(3)
     :           + D12(2)*X12(2) - D12(3)*X12(3)
     :           + D14(2)*X14(2) - D14(3)*X14(3)
     :           + D16(2)*X16(2) - D16(3)*X16(3)
     :           + D20(2)*X20(2) - D20(3)*X20(3)
     :           + D3(2)*X3(2) - D3(3)*X3(3)
            EQU(8) = (SG2 + 0.5*(DA16(2)+DA16(3)))*(X16(3)-X16(2))
     :           - (SG1+0.5*(DA16(1)+DA16(2)))*(X16(2)-X16(1)) - X16T(2)
     :           - D16(2)*X16(2) + D16(3)*X16(3)
            EQU(9) = (SG2 + 0.5*(DA4(2)+DA4(3)))*(X4(3)-X4(2))
     :           - (SG1+0.5*(DA4(1)+DA4(2)))*(X4(2)-X4(1)) - X4T(2)
     :           - D4(2)*X4(2) + D4(3)*X4(3)
            EQU(10) = (SG2 + 0.5*(DA12(2)+DA12(3)))*(X12(3)-X12(2))
     :           - (SG1+0.5*(DA12(1)+DA12(2)))*(X12(2)-X12(1)) - X12T(2)
     :           - D12(2)*X12(2) + D12(3)*X12(3)
            EQU(11)= (SG2 + 0.5*(DA20(2)+DA20(3)))*(X20(3)-X20(2))
     :           - (SG1+0.5*(DA20(1)+DA20(2)))*(X20(2)-X20(1)) - X20T(2)
     :           - D20(2)*X20(2) + D20(3)*X20(3)
            EQU(12)= (SG2 + 0.5*(DA14(2)+DA14(3)))*(X14(3)-X14(2))
     :           - (SG1+0.5*(DA14(1)+DA14(2)))*(X14(2)-X14(1)) - X14T(2)
     :           - D14(2)*X14(2) + D14(3)*X14(3)
            EQU(14)= (SG2 + 0.5*(DA3(2)+DA3(3)))*(X3(3)-X3(2))
     :           - (SG1+0.5*(DA3(1)+DA3(2)))*(X3(2)-X3(1)) - X3T(2)
     :           - D3(2)*X3(2) + D3(3)*X3(3)
            RETURN
         END IF
      END IF
C central boundary conditions for first-order equations
      IF (WTM.GT.1.0) WTM = 1.0
      IF (WTM.LT.0.0) WTM = 0.0
C     Normal BCs
      IF (IQSM_FLAG.EQ.0) THEN      
          IF (WTM.NE.0.0) THEN
             EQU(1) = MT(3)
          ELSE
             EQU(1) = VM(3) + 1.5D0*VMK(3) - 0.5D0*VMK(2)
          END IF
C This was the original central L boundary condition -- PPE said there
C was a good reason for it but he couldn't remember it.
C      EQU(2) = L(3) + 0.93333D0*LK(3) - 0.1885D0*LK(2) - LQ(3)*GTA(3)
C     :         *MT(3)
          EQU(2) = L(3) + 1.5D0*LK(3) - 0.5D0*LK(2) - LQ(3)*GTA(3)
     :              *MT(3)
          EQU(3) = R2(3) + 1.5D0*R2K(3) - 0.5D0*R2K(2)
          
          
      ENDIF
      
      ! Quasi-Star BCs
      IF (IQSM_FLAG.EQ.11) THEN
      ! Compute the Radius of the NS from a simple TOV EOS
        QSM_CORERAD = (QSM_COREMASS / (2.08D-6))**(-1.D0 / 3.D0)
      
      ! Compute the mass flux from Bondi-Hoyle approach
      
        SOUND_SPEED = SQRT((4.D0 / 3.D0) * (P/RHO))
        
        !ACCRETE_RATE = DPI * (1.D0 / SQRT(2.D0)) * (QSM_CORERAD**2)
      !:  !    * RHO * sound_speed
      ! New accretion rate from Ball 2012, Convective luminosity 
      ! Argument for cutting the max flux by a factor of 4cs^2 / gammalambdaepsilonprimec^2
      ! Awful G, I'm very sorry
      ! will get a cgs style g/s accretion rate from this, (hopefully)
      ! Untimately want one in Msun yrs^-1
        ODD_BIG_G = CG
        ! Fiducial run values for these from Warrick
        ETA = 0.1D0
        EPSILONPRIME = 0.1D0
        ADIABAT = 4.D0 / 3.D0
        ACCRETE_RATE_CGS = (16.D0 * DPI) * (ETA / (EPSILONPRIME * ADIABAT))
     :      * (((ODD_BIG_G * (QSM_COREMASS * 1.989D33))**2)/
     :      (sound_speed * (CL)**2)) * RHO 
      ! Convert accretion rate from g/s to Msun yr^-1
        ACCRETE_RATE = ACCRETE_RATE_CGS * (5.0279D-34) * CSECYR
        
        
      ! Compute the accretion luminosity from the accretion rate
        ! I'm unsure what Warrick was desribing exactly as his
        ! modified accretion luminosity, so I'm including both, for
        ! posterity's sake
        
        !Super simple black hole style accretion
        !QSM_CORELUM_CGS=EPSILONPRIME * ACCRETE_RATE_CGS * CL**2
        
        
        !Maximum convective flux based luminosity calculation from
        ! Warrick 2012, that I don't quite understand, but I feel like
        ! I should try this out
        QSM_CORELUM_CGS = ((4.D0)/(ADIABAT * (1.D0 / SQRT(2.D0))))
     :      * ACCRETE_RATE_CGS * (sound_speed**2)
        
      ! Convert to solar units
        QSM_CORELUM = QSM_CORELUM_CGS / 3.839D33
      ! Convert from solar input units to code units
        QSM_COREMASS_EGG = QSM_COREMASS * 1.989D0
        QSM_CORERAD_EGG = QSM_CORERAD * 0.696D0
        QSM_CORELUM_EGG = QSM_CORELUM * 3.844D0
        
        EQU(1) = (QSM_COREMASS_EGG)
        EQU(2) = QSM_CORELUM_EGG
        EQU(3) = (QSM_CORERAD_EGG)
        
      !Update the core mass based on the accretion rate
        QSM_COREMASS = QSM_COREMASS + (ACCRETE_RATE * (DT/CSECYR))
        
      ENDIF
      
      
      ! Try grow core radius again
      IF (IQSM_FLAG.EQ.1) THEN
      
      EQU(1) = QSM_COREMASS
      EQU(2) = QSM_CORELUM
      EQU(3) = QSM_CORERAD
      
      !EQU(1) = 0.388675252d0!2.695766905533690D-10
      !EQU(2) = 3.89051216d02!-4.82362145D-06
      !EQU(3) = 13.6056164d0!0.016404132d0
!      IF (QSM_COREMASS.EQ.0) THEN
!      EQU(1) = (2.695766905533690D-10)
!      EQU(2) = (01.82362145D-06)
!      EQU(3) = (0.016404132d0)
!      QSM_COREMASS = EQU(1)
!      QSM_CORELUM = EQU(2)
!      QSM_CORERAD = EQU(3)
!      ELSEIF (QSM_COREMASS.NE.0) THEN
!      EQU(1) = QSM_COREMASS
!      EQU(2) = QSM_CORELUM
!      EQU(3) = QSM_CORERAD
!     
!      IF (QSM_COREMASS.LT.0.388675252d0) QSM_COREMASS=QSM_COREMASS+(1.d-4)
!      IF (QSM_CORELUM.LT.3.89051216d02) QSM_CORELUM=QSM_CORELUM+(1.d-4)
!      IF (QSM_CORERAD.LT.13.6056164d0) QSM_CORERAD=QSM_CORERAD+(1.d-4)
!      
!      IF (QSM_COREMASS.GE.0.388675252d0) WRITE (*,*) '!!MASS BC REACHED!!'
!      IF (QSM_CORELUM.GE.3.89051216d02) WRITE (*,*) '!!LUMINOSITY BC REACHED!!'
!      IF (QSM_CORERAD.GE.13.6056164d0) WRITE (*,*) '!!RADIUS BC REACHED!!'
  
      !ENDIF
      
      
      !WRITE (*,*) 'EQU(1), EQU(2), EQU(3)', EQU(1), EQU(2), EQU(3)
      !EQU(1) = 0.D0
      !EQU(2) = 0.D0
      !EQU(3) = 0.D0
      ENDIF
      
      
*      ! Grow core radius slowly
*      IF (IQSM_FLAG.EQ.2) THEN
*      ! Start with normal BCs, slowly increase core radius
*        !Hack using the radius to, well, store the radius
*        IF (WTM.NE.0.0) THEN
*                 EQU(1) = MT(3)
*              ELSE
*                 EQU(1) = VM(3) + 1.5D0*VMK(3) - 0.5D0*VMK(2)
*              END IF
*              EQU(2) = L(3) + 1.5D0*LK(3) - 0.5D0*LK(2) - LQ(3)*GTA(3)
*     :                *MT(3)
*     
*        IF (QSM_CORERAD.EQ.0.D0) THEN
*              EQU(3) = R2(3) + 1.5D0*R2K(3) - 0.5D0*R2K(2)
*              QSM_CORERAD = EQU(3)
*        ELSE
*            QSM_CORELUM = 1.1D0
*            
*            EQU(3) = EQU(3) + ((QSM_CORELUM * EQU(3)) * (DT/CSECYR))
*            WRITE (*,*) 'EQU(3)', EQU(3)
*            QSM_CORERAD = EQU(3)
*        ENDIF
*      ENDIF
          
      
C central boundary conditions for second-order equations
      SG2 = 0.5D0*(SG(2)+SG(3))
C Plus thermohaline mixing
      SG2 = SG2 + 0.5*(SGTH(2)+SGTH(3))*PS((MU(2)-MU(3)))
      EQU(4) = SG2*(X1(3)-X1(2)) + X1T(3)
      EQU(5) = SG2*(X16(3)-X16(2)) + X16T(3)
      EQU(6) = SG2*(X4(3)-X4(2)) + X4T(3)
      EQU(7) = SG2*(X12(3)-X12(2)) + X12T(3)
      EQU(8) = SG2*(X20(3)-X20(2)) + X20T(3)
      EQU(9) = SG2*(X14(3)-X14(2)) + X14T(3)
      EQU(10) = SG2*(X3(3)-X3(2)) + X3T(3)
      RETURN
      END
