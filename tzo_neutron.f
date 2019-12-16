**==.FOR     
      SUBROUTINE TZO_NEUTRON(rho,T,np,nn,conde,kappa,en)
*      IMPLICIT NONE
      COMMON /CNSTS / CPI, PI4, CLN10, CA, CB, CL, CD, CG, CR, CT, CEVB, 
     &     CEN, CPL, CMEVMU, CSECYR, LSUN, MSUN, RSUN, TSUNYR
     
      REAL*8 rho,T,kappa,en
      INTEGER ITMAX
      PARAMETER(ITMAX=100)

c Constants
      REAL*8 me,mn,Q,lambdan,lambdae,n0,tp2,third
      DATA me,mn,Q,lambdan,lambdae,n0,tp2,third /
     &     9.10938188d-28,       ! Electon mass (g)
     &     1.67492716d-24,       ! Neutron mass (g)
     &     2.30558d-27,          ! Q=m_n-m_p (g)
     &     2.1001942d-14,        ! lambda_n = hbar/(m_n c) (cm)
     &     3.86159276e-11,       ! lambda_e = hbar/(m_e c) (cm)
     &     1.6e+38,              ! std. num. dens. n0 = .16 fm^-3 
     &     19.739209d0,          ! 2 pi^2
     &     .3333333333333333d0/  ! 1/3

      INTEGER i
      REAL*8 nn,np,ne,ncrit,rpn,rpn0 ! Neut., prot. & elect. num. dens., critical nn, ratio
      REAL*8 chin,chin2,chip,chie    ! Chemical potentials of n(+sq),p,e
      REAL*8 mnstar,mpstar,mestar    ! Effective masses
      REAL*8 T8,conde
      real*8 alphan,betan,Qp0,Qn0

      REAL*8 mp

c If rho < the critical density (eq. 2.5.12 of Shapiro & Teukolski 83 [ST83]) no neutrons
      ncrit = ((Q/me)**2-1.d0)**1.5/(tp2*lambdae**3)
      IF (rho.LT.ncrit*mn) THEN
         np = rho/mn
         ne = np
         nn = 0.d0
         conde = 1.d0
         kappa = 1.d0
         Qn0 = 0.d0
         RETURN
      END IF

c Solve iteratively for the proton and neutron number densities
      rpn = 1.d0/8.d0
      rpn0 = 1.d0
      i = 0
      DO WHILE (dabs(rpn-rpn0)/(rpn+rpn0) .GT. 1.d-3 .AND. i.LT.itmax)
         rpn0 = rpn
         nn = rho/(mn*(1+rpn))
         chin = (tp2*nn)**third*lambdan
         chin2 = chin*chin
         rpn = .125d0*((1.d0 + 4.d0*Q/(mn*chin2) + 4.d0*(Q**2-me**2)/(mn**2*chin**4))/
     &                 (1.d0+1.d0/chin2))**1.5        ! ST83 eq. 2.5.17
         i = i + 1
      END DO
      np = nn * rpn
      ne = np

c Calculate effective masses (=chemical potentials)
      mp = mn - Q
      chip = chin * mn / mp * rpn**third       ! ST83 eq. 2.5.16
      chie = mp / me * chip                    ! ST83 eq. 2.5.7
      mnstar = mn * sqrt(1.d0 + chin2)
      mpstar = mp * sqrt(1.d0 + chip**2)
      mestar = me * sqrt(1.d0 + chie**2)

c nue is the electron collision frequency, expressed as a sum over species
      

c Electron conductivities when scattering off protons dominates: Eq. 77 of Gnedin & Yakovlev (95) Nu.Ph.A 582 697
      T8 = T/1.e8
      conde = 6.3d24 * 1.2d0 / T8 * sqrt(mp/mpstar) * (np/n0)**(7.d0/6.d0)
      kappa = 4.d0*CA*CL*T**3/(3.d0*rho*conde)
c      conde = 1.7d24 * 1.2d0 * T8 * (np/n0)**(2.d0/3.d0) / nue

c Neutrino loss rates in erg/(cm^3 s) from Yakovlev & Levenfish (1995) A&A 297 717 (ignoring surpression)
      RETURN
      alphan = max(0.d0, 1.76d0 - 0.63d0 * (n0/nn)**(2.d0/3.d0))
      betan = 0.68d0
      Qn0 = 8.55d21 * (mnstar/mn)**3 * (mpstar/mp) * (ne/n0)**third * T8**9 * alphan * betan
      if (chin*mn .lt. 3.d0*chip*mp + chie*me) then
         Qp0 = Qn0 * 0.75d0 * (mpstar/mnstar)**2
      else
         Qp0 = 0.d0
      end if
      EN = (Qn0 + Qp0)/rho
      END