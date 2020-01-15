**==OPACTY.FOR
      SUBROUTINE OPACTY (JX,TF,FR,FKL,FKH)
      IMPLICIT REAL*8 (A-H,N-Z)
      PARAMETER (MT = 191,MR = 121)
      COMMON /STAT1 / CSX(10), CS(121,191,10), CNIU(60,41,2), W(18400)
      COMMON /STAT3 / F(4,4,191,121,10),TFM(191),FRM(121),
     :                FKLM(6),FKHM(6)
      SAVE
*
* Calculate a bicubic spline interpolation fit for the temperature
* and density opacity fit.  Do not stop if the input lies outside the
* array but rather use the previous result.
*
* Check that we are in the table.
*
      IF((TF.LT.TFM(1)).OR.(TF.GE.TFM(MT)).OR.
     &   (FR.LT.FRM(1)).OR.(FR.GE.FRM(MR))) THEN
         FKL = FKLO
         FKH = FKHO
         WRITE (6,100) TF,FR,FKL,FKH
!         WRITE (*,*) 'MIN T: ', TFM(1)
!         WRITE (*,*) 'MAX T: ', TFM(MT)
!         WRITE (*,*) 'MIN RHO: ', FRM(1)
!         WRITE (*,*) 'MAX RHO: ', FRM(MR)
         WRITE (333,100) TF,FR,FKL,FKH
!         ICHEATFLAG = 0
!         IF (TF.LT.TFM(1)) THEN
!            TF = TFM(1)
!            WRITE (*,*) '!!!SETTING TF TO TMIN!!!'
!            ICHEATFLAG = 1
!         ENDIF
!         IF (FR.LT.FRM(1)) THEN
!            FR = FRM(1)
!            WRITE (*,*) '!!!SETTING FR TO RMIN!!!'
!            ICHEATFLAG = 1
!         ENDIF
!         IF (ICHEATFLAG.EQ.1) GOTO 5
      ELSE
!    5    CONTINUE !AWFUL, AWFUL TMIN RHOMIN CHEAT
         FKLO = FKL
         FKHO = FKH
*
* Find interval in which target point lies.
*
C Debug, check for integer overflow
        
         I = 1 + (MT-1)*(TF-TFM(1))/(TFM(MT)-TFM(1))
         J = 1 + (MR-1)*(FR-FRM(1))/(FRM(MR)-FRM(1))
         !WRITE (*,*) MT, TF, TFM(1), TFM(MT)

         DT = TF-TFM(I)
         DR = FR-FRM(J)
*         WRITE (*,*) 'FRM(J), FR, DR', FRM(J), FR, DR         
*
* Evaluate the splines.
*
         FKL = F(1,1,I,J,JX) + DR*(F(1,2,I,J,JX)
     &   + DR*(F(1,3,I,J,JX) + DR*F(1,4,I,J,JX)))
     &   + DT*(F(2,1,I,J,JX) + DR*(F(2,2,I,J,JX)
     &   + DR*(F(2,3,I,J,JX) + DR*F(2,4,I,J,JX)))
     &   + DT*(F(3,1,I,J,JX) + DR*(F(3,2,I,J,JX)
     &   + DR*(F(3,3,I,J,JX) + DR*F(3,4,I,J,JX)))
     &   + DT*(F(4,1,I,J,JX) + DR*(F(4,2,I,J,JX)
     &   + DR*(F(4,3,I,J,JX) + DR*F(4,4,I,J,JX))))))
*
         FKH = F(1,1,I,J,JX+1) + DR*(F(1,2,I,J,JX+1)
     &   + DR*(F(1,3,I,J,JX+1) + DR*F(1,4,I,J,JX+1)))
     &   + DT*(F(2,1,I,J,JX+1) + DR*(F(2,2,I,J,JX+1)
     &   + DR*(F(2,3,I,J,JX+1) + DR*F(2,4,I,J,JX+1)))
     &   + DT*(F(3,1,I,J,JX+1) + DR*(F(3,2,I,J,JX+1)
     &   + DR*(F(3,3,I,J,JX+1) + DR*F(3,4,I,J,JX+1)))
     &   + DT*(F(4,1,I,J,JX+1) + DR*(F(4,2,I,J,JX+1)
     &   + DR*(F(4,3,I,J,JX+1) + DR*F(4,4,I,J,JX+1))))))
      ENDIF
  100 FORMAT ('OPACITY OUT OF RANGE',' TF FR FKL FKH ',4F9.4)
      RETURN
      END


