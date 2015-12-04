C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C& SHLI    
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      SUBROUTINE SHLI(ISTEP,ID)
C  
      INCLUDE 'bemlayer.fin'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      VLO1=RMATE(ID,3) 
      VLO2=VLO1
      IF(ID.NE.ND)VLO2=RMATE(ID+1,3)
      WAFF=ISTEP/(NSAMP*DELTA)
      XLEN1=VLO1/(NGQ*WAFF)
      XLEN2=VLO2/(NGQ*WAFF)
      IF(XLEN2.GT.XLEN1)XLEN2=XLEN1
      IF(XLEN1.LT.XMIN)XLEN1=XMIN
      IF(XLEN2.LT.XMIN)XLEN2=XMIN
      IF(XLEN1.GT.XMAX)XLEN1=XMAX
      IF(XLEN2.GT.XMAX)XLEN2=XMAX
      NGP=3
      IF(XLEN1.GT.1000.)NGP=4
      IC=0
      JC=0
      KC=NS
C---------------------------------------------
      IF(ID.EQ.1)THEN
      J1=MSTA(ID,1,1)+1   
      J2=MSTA(ID,1,1)
      X1=XCOORD(ID,J1,1)
      Y1=XCOORD(ID,J1,2)
      X2=XCOORD(ID,J2,1)
      Y2=XCOORD(ID,J2,2)
      RLEN=SQRT((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2))
      BM=KC*XLEN1
      AMD=-(BM+RLEN)/BM
      X0=(X1+AMD*X2)/(1.0+AMD)
      Y0=(Y1+AMD*Y2)/(1.0+AMD)
      IC=IC+1
      COORD(ID,IC,1)=X0
      COORD(ID,IC,2)=Y0
      COORD(ID,IC,3)=0.0
      COORD(ID,IC,4)=PI
      DO INL=1,KC-1
      IC=IC+1
      RATO=FLOAT(INL)/(FLOAT(KC)-FLOAT(INL))
      COORD(ID,IC,1)=(X0+RATO*X2)/(1.0+RATO)
      COORD(ID,IC,2)=(Y0+RATO*Y2)/(1.0+RATO)
      COORD(ID,IC,3)=0.0
      COORD(ID,IC,4)=PI
      END DO
c
      DO 10 J=MSTA(ID,1,1),MSTA(ID,1,2)-1
      XFT=XCOORD(ID,J,1)
      YFT=XCOORD(ID,J,2)
      XSE=XCOORD(ID,J+1,1)
      YSE=XCOORD(ID,J+1,2)
      RLEN=SQRT((XFT-XSE)*(XFT-XSE)+(YFT-YSE)*(YFT-YSE))
      INCR=IFIX(RLEN/XLEN1+0.75)
      IC=IC+1
      COORD(ID,IC,1)=XFT
      COORD(ID,IC,2)=YFT
      COORD(ID,IC,3)=BLPH(ID,J,1)
      COORD(ID,IC,4)=BLPH(ID,J,2)
      IF(INCR.GT.1)THEN
        DO INL=1,INCR-1
         IC=IC+1
         RATO=FLOAT(INL)/(FLOAT(INCR)-FLOAT(INL))
         COORD(ID,IC,1)=(XFT+RATO*XSE)/(1.0+RATO)
         COORD(ID,IC,2)=(YFT+RATO*YSE)/(1.0+RATO)
         COORD(ID,IC,3)=COORD(ID,IC-1,4)-PI
         COORD(ID,IC,4)=COORD(ID,IC-1,4)
        END DO
      ENDIF
 10   CONTINUE
c
      J1=MSTA(ID,1,2)-1   
      J2=MSTA(ID,1,2)
      X1=XCOORD(ID,J1,1)
      Y1=XCOORD(ID,J1,2)
      X2=XCOORD(ID,J2,1)
      Y2=XCOORD(ID,J2,2)
      RLEN=SQRT((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2))
      BM=KC*XLEN1
      AMD=-(BM+RLEN)/BM
      X3=(X1+AMD*X2)/(1.0+AMD)
      Y3=(Y1+AMD*Y2)/(1.0+AMD)
      IC=IC+1
      COORD(ID,IC,1)=X2
      COORD(ID,IC,2)=Y2
      COORD(ID,IC,3)=COORD(ID,IC-1,4)-PI
      COORD(ID,IC,4)=COORD(ID,IC-1,4)
      DO INL=1,KC-1   
      IC=IC+1
      RATO=FLOAT(INL)/(FLOAT(KC)-FLOAT(INL))
      COORD(ID,IC,1)=(X2+RATO*X3)/(1.0+RATO)
      COORD(ID,IC,2)=(Y2+RATO*Y3)/(1.0+RATO)
      COORD(ID,IC,3)=COORD(ID,IC-1,4)-PI
      COORD(ID,IC,4)=COORD(ID,IC-1,4)
      END DO
c
      IC=IC+1
      COORD(ID,IC,1)=X3
      COORD(ID,IC,2)=Y3
      COORD(ID,IC,3)=COORD(ID,IC-1,4)-PI
      COORD(ID,IC,4)=COORD(ID,IC-1,4)
      NSTA(ID,1,1)=1
      NSTA(ID,1,2)=IC
c
      DO I=1,IC-1
      JC=JC+1
      NEL(ID,JC,1)=I
      NEL(ID,JC,2)=I+1
      END DO
c
       ELSE
      NSTA(ID,1,1)=1
      N1=NSTA(ID-1,2,1)
      N2=NSTA(ID-1,2,2)
      DO INODE=N2,N1,-1
      IC=IC+1
      COORD(ID,IC,1)=COORD(ID-1,INODE,1)
      COORD(ID,IC,2)=COORD(ID-1,INODE,2)
      COORD(ID,IC,3)=COORD(ID-1,INODE,4)-2.0*PI
      COORD(ID,IC,4)=COORD(ID-1,INODE,3)
      END DO
      NSTA(ID,1,2)=IC
      DO I=1,IC-1
      JC=JC+1
      NEL(ID,JC,1)=I
      NEL(ID,JC,2)=I+1
      END DO
       ENDIF
C-----------------------------------------------
       IF(ID.EQ.ND)THEN
      NSTA(ID,2,1)=0
      NSTA(ID,2,2)=-1
       ELSE
      NSTA(ID,2,1)=IC+1
      J1=MSTA(ID,2,1)+1    
      J2=MSTA(ID,2,1)
      X1=XCOORD(ID,J1,1)
      Y1=XCOORD(ID,J1,2)
      X2=XCOORD(ID,J2,1)
      Y2=XCOORD(ID,J2,2)
      RLEN=SQRT((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2))
      BM=KC*XLEN2
      AMD=-(BM+RLEN)/BM
      X0=(X1+AMD*X2)/(1.0+AMD)
      Y0=(Y1+AMD*Y2)/(1.0+AMD)
      IC=IC+1
      COORD(ID,IC,1)=X0
      COORD(ID,IC,2)=Y0
      COORD(ID,IC,3)=PI
      COORD(ID,IC,4)=2.*PI
      DO INL=1,KC-1  
      IC=IC+1
      RATO=FLOAT(INL)/(FLOAT(KC)-FLOAT(INL))
      COORD(ID,IC,1)=(X0+RATO*X2)/(1.0+RATO)
      COORD(ID,IC,2)=(Y0+RATO*Y2)/(1.0+RATO)
      COORD(ID,IC,3)=PI
      COORD(ID,IC,4)=2.*PI
        END DO
C
      DO 20 J=MSTA(ID,2,1),MSTA(ID,2,2)-1
      XFT=XCOORD(ID,J,1)
      YFT=XCOORD(ID,J,2)
      XSE=XCOORD(ID,J+1,1)
      YSE=XCOORD(ID,J+1,2)
      RLEN=SQRT((XFT-XSE)*(XFT-XSE)+(YFT-YSE)*(YFT-YSE))
      INCR=IFIX(RLEN/XLEN2+0.75)
      IC=IC+1
      COORD(ID,IC,1)=XFT
      COORD(ID,IC,2)=YFT
      COORD(ID,IC,3)=BLPH(ID,J,1)
      COORD(ID,IC,4)=BLPH(ID,J,2)
      IF(INCR.GT.1)THEN
      DO INL=1,INCR-1
      RATO=FLOAT(INL)/(FLOAT(INCR)-FLOAT(INL))
      IC=IC+1
      COORD(ID,IC,1)=(XFT+RATO*XSE)/(1.0+RATO)
      COORD(ID,IC,2)=(YFT+RATO*YSE)/(1.0+RATO)
      COORD(ID,IC,3)=COORD(ID,IC-1,4)-PI
      COORD(ID,IC,4)=COORD(ID,IC-1,4)
        END DO
      ENDIF
 20     CONTINUE
      IC=IC+1
      COORD(ID,IC,1)=XSE
      COORD(ID,IC,2)=YSE
      COORD(ID,IC,3)=PI
      COORD(ID,IC,4)=2.0*PI
C
      J1=MSTA(ID,2,2)-1  !RIGHT POINT
      J2=MSTA(ID,2,2)
      X1=XCOORD(ID,J1,1)
      Y1=XCOORD(ID,J1,2)
      X2=XCOORD(ID,J2,1)
      Y2=XCOORD(ID,J2,2)
      RLEN=SQRT((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2))
      BM=KC*XLEN2
      AMD=-(BM+RLEN)/BM
      X3=(X1+AMD*X2)/(1.0+AMD)
      Y3=(Y1+AMD*Y2)/(1.0+AMD)
      DO INL=1,KC-1  
      IC=IC+1
      RATO=FLOAT(INL)/(FLOAT(KC)-FLOAT(INL))
      COORD(ID,IC,1)=(X2+RATO*X3)/(1.0+RATO)
      COORD(ID,IC,2)=(Y2+RATO*Y3)/(1.0+RATO)
      COORD(ID,IC,3)=PI
      COORD(ID,IC,4)=2.0*PI
      END DO
      IC=IC+1
      COORD(ID,IC,1)=X3
      COORD(ID,IC,2)=Y3
      COORD(ID,IC,3)=PI
      COORD(ID,IC,4)=2.0*PI
C
      NSTA(ID,2,2)=IC
      DO I=NSTA(ID,2,1),NSTA(ID,2,2)-1
      JC=JC+1
      NEL(ID,JC,1)=I
      NEL(ID,JC,2)=I+1
      END DO
      ENDIF
C-------------------------------------
      NBDATA(ID,1)=IC
      NBDATA(ID,2)=JC
C---------------------------------------
      RETURN
      END
