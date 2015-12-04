C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C& FRHM               
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      SUBROUTINE FRHM(ID,ISTEP)
        COMPLEX FFFF,Y1,FFUU,FFWW,GL(2,2),P0,P1,P2
        COMPLEX S0,S1,S2,F1,F2,P3,S3
        COMPLEX TTTT,AICON,ttttp,tttts
        complex G111P,G111S,G112P,G112S
        complex G121P,G121S,G122P,G122S
        complex G221P,G221S,G222P,G222S
        REAL  M11,M12,M21,M22
        REAL DR(2),D(2,2)
C
        INCLUDE 'bemlayer.fin'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        DE=RMATE(ID,1)
        VP=RMATE(ID,2)
        VS=RMATE(ID,3)
        Y1=-(0.0,1.0)
        AICON=-CMPLX(0.0,0.25)
        FFFF=-CMPLX(REPA(ISTEP),AIMPA(ISTEP))
        W=2.0*PI*ISTEP/(NSAMP*DELTA)
        WP=W/VP
        WS=W/VS
        C1=1.0/(DE*VS*VS)
        C2=1.0/(DE*VP*VP)
        C3=VS/VP
        C4=C3*C3
        D(1,1)=1.0
        D(1,2)=0.0
        D(2,1)=0.0
        D(2,2)=1.0
        NBBP=NBDATA(ID,1)
C------------------------------
      DO IFS=1,NFS
      DO 300 INODE=1,NBBP
       IF(SFLAG(IFS).EQ.0) THEN
C-----------------------------------
C PLANE WAVE
C-----------------------------------
        FFUU=FFFF*SIN(THETA(IFS))
        FFWW=FFFF*COS(THETA(IFS))
        XF=COORD(1,INODE,1)
        YF=COORD(1,INODE,2)
        TTTT=cexp(-Y1*WP*(XF*SIN(THETA(IFS))-YF*COS(THETA(IFS))))/AICON
        IO=2*INODE-1
        IE=IO+1
        RESU(IO,IFS)=FFUU*TTTT
        RESU(IE,IFS)=-FFWW*TTTT
       ELSEIF(SFLAG(IFS).EQ.1) THEN
C-----------------------------------
C PLANE WAVE
C-----------------------------------
        FFUU=FFFF*SIN(THETA(IFS))
        FFWW=FFFF*COS(THETA(IFS))
        XF=COORD(1,INODE,1)
        YF=COORD(1,INODE,2)
        TTTT=exp(-Y1*WS*(XF*SIN(THETA(IFS))-YF*COS(THETA(IFS))))/AICON
        IO=2*INODE-1
        IE=IO+1
        RESU(IO,IFS)=-Y1*WS*FFFF*TTTT*cos(theta(IFS))
        RESU(IE,IFS)=-Y1*WS*FFFF*TTTT*sin(theta(IFS))
       ELSEIF(SFLAG(IFS).EQ.2) THEN
C-----------------------------------
CSINGELE FORCE
C-----------------------------------
        FFUU=FFFF*SIN(THETA(IFS))
        FFWW=FFFF*COS(THETA(IFS))
        XF=COORD(1,INODE,1)
        YF=COORD(1,INODE,2)
        XMX=XF-XSEC(IFS) 
        YMY=YF-YSEC(IFS)
        R=SQRT(XMX*XMX+YMY*YMY)
        IF(R.LT.0.1)R=0.1
        DR(1)=XMX/R
        DR(2)=YMY/R
        ZP=WP*R
        ZS=WS*R
        CALL HANKEX(ZP,P0,P1)
        CALL HANKEX(ZS,S0,S1)
        P2=2.0*P1/ZP-P0   
        S2=2.0*S1/ZS-S0
        F1=S0+(C3*P1-S1)/ZS
        F2=C4*P2-S2
        DO  K=1,2
        DO  L=1,2
        GL(K,L)=C1*(F1*D(K,L)-F2*DR(K)*DR(L))
        ENDDO
        ENDDO
        IO=2*INODE-1
        IE=IO+1
        RESU(IO,IFS)=GL(1,1)*FFUU+GL(1,2)*FFWW
        RESU(IE,IFS)=GL(2,1)*FFUU+GL(2,2)*FFWW
       ELSEIF(SFLAG(IFS).EQ.3) THEN
C-----------------------------------
C EXPLOSIVE 
C-----------------------------------
        XF=COORD(ID,INODE,1)
        YF=COORD(ID,INODE,2)
        XMX=XF-XSEC(IFS)
        YMY=YF-YSEC(IFS)
        R=SQRT(XMX*XMX+YMY*YMY)
        IF(R.LT.0.1)R=0.1
        DR(1)=XMX/R
        DR(2)=YMY/R
        ZP=WP*R
        CALL HANKEX(ZP,P0,P1)
        IO=2*INODE-1
        IE=IO+1
        RESU(IO,IFS)=-FFFF*P1*WP*DR(1)*C2
        RESU(IE,IFS)=-FFFF*P1*WP*DR(2)*C2
       ELSEIF(SFLAG(IFS).EQ.4) THEN
C-------------------------------------
C DOUBLE COUPLE
C-------------------------------------
            M11=-sin(theta1(IFS))*cos(theta2(IFS))*sin(2*theta(IFS))
     &          -sin(2*theta1(IFS))*sin(theta2(IFS))*sin(theta(IFS))**2
            M12=-cos(theta1(IFS))*cos(theta2(IFS))*cos(theta(IFS))
     &          -cos(2*theta1(IFS))*sin(theta2(IFS))*sin(theta(IFS))
            M21=M12
            M22=sin(2*theta1(IFS))*sin(theta2(IFS))
            XF=COORD(1,INODE,1)
            YF=COORD(1,INODE,2)
            XMX=XF-XSEC(IFS)
            YMY=YF-YSEC(IFS)
            R=SQRT(XMX*XMX+YMY*YMY)
            IF(R.LT.0.1)R=0.1
            DR(1)=XMX/R
            DR(2)=YMY/R
            ZP=WP*R
            ZS=WS*R
            CALL HANKEX(ZP,P0,P1)
            CALL HANKEX(ZS,S0,S1)
            P2=2.0*P1/ZP-P0
            S2=2.0*S1/ZS-S0
            P3=4.0*P2/ZP-P1
            S3=4.0*S2/ZS-S1

            G111P=FFFF*C2*((-P1+(P1-P3)/2*(1-2*DR(1)*DR(1)))*WP
     &           -4.0*DR(2)*DR(2)*P2/R)*DR(1)/2.0
            G111S=FFFF*C1*((-S1-(S1-S3)/2*(1-2*DR(1)*DR(1)))*WS
     &           +4.0*DR(2)*DR(2)*S2/R)*DR(1)/2.0
            G112P=FFFF*C2*((-P1+(P1-P3)/2*(1-2*DR(1)*DR(1)))*WP
     &           +4.0*DR(1)*DR(1)*P2/R)*DR(2)/2.0
            G112S=FFFF*C1*((-S1-(S1-S3)/2*(1-2*DR(1)*DR(1)))*WS
     &           -4.0*DR(1)*DR(1)*S2/R)*DR(2)/2.0

            G121P=FFFF*C2*(-DR(1)*DR(1)*(P1-P3)*WP
     &           +(4.0*DR(1)*DR(1)-2.0)*P2/R)*DR(2)/2.0
            G121S=FFFF*C1*((S1-S3)*WS*DR(1)*DR(1)
     &           +(2.0-4.0*DR(1)*DR(1))*S2/R)*DR(2)/2.0
            G122P=FFFF*C2*(-DR(2)*DR(2)*(P1-P3)*WP
     &           +(4.0*DR(2)*DR(2)-2.0)*P2/R)*DR(1)/2.0
            G122S=FFFF*C1*(DR(2)*DR(2)*(S1-S3)*WS
     &           +(2.0-4.0*DR(2)*DR(2))*S2/R)*DR(1)/2.0

            G221P=FFFF*C2*((-P1+(P1-P3)/2*(1-2*DR(2)*DR(2)))*WP
     &           +4.0*DR(2)*DR(2)*P2/R)*DR(1)/2.0
            G221S=FFFF*C1*((-S1-(S1-S3)/2*(1-2*DR(2)*DR(2)))*WS
     &           -4.0*DR(2)*DR(2)*S2/R)*DR(1)/2.0
            G222P=FFFF*C2*((-P1+(P1-P3)/2*(1-2*DR(2)*DR(2)))*WP
     &           -4.0*DR(1)*DR(1)*P2/R)*DR(2)/2.0
            G222S=FFFF*C1*((-S1-(S1-S3)/2*(1-2*DR(2)*DR(2)))*WS
     &           +4.0*DR(1)*DR(1)*S2/R)*DR(2)/2.0
            IO=2*INODE-1
            IE=IO+1
            RESU(IO,IFS)=M11*(G111P+G111S)+M12*(G112P+G112S)
     &                  +M21*(G121P+G121S)+M22*(G122P+G122S)
            RESU(IE,IFS)=M11*(G121P+G121S)+M12*(G122P+G122S)
     &                  +M21*(G221P+G221S)+M22*(G222P+G222S)

        ENDIF
300     CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-----------------------------
        IN1=NSTA(ID,1,1)
        AA=-NS
        DAMP=0.5/AA
        DO INODE=IN1,IN1+NS-1
         DECAY=COS(AA*DAMP*PI)**3
cc         print*,decay
         IO=2*INODE-1
         IE=IO+1
         RESU(IO,IFS)=RESU(IO,IFS)*DECAY
         RESU(IE,IFS)=RESU(IE,IFS)*DECAY
         AA=AA+1.0
        END DO
C-------------------------------------
        IN2=NSTA(ID,1,2)
        AA=-NS
        DO INODE=IN2,IN2-NS+1,-1
c          DECAY=EXP(AA*DAMP)
          DECAY=COS(AA*DAMP*PI)**3
          IO=2*INODE-1
          IE=IO+1
          RESU(IO,IFS)=RESU(IO,IFS)*DECAY
          RESU(IE,IFS)=RESU(IE,IFS)*DECAY
          AA=AA+1.0
        END DO
C-------------------------------------
        IF(ND.NE.ID) THEN
          IN1=NSTA(ID,2,1)
          AA=-NS
          DO INODE=IN1,IN1+NS-1
           IO=2*INODE-1
           IE=IO+1
           DECAY=EXP(AA*DAMP)
           RESU(IO,IFS)=RESU(IO,IFS)*DECAY
           RESU(IE,IFS)=RESU(IE,IFS)*DECAY
           AA=AA+1.0
          END DO
C-------------------------------------
        IN2=NSTA(ID,2,2)
        AA=-NS
        DO INODE=IN2,IN2-NS+1,-1
          IO=2*INODE-1
          IE=IO+1
          DECAY=EXP(AA*DAMP)
          RESU(IO,IFS)=RESU(IO,IFS)*DECAY
          RESU(IE,IFS)=RESU(IE,IFS)*DECAY
          AA=AA+1.0
        END DO
       ENDIF
ccccc enddo for ifs
      ENDDO 
ccccccccccccccccccccccccccccccccccccccccccccccc

C--------------------------------------            
      RETURN
      END
