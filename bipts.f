C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C& BIPTS                
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      SUBROUTINE BIPTS(ISTEP,MAB)
      COMPLEX FFFF,Y1,FFUU,FFWW,AICON,SUM(2)
      COMPLEX P0,P1,P2,P3,S0,S1,S2,S3,F1,F2,GL(2,2)
      complex G111P,G111S,G112P,G112S
      complex G121P,G121S,G122P,G122S
      complex G221P,G221S,G222P,G222S
      real m11,m12,m21,m22
      REAL D(2,2),DR(2),YY1(2,2)
      complex, allocatable :: GM(:,:)
C
      INCLUDE 'bemlayer.fin'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        TOL=1.0E-4
        AICON=CMPLX(0.0,0.25)
        Y1=-CMPLX(0.0,1.0)
        D(1,1)=1.0
        D(1,2)=0.0
        D(2,1)=0.0
        D(2,2)=1.0
        YY1(1,1)=1.0
        YY1(1,2)=0.0
        YY1(2,1)=0.0
        YY1(2,2)=1.0
        FFFF=CMPLX(REPA(ISTEP),AIMPA(ISTEP))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcccccc
        DO IN=1,NG
        ID=RECFLAG(IN,1)
        IFLAG=RECFLAG(IN,2)
        DE=RMATE(ID,1)
        VP=RMATE(ID,2)
        VS=RMATE(ID,3)
        W=2.0*PI*ISTEP/(NSAMP*DELTA)
        WP=W/VP
        WS=W/VS
        C1=1./(DE*VS*VS)
        C2=1./(DE*VP*VP)
        C3=VS/VP
        C4=C3*C3
        C7=1./C4
        IF(IFLAG.EQ.0) THEN
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc When the receiver is on the freesurface, we directly read
cc the results from the boundary nodes.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          XIN=ROCOOR(IN,1)
          YIN=ROCOOR(IN,2)     
          II=0  
          DO IC=NSTA(ID,1,2),NSTA(ID,1,1),-1
           II=II+1 
           IF(abs(COORD(ID,IC,1)-xin).lt.tol
     &         .and.abs(coord(id,ic,2)-yin).lt.tol) THEN
            DO IFS=1,NFS
              IE=2*IN-1
              RHBA(IE,IFS)=DISP(1,2*II-1,ifs)
              RHBA(IE+1,IFS)=DISP(1,2*II,ifs)
            ENDDO
           ENDIF
          ENDDO
        ELSE
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C ONLY WHEN THE RECIEVER AND SOURCE ARE in the same region
c will we calculate the direct wave from source
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            IO=2*IN-1
            IE=IO+1
            XIN=ROCOOR(IN,1)
            YIN=ROCOOR(IN,2)
            DO IFS=1,NFS
            IF(ID.EQ.ISOURCE(IFS)) THEN     
            XMX=XIN-XSEC(IFS)
            YMY=YIN-YSEC(IFS)
            R=SQRT(XMX*XMX+YMY*YMY)
            IF(R.LT.0.0001)R=0.0001
            DR(1)=XMX/R
            DR(2)=YMY/R
            IF(SFLAG(IFS).EQ.1) THEN

            ELSEIF(SFLAG(IFS).EQ.2) THEN
C----------------------------------------
c SINGLE FORCE
C----------------------------------------
            ZP=WP*R
            ZS=WS*R
            FFUU=FFFF*SIN(THETA(IFS))
            FFWW=FFFF*COS(THETA(IFS))
            CALL HANKEX(ZP,P0,P1)
            CALL HANKEX(ZS,S0,S1)
            P2=2.0*P1/ZP-P0
            S2=2.0*S1/ZS-S0
            F1=S0+(C3*P1-S1)/ZS
            F2=C4*P2-S2
            DO K=1,2 
            DO L=1,2
            GL(K,L)=C1*(F1*D(K,L)-F2*DR(K)*DR(L))
            ENDDO
            ENDDO
            RHBA(IO,IFS)=(GL(1,1)*FFUU+GL(1,2)*FFWW)
            RHBA(IE,IFS)=(GL(2,1)*FFUU+GL(2,2)*FFWW)
            ELSEIF(SFLAG(IFS).EQ.3) THEN
C-----------------------------------------
C           EXPLOSIVE SOURCe
C-----------------------------------------
            ZP=WP*R
            CALL HANKEX(ZP,P0,P1)
            RHBA(IO,IFS)=-FFFF*DR(1)*WP*P1*C2
            RHBA(IE,IFS)=-FFFF*DR(2)*WP*P1*C2 

            ELSEIF(SFLAG(IFS).EQ.4) THEN
C-----------------------------------------
C           SOURCE with a moment
C-----------------------------------------
            ZP=WP*R
            ZS=WS*R
            M11=-sin(theta1(IFS))*cos(theta2(IFS))*sin(2*theta(IFS))
     &          -sin(2*theta1(IFS))*sin(theta2(IFS))*sin(theta(IFS))**2
            M12=-cos(theta1(IFS))*cos(theta2(IFS))*cos(theta(IFS))
     &          -cos(2*theta1(IFS))*sin(theta2(IFS))*sin(theta(IFS))
            M21=M12
            M22=sin(2*theta1(IFS))*sin(theta2(IFS))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
            RHBA(IO,IFS)=M11*(G111P+G111S)+M12*(G112P+G112S)
     &                  +M21*(G121P+G121S)+M22*(G122P+G122S)
            RHBA(IE,IFS)=M11*(G121P+G121S)+M12*(G122P+G122S)
     &                  +M21*(G221P+G221S)+M22*(G222P+G222S)
            ENDIF
cccccccccccc* endif SFLAG
           ELSE
              RHBA(IO,IFS)=(0.0,0.0)
              RHBA(IE,IFS)=(0.0,0.0)
           ENDIF
           ENDDO
cccccccccccc* enddo IFS
C-----------------------------------------------
           NM1=MSTA(ID,1,1)
           NM2=MSTA(ID,2,2)
           yl1=XCOORD(ID,NM1,2)
           yl2=XCOORD(ID,NM2,2)
           XIN=ROCOOR(IN,1)
           YIN=ROCOOR(IN,2)
           IF(IFLAG.EQ.2) then
                YY1(1,1)=0.5
                YY1(2,2)=0.5
           else
                YY1(1,1)=1.0
                YY1(2,2)=1.0
           endif
           CALL INTAL(IN,XIN,YIN,ID,W,WP,WS,C1,C3,C4,C7,D)
            I=0
            IF(ID.EQ.1) THEN
            DO  IC=NSTA(ID,1,2),NSTA(ID,1,1),-1
            ICA=2*IC-2
            I=I+1
            IX=2*I-2
            DO  IH=1,2
            DO  JH=1,2
               A(IH,IX+JH)=BMT(IH,ICA+JH)
            ENDDO
            ENDDO
            ENDDO
           IF(ND.NE.ID) THEN
            DO IC=NSTA(ID,2,1),NSTA(ID,2,2)
            ICA=2*IC-2
            I=I+1
            IX=2*I-2
            DO  IH=1,2
            DO  JH=1,2
              A(IH,IX+JH)=BMT(IH,ICA+JH)
            ENDDO
            ENDDO
            ENDDO
            DO  IC=NSTA(ID,2,1),NSTA(ID,2,2)
            ICA=2*IC-2
            I=I+1
            IX=2*I-2
            DO  IH=1,2
            DO  JH=1,2
              A(IH,IX+JH)=-AMT(IH,ICA+JH)
            ENDDO
            ENDDO
            ENDDO
           ENDIF
           ELSE
            DO  IC=NSTA(ID,1,2),NSTA(ID,1,1),-1
            ICA=2*IC-2
            I=I+1
            IX=2*I-2
            DO  IH=1,2
            DO  JH=1,2
               A(IH,IX+JH)=BMT(IH,ICA+JH)
            ENDDO
            ENDDO
            ENDDO
            DO  IC=NSTA(ID,1,2),NSTA(ID,1,1),-1
            ICA=2*IC-2
            I=I+1
            IX=2*I-2
            DO  IH=1,2
            DO  JH=1,2
               A(IH,IX+JH)=AMT(IH,ICA+JH)
            ENDDO
            ENDDO
            ENDDO
           IF(ND.NE.ID) THEN
            DO IC=NSTA(ID,2,1),NSTA(ID,2,2)
            ICA=2*IC-2
            I=I+1
            IX=2*I-2
            DO  IH=1,2
            DO  JH=1,2
              A(IH,IX+JH)=BMT(IH,ICA+JH)
            ENDDO
            ENDDO
            ENDDO
            DO  IC=NSTA(ID,2,1),NSTA(ID,2,2)
            ICA=2*IC-2
            I=I+1
            IX=2*I-2
            DO  IH=1,2
            DO  JH=1,2
              A(IH,IX+JH)=-AMT(IH,ICA+JH)
            ENDDO
            ENDDO
            ENDDO
           ENDIF
           ENDIF 
           NA=I
C-------------------------------
          N1=2*(NSTA(ID,1,2)-NSTA(ID,1,1)+1)
          N2=2*(NSTA(ID,2,2)-NSTA(ID,2,1)+1)
          IF(ID.EQ.1.and.ID.EQ.ND) THEN
            ALLOCATE(GM(N1,NFS))
            GM(1:N1,1:NFS)=DISP(id,1:N1,1:NFS)
          ELSEIF(ID.EQ.1) THEN
            ALLOCATE(GM(N1+2*N2,NFS))
            GM(1:N1,1:NFS)=DISP(ID,1:N1,1:NFS)
            GM(N1+1:N1+N2,1:NFS)=DISP(ID+1,1:N2,1:NFS)
            GM(N1+N2+1:N1+2*N2,1:NFS)=TRACT(ID+1,1:N2,1:NFS)
          ELSEIF(ID.EQ.ND) THEN
            ALLOCATE(GM(2*N1,NFS))
            GM(1:N1,1:NFS)=DISP(ID,1:N1,1:NFS)
            GM(N1+1:2*N1,1:NFS)=TRACT(ID,1:N2,1:NFS)
          ELSE
            ALLOCATE(GM(2*N1+2*N2,NFS))
            GM(1:N1,1:NFS)=DISP(ID,1:N1,1:NFS)
            GM(N1+1:2*N1,1:NFS)=TRACT(ID,1:N1,1:NFS)
            GM(2*N1+1:2*N1+N2,1:NFS)=DISP(ID+1,1:N2,1:NFS)
            GM(2*N1+N2+1:2*N1+2*N2,1:NFS)=TRACT(ID+1,1:N2,1:NFS)
          ENDIF
          DO IFS=1,NFS
            SUM(1)=CMPLX(0.0,0.0)
            SUM(2)=CMPLX(0.0,0.0)
            DO I=1,2
             DO K=1,2*NA
                SUM(I)=SUM(I)+A(I,K)*GM(K,IFS)
             END DO
             IE=2*IN+I-2
             RHBA(IE,IFS)=RHBA(IE,IFS)+SUM(I)
cc             RHBA(IE,IFS)=RHBA(IE,IFS)
             RHBA(IE,IFS)=RHBA(IE,IFS)*AICON/YY1(I,I)
            END DO
          END DO
          deallocate(GM)
c-------------------------------------------------------
       ENDIF
      ENDDO
      DO IFS=1,NFS
         M1=(MAB-1)*NFS+1
         M2=MAB*NFS
         IG=0
         DO IREC=M1,M2
            IG=IG+1
            WRITE(10,REC=IREC) (RHBA(K,IG),K=1,2*NG)
         END DO
      ENDDO
C--------------------------------------
      RETURN
      END
