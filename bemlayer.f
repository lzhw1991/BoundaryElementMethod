CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C--Elastic modeling for layering media
C  with irregular interfaces
C--Forward method to assemble matrix
C--Point source in the surface domain ID=1
C--Using linear elements
C-------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM EBEM
      INCLUDE 'bemlayer.fin'
      complex, allocatable:: gatd(:,:),gadd(:,:) 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL EDITP
c******************************************
      LECL10=2*(2*NG) 
      OPEN(10,FILE=TEMNAM,ACCESS='DIRECT',RECL=4*LECL10)
      MAB=0
      DO 999 ISTEP=KSTEP,NSTEP
      MAB=MAB+1
      WRITE(*,*) ISTEP
      MXX=0
      MYY=0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO ID=1,ND
         CALL SHLI(ISTEP,ID)
      ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Compute the coeffient matrix above the source layer
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO  ID=1,LS-1
         CALL FMAT1(ID,ISTEP)
      ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Compute the coeffient matrix below the source layer
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     
      DO ID=ND,LS+1,-1
         CALL FMAT2(ID,ISTEP)
      ENDDO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc when ID=LS SOLVE the matrix
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CALL FMATS(LS,ISTEP)
      RESU=(0.0,0.0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Complute the displacement and traction along the boundries
c 
      DO ID=LS-1,1,-1
         N1=2*(NSTA(ID,1,2)-NSTA(ID,1,1)+1)
         N2=2*(NSTA(ID,2,2)-NSTA(ID,2,1)+1)
         allocate(gadd(N1,N2),gatd(n2,n2))
         call gaload(istep,id,gatd,n2)
         call gaload1(istep,id,gadd,n1,n2)
         DISP(ID,1:N1,1:NFS)
     *         =MATMUL(GADD,DISP(ID+1,1:N2,1:NFS))
         TRACT(ID+1,1:N2,1:NFS)
     *         =MATMUL(GATD,DISP(ID+1,1:N2,1:NFS))
         deallocate(gadd,gatd)
      ENDDO
      IF(.FALSE.) THEN
        DO ID=LS+1,ND-1
           N1=2*(NSTA(ID,1,2)-NSTA(ID,1,1)+1)
           N2=2*(NSTA(ID,2,2)-NSTA(ID,2,1)+1)
           allocate(gadd(n2,N1),gatd(n1,n1))
           call gaload(istep,id,gatd,n1)
           call gaload1(istep,id,gadd,n2,n1)
           DISP(ID+1,1:N2,1:NFS)=
     *        MATMUL(GADD,DISP(ID,1:N1,1:NFS))
           TRACT(ID,1:N1,1:NFS)=
     *        MATMUL(GATD,DISP(ID,1:N1,1:NFS))
           deallocate(gadd,gatd)
        ENDDO
        IF(ND.NE.LS) THEN
        ID=ND
        N1=2*(NSTA(ID,1,2)-NSTA(ID,1,1)+1)
        allocate(gatd(n1,n1))
        call gaload(istep,id,gatd,n1)
        TRACT(ID,1:N1,1:NFS)=
     *      MATMUL(GATD,DISP(ID,1:N1,1:NFS))
        deallocate(gatd)
        ENDIF
      ENDIF
C-------------------------------------------
      CALL BIPTS(ISTEP,MAB)
C---------------------------------------------
999   CONTINUE
      CLOSE(10)
C---------------------------------------------
      OPEN(10,FILE=TEMNAM,ACCESS='DIRECT',RECL=4*LECL10)
      LECL11=LENGTH 
      PRINT*,'LECL11=',LECL11
      OPEN(11,FILE=OUTNAMU,ACCESS='DIRECT',RECL=4*LECL11)
      OPEN(12,FILE=OUTNAMW,ACCESS='DIRECT',RECL=4*LECL11)
      CALL FRES
      CLOSE(10)
      CLOSE(12)
      CLOSE(11)
C----------------------------------
      END
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C& EDITP 
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        SUBROUTINE EDITP
        CHARACTER(77) DFILE
        CHARACTER(6) DOMAIN
        REAL TEMPX(2000),TEMPY(2000),TEMPR(2000)
        INTEGER TTR
C       
        INCLUDE 'bemlayer.fin'             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
        WRITE(*,*) 'Input the PARfile name.'
C        READ(*,'(a77)') DFILE
        OPEN(2,FILE='ebem.dat',ACCESS='SEQUENTIAL')
        READ(2,'(64X)')
C---------SET GLOBAL PARAMETERS----------- 
        READ(2,*) ND,NFS,NGQ,NREQ,TLEN,DELTA,NS,DAMP,TSHIFT
        WRITE(*,*) ND,NFS,NGQ,NREQ,TLEN,DELTA,NS,DAMP,TSHIFT
C        NFS=1        
        LENGTH=ANINT(TLEN/DELTA)+1
C--------------------------------------------------
        READ(2,*) WLEN,F0,FMIN,FMAX,XMIN,XMAX
        WRITE(*,*) WLEN,F0,FMIN,FMAX,XMIN,XMAX
        XMIN=XMIN
        XMAX=XMAX
C-------------------------------------------------
        ANG2RAD=PI/180.0
        DO IFS=1,NFS
          READ(2,*) ISOURCE(IFS),XSEC(IFS),YSEC(IFS),SFLAG(IFS),
     &            THETA(IFS),THETA1(IFS),THETA2(IFS)
          write(*,*) ISOURCE(IFS),XSEC(IFS),YSEC(IFS),SFLAG(IFS),
     &            THETA(IFS),THETA1(IFS),THETA2(IFS)
          LS=ISOURCE(IFS)
          XSEC(IFS)=XSEC(IFS)
          YSEC(IFS)=YSEC(IFS)
          THETA(IFS)=THETA(IFS)*ANG2RAD
          THETA1(IFS)=THETA1(IFS)*ANG2RAD
          THETA2(IFS)=THETA2(IFS)*ANG2RAD
        ENDDO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Change the input format of receiver positions, 09/09/2005
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        READ(2,*) NREC
        open(55, file='receiver.dat',status='old')
        NG=NREC
        K=0
        DO IREC=1,NREC
           READ(55,*), RECFLAG(IREC,1),RECFLAG(IREC,2),XREC,YREC
           ROCOOR(IREC,1)=XREC
           ROCOOR(IREC,2)=YREC
        ENDDO
        CLOSE(55)        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         READ(2,'(a77)') TEMNAM
         WRITE(*,'(a77)') TEMNAM
         READ(2,'(a77)') OUTNAMU
         WRITE(*,'(a77)') OUTNAMU
         READ(2,'(a77)') OUTNAMW
         WRITE(*,'(a77)') OUTNAMW
C-----------------------------------------
        DO I=1,ND
         READ(2,*) DOMAIN,ID,DEN,PVEL,SVEL,MBN
         WRITE(*,*) DOMAIN,ID,DEN,PVEL,SVEL,MBN
         RMATE(I,1)=DEN
         RMATE(I,2)=PVEL
         RMATE(I,3)=SVEL
         MBDATA(I)=MBN
         IF(I.eq.1) then
            READ(2,*) N1S,N1E,N2S,N2E
            WRITE(*,*) N1S,N1E,N2S,N2E
         else
            READ(2,*) N2S,N2E
            WRITE(*,*)N2S,N2E
         endif
         IF(N2E.NE.MBN) THEN
          PRINT*,'MBN and N is not same'
          stop
         endif
            MSTA(I,1,1)=N1S
            MSTA(I,1,2)=N1E
            MSTA(I,2,1)=N2S
            MSTA(I,2,2)=N2E
C-----------SET REPEAT PARAMETERS--------
       if(id.eq.1) then
           DO J=1,MBN
             READ(2,*) IP,XD,YD
             WRITE(*,*) IP,XD,YD
             XCOORD(I,J,1)=XD
             XCOORD(I,J,2)=YD
           ENDDO
        else
           DO J=N2S,N2E
             READ(2,*) IP,XD,YD
             WRITE(*,*) IP,XD,YD
             XCOORD(I,J,1)=XD
             XCOORD(I,J,2)=YD
           ENDDO
         endif
        ENDDO
        CLOSE(2)
C----------------------------------------
      IF(NREQ.EQ.0)THEN 
      IDIV=2*ANINT(WLEN/DELTA)+1
      ELSE
      IDIV=ANINT(WLEN/DELTA)+1      
      ENDIF
      IDIV=(IDIV/2)*2+1

      NSAMP=512
      IF(LENGTH.GT.512)NSAMP=1024
      IF(LENGTH.GT.1024)NSAMP=2048
      IF(LENGTH.GT.2048)NSAMP=4096
      IF(LENGTH.GT.4096)NSAMP=8192

      I=0
      N=1
 5    I=I+1
      N=2*N
      IF(N.EQ.NSAMP)THEN
        NN=I
      ELSE
      GO TO 5
      ENDIF

      KSTEP=NINT(FMIN*NSAMP*DELTA)
      IF(KSTEP.EQ.0)KSTEP=1
      NSTEP=NINT(FMAX*NSAMP*DELTA+0.5)
c      print*,NSTEP,FMAX,'xixix',FMAX*NSAMP*DELTA+0.5
C-----------DESIGN WAVELET-----------
      DO I=0,NSAMP-1
      REPA(I)=0.0
      AIMPA(I)=0.0
      END DO
      DF=1.0/(NSAMP*delta)
      print*,'df=',df
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Source time function
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SELECT CASE (NREQ)
      CASE(0)
         CALL FSOU(WLEN,IDIV,DELTA,REPA,NSAMP0)
         CALL FFT(NN,1,NSAMP,REPA,AIMPA,WK1,WK2,NSAMP0)
      CASE(1)
          CALL RICKER(kstep,nstep,F0,df,REPA,aimpa,NSAMP0, TSHIFT)
      CASE(2)
          CALL GAUSS(IDIV,REPA,NSAMP0)
          CALL FFT(NN,1,NSAMP,REPA,AIMPA,WK1,WK2,NSAMP0)
      END SELECT
C----------------------------------------
C-----------------------------------------
      WGP=0.0
      SHN=0.0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Shape function and weigh for Gaussain quad
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CALL SHAP(SHN,WGP)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CALL ANGLA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RETURN
      END
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C& ANGLA   
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      SUBROUTINE ANGLA
C  
      INCLUDE 'bemlayer.fin'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 999 ID=1,ND
      DO 888 IC=1,2  
      N1=MSTA(ID,IC,1) 
      N2=MSTA(ID,IC,2)
      IF(IC.EQ.1)THEN 
      BLPH(ID,N1,1)=0.0
      BLPH(ID,N1,2)=PI
      BLPH(ID,N2,1)=0.0
      BLPH(ID,N2,2)=PI
      ELSE
      BLPH(ID,N1,1)=PI
      BLPH(ID,N1,2)=2.0*PI        
      BLPH(ID,N2,1)=PI
      BLPH(ID,N2,2)=2.0*PI
      ENDIF
      DO 100 IQ=N1+1,N2-1       
      X0=XCOORD(ID,IQ-1,1)
      Y0=XCOORD(ID,IQ-1,2)
      X1=XCOORD(ID,IQ,1)
      Y1=XCOORD(ID,IQ,2)
      X2=XCOORD(ID,IQ+1,1)
      Y2=XCOORD(ID,IQ+1,2)
      ALFA1=ATAN2(Y0-Y1,X0-X1)
      ALFA2=ATAN2(Y2-Y1,X2-X1)
      IF(ALFA2.lt.ALFA1) ALFA2=ALFA2+2*PI
      BLPH(ID,IQ,1)=ALFA1
      BLPH(ID,IQ,2)=ALFA2
100   CONTINUE
888   CONTINUE
999   CONTINUE
C-----------------------------------------
      RETURN
      END
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C& FREQUENCY TO TIME DOMAIN                
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      SUBROUTINE FRES
C
      INCLUDE 'bemlayer.fin'
      complex, allocatable :: temp(:,:)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      MAX=(NSTEP-KSTEP+1)*NFS
      IMTY=NFS
      LSAMP=NSAMP-1
      NENGTH=LENGTH-1
      allocate(temp(2*NG,NFS*NSTEP))
C-------------------------------

      DO 200 IFS=1,NFS
         NCC=(IFS-1)*2*NG
         NCUT=0
         DO IREC=IFS,MAX,IMTY
            NCUT=NCUT+1
            READ(10,REC=IREC) (temp(K,NCUT),K=1,2*NG)
         END DO
         IK=(IFS-1)*NG
         DO 400 IGT=1,2*NG,2
            IK=IK+1 
            DO K=0,LSAMP
               REPA(K)=0.0
               AIMPA(K)=0.0
            END DO
            ICOUNT=0
            DO I=KSTEP,NSTEP
                ICOUNT=ICOUNT+1
                REPA(I)=REAL(temp(IGT,ICOUNT))
                REPA(NSAMP-I)=REPA(I)
                AIMPA(I)=AIMAG(temp(IGT,ICOUNT))
                AIMPA(NSAMP-I)=-AIMPA(I)
            END DO
            CALL FFT(NN,-1,1,REPA,AIMPA,WK1,WK2,NSAMP0)  
            IREC=IK
            WRITE(11,REC=IREC) (REPA(K),K=0,NENGTH)
400      CONTINUE
C-------------------------------------------
         JK=(IFS-1)*NG
         DO 500 IGT=2,2*NG,2
            JK=JK+1 
            DO K=0,LSAMP
               REPA(K)=0.0
               AIMPA(K)=0.0
            END DO
            ICOUNT=0
            DO I=KSTEP,NSTEP
               ICOUNT=ICOUNT+1
               REPA(I)=REAL(temp(IGT,ICOUNT))
               REPA(NSAMP-I)=REPA(I)
               AIMPA(I)=AIMAG(temp(IGT,ICOUNT))
               AIMPA(NSAMP-I)=-AIMPA(I)
            END DO
            CALL FFT(NN,-1,1,REPA,AIMPA,WK1,WK2,NSAMP0)  
            IREC=JK
            WRITE(12,REC=IREC) (REPA(K),K=0,NENGTH)
500      CONTINUE
200   CONTINUE
      deallocate(temp)
C-------------------------------------
      RETURN
      END
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c& INTAL                
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        SUBROUTINE INTAL(ISEG,XPT,YPT,ID,W,WP,WS,C1,C3,C4,C7,D)
        COMPLEX H(2,4),G(2,4),UL(2,2),PL(2,2)
        COMPLEX P0,P1,P2,P3,S0,S1,S2,S3
        COMPLEX F1,F2,DF1,DF2,F3,F4,F5,F6
        REAL D(2,2),DXY(2),BN(2),DR(2),B(2)
C
        INCLUDE 'bemlayer.fin'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        NBBP=NBDATA(ID,1)
        NBBE=NBDATA(ID,2)      
        DO 10 I=1,2
        DO  J=1,2*(NBBP+1)
        AMT(I,J)=(0.0,0.0)
        BMT(I,J)=(0.0,0.0)
        ENDDO
   10   CONTINUE
      DO 900 IS=1,NBBE
        ISA=NEL(ID,IS,1)
        ISB=NEL(ID,IS,2)
        XA=COORD(ID,ISA,1)
        YA=COORD(ID,ISA,2)
        XB=COORD(ID,ISB,1)
        YB=COORD(ID,ISB,2)
        DO  K=1,2
        DO  L=1,4
        H(K,L)=(0.0,0.0)
        G(K,L)=(0.0,0.0)
        ENDDO
        ENDDO
        DXY(1)=XB-XA
        DXY(2)=YB-YA
        RAB=SQRT(DXY(1)*DXY(1)+DXY(2)*DXY(2))
        BN(1)=-DXY(2)/RAB
        BN(2)=DXY(1)/RAB
c-----------------------------
        JGG=NGP-1
      DO 50 IG=1,NGP        
        XS1=SHN(JGG,IG,1)
        XS2=SHN(JGG,IG,2)
        XGP=XS1*XA+XS2*XB
        YGP=XS1*YA+XS2*YB
        XMX=XGP-XPT
        YMY=YGP-YPT
        R=SQRT(XMX*XMX+YMY*YMY)
        B(1)=0.5*XS1*RAB
        B(2)=0.5*XS2*RAB
        IF(R.LT.0.0001) R=0.0001
        DR(1)=XMX/R
        DR(2)=YMY/R
        DRDN=DR(1)*BN(1)+DR(2)*BN(2)
        ZP=WP*R
        ZS=WS*R
        CALL HANKEX(ZP,P0,P1)
        CALL HANKEX(ZS,S0,S1)
        P2=2.0*P1/ZP-P0
        P3=4.0*P2/ZP-P1
        S2=2.0*S1/ZS-S0
        S3=4.0*S2/ZS-S1
        F1=S0+(C3*P1-S1)/ZS
        F2=C4*P2-S2
        ZSD=ZS*ZS
        DF1=-2.*WP*P1/ZSD-S0/R-WS*S1+2.*S1/(ZS*R)+C4*P0/R
        DF2=WP*C4*(P1-P3)/2.-WS*(S1-S3)/2.
        DO  I=1,2
        DO  J=1,2
        UL(I,J)=C1*(F1*D(I,J)-F2*DR(I)*DR(J))
        F3=(DF1-F2/R)*(D(I,J)*DRDN+DR(J)*BN(I))
        F4=2.*F2*(DR(I)*BN(J)-2.*DR(I)*DR(J)*DRDN)/R
        F5=2.*DF2*DR(I)*DR(J)*DRDN
        F6=(C7-2.)*(DF1-DF2-F2/R)*DR(I)*BN(J)
        PL(I,J)=-F3-F6+F5+F4
        ENDDO
        ENDDO
        DO  LA=1,2
        IC=0
        DO  LL=1,2
        DO  JJ=1,2
        IC=IC+1
        G(LA,IC)=G(LA,IC)+UL(LA,JJ)*B(LL)*WGP(JGG,IG)
        H(LA,IC)=H(LA,IC)+PL(LA,JJ)*B(LL)*WGP(JGG,IG)
        ENDDO
        ENDDO
        ENDDO
   50   CONTINUE
        DO 60 IM=1,2
        AMT(IM,2*ISA-1)=G(IM,1)+AMT(IM,2*ISA-1)
        AMT(IM,2*ISA)=G(IM,2)+AMT(IM,2*ISA)
        AMT(IM,2*ISB-1)=G(IM,3)+AMT(IM,2*ISB-1)
        AMT(IM,2*ISB)=G(IM,4)+AMT(IM,2*ISB)
        BMT(IM,2*ISA-1)=H(IM,1)+BMT(IM,2*ISA-1)
        BMT(IM,2*ISA)=H(IM,2)+BMT(IM,2*ISA)
        BMT(IM,2*ISB-1)=H(IM,3)+BMT(IM,2*ISB-1)
        BMT(IM,2*ISB)=H(IM,4)+BMT(IM,2*ISB)
   60   CONTINUE
 900    CONTINUE
c-----------------------------------
      RETURN
      END
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C& FSOU 
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      SUBROUTINE FSOU(WLEN,IDIV,DT,REPA,NSAMP1)
      REAL REPA(0:NSAMP1)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C A:0~1.6HZ; B:0~2HZ; C:0~2.5HZ; D:0~5HZ
C WLEN=4000MS=4S
C DT=100MS=0.1S      
C-----------------------------------      
        SGM=1.0
          DO II=1,IDIV
        CC1=DT*II-WLEN
        CC=-2.0*CC1*CC1/SGM/SGM
        CC2=WLEN/SGM
        CC=(EXP(CC)-EXP(-2.0*CC2*CC2))/1.2533/SGM
        CC=-4.0/SGM/SGM*CC*CC1
        REPA(II-1)=CC
        END DO
C------------------------------    
              RETURN
        END
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C& GAUSS 
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        SUBROUTINE GAUSS(IDIV,REPA,NSAMP1)
        REAL REPA(0:NSAMP1)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-----------------------------------      
        C8=4.0
        C9=0.1
        SGM=1.0
          DO II=1,IDIV
        CC1=C9*II-C8
        CC=-2.0*CC1*CC1/SGM/SGM
        CC2=C8/SGM
        CC=(EXP(CC)-EXP(-2.0*CC2*CC2))/1.2533/SGM
        CC=-4.0*CC*CC1
        REPA(II-1)=-CC
        END DO
C------------------------------    
        RETURN
        END
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C RICKER
Cccccccccccccccccccccccccccccc
      subroutine ricker(kstep,nstep,f0,df,repa,aimpa,nsamp1,TSHIFT)
c subroutine to computer ricker wavelet in frequency domain
      real repa(0:nsamp1),aimpa(0:nsamp1)
      real f0,ts
      pi=3.14159265399
      !ts=3.0/f0
      ts=-TSHIFT
      do i=kstep,nstep
         f=i*df
         temp=-(f/f0)**2.0*exp(-(f/f0)**2.0)
         repa(i)=temp*cos(2*pi*f*ts)
         aimpa(i)=temp*sin(2*pi*f*ts)
      enddo   
      return
      end
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C& FFT             
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      SUBROUTINE FFT(NN,INT,NCONS,A1,A2,A3,A4,NSAMP1)  
      REAL A1(0:NSAMP1),A2(0:NSAMP1)
      REAL W1(0:NSAMP1),W2(0:NSAMP1)
      REAL A3(0:NSAMP1),A4(0:NSAMP1)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      N=2**NN
      CONS=1.0/FLOAT(NCONS)
      YY=FLOAT(INT)/FLOAT(N)
      DO 40 I=0,N/2-1
      Y=6.2832*I*YY
      Z1=COS(Y)
      Z2=SIN(Y)
      W1(I)=Z1
      W2(I)=Z2
 40   CONTINUE
      N2=2**(NN-1)
      DO 50 I=1,NN
      LK=2**(NN-I)-1
      LJ=2**(I-1)-1
      N1=2**(I-1)
      N3=2**I
      MX=(I/2)*2-I
      IF(MX.NE.0) THEN
      DO K=0,LK
      N4=N1*K
      N5=N3*K
      DO J=0,LJ
      I1=N5+J
      I2=N4+J
      I3=I2+N2
      I4=I1+N1
      I5=N4
      A3(I1)=A1(I2)+A1(I3)
      A4(I1)=A2(I2)+A2(I3)
      X1=A1(I2)-A1(I3)
      X2=A2(I2)-A2(I3)
      A3(I4)=X1*W1(I5)-X2*W2(I5)
      A4(I4)=X1*W2(I5)+X2*W1(I5)
      ENDDO
      ENDDO
      ELSE
      DO K=0,LK
      N4=N1*K
      N5=N3*K
      DO J=0,LJ
      I1=N5+J
      I2=N4+J
      I3=I2+N2
      I4=I1+N1
      I5=N4
      A1(I1)=A3(I2)+A3(I3)
      A2(I1)=A4(I2)+A4(I3)
      X1=A3(I2)-A3(I3)
      X2=A4(I2)-A4(I3)
      A1(I4)=X1*W1(I5)-X2*W2(I5)
      A2(I4)=X1*W2(I5)+X2*W1(I5)
      ENDDO
      ENDDO
      ENDIF
 50   CONTINUE
      MX=(NN/2)*2-NN
      IF(MX.NE.0) THEN
      DO I=0,N-1
      A1(I)=A3(I)
      A2(I)=A4(I)
      ENDDO
      ELSE
      ENDIF
      DO I=0,N-1
      A1(I)=A1(I)*CONS
      A2(I)=A2(I)*CONS
      ENDDO
C----------------------------
      RETURN
      END
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C& SHAP            
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      SUBROUTINE SHAP(SHN,WGP)  
      REAL SGP(7,8),WEG(7,8),SHN(7,8,2),WGP(7,8)
C
C-----------FOR LINEAR ELEMENT-------------------------
C------------------------------------------------------------
      DATA(SGP(1,I),I=1,2)/-0.57735027,0.57735027/
      DATA(WEG(1,I),I=1,2)/1.000000,1.000000/
      DATA(SGP(2,I),I=1,3)/-0.774597,0.000000,0.774597/
      DATA(WEG(2,I),I=1,3)/0.555556,0.888889,0.555556/
      DATA(SGP(3,I),I=1,4)/-0.861136,-0.339981,0.339981,0.861136/
      DATA(WEG(3,I),I=1,4)/0.347855,0.652145,0.652145,0.347855/
      DATA(SGP(4,I),I=1,5)/-0.906179,-0.538469,0.0,0.538469,0.906179/
      DATA(WEG(4,I),I=1,5)/0.236927,0.478629,0.568889,0.478629,
     &  0.236927/
      DATA(SGP(5,I),I=1,6)/-0.932470,-0.661210,-0.238619,0.238619,
     &  0.661210,0.932470/
      DATA (WEG(5,I),I=1,6)/0.171324,0.360761,0.467914,0.467914,
     &  0.360761,0.171324/
      DATA(SGP(6,I),I=1,7)/-0.949110,-0.741531,-0.405845,0.000000,
     &  0.405845,0.741531,0.949110/
      DATA(WEG(6,I),I=1,7)/0.129485,0.279710,0.381830,0.417960,0.381830,
     &  0.279710,0.129485/
      DATA(SGP(7,I),I=1,8)/-0.960290,-0.796670,-0.525530,-0.183435,
     &  0.183435,0.525530,0.796670,0.960290/
      DATA(WEG(7,I),I=1,8)/0.101220,0.222380,0.313710,0.362680,0.362680,
     &  0.3213710,0.222380,0.101220/
      DO I=1,7
      DO J=1,I+1
      WGP(I,J)=WEG(I,J)
      ENDDO
      ENDDO
      DO I=1,7
      NT=I+1
      DO  J=1,NT
      XSS=SGP(I,J)
      SHN(I,J,1)=0.5*(1.0-XSS)
      SHN(I,J,2)=0.5*(1.0+XSS)
      ENDDO
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCC
      RETURN
      END
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C& HANKEX                
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      SUBROUTINE HANKEX(Z,H0,H1)
      COMPLEX H0,H1
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FR=ALOG(Z/2.0)
      IF(Z.LE.3.0)THEN
        T=(Z/3.0)**2
        BJ0=(((((0.00021*T-0.003944)*T+0.0444479)*T-
     &      0.3163866)*T+1.2656208)*T-2.2499997)*T+1.0
        YJ0=(((((-0.00024846*T+0.00427916)*T-0.04261214)*T+0.25300117)*
     &      T-0.74350384)*T+0.60559366)*T+0.36746691+0.6366197*FR*BJ0
        H0=CMPLX(BJ0,YJ0)
        BJ1=((((((0.1109E-4*T-0.31761E-3)*T+4.43319E-3)*
     &      T-0.03954289)*T+0.21093573)*T-0.56249985)*T+0.5)*Z
        YJ1=((((((0.0027873*T-0.0400976)*T+0.3123951)*T-1.3164827)*T+
     &      2.1682709)*T+0.2212091)*T-0.6366198+0.6366197*Z*FR*BJ1)/Z
        H1=CMPLX(BJ1,YJ1)
      ELSE
        T=3.0/Z
        F0=(((((0.00014476*T-0.00072805)*T+0.00137237)*T-0.9512E-4)*T-
     &      0.55274E-2)*T-0.77E-6)*T+0.79788456
        FA0=(((((-0.13558E-3*T+0.29333E-3)*T+0.54125E-3)*T-0.00262573)*
     &  T+0.3954E-4)*T+0.04166397)*T+0.78539816
        BJ0=F0*COS(Z-FA0)/SQRT(Z)
        YJ0=F0*SIN(Z-FA0)/SQRT(Z)
        H0=CMPLX(BJ0,YJ0)
        F1=(((((-2.0033E-4*T+1.13653E-3)*T-2.49511E-3)*T+1.7105E-4)*T+
     &      0.01659667)*T+1.56E-6)*T+0.79788456
        FA1=(((((2.9166E-4*T-7.9824E-4)*T-7.4348E-4)*T+6.37879E-3)*T-
     &      5.65E-5)*T-0.12499612)*T+2.35619449
        BJ1=F1*COS(Z-FA1)/SQRT(Z)
        YJ1=F1*SIN(Z-FA1)/SQRT(Z)
        H1=CMPLX(BJ1,YJ1)
      ENDIF 
C--------------------------------------
        RETURN
      END

























