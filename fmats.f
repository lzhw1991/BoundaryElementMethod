C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C& FMATS               
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      SUBROUTINE FMATS(ID,ISTEP)
        COMPLEX AICON
        COMPLEX, ALLOCATABLE :: MAT1(:,:),MAT2(:,:),MAT3(:,:),MAT4(:,:)
     &                        ,GM(:,:)
        complex,allocatable :: gatdm1(:,:),gatdp1(:,:),temp(:,:)
        INTEGER N1,N2,N12 
        REAL D(2,2),YY1(2,2)
C
      INCLUDE 'bemlayer.fin'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        AICON=-CMPLX(0.0,0.25)
        DE=RMATE(ID,1)
        VP=RMATE(ID,2)
        VS=RMATE(ID,3)
        N1=2*(NSTA(ID,1,2)-NSTA(ID,1,1)+1)
        N2=2*(NSTA(ID,2,2)-NSTA(ID,2,1)+1)
        NBBP=NBDATA(ID,1)
        N12=2*NBBP
        ALLOCATE(MAT1(N12,N1))
        IF(ID.NE.1) THEN 
          ALLOCATE(MAT2(N12,N1))
        ENDIF
        IF(ID.NE.ND) THEN 
          ALLOCATE(MAT3(N12,N2),MAT4(N12,N2))
        ENDIF
        W=2.0*PI*ISTEP/(NSAMP*DELTA)
        WP=W/VP
        WS=W/VS
        C1=1.0/(DE*VS*VS)
        C3=VS/VP
        C4=C3*C3
        C7=1.0/C4
        PO=0.5*(2.0*C4-1.0)/(C4-1.0)
        B0=2.0*(1.0-PO)
        C9=2.0*PI*B0
        D(1,1)=1.0
        D(1,2)=0.0
        D(2,1)=0.0
        D(2,2)=1.0
C---------------------------
        PRINT*,'NBBP=',NBBP
        DO 100 ISEG=1,NBBP
          XPT=COORD(ID,ISEG,1)
          YPT=COORD(ID,ISEG,2)
          ALPF1=COORD(ID,ISEG,3)
          ALPF2=COORD(ID,ISEG,4)
          SALP=SIN(ALPF2)**2-SIN(ALPF1)**2
          SSALP=SIN(2.0*ALPF2)-SIN(2.0*ALPF1)
          YY1(1,1)=((ALPF2-ALPF1)*B0+0.5*SSALP)/C9
          YY1(2,2)=((ALPF2-ALPF1)*B0-0.5*SSALP)/C9
          YY1(1,2)=SALP/C9
          YY1(2,1)=YY1(1,2)
          CALL INTAL(ISEG,XPT,YPT,ID,W,WP,WS,C1,C3,C4,C7,D)
C*----------------------------------------------------
          DO  I=1,2
          DO  J=1,2
            K=ISEG*2-2+J
            BMT(I,K)=BMT(I,K)+YY1(I,J)/AICON
          enddo
          enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          IRA=2*ISEG-2
          I=0
          DO  IC=NSTA(ID,1,2),NSTA(ID,1,1),-1
              ICA=2*IC-2
              I=I+1
              IX=2*I-2
              DO  IH=1,2
              DO  JH=1,2
                 MAT1(IRA+IH,IX+JH)=BMT(IH,ICA+JH)
              enddo
              enddo
          enddo
          IF(ID.NE.1) THEN
          I=0
          DO  IC=NSTA(ID,1,2),NSTA(ID,1,1),-1
             ICA=2*IC-2
             I=I+1
             IX=2*I-2
             DO  IH=1,2
             DO  JH=1,2
                 MAT2(IRA+IH,IX+JH)=AMT(IH,ICA+JH)
             enddo
             enddo
          enddo
          ENDIF
          IF(ID.NE.ND) THEN
          I=0
          DO IC=NSTA(ID,2,1),NSTA(ID,2,2)
             ICA=2*IC-2
             I=I+1
             IX=2*I-2
             DO  IH=1,2
             DO  JH=1,2
                 MAT3(IRA+IH,IX+JH)=BMT(IH,ICA+JH)
             enddo
             enddo
          enddo
          I=0
          DO  IC=NSTA(ID,2,1),NSTA(ID,2,2)
              ICA=2*IC-2
              I=I+1
              IX=2*I-2
              DO  IH=1,2
              DO  JH=1,2
                  MAT4(IRA+IH,IX+JH)=-AMT(IH,ICA+JH)
              enddo
              enddo
          enddo
          endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  100    CONTINUE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ALLOCATE(GM(N12,N12))
         IF(ID.EQ.1) THEN
                GM(1:N12,1:N1)=MAT1
                DEALLOCATE(MAT1)
         ELSE 
                allocate(gatdm1(n1,n1),TEMP(n12,n1))
                call gaload(istep,id-1,gatdm1,n1)
                call MATXMAT(MAT2,GATDm1,TEMP,n12,n1,n1)
                GM(1:N12,1:N1)=MAT1+TEMP
                DEALLOCATE(MAT1,MAT2,gatdm1,temp)
         ENDIF
c---------------------------------
         IF(ID.NE.ND) THEN
            allocate(gatdp1(n2,n2),temp(n12,n2))
            call gaload(istep,id+1,gatdp1,n2)
            CALL MATXMAT(MAT4,GATDp1,TEMP,n12,n2,n2)
            GM(1:N12,N1+1:N12)=MAT3+TEMP
            DEALLOCATE(MAT3,MAT4,TEMP)
         ENDIF
         CALL FRHM(ID,ISTEP)
         call invaxb(gm,resu(1:N12,1:NFS),N12,NFS)
         DISP(ID,1:N1,1:NFS)=resu(1:N1,1:NFS)
         if(id.ne.nd) then 
            DISP(ID+1,1:N2,1:NFS)=resu(N1+1:N12,1:NFS)
            TRACT(ID+1,1:N2,1:NFS)=MATMUL(gatdp1,resu(N1+1:N12,1:NFS))
            deallocate(gatdp1)
         endif
         DEALLOCATE(GM) 
         print*,'endof fmats'
         RETURN
       END
