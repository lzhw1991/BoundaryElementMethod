C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C& FMAT1                
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      SUBROUTINE FMAT1(ID,ISTEP)
      COMPLEX AICON
      COMPLEX,ALLOCATABLE :: HT1(:,:),HT2(:,:),GT1(:,:),GT2(:,:)
      COMPLEX,ALLOCATABLE :: MAT1(:,:),MAT2(:,:)
      complex,allocatable:: gadd(:,:),gatd(:,:),gatdm1(:,:),temp(:,:)
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
        allocate(HT1(N12,N1),HT2(N12,N2),GT2(N12,N2))
        allocate(MAT1(N12,N12))
        IF(ID.NE.1) THEN
          allocate(GT1(N12,N1))
        endif
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
        PRINT*,'NBBP=',NBBP,',   ID=',ID
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

          I=0
          IRA=2*ISEG-2
          DO  IC=NSTA(ID,1,2),NSTA(ID,1,1),-1
              ICA=2*IC-2
              I=I+1
              IX=2*I-2
              DO  IH=1,2
              DO  JH=1,2
                 HT1(IRA+IH,IX+JH)=BMT(IH,ICA+JH)
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
                 GT1(IRA+IH,IX+JH)=-AMT(IH,ICA+JH)
             enddo
             enddo
            enddo
          ENDIF
          I=0
          DO  IC=NSTA(ID,2,1),NSTA(ID,2,2)
              ICA=2*IC-2
              I=I+1
              IX=2*I-2
              DO  IH=1,2
              DO  JH=1,2
                  GT2(IRA+IH,IX+JH)=-AMT(IH,ICA+JH)
              enddo
              enddo
          enddo
          I=0
          DO IC=NSTA(ID,2,1),NSTA(ID,2,2)
             ICA=2*IC-2
             I=I+1
             IX=2*I-2
             DO  IH=1,2
             DO  JH=1,2
                 HT2(IRA+IH,IX+JH)=-BMT(IH,ICA+JH)
             enddo
             enddo
          enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  100    CONTINUE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         allocate(gadd(n1,n2),gatd(n2,n2))
         if(id.eq.1) then
	      MAT1(1:N12,1:N1)=HT1
	   ELSE
	      allocate(gatdm1(n1,n1),temp(n12,n1))
            CALL GALOAD(istep,id-1,GATDm1,N1) 
            call matxmat(GT1,GATDM1,TEMP,N12,N1,N1)
	      MAT1(1:N12,1:N1)=HT1-TEMP
	      deallocate(GATDM1,GT1,TEMP)
	   ENDIF
         deallocate(HT1)
	   MAT1(1:N12,N1+1:N12)=GT2
	   CALL INVAXB(MAT1,HT2,n12,n2)
         GADD=HT2(1:N1,1:N2)
         GATD=HT2(N1+1:N12,1:N2)
         CALL GASAVE(istep,id,GADD,GATD,N1,N2)
         deallocate(HT2,GT2,mat1) 
         deallocate(gadd,gatd)
         MYY=MYY0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCC            
         RETURN
       END
