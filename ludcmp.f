      SUBROUTINE INVAXB(A,B,N,M)
      INTEGER N,indx(n),np
      complex A(n,n),B(n,m)
 
      np=n
      call ludcmp(A,n,np,indx,d)
      do j=1,m 
         call lubksb(a,n,np,indx,B(1,j))
      enddo  
      return 
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      PARAMETER (NMAX=5000,TINY=1.0E-20)
	COMPLEX A(NP,NP),SUM,DUMC
      DIMENSION INDX(N),VV(NMAX)
      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
          !IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        IF (AAMAX.EQ.0.) then
            write(*,*)  'Singular matrix.'
            stop
        ENDIF
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUMC=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUMC
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUMC=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUMC
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      COMPlEX A(NP,NP),B(N)
	complex sum
      DIMENSION INDX(N)
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc  Matrix multiplication
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
        SUBROUTINE MATXMAT(A,B,C,N1,N2,N3)
	  integer N1,N2,N3
	  complex :: A(N1,N2),B(N2,N3),C(N1,N3)
      
        DO I=1,N1
	   DO J=1,N3
	      C(I,J)=0.0
	      DO K=1,N2
	         C(I,J)=C(I,J)+A(I,K)*B(K,J)
	      ENDDO
	   ENDDO
	  ENDDO
	  RETURN
	END 
	
	
