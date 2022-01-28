c     SUBROUTINE simplx(A,M,N,MP,NP,M1,M2,M3,ICASE,IZROV,IPOSV)
      SUBROUTINE simplx(A,M,N,MP,NP,M1,M2,M3,ICASE,IZROV,IPOSV,
     *                  NMAX,L1,MMAX,L2,L3)
      implicit real*8 (a-h,o-z)
c_old PARAMETER(MMAX=100,EPS=1.E-6)
      PARAMETER(EPS=1.E-6)
      DIMENSION A(MP,NP),IZROV(N),IPOSV(M)
      DIMENSION L1(NMAX),L2(MMAX),L3(MMAX)
c_old DIMENSION L1(MMAX),L2(MMAX),L3(MMAX)
      i10 = 0
      i20 = 0
      IF(M.NE.M1+M2+M3)PAUSE 'Bad input constraint counts.'
      NL1=N
      DO 11 K=1,N
        L1(K)=K
        IZROV(K)=K
11    CONTINUE
      NL2=M
      DO 12 I=1,M
        IF(A(I+1,1).LT.0.)PAUSE 'Bad input tableau.'
        L2(I)=I
        IPOSV(I)=N+I
12    CONTINUE
      DO 13 I=1,M2
        L3(I)=1
13    CONTINUE
      IR=0
      IF(M2+M3 .EQ. 0) GO TO 30
      IR=1
      DO 15 K=1,N+1
        Q1=0.
        DO 14 I=M1+1,M
          Q1=Q1+A(I+1,K)
14      CONTINUE
        A(M+2,K)=-Q1
15    CONTINUE
c      write(*,*) "SIMP1"
10    CALL SIMP1(A,MP,NP,M+1,L1,NL1,0,KP,BMAX)
      IF(BMAX.LE.EPS.AND.A(M+2,1).LT.-EPS)THEN
        ICASE=-1
        RETURN
      ELSE IF(BMAX.LE.EPS.AND.A(M+2,1).LE.EPS)THEN
        M12=M1+M2+1
        IF(M12.LE.M)THEN
          DO 16 IP=M12,M
            IF(IPOSV(IP).EQ.IP+N)THEN
              CALL SIMP1(A,MP,NP,IP,L1,NL1,1,KP,BMAX)
              IF(BMAX.GT.0.)GO TO 1
            ENDIF
16        CONTINUE
        ENDIF
        IR=0
        M12=M12-1
        IF(M1+1.GT.M12)GO TO 30
        DO 18 I=M1+1,M12
          IF(L3(I-M1).EQ.1)THEN
            DO 17 K=1,N+1
              A(I+1,K)=-A(I+1,K)
17          CONTINUE
          ENDIF
18      CONTINUE
        GO TO 30
      ENDIF
c      write(*,*) "SIMP2"
      CALL SIMP2(A,M,N,MP,NP,L2,NL2,IP,KP,Q1)
      IF(IP .EQ. 0)THEN
        ICASE = -1
        RETURN
      ENDIF
1     CALL SIMP3(A,MP,NP,M+1,N,IP,KP)
      IF(IPOSV(IP).GE.N+M1+M2+1)THEN
        DO 19 K=1,NL1
          IF(L1(K).EQ.KP)GO TO 2
19      CONTINUE
2       NL1=NL1-1
        DO 21 IS=K,NL1
          L1(IS)=L1(IS+1)
21      CONTINUE
      ELSE
        IF(IPOSV(IP).LT.N+M1+1)GO TO 20
        KH=IPOSV(IP)-M1-N
        IF(L3(KH).EQ.0)GO TO 20
        L3(KH)=0
      ENDIF
      A(M+2,KP+1)=A(M+2,KP+1)+1.
      DO 22 I=1,M+2
        A(I,KP+1)=-A(I,KP+1)
22    CONTINUE
20    IS=IZROV(KP)
      IZROV(KP)=IPOSV(IP)
      IPOSV(IP)=IS
      i10 = i10 + 1
      write(*,*) "i10 =", i10
      IF(IR .NE. 0) GO TO 10
30    CALL SIMP1(A,MP,NP,0,L1,NL1,0,KP,BMAX)
      IF(BMAX.LE.0.)THEN
        ICASE=0
        RETURN
      ENDIF
      CALL SIMP2(A,M,N,MP,NP,L2,NL2,IP,KP,Q1)
      IF(IP.EQ.0)THEN
        ICASE = 1
        RETURN
      ENDIF
      CALL SIMP3(A,MP,NP,M,N,IP,KP)
      i20 = i20 + 1
      write(*,*) "i20 =", i20
      GO TO 20

      END
