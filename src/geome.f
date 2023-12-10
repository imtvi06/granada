      SUBROUTINE GEOME (NOAT, coord, ncg, nbg, nag)
      include 'include.inc'

C     LECTURA DE LOS INDICADORES Y DE LAS COORDENADAS INTERNAS
C     CALCULO DE LAS COORDENADAS CARTESIANAS DE LOS ATOMOS EN EL SISTEMA
C
      dimension rcdg(NATMAX), tbcdg(NATMAX), pabcg(NATMAX),
     &          x(NATMAX), y(NATMAX), z(NATMAX),
     &          nag(*), nbg(*), ncg(*), coord(3,*)
      common /geos/ geo(3,NATMAX)

C SEXRAD ES EL FACTOR DE CONVERSION DE UNIDADES SEXAGESIMALES A RADIANES

      real*8 SEXRAD /.0174532d0/

      do ii=1,NOAT
        rcdg(ii)=geo(1,ii)
        tbcdg(ii)=geo(2,ii)
        pabcg(ii)=geo(3,ii)
      enddo

      THETA=tbcdg(3)*SEXRAD
      CCOS=COS(THETA)
      SSIN=SIN(THETA)


4     DO I=1,3
        X(I)=CERO
        Y(I)=CERO
        Z(I)=CERO
      enddo
      X(2)=rcdg(2)
      X(3)=rcdg(2)-rcdg(3)*CCOS
      Y(3)=rcdg(3)*SSIN

      DO I=4,NOAT
        X(I)=10000.d0
        y(I)=CERO
        z(I)=CERO
      enddo
      DO 52 I=4,NOAT
        NA = NAG(I)
        NB = NBG(I)
        NC = NCG(I)
        ND = I
        RCD = RCDG(I)
        THBCD = TBCDG(I)
        PHABCD = PABCG(I)
        TAF=X(NA)+X(NB)+X(NC)
        IF (TAF.GE.7000.D0) stop 

C MOVIMIENTO DEL ATOMO C AL ORIGEN

79      XA=X(NA)-X(NC)
        YA=Y(NA)-Y(NC)
        ZA=Z(NA)-Z(NC)
        XB=X(NB)-X(NC)
        YB=Y(NB)-Y(NC) 
        ZB=Z(NB)-Z(NC)
C ROTACION ALREDEDOR DEL EJE Z PARA HACER YB=0,XB+VE. SI XYB ES
C PEQUENO,ROTACION PRIMERO DE 90 GRADOS ALREDEDOR DEL EJE Y
        XYB=SQRT(XB**2+YB**2)
        K=1
        IF (XYB.LT.0.1D0) then
          K=0
          XPA=ZA
          ZPA=-XA
          XA=XPA
          ZA=ZPA
          XPB=ZB
          ZPB=-XB
          XB=XPB
          ZB=ZPB
          XYB=SQRT(XB**2+YB**2)
        endif
        COSTH=XB/XYB
        SINTH=YB/XYB
        XPA=XA*COSTH+YA*SINTH
        YPA=YA*COSTH-XA*SINTH
C ROTACION ALREDEDOR DEL EJE Y PARA ANULAR ZB
11      RBC=SQRT(XB**2+YB**2+ZB**2)
        SINPH=ZB/RBC
        COSPH=SQRT(1.0-SINPH**2)
        XQA=XPA*COSPH+ZA*SINPH
        ZQA=ZA*COSPH-XPA*SINPH
C ROTACION ALREDEDOR DEL EJE X PARA HACER Z=0,YA+VE
12      YZA=SQRT(YPA**2+ZQA**2)
        COSKH=YPA/YZA
        SINKH=ZQA/YZA
C COORDENADAS A,(XQA,YZA,0),B,(RBC,0,0),C,(0,0,0),NO -VE
C COORDENADAS DE D CALCULADAS AHORA EN UN NUEVO MARCO

        COSD=1.d0
        SIND=CERO
28      THBCD=THBCD*SEXRAD
        PHABCD=PHABCD*SEXRAD
        SINA=SIN(THBCD)
        COSA=COS(THBCD)
        SIND=SIN(PHABCD)
        COSD=COS(PHABCD)
29      CONTINUE
        XD=RCD*COSA
        YD=RCD*SINA*COSD
        ZD=RCD*SINA*SIND
C TRANSFORMACION DE LAS COORDENADAS DE D PARA REGRESAR AL SISTEMA
C ORIGINAL
30      YPD=YD*COSKH-ZD*SINKH
        ZPD=ZD*COSKH+YD*SINKH
        XPD=XD*COSPH-ZPD*SINPH
        ZQD=ZPD*COSPH+XD*SINPH
        XQD=XPD*COSTH-YPD*SINTH
        YQD=YPD*COSTH+XPD*SINTH
        IF (K.NE.10) THEN
          XRD=-ZQD
          ZRD=XQD
          XQD=XRD
          ZQD=ZRD
        ENDIF
        X(ND)=XQD+X(NC)
        Y(ND)=YQD+Y(NC)
        Z(ND)=ZQD+Z(NC)
52    CONTINUE

      do ii=1,NOAT
        coord(1,ii)=x(ii)
        coord(2,ii)=y(ii)
        coord(3,ii)=z(ii)
      enddo


      RETURN
      END

