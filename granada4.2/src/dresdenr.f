      subroutine DRESDENR 
     &  (ii)
      include 'include.inc'
      common /indices/    ninput, qsolv, ncatal, natop, ncard
      common /options/ qroom, qgau, qorca, qmopac, qsolid1, qsolid2
      common /strings/ cadenas(6), title1, title3
      common /datos/      nlim, dimen, nconf,
     &                    nfsolv, fracm(5), ncavid,
     &                    h2o(3,3), nath2o(3), nh2o
      common /dvariables/ aa(3,NATMAX), as(3,NATSMAX), ap(3,NATPMAX),
     &                    asl(3,NATSMAX), d(3), aqpm(3), aqpp(3)
      common /ivariables/ nat(NATMAX), natp(NATPMAX), nats(NATSMAX),
     &                    ntot, nprin, nsolv, nl(9)
      common /filenames/ files, exts
      real*4 Rand
      character*80 cadenas, title1, title3
      common /solid2/ apzmax, apzmin
      common /elements/ elemnt(107)
      character*2 elemnt
      COMMON /simbolos/ isymb
      character*2 ISYMB(107)
      dimension nsrest(9)
      integer*2 dummy, Rseed, sec, csec
      character*128 files(5)
      character*81 titl1s, titl2s
      character*4 exts
      common /solvated/ smop, sxyz, sinp, sgjf, slog
      character(len=128)::smop, sxyz, sinp, sgjf, slog
      

C Inicializacion de contadores y variables
c Para las moleculas de solvente
      
      if (.not.qsolv) then
        nsolv = nh2o
      else
        do i=1,nfsolv
          nsrest(i) = nl(i)
        enddo
      endif

C ntot es el numero de atomos de la molecula o el cluster central en
C este ciclo
      ntot = nprin
      nn = 0
C d(i) contiene ahora las dimensiones de entrada de la caja
      do i=1,3
        d(i) = dimen
      enddo

C EL CONVENIO DE SUBINDICES ESCOGIDO HACE I COMO EL NUMERADOR DE ATOMOS
C Y J COMO EL DE COORDENADAS CARTESIANAS DE CADA ATOMO. Aqui se crean los
C valores locales de este ciclo para los numeros de los atomos nat.

      DO i=1,nprin
        nat(i) = natp(i)
      ENDDO

C Inicializacion del generador de numeros aleatorios con una semilla
C dependiente del reloj de la computadora para este ciclo

      if (ii.eq.1) then
        call init_random_seed()
      endif

C Se determinan las dimensiones del elipsoide central a partir del maximo
C valor de cada coordenada en la supermolecula ya creada en este
C ciclo.
C En el caso de solutos solidos la coordenada Z se dimensiona segun
C el dato "dimen" de entrada

40    if (qsolid1.or.qsolid2) then

C aqpm(j) es la coordenada j de menor valor de la molecula de soluto solido
C aqpp(j) es la coordenada j de mayor valor de la molecula de soluto solido
C d(j) es el tamanyo maximo en la coordenada j 
C Caso de solidos con el solvente en una cara del plano XY del soluto
C DIMEN da la distancia Z fija del cubo hacia el plano XY (SOLID1)

        do j=1,2
          aap = aa(j,1)
          aam = aa(j,1)
          do i=2,nprin
            if (aa(j,i).gt.aap) aap = aa(j,i)
            if (aa(j,i).le.aam) aam = aa(j,i)
          enddo
          aqpm(j) = aam
          aqpp(j) = aap
          d(j) = abs(aap - aam)
        enddo
        if (qsolid1) then
          d(3) = dimen
        else
          d(3) = 2.d0*dimen
        endif
      else

C Caso de soluciones
C Creacion de un elipsoide de dimensiones d(j) con el agregado
C molecular central
C El centro espacial de la molecula se situa al centro de coordenadas 

        do j=1,3
          aaa = aa(j,1)
          aap = aaa
          aam = aaa
          do i=2,nprin
            if (aa(j,i).gt.aap) aap = aa(j,i) 
            if (aa(j,i).le.aam) aam = aa(j,i)
          enddo
          aqpm(j) = aam
          aqpp(j) = aap
          dj = abs(aap - aam)
          do i=1,nprin
            aa(j,i) = aa(j,i) - (aqpm(j) + 0.5d0*dj)
          enddo
          if (.not.qroom) then
            d(j) = dj
          endif
        enddo
      endif

C Rotacion aleatoria del solvente

10    if (.not.qsolv) then
        do i=1,nsolv
          nats(i) = nath2o(i)
        enddo
        CALL ROTAT (nsolv,h2o,as)
      else
15      call RANDOM_NUMBER(Rand)
        ns = Rand*10.
        if (ns.eq.0) goto 15
        if ((ns.gt.nfsolv) .or.
     &     (nsrest(ns).le.0)) goto 15
        call OPENF (nsolv,nats,files(ns),titl1s)
        call ROTAT (nsolv,asl,as)
      endif

C Desplazamiento de las coordenadas de la molecula
C del solvente ya rotada a un origen aleatorio dentro del elipsoide
C y no mayor que la caja

      if (qsolid1) then
        DO j=1,2
          CALL RANDOM_NUMBER (Rand)
          ass = signo()*Rand*abs(aqpp(j) - aqpm(j))
          DO i=1,nsolv
            as(j,i) = as(j,i) - ass
            if (as(j,i).gt.aqpp(j)) goto 40
            if (as(j,i).lt.aqpm(j)) goto 40
          ENDDO
        ENDDO
        CALL RANDOM_NUMBER (Rand)
        ass = signo()*Rand*dimen
        DO i=1,nsolv
          as(3,i) = as(3,i) - ass
          if (as(3,i).lt.0) goto 40
        enddo
      elseif (qsolid2) then
        DO j=1,3
          CALL RANDOM_NUMBER (Rand)
          if (j.eq.3) then
            ass = apzmin + Rand*(apzmax - apzmin)
            DO i=1,nsolv
              as(3,i) = as(3,1) - ass
            if (as(3,i).gt.apzmax) goto 40
            if (as(3,i).lt.apzmin) goto 40
            ENDDO
          else
            ass = signo()*Rand*d(j)
            DO i=1,nsolv
              as(j,i) = as(j,i) - ass
            ENDDO
          endif
          DO i=1,nsolv
            as(j,i) = as(j,i) - ass
            if (as(3,i).gt.apzmax) goto 40
          ENDDO
        ENDDO
      else      
        DO j=1,3
          CALL RANDOM_NUMBER (Rand)
*          if (qroom) then
             ass = signo()*(Rand*d(j))
*          else
*            ass = signo()*(aqpm(j)+Rand*d(j))
*          endif
          DO i=1,nsolv
            as(j,i) = as(j,i) - ass
            if (abs(as(j,i)).gt.dimen) goto 40
          ENDDO
        ENDDO
      endif

C Esta posicion es aceptada si no hay atomos del soluto y el solvente
C mas proximos entre si que los parametros de van der Waals

      CALL TESTVDW (ntot, qroom, aa, nat, nsolv, as, nats, *10)
      if (qsolv) nsrest(ns) = nsrest(ns) - 1

C Adicion de la conformacion generada al sistema

      do i=1,nsolv
        nat(ntot+i) = nats(i)
        do j=1,3
          aa(j,ntot+i) = as(j,i)
        enddo
      enddo
      ntot = ntot + nsolv
      nn = nn + 1
      write (7,*) ' Solvent molecule:',nn,' is added ...'
      if (nn.lt.nlim) goto 10
C
      open (6,file=trim(sxyz),status='unknown')
      if (ii.gt.1) call appe (6,1)
      write (6,'(i5)') ntot
*      write (6,'(a)') title1
      write (title3,'(a,i5,a,i3,a)')
     &' CELL:',ii,'. Central system with',nlim,' distributed molecules.'
      write (6,'(a)') title3
      do i=1,ntot
        write (6,'(i4,3f12.5)') nat(i), (aa(j,i),j=1,3)
      enddo
      close (6)
      
C Creacion eventual de la entrada para ORCA

      if (qorca) then
        open (9,file=trim(sinp),status='unknown')
        if (ii.gt.1) then
          call appe (9,1)
          write (9,'(a)') '$new_job'
        endif
        write (9,'(a,a)') '# ', title1
        write (title3,'(a,i5,a,i3,a)')
     &' CELL:',ii,'. Central system with',nlim,' distributed molecules.'
        write (9,'(a,a)') '# ', title3
        do icard=1,ncard
          write (9,'(a)') cadenas(icard)
        enddo
        write (9,'(a)') '* xyz 0  1' 
        do i=1,ntot
          write (9,'(a2,3(f12.5))')
     &           ISYMB(nat(i)),(aa(j,i),j=1,3)
        enddo
        write (9,'(a)') '*' 
        close (9)
      endif
      
C Creacion eventual de la entrada para Gaussian

      if (qgau) then
        open (8,file=trim(sgjf),status='unknown')
        if (ii.gt.1) then
          call appe (8,1)
          write (8,'(/a)') '--Link1--'
        endif
        do icard=1,ncard
          write (8,'(a)') cadenas(icard)
        enddo
        write (8,'(/a)') title1
        write (title3,'(a,i5,a,i3,a)')
     &' CELL:',ii,'. Central system with',nlim,' distributed molecules.'
        write (8,'(a)') title3
        write (8,'(/a)') ' 0  1' 
        do i=1,ntot
          write (8,'(a2,3(f12.5))')
     &           ISYMB(nat(i)),(aa(j,i),j=1,3)
        enddo
        close (8)
      endif
      return
      END

      REAL*8 FUNCTION signo()
C Generador aleatorio de signos positivos y negativos
      CALL RANDOM_NUMBER( Rand )
      iii=int(100.*Rand)
      signo = (-1.d0)**iii
      RETURN
      END

      subroutine ROTAT (n, ai, ao)
      implicit real*8 (a-h,o,p,r-z)
      parameter (PI=3.141596d0, DOS=2.d0)
      real*8 ai(3,*), ao(3,*), t(3,3)
      real*4 Rand

      call RANDOM_NUMBER (Rand)
      theta = Rand*PI
      call RANDOM_NUMBER (Rand)
      beta = DOS*Rand*PI
      call RANDOM_NUMBER (Rand)
      phi = DOS*Rand*PI
      cosb = cos(beta)
      cosp = cos(phi)
      cost = cos(theta)
      sinb = sin(beta)
      sinp = sin(phi)
      sint = sin(theta)
      t(1,1) = cosb*cosp - cost*sinb*sinp
      t(1,2) = -cosb*sinp - cost*sinb*cosp
      t(1,3) = sint*sinb
      t(2,1) = sinb*cosp + cost*cosb*sinp
      t(2,2) = -sinb*sinp + cost*cosb*cosp
      t(2,3) = -sint*cosb
      t(3,1) = sint*sinp
      t(3,2) = sint*cosp
      t(3,3) = cost
      do i=1,n
          ao(1,i) = t(1,1)*ai(1,i) + t(1,2)*ai(2,i) + t(1,3)*ai(3,i)
          ao(2,i) = t(2,1)*ai(1,i) + t(2,2)*ai(2,i) + t(2,3)*ai(3,i)
          ao(3,i) = t(3,1)*ai(1,i) + t(3,2)*ai(2,i) + t(3,3)*ai(3,i)
      enddo
      return
      end

      SUBROUTINE TESTVDW (np, qroom, aa, natp, ns, as, nats, *)
      implicit real*8 (a-h,o,p,r-z)
      integer*4 natp(np), nats(ns)
      logical*2 qroom
      real*8
     & aa(3,*), as(3,*),
     & s(107)  /1.447d0,1.48d0,3*1.8d0,1.87d0,2*1.54d0,1.6d0,1.54d0,
     &          4*2.1d0,2.05d0,2.d0,1.95d0,1.92d0,14*2.25d0,2.2d0,
     &          2.15d0,2.1d0,2.07d0,14*2.4d0,57*2.25d0/
      if (qroom) then
        fff = 0.5d0
      else
        fff = 0.75d0
      endif
      DO ip=1,np
        DO is=1,ns
          r = sqrt((aa(1,ip)-as(1,is))**2 +
     &             (aa(2,ip)-as(2,is))**2 +
     &             (aa(3,ip)-as(3,is))**2)
          rmin = fff*(s(natp(ip)) + s(nats(is)))
          if (r.lt.rmin) then
            write (7,*)
     &      ' One random configuration rejected because overlaping ...'
            return 1
          endif
        ENDDO
      ENDDO
      RETURN
      END

      Subroutine OPENF (nsinp,nats,filei,tit1)
      include 'include.inc'
      common /dvariables/ aa(3,NATMAX), as(3,NATSMAX), ap(3,NATPMAX),
     &                    asl(3,NATSMAX), d(3), aqpm(3), aqpp(3) 
      character*128 filei
*      character*4 exts
      character*81 tit1
      dimension nats(*)
      common /elements/ elemnt(107)
      character*2 elemnt
      COMMON /simbolos/ isymb
      character*2 ISYMB(107), symbol
      open (14,file=filei,status='old')
      read (14,*) nsinp
      read (14,'(80a1)') tit1
      do 400 i=1,nsinp
        read (14,*) symbol, (asl(j,i),j=1,3) 
        nats(i) = reada(symbol,1)
        if (nats(i).eq.0) then
          do jj=1,107
            if (symbol.eq.elemnt(jj) .or. symbol.eq.ISYMB(jj)) then
              nats(i) = jj
                goto 400
            endif
          enddo
        endif
400   continue      
      close (14)
c
      do j=1,3
        aqi = asl(j,1)
        aqp = aqi
        aqm = aqi
        do i=2,nsinp
          aqi = asl(j,i)
          if (aqi.gt.aqp) aqp = aqi
          if (aqi.le.aqm) aqm = aqi
        enddo
        dd = aqp - aqm
        do i=1,nprin
          asl(j,i) = asl(j,i) - HALF*dd
        enddo
      enddo
      return
      end

      SUBROUTINE init_random_seed()
      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))
      CALL SYSTEM_CLOCK(COUNT=clock)
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)
      DEALLOCATE(seed)
      END SUBROUTINE
