      subroutine input
      include 'include.inc'

C Entrada de datos para todas las configuraciones

      common /filenames/ files, exts
      common /options/ qroom, qgau, qorca, qmopac, qsolid1, qsolid2
      common /strings/ cadenas(6), title1, title3
      common /datos/      nlim, dimen, nconf,
     &                    nfsolv, fracm(5), ncavid,
     &                    h2o(3,3), nath2o(3), nh2o
      common /indices/    ninput, qsolv, ncatal, natop, ncard
      common /dvariables/ aa(3,NATMAX), as(3,NATSMAX), ap(3,NATPMAX),
     &                    asl(3,NATSMAX), d(3), aqpm(3), aqpp(3), fvdw
      common /ivariables/ nat(NATMAX), natp(NATPMAX), nats(NATSMAX),
     &                    ntot, nprin, nsolv, nl(9)
!      save /indices/,/dvariables/,/ivariables/,/strings/,/datos/
      common /solid/ apzmax, apzmin, dx
      common /simbolos/ ISYMB
      character*2 ISYMB(107)
      common /elements/ elemnt(107)
      character*2 elemnt
      character*2 symbol
      COMMON /KEYWRD/ KEYWRD
      CHARACTER KEYWRD*241
      character*128 filep, filesol, files(5)
      character*80 cadenas, title1, title3
      character*4 extcar /'.xyz'/, extmop /'.mop'/, ext
      integer*2 nuno
      common /solvated/ smop, sxyz, sinp, sgjf, slog
      character(len=128)::smop, sxyz, sinp, sgjf, slog
      character(len=10):: cvalue(10)
      real*8 value(10)

      open (4,file='input.mmh',status='old')

      qsolid1 = .false.
      qsolid2 = .false.
      qgau = .false.
      qorca = .false.
      qmopac = .false.
      dimen = CERO
      fvdw = CERO
      ncatal = 0
      ncard = 1

C Lectura de:
C Linea 1
C - numero de moleculas de solvente NLIM (imp. 10),
C - semi-lado de la caja (en Angstroms) DIMEN (imp. 5.),
C - numero de configuraciones a calcular NCONF (imp. 1)
C - numero eventual de solventes diferentes NFSOLV (hasta 5),
C - fraccion molar de cada solvente FRACM(5)
C Linea 2
C - palabras clave con opciones

      read (4,'(i5,f10.0,2i5,5f10.0)') nlim, dimen, nconf, nfsolv, fracm
      KEYWRD=' '
      READ(4,'(A)')KEYWRD
      CALL UPCASE(KEYWRD)
c qroom es .true. cuando las moléculas del solvente pueden estar
c situadas dentro del elipsoide del cluster central
      if (INDEX(KEYWRD,'EXCLUDE').ne.0) then
        qroom = .false.
      else
        qroom = .true.
      endif
      if (nlim.eq.0) nlim = 10
c qgau es .true. cuando se desea una salida para Gaussian del tipo gjf
c qmopac es .true. cuando se desea una salida para Mopac del tipo mop
c qorca es .true. cuando se desea una salida para Orca del tipo inp
      qmopac = (INDEX(KEYWRD,'MOPAC').ne.0)
      qgau = (INDEX(KEYWRD,'GAUSSIAN').ne.0)
      if (INDEX(KEYWRD,'GAUSSIAN1').ne.0) ncard =1
      if (INDEX(KEYWRD,'GAUSSIAN2').ne.0) ncard =2
      if (INDEX(KEYWRD,'GAUSSIAN3').ne.0) ncard =3
      if (INDEX(KEYWRD,'GAUSSIAN4').ne.0) ncard =4
      if (INDEX(KEYWRD,'GAUSSIAN5').ne.0) ncard =5
      if (INDEX(KEYWRD,'GAUSSIAN6').ne.0) ncard =6
      if (qgau) then
        call valstr (INDEX(KEYWRD,'GAUSSIAN'),KEYWRD,cvalue)
        read (cvalue(1),'(f10.0)') value(1)
        ncard = value(1)
      endif
!
      qorca = (INDEX(KEYWRD,'ORCA').ne.0)
      if (INDEX(KEYWRD,'ORCA1').ne.0) ncard =1
      if (INDEX(KEYWRD,'ORCA2').ne.0) ncard =2
      if (INDEX(KEYWRD,'ORCA3').ne.0) ncard =3
      if (INDEX(KEYWRD,'ORCA4').ne.0) ncard =4
      if (INDEX(KEYWRD,'ORCA5').ne.0) ncard =5
      if (INDEX(KEYWRD,'ORCA6').ne.0) ncard =6
      if (qorca) then
        call valstr (INDEX(KEYWRD,'ORCA'),KEYWRD,cvalue)
        read (cvalue(1),'(f10.0)') value(1)
        ncard = value(1)
      endif
!
      if (nconf.eq.0) nconf = 1
      qsolid1 = (INDEX(KEYWRD,'SOLID1').ne.0)
      if (INDEX(KEYWRD,'SOLID1_0.5').ne.0) then
        dx = 0.5d0
      elseif (INDEX(KEYWRD,'SOLID1_1.0').ne.0) then
        dx = 1.d0
      elseif (INDEX(KEYWRD,'SOLID1_1.5').ne.0) then
        dx = 1.5d0
      endif  
      if (qsolid1) then
        call valstr (INDEX(KEYWRD,'SOLID1'),KEYWRD,cvalue)
        read (cvalue(1),'(f10.0)') value(1)
        dx = value(1)
      endif

      qsolid2 = (INDEX(KEYWRD,'SOLID2').ne.0)
      if (INDEX(KEYWRD,'SOLID2_0.5').ne.0) then
        dx = 0.5d0
      elseif (INDEX(KEYWRD,'SOLID2_1.0').ne.0) then
        dx = 1.d0
      elseif (INDEX(KEYWRD,'SOLID2_1.5').ne.0) then
        dx = 1.5d0
      endif  
      if (qsolid2) then
        call valstr (INDEX(KEYWRD,'SOLID2'),KEYWRD,cvalue)
        read (cvalue(1),'(f10.0)') value(1)
        dx = value(1)
      endif

      qfvdw = (INDEX(KEYWRD,'FVDW').ne.0)
      if (qfvdw) then
        call valstr (INDEX(KEYWRD,'FVDW'),KEYWRD,cvalue)
        read (cvalue(1),'(f10.0)') value(1)
        fvdw = value(1)
      endif
      if (qroom .and. fvdw.eq.CERO) then
        fvdw = 0.5d0
      elseif (fvdw.eq.CERO) then
        fvdw = 0.75d0
      endif

      
      if (dimen.eq.CERO) dimen=5.d0
      do i=1,ncard
        read (4,'(a)') cadenas (i)
      enddo
      close (4)

C Deteccion de cadenas en la linea de comandos y lectura del
C nombre del fichero del soluto

      ncom = COMMAND_ARGUMENT_COUNT()
      if (ncom-1.ne.nfsolv) then
        write (*,*)
     & ' Non congruent command line parameters with the amount of reques
     &ted solvent'
        write (*,*) ' input files'
        stop
      endif
      nuno = 1
      call GETARG (nuno,filep)
      filesol = filep
      kk = 0
      call filen1 (kk,filep,ext)
      smop = filep(1:kk)//'_solvated.mop'
      sxyz = filep(1:kk)//'_solvated.xyz'
      sinp = filep(1:kk)//'_solvated.inp'
      sgjf = filep(1:kk)//'_solvated.gjf'
      slog = filep(1:kk)//'_solvated.log'
      
      open (5,file=filesol,status='old')
C Lectura del fichero de cartesianas del soluto dado en la linea de comandos
C El formato de los ficheros de cartesianas es el universal xyz
* Lectura de ficheros .XYZ      
*
* LECTURA DEL NUMERO DE ATOMOS "NA" EN EL FICHERO DE ENTRADA DE CARTESIANAS
      read (5,*) ninput
* ENCABEZAMIENTO DE CADA JUEGO DE DATOS A CALCULAR
* LECTURA Y ESCRITURA DEL TEXTO DE IDENTIFICACION IDENT DE ESTA CORRIDA
      read (5,'(a)') title1
* LECTURA DE LAS COORDENADAS CARTESIANAS EN ANGTROMS Y DE SUS NUMEROS
* ATOMICOS
      do 400 I=1,ninput
        read (5,*) symbol, (ap(j,i),j=1,3) 
          natp(i) = reada(symbol,1)
          if (natp(i).eq.0) then
            do jj=1,107
              if (symbol.eq.elemnt(jj) .or. symbol.eq.ISYMB(jj)) then
                natp(i) = jj
                goto 400
              endif
            enddo
          endif
400   continue
      close (5)

      if (qsolid2) then
c Formacion de la caja con doble soluto
      do i=1,ninput
          ap(3,i) = ap(3,i) - dimen
        enddo
        do i=1,ninput
          natp(ninput+i) = natp(i)
          ap(1,ninput+i) = ap(1,i)
          ap(2,ninput+i) = ap(2,i)
          ap(3,ninput+i) = -ap(3,i)
        enddo
        ninput = ninput + ninput
c Determinación de la máxima y mínima z
        apzmax = ap(3,1)
        apzmin = ap(3,1)
        do i=2,ninput
          if (ap(3,i).gt.apzmax) apzmax = ap(3,i)
          if (ap(3,i).lt.apzmin) apzmin = ap(3,i)
        enddo
      endif
      nprin = ninput

C Aqui se crean los valores iniciales para las coordenadas aa.
C Ubicacion del sistema con el centro de coordenadas en el
C centroide de la molecula central. aqp es la coordenada positiva
C mayor del eje j y aqm es la coordenada negativa menor del eje j

      if (qsolid2) then
        do j=1,2
          aqi = ap(j,1)
          aqp = aqi
          aqm = aqi
          do i=2,nprin
            aqi = ap(j,i)
            if (aqi.gt.aqp) aqp = aqi
            if (aqi.lt.aqm) aqm = aqi
          enddo
          dd = aqp - aqm
          do i=1,nprin
            aa(j,i) = ap(j,i) - HALF*dd
          enddo
        enddo
        do i=1,nprin
          aa(3,i) = ap(3,i)
        enddo
      else
        do j=1,3
          aqi = ap(j,1)
          aqp = aqi
          aqm = aqi
          do i=2,nprin
            aqi = ap(j,i)
            if (aqi.gt.aqp) aqp = aqi
            if (aqi.lt.aqm) aqm = aqi
          enddo
          dd = aqp - aqm
          do i=1,nprin
            aa(j,i) = ap(j,i) - HALF*dd
          enddo
        enddo
      endif
C Ahora entran los nombres de ficheros de solventes, si es que
C no se usa el agua implicita. Tambien se determina el numero nl(i)
C de moleculas de cada solvente que entrara, de acuerdo con sus frac-
C ciones molares de entrada respectivas. La suma de todos los
C nl(i), donde i es el numero de orden del solvente, debe dar
C nlim, o el numero total de moleculas de solvente.

c qsolv es .true. cuando el solvente es diferente del agua implicita

      qsolv = nfsolv.gt.0
      if (qsolv) then
        nntot = 0
        do nuno=2,nfsolv+1
          indice = nuno-1
          call GETARG(nuno,files(indice))
          nl(indice) = fracm(indice)*dble(nlim) + .5d0
          nntot = nntot + nl(indice)
        enddo
        ndif = nlim - nntot
        if (ndif.ne.0) then

c Si la normalizacion no conduce a que la suma de los nl sea
c nlim entonces se renormaliza

          nmax = 0
          do j=1,nfsolv
            if (nl(j).gt.nmax) then
              nmax = nl(j)
              jmax = j
            endif
          enddo
          nl(jmax) = nl(jmax) + ndif
          write (7,*) ' Molar fractions have been adjusted to real propo
     &rtions of solvent molecules'
          write (*,*) ' Molar fractions have been adjusted to real propo
     &rtions of solvent molecules'
          do j=1,nfsolv
            fracm(j) = dble(nl(j)) / dble(nlim)
          enddo
          write (7,*) ' New molar fractions are:', fracm
          write (*,*) ' New molar fractions are:', fracm
        endif
      endif
      return
      end
