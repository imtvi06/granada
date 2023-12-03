      program GRANADA
c Version 3.0
c Luis A. Montero Cabrera, Universidad de La Habana, October 27, 2019
c Previous versions were called al GRANADA xxxx where xxxx is was
c the year of production
      include 'include.inc'
      character*80 cadenas
      character*80 title1, title3
      character*4 exts
      character*32 files(5)
      common /filenames/ files, exts
      common /options/ qroom, qgau, qorca, qmopac, qsolid1, qsolid2
      common /strings/ cadenas(6), title1, title3
      common /indices/    ninput, qsolv, ncatal, natop, ncard
      common /datos/      nlim, dimen, nconf,
     &                    nfsolv, fracm(5), ncavid,
     &                    h2o(3,3), nath2o(3), nh2o
      common /dvariables/ aa(3,NATMAX), as(3,NATSMAX), ap(3,NATPMAX),
     &                    asl(3,NATSMAX), d(3), aqpm(3), aqpp(3), fvdw
      common /ivariables/ nat(NATMAX), natp(NATPMAX), nats(NATSMAX),
     &                    ntot, nprin, nsolv, nl(9)
      common /solvated/ smop, sxyz, sinp, sgjf, slog
      character(len=128)::smop, sxyz, sinp, sgjf, slog
      
      call input
      open (7,file=trim(slog),status='unknown')

! Salida de los datos de entrada
 
      write (7,*) ' Granada MMH v.4.5 main log information'
      write (7,*) ' Released in Havana, October 2020'
      write (*,*) ' Granada MMH v.4.5 main log information'
      write (*,*) ' Released in Havana, October 2020'
      write (7,*) ' Maximum number of atoms:',NATMAX
      write (7,*) ' Maximum number of atoms in the solute:',NATPMAX
      write (7,*) ' Maximum number of atoms in the solvent:',NATSMAX
      write (7,*) ' *** Summary of input data:'
      write (*,*) ' SUMMARY OF INPUT DATA:'
      write (7,*) ' Solvent molecules:', nlim
      write (*,*) ' Solvent molecules:', nlim
      write (7,1000) ' Box dimension input (Angst.):', dimen
      write (*,1000) ' Box dimension input (Angst.):', dimen
      write (7,1000) ' Factor for rejection by VDW radii sums:', fvdw
      write (*,1000) ' Factor for rejection by VDW radii sums:', fvdw
      write (7,*) ' Requested configurations:', nconf
      write (*,*) ' Requested configurations:', nconf
      if (nfsolv.ne.CERO) then
        write (7,*) ' Different kinds of solvents:', nfsolv
        write (7,1000) ' Input molar fractions:', fracm
        write (*,*) ' Different kinds of solvents:', nfsolv
        write (*,1000) ' Input molar fractions:', fracm
      endif
1000  format (a,f10.4)
      if (qsolid1) write (7,'(a)') 
     &   ' The solute is a solid layer in the XY plane'
      if (qsolid2) write (7,'(a,g6.2,a)') 
     &   ' The solute is a mirrored bilayer separated by',dimen*2.d0,
     &   ' A in the Z dimension'

      do ii=1,nconf
        write (7,'(a,i5.5)') ' Cell No. ', ii
        call dresdenr (ii) 
        if (qmopac) call dresdenc (ii)
      enddo
      close (7)
      stop
      end

      subroutine appe (ifile,iapp4)
      implicit integer*2 (i-n)
      integer*4 ifile, iapp4
      character*1 a
      character*11 ans
*      iff = ifile
*      inquire (iff,form=ans)
*      if (ans.eq.'unformatted') then
*        write(*,'(/a/a/)') ' ERROR: file is sequential and unformatted',
*     &                     ' it can not be appended by appe...'
*        stop
*      endif
      if (iapp4.ne.0) then
1       read (ifile,'(a)',end=2) a
        go to 1
      endif
      return
2     backspace ifile
      return
      end
      
      subroutine filen1 (kk,fname,ext)
      implicit integer*2 (i-n)
      integer*4 kk
c     default extension for dos filenames
c     kk=0 returns fname without extension and kk becomes the
c          number of characters of the filename without extension
c       =1 returns fname with the new extension ext
      character*1 fname(128),ext(4)
      do 10 i=1,128
      if (fname(i).eq.'.') then
          if (kk.eq.0) then
            kk = i - 1
            go to 50
          endif
          go to 20
      endif
      if (fname(i).eq.' ') go to 20
10    continue
20    do  j=1,4
        fname(i-1+j) = ext(j)
      enddo
50    return
      end


      block data
      include 'include.inc'
      common /datos/      nlim, dimen, nconf,
     &                    nfsolv, fracm(5), ncavid,
     &                    h2o(3,3), nath2o(3), nh2o
      common /simbolos/ isymb
      character*2 ISYMB(107)
      common /elements/ elemnt(107)
      character*2 elemnt
      DATA ISYMB /' H','He','Li','Be',' B',' C',' N',' O',' F','Ne','Na'
     .            ,'Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc',
     .            'Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge'
     .            ,'As','Se','Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo',
     .            'Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I'
     .            ,'Xe','Cs','Ba','La',14*' *','Hf','Ta',' W','Re','Os',
     .            'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr'
     .            ,'Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf'
     .            ,'XX','Fm','Md','Cb','++','+','--','-','TV'/
* ELEMNT es el arreglo de simbolos atomicos para las comparaciones
* de caracteres en mayúsculas
      DATA ELEMNT/'H','HE',
     1 'LI','BE','B','C','N','O','F','NE',
     2 'NA','MG','AL','SI','P','S','CL','AR',
     3 'K','CA','SC','TI','V','CR','MN','FE','CO','NI','CU',
     4 'ZN','GA','GE','AS','SE','BR','KR',
     5 'RB','SR','Y','ZR','NB','MO','TC','RU','RH','PD','AG',
     6 'CD','IN','SN','SB','TE','I','XE',
     7 'CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB','DY',
     8 'HO','ER','TM','YB','LU','HF','TA','W','RE','OS','IR','PT',
     9 'AU','HG','TL','PB','BI','PO','AT','RN',
     1 'FR','RA','AC','TH','PA','U','NP','PU','AM','CM','BK','CF','XX',
     2 'FM','MD','CB','++','+','--','-','TV'/
 
      data
     &    h2o             /0.d0,0.d0,1.d0,
     &                    .95097d0,0.d0,1.d0,
     &                    1.23998d0,.90599d0,1.d0/,
     &    nath2o          /1,8,1/,
     &    nh2o            /3/,
     &    nlim            /0/,
     &    nfsolv          /0/
      end
