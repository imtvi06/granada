      subroutine readint
     &  (iunit, NATOMS, coord,natdum, title1, title2)
      include 'include.inc'
c Los siguientes commons son para casos especiales en la subrutina GMETRY
c Si el numero atomico del ultimo atomo es 107, se procesa como un solido y
c es preciso evaluar NUMCAL.Ver la subrutina GMETRY.
      COMMON /KEYWRD/ KEYWRD
      CHARACTER KEYWRD*241, TITLE1*81, TITLE2*81
      common /geos/ geo(3,NATMAX)
      DIMENSION natdum(NATMAX),coord(3,*),
     &          NA(NATMAX),NB(NATMAX),NC(NATMAX),LOPT(3,NATMAX)
      step = CERO
      call gettxt (iunit, title1, title2)
      call getgeo (iunit,natdum,geo,lopt,na,nb,nc,NATOMS,qaint)
      call geome (NATOMS, coord, na, nb, nc)
      
      return
      end