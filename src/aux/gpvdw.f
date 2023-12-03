      PROGRAM gpvdw
      CHARACTER (len=241) :: cadena
      CHARACTER (len=10) :: value(10)
      REAL*8 v(10)
      cadena = 'Esta es una prueba de VDW(3.487,35.8),test(45.6)'
      print *, cadena
      nt = 23
      call valstr (nt, cadena, value)
      print *, value(1), value(2)
      do i=1,10
        read(value(i),'(f10.0)') v(i)
        write(*,'(f10.5)') v(i)
      enddo
      nt = 39
      call valstr (nt, cadena, value)
      print *, value(1)
      do i=1,10
        read(value(i),'(f10.0)') v(i)
        write(*,'(f10.5)') v(i)
      enddo
      END PROGRAM

      subroutine valstr (nt, c, cvalue)
      CHARACTER :: c(241), cvalue(10,10)
      LOGICAL :: qparentesis, qcoma
      do i = 1,10
        do j = 1,10
          cvalue(i,j) = ' '
        enddo
      enddo
      qparentesis = .FALSE.
      qcoma = .FALSE.
      j = 1
      k = 1
      print *, nt
      do i=nt,241
        if (c(i).eq."(") then
          qparentesis = .TRUE.
          cycle
        elseif (c(i).eq.")") then
          exit
        elseif (i.eq.(nt+20) .and. .not.qparentesis) then
          return
        endif
        if (c(i).eq.",") then
          qcoma = .TRUE.
          j = 1
          k = k + 1
          cycle
        endif   
        if (qparentesis .or. qcoma) then
          cvalue(j,k) = c(i)
          j = j + 1
          print *, i,j-1,k, c(i)
        endif
      enddo
      end
