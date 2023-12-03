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
!      print *, nt
      do i=nt,241
        if (c(i).eq."(") then
          qparentesis = .TRUE.
          cycle
        elseif (c(i).eq.")") then
          exit
        elseif (i.eq.(nt+10) .and. .not.qparentesis) then
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
!          print *, i,j-1,k, c(i)
        endif
      enddo
      end
