      subroutine DRESDENC (ii)
      include 'include.inc'
      common /indices/ ninput, qsolv, ncatal, natop, ncard
      common /strings/ cadenas(6), title1, title3
      common /dvariables/ aa(3,NATMAX), as(3,NATSMAX), ap(3,NATPMAX),
     &                    asl(3,NATSMAX), d(3), aqpm(3)
      common /ivariables/ nat(NATMAX), natp(NATPMAX), nats(NATSMAX),
     &                    ntot, nprin, nsolv, nl(9)
      character*80 cadenas,title1, title3
      real*8 geo(3,NATMAX)
      integer*4
     & na(NATMAX) /NATMAX*0/,
     & nb(NATMAX) /NATMAX*0/,
     & nc(NATMAX) /NATMAX*0/
      COMMON /simbolos/ isymb
      character*2 ISYMB(106)
      common /solvated/ smop, sxyz, sinp, sgjf, slog
      character(len=128)::smop, sxyz, sinp, sgjf, slog
      
      open (6,file=trim(smop),status='unknown')
      if (ii.gt.1) call appe (6,1)
      noutput = ntot
      na(2) = 1
      na(3) = 2
      nb(3) = 1
      if (ncatal.eq.0) then
        do i=4,noutput
          na(i) = i - 1
          nb(i) = i - 2
          nc(i) = i - 3
        enddo
      else
        do i=4,natop
          na(i) = i - 1
          nb(i) = i - 2
          nc(i) = i - 3
        enddo
        do i=natop+1,noutput
          na(i) = 3
          nb(i) = 2
          nc(i) = 1
        enddo
      endif
      degree = 180.d0/3.141569d0         
      call xyzgeo (aa,noutput,na,nb,nc,degree,geo)
      write (6,'(1x,a)') cadenas(1)
      write (6,'(1x,a)') title1
      write (6,'(1x,a)') title3
      write (6,'(1x,a2)') ISYMB(nat(1))
      if (ncatal.eq.0) then
        write (6,'(1x,a2,f10.5,i3)')
     &    ISYMB(nat(2)),geo(1,2),1
        write (6,'(1x,a2,2(f10.5,i3))')
     &    ISYMB(nat(3)),geo(1,3),1,geo(2,3),1
        do i=4,noutput
        write (6,'(1x,a2,3(f10.5,i3),3i5)')
     &    ISYMB(nat(i)),geo(1,i),1,geo(2,i),1,geo(3,i),1,
     &    na(i),nb(i),nc(i)
        enddo
      else
        if (natop.ge.2) then
          iop2 = 0
        else
          iop2 = 1
        endif    
        write (6,'(1x,a2,f10.5,i3)')
     &    ISYMB(nat(2)),geo(1,2),iop2
        if (natop.ge.3) then
          iop3 = 0
        else
          iop3 = 1
        endif    
        write (6,'(1x,a2,2(f10.5,i3))')
     &    ISYMB(nat(3)),geo(1,3),iop3,geo(2,3),iop3
        if (natop.ge.4) then
          do i=4,natop
            write (6,'(1x,a2,3(f10.5,i3),3i5)')
     &        ISYMB(nat(i)),geo(1,i),0,geo(2,i),0,geo(3,i),0,
     &        na(i),nb(i),nc(i)
          enddo
          do i=natop+1,noutput
            write (6,'(1x,a2,3(f10.5,i3),3i5)')
     &        ISYMB(nat(i)),geo(1,i),1,geo(2,i),1,geo(3,i),1,
     &        na(i),nb(i),nc(i)
          enddo
        else
          do i=4,noutput
            write (6,'(1x,a2,3(f10.5,i3),3i5)')
     &        ISYMB(nat(i)),geo(1,i),1,geo(2,i),1,geo(3,i),1,
     &        na(i),nb(i),nc(i)
          enddo
        endif    
      endif  
      write (6,'(a)') ' 0'
      close (6)
      return
      end
      
