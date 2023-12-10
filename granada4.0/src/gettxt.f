      SUBROUTINE GETTXT (iunit, KOMENT, TITLE)
      COMMON /KEYWRD/ KEYWRD
      DIMENSION IS(3)
      CHARACTER KEYWRD*241, KOMENT*81, TITLE*81, CH*1, CH2*1,
     +  OLDKEY*80

      IS(1)=161
      IS(2)=81
      IS(3)=1
      KEYWRD=' '
      KOMENT='    NULL  '
      TITLE ='    NULL  '
      READ(IUNIT,'(A)',END=100,ERR=90)KEYWRD(:80)
      OLDKEY=KEYWRD
      CALL UPCASE(KEYWRD(1:80))
      IF(INDEX(KEYWRD(1:80),' +') .NE.0)THEN
C
C  READ SECOND KEYWORD LINE
C
         READ(iunit,'(A)',END=100,ERR=90)KEYWRD(81:160)
         OLDKEY=KEYWRD(81:160)
         CALL UPCASE(KEYWRD(81:160))
         IF(INDEX(KEYWRD(81:160),' +') .NE.0)THEN
C
C  READ THIRD KEYWORD LINE
C
            READ(iunit,'(A)',END=100,ERR=90)KEYWRD(161:240)
            CALL UPCASE(KEYWRD(161:240))
         ENDIF
C
C  READ TITLE LINE
C
         READ(iunit,'(A)',END=100,ERR=90)KOMENT,TITLE
      ELSEIF(INDEX(KEYWRD(:80),'&').NE.0)THEN
         READ(iunit,'(A)',END=100,ERR=90)KEYWRD(81:160)
         OLDKEY=KEYWRD(81:160)
         CALL UPCASE(KEYWRD(81:160))
         IF(INDEX(KEYWRD(81:160),'&').NE.0)THEN
            READ(iunit,'(A)',END=100,ERR=90)KEYWRD(161:240)
         ELSE
            READ(iunit,'(A)',END=100,ERR=90)TITLE
         ENDIF
      ELSE
         READ(iunit,'(A)',END=100,ERR=90)KOMENT,TITLE
      ENDIF
   50 DO 80 J=1,3
         IF(KEYWRD(IS(J):IS(J)) .NE. ' ') THEN
            CH=KEYWRD(IS(J):IS(J))
            KEYWRD(IS(J):IS(J))=' '
            DO 60 I=IS(J)+1,239
               CH2=KEYWRD(I:I)
               KEYWRD(I:I)=CH
               CH=CH2
               IF(KEYWRD(I+1:I+2) .EQ. '  ') THEN
                  KEYWRD(I+1:I+1)=CH
                  GOTO 70
               ENDIF
   60       CONTINUE
            WRITE(6,'(A,I2,A)')' LINE',J,' OF KEYWORDS DOES NOT HAVE ENO
     1UGH'
            WRITE(6,'(A)')' SPACES FOR PARSING.  PLEASE CORRECT LINE.'
            STOP
   70       CONTINUE
         ENDIF
   80 CONTINUE
      RETURN
   90 WRITE(6,'(A)')' ERROR IN READ OF FIRST THREE LINES'
  100 STOP
      END
