!-----------------------------------------------------------------------------------------
SUBROUTINE configinfo(num,wattyp,nwat,nbet0,nbet1,nbet2,nnwata)
USE coords
USE atoms
IMPLICIT NONE
INTEGER :: wattyp
INTEGER :: i, num, nwata, nnwata, nwat, nbet0, nbet1, nbet2
DOUBLE PRECISION, DIMENSION(1:3) :: rbeta1
   nwata  = 0
   nnwata = 0
   nwat   = 0
   nbet0  = 0
   nbet1  = 0
   nbet2  = 0
   IF (wattyp.EQ.0) THEN
      DO i=1,num
         IF (TRIM(resname(i)).EQ.'TIP3') THEN
            nwata  = nwata + 1
         END IF
         IF (TRIM(resname(i)).NE.'TIP3') THEN
            nnwata = nnwata + 1
         END IF
         IF (INT(beta(i)).EQ.0) THEN
            nbet0  = nbet0 + 1
         END IF
         IF (INT(beta(i)).EQ.1) THEN
            nbet1  = nbet1 + 1
         END IF
         IF (INT(beta(i)).EQ.2) THEN
            nbet2  = nbet2 + 1
         END IF
      END DO
      IF (MOD(nwata,3).NE.0) THEN
         WRITE(*,*) 'TIP3P : COUNTED ', nwata,' NUMBER OF WATER ATOMS WHICH'
         WRITE(*,*) 'IS NOT DIVISIBLE BY 3 ...'
         WRITE(*,*) 'SOMETHING IS SERIOUSLY WRONG.'
         STOP
      END IF
      nwat = nwata/3
   ELSE IF (wattyp.EQ.1) THEN
      DO i=1,num
         IF (TRIM(resname(i)).EQ.'TIP4') THEN
            nwata  = nwata + 1
         END IF
         IF (TRIM(resname(i)).NE.'TIP4') THEN
            nnwata = nnwata + 1
         END IF
         IF (INT(beta(i)).EQ.0) THEN
            nbet0  = nbet0 + 1
         END IF
         IF (INT(beta(i)).EQ.1) THEN
            nbet1  = nbet1 + 1
         END IF
         IF (INT(beta(i)).EQ.2) THEN
            nbet2  = nbet2 + 1
         END IF
      END DO
      IF (MOD(nwata,4).NE.0) THEN
         WRITE(*,*) 'TIP4P : COUNTED ', nwata,' NUMBER OF WATER ATOMS WHICH '
         WRITE(*,*) 'IS NOT DIVISIBLE BY 4 ...'
         WRITE(*,*) 'SOMETHING IS SERIOUSLY WRONG.'
         STOP
      END IF
      nwat = nwata/4
   ELSE IF (wattyp.EQ.2) THEN
      DO i=1,num
         ! total waters with beta equal to 2
         !IF (TRIM(resname(i)).EQ.'SPCE' .AND. INT(beta(i)).EQ.2) THEN
         IF (TRIM(resname(i)).EQ.'SPCE') THEN
            nwata  = nwata + 1
         END IF
         IF (TRIM(resname(i)).NE.'SPCE') THEN
            nnwata = nnwata + 1
         END IF
         IF (INT(beta(i)).EQ.0) THEN
            nbet0  = nbet0 + 1
         END IF
         IF (INT(beta(i)).EQ.1) THEN
            nbet1  = nbet1 + 1
         END IF
         IF (INT(beta(i)).EQ.2) THEN
            nbet2  = nbet2 + 1
         END IF
      END DO
      IF (MOD(nwata,3).NE.0) THEN
         WRITE(*,*) 'SPCE : COUNTED ', nwata,' NUMBER OF WATER ATOMS WHICH'
         WRITE(*,*) 'IS NOT DIVISIBLE BY 3 ...'
         WRITE(*,*) 'SOMETHING IS SERIOUSLY WRONG.'
         STOP
      END IF
      nwat = nwata/3
   ELSE
      WRITE(*,*) 'UNRECOGNIZED WATER TYPE ...'
   END IF
RETURN
END SUBROUTINE configinfo
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
