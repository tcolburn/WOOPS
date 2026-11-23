      SUBROUTINE copyright
      IMPLICIT NONE
! ********************************************************************
! * This routine prints the copyright information.
! ********************************************************************
      ! Args
      CHARACTER(LEN=57) :: strTmp
! *******************************************************************

      WRITE(*,5008) '' 
      strTmp='+-------------------------------------------------------+'
      WRITE(*,5008) strTmp
      strTmp='| (c) 2015 Arizona Board of Regents acting for and on   |'
      WRITE(*,5008) strTmp
      strTmp='| behalf of Arizona State University ("ABOR").  All     |'
      WRITE(*,5008) strTmp
      strTmp='| rights reserved.   This software is protected by the  |'
      WRITE(*,5008) strTmp
      strTmp='| copyright laws of the United States and international |'
      WRITE(*,5008) strTmp
      strTmp='| treaties, may not be copied, reproduced, distributed, |'
      WRITE(*,5008) strTmp
      strTmp='| sold or offered for sale without the express written  |'
      WRITE(*,5008) strTmp
      strTmp='| permission of the copyright owner, and is only for    |'
      WRITE(*,5008) strTmp
      strTmp='| academic, non-commercial use.  No modifications or    |'
      WRITE(*,5008) strTmp
      strTmp='| derivations of this software are permitted under any  |'
      WRITE(*,5008) strTmp
      strTmp='| circumstances, and this software may not be incorp-   |'
      WRITE(*,5008) strTmp
      strTmp='| orated into any other software or other work.  ABOR   |'
      WRITE(*,5008) strTmp
      strTmp='| shall not be liable for any actions taken in reliance |'
      WRITE(*,5008) strTmp
      strTmp='| on use of this software.                              |'
      WRITE(*,5008) strTmp
      strTmp='+-------------------------------------------------------+'
      WRITE(*,5008) strTmp
      WRITE(*,5008) '' 

      RETURN

5008  FORMAT(57a)
      END SUBROUTINE copyright
