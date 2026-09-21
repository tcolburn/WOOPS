!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
SUBROUTINE alStructArrays
USE atoms
IMPLICIT NONE
   ALLOCATE(resid(1:natom),resname(1:natom))
   ALLOCATE(aindex(1:natom),atomname(1:natom),atomtype(1:natom))
   ALLOCATE(charge(1:natom),mass(1:natom))
   ALLOCATE(pdbresid(1:natom),pdbresname(1:natom))
   ALLOCATE(pdbaindex(1:natom),pdbatomname(1:natom))
   ALLOCATE(pdbaindexC(1:natom))
   ALLOCATE(coor(1:natom),beta(1:natom),segname(1:natom))
   ALLOCATE(xini(1:natom),yini(1:natom),zini(1:natom))
   ALLOCATE(atomname2(1:natom))
RETURN
END SUBROUTINE alStructArrays
!-------------------------------------------------------------------------------
SUBROUTINE deStructArrays
USE atoms
IMPLICIT NONE
   DEALLOCATE(resid,resname)
   DEALLOCATE(aindex,atomname,atomtype)
   DEALLOCATE(charge,mass)
   DEALLOCATE(pdbresid,pdbresname)
   DEALLOCATE(pdbaindex,pdbatomname)
   DEALLOCATE(pdbaindexC)
   DEALLOCATE(coor,beta,segname)
   DEALLOCATE(xini,yini,zini)
   DEALLOCATE(atomname2)
RETURN
END SUBROUTINE deStructArrays
!-------------------------------------------------------------------------------
SUBROUTINE alCoordArrays(n)
USE coords
USE atoms, ONLY : natom
IMPLICIT NONE
INTEGER :: n
! Streaming implementation: only the current frame is kept in memory.
! The argument n is retained for interface compatibility.
  ALLOCATE(rx(1:natom),ry(1:natom),rz(1:natom))
RETURN
END SUBROUTINE alCoordArrays
!-------------------------------------------------------------------------------
SUBROUTINE deCoordArrays
USE coords
IMPLICIT NONE
  IF (ALLOCATED(rxf)) DEALLOCATE(rxf)
  IF (ALLOCATED(ryf)) DEALLOCATE(ryf)
  IF (ALLOCATED(rzf)) DEALLOCATE(rzf)
  IF (ALLOCATED(box)) DEALLOCATE(box)
  IF (ALLOCATED(rx)) DEALLOCATE(rx)
  IF (ALLOCATED(ry)) DEALLOCATE(ry)
  IF (ALLOCATED(rz)) DEALLOCATE(rz)
RETURN
END SUBROUTINE deCoordArrays
!-------------------------------------------------------------------------------
!SUBROUTINE alWatArrays(num)
!USE water
!IMPLICIT NONE
!INTEGER :: num, num2
!   num2 = 2*num
!   ALLOCATE(rxw(1:num),ryw(1:num),rzw(1:num))
!   ALLOCATE(exw(1:num),eyw(1:num),ezw(1:num))
!   ALLOCATE(rxOx(1:num),ryOx(1:num),rzOx(1:num))
!   ALLOCATE(rxHy(1:num2),ryHy(1:num2),rzHy(1:num2))
!RETURN
!END SUBROUTINE alWatArrays
!-------------------------------------------------------------------------------
!SUBROUTINE deWatArrays
!USE water
!IMPLICIT NONE
!   DEALLOCATE(rxw,ryw,rzw)
!   DEALLOCATE(exw,eyw,ezw)
!   DEALLOCATE(rxOx,ryOx,rzOx)
!   DEALLOCATE(rxHy,ryHy,rzHy)
!RETURN
!END SUBROUTINE deWatArrays
!-------------------------------------------------------------------------------
!SUBROUTINE alTip4pWatArrays(num)
!USE water
!IMPLICIT NONE
!INTEGER :: num, num2
!   num2 = 2*num
!   ALLOCATE(rxw(1:num),ryw(1:num),rzw(1:num))
!   ALLOCATE(exw(1:num),eyw(1:num),ezw(1:num))
!   ALLOCATE(rxOx(1:num),ryOx(1:num),rzOx(1:num))
!   ALLOCATE(rxOm(1:num),ryOm(1:num),rzOm(1:num))
!   ALLOCATE(rxHy(1:num2),ryHy(1:num2),rzHy(1:num2))
!RETURN
!END SUBROUTINE alTip4pWatArrays
!-------------------------------------------------------------------------------
!SUBROUTINE deTip4pWatArrays
!USE water
!IMPLICIT NONE
!   DEALLOCATE(rxw,ryw,rzw)
!   DEALLOCATE(exw,eyw,ezw)
!   DEALLOCATE(rxOx,ryOx,rzOx)
!   DEALLOCATE(rxOm,ryOm,rzOm)
!   DEALLOCATE(rxHy,ryHy,rzHy)
!RETURN
!END SUBROUTINE deTip4pWatArrays
!-------------------------------------------------------------------------------
SUBROUTINE alBeta1Arrays(num)
USE beta1
IMPLICIT NONE
INTEGER :: num
   ALLOCATE(rxbt1(1:num),rybt1(1:num),rzbt1(1:num))
   ALLOCATE(chrgbt1(1:num),massbt1(1:num),rnamebt1(1:num),ridbt1(1:num))
   ALLOCATE(aidbt1(1:num),anamebt1(1:num),atypbt1(1:num),segbt1(1:num),atmnm2bt1(1:num))
   aidbt1(:) = 0
RETURN
END SUBROUTINE alBeta1Arrays
!-------------------------------------------------------------------------------
SUBROUTINE deBeta1Arrays
USE beta1
IMPLICIT NONE
   DEALLOCATE(rxbt1,rybt1,rzbt1)
   DEALLOCATE(chrgbt1,massbt1,rnamebt1,ridbt1)
   DEALLOCATE(aidbt1,anamebt1,atypbt1,segbt1,atmnm2bt1)
RETURN
END SUBROUTINE deBeta1Arrays
!-------------------------------------------------------------------------------
SUBROUTINE alBeta2Arrays(num)
USE beta2
IMPLICIT NONE
INTEGER :: num
   ALLOCATE(rxbt2(1:num),rybt2(1:num),rzbt2(1:num))
   ALLOCATE(chrgbt2(1:num),massbt2(1:num),rnamebt2(1:num),ridbt2(1:num))
   ALLOCATE(aidbt2(1:num),anamebt2(1:num),atypbt2(1:num),segbt2(1:num),atmnm2bt2(1:num))
   aidbt2(:) = 0
RETURN
END SUBROUTINE alBeta2Arrays
!-------------------------------------------------------------------------------
SUBROUTINE deBeta2Arrays
USE beta2
IMPLICIT NONE
   DEALLOCATE(rxbt2,rybt2,rzbt2)
   DEALLOCATE(chrgbt2,massbt2,rnamebt2,ridbt2)
   DEALLOCATE(aidbt2,anamebt2,atypbt2,segbt2,atmnm2bt2)
RETURN
END SUBROUTINE deBeta2Arrays
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
