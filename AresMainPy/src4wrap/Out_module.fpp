# 1 "Out_module.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "Out_module.f90"
MODULE out_module
   IMPLICIT NONE
CONTAINS
!------------------------list--------------------
   SUBROUTINE KSout_list()
      USE parameters , ONLY: Lpbc
      IMPLICIT NONE

      if(Lpbc)then
         call KSout_list_per()
      else
         call KSout_list_iso()
      endif
   ENDSUBROUTINE KSout_list
   !------------------- divide line ---------------------

   SUBROUTINE KSout_list_per()
      USE pspot_module , ONLY : max_nproj
      USE struct_module , ONLY : natom,energy
      USE energy_module
      USE grid_module , ONLY : n1,n2,n3,grid
      USE potential_module , ONLY : sumrhoS
      USE m_time_evaluate, ONLY: memory_sum,memory_free

      USE smpi_math_module , ONLY : parallel

      IMPLICIT NONE
      REAL(DP) :: rho(n1,n2,n3)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      if(parallel%isroot)then

      call memory_sum('KS_out_list_local',real(size(rho),DP)*DP)
      !ground data
         WRITE(6,*)'*******************OUTPUT DATA***********************'

         WRITE(6,10)'*Fermi level     ->',Efm*hart2ev,'eV'
         WRITE(6,*)'-----------------------------------------------------'
         WRITE(6,10)'*Max Last band w ->',WmaxL,'eV'
         WRITE(6,10)'*Min Last band w ->',WminL,'eV'
         WRITE(6,*)'-----------------------------------------------------'
         WRITE(6,10)'*Band energy     ->',Eband*hart2ev,'eV'
         !WRITE(6,10)'*Kinetic energy  ->',Ekine*hart2ev,'eV'
         WRITE(6,10)'*Hartree energy  ->',Ehart*hart2ev,'eV'
         WRITE(6,10)'*XC energy       ->',Exc*hart2ev,'eV'
         WRITE(6,10)'*I-E energy(loc) ->',Eie*hart2ev,'eV'
         !IF(max_nproj>0)THEN
         !   WRITE(6,10)'*I-E energy(nloc)->',Eienl*hart2ev,'eV'
         !ENDIF
         WRITE(6,10)'*Ewald energy    ->',Eewald*hart2ev,'eV'
         WRITE(6,*)'-----------------------------------------------------'
         WRITE(6,10)'*Total energy    ->',Etot*hart2eV,'eV'
         WRITE(6,10)'*Free energy     ->',FE*hart2eV,'eV'
         WRITE(6,10)'*0K Total energy ->',FE0*hart2eV,'eV'
         WRITE(6,*)'-----------------------------------------------------'
         WRITE(6,10)'*Energy per atom ->',FE0*hart2eV/natom,'eV'
         WRITE(6,*)'-----------------------------------------------------'
         WRITE(6,*)'*******************OUTPUT DATA***********************'

      endif

10    FORMAT(A20,1X,' ',1X,SP,F25.7,A3)
      !> For CG
      energy(1)=Etot
      !ground density
      CALL sumrhoS(grid%rhoS,rho)

      if(parallel%isroot)then

      OPEN(1110,FILE='ares.chg')
         WRITE(1110,*) n1,n2,n3
         WRITE(1110,*) rho
      CLOSE(1110)

      endif

      call memory_free('KS_out_list_local',real(size(rho),DP)*DP)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE KSout_list_per
   !------------------- divide line ---------------------
   SUBROUTINE KSout_list_iso()
      USE pspot_module , ONLY : max_nproj
      USE struct_module , ONLY : natom,energy
      USE energy_module



      USE grid_module , ONLY : n1,n2,n3,grid,rho_calc,parallel_s2g

      USE potential_module , ONLY : sumrhoS

      USE smpi_math_module , ONLY : parallel
      USE parameters , ONLY: Nspin

      USE m_time_evaluate, ONLY: memory_sum,memory_free
      IMPLICIT NONE



      REAL(DP) :: rho(n1,n2,n3)
      INTEGER(I4B) :: i

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call memory_sum('KS_out_list_local',real(size(rho),DP)*DP)
      !ground data

      if(parallel%isroot)then

         WRITE(6,*)'*******************OUTPUT DATA***********************'

         WRITE(6,10)'*Fermi level     ->',Efm*hart2ev,'eV'
         WRITE(6,*)'-----------------------------------------------------'
         WRITE(6,10)'*Max Last band w ->',WmaxL,'eV'
         WRITE(6,10)'*Min Last band w ->',WminL,'eV'
         WRITE(6,*)'-----------------------------------------------------'
         WRITE(6,10)'*Band energy     ->',Eband*hart2ev,'eV'
         !WRITE(6,10)'*Kinetic energy  ->',Ekine*hart2ev,'eV'
         WRITE(6,10)'*Hartree energy  ->',Ehart*hart2ev,'eV'
         WRITE(6,10)'*XC energy       ->',Exc*hart2ev,'eV'
         WRITE(6,10)'*I-E energy(loc) ->',Eie*hart2ev,'eV'
         !IF(max_nproj>0)THEN
         !   WRITE(6,10)'*I-E energy(nloc)->',Eienl*hart2ev,'eV'
         !ENDIF
         WRITE(6,10)'*Ewald energy    ->',Eewald*hart2ev,'eV'
         WRITE(6,*)'-----------------------------------------------------'
         WRITE(6,10)'*Total energy    ->',Etot*hart2eV,'eV'
         WRITE(6,10)'*Free energy     ->',FE*hart2eV,'eV'
         WRITE(6,10)'*0K Total energy ->',FE0*hart2eV,'eV'
         WRITE(6,*)'-----------------------------------------------------'
         WRITE(6,10)'*Energy per atom ->',FE0*hart2eV/natom,'eV'
         WRITE(6,*)'-----------------------------------------------------'
         WRITE(6,*)'*******************OUTPUT DATA***********************'

      endif

10    FORMAT(A20,1X,' ',1X,SP,F25.7,A3)
      !> For CG
      energy(1)=Etot
      !ground density



      do i =1,Nspin
         CALL parallel_s2g(rho_calc%OneDSphere(:,i),grid%rhoS(:,:,:,i))
      enddo
      rho(:,:,:)=0.d0
      DO i=1,NSPIN
         rho(:,:,:)=rho(:,:,:) + grid%rhoS(:,:,:,i)
      ENDDO


      if(parallel%isroot)then

      OPEN(1110,FILE='ares.chg')
         WRITE(1110,*) n1,n2,n3
         WRITE(1110,*) rho
      CLOSE(1110)

      endif

      call memory_free('KS_out_list_local',real(size(rho),DP)*DP)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE KSout_list_iso


ENDMODULE out_module
