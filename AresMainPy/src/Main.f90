PROGRAM ARES
   !##########################BEAST WAR#############################!
   !*BEAST   : real space Kohn-Sham DFT Software                    !
   !*Author  : Qiang Xu                                             !
   !*E-mail  : xq@calypso.cn                                        !
   !*Date    : 2017/07/01                                           !
   !*Diag H  : ARPACK                                               !
   !################################################################!
  !USE omp_lib
   USE constants
   USE read_module
   USE parameters  , ONLY : ISTART,Nssp,outfile,Lpbc,isp
   USE End_module  , ONLY : destroy_beast
   USE forcestress_module , ONLY : cal_force_stress
   USE relax_module , ONLY : initialize_relax,relaxer,Ldone
#ifdef MPI
   USE SMPI_MATH_MODULE
   USE ScaLapack_module, ONLY: Init_scala,L_useless
#endif
   use m_time_evaluate ,ONLY: memory_sum
   USE band_structure, ONLY: build_bands
   IMPLICIT NONE
   INTEGER(I4B) :: it1,it2
   LOGICAL :: l1=.true.,l2=.false.
   LOGICAL :: Lares_p
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef MPI
   !> Initialize the mpi environment {comm, size, rank}
   CALL smpi_init()
#endif

#include "my_macro.h"
#ifdef MPI
   IF (PARALLEL%ISROOT) THEN
#endif
   WRITE(6,*) '=========================== ARES =============================='
   WRITE(6,*) '>>>>>>> Ab initio Real space Electronic Structer solver <<<<<<<'
   VERSION_GIT
#ifdef MPI
   END IF
#endif

   !> Initialize the calculating information {parameters,
   !> position, pseudopotential}
   !Read file
   CALL read_file('ares.in')
#ifdef MPI
   CALL Init_scala()
   if(L_useless)then
      print*,'error,rank',parallel%myid,'is not be used for scala setting'
      call SMPI_EXIT()
   endif
#endif
   IF(trim(adjustl(outfile))/='screen')then
      open(6,file=trim(adjustl(outfile)))
   ENDIF
   SELECTCASE(ISTART)
   case(0)
#ifndef MPI
       WRITE(6,*)'[TASK] Structural relaxtion'
#else
       IF(parallel%isroot) WRITE(6,*)'[TASK] Structural relaxtion'
#endif
       IF(Nssp>0)THEN
          CALL initialize_relax()
       ENDIF
       DO Isp=1,Nssp
#ifndef MPI
          WRITE(6,*)'+++++++++++++++++++++++++++++++++++++++++++++'
          WRITE(6,*)'Present Ionic Step',Isp
          WRITE(6,*)'+++++++++++++++++++++++++++++++++++++++++++++'
#else
          IF(parallel%isroot) &
               & WRITE(6,*)'+++++++++++++++++++++++++++++++++++++++++++++'
          IF(parallel%isroot) &
               & WRITE(6,*)'Present Ionic Step',Isp
          IF(parallel%isroot) &
               & WRITE(6,*)'+++++++++++++++++++++++++++++++++++++++++++++'
#endif
          !count the SCF time
#ifndef MPI
          CALL system_clock(it1)
#else
          CALL start_time('[Total scf time]',l1)
          CALL start_time('Init_data',l2)
          if(Isp ==1 )call memory_sum("time_counts",real(100,DP)*(100+4*DP+I4B))
#endif
          !SCF,Force,Stress
          CALL cal_force_stress()
          !count the SCF time
#ifndef MPI
          CALL system_clock(it2)
          print*,'[Total SCF time]',(it2-it1)/10000.d0
#else
          CALL end_time('[Total scf time]',l1)
          ! if(parallel%isroot)CALL write_time('[Total scf time]',l1)
          ! if(parallel%isroot)call write_sum_time('comm_ata',.true.)
          ! if(parallel%isroot)call write_sum_time('comm_assign1',.true.)
          ! if(parallel%isroot)call write_sum_time('comm_assign2',.true.)
          ! if(parallel%isroot)call write_sum_time('comm_assign3',.true.)
          ! if(parallel%isroot)call write_sum_time('nabla1',.true.)
          ! if(parallel%isroot)call write_sum_time('nabla2',.true.)
#endif
          !============================================
          IF(Ldone.OR.Isp>=Nssp)  EXIT
          !============================================
          !move the ion and cell accroding force stress
          CALL relaxer(Isp)
       ENDDO
   CASE default
#ifndef MPI
        PRINT*,'[TASK] Construct the band structures'
        call build_bands()
#else
        IF(parallel%isroot)THEN
           print*,'[TASK] Construct the band structures'
        ENDIF
        call build_bands()
#endif
   ENDSELECT
   !> destroy grid
#ifndef MPI
   CALL destroy_beast()
   WRITE(6,*) '=================================Well Done==================================='
#else
   ! CALL MPI_Barrier(parallel%comm,mpinfo)
   IF(parallel%isroot)THEN
      WRITE(6,*) '=================================Well Done==================================='
      IF(trim(adjustl(outfile))/='screen')then
         close(6)
      ENDIF
   ENDIF
   CALL MPI_FINALIZE(mpinfo)
   !stop 'ares done'
#endif
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDPROGRAM ARES
