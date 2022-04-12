MODULE end_module
   USE struct_module, ONLY: destroy_struct
   USE grid_module, ONLY: destroy_rgrid
   USE mixer_module , ONLY : destroy_mixer
   USE Fourier,         only : CleanFFT
   USE smearing_module , ONLY :destroy_smear
CONTAINS
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   SUBROUTINE destroy_beast()
      IMPLICIT NONE
      !
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !print*,'111'
      CALL destroy_struct()
      !print*,'112'
      CALL destroy_rgrid()
      !print*,'113'
      !CALL CleanFFT()
      !print*,'114'
      CALL destroy_mixer()
      !print*,'115'
      CALL destroy_smear()
      !print*,'116'
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE destroy_beast
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE end_module
