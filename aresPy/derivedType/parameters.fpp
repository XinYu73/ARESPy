# 1 "parameters.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "parameters.f90"
module parameters

    implicit none
    private
    public :: idp, isp
    
    INTEGER, PARAMETER :: idp=kind(1d0) ! double precision, np.float64 equivalent
    INTEGER, PARAMETER :: isp=kind(1e0) ! single precision, np.float32 equivalent
    
    contains
    
    end module parameters
