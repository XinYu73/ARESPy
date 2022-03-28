subroutine matrixMult(n1, A, D)

    INTEGER, intent(in) :: n1
    DOUBLE PRECISION, intent(in) :: A(n1,n1)
    DOUBLE PRECISION, intent(out) :: D(n1)
    CALL EVLRG(A,D)
    RETURN

end subroutine matrixMult