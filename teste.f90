PROGRAM TESTE
    IMPLICIT REAL*8(a-h,o-z)

    S = rand()*(20*5)
    T = rand()*(20*5)
    U = rand()*(20*5)

    I = IDNINT(S)
    J = IDNINT(T)
    K = IDNINT(U)

    PRINT*, S, T, U
    PRINT*, I, J, K

END PROGRAM