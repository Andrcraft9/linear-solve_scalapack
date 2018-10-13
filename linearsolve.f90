! Solve linear system A X = B
! Use scalapack
!
program linearsolve
    use mpi
    implicit none

    integer, parameter :: DLEN_ = 9, DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,  &
                          MB_ = 5, NB_ = 6, RSRC_ = 7, CSRC_ = 8, LLD_ = 9 
    integer, parameter :: M = 100, N = 100
    integer, parameter :: MB = 2, NB = 2
    integer, parameter :: MXLLDA = 50, MXLOCC = 50, MXLOCR = 50
    integer, parameter :: RSRC = 0, CSRC = 0
    
    double precision :: A( MXLLDA, MXLOCC ), A0( MXLLDA, MXLOCC )
    double precision :: B( MXLLDA, 1 ), B0( MXLLDA, 1 )
    integer :: DESCA( DLEN_ ), DESCB( DLEN_ )
    integer :: IPIV( MXLOCR + NB ), WORK( MXLOCR )
    
    integer :: ICTXT, INFO, IASEED
    integer :: MYCOL, MYROW, NPCOL, NPROW
    integer :: IAM, NPROCS, k

    double precision :: ANORM, BNORM, EPS, RESID, XNORM
    ! Functions
    double precision :: PDLAMCH, PDLANGE
    real*8 :: stime, etime

    integer :: ierr
    integer, dimension(2) :: p_size

    p_size = (/0,0/)
    ierr = 0

    ! Initialize an NPROW x NPCOL process grid using a row-major
    ! ordering  of  the  processes. This routine retrieves a default system
    ! context  which  will  include all available processes. In addition it
    ! spawns the processes if needed.
    call BLACS_PINFO( IAM, NPROCS )
    call MPI_DIMS_CREATE(NPROCS, 2, p_size, ierr)
    NPROW = p_size(1)
    NPCOL = p_size(2)
    if (IAM == 0) then
        print *, "NPROW = ", NPROW, " NPCOL = ", NPCOL
    endif

    ! If machine needs additional set up, do it now
    if ( NPROCS.LT.1 ) then
        if ( IAM.EQ.0 ) NPROCS = NPROW*NPCOL
        call BLACS_SETUP( IAM, NPROCS )
    endif
    ! Define process grid
    call BLACS_GET( -1, 0, ICTXT )
    call BLACS_GRIDINIT( ICTXT, 'Row-major', NPROW, NPCOL )
    call BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

    print *, IAM, MYROW, MYCOL, NPROW, NPCOL

    ! Initialize the array descriptors for the matrices A
    call DESCINIT( DESCA, M, N, MB, NB, RSRC, CSRC, ICTXT, MXLLDA, INFO )
    call DESCINIT( DESCB, N, 1, NB, 1, RSRC, CSRC, ICTXT, MXLLDA, INFO )

    do k = 1, MXLLDA
        B(k, :) = 1.0d0
    enddo

    ! Generate random matrix
    !call PDMATGEN(ICTXT, 'R', 'D', M, N, MB, NB, A, LDA,  &
    !              IAROW, IACOL, ISEED, IROFF, IRNUM, ICOFF,  &
    !              ICNUM, MYROW, MYCOL, NPROW, NPCOL )
    IASEED = 10
    call PDMATGEN( ICTXT, 'N', 'N', DESCA( M_ ),         &
                   DESCA( N_ ), DESCA( MB_ ),            &
                   DESCA( NB_ ), A,                      &
                   DESCA( LLD_ ), DESCA( RSRC_ ),        &
                   DESCA( CSRC_ ), IASEED, 0, MXLOCR, 0, MXLOCC, &
                   MYROW, MYCOL, NPROW, NPCOL )


    ! Make a copy of A and B for checking purposes
    call PDLACPY( 'All', N, N, A, 1, 1, DESCA, A0, 1, 1, DESCA )
    call PDLACPY( 'All', N, 1, B, 1, 1, DESCB, B0, 1, 1, DESCB )

    ! Linear solve
    stime = mpi_wtime()
    call PDGESV( N, 1, A, 1, 1, DESCA, IPIV, B, 1, 1, DESCB, INFO )
    etime = mpi_wtime()

    ! Resudial computational
    !EPS = PDLAMCH( ICTXT, 'Epsilon' )
    !ANORM = PDLANGE( 'I', N, N, A, 1, 1, DESCA, WORK )
    !BNORM = PDLANGE( 'I', N, 1, B, 1, 1, DESCB, WORK )
    call PDGEMM( 'N', 'N', N, 1, N, 1.0d0, A0, 1, 1, DESCA, B, 1, 1,   &
                 DESCB, -1.0d0, B0, 1, 1, DESCB )
    XNORM = PDLANGE( 'I', N, 1, B0, 1, 1, DESCB, WORK )
    !RESID = XNORM / ( ANORM*BNORM*EPS*DBLE( N ) )
    RESID = XNORM

    if (IAM == 0) then
        print *, "resid = ", RESID, " time for linear solver = ", etime-stime
    endif

    call BLACS_GRIDEXIT( ICTXT )

    ! Exit the BLACS
    call BLACS_EXIT( 0 )

end program
