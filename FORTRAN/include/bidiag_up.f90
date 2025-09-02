!>
!> bidiag_up
!>
!> Algorithm for updating the bidiagonal factorization  
!> 
!> B = Q*(B0 + w*p')*P'
!>
!> Here w,p are vectors of size m,n respectively. The algorithm uses two phases.
!> In phase 1 the bidiagonal is transformed into a quaddiagonal (i.e., four diagonal bands).
!> In phase 2 the quaddiagonal is restored to bidiagonal.
!>
!> B1 = Q2*Q1*(B0 + w*p')*P1'               (phase 1)
!> B  = Q3*B1*P2'                           (phase 2)
!>
!> The algorithms do not form orthogonal matrices Q1,Q2,Q3 and P1,P2 explicitly. Instead 
!> the orthogonal transformations are implicitly represented by Givens rotations.
!>
!> J.J. Brust
!> Initial version Aug. 2025
!> johannesbrust@yahoo.com
!> 
!>----------------------------------------------------------------------------------------
!> 01/07/25, J.B., Initial version 
!> 01/09/25, J.B., Buldge chasing, implementation of phase 1,
!> 01/24/25, J.B., Addition of a C interface
!> 01/25/25, J.B., C interface for phase 2
!> 02/28/25, J.B., Givens apply, column-wise storage of rotations in the c-interface
!> 03/03/25, J.B., C interface for givens apply
!> 08/28/25, J.B., Machine precision tolerance in givel, used by 
!>                 bdup_mulp and bdup_mulq
                
module bidiag_up

    use iso_c_binding
    use kind_parameter, only: dp
    use givens

    implicit none

    contains

    !>
    !> phase1
    !> A function description in Matlab notation:
    !
    ! function [Q,QB,B,P] = givphase1_SPMF(B,w,p)
    !      GIVPHASE1_SPMF Sparse matrix free Phase 1 of a Given's factorization of 
    !         bidiag plus rank-1. Computing the factors of B4 = Q'*(B+w*p')*P given 
    !         that B is bidiagonal and B4 is quaddiagonal. This algorithm uses 
    !         permutations via row/column swaps and a sequence 
    !         of Given's rotations so that e.g., p'*P1 = [0,*,*,...,*]. Similarly 
    !         Q1'*w = [0,*;*;...;*]. When fill-ins occur this method restores 
    !         the quaddiagonal form. This "matrix-free" version represents the orthogonal 
    !         matrices implictly.
    !      
    !         Input:
    !         B = Bidiagonal
    !         w = sparse rank one vector (mx1)
    !         p = sparse rank one vector (nx1)
    !      
    !         Output:
    !         Q   = Givens information of Orthogonal matrix ((n*n)/2x4)
    !         QB  = Givens information of Orthogonal matrix ((m-n)x4)
    !         P   = Givens information of Orthogonal matrix ((n*n)/2x4)
    !         B   = Quaddiagonal (mxn)
    ! --------------------------------------------------------------------------
    subroutine phase1(B,w,p,Q,QB,PB)

        !>
        !> Initializations
        !>
        real(dp), intent(out)                   :: B(:,:)               ! bidiagonal
        real(dp), intent(out)                   :: w(:)
        real(dp), intent(out)                   :: p(:)
        real(dp), allocatable, intent(out)      :: Q(:,:)
        real(dp), allocatable, intent(out)      :: QB(:,:)
        real(dp), allocatable, intent(out)      :: PB(:,:)

        integer                                 :: qbnd = 2
        integer                                 :: pbnd = 1
        integer                                 :: m
        integer                                 :: n
        integer                                 :: nit
        integer                                 :: nitp
        integer                                 :: nitq
        integer                                 :: nitqb
        integer                                 :: i
        integer                                 :: i1
        integer                                 :: i2
        integer                                 :: b1
        integer                                 :: b2
        integer                                 :: ck
        integer                                 :: ii

        real(dp)                                :: btol = 1.0e-13_dp
        real(dp), allocatable                   :: buff(:)
        real(dp)                                :: G(2,2)

        m = size(w)
        n = size(p)

        allocate( Q(ceiling( n*n/2.0 )+1, 4))
        allocate(QB(( m-n )+1, 4))
        allocate(PB(ceiling( n*n/2.0 )+1, 4))        
        allocate(buff(2*qbnd+pbnd+1))

        nit     = n-1
        nitp    = 1
        nitq    = 1
        nitqb   = 1 

        !>
        !> Main loop
        !>
        do i = 1, nit

            !> Zeroing an element in p
            if ( abs(p(i)) > btol ) then 

                i1 = max(i - qbnd, 1)
                i2 = min(i1 + qbnd + pbnd, m)
            
                ! Check if a permutation is sufficient
                if ( abs(p(i+1)) < btol ) then

                    buff(1)     = p(i)
                    p(i)        = p(i+1)
                    p(i+1)      = buff(1)
                    !p(i+1)      = p(i)
                    !p(i)        = 0_dp

                    !p(i:(i+1)) = [0_dp,buff(1)]
                    !p(i:(i+1)) = [0,p(i)]
        
                    ! Swap columns or their indices 
                    buff(1:((i2-i1)+1)) = B(i1:i2,i)                           ! 1:length(i1:i2)
                    B(i1:i2,i)          = B(i1:i2,(i+1))
                    B(i1:i2,i+1)        = buff(1:((i2-i1)+1))                  ! 1:length(i1:i2)
                    
                    nitp            = nitp + 1
                    PB(nitp,1:2)    = [0_dp,1_dp]
                    PB(nitp,3:4)    = [i,i+1]  

                    !PB(nitp,:)  = [0,1,i,i+1]

                else

                    call giv(G,p(i:(i+1)),2,1)                    
                    
                    ![G,y] = giv(p(i:(i+1))');
                    !p(i:(i+1)) = y(:);
                    
                    !B(i1:i2,i:(i+1)) = B(i1:i2,i:(i+1))*G;
                    
                    B(i1:i2,i:(i+1)) = matmul(B(i1:i2,i:(i+1)),G)

                    ! Orthogonal matrix
                    !PP(:,pix(i:(i+1))) = PP(:,pix(i:(i+1)))*G;

                    nitp            = nitp + 1;
                    PB(nitp,1:2)    = [G(1,1),G(1,2)]
                    PB(nitp,3:4)    = [i,i+1]
                    !PB(nitp,:)  = [G(1,1),G(1,2),i,i+1]

                end if
                
                b1 = i - qbnd
                b2 = i + 1
            
                ! Chasing a buldge (after column reduction)
                ck = 1

                do while ( (b1 > 0) .and. (b2 > 0) )

                    ! Break-out early if buldge value is small
                    if ( abs(B(b1,b2)) < btol ) then 
                        exit
                    end if

                    ! remove the buldge using row or column transformations
                    if ( mod(ck,2) == 0 ) then 

                        ii = b2

                        call rmbp1(B,PB,G,buff,nitp,b1,b2,qbnd,pbnd,m,btol,1)

                        ! update buldge indices
                        b1 = ii - qbnd
                        b2 = ii + 1
                        ck = ck + 1

                    else

                        ii = b1

                        call rmbp1(B,Q,G,buff,nitq,b1,b2,qbnd,pbnd,n,btol,2)

                        ! update buldge indices
                        b1 = ii + 1
                        b2 = ii - pbnd
                        ck = ck + 1

                    end if

                end do

            end if

            ! Skipping the last elimination of w
            if ( i == nit ) then 
                cycle   ! exit
            end if

            ! Zeroing an element in w 
            if ( abs(w(i)) > btol ) then 

                i1 = max(i - pbnd, 1)
                i2 = min(i1 + qbnd + pbnd, n)

                ! Check if a permutation is sufficient
                if ( abs(w(i+1)) < btol ) then 

                    buff(1)     = w(i)
                    w(i)        = w(i+1)
                    w(i+1)      = buff(1)
                    !w(i+1)      = w(i)
                    !w(i)        = 0_dp

                    ! Swap rows or their indices 
                    buff(1:((i2-i1)+1)) = B(i,i1:i2)                           
                    B(i,i1:i2)          = B((i+1),i1:i2)
                    B(i+1,i1:i2)        = buff(1:((i2-i1)+1))                  
                    
                    nitq            = nitq + 1
                    Q(nitq,1:2)     = [0_dp,1_dp]
                    Q(nitq,3:4)     = [i,i+1]

                else ! Givens rotation

                    call giv(G,w(i:(i+1)),1,1)

                    B(i:(i+1),i1:i2)    = matmul(G,B(i:(i+1),i1:i2))

                    nitq                = nitq + 1;
                    Q(nitq,1:2)         = [G(1,1),G(1,2)]
                    Q(nitq,3:4)         = [i,i+1]

                end if

                b1 = i + 1
                b2 = i - pbnd

                ! Chasing a buldge (after row reduction)
                ck = 2

                do while ( (b1 > 0) .and. (b2 > 0) )

                    ! Break-out early if buldge value is small
                    if ( abs(B(b1,b2)) < btol ) then 
                        exit
                    end if

                    ! remove the buldge using row or column transformations
                    if ( mod(ck,2) == 0 ) then 

                        ii = b2

                        call rmbp1(B,PB,G,buff,nitp,b1,b2,qbnd,pbnd,m,btol,1)

                        ! update buldge indices
                        b1 = ii - qbnd
                        b2 = ii + 1
                        ck = ck + 1

                    else

                        ii = b1

                        call rmbp1(B,Q,G,buff,nitq,b1,b2,qbnd,pbnd,n,btol,2)

                        ! update buldge indices
                        b1 = ii + 1
                        b2 = ii - pbnd
                        ck = ck + 1

                    end if

                end do

            end if

        end do

        !> Final updates (remove a spike if needed)

        B((n-1):m,n) = B((n-1):m,n) + p(n)*w(n-1:m)

        do i = m, (n+2), -1

            ii = i - 1 

            ! Check for zero buldge
            if ( abs( B(ii+1,n) ) > btol ) then 

                ! Check if a permutation is possible
                if ( abs(B(ii,n)) < btol ) then 

                    buff(1)         = B(ii+1,n)
                    B(ii+1,n)       = B(ii,n)
                    B(ii,n)         = buff(1)
                    nitqb           = nitqb + 1

                    QB(nitqb,1:2)   = [0_dp,1_dp]
                    QB(nitqb,3:4)   = [ii,ii+1]

                else 

                    call giv(G,B(ii:(ii+1),n),1,2)

                    nitqb           = nitqb + 1  
                    QB(nitqb,1:2)   = [G(1,1),G(1,2)]
                    QB(nitqb,3:4)   = [ii,ii+1]

                end if

            end if 

        end do

        Q(1,1)  = nitq - 1
        PB(1,1) = nitp - 1
        QB(1,1) = nitqb - 1 

    end subroutine

    !>
    !> C interface
    !> 
    subroutine phase1_c(B,w,p,Q,QB,PB,m,n,m1,m2) bind(C,name="phase1_c")

        !>
        !> Initializations
        !>
        real(c_double), intent(out)                     :: B(m,n)       ! bidiagonal
        real(c_double), intent(out)                     :: w(m)
        real(c_double), intent(out)                     :: p(n)
        real(c_double), intent(out)                     :: Q(4,m1)      ! Q(m1,4)
        real(c_double), intent(out)                     :: QB(4,m2)     ! QB(m2,4)
        real(c_double), intent(out)                     :: PB(4,m1)     ! PB(m1,4)   

        integer                                         :: qbnd = 2
        integer                                         :: pbnd = 1
        integer(c_int), intent(in)                      :: m
        integer(c_int), intent(in)                      :: n
        integer(c_int), intent(in)                      :: m1
        integer(c_int), intent(in)                      :: m2
        integer                                         :: nit
        integer                                         :: nitp
        integer                                         :: nitq
        integer                                         :: nitqb
        integer                                         :: i
        integer                                         :: i1
        integer                                         :: i2
        integer                                         :: b1
        integer                                         :: b2
        integer                                         :: ck
        integer                                         :: ii

        real(dp)                                        :: btol = 1.0e-13_dp
        real(dp), allocatable                           :: buff(:)
        real(dp)                                        :: G(2,2)

        !m = size(w)
        !n = size(p)
               
        allocate(buff(2*qbnd+pbnd+1))

        nit     = n-1
        nitp    = 1
        nitq    = 1
        nitqb   = 1 

        !>
        !> Main loop
        !>

        do i = 1, nit

            !> Zeroing an element in p
            if ( abs(p(i)) > btol ) then 

                i1 = max(i - qbnd, 1)
                i2 = min(i1 + qbnd + pbnd, m)
            
                ! Check if a permutation is sufficient
                if ( abs(p(i+1)) < btol ) then

                    buff(1)     = p(i)
                    p(i)        = p(i+1)
                    p(i+1)      = buff(1)
                    !p(i+1)      = p(i)
                    !p(i)        = 0_dp

                    !p(i:(i+1)) = [0_dp,buff(1)]
                    !p(i:(i+1)) = [0,p(i)]
        
                    ! Swap columns or their indices 
                    buff(1:((i2-i1)+1)) = B(i1:i2,i)                           ! 1:length(i1:i2)
                    B(i1:i2,i)          = B(i1:i2,(i+1))
                    B(i1:i2,i+1)        = buff(1:((i2-i1)+1))                  ! 1:length(i1:i2)
                    
                    nitp            = nitp + 1
                                        
                    PB(1:2,nitp)    = [0_dp,1_dp]
                    PB(3:4,nitp)    = [i,i+1]

                    !PB(nitp,1:2)    = [0_dp,1_dp]
                    !PB(nitp,3:4)    = [i,i+1]

                    !PB(nitp,:)  = [0,1,i,i+1]

                else

                    call giv_c(G,p(i:(i+1)),2,1)                    
                    
                    ![G,y] = giv(p(i:(i+1))');
                    !p(i:(i+1)) = y(:);
                    
                    !B(i1:i2,i:(i+1)) = B(i1:i2,i:(i+1))*G;
                    
                    B(i1:i2,i:(i+1)) = matmul(B(i1:i2,i:(i+1)),G)

                    ! Orthogonal matrix
                    !PP(:,pix(i:(i+1))) = PP(:,pix(i:(i+1)))*G;

                    nitp            = nitp + 1;                    

                    PB(1:2,nitp)    = [G(1,1),G(1,2)]
                    PB(3:4,nitp)    = [i,i+1]

                    ! PB(nitp,1:2)    = [G(1,1),G(1,2)]
                    ! PB(nitp,3:4)    = [i,i+1]

                    !PB(nitp,:)  = [G(1,1),G(1,2),i,i+1]

                end if
                
                b1 = i - qbnd
                b2 = i + 1
            
                ! Chasing a buldge (after column reduction)
                ck = 1

                do while ( (b1 > 0) .and. (b2 > 0) )

                    ! Break-out early if buldge value is small
                    if ( abs(B(b1,b2)) < btol ) then 
                        exit
                    end if

                    ! remove the buldge using row or column transformations
                    if ( mod(ck,2) == 0 ) then 

                        ii = b2

                        call rmbp1_c(B,PB,G,buff,nitp,b1,b2,qbnd,pbnd,m,btol,1)
                        !call rmbp1(B,PB,G,buff,nitp,b1,b2,qbnd,pbnd,m,btol,1)

                        ! update buldge indices
                        b1 = ii - qbnd
                        b2 = ii + 1
                        ck = ck + 1

                    else

                        ii = b1

                        call rmbp1_c(B,Q,G,buff,nitq,b1,b2,qbnd,pbnd,n,btol,2)
                        !call rmbp1(B,Q,G,buff,nitq,b1,b2,qbnd,pbnd,n,btol,2)

                        ! update buldge indices
                        b1 = ii + 1
                        b2 = ii - pbnd
                        ck = ck + 1

                    end if

                end do

            end if

            ! Skipping the last elimination of w
            if ( i == nit ) then 
                cycle   ! exit
            end if

            ! Zeroing an element in w 
            if ( abs(w(i)) > btol ) then 

                i1 = max(i - pbnd, 1)
                i2 = min(i1 + qbnd + pbnd, n)

                ! Check if a permutation is sufficient
                if ( abs(w(i+1)) < btol ) then 

                    buff(1)     = w(i)
                    w(i)        = w(i+1)
                    w(i+1)      = buff(1)
                    !w(i+1)      = w(i)
                    !w(i)        = 0_dp

                    ! Swap rows or their indices 
                    buff(1:((i2-i1)+1)) = B(i,i1:i2)                           
                    B(i,i1:i2)          = B((i+1),i1:i2)
                    B(i+1,i1:i2)        = buff(1:((i2-i1)+1))                  
                    
                    nitq            = nitq + 1
                    !Q(nitq,1:2)     = [0_dp,1_dp]
                    !Q(nitq,3:4)     = [i,i+1]
                    Q(1:2,nitq)     = [0_dp,1_dp]
                    Q(3:4,nitq)     = [i,i+1]

                else ! Givens rotation

                    call giv_c(G,w(i:(i+1)),1,1)

                    B(i:(i+1),i1:i2)    = matmul(G,B(i:(i+1),i1:i2))

                    nitq                = nitq + 1;
                    ! Q(nitq,1:2)         = [G(1,1),G(1,2)]
                    ! Q(nitq,3:4)         = [i,i+1]
                    Q(1:2,nitq)         = [G(1,1),G(1,2)]
                    Q(3:4,nitq)         = [i,i+1]

                end if

                b1 = i + 1
                b2 = i - pbnd

                ! Chasing a buldge (after row reduction)
                ck = 2

                do while ( (b1 > 0) .and. (b2 > 0) )

                    ! Break-out early if buldge value is small
                    if ( abs(B(b1,b2)) < btol ) then 
                        exit
                    end if

                    ! remove the buldge using row or column transformations
                    if ( mod(ck,2) == 0 ) then 

                        ii = b2

                        call rmbp1_c(B,PB,G,buff,nitp,b1,b2,qbnd,pbnd,m,btol,1)
                        !call rmbp1(B,PB,G,buff,nitp,b1,b2,qbnd,pbnd,m,btol,1)

                        ! update buldge indices
                        b1 = ii - qbnd
                        b2 = ii + 1
                        ck = ck + 1

                    else

                        ii = b1

                        call rmbp1_c(B,Q,G,buff,nitq,b1,b2,qbnd,pbnd,n,btol,2)
                        !call rmbp1(B,Q,G,buff,nitq,b1,b2,qbnd,pbnd,n,btol,2)

                        ! update buldge indices
                        b1 = ii + 1
                        b2 = ii - pbnd
                        ck = ck + 1

                    end if

                end do

            end if

        end do

        !> Final updates (remove a spike if needed)

        B((n-1):m,n) = B((n-1):m,n) + p(n)*w(n-1:m)

        do i = m, (n+2), -1

            ii = i - 1 

            ! Check for zero buldge
            if ( abs( B(ii+1,n) ) > btol ) then 

                ! Check if a permutation is possible
                if ( abs(B(ii,n)) < btol ) then 

                    buff(1)         = B(ii+1,n)
                    B(ii+1,n)       = B(ii,n)
                    B(ii,n)         = buff(1)
                    nitqb           = nitqb + 1

                    !QB(nitqb,1:2)   = [0_dp,1_dp]
                    !QB(nitqb,3:4)   = [ii,ii+1]

                    QB(1:2,nitqb)   = [0_dp,1_dp]
                    QB(3:4,nitqb)   = [ii,ii+1]

                else 

                    call giv_c(G,B(ii:(ii+1),n),1,2)

                    nitqb           = nitqb + 1  
                    ! QB(nitqb,1:2)   = [G(1,1),G(1,2)]
                    ! QB(nitqb,3:4)   = [ii,ii+1]

                    QB(1:2,nitqb)   = [G(1,1),G(1,2)]
                    QB(3:4,nitqb)   = [ii,ii+1]

                end if

            end if 

        end do

        Q(1,1)  = nitq - 1
        PB(1,1) = nitp - 1
        QB(1,1) = nitqb - 1 

    end subroutine

    !> second phase of the algorithm
    !>
    !> phase2
    !>
    !> Function description in Matlab notation:
    !
    !   function [Q,B,P] = givphase2_SPMF(B)
    !     GIVPHASE2_SPMF Phase 2 of a Given's factorization of bidiagonal plus rank-1
    !     Computing the factors of B4 = Q'*Q4'*B4*P4*P given that B is banded
    !     with bandwidth 4, (pb=1,qb=2). This algorithm uses a sequence of Given's  
    !     rotations so that the banded matrix is reduced to bidiagonal. 
    !     The offending buldge elements are chased off the matrix. This version
    !     does not explicitly compute the orthogonal matrices but only stores
    !     the Givens coefficients and indicies of transformed rows/columns.
    !  
    !     Input:
    !     B = Banded (4 diagonals)
    !  
    !     Output:
    !     B   = Bidiagonal
    !     Q   = Givens information of Orthogonal matrix ((n*n)/2x4)
    !     P   = Givens information of Orthogonal matrix ((n*n)/2x4)
    ! ---------------------------------------------------------------------------
    subroutine phase2(B,Q,P)

        !>
        !> Initializations
        !>
        real(dp), intent(out)                   :: B(:,:)               ! bidiagonal
        real(dp), allocatable, intent(out)      :: Q(:,:)
        real(dp), allocatable, intent(out)      :: P(:,:)

        integer                                 :: qbnd = 2
        integer                                 :: pbnd = 1
        integer                                 :: m
        integer                                 :: n
        integer                                 :: nit
        integer                                 :: nitp
        integer                                 :: nitq        
        integer                                 :: i
        integer                                 :: b1
        integer                                 :: b2
        integer                                 :: ck
       
        real(dp)                                :: btol = 1.0e-13_dp
        real(dp), allocatable                   :: buff(:)
        real(dp)                                :: G(2,2)

        m = size(B,dim=1)
        n = size(B,dim=2)

        allocate( Q(ceiling( n*n/2.0 )+1, 4))
        allocate( P(ceiling( n*n/2.0 )+1, 4))
        allocate(buff(2*qbnd+pbnd+1))

        if ( m == n ) then 
            nit = n
        else
            nit = n + 1
        end if
        
        nitp    = 1
        nitq    = 1        

        !> main loop 
        do i = 1, nit

            ! subdiagonal buldge
            b1 = i + pbnd
            b2 = i

            if ( b1 <= min(n+1,m) ) then 

                if ( abs(B(b1,b2)) > btol ) then 

                    ! debugging
                    !print *, B

                    ! buff(1:2) = B((b1-1):b1,b2)
                    ! call giv(G,buff(1:2),1,2)
                    ! print *, buff(1:2)

                    call rmbp2(B,Q,G,buff,nitq,b1,b2,qbnd,pbnd,n,btol,2)

                    !print *, B

                    b1 = b1 - pbnd
                    b2 = b2 + pbnd + qbnd

                    ! Chasing a buldge (after row reduction)
                    ck = 1
                    do while ( ( b1 <= min(n+1,m) ) .and. ( b2 <= n ) )

                        ! Break out early for a small buldge
                        if ( abs(B(b1,b2)) < btol ) then 
                            exit
                        end if

                        ! Elimitate buldge
                        if ( mod(ck,2) .ne. 0 ) then
                            
                            call rmbp2(B,P,G,buff,nitp,b1,b2,qbnd,pbnd,m,btol,1)

                            b1 = b1 + pbnd + qbnd + 1
                            b2 = b2 - pbnd
                            ck = ck + 1

                        else

                            call rmbp2(B,Q,G,buff,nitq,b1,b2,qbnd,pbnd,n,btol,2)

                            b1 = b1 - pbnd
                            b2 = b2 + pbnd + qbnd + 1
                            ck = ck + 1

                        end if

                    end do

                end if

            end if

            ! Debugging
            ! print *, 'After first subdiagonal'
            ! print *, B

            ! Superdiagonal 
            b1 = i
            b2 = i + qbnd

            if ( b2 <= n ) then 

                if ( abs(B(b1,b2)) > btol ) then 

                    call rmbp2(B,P,G,buff,nitp,b1,b2,qbnd,pbnd,m,btol,1)

                    b1 = b1 + pbnd + qbnd
                    b2 = b2 - pbnd

                    ! Chasing buldge (after column reduction)
                    ck = 2

                    ! Debugging
                    ! print *, 'Before 1st chase after column reduction'
                    ! print *, B
                    ! print *, b1
                    ! print *, b2 
                    ! print *, n
                    ! print *, m

                    do while ( ( b1 <= min(n+1,m) ) .and. ( b2 <= n ) )

                        

                        ! Break out early for a small buldge
                        if ( abs(B(b1,b2)) < btol ) then
                            !print *, 'Early exit' 
                            exit
                        end if

                        ! Debugging
                        ! print *, b1
                        ! print *, b2 
                        ! print *, B

                        ! Elimitate buldge
                        if ( mod(ck,2) .ne. 0 ) then
                            
                            call rmbp2(B,P,G,buff,nitp,b1,b2,qbnd,pbnd,m,btol,1)

                            b1 = b1 + pbnd + qbnd + 1
                            b2 = b2 - pbnd
                            ck = ck + 1

                        else

                            call rmbp2(B,Q,G,buff,nitq,b1,b2,qbnd,pbnd,n,btol,2)

                            b1 = b1 - pbnd
                            b2 = b2 + pbnd + qbnd + 1
                            ck = ck + 1

                        end if

                    end do

                end if

            end if

        end do

        ! Store transformation count
        P(1,1) = nitp - 1
        Q(1,1) = nitq - 1

    end subroutine

    !>
    !> C interface to phase 2
    ! ---------------------------------------------------------------------------
    subroutine phase2_c(B,Q,P,m,n,m1) bind(C,name="phase2_c")

        !>
        !> Initializations
        !>
        real(c_double), intent(out)             :: B(m,n)               ! bidiagonal
        real(c_double), intent(out)             :: Q(4,m1)              ! Q(m1,4)
        real(c_double), intent(out)             :: P(4,m1)              ! P(m1,4)

        integer                                 :: qbnd = 2
        integer                                 :: pbnd = 1
        integer(c_int), intent(in)              :: m
        integer(c_int), intent(in)              :: n
        integer(c_int), intent(in)              :: m1
        integer                                 :: nit
        integer                                 :: nitp
        integer                                 :: nitq        
        integer                                 :: i
        integer                                 :: b1
        integer                                 :: b2
        integer                                 :: ck
       
        real(dp)                                :: btol = 1.0e-13_dp
        real(dp), allocatable                   :: buff(:)
        real(dp)                                :: G(2,2)

        !m = size(B,dim=1)
        !n = size(B,dim=2)

        !allocate( Q(ceiling( n*n/2.0 )+1, 4))
        !allocate( P(ceiling( n*n/2.0 )+1, 4))
        allocate(buff(2*qbnd+pbnd+1))

        if ( m == n ) then 
            nit = n
        else
            nit = n + 1
        end if
        
        nitp    = 1
        nitq    = 1        

        !> main loop 
        do i = 1, nit

            ! subdiagonal buldge
            b1 = i + pbnd
            b2 = i

            if ( b1 <= min(n+1,m) ) then 

                if ( abs(B(b1,b2)) > btol ) then 

                    ! debugging
                    !print *, B

                    ! buff(1:2) = B((b1-1):b1,b2)
                    ! call giv(G,buff(1:2),1,2)
                    ! print *, buff(1:2)

                    call rmbp2_c(B,Q,G,buff,nitq,b1,b2,qbnd,pbnd,n,btol,2)

                    !print *, B

                    b1 = b1 - pbnd
                    b2 = b2 + pbnd + qbnd

                    ! Chasing a buldge (after row reduction)
                    ck = 1
                    do while ( ( b1 <= min(n+1,m) ) .and. ( b2 <= n ) )

                        ! Break out early for a small buldge
                        if ( abs(B(b1,b2)) < btol ) then 
                            exit
                        end if

                        ! Elimitate buldge
                        if ( mod(ck,2) .ne. 0 ) then
                            
                            call rmbp2_c(B,P,G,buff,nitp,b1,b2,qbnd,pbnd,m,btol,1)

                            b1 = b1 + pbnd + qbnd + 1
                            b2 = b2 - pbnd
                            ck = ck + 1

                        else

                            call rmbp2_c(B,Q,G,buff,nitq,b1,b2,qbnd,pbnd,n,btol,2)

                            b1 = b1 - pbnd
                            b2 = b2 + pbnd + qbnd + 1
                            ck = ck + 1

                        end if

                    end do

                end if

            end if

            ! Debugging
            ! print *, 'After first subdiagonal'
            ! print *, B

            ! Superdiagonal 
            b1 = i
            b2 = i + qbnd

            if ( b2 <= n ) then 

                if ( abs(B(b1,b2)) > btol ) then 

                    call rmbp2_c(B,P,G,buff,nitp,b1,b2,qbnd,pbnd,m,btol,1)

                    b1 = b1 + pbnd + qbnd
                    b2 = b2 - pbnd

                    ! Chasing buldge (after column reduction)
                    ck = 2

                    ! Debugging
                    ! print *, 'Before 1st chase after column reduction'
                    ! print *, B
                    ! print *, b1
                    ! print *, b2 
                    ! print *, n
                    ! print *, m

                    do while ( ( b1 <= min(n+1,m) ) .and. ( b2 <= n ) )

                        

                        ! Break out early for a small buldge
                        if ( abs(B(b1,b2)) < btol ) then
                            !print *, 'Early exit' 
                            exit
                        end if

                        ! Debugging
                        ! print *, b1
                        ! print *, b2 
                        ! print *, B

                        ! Elimitate buldge
                        if ( mod(ck,2) .ne. 0 ) then
                            
                            call rmbp2_c(B,P,G,buff,nitp,b1,b2,qbnd,pbnd,m,btol,1)

                            b1 = b1 + pbnd + qbnd + 1
                            b2 = b2 - pbnd
                            ck = ck + 1

                        else

                            call rmbp2_c(B,Q,G,buff,nitq,b1,b2,qbnd,pbnd,n,btol,2)

                            b1 = b1 - pbnd
                            b2 = b2 + pbnd + qbnd + 1
                            ck = ck + 1

                        end if

                    end do

                end if

            end if

        end do

        ! Store transformation count
        P(1,1) = nitp - 1
        Q(1,1) = nitq - 1

    end subroutine


    !> remove buldge phase1 (called from within phase1)
    subroutine rmbp1(B,PQ,G,buff,nit,b1,b2,qb,pb,mn,btol,mode)

        real(dp), intent(out)           :: B(:,:)
        real(dp), intent(out)           :: PQ(:,:)
        real(dp), intent(out)           :: G(:,:)
        real(dp), intent(out)           :: buff(:)
        real(dp), intent(in)            :: btol
        
        integer, intent(in)             :: b1
        integer, intent(in)             :: b2
        integer, intent(out)            :: nit        
        integer, intent(in)             :: qb
        integer, intent(in)             :: pb
        integer, intent(in)             :: mn
        integer, intent(in)             :: mode     ! mode=1 buldge row, mode=2 buldge col  

        integer                         :: ii
        integer                         :: ii1 
        integer                         :: ii2
        
        ! buldge row index fixed
        if ( mode==1 ) then
            
            ii  = b2
            ii1 = max(ii - qb,1)
            ii2 = min(ii1 + qb + pb + 1, mn)

            !> Check for a permutation
            if ( abs(B(b1,ii+1)) < btol ) then 

                ! Swap columns or their indices 
                buff(1:((ii2-ii1)+1))   = B(ii1:ii2,ii);
                B(ii1:ii2,ii)           = B(ii1:ii2,(ii+1));
                B(ii1:ii2,ii+1)         = buff(1:((ii2-ii1)+1));
                
                nit                     = nit + 1
                PQ(nit,1:2)             = [0_dp,1_dp]
                PQ(nit,3:4)             = [ii,ii+1]

            else ! zero elements using Givens

                call giv(G,B(b1,ii:(ii+1)),2,1)
                ![G,y] = giv(B(b1,ii:(ii+1)));            

                buff(1:2)               = B(b1,ii:(ii+1))
                
                B(ii1:ii2,ii:(ii+1))    = matmul(B(ii1:ii2,ii:(ii+1)),G);                
                
                B(b1,ii:(ii+1))         = buff(1:2)

                nit                     = nit + 1;
                PQ(nit,1:2)             = [G(1,1),G(1,2)]
                PQ(nit,3:4)             = [ii,ii+1]
                
                !P(nit,:) = [G(1,1),G(1,2),ii,ii+1];

            end if

        else ! buldge column index fixed

            ii  = b1
            ii1 = max(ii - pb,1)
            ii2 = min(ii1 + qb + pb + 1, mn)
            
            if ( abs(B(ii+1,b2)) < btol ) then 

                ! Swap rows or their indices 
                buff(1:((ii2-ii1)+1))   = B(ii,ii1:ii2);
                B(ii,ii1:ii2)           = B((ii+1),ii1:ii2);
                B(ii+1,ii1:ii2)         = buff(1:((ii2-ii1)+1));
                
                nit                     = nit + 1
                PQ(nit,1:2)             = [0_dp,1_dp]
                PQ(nit,3:4)             = [ii,ii+1]
                
            else 

                call giv(G,B(ii:(ii+1),b2),1,1)

                buff(1:2)               = B(ii:(ii+1),b2)

                B(ii:(ii+1),ii1:ii2)    = matmul(G,B(ii:(ii+1),ii1:ii2))

                B(ii:(ii+1),b2)         = buff(1:2)

                nit                     = nit + 1;
                PQ(nit,1:2)             = [G(1,1),G(1,2)]
                PQ(nit,3:4)             = [ii,ii+1]

            end if

        end if

    end subroutine

    !> remove buldge phase1 (called from within phase1)
    !> c interface
    subroutine rmbp1_c(B,PQ,G,buff,nit,b1,b2,qb,pb,mn,btol,mode)

        real(c_double), intent(out)             :: B(:,:)
        real(c_double), intent(out)             :: PQ(:,:)
        real(dp), intent(out)                   :: G(:,:)
        real(dp), intent(out)                   :: buff(:)
        real(dp), intent(in)                    :: btol
        
        integer, intent(in)                     :: b1
        integer, intent(in)                     :: b2
        integer, intent(out)                    :: nit        
        integer, intent(in)                     :: qb
        integer, intent(in)                     :: pb
        integer(c_int), intent(in)              :: mn
        integer, intent(in)                     :: mode     ! mode=1 buldge row, mode=2 buldge col  

        integer                                 :: ii
        integer                                 :: ii1 
        integer                                 :: ii2
        
        ! buldge row index fixed
        if ( mode==1 ) then
            
            ii  = b2
            ii1 = max(ii - qb,1)
            ii2 = min(ii1 + qb + pb + 1, mn)

            !> Check for a permutation
            if ( abs(B(b1,ii+1)) < btol ) then 

                ! Swap columns or their indices 
                buff(1:((ii2-ii1)+1))   = B(ii1:ii2,ii);
                B(ii1:ii2,ii)           = B(ii1:ii2,(ii+1));
                B(ii1:ii2,ii+1)         = buff(1:((ii2-ii1)+1));
                
                nit                     = nit + 1
                ! PQ(nit,1:2)             = [0_dp,1_dp]
                ! PQ(nit,3:4)             = [ii,ii+1]

                PQ(1:2,nit)             = [0_dp,1_dp]
                PQ(3:4,nit)             = [ii,ii+1]

            else ! zero elements using Givens

                call giv_c(G,B(b1,ii:(ii+1)),2,1)
                ![G,y] = giv(B(b1,ii:(ii+1)));         
                
                buff(1:2)               = B(b1,ii:(ii+1))
                
                B(ii1:ii2,ii:(ii+1))    = matmul(B(ii1:ii2,ii:(ii+1)),G);                
                
                B(b1,ii:(ii+1))         = buff(1:2)

                nit                     = nit + 1;
                ! PQ(nit,1:2)             = [G(1,1),G(1,2)]
                ! PQ(nit,3:4)             = [ii,ii+1]

                PQ(1:2,nit)             = [G(1,1),G(1,2)]
                PQ(3:4,nit)             = [ii,ii+1]
                
                !P(nit,:) = [G(1,1),G(1,2),ii,ii+1];

            end if

        else ! buldge column index fixed

            ii  = b1
            ii1 = max(ii - pb,1)
            ii2 = min(ii1 + qb + pb + 1, mn)
            
            if ( abs(B(ii+1,b2)) < btol ) then 

                ! Swap rows or their indices 
                buff(1:((ii2-ii1)+1))   = B(ii,ii1:ii2);
                B(ii,ii1:ii2)           = B((ii+1),ii1:ii2);
                B(ii+1,ii1:ii2)         = buff(1:((ii2-ii1)+1));
                
                nit                     = nit + 1
                ! PQ(nit,1:2)             = [0_dp,1_dp]
                ! PQ(nit,3:4)             = [ii,ii+1]

                PQ(1:2,nit)             = [0_dp,1_dp]
                PQ(3:4,nit)             = [ii,ii+1]
                
            else 

                call giv_c(G,B(ii:(ii+1),b2),1,1)

                buff(1:2)               = B(ii:(ii+1),b2)

                B(ii:(ii+1),ii1:ii2)    = matmul(G,B(ii:(ii+1),ii1:ii2))

                B(ii:(ii+1),b2)         = buff(1:2)

                nit                     = nit + 1;
                ! PQ(nit,1:2)             = [G(1,1),G(1,2)]
                ! PQ(nit,3:4)             = [ii,ii+1]

                PQ(1:2,nit)             = [G(1,1),G(1,2)]
                PQ(3:4,nit)             = [ii,ii+1]

            end if

        end if

    end subroutine


    !> remove buldge phase2 (called from within phase2)
    subroutine rmbp2(B,PQ,G,buff,nit,b1,b2,qb,pb,mn,btol,mode)

        real(dp), intent(out)           :: B(:,:)
        real(dp), intent(out)           :: PQ(:,:)
        real(dp), intent(out)           :: G(:,:)
        real(dp), intent(out)           :: buff(:)
        real(dp), intent(in)            :: btol
        
        integer, intent(in)             :: b1
        integer, intent(in)             :: b2
        integer, intent(out)            :: nit        
        integer, intent(in)             :: qb
        integer, intent(in)             :: pb
        integer, intent(in)             :: mn
        integer, intent(in)             :: mode     ! mode=1 buldge row, mode=2 buldge col  

        integer                         :: ii
        integer                         :: ii1 
        integer                         :: ii2
        
        ! buldge row index fixed
        if ( mode==1 ) then
            
            ii  = b2 - 1
            ii1 = max(ii + 1 - qb - pb,1)
            ii2 = min(ii1 + qb + pb + 1, mn)

            !> Check for a permutation
            if ( abs(B(b1,ii)) < btol ) then 

                ! Swap columns or their indices 
                buff(1:((ii2-ii1)+1))   = B(ii1:ii2,ii);
                B(ii1:ii2,ii)           = B(ii1:ii2,(ii+1));
                B(ii1:ii2,ii+1)         = buff(1:((ii2-ii1)+1));
                
                nit                     = nit + 1
                PQ(nit,1:2)             = [0_dp,1_dp]
                PQ(nit,3:4)             = [ii,ii+1]

            else ! zero elements using Givens

                call giv(G,B(b1,ii:(ii+1)),2,2)
                ![G,y] = giv(B(b1,ii:(ii+1)));       
                
                buff(1:2) = B(b1,ii:(ii+1))
                
                B(ii1:ii2,ii:(ii+1))    = matmul(B(ii1:ii2,ii:(ii+1)),G);                
                
                B(b1,ii:(ii+1)) = buff(1:2)

                nit                     = nit + 1;
                PQ(nit,1:2)             = [G(1,1),G(1,2)]
                PQ(nit,3:4)             = [ii,ii+1]
                
                !P(nit,:) = [G(1,1),G(1,2),ii,ii+1];

            end if

        else ! buldge column index fixed

            ii  = b1 - 1
            ii1 = max(ii - pb,1)
            ii2 = min(ii1 + qb + pb + 1, mn)
            
            if ( abs(B(ii,b2)) < btol ) then 

                ! Swap rows or their indices 
                buff(1:((ii2-ii1)+1))   = B(ii,ii1:ii2);
                B(ii,ii1:ii2)           = B((ii+1),ii1:ii2);
                B(ii+1,ii1:ii2)         = buff(1:((ii2-ii1)+1));
                
                nit                     = nit + 1
                PQ(nit,1:2)             = [0_dp,1_dp]
                PQ(nit,3:4)             = [ii,ii+1]
                
            else 

                call giv(G,B(ii:(ii+1),b2),1,2)

                buff(1:2) = B(ii:(ii+1),b2)

                B(ii:(ii+1),ii1:ii2)    = matmul(G,B(ii:(ii+1),ii1:ii2))

                B(ii:(ii+1),b2) = buff(1:2)

                nit                     = nit + 1;
                PQ(nit,1:2)             = [G(1,1),G(1,2)]
                PQ(nit,3:4)             = [ii,ii+1]

            end if

        end if

    end subroutine

    !> C interface to remove buldge in phase 2
    subroutine rmbp2_c(B,PQ,G,buff,nit,b1,b2,qb,pb,mn,btol,mode)

        real(c_double), intent(out)             :: B(:,:)
        real(c_double), intent(out)             :: PQ(:,:)
        real(dp), intent(out)                   :: G(:,:)
        real(dp), intent(out)                   :: buff(:)
        real(dp), intent(in)                    :: btol
        
        integer, intent(in)                     :: b1
        integer, intent(in)                     :: b2
        integer, intent(out)                    :: nit        
        integer, intent(in)                     :: qb
        integer, intent(in)                     :: pb
        integer(c_int), intent(in)              :: mn
        integer, intent(in)                     :: mode     ! mode=1 buldge row, mode=2 buldge col  

        integer                                 :: ii
        integer                                 :: ii1 
        integer                                 :: ii2
        
        ! buldge row index fixed
        if ( mode==1 ) then
            
            ii  = b2 - 1
            ii1 = max(ii + 1 - qb - pb,1)
            ii2 = min(ii1 + qb + pb + 1, mn)

            !> Check for a permutation
            if ( abs(B(b1,ii)) < btol ) then 

                ! Swap columns or their indices 
                buff(1:((ii2-ii1)+1))   = B(ii1:ii2,ii);
                B(ii1:ii2,ii)           = B(ii1:ii2,(ii+1));
                B(ii1:ii2,ii+1)         = buff(1:((ii2-ii1)+1));
                
                nit                     = nit + 1
                
                ! PQ(nit,1:2)             = [0_dp,1_dp]
                ! PQ(nit,3:4)             = [ii,ii+1]

                PQ(1:2,nit)             = [0_dp,1_dp]
                PQ(3:4,nit)             = [ii,ii+1]

            else ! zero elements using Givens

                call giv_c(G,B(b1,ii:(ii+1)),2,2)
                ![G,y] = giv(B(b1,ii:(ii+1)));       
                
                buff(1:2) = B(b1,ii:(ii+1))
                
                B(ii1:ii2,ii:(ii+1))    = matmul(B(ii1:ii2,ii:(ii+1)),G);                
                
                B(b1,ii:(ii+1)) = buff(1:2)

                nit                     = nit + 1;
                ! PQ(nit,1:2)             = [G(1,1),G(1,2)]
                ! PQ(nit,3:4)             = [ii,ii+1]

                PQ(1:2,nit)             = [G(1,1),G(1,2)]
                PQ(3:4,nit)             = [ii,ii+1]
                
                !P(nit,:) = [G(1,1),G(1,2),ii,ii+1];

            end if

        else ! buldge column index fixed

            ii  = b1 - 1
            ii1 = max(ii - pb,1)
            ii2 = min(ii1 + qb + pb + 1, mn)
            
            if ( abs(B(ii,b2)) < btol ) then 

                ! Swap rows or their indices 
                buff(1:((ii2-ii1)+1))   = B(ii,ii1:ii2);
                B(ii,ii1:ii2)           = B((ii+1),ii1:ii2);
                B(ii+1,ii1:ii2)         = buff(1:((ii2-ii1)+1));
                
                nit                     = nit + 1
                ! PQ(nit,1:2)             = [0_dp,1_dp]
                ! PQ(nit,3:4)             = [ii,ii+1]

                PQ(1:2,nit)             = [0_dp,1_dp]
                PQ(3:4,nit)             = [ii,ii+1]
                
            else 

                call giv_c(G,B(ii:(ii+1),b2),1,2)

                buff(1:2) = B(ii:(ii+1),b2)

                B(ii:(ii+1),ii1:ii2)    = matmul(G,B(ii:(ii+1),ii1:ii2))

                B(ii:(ii+1),b2) = buff(1:2)

                nit                     = nit + 1;
                ! PQ(nit,1:2)             = [G(1,1),G(1,2)]
                ! PQ(nit,3:4)             = [ii,ii+1]

                PQ(1:2,nit)             = [G(1,1),G(1,2)]
                PQ(3:4,nit)             = [ii,ii+1]

            end if

        end if

    end subroutine


    !>
    !> givapply 
    !> 
    !> Subroutine to apply givens rotations for 4 different cases.
    !> A Matlab style description of the function
    !>
    ! function [X,idx] = givApply2_tst(X,Q,isL,isT,mode)
    !     %givApply2_tst Applying compact Given's rotations to matrix X
    !     %
    !     % [X] = givApply2(X,Q,isL,isT,mode) applies the compactly stored rotations
    !     %   in Q to X. Q is a m x 4 matrix where m represents the number of
    !     %   rotations. The ith row of Q(i,:) = [c,s,i1,i2] contains the rotation
    !     %   coefficients c,s as computed from giv(y) and i1, i2 correspond to the 
    !     %   effected rows or columns. This method permutes row/column indices
    !     %   instead of copys when swaps are possible.
    !     %
    !     %   Eg., when isL=1,isT=1 and mode=1 the result is computing
    !     %
    !     %    (Qm Qm-1 ... Q1)' X   
    !     %
    !     %   INPUTS:  
    !     %   isL     = 0,1 a flag for left transformations,
    !     %   isT     = 0,1 a flag for transpose,
    !     %   mode    = 0,1 a flag for whether zeroes are top-down or bottom-up
    !     %
    !     %   OUTPUS:
    !     %   X       = Product of rotations with original matrix in permuted form
    !     %   idx     = X(idx,:) or X(:,idx) is the result of the orginal matrix
    !     %           transformed
    !     %
    !     %--------------------------------------------------------------------------   
    subroutine givapply(X,Q,idx,isL,isT,mode)

        real(dp),   intent(out)                 :: X(:,:)
        real(dp),   intent(in)                  :: Q(:,:)
        integer,    intent(out)                 :: idx(:)

        integer,    intent(in)                  :: isL
        integer,    intent(in)                  :: isT
        integer,    intent(in)                  :: mode
        
        real(dp),   allocatable                 :: buff(:)

        integer                                 :: i
        integer                                 :: m
        integer                                 :: n
        integer                                 :: ssgn
        integer                                 :: te
                
        te  = int(Q(1,1)) + 1
        m   = size( X, dim=1 )
        n   = size( X, dim=2 )

        if ( mode .eq. 1 ) then 
            ssgn = 1
        else
            ssgn = -1
        endif  

        if ( isL .eq. 1 ) then ! left apply

            allocate( buff(n) )
            !idx(:)  = (/ (i, i=1, m, 1) /)

            if ( isT .eq. 0 ) then ! no transpose
                
                do i = 2, te, 1

                    call givel(i, Q, X, buff, idx, 1, ssgn, ssgn, 0)

                end do

            else ! transpose

                do i = te, 2, -1

                    call givel(i, Q, X, buff, idx, ssgn, 1, ssgn, 0)

                    ! DEBUG
                    ! if ( i .eq. te ) then
                    !     print '("i=", i0)', i
                    !     print *, 'X'
                    !     print *, X
                    !     print *, 'buff'
                    !     print *, buff
                    ! endif

                end do

            endif

        else ! right apply

            allocate( buff(m) )
            !idx(:)  = (/ (i, i=1, n, 1) /)

            if ( isT .eq. 0 ) then ! no transpose
                
                do i = 2, te, 1

                    call givel(i, Q, X, buff, idx, ssgn, 1, ssgn, 1)

                end do

            else ! transpose

                do i = te, 2, -1

                    call givel(i, Q, X, buff, idx, 1, ssgn, ssgn, 1)

                end do

            endif

        endif

    end subroutine givapply

    !>
    !> givel
    !> 
    !> subroutine to apply each givens rotation
    !>
    !> 08/28/25, J.B., modify ztol to machine precision
    subroutine givel(i,Q,X,buff,idx,c1,c2,ssgn,colrow)

        real(dp), intent(in)                    :: Q(:,:)
        real(dp), intent(out)                   :: X(:,:)
        real(dp), intent(out)                   :: buff(:)

        integer, intent(in)                     :: i
        integer, intent(out)                    :: idx(:)

        integer                                 :: c1
        integer                                 :: c2
        integer                                 :: ssgn     ! ssgn = 1 or ssgn = -1
        integer                                 :: colrow   ! colrow = 1/0 for column/row
        integer                                 :: ibuff

        !real(dp)                                :: ztol = 1.0e-14_dp     ! zero tolerance
        real(dp)                                :: ztol = epsilon(1.0_dp) ! zero tolerance

        ! If a permutation is possible then perform an index swap, otherwise
        ! apply the rotation
        if ( abs( Q(1,i) ) > ztol ) then

            ! row givens apply
            if ( colrow .eq. 0 ) then 
            
                buff(:) = Q( 1, i )*X( idx( int(Q(3,i))), : ) + c1*Q( 2, i )*X( idx(int(Q(4,i))), : )

                X( idx(int(Q(4,i))), : ) = c2*Q( 2, i )*X( idx(int(Q(3,i))), :) - ssgn*Q( 1, i )*X( idx(int(Q(4,i))), :)
                
                X( idx(int(Q(3,i))), : ) = buff(:)

                ! DEBUG
                !print *, 'idx'
                !print *, idx
                
            else ! column apply

                buff(:) = Q( 1, i )*X( :, idx( int(Q(3,i))) ) + c1*Q( 2, i )*X( :, idx(int(Q(4,i))) )

                X( :, idx(int(Q(4,i))) ) = c2*Q( 2, i )*X( :, idx(int(Q(3,i)))) - ssgn*Q( 1, i )*X( :, idx(int(Q(4,i))) )
                
                X( :, idx(int(Q(3,i))) ) = buff(:)

            endif

        else

            ibuff               = idx( int(Q(4,i)) )
            idx( int(Q(4,i)) )  = idx( int(Q(3,i)) )
            idx( int(Q(3,i)) )  = ibuff

        endif

    end subroutine givel

    !>
    !> Left multiplications with the compact Qs, i.e.,
    !>
    !> Q1  * Q2  * Q3 * X,  (trans==1)
    !> Q3' * Q2' * Q1 * X,  (trans==0)
    !>
    subroutine bidiag_up_mulq(X,Q1,Q2,Q3,idx,m,n,n1,n2,n3,trans)

        real(dp), intent(out)                   :: X(m,n)
        real(dp), intent(in)                    :: Q1(4,n1)
        real(dp), intent(in)                    :: Q2(4,n2)
        real(dp), intent(in)                    :: Q3(4,n3)

        integer, intent(out)                    :: idx(m)
        integer, intent(in)                     :: m
        integer, intent(in)                     :: n
        integer, intent(in)                     :: n1
        integer, intent(in)                     :: n2
        integer, intent(in)                     :: n3
        integer, intent(in)                     :: trans    ! trans=1/0, transpose/notrans 
        integer                                 :: i

        idx(:)  = (/ (i, i=1, m, 1) /)

        if ( trans .eq. 1 ) then 

            ! DEBUG
            !print *, 'idx'
            !print *, idx

            call givapply(X,Q3,idx,1,1,2)

            ! DEBUG
            !print *, 'idx'
            !print *, idx

            call givapply(X,Q2,idx,1,1,2)

              ! DEBUG
            !print *, 'idx'
            !print *, idx

            call givapply(X,Q1,idx,1,1,1)

        else

            call givapply(X,Q1,idx,1,0,1)
            call givapply(X,Q2,idx,1,0,2)
            call givapply(X,Q3,idx,1,0,2)

        endif
            
    end subroutine bidiag_up_mulq


    !>
    !> C interface
    !>
    subroutine bidiag_up_mulq_c(X,Q1,Q2,Q3,idx,m,n,n1,n2,n3,trans) bind(C,name="bidiag_up_mulq_c") 

        real(c_double), intent(out)                     :: X(m,n)
        real(c_double), intent(in)                      :: Q1(4,n1)
        real(c_double), intent(in)                      :: Q2(4,n2)
        real(c_double), intent(in)                      :: Q3(4,n3)

        integer(c_int), intent(out)                     :: idx(m)
        integer(c_int), intent(in)                      :: m
        integer(c_int), intent(in)                      :: n
        integer(c_int), intent(in)                      :: n1
        integer(c_int), intent(in)                      :: n2
        integer(c_int), intent(in)                      :: n3
        integer(c_int), intent(in)                      :: trans    ! trans=1/0, transpose/notrans 
        !integer                                         :: i

        call bidiag_up_mulq(X,Q1,Q2,Q3,idx,m,n,n1,n2,n3,trans)

        ! idx(:)  = (/ (i, i=1, m, 1) /)

        ! if ( trans .eq. 1 ) then 

        !     ! DEBUG
        !     !print *, 'idx'
        !     !print *, idx

        !     call givapply(X,Q3,idx,1,1,2)

        !     ! DEBUG
        !     !print *, 'idx'
        !     !print *, idx

        !     call givapply(X,Q2,idx,1,1,2)

        !       ! DEBUG
        !     !print *, 'idx'
        !     !print *, idx

        !     call givapply(X,Q1,idx,1,1,1)

        ! else

        !     call givapply(X,Q1,idx,1,0,1)
        !     call givapply(X,Q2,idx,1,0,2)
        !     call givapply(X,Q3,idx,1,0,2)

        ! endif
            
    end subroutine bidiag_up_mulq_c

    !>
    !> Right multiplications with the compact Ps, i.e.,
    !>
    !> X * P2  * P1 ,  (trans==1)
    !> X * P1' * P2',  (trans==0)
    !>
    subroutine bidiag_up_mulp(X,Q1,Q2,idx,m,n,n1,n2,trans)

        real(dp), intent(out)                   :: X(m,n)
        real(dp), intent(in)                    :: Q1(4,n1)
        real(dp), intent(in)                    :: Q2(4,n2)        

        integer, intent(out)                    :: idx(n)
        integer, intent(in)                     :: m
        integer, intent(in)                     :: n
        integer, intent(in)                     :: n1
        integer, intent(in)                     :: n2
        integer, intent(in)                     :: trans    ! trans=1/0, transpose/notrans 
        integer                                 :: i

        idx(:)  = (/ (i, i=1, n, 1) /)

        if ( trans .eq. 1 ) then 
            
            call givapply(X,Q2,idx,0,1,2)
            
            ! DEBUG
            !print *, 'idx'
            !print *, idx
            
            call givapply(X,Q1,idx,0,1,1)

            ! DEBUG
            !print *, 'idx'
            !print *, idx

        else

            call givapply(X,Q1,idx,0,0,1)
            call givapply(X,Q2,idx,0,0,2)            

        endif
            
    end subroutine bidiag_up_mulp

    !>
    !> C interface
    !>    
    subroutine bidiag_up_mulp_c(X,Q1,Q2,idx,m,n,n1,n2,trans) bind(C,name="bidiag_up_mulp_c")

        real(c_double), intent(out)                    :: X(m,n)
        real(c_double), intent(in)                     :: Q1(4,n1)
        real(c_double), intent(in)                     :: Q2(4,n2)        

        integer(c_int), intent(out)                    :: idx(n)
        integer(c_int), intent(in)                     :: m
        integer(c_int), intent(in)                     :: n
        integer(c_int), intent(in)                     :: n1
        integer(c_int), intent(in)                     :: n2
        integer(c_int), intent(in)                     :: trans    ! trans=1/0, transpose/notrans 
        !integer                                 :: i

        call bidiag_up_mulp(X,Q1,Q2,idx,m,n,n1,n2,trans)

        ! idx(:)  = (/ (i, i=1, n, 1) /)

        ! if ( trans .eq. 1 ) then 
            
        !     call givapply(X,Q2,idx,0,1,2)
            
        !     ! DEBUG
        !     !print *, 'idx'
        !     !print *, idx
            
        !     call givapply(X,Q1,idx,0,1,1)

        !     ! DEBUG
        !     !print *, 'idx'
        !     !print *, idx

        ! else

        !     call givapply(X,Q1,idx,0,0,1)
        !     call givapply(X,Q2,idx,0,0,2)            

        ! endif
            
    end subroutine bidiag_up_mulp_c

end module