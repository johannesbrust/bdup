!>
!> bidiagh_up
!>
!> Algorithm for updating a bidiagonal factorization using Householder reflectors
!> This is for overdetermined systems, i.e., numrows >= numcols
!>
!> Given B + wh*p' this method finds sequences of Householder transformations so that
!> Bk = Qk'*(B + wh*p')*Pk at k=n is a new bidiagonal.
!>
!> The orthogonal matrices Q and P are not explicitly computed, only the Householder
!> reflectors are stored. This implementation uses a compact representation.
!> In particular, denote the unit norm left and right Householder vectors as yk and wk
!>
!> Yk = [y1,...,yk]
!> Wk = [w1,...,wk]
!> Lk = I + 2 tril(Yk'*Yk,-1)
!> Tk = I + 2 triu(Wk'*Wk,1)
!>
!> The compact representation of Bk is 
!>
!>                             [ -2        0   2p'*W   ]     [ p'   ]
!> Bk = B -2 [wh Y B*W] * inv( [ 2Y'*wh    L   2Y'*B*W ] ) * [ Y'*B ]
!>                             [0          0       T   ]     [ W'   ]  
!> 
!> This module contains functions:
!>
!> trapmul
!> bidagh_up
!> triangsolve
!>----------------------------------------------------------------------------------------
!> 05/02/25, J.B., Initial implementation
!> 06/02/25, J.B., C interface functions
!> 02/12/25, J.B., Reordering loops for more data congruency

module bidiagh_up

    use iso_c_binding
    use kind_parameter, only: dp

    implicit none 

    contains

    ! Main algorithm : bhu (bidiagonal-householder-update)
    subroutine bhu(Y,W,B,BO,wh,p,m,n,kk)
        
        real(dp), intent(out)           :: Y(m,kk)  ! left householder reflectors
        real(dp), intent(out)           :: W(n,kk)  ! right householder reflectors
        real(dp), intent(in)            :: B(n,2)   ! bidiagonal in
        real(dp), intent(out)           :: BO(kk,2) ! bidiagonal out
        
        real(dp), intent(in)            :: wh(m)    ! left rank-1 update
        real(dp), intent(in)            :: p(n)     ! right rank-1 update
        
        integer, intent(in)             :: m        
        integer, intent(in)             :: n
        integer, intent(in)             :: kk       ! number of iterations

        integer                         :: k

        real(dp)                        :: a1(2*n+1) ! buffer 1
        real(dp)                        :: a2(2*n+1) ! buffer 2 
        real(dp)                        :: a3(n)     ! buffer 3
        real(dp)                        :: Yw(kk)
        real(dp)                        :: Wp(kk)
        real(dp)                        :: alpbet
        real(dp)                        :: y1(m)
        real(dp)                        :: w1(n)

        a1      = 0        
        BO      = 0

        ! first update
        k       = 1

        !>
        !> left householder vector
        !>
        y1(k:m) = p(k) * wh(k:m)

        call hvec(y1,alpbet,B(k,1),k,m)
        
        ! updates
        BO(k,1) = alpbet
        Y(k:m,k)= y1(k:m)
        Yw(k)   = dot_product(y1(k:m),wh(k:m))

        !>
        !> right householder vector
        !>
        a1(2)       = Y(k,k)
        a1(1)       = (2*Yw(k)*a1(2) -  wh(k)) / 2
        
        w1(k+1:n)   = -2*a1(1) * p(k+1:n) - 2*a1(2)*( B(k+1:n,1) * Y(k+1:n,k) + B(k:n,2) * Y(1:n-1,k) )
        
        ! Debug
        !print *, "w1(k+1:n)"
        !print *, w1(k+1:n)

        call hvec(w1,alpbet,B(k,2),k+1,n)
        
        BO(k,2)     = alpbet
        W(k+1:n,k)  = w1(k+1:n)
        Wp(k)       = dot_product(w1(k+1:n),p(k+1:n))

        ! start loop here
        do k = 2, kk 

            !>
            !> [ p'  ] * ek
            !> [ Y'B ]
            !>
            a1(1)   = p(k)
            a1(2:k) = B(k-1,2)*Y(k-1,1:k-1) + B(k,1)*Y(k,1:k-1)

            !> Solution to block triangular system
            !>
            !>  [ -2        0   2p'*W   ] [ a2(1)        ] 
            !>  [ 2Y'*wh    L   2Y'*B*W ] [ a2(2:k)      ] 
            !>  [0          0       T   ] [ a2(k+1:2k-1) ]
            !>
            !>

            !> triangular solve for a2(k+1:2k-1) 
            call trisolve(W,W(k,1:k-1),a2(k+1:(2*k-1)),k-1,0)

            ! solve for a2(1)
            a2(1)   = (a1(1) - 2*dot_product(Wp(1:k-1),a2(k+1:(2*k-1)))) / (-2.0)
            !a2(1)   = a1(1)

            ! solve for a2(2:k)
            call trapmul(W,a2(k+1:(2*k-1)),a1(n+1:2*n),n,k-1,1,0)

            !a1(n+1:2*n) = B(1:n,1)*a1(n+1:2*n) + B(1:n,2)*a1(n+2:(2*n+1))
            a3(1:n) = B(1:n,1)*a1(n+1:2*n) + B(1:n,2)*a1(n+2:(2*n+1))

            ! DEBUG
            !print *, "a3(1:n)"
            !print *, a3(1:n)
            
            !print *, "Y"
            !print *, Y

            !call trapmul(Y,a1(n+1:2*n),a1(k+1:(2*k-1)),n,k-1,0,1)
            call trapmul(Y,a3(1:n),a1(k+1:(2*k-1)),n,k-1,0,1)
            
            !print *, "a1(k+1:(2*k-1))"
            !print *, a1(k+1:(2*k-1))

            a1(2:k) = a1(2:k) - 2 * a1(k+1:(2*k-1)) - 2 * a2(1) * Yw(1:k-1)
            
            ! DEBUG
            !print *, "a1(2:k)"
            !print *, a1(2:k)

            call trisolve(Y,a1(2:k),a2(2:k),k-1,1)

            ! DEBUG
            !print *, "a2(1:2k-1)"
            !print *, a2(1:2*k-1)

            !>
            !> -2 [ wh Y B*W ] * [ a2(1) a2(2:k) a2(k+1:2k-1) ]'
            !>
            a1(n+1:2*n-k+1) = matmul(W(k:n,1:(k-1)),-2*a2(k+1:(2*k-1)))

            a1(n+1:2*n-k+1) = B(k:n,1) * a1(n+1:2*n-k+1) + B(k:n,2) * a1(n+2:2*n-k+2)

            y1(k:m)         = wh(k:m)*(-2*a2(1)) + matmul(Y(k:m,1:k-1),-2*a2(2:k))

            y1(k:n)         = y1(k:n) + a1(n+1:2*n-k+1)

            ! householder vector
            call hvec(y1,alpbet,B(k,1),k,m)

            ! updates
            BO(k,1)     = alpbet
            Y(k:m,k)    = y1(k:m)
            Y(1:k-1,k)  = 2*matmul(y1(k:m),Y(k:m,1:k-1)) 
            Yw(k)       = dot_product(y1(k:m),wh(k:m))

            ! early exit
            if ( (kk == 2) .or. (k == n) ) then
                exit
            end if

            ! right transform

            a1(1)   = wh(k)
            a1(2:k) = B(k,2)*W(k+1,1:k-1) + B(k,1)*W(k,1:k-1)

            !> Solution to transposed block triangular system
            !>
            !>         [ -2        0   2p'*W   ]   [ a2(1)        ] 
            !>  trans( [ 2Y'*wh    L   2Y'*B*W ] ) [ a2(2:k)      ] 
            !>         [0          0       T   ]   [ a2(k+1:2k-1) ]
            !>
            !>

            !> triangular solve for a2(k+1:2k-1) 
            call trisolve(Y,Y(k,1:k),a2(k+1:2*k),k,0)

            !print *, "Y(k,1:k)"
            !print *, Y(k,1:k)

            !print *, "a2(k+1:2*k)"
            !print *, a2(k+1:2*k)

            !print *, "Y"
            !print *, Y

            ! solve for a2(1)
            a2(1)   = (a1(1) - 2*dot_product(Yw(1:k),a2(k+1:2*k))) / (-2.0)

            !print *, "a1(1)"
            !print *, a1(1)

            !print *, "a2(1)"
            !print *, a2(1)

            !print *, "dot_product(Yw(1:k),a2(k+1:2*k))"
            !print *, dot_product(Yw(1:k),a2(k+1:2*k))

            call trapmul(Y,a2(k+1:2*k),a3(1:n),n,k,0,0)

            a3(2:n) = B(2:n,1) * a3(2:n) + B(1:n-1,2) * a3(1:n-1)
            a3(1)   = B(1,1) * a3(1)

            !print *, "a3(1:n)"
            !print *, a3(1:n)

            !a3(1:n) = B(1:n,1)*a3(1:n)
            !a3(2:n) = a3(2:n) + B(1:n-1,2) 

            call trapmul(W,a3(1:n),a1(k+1:(2*k-1)),n,k-1,1,1)

            a1(2:k) = a1(2:k) - 2 * a1(k+1:(2*k-1)) - 2 * a2(1) * Wp(1:k-1)

            call trisolve(W,a1(2:k),a2(2:k),k-1,1)

            !print *, "a2(1:2k)"
            !print *, a2(1:2*k)


            !>
            !> copied from left hand update
            !>

            a1(n+1:2*n-k+1) = matmul(Y(k:n,1:k),-2*a2(k+1:(2*k)))

            a1(n+1:2*n-k)   = B(k+1:n,1) * a1(n+2:2*n-k+1) + B(k:n-1,2) * a1(n+1:2*n-k)

            w1(k+1:n)       = p(k+1:n)*(-2*a2(1)) + matmul(W(k+1:n,1:k-1),-2*a2(2:k)) 

            w1(k+1:n)       = w1(k+1:n) + a1(n+1:2*n-k)

            ! householder vector
            call hvec(w1,alpbet,B(k,2),k+1,n)

            ! updates
            BO(k,2)     = alpbet
            w(k+1:n,k)    = w1(k+1:n)
            W(1:k-1,k)  = 2*matmul(w1(k+1:n),W(k+1:n,1:k-1)) 
            Wp(k)       = dot_product(w1(k+1:n),p(k+1:n))

        end do
        
    end subroutine

    ! C interface
    subroutine bhu_c(Y,W,B,BO,wh,p,m,n,kk) bind(C,name="bhu_c")
        
        real(c_double), intent(out)             :: Y(m,kk)  ! left householder reflectors
        real(c_double), intent(out)             :: W(n,kk)  ! right householder reflectors
        real(c_double), intent(in)              :: B(n,2)   ! bidiagonal in
        real(c_double), intent(out)             :: BO(kk,2) ! bidiagonal out
        
        real(c_double), intent(in)              :: wh(m)    ! left rank-1 update
        real(c_double), intent(in)              :: p(n)     ! right rank-1 update
        
        integer(c_int), intent(in)              :: m        
        integer(c_int), intent(in)              :: n
        integer(c_int), intent(in)              :: kk       ! number of iterations

        integer                                 :: k

        real(c_double)                          :: a1(2*n+1) ! buffer 1
        real(c_double)                          :: a2(2*n+1) ! buffer 2 
        real(c_double)                          :: a3(n)     ! buffer 3
        real(c_double)                          :: Yw(kk)
        real(c_double)                          :: Wp(kk)
        real(c_double)                          :: alpbet
        real(c_double)                          :: y1(m)
        real(c_double)                          :: w1(n)

        a1      = 0        
        BO      = 0

        ! first update
        k       = 1

        !>
        !> left householder vector
        !>
        y1(k:m) = p(k) * wh(k:m)

        call hvec_c(y1,alpbet,B(k,1),k,m)
        
        ! updates
        BO(k,1) = alpbet
        Y(k:m,k)= y1(k:m)
        Yw(k)   = dot_product(y1(k:m),wh(k:m))

        !>
        !> right householder vector
        !>
        a1(2)       = Y(k,k)
        a1(1)       = (2*Yw(k)*a1(2) -  wh(k)) / 2
        
        w1(k+1:n)   = -2*a1(1) * p(k+1:n) - 2*a1(2)*( B(k+1:n,1) * Y(k+1:n,k) + B(k:n,2) * Y(1:n-1,k) )
        
        ! Debug
        !print *, "w1(k+1:n)"
        !print *, w1(k+1:n)

        call hvec_c(w1,alpbet,B(k,2),k+1,n)
        
        BO(k,2)     = alpbet
        W(k+1:n,k)  = w1(k+1:n)
        Wp(k)       = dot_product(w1(k+1:n),p(k+1:n))

        ! start loop here
        do k = 2, kk 

            !>
            !> [ p'  ] * ek
            !> [ Y'B ]
            !>
            a1(1)   = p(k)
            a1(2:k) = B(k-1,2)*Y(k-1,1:k-1) + B(k,1)*Y(k,1:k-1)

            !> Solution to block triangular system
            !>
            !>  [ -2        0   2p'*W   ] [ a2(1)        ] 
            !>  [ 2Y'*wh    L   2Y'*B*W ] [ a2(2:k)      ] 
            !>  [0          0       T   ] [ a2(k+1:2k-1) ]
            !>
            !>

            !> triangular solve for a2(k+1:2k-1) 
            call trisolve_c(W,W(k,1:k-1),a2(k+1:(2*k-1)),k-1,0)

            ! solve for a2(1)
            a2(1)   = (a1(1) - 2*dot_product(Wp(1:k-1),a2(k+1:(2*k-1)))) / (-2.0)
            !a2(1)   = a1(1)

            ! solve for a2(2:k)
            call trapmul_c(W,a2(k+1:(2*k-1)),a1(n+1:2*n),n,k-1,1,0)

            !a1(n+1:2*n) = B(1:n,1)*a1(n+1:2*n) + B(1:n,2)*a1(n+2:(2*n+1))
            a3(1:n) = B(1:n,1)*a1(n+1:2*n) + B(1:n,2)*a1(n+2:(2*n+1))

            ! DEBUG
            !print *, "a3(1:n)"
            !print *, a3(1:n)
            
            !print *, "Y"
            !print *, Y

            !call trapmul(Y,a1(n+1:2*n),a1(k+1:(2*k-1)),n,k-1,0,1)
            call trapmul_c(Y,a3(1:n),a1(k+1:(2*k-1)),n,k-1,0,1)
            
            !print *, "a1(k+1:(2*k-1))"
            !print *, a1(k+1:(2*k-1))

            a1(2:k) = a1(2:k) - 2 * a1(k+1:(2*k-1)) - 2 * a2(1) * Yw(1:k-1)
            
            ! DEBUG
            !print *, "a1(2:k)"
            !print *, a1(2:k)

            call trisolve_c(Y,a1(2:k),a2(2:k),k-1,1)

            ! DEBUG
            !print *, "a2(1:2k-1)"
            !print *, a2(1:2*k-1)

            !>
            !> -2 [ wh Y B*W ] * [ a2(1) a2(2:k) a2(k+1:2k-1) ]'
            !>
            a1(n+1:2*n-k+1) = matmul(W(k:n,1:(k-1)),-2*a2(k+1:(2*k-1)))

            a1(n+1:2*n-k+1) = B(k:n,1) * a1(n+1:2*n-k+1) + B(k:n,2) * a1(n+2:2*n-k+2)

            y1(k:m)         = wh(k:m)*(-2*a2(1)) + matmul(Y(k:m,1:k-1),-2*a2(2:k))

            y1(k:n)         = y1(k:n) + a1(n+1:2*n-k+1)

            ! householder vector
            call hvec_c(y1,alpbet,B(k,1),k,m)

            ! updates
            BO(k,1)     = alpbet
            Y(k:m,k)    = y1(k:m)
            Y(1:k-1,k)  = 2*matmul(y1(k:m),Y(k:m,1:k-1)) 
            Yw(k)       = dot_product(y1(k:m),wh(k:m))

            ! early exit
            if ( (kk == 2) .or. (k == n) ) then
                exit
            end if

            ! right transform

            a1(1)   = wh(k)
            a1(2:k) = B(k,2)*W(k+1,1:k-1) + B(k,1)*W(k,1:k-1)

            !> Solution to transposed block triangular system
            !>
            !>         [ -2        0   2p'*W   ]   [ a2(1)        ] 
            !>  trans( [ 2Y'*wh    L   2Y'*B*W ] ) [ a2(2:k)      ] 
            !>         [0          0       T   ]   [ a2(k+1:2k-1) ]
            !>
            !>

            !> triangular solve for a2(k+1:2k-1) 
            call trisolve_c(Y,Y(k,1:k),a2(k+1:2*k),k,0)

            !print *, "Y(k,1:k)"
            !print *, Y(k,1:k)

            !print *, "a2(k+1:2*k)"
            !print *, a2(k+1:2*k)

            !print *, "Y"
            !print *, Y

            ! solve for a2(1)
            a2(1)   = (a1(1) - 2*dot_product(Yw(1:k),a2(k+1:2*k))) / (-2.0)

            !print *, "a1(1)"
            !print *, a1(1)

            !print *, "a2(1)"
            !print *, a2(1)

            !print *, "dot_product(Yw(1:k),a2(k+1:2*k))"
            !print *, dot_product(Yw(1:k),a2(k+1:2*k))

            call trapmul_c(Y,a2(k+1:2*k),a3(1:n),n,k,0,0)

            a3(2:n) = B(2:n,1) * a3(2:n) + B(1:n-1,2) * a3(1:n-1)
            a3(1)   = B(1,1) * a3(1)

            !print *, "a3(1:n)"
            !print *, a3(1:n)

            !a3(1:n) = B(1:n,1)*a3(1:n)
            !a3(2:n) = a3(2:n) + B(1:n-1,2) 

            call trapmul_c(W,a3(1:n),a1(k+1:(2*k-1)),n,k-1,1,1)

            a1(2:k) = a1(2:k) - 2 * a1(k+1:(2*k-1)) - 2 * a2(1) * Wp(1:k-1)

            call trisolve_c(W,a1(2:k),a2(2:k),k-1,1)

            !print *, "a2(1:2k)"
            !print *, a2(1:2*k)


            !>
            !> copied from left hand update
            !>

            a1(n+1:2*n-k+1) = matmul(Y(k:n,1:k),-2*a2(k+1:(2*k)))

            a1(n+1:2*n-k)   = B(k+1:n,1) * a1(n+2:2*n-k+1) + B(k:n-1,2) * a1(n+1:2*n-k)

            w1(k+1:n)       = p(k+1:n)*(-2*a2(1)) + matmul(W(k+1:n,1:k-1),-2*a2(2:k)) 

            w1(k+1:n)       = w1(k+1:n) + a1(n+1:2*n-k)

            ! householder vector
            call hvec_c(w1,alpbet,B(k,2),k+1,n)

            ! updates
            BO(k,2)     = alpbet
            w(k+1:n,k)    = w1(k+1:n)
            W(1:k-1,k)  = 2*matmul(w1(k+1:n),W(k+1:n,1:k-1)) 
            Wp(k)       = dot_product(w1(k+1:n),p(k+1:n))

        end do
        
    end subroutine

    ! computes a householder vector
    subroutine hvec(wy,alpbet,Bk,k,mn)

        real(dp), intent(out)       :: wy(:)
        real(dp), intent(out)       :: alpbet
        real(dp), intent(in)        :: Bk
        integer, intent(in)         :: k
        integer, intent(in)         :: mn

        real(dp)                    :: a1

        real(dp)                    :: btol = 1.0e-13_dp

        wy(k)       = wy(k) + Bk
        a1          = wy(k)
        alpbet      = norm2(wy(k:mn))

        if (alpbet > btol) then

            wy(k)       = wy(k) + sign(1.0d0,a1) * alpbet        
            wy(k:mn)    = wy(k:mn) / sqrt(2 * alpbet * (alpbet + abs(a1)))
            alpbet      = sign(alpbet,-a1)
        
        
        else
            wy(k:mn)    = 0.0
        
        end if

    end subroutine

    ! C interface
    subroutine hvec_c(wy,alpbet,Bk,k,mn)

        real(c_double), intent(out)         :: wy(:)
        real(c_double), intent(out)         :: alpbet
        real(c_double), intent(in)          :: Bk
        integer, intent(in)                 :: k
        integer(c_int), intent(in)          :: mn

        real(c_double)                      :: a1

        real(c_double)                      :: btol = 1.0e-13

        wy(k)       = wy(k) + Bk
        a1          = wy(k)
        alpbet      = norm2(wy(k:mn))
        
        if (alpbet > btol) then

            wy(k)       = wy(k) + sign(1.0d0,a1) * alpbet        
            wy(k:mn)    = wy(k:mn) / sqrt(2 * alpbet * (alpbet + abs(a1)))
            alpbet      = sign(alpbet,-a1)
        
        
        else
            wy(k:mn)    = 0.0
        
        end if

    end subroutine

    !> Forming products with trapezoidal matrices
    !>
    !> Yk = [ y11         ]         Wk = [ 0            ] 
    !>      [ y21 y22     ]              [ w21          ]
    !>      [ y31 y32 y33 ]              [ w31 w32      ]
    !>      [ y41 y42 y43 ]              [ w41 w42  w43 ]
    !>      [ y51 y52 y53 ]              [ w51 w52  w53 ]
    !>      [ y61 y62 y63 ]
    subroutine trapmul(WY,a,b,n,km1,isWY,trans)

        real(dp), intent(in)    :: WY(:,:)
        real(dp), intent(in)    :: a(:)
        real(dp), intent(out)   :: b(:)
        integer,  intent(in)    :: n
        integer,  intent(in)    :: km1
        integer,  intent(in)    :: isWY     ! isWY=1 then WY=W, isWY=0 then WY=Y
        integer,  intent(in)    :: trans    ! flag to transpose

        ! offset and indices 
        integer                 :: i1
        integer                 :: i
        integer                 :: j  

        ! check for offset, and initial zero
        if (isWY == 1) then 
            i1 = 1            
        else
            i1 = 0
        end if
        b(:) = 0.0

        ! main method
        if (trans == 0) then

            do j = 1, km1

                b((j+i1):n) = b((j+i1):n) + WY((j+i1):n,j)*a(j)

            end do

            !>
            !> Alternative loop
            !>

            ! do j = 1, km1
            !     do i = j, km1 
            !         b(i+i1) = b(i+i1) + WY(i+i1,j)*a(j)
            !     end do
            ! end do

            ! if (km1 < n) then 
            !     b((km1+1+i1):n) = matmul(WY((km1+1+i1):n,1:km1),a(1:km1))                
            ! end if

        else ! transposed loop

            do j = 1, km1 

                do i = j, (n - i1) 

                    b(j) = b(j) + WY(i+i1,j)*a(i+i1)

                end do

            end do

            !>
            !> Alternative loop
            !>

            ! do j = km1, 1, -1 
            !     do i = j, 1, -1 
            !         b(i) = b(i) + WY(j+i1,i)*a(j+i1)
            !     end do
            ! end do

            ! if (km1 < n) then                 
            !     b(1:km1) = b(1:km1) + matmul(transpose(WY((km1+1+i1):n,1:km1)),a(km1+1+i1:n))
            ! end if

        end if

    end subroutine

    ! Interface for C data
    subroutine trapmul_c(WY,a,b,n,km1,isWY,trans)

        real(c_double), intent(in)      :: WY(:,:)
        real(c_double), intent(in)      :: a(:)
        real(c_double), intent(out)     :: b(:)
        integer(c_int), intent(in)      :: n
        integer,        intent(in)      :: km1
        integer,        intent(in)      :: isWY     ! isWY=1 then WY=W, isWY=0 then WY=Y
        integer,        intent(in)      :: trans    ! flag to transpose

        ! offset and indices 
        integer                         :: i1
        integer                         :: i
        integer                         :: j  

        ! check for offset, and initial zero
        if (isWY == 1) then 
            i1 = 1            
        else
            i1 = 0
        end if
        b(:) = 0.0

        ! main method
        if (trans == 0) then

            do j = 1, km1

                b((j+i1):n) = b((j+i1):n) + WY((j+i1):n,j)*a(j)

            end do

            !>
            !> Alternative loop
            !>

            ! do j = 1, km1
            !     do i = j, km1 
            !         b(i+i1) = b(i+i1) + WY(i+i1,j)*a(j)
            !     end do
            ! end do

            ! if (km1 < n) then 
            !     b((km1+1+i1):n) = matmul(WY((km1+1+i1):n,1:km1),a(1:km1))                
            ! end if

        else ! transposed loop

            do j = 1, km1 

                do i = j, (n - i1) 

                    b(j) = b(j) + WY(i+i1,j)*a(i+i1)
  
                end do

            end do

            !>
            !> Alternative loop
            !>
            ! do j = km1, 1, -1 
            !     do i = j, 1, -1 
            !         b(i) = b(i) + WY(j+i1,i)*a(j+i1)
            !     end do
            ! end do

            ! if (km1 < n) then 
            !     b(1:km1) = b(1:km1) + matmul(transpose(WY((km1+1+i1):n,1:km1)),a(km1+1+i1:n))
            ! end if

        end if

    end subroutine

    !> triangular solves, when non unit element of the triangular matrices are 
    !> stored in the zeros of upper trapezoidal matrices
    !>
    !> Yk = [ y11 l21 l31 ]         Wk = [ 0   r12  r13 ] 
    !>      [ y21 y22 l32 ]              [ w21 0    r23 ]
    !>      [ y31 y32 y33 ]              [ w31 w32  0   ]
    !>      [ y41 y42 y43 ]              [ w41 w42  w43 ]
    !>      [ y51 y52 y53 ]              [ w51 w52  w53 ]
    !>      [ y61 y62 y63 ]
    !> 
    !> Yk stores a lower triangular matrix in its upper triangle. The diagonals are 1.
    !> Wk stores a upper triangular matrix in its upper triangle. The diagonals are 1.
    !>
    subroutine trisolve(WY,a,b,k,trans)

        real(dp), intent(in)        :: WY(:,:)
        real(dp), intent(in)        :: a(:)
        real(dp), intent(out)       :: b(:)
        integer, intent(in)         :: trans ! flag for transpose

        integer, intent(in)         :: k
        integer                     :: i,j

        b(:) = a(:)

        if (trans == 0) then

            do j = k, 1, -1 

                do i = j-1, 1, -1

                    b(i) = b(i) - b(j)*WY(i,j)

                end do

            end do

        else

            do j = 1, k, 1

                do i = 1, (j-1), 1

                    b(j) = b(j) - b(i)*WY(i,j)

                end do

            end do

            !>
            !> Alternative loop
            !>
            ! do j = 1, k, 1
            !     do i = j+1, k, 1
            !         b(i) = b(i) - b(j)*WY(j,i)
            !     end do
            ! end do

        end if

    end subroutine

    ! Interface for C data
    subroutine trisolve_c(WY,a,b,k,trans)

        real(c_double), intent(in)          :: WY(:,:)
        real(c_double), intent(in)          :: a(:)
        real(c_double), intent(out)         :: b(:)
        integer, intent(in)                 :: trans ! flag for transpose

        integer, intent(in)                 :: k
        integer                             :: i,j

        b(:) = a(:)

        if (trans == 0) then

            do j = k, 1, -1 

                do i = j-1, 1, -1

                    b(i) = b(i) - b(j)*WY(i,j)

                end do

            end do

        else

            do j = 1, k, 1

                do i = 1, (j-1), 1

                    b(j) = b(j) - b(i)*WY(i,j)

                end do

            end do

            !>
            !> Alternative loop
            !>
            ! do j = 1, k, 1
            !     do i = j+1, k, 1
            !         b(i) = b(i) - b(j)*WY(j,i)
            !     end do
            ! end do

        end if

    end subroutine


end module
