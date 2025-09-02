!>
!> test_bdup
!>
!> First test for updating a bidiagonal factorization
!> Extended to include the givapply method
!>
!>----------------------------------------------------------------------------------------
!> 07/01/25, J.B., Initial version
!> 09/01/25, J.B., 
!> 10/01/25, J.B., testing both phases
!> 02/28/25, J.B., testing givapply

program test_bdup

    use kind_parameter, only: dp
    use bidiag_up

    implicit none

    integer, parameter      :: m = 6
    integer, parameter      :: n = 3
    
    real(dp)                :: B(m,n)
    real(dp)                :: A(m,n)
    real(dp)                :: C(m,n)
    real(dp)                :: w(m)
    real(dp)                :: p(n)
    !real(dp)    :: p1(2)
    real(dp)                :: In(n,n)
    real(dp)                :: Im(m,m)
    real(dp)                :: In2(n,n)
    real(dp)                :: Im2(m,m)
    integer                 :: idxm(m)
    integer                 :: idxn(n)    

    real(dp), allocatable   :: Q1(:,:)
    real(dp), allocatable   :: Q2(:,:)
    real(dp), allocatable   :: Q3(:,:)
    real(dp), allocatable   :: P1(:,:)
    real(dp), allocatable   :: P2(:,:)

    real(dp)                :: errm
    real(dp)                :: errn
    integer                 :: i
    integer                 :: j

    B(:,:)  = 0
    B(1,1)  = 1
    B(1,2)  = 2
    B(2,2)  = 2
    B(2,3)  = 3
    B(3,3)  = 3

    w(:)    = [2_dp,1_dp,1_dp,1_dp,1_dp,1_dp]
    p(:)    = [2_dp,0_dp,1_dp]
    !p1(:)   = [1_dp,0_dp]

    C(:,:)  = B + matmul(reshape(w,[size(w),1]),reshape(p,[1,size(p)]))

    print *, "C"
    print *, C

    In = 0
    Im = 0
    In2 = 0
    Im2 = 0

    do i = 1, n
        In(i,i)     = 1.0
        In2(i,i)    = 1.0
    enddo
    do i = 1, m
        Im(i,i)     = 1.0
        Im2(i,i)    = 1.0
    enddo

    !>
    !> Call to the functions
    !>
    call phase1(B,w,p,Q1,Q2,P1)

    call phase2(B,Q3,P2)

    A(:,:) = B(:,:)
    
    print *, 'A'
    print *, A

    ! Applying transposed reflectors
    call bidiag_up_mulq(A,transpose(Q1),transpose(Q2),transpose(Q3),idxm,m,n,int(Q1(1,1)),int(Q2(1,1)),int(Q3(1,1)),1)
    call bidiag_up_mulp(A,transpose(P1),transpose(P2),idxn,m,n,int(P1(1,1)),int(P2(1,1)),1)

    print *, "A(idxm,idxn)"
    print *, A(idxm,idxn)

    ! Applying reflectors
    call bidiag_up_mulq(C,transpose(Q1),transpose(Q2),transpose(Q3),idxm,m,n,int(Q1(1,1)),int(Q2(1,1)),int(Q3(1,1)),0)
    call bidiag_up_mulp(C,transpose(P1),transpose(P2),idxn,m,n,int(P1(1,1)),int(P2(1,1)),0)

    print *, "C(idxm,idxn)"
    print *, C(idxm,idxn)

    ! Computing Q
    call bidiag_up_mulq(Im2,transpose(Q1),transpose(Q2),transpose(Q3),idxm,m,m,int(Q1(1,1)),int(Q2(1,1)),int(Q3(1,1)),1)
    print *, "Im2(idxm,:)"
    print *, Im2(idxm,:)

    ! Computing P
    call bidiag_up_mulp(In,transpose(P1),transpose(P2),idxn,n,n,int(P1(1,1)),int(P2(1,1)),1)
    print *, "In(:,idxn)"
    print *, In(:,idxn)

end program