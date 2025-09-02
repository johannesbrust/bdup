!>
!> givens
!>
!>
!> Module to apply givens rotations. There are two subroutines, planerot and
!> giv. giv is a wrapper around planerot, and is called from the main algorithms
!>
!>

module givens

    use iso_c_binding, only: c_double, c_int
    use kind_parameter, only: dp

    implicit none
    
    private
    
    public planerot, giv, planerot_c, giv_c

    contains

    !>
    !> This subroutine is a translation of the corresponding build-in MATLAB function
    subroutine planerot( G, x )

        real(dp), intent(out)   :: G(:,:)
        real(dp), intent(out)   :: x(:)

        real(dp)                :: zerotol      = 1.0e-14_dp
        real(dp)                :: r

        if ( abs(x(2)) > zerotol ) then

            r = norm2(x)
            G(1,:) = x(:)
            G(2,1) = -x(2)
            G(2,2) = x(1)
            G = G / r
            x(1) = r
            x(2) = 0.0_dp

        else

            G(:,:) = 0.0_dp
            G(1,1) = 1.0_dp
            G(2,2) = 1.0_dp

        end if

    end

    !>
    !> C interface
    !> This subroutine is a translation of the corresponding build-in MATLAB function
    subroutine planerot_c( G, x )

        real(c_double), intent(out)   :: G(:,:)
        real(c_double), intent(out)   :: x(:)

        real(c_double)                :: zerotol      = 1.0e-14
        real(c_double)                :: r

        if ( abs(x(2)) > zerotol ) then

            r = norm2(x)
            G(1,:) = x(:)
            G(2,1) = -x(2)
            G(2,2) = x(1)
            G = G / r
            x(1) = r
            x(2) = 0.0

        else

            G(:,:) = 0.0
            G(1,1) = 1.0
            G(2,2) = 1.0

        end if

    end

    !>
    !> Wrapper around the planerot 
    subroutine giv(G, a, colrow, mode)

        real(dp), intent(out)   :: G(:,:)
        real(dp), intent(out)   :: a(:)
        integer, intent(in)     :: colrow   ! 1 for column and 2 for row 
        integer, intent(in)     :: mode

        real(dp)                :: buff(2)

        !real(dp)                :: zerotol      = 1.0e-14_dp
        !real(dp)                :: r

        call planerot(G, a)

        ! if ( abs(a(2)) > zerotol ) then

        !     r = norm2(a)
        !     G(1,:) = a(:)
        !     G(2,1) = -a(2)
        !     G(2,2) = a(1)
        !     G = G / r
        !     a(1) = r
        !     a(2) = 0.0_dp

        ! else

        !     G(:,:) = 0.0_dp
        !     G(1,1) = 1.0_dp
        !     G(2,2) = 1.0_dp

        ! end if

        if ( mode==1 ) then

            buff(:) = G(1,:)                
            G(1,:)  = G(2,:)
            G(2,:)  = buff(:)
            buff(1) = a(1)
            a(1)    = a(2)
            a(2)    = buff(1)

        end if

        if ( colrow .ne. 1 ) then

            G(:,:) = transpose(G(:,:))

        end if 

    end

    !>
    !> C interface
    !> Wrapper around the planerot 
    subroutine giv_c(G, a, colrow, mode)

        real(c_double), intent(out)     :: G(:,:)
        real(c_double), intent(out)     :: a(:)
        integer, intent(in)             :: colrow   ! 1 for column and 2 for row 
        integer, intent(in)             :: mode

        real(c_double)                  :: buff(2)

        !real(dp)                :: zerotol      = 1.0e-14_dp
        !real(dp)                :: r

        call planerot(G, a)

        ! if ( abs(a(2)) > zerotol ) then

        !     r = norm2(a)
        !     G(1,:) = a(:)
        !     G(2,1) = -a(2)
        !     G(2,2) = a(1)
        !     G = G / r
        !     a(1) = r
        !     a(2) = 0.0_dp

        ! else

        !     G(:,:) = 0.0_dp
        !     G(1,1) = 1.0_dp
        !     G(2,2) = 1.0_dp

        ! end if

        if ( mode==1 ) then

            buff(:) = G(1,:)                
            G(1,:)  = G(2,:)
            G(2,:)  = buff(:)
            buff(1) = a(1)
            a(1)    = a(2)
            a(2)    = buff(1)

        end if

        if ( colrow .ne. 1 ) then

            G(:,:) = transpose(G(:,:))

        end if 

    end

end module