module gaussian_filter
  use params
  implicit none
  real(8),parameter :: GF_SIGMA = 0.5d0  ! GF_SIGMA > 0.d0
  !real(8),parameter :: GF_SIGMA = 1.00d0  ! GF_SIGMA > 0.d0
  real(8),allocatable :: gaussian_filter_coef_1D(:)
  real(8),allocatable :: gaussian_filter_coef_2D(:,:)

  contains

    subroutine set_gaussian_filter_1D
      integer            :: i
      real(8)            :: r2
      real(8)            :: weight
      real(8)            :: sum_weight
      ! ------------------------------------------------------------------------
      if( .NOT. allocated( gaussian_filter_coef_1D ) ) then
        allocate( gaussian_filter_coef_1D(-ng:ng) )
        sum_weight = 0.d0
        do i = -ng, ng
          r2 = dble(i*i)
          weight = exp( -1.d0 * r2/(2.d0*GF_SIGMA*GF_SIGMA) )
          gaussian_filter_coef_1D(i) = weight
          sum_weight = sum_weight + weight
        end do
        gaussian_filter_coef_1D(:) = gaussian_filter_coef_1D(:) / sum_weight
      end if
    end subroutine set_gaussian_filter_1D

    subroutine gaussian_filtering_1D( var_in, var_out )
      real(8),intent(in)    :: var_in(-ng:ng)
      real(8),intent(inout) :: var_out
      ! +++++ Local variables +++++
      integer               :: i, j, k
      real(8)               :: var_filtered
      ! ------------------------------------------------------------------------
      var_filtered = 0.d0
      do i = -ng, ng
        var_filtered = var_filtered + gaussian_filter_coef_1D(i) * var_in(i)
      end do
      var_out = var_filtered
    end subroutine gaussian_filtering_1D

  subroutine set_gaussian_filter_2D
      integer            :: i, j
      real(8)            :: r2
      real(8)            :: weight
      real(8)            :: sum_weight
      ! ------------------------------------------------------------------------
      if( .NOT. allocated( gaussian_filter_coef_2D ) ) then
        allocate( gaussian_filter_coef_2D(-ng:ng,-ng:ng) )
        sum_weight = 0.d0
        do j = -1, 1
          do i = -1, 1
        !do j = -ng, ng
        !  do i = -ng, ng
            r2 = dble(i*i+j*j)
            weight = exp( -1.d0 * r2/(2.d0*GF_SIGMA*GF_SIGMA) )
            gaussian_filter_coef_2D(i,j) = weight
            sum_weight = sum_weight + weight
          end do
        enddo
        gaussian_filter_coef_2D(:,:) = gaussian_filter_coef_2D(:,:) / sum_weight
      end if
    end subroutine set_gaussian_filter_2D

    subroutine gaussian_filtering_2D( var_in, var_out )
      real(8),intent(in)    :: var_in(-ng:ng,-ng:ng)
      real(8),intent(inout) :: var_out
      ! +++++ Local variables +++++
      integer               :: i, j, k
      real(8)               :: var_filtered
      ! ------------------------------------------------------------------------
      var_filtered = 0.d0
      do j = -1, 1
        do i = -1, 1
      !do j = -ng, ng
      !  do i = -ng, ng
          var_filtered = var_filtered + gaussian_filter_coef_2D(i,j) * var_in(i,j)
        end do
      enddo
      var_out = var_filtered
    end subroutine gaussian_filtering_2D

end module gaussian_filter