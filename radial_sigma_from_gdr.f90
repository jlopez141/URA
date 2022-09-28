program kk
implicit none

integer :: bins, i, N_neigh, numions(2)
real*8  :: r, rmax, V_mean
real*8, allocatable :: total(:), partial(:,:), mean(:), sigma(:)

bins    = 400
N_neigh = 6
rmax    = 5.0

numions = (/192, 384/)
V_mean  = 7480.6997178475131

allocate(total(bins), partial(N_neigh,bins), mean(N_neigh), sigma(N_neigh))

open(unit = 1, action="read", status = "old", file = "contr_tot_gdr12.dat")
do i = 1, bins
    read(1,*) r, total(i), partial(:,i)
enddo
close(unit = 1)


call get_mean_sigma_dist(N_neigh,partial, mean, sigma, rmax, numions(2)/V_mean, .true.)

do i = 1, N_neigh
    print*, i, mean(i), sigma(i), sigma(i)/mean(i)*100
enddo


contains

subroutine get_mean_sigma_dist( N_pair, dist_distr, mean, sigma, rmax, norm, normalize )
    implicit none
    integer, intent(in)   :: N_pair
    real*8, intent(in)    :: rmax, norm
    real*8, intent(inout) :: dist_distr(:,:)
    real*8, intent(out)   :: mean(:), sigma(:)
    logical, intent(in)   :: normalize
    real*8, parameter     :: pi = acos(-1.0d0)
    real*8                :: x, dx, aux
    integer               :: ipair, i, N   

    N = size(dist_distr,2)
    dx = rmax/N

    mean = 0.0d0
    sigma = 0.0d0


    do ipair = 1, N_pair

        if (normalize) then
            aux = 0.0d0
            do i = 1, N
                x = (i-0.5)*dx
                aux = aux+  dist_distr(ipair,i) * dx*(4.0*pi*(dx*(i-0.5d0))**2)*norm  !norm = numions(iesp)/V_mean
            enddo
        else
            aux = 1.0d0
        endif

        print*, aux

        do i = 1, N
            x = (i-0.5)*dx
            mean(ipair) = mean(ipair) +  dist_distr(ipair,i) * x/aux * dx *(4.0*pi*(dx*(i-0.5d0))**2)*norm  !norm = numions(iesp)/V_mean
        enddo

        do i = 1, N
            x = (i-0.5)*dx
            sigma(ipair) = sigma(ipair) + ( x/aux-mean(ipair) )**2 * dist_distr(ipair,i) * dx *(4.0*pi*(dx*(i-0.5d0))**2)*norm 
        enddo

    enddo
    sigma = sqrt(sigma)



end subroutine get_mean_sigma_dist

endprogram kk