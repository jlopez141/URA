program analize_constr

    implicit none
    real*8 :: max_sigma_angle, max_sigma_dist, mat_mean_angle(2,2), mat_mean_dist(2,2)


    max_sigma_dist  = 3.0
    max_sigma_angle = 10.0

    mat_mean_angle = reshape( (/   0.0, 152.3, &
                                 108.0, 0.0   /), (/2,2/))
    mat_mean_dist  = reshape( (/   0.0, 1.6, &
                                   1.6, 0.0   /), (/2,2/))

    call compute_number_constraints(max_sigma_angle, max_sigma_dist, mat_mean_angle, mat_mean_dist)

contains

subroutine compute_number_constraints(max_sigma_angle, max_sigma_dist, mat_mean_angle, mat_mean_dist)
    implicit none
    real*8, intent(in) :: max_sigma_angle, max_sigma_dist, mat_mean_angle(:,:), mat_mean_dist(:,:)
    integer              :: iat, iesp, jesp, natoms, N_species, max_neigh, &
                            max_pair, N_neigh, N_pair, d1, d2, i, kkk
    integer, allocatable :: mat_neighbor(:,:), atomtype(:), numions(:), &
                            N_const_angle(:,:), N_const_dist(:,:)
    real*8, allocatable  :: sigma_angle(:), sigma_dist(:), mean_angle(:), mean_dist(:)

    open(unit = 1, action = "read", status = "old", file = "deviations.log")

    read(1,*) natoms, N_species

    allocate(mat_neighbor(N_species, N_species), atomtype(natoms), numions(N_species), &
             N_const_angle(N_species, N_species), N_const_dist(N_species, N_species))
    read(1,*) atomtype(:)
    do iesp = 1, N_species
        read(1,*) mat_neighbor(iesp,:)
    enddo
    read(1,*)

    max_neigh = maxval(mat_neighbor)
    max_pair  = max_neigh*(max_neigh-1)/2
    allocate(sigma_angle(max_pair), sigma_dist(max_neigh), &
             mean_angle(max_pair), mean_dist(max_neigh))

    numions = 0
    N_const_angle = 0
    N_const_dist  = 0
    do iat = 1, natoms
        numions(atomtype(iat)) = numions(atomtype(iat)) + 1
        do iesp = 1, N_species

            if ( atomtype(iat) == iesp ) cycle

            N_neigh = mat_neighbor(atomtype(iat),iesp)
            N_pair = N_neigh*(N_neigh-1)/2

            read(1,*) d1, d2, mean_angle(:)
            read(1,*) d1, d2, sigma_angle(:)
            read(1,*) d1, d2, mean_dist(:)
            read(1,*) d1, d2, sigma_dist(:)
            read(1,*) 

            kkk = 0
            do i = 1, N_pair
                !if ( abs(mean_angle(i) - mat_mean_angle(atomtype(iat),iesp) ) < 10 ) then
                    if ( sigma_angle(i) < max_sigma_angle) then
                    KKK = kkk + 1
                        N_const_angle(atomtype(iat),iesp) = N_const_angle(atomtype(iat),iesp) + 1
                    !else
                        !if (atomtype(iat) == 2) cycle
                        !print*, "angle", iat, iesp, mean_angle(i), sigma_angle(i)
                    endif
                !endif
            enddo
            !print*, atomtype(iat), iesp, kkk
            do i = 1, N_neigh
                !if ( abs(mean_angle(i) - mat_mean_angle(atomtype(iat),iesp) ) < 0.6 ) then
                    if (sigma_dist(i) < max_sigma_dist) then
                        N_const_dist(atomtype(iat),iesp) = N_const_dist(atomtype(iat),iesp) + 1
                    !else
                       !print*, "dist", iat, iesp, sigma_dist(i)
                    endif
                !endif
            enddo

            print*, iat, N_neigh, mean_angle(:N_pair)

        enddo
    enddo

    do iesp = 1, N_species
        do jesp = 1, N_species
            if (iesp == jesp) cycle
            N_neigh = mat_neighbor(iesp,jesp)
            N_pair = N_neigh*(N_neigh-1)/2

            print*, iesp, jesp, N_const_angle(iesp,jesp)/numions(iesp), N_const_dist(iesp,jesp)/numions(iesp)


            ! print*, N_const_angle(iesp,jesp), real(numions(iesp),8), real(N_pair,8)
            ! print*, N_const_angle(iesp,jesp) / real(numions(iesp),8)! / real(N_pair,8)

            ! print*, N_const_dist(iesp,jesp), real(numions(iesp),8), real(N_pair,8)
            ! print*, N_const_angle(iesp,jesp) / real(numions(iesp),8)! / real(N_pair,8)
            ! print*, ""
        enddo
    enddo

    !print*, N_const_angle!, numions

    close(unit = 1)

end subroutine compute_number_constraints


end program analize_constr