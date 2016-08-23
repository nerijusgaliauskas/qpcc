!
! October 19, 2015
!
! This driver demonstrates the usage of procedure "qpcc_global". The
! procedure is defined in module "qpcc" (qpcc_mod.f90).
!
program qpcc_pglobal_dr
	use mpi
	use qpcc
	implicit none
	character(len = string_length) :: data_name, string
	double precision, allocatable :: delta(:, :), image(:, :), image_tmp(:, :)
	double precision :: re, re_tmp
	integer, parameter :: tag_m = 1001, tag_n = 1002, tag_delta = 1003,&
		tag_image = 1004
	integer :: m, n, i, j, k, iter, info, values0(8), values1(8), ierr, size,&
		rank, status(mpi_status_size), level, level2, rest, num, lev
	logical, allocatable :: snode(:)

	call mpi_init(ierr)
	call mpi_comm_size(mpi_comm_world, size, ierr)
	call mpi_comm_rank(mpi_comm_world, rank, ierr)
	num = int(size / 2.0d+0)
	level = 0
	do while (num > 0)
		num = int(num / 2.0d+0)
		level = level + 1
	end do
	if (rank == 0) then
		read(*, *) data_name
		read(*, *) m, n
		allocate(delta(n, n), image(m, n), image_tmp(m, n),&
			snode(m * n * (n - 1)))
		read(*, *) ((delta(i, j), j = 1, n), i = 1, n)
		delta = sqrt((n * (n - 1) / 2) / sum((/((delta(i, j)**2,&
			j = (i + 1), n), i = 1, (n - 1))/))) * delta
		do i = 1, size - 1
			call mpi_send(m, 1, mpi_integer, i, tag_m, mpi_comm_world, ierr)
			call mpi_send(n, 1, mpi_integer, i, tag_n, mpi_comm_world, ierr)
			call mpi_send(delta, n * n, mpi_double_precision, i, tag_delta,&
				mpi_comm_world, ierr)
		end do
		snode = .false.
		snode((/(2 * i - 1, i = 1, m + level)/)) = .true.
		call date_and_time(VALUES = values0)
		call qpcc_global(m, n, delta, image, snode, iter, info)
		re = sqrt( sum( (/( ((&
			sum( (/( abs(image(k, i) - image(k, j)),&
			k = 1, m )/) ) - delta(i, j))**2, j = (i + 1), n),&
			i = 1, (n - 1) )/) ) / sum( (/( (delta(i, j)**2,&
			j = (i + 1), n), i = 1, (n - 1) )/) ) )
		do i = 1, size - 1
			call mpi_recv(image_tmp, m * n, mpi_double_precision, i, tag_image,&
				mpi_comm_world, status, ierr)
			re_tmp = sqrt( sum( (/( ((&
				sum( (/( abs(image_tmp(k, i) - image_tmp(k, j)),&
				k = 1, m )/) ) - delta(i, j))**2, j = (i + 1), n),&
				i = 1, (n - 1) )/) ) / sum( (/( (delta(i, j)**2,&
				j = (i + 1), n), i = 1, (n - 1) )/) ) )
			if (re_tmp < re) then
				image = image_tmp
				re = re_tmp
			end if
		end do
		call date_and_time(VALUES = values1)
		write(string, '("image (", I0, " x ", I0, ")")') m, n
		call print_matrix(m, n, image, string)
		call print_matrix(1, 1, (/re/), 'relative error')
		write(*, '(" calculations started")')
		write(*, '("    ", I0, "/", I0, "/", I0, ", ", I0, ":", I0, ":", I0, ".",&
			I0)') values0(2), values0(3), values0(1), values0(5), values0(6),&
			values0(7), values0(8)
		write(*, *)
		write(*, '(" calculations finished")')
		write(*, '("    ", I0, "/", I0, "/", I0, ", ", I0, ":", I0, ":", I0, ".",&
			I0)') values1(2), values1(3), values1(1), values1(5), values1(6),&
			values1(7), values1(8)
		deallocate(delta, image, image_tmp, snode)
	else
		call mpi_recv(m, 1, mpi_integer, 0, tag_m, mpi_comm_world, status, ierr)
		call mpi_recv(n, 1, mpi_integer, 0, tag_n, mpi_comm_world, status, ierr)
		allocate(delta(n, n), image(m, n), snode(m * n * (n - 1)))
		call mpi_recv(delta, n * n, mpi_double_precision, 0, tag_delta,&
			mpi_comm_world, status, ierr)
		level2 = 2**level
		rest = size - level2
		if (rank < size - 2 * rest) then
			num = rank
			lev = level
		else
			num = rank + 2 * level2 - size
			lev = level + 1
		end if
		snode = .false.
		snode((/(2 * i - 1, i = 1, m)/)) = .true.
		do i = m + 1, m + lev
			if (mod(num, 2) == 1) then
				snode(2 * i - 1) = .true.
			else
				snode(2 * i) = .true.
			end if
			num = int(num / 2.0d+0)
		end do
		call qpcc_global(m, n, delta, image, snode, iter, info)
		call mpi_send(image, m * n, mpi_double_precision, 0, tag_image,&
			mpi_comm_world, ierr)
		deallocate(delta, image, snode)
	end if
	call mpi_finalize(ierr)
end program qpcc_pglobal_dr
