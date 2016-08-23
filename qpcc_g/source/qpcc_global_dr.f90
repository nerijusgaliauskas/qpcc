!
! December 27, 2014
!
! This driver demonstrates the usage of procedure "qpcc_global". The
! procedure is defined in module "qpcc" (qpcc_mod.f90).
!
program qpcc_global_dr
	use qpcc
	implicit none
	character(len = string_length) :: data_name, string
	double precision, allocatable :: delta(:, :), image(:, :)
	double precision :: re
	integer :: m, n, i, j, k, iter, info, values0(8), values1(8)

	read(*, *) data_name
	read(*, *) m, n
	allocate(delta(n, n), image(m, n))
	read(*, *) ((delta(i, j), j = 1, n), i = 1, n)
	delta = sqrt((n * (n - 1) / 2) / sum((/((delta(i, j)**2,&
		j = (i + 1), n), i = 1, (n - 1))/))) * delta
	call date_and_time(VALUES = values0)
	call qpcc_global(m, n, delta, image, iter, info)
	call date_and_time(VALUES = values1)
	re = sqrt( sum( (/( ((&
		sum( (/( abs(image(k, i) - image(k, j)),&
		k = 1, m )/) ) - delta(i, j))**2, j = (i + 1), n),&
		i = 1, (n - 1) )/) ) / sum( (/( (delta(i, j)**2,&
		j = (i + 1), n), i = 1, (n - 1) )/) ) )
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
	deallocate(delta, image)
end program qpcc_global_dr
