!
! November 12, 2014
!
! This module defines the following functions and procedures:
!  qpcc_global
!  qpcc_local
!  qpcc_crlarge
!  qpcc_crsmall
!  qpcc_qrupdaddc12
!  qpcc_qcrtp
!  qpcc_start0w
!  qpcc_startxw
!  qpcc_mas
!  qpcc_xdota1d
!  qpcc_kktns
!  qpcc_xdota2d
!  qpcc_lls (uses LAPACK subroutine "dgelsd")
!  qpcc_lslu (uses LAPACK subroutine "dgetrs")
!  qpcc_qrupdremc3
!  qpcc_qrupdaddc3
!  qpcc_as
!  rseed (helps to generate random numbers; in Unix-like OS uses a
!   special file "/dev/random", in other OS uses "date and time"
!   method)
!  print_matrix
!
module qpcc

	integer, parameter :: string_length = 128
	double precision, parameter :: zero = 1.0d-5
	double precision, parameter :: infinity = 1.0d+100
	integer, parameter :: max_iter = 1000

contains

subroutine qpcc_global(m, n, delta, image, iter, info)
	implicit none

	double precision :: delta(n, n), image(m, n), b(m * n**2),&
		q(m * n**2, m * n**2), r(m * n**2, m * n**2),&
		q0(m * n**2, m * n**2), r0(m * n**2, m * n**2),&
		q1(m * n**2, m * n**2), r1(m * n**2, m * n**2), x(m * n**2),&
		obj, obj_tmp
	integer :: m, n, iter, info, step, nn, mm, mm_eq, step2, mn,&
		c2(4, m * n * (n - 1) / 2), i, j, k, l, prob_size, level, w0_t,&
		w1_t, w(m + 3 * m * n * (n - 1) / 2), indx,&
		w0(m + 3 * m * n * (n - 1) / 2),&
		w1(m + 3 * m * n * (n - 1) / 2), pos(m, n, n),&
		ind(3, m * n * (n - 1) / 2), s, const, v(m), from, to
	logical :: prob(m * n * (n - 1) / 2 + 1 - m, m * n * (n - 1)),&
		p(m * n * (n - 1)), feasible, next(2), p_s(m * n * (n - 1))

	info = 0
	iter = 0
	image = 0.0d+0
	const = n * (n - 1) / 2
	step = m * const
	nn = m * n**2
	mm = m + 3 * step
	mm_eq = m + step
	step2 = step * 2
	mn = m * n
	obj = infinity
	l = 0
	do i = 1, (n - 1)
		do j = (i + 1), n
			do k = 1, m
				l = l + 1
				pos(k, i, j) = l
				ind(:, l) = (/k, i, j/)
			end do
		end do
	end do
	call qpcc_crsmall(m, n, step, delta, nn, mm, b, c2, info)
	call qpcc_qrupdaddc12(m, n, step, c2, nn, q, r, info)
	w = (/(i, i = 1, mm)/)
	prob_size = 1
	prob(prob_size, :) = .false.
	prob(prob_size, (/(2 * i - 1, i = 1, m)/)) = .true.
	do while ((prob_size > 0) .and. (info == 0))
		p = prob(prob_size, :)
		prob_size = prob_size - 1
		level = count(p)
		next = .true.
		do s = 1, 2
			p_s = p
			p_s(2 * level + s) = .true.
			if (next(s)) then
				k = ind(1, level + 1)
				j = ind(2, level + 1)
				l = ind(3, level + 1)
				i = 1
				do while ((i < j) .and. next(s))
					if (p_s(2 * (level + 1) - 1)) then
						next(s) = .not. (p_s(2 * pos(k, i, j) - 1) .and.&
							p_s(2 * pos(k, i, l)))
					else
						next(s) = .not. (p_s(2 * pos(k, i, j)) .and.&
							p_s(2 * pos(k, i, l) - 1))
					end if
					i = i + 1
				end do
			end if
			if (next(s)) then
				do k = 1, m
					l = 0
					do i = 1, const
						j = 0
						if (p_s(2 * (k + (i - 1) * m) - 1)) j = 1
						l = l + j * 2**(const - i)
					end do
					v(k) = l
				end do
				i = 1
				do while ((i < m) .and. next(s))
					j = i + 1
					do while ((j < m + 1) .and. next(s))
						if (v(i) > v(j)) next(s) = .false.
						j = j + 1
					end do
					i = i + 1
				end do
			end if
		end do
		from = 1
		to = 2
		if (.not. (next(1) .and. next(2))) then
			to = 0
			if (next(1)) to = 1
			if (next(2)) to = 2
			from = to + 1
			if (to > 0) then
				if (level + 1 == step) then
					from = to
				else
					prob_size = prob_size + 1
					prob(prob_size, :) = p
					prob(prob_size, 2 * level + to) = .true.
				end if
			end if
		end if
		if ((from < to) .or. (from == to)) then
			w0_t = mm_eq + level
			w0 = (/w(1:mm_eq), pack(w((mm_eq + 1):mm), p),&
				pack(w((mm_eq + 1):mm), .not. p)/)
			q0 = q
			r0 = r
			do i = (mm_eq + 1), w0_t
				call qpcc_qrupdaddc3(nn, w0(i) - mm_eq + mn, i,	q0, r0, info)
			end do
		end if
		do l = from, to
			w1_t = w0_t + 1
			w1 = w0
			indx = mm_eq + 2 * level + l
			w1(w1_t) = w0(indx)
			w1(indx) = w0(w1_t)
			q1 = q0
			r1 = r0
			call qpcc_qrupdaddc3(nn, w1(w1_t) - mm_eq + mn,	w1_t, q1,&
				r1, info)
			x = 0.0d+0
			call qpcc_as(m, n, step, nn, mm, mm_eq, b, w1, w1_t, q1, r1,&
				x, info)
			iter = iter + 1
			obj_tmp = 5.0d-1 * dot_product(qpcc_xdota1d(m, n, nn, x),&
				x) + dot_product(b, x)
			if (obj_tmp < obj) then
				feasible = .true.
				i = 1
				do while ((i < step + 1) .and. feasible)
					if (abs(x(mn + 2 * i - 1) * x(mn + 2 * i)) > zero)&
						feasible = .false.
					i = i + 1
				end do
				if (feasible) then
					image = reshape(x(1:mn), (/m, n/))
					obj = obj_tmp
				else
					prob_size = prob_size + 1
					prob(prob_size, :) = p
					prob(prob_size, 2 * level + l) = .true.
				end if
			end if
		end do
	end do
end subroutine qpcc_global

subroutine qpcc_local(m, n, delta, nran, image, info)
	implicit none
	double precision :: delta(n, n), image(m, n),&
		b(m * n**2), x(m * n**2), obj, obj_tmp,&
		q(m * n**2, m * n**2), q_tmp(m * n**2, m * n**2),&
		r(m * n**2, m * n**2), r_tmp(m * n**2, m * n**2), t0, t1
	integer :: m, n, nran, info, desc, c2(4, m * n * (n - 1) / 2), nn,&
		step, mm, mm_eq, mn, w(m + 3 * m * n * (n - 1) / 2), w_t, i, j, k

	step = m * n * (n - 1) / 2
	nn = m * n**2
	mm = m + 3 * step
	mm_eq = m + step
	mn = m * n
	obj = infinity
	call qpcc_crsmall(m, n, step, delta, nn, mm, b, c2, info)
	call qpcc_qrupdaddc12(m, n, step, c2, nn, q, r, info)
	do desc = 1, nran
		if (mod(desc, 2) == 1) then
			call qpcc_startxw(m, n, step, nn, mm, mm_eq, x, w, w_t, info)
		else
			call qpcc_start0w(step, nn, mm, mm_eq, x, w, w_t, info)
		end if
		q_tmp = q
		r_tmp = r
		do i = (mm_eq + 1), (mm_eq + step)
			call qpcc_qrupdaddc3(nn, w(i) - mm_eq + mn, i, q_tmp, r_tmp,&
				info)
		end do
		call qpcc_mas(m, n, step, nn, mm, mm_eq, b, w, w_t, q_tmp,&
			r_tmp, x, info)
		obj_tmp = 5.0d-1 * dot_product(qpcc_xdota1d(m, n, nn, x),&
			x) + dot_product(b, x)
		if (obj_tmp < obj) then
			image = reshape(x(1:(m * n)), (/m, n/))
			obj = obj_tmp
		end if
	end do
end subroutine qpcc_local

subroutine qpcc_crlarge(m, n, step, delta, nn, mm, a, b, c, info)
	implicit none
	double precision :: delta(n, n), a(nn, nn), b(nn), c(nn, mm)
	integer :: m, n, step, nn, mm, info, m2, mn, i, j, row_start,&
		row_end, const, ind, c2(4, step)

	a = 0.0d+0
	m2 = m * 2
	mn = m * n
	const = n * (n - 1) / 2
	do i = 1, const
		row_start = mn + 1 + (i - 1) * m2;
		row_end = row_start + m2 - 1;
		a(row_start:row_end, row_start:row_end) = 1.0d+0;
	end do
	call qpcc_crsmall(m, n, step, delta, nn, mm, b, c2, info)
	c = 0.0d+0
	do i = 1, mm
		if (i < m + 1) then
			do j = 1, n
				ind = i + (j - 1) * m
				c(ind, i) = 1.0d+0
			end do
		else if ((m < i) .and. (i < m + step + 1)) then
			c(c2(:, i - m), i) = (/1.0d+0, -1.0d+0, -1.0d+0, 1.0d+0/)
		else
			ind = i - m - step + mn
			c(ind, i) = 1.0d+0
		end if
	end do
end subroutine qpcc_crlarge

subroutine qpcc_crsmall(m, n, step, delta, nn, mm, b, c2, info)
	implicit none
	double precision :: delta(n, n), a(nn, nn), b(nn)
	integer :: m, n, step, nn, mm, c2(4, step), info, mn, m2, p, q, i,&
		j, k, from, to

	info = 0
	mn = m * n
	m2 = m * 2
	p = 0
	q = 0
	b(1:mn) = 0.0d+0
	do i = 1, (n - 1)
		do j = (i + 1), n
			p = p + 1
			b((mn + (p - 1) * m2 + 1):(mn + p * m2)) = -delta(i, j)
			do k = 1, m
				q = q + 1
				c2(:, q) = (/(i - 1) * m + k, (j - 1) * m + k,&
					(mn + 2 * q - 1), (mn + 2 * q)/)
			end do
		end do
	end do
end subroutine qpcc_crsmall

subroutine qpcc_qrupdaddc12(mm, nn, step, c2, n, q, r, info)
	implicit none
	double precision :: q(n, n), r(n, n), z(2), p(2, 2)
	integer :: mm, nn, step, c2(4, step), n, k, info, indices(nn), i,&
		k_mm

	info = 0
	indices = (/(mm * (i - 1), i = 1, nn)/)
	do i = 1, n
		q(1:(i - 1), i) = 0.0d+0
		q(i, i) = 1.0d+0
		q((i + 1):n, i) = 0.0d+0
		r(:, i) = 0.0d+0
	end do
	do k = 1, mm
		do i = n, (k + 1), -1
			z = sum(q(k + indices, (i - 1):i), dim = 1)
			if (abs(z(2)) > 0.0d+0) then
				call qpcc_qcrtp(z, p, info)
				q(:, (i - 1):i) = matmul(q(:, (i - 1):i), p)
			end if
		end do
		r(1:k, k) = sum(q(k + indices, 1:k), dim = 1)
	end do
	do k = (mm + 1), (mm + step)
		k_mm = k - mm
		do i = n, (k + 1), -1
			z(1) = q(c2(1, k_mm), i - 1) - q(c2(2, k_mm), i - 1)-&
				q(c2(3, k_mm), i - 1) + q(c2(4, k_mm), i - 1)
			z(2) = q(c2(1, k_mm), i) - q(c2(2, k_mm), i) -&
				q(c2(3, k_mm), i) + q(c2(4, k_mm), i)
			if (abs(z(2)) > 0.0d+0) then
				call qpcc_qcrtp(z, p, info)
				q(:, (i - 1):i) = matmul(q(:, (i - 1):i), p)
			end if
		end do
		do i = 1, k
			r(i, k) = q(c2(1, k_mm), i) - q(c2(2, k_mm), i) -&
				q(c2(3, k_mm), i) + q(c2(4, k_mm), i)
		end do
	end do
end subroutine qpcc_qrupdaddc12

subroutine qpcc_qcrtp(z, p, info)
	implicit none
	double precision :: z(2), p(2, 2), ro, c, s, dnrm2
	integer :: info

	info = 0
	ro = dnrm2(2, z, 1)
	if (z(1) < 0.0d+0) ro = -ro
	c = z(1) / ro
	s = z(2) / ro
	p(1, 1) = c
	p(1, 2) = s
	p(2, 1) = s
	p(2, 2) = -c
end subroutine qpcc_qcrtp

subroutine qpcc_start0w(step, nn, mm, mm_eq, x, w, w_t, info)
	implicit none
	double precision :: x(nn), rn(step)
	integer :: step, nn, mm, mm_eq, w(mm), w_t, info, i, ind, tmp, left

	info = 0
	x = 0.0d+0
	w = (/(i, i = 1, mm)/)
	w_t = mm_eq + step
	call rseed()
	call random_number(rn)
	do i = 1, step
		left = 0
		if (rn(i) > 0.5d+0) then
			left = 1
		end if
		ind = mm_eq + 2 * i - left
		tmp = w(ind)
		w(ind) = w(mm_eq + i)
		w(mm_eq + i) = tmp
	end do
end subroutine qpcc_start0w

subroutine qpcc_startxw(m, n, step, nn, mm, mm_eq, x, w, w_t, info)
	implicit none
	double precision :: x(nn), val
	integer, parameter :: from = -10, to = 10
	integer :: m, n, step, nn, mm, mm_eq, w(mm), w_t, info, mn,&
		set(n), i, j, k, p, q, tmp, left

	info = 0
	x = 0.0d+0
	w = (/(i, i = 1, mm)/)
	w_t = mm_eq + step
	mn = m * n
	set = (/(m * (i - 1), i = 1, n)/)
	call rseed()
	call random_number(x)
	x = from + (to - from) * x
	do k = 1, m
		x(k + set) = x(k + set) - sum(x(k + set)) / n
	end do
	p = 1
	q = 1
	do i = 1, (n - 1)
		do j = (i + 1), n
			do k = 1, m
				val = x(k + set(i)) - x(k + set(j))
				left = 0
				if (val > 0.0d+0) left = 1
				x(mn + p + 1 - left) = abs(val)
				x(mn + p + left) = 0.0d+0
				tmp = w(mm_eq + p + left)
				w(mm_eq + p + left) = w(mm_eq + q)
				w(mm_eq + q) = tmp
				p = p + 2
				q = q + 1
			end do
		end do
	end do
end subroutine qpcc_startxw

subroutine qpcc_mas(mm, nn, step, n, m, m_eq, b, w, w_t, q, r, x,&
	info)
	implicit none
	double precision :: b(n), x(n), p(n), lambda(m), val1, val2,&
		alpha, q(n, n), r(n, n)
	integer :: mm, nn, step, n, m, m_eq, w(m), w_t, info, i, j, k, l,&
		tmp, mmnn, iter
	logical :: next, remove, odd

	info = 0
	iter = 0
	mmnn = mm * nn
	next = .false.
	if (info == 0) next = .true.
	do while (next)
		iter = iter + 1
		if (iter > max_iter)then
			print *, 'iter > ', max_iter
			return
		end if
		p = qpcc_xdota1d(mm, nn, n, x) + b
		call qpcc_kktns(mm, nn, step, n, m, w, w_t, q, r, p, lambda,&
			info)
		if (info == 0) then
			if (all(abs(p) < zero)) then
				remove = .false.
				do while (next .and. (.not. remove))
					j = 0
					val1 = -zero
					do i = (m_eq + 1), w_t
						if ((w(i) > m_eq) .and.&
							(lambda(w(i)) < val1)) then
							val1 = lambda(w(i))
							j = i
						end if
					end do
					if (j > 0) then
						i = m_eq
						odd = mod(w(j) - m_eq, 2) == 1
						do while ((i < w_t) .and. (.not. remove))
							i = i + 1
							if ((i /= j) .and.&
								(odd .and. (w(i) == w(j) + 1)) .or.&
								((.not. odd) .and. (w(i) == w(j) - 1)))&
								remove = .true.
						end do
						if (remove) then
							if (j < w_t) call qpcc_qrupdremc3(n, w_t,&
								w((j + 1):w_t) - m_eq + mmnn, j, q, r, info)
							tmp = w(j)
							do i = j, (w_t - 1)
								w(i) = w(i + 1)
							end do
							w(w_t) = tmp
							w_t = w_t - 1
						else
							lambda(w(j)) = 0.0d+0
						end if
					else
						next = .false.
					end if
				end do
			else
				j = 0
				alpha = 1.0d+0
				do i = (w_t + 1), m
					val1 = p(w(i) - m_eq + mmnn)
					if (val1 < -zero) then
						val2 = (-x(w(i) - m_eq + mmnn)) / val1
						if (val2 - alpha < zero) then
							alpha = val2
							j = i
						end if
					end if
				end do
				x = x + alpha * p
				if ((j > 0) .and. (w_t < n)) then
					w_t = w_t + 1
					tmp = w(j)
					w(j) = w(w_t)
					w(w_t) = tmp
					if (w_t < (n + 1)) call qpcc_qrupdaddc3(n, w(w_t) -&
						m_eq + mmnn, w_t, q, r, info)
				end if
			end if
		else
			next = .false.
		end if
	end do
end subroutine qpcc_mas

function qpcc_xdota1d(mm, nn, n, x)
	implicit none
	double precision :: x(n), qpcc_xdota1d(n)
	integer :: mm, nn, n, mmnn, mm2, finish, j, from, to

	mmnn = mm * nn
	mm2 = mm * 2
	finish = nn * (nn - 1) / 2
	qpcc_xdota1d(1:mmnn) = 0.0d+0
	do j = 1, finish
		from = mmnn + 1 + (j - 1) * mm2
		to = mmnn + j * mm2
		qpcc_xdota1d(from:to) = sum(x(from:to))
	end do
	return
end function qpcc_xdota1d

subroutine qpcc_kktns(mm, nn, step, n, m, w, w_t, q, r, p, lambda,&
	info)
	implicit none
	double precision :: q(n, n), r(n, n), p(n), lambda(m),&
		red_hess(n - w_t, n - w_t), red_grad(n - w_t), ls_m(w_t, w_t),&
		ls_rhs(w_t)
	integer :: mm, nn, step, n, m, w(m), w_t, info,&
		n_w_t, i

	info = 0
	n_w_t = n - w_t
	if (w_t < n) then
		red_hess = matmul(qpcc_xdota2d(mm, nn, n, n_w_t,&
			q(:, (w_t + 1):n)), q(:, (w_t + 1):n))
		red_grad = -matmul(p, q(:, (w_t + 1):n))
		call qpcc_lls(n_w_t, red_hess, red_grad, info)
	end if
	if (info == 0) then
		if (((w_t < n) .and. (all(abs(red_grad) < zero))) .or.&
			(w_t == n)) then
			do i = 1, w_t
				ls_m(1:i, i) = r(1:i, i)
				ls_m((i + 1):w_t, i) = 0.0d+0
			end do
			ls_rhs = matmul(p, q(:, 1:w_t))
			call qpcc_lslu(w_t, ls_m, ls_rhs, info)
			if (info == 0) then
				p = 0.0d+0
				lambda = 0.0d+0
				lambda(w(1:w_t)) = ls_rhs
			end if
		else
			p = matmul(q(:, (w_t + 1):n), red_grad)
		end if
	end if
end subroutine qpcc_kktns

function qpcc_xdota2d(mm, nn, n, n_w_t, x)
	implicit none
	double precision :: x(n, n_w_t), qpcc_xdota2d(n_w_t, n)
	integer :: mm, nn, n, n_w_t, mmnn, mm2, finish, i, j, from, to

	mmnn = mm * nn
	mm2 = mm * 2
	finish = nn * (nn - 1) / 2
	do i = 1, n_w_t
		qpcc_xdota2d(i, 1:mmnn) = 0.0d+0
		do j = 1, finish
			from = mmnn + 1 + (j - 1) * mm2
			to = mmnn + j * mm2
			qpcc_xdota2d(i, from:to) = sum(x(from:to, i))
		end do
	end do
	return
end function qpcc_xdota2d

subroutine qpcc_lls(n, a, b, info)
	implicit none
	double precision, allocatable :: work(:)
	integer, allocatable :: iwork(:)
	double precision :: a(n, n), b(n), s(n), rcond
	integer :: n, info, rank, lwork, liwork

	rcond = -1.0d+0
	allocate(work(1), iwork(1))
	lwork = -1
	call dgelsd(n, n, 1, a, n, b, n, s, rcond, rank, work, lwork,&
		iwork, info)
	if (info == 0) then
		lwork = work(1)
		liwork = iwork(1)
		deallocate(work, iwork)
		allocate(work(max(1, lwork)), iwork(max(1, liwork)))
		call dgelsd(n, n, 1, a, n, b, n, s, rcond, rank, work, lwork,&
			iwork, info)
	end if
	deallocate(work, iwork)
end subroutine qpcc_lls

subroutine qpcc_lslu(n, a, b, info)
	implicit none
	double precision :: a(n, n), b(n)
	integer :: n, info, i, ipiv(n)

	ipiv = (/(i, i = 1, n)/)
	call dgetrs('n', n, 1, a, n, ipiv, b, n, info)
end subroutine qpcc_lslu

subroutine qpcc_qrupdremc3(n, m, nr, k, q, r, info)
	implicit none
	double precision :: q(n, n), r(n, n), z(2), p(2, 2)
	integer :: n, m, k, nr(m - k), info, i

	info = 0
	do i = (k + 1), m
		z = q(nr(i - k), (i - 1):i)
		if (abs(z(2)) > 0.0d+0) then
			call qpcc_qcrtp(z, p, info)
			q(:, (i - 1):i) = matmul(q(:, (i - 1):i), p)
			r((i - 1):i, i:m) = matmul(p, r((i - 1):i, i:m))
		end if
		r(1:(i - 1), i - 1) = r(1:(i - 1), i)
	end do
	r(:, m) = 0.0d+0
end subroutine qpcc_qrupdremc3

subroutine qpcc_qrupdaddc3(n, nr, k, q, r, info)
	implicit none
	double precision :: q(n, n), r(n, n), z(2), p(2, 2)
	integer :: n, nr, k, info, i

	info = 0
	do i = n, (k + 1), -1
		z = q(nr, (i - 1):i)
		if (abs(z(2)) > 0.0d+0) then
			call qpcc_qcrtp(z, p, info)
			q(:, (i - 1):i) = matmul(q(:, (i - 1):i), p)
		end if
	end do
	r(1:k, k) = q(nr, 1:k)
end subroutine qpcc_qrupdaddc3

subroutine qpcc_as(mm, nn, step, n, m, m_eq, b, w, w_eq, q, r, x,&
	info)
	implicit none
	double precision :: b(n), x(n), p(n), lambda(m), val1, val2,&
		alpha, q(n, n), r(n, n)
	integer :: mm, nn, step, n, m, m_eq, w(m), w_eq, info, w_t, i, j,&
		k, l, tmp, mmnn, iter
	logical :: next

	info = 0
	iter = 0
	w_t = w_eq
	mmnn = mm * nn
	next = .false.
	if (info == 0) next = .true.
	do while (next)
		iter = iter + 1
		if (iter > max_iter) then
			return
		end if
		p = qpcc_xdota1d(mm, nn, n, x) + b
		call qpcc_kktns(mm, nn, step, n, m, w, w_t, q, r, p, lambda,&
			info)
		if (info == 0) then
			j = 0
			if (all(abs(p) < zero)) then
				val1 = -zero
				do i = (w_eq + 1), w_t
					if (lambda(w(i)) < val1) then
						val1 = lambda(w(i))
						j = i
					end if
				end do
				if (j > 0) then
					if (j < w_t) call qpcc_qrupdremc3(n, w_t,&
						w((j + 1):w_t) - m_eq + mmnn, j, q, r, info)
					tmp = w(j)
					do i = j, (w_t - 1)
						w(i) = w(i + 1)
					end do
					w(w_t) = tmp
					w_t = w_t - 1
				else
					next = .false.
				end if
			else
				alpha = 1.0d+0
				do i = (w_t + 1), m
					val1 = p(w(i) - m_eq + mmnn)
					if (val1 < -zero) then
						val2 = (-x(w(i) - m_eq + mmnn)) / val1
						if (val2 - alpha < zero) then
							alpha = val2
							j = i
						end if
					end if
				end do
				x = x + alpha * p
				if ((j > 0) .and. (w_t < n)) then
					w_t = w_t + 1
					tmp = w(j)
					w(j) = w(w_t)
					w(w_t) = tmp
					if (w_t < (n + 1)) call qpcc_qrupdaddc3(n, w(w_t) -&
						m_eq + mmnn, w_t, q, r, info)
				end if
			end if
		else
			next = .false.
		end if
	end do
end subroutine qpcc_as

subroutine rseed()
	implicit none
	integer, allocatable :: seed(:)
	integer :: n, stat, dt(8), i
		       
	call random_seed(size = n)
	allocate(seed(n))
	open(1, file = '/dev/urandom', access = 'stream',&
		form = 'unformatted', iostat = stat)
	if (stat == 0) then
		read(1) seed
		close(1)
	else
		call date_and_time(values = dt)
		seed = (dt(8) + dt(7) * 1000 + dt(6) * 1000 * 60) +&
			37 * (/(i, i = 0, (n - 1))/)
	end if
	call random_seed(put = seed)
	deallocate(seed)
end subroutine rseed

subroutine print_matrix(n, m, a, a_name)
	implicit none

	character(len = *) :: a_name
	character(len = string_length) :: output_format
	double precision :: a(n, m)
	integer :: n, m, i, j

	output_format = '(f10.4)'
	write(*, *) trim(a_name)
	do i = 1, n
		do j = 1, m
			write(*, output_format, advance = 'no') a(i, j)
		end do
		write(*, *)
	end do
	write(*, *)
end subroutine print_matrix

end module qpcc
