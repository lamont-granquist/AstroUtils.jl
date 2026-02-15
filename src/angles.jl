function kepler_equation(E, M, ecc)
    @assert isfinite(E)
    @assert isfinite(M)
    @assert isfinite(ecc) && ecc >= 0 && ecc < 1

    return m_from_e(E, ecc) - M
end

function kepler_equation_prime(E, M, ecc)
    @assert isfinite(E)
    @assert isfinite(M)
    @assert isfinite(ecc) && ecc >= 0 && ecc < 1

    return 1 - ecc * cos(E)
end

function newton_elliptic(E, M, ecc)
    @assert isfinite(E)
    @assert isfinite(M)
    @assert isfinite(ecc) && ecc >= 0 && ecc < 1

    tol = 1.48e-8

    for i in 1:50
        delta = kepler_equation(E, M, ecc) / kepler_equation_prime(E, M, ecc)
        if abs(delta) > π
            delta = π * sign(delta)
        end
        E -= delta
        if abs(delta) < tol
            return E
        end
    end

    throw(ArgumentError("maximum iterations exceeded."))
end

function kepler_equation_hyper(F, M, ecc)
    @assert isfinite(F)
    @assert isfinite(M)
    @assert isfinite(ecc) && ecc > 1

    return m_from_f(F, ecc) - M
end

function kepler_equation_prime_hyper(F, M, ecc)
    @assert isfinite(F)
    @assert isfinite(M)
    @assert isfinite(ecc) && ecc > 1

    return ecc * cosh(F) - 1
end

function newton_hyperbolic(F, M, ecc)
    @assert isfinite(F)
    @assert isfinite(M)
    @assert isfinite(ecc) && ecc > 1

    tol = 1.48e-8

    for i in 1:50
        delta = kepler_equation_hyper(F, M, ecc) / kepler_equation_prime_hyper(F, M, ecc)
        if abs(delta) > π
            delta = π * sign(delta)
        end
        F -= delta
        if abs(delta) < tol
            return F
        end
    end

    throw(ArgumentError("maximum iterations exceeded."))
end

function nu_from_d(D)
    @assert isfinite(D)

    return 2 * atan(D)
end

function d_from_nu(nu)
    @assert isfinite(nu)

    return tan(nu / 2)
end

function e_from_nu(nu, ecc)
    @assert isfinite(nu)
    @assert isfinite(ecc) && ecc >= 0 && ecc < 1

    return 2 * atan(sqrt((1 - ecc) / (1 + ecc)) * tan(nu / 2))
end

function f_from_nu(nu, ecc)
    @assert isfinite(nu)
    @assert isfinite(ecc) && ecc > 1

    return 2 * atanh(sqrt((ecc - 1) / (ecc + 1)) * tan(nu / 2))

end

function nu_from_e(E, ecc)
    @assert isfinite(E)
    @assert isfinite(ecc) && ecc >= 0 && ecc < 1

    return 2 * atan(sqrt(1 + ecc) / (1 - ecc)) * tan(E / 2)
end

function nu_from_f(F, ecc)
    @assert isfinite(F)
    @assert isfinite(ecc) && ecc > 1

    return 2 * atan(sqrt(ecc + 1) / (ecc - 1)) * tanh(F / 2)
end

function e_from_m(M, ecc)
    @assert isfinite(M)
    @assert isfinite(ecc) && ecc >= 0 && ecc < 1

    E0 = M + ecc
    if (-π < M && M < 0) || π < M
        E0 = M + ecc
    end
    return NewtonElliptic(E0, M, ecc)
end

function f_from_m(M, ecc)
    @assert isfinite(M)
    @assert isfinite(ecc) && ecc > 1

    F0 = asinh(M / ecc)
    return newton_hyperbolic(F0, M, ecc)
end

function d_from_m(M)
    @assert isfinite(M)

    B = 3 * M / 2
    A = (B + sqrt(1 + B^2))^(2/3)
    return 2 * A * B / (1 + A + A^2)
end

function m_from_e(E, ecc)
    @assert isfinite(E)
    @assert isfinite(ecc) && ecc >= 0 && ecc < 1

    return E - ecc * sin(E)
end

function m_from_f(F, ecc)
    @assert isfinite(F)
    @assert isfinite(ecc) && ecc > 1

    return ecc * sinh(F) - F
end

function m_from_d(D)
    @assert isfinite(D)

    return D + D^3/3
end

function m_from_nu(nu, ecc)
    @assert isfinite(nu)
    @assert isfinite(ecc) && ecc >= 0

    if ecc < 1
        return m_from_e(e_from_nu(nu, ecc), ecc)
    end

    if ecc > 1
        return m_from_f(f_from_nu(nu, ecc), ecc)
    end

    return m_from_d(d_from_nu(nu))
end
