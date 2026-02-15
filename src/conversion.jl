using LinearAlgebra

"""
    keplerian_from_state_vectors(μ::Number,
                                  r::AbstractVector{<:Number},
                                  v::AbstractVector{<:Number})
                                  -> Tuple{AbstractFloat, AbstractFloat, AbstractFloat,
                                           AbstractFloat, AbstractFloat, AbstractFloat, AbstractFloat}

Converts position and velocity state vectors into Keplerian orbital elements.

# Arguments
- `μ`: Standard gravitational parameter
- `r`: Position vector (3-element)
- `v`: Velocity vector (3-element)

# Returns
A tuple containing the following Keplerian elements (all scalars):

- `sma`: Semi-major axis (same units as position)
- `ecc`: Eccentricity (unitless)
- `inc`: Inclination (radians)
- `lan`: Longitude of ascending node (Ω, radians)
- `argp`: Argument of periapsis (ω, radians)
- `nu`: True anomaly (ν, radians)
- `l`: Semilatus rectum (same units as position)

All angular elements are in radians and lie within the interval [0, 2π) when appropriate.
"""
function keplerian_from_state_vectors(
    μ::Number,
    r::AbstractVector{<:Number},
    v::AbstractVector{<:Number}
)::Tuple{AbstractFloat, AbstractFloat, AbstractFloat, AbstractFloat, AbstractFloat, AbstractFloat, AbstractFloat}
    @assert μ > 0
    @assert length(r) == 3 && norm(r) > 0
    @assert length(v) == 3

    rmag = sqrt(sum(i^2 for i in r))         # Magnitude of position vector
    vmag = sqrt(sum(i^2 for i in v))         # Magnitude of velocity vector
    rhat = r ./ rmag                         # Unit vector in the direction of r
    hv = cross(r, v)                         # Specific angular momentum vector
    hhat = hv ./ sqrt(sum(i^2 for i in hv))  # Unit vector of angular momentum
    eccvec = cross(v / μ, hv) - rhat         # Eccentricity vector
    sma = 1.0 / (2.0 / rmag - vmag^2 / μ)    # Semi-major axis
    slr = sum(i^2 for i in hv) / μ           # Semi-latus rectum

    # Parameters for frame transformation
    d = 1.0 + hhat[3]
    p = d == 0 ? 0 : hhat[1] / d
    q = d == 0 ? 0 : -hhat[2] / d
    const1 = 1.0 / (1.0 + p^2 + q^2)

    fhat = [
        const1 * (1.0 - p^2 + q^2),
        const1 * 2.0 * p * q,
        -const1 * 2.0 * p
    ]

    ghat = [
        const1 * 2.0 * p * q,
        const1 * (1.0 + p^2 - q^2),
        const1 * 2.0 * q
    ]

    # Calculate Keplerian elements
    h = dot(eccvec, ghat)
    xk = dot(eccvec, fhat)
    x1 = dot(r, fhat)
    y1 = dot(r, ghat)
    λt = atan(y1, x1)                             # True longitude
    ecc = sqrt(h^2 + xk^2)                        # Eccentricity
    inc = 2.0 * atan(sqrt(p^2 + q^2))             # Inclination
    lan = inc > eps() ? atan(p, q) : 0.0          # Longitude of ascending node
    argp = ecc > eps() ? atan(h, xk) - lan : 0.0  # Argument of periapsis
    nu = λt - lan - argp                          # True anomaly

    # Normalize angles to [0, 2π]
    lan = mod2pi(lan)
    argp = mod2pi(argp)
    nu = mod2pi(nu)

    return sma, ecc, inc, lan, argp, nu, slr
end

function state_vectors_from_keplerian(μ, sma, ecc, inc, lan, argp, tanom)
    # semi-latus rectum
    slr = sma * (1 - ecc * ecc)

    # magnitude of radius
    rm = slr / (1 + ecc * cos(tanom))

    # argument of latitude
    arglat = argp + tanom

    sarglat = sin(arglat)
    carglat = cos(arglat)

    c4 = sqrt(μ / slr)
    c5 = ecc * cos(argp) + carglat
    c6 = ecc * sin(argp) + sarglat

    sinc = sin(inc)
    cinc = cos(inc)

    slan = sin(lan)
    clan = cos(lan)

    r = [
         rm * (clan * carglat - slan * cinc * sarglat),
         rm * (slan * carglat + cinc * sarglat * clan),
         rm * sinc * sarglat,
        ]

    v = [
         -c4 * (clan * c6 + slan * cinc * c5),
         -c4 * (slan * c6 - clan * cinc * c5),
         c4 * c5 * sinc,
        ]

    return r, v
end
