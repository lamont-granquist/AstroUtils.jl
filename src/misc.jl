using LinearAlgebra

include("angles.jl")

"""
    time_to_next_radius(mu::AbstractFloat, 
                        r::AbstractVector{<:AbstractFloat}, 
                        v::AbstractVector{<:AbstractFloat}, 
                        radius::AbstractFloat) -> AbstractFloat

Calculates the time to the next radial distance from state vectors.

# Arguments
- `mu`: Gravitational parameter
- `r`: Current position vector (3-element)
- `v`: Current velocity vector (3-element)
- `radius`: Target radial distance

# Returns
- Time in seconds until the body next reaches the given `radius`.  May return a negative
  value for a hyperbolic orbit which has no future encounter.
"""
function time_to_next_radius(
    mu::Number,
    r::AbstractVector{<:Number},
    v::AbstractVector{<:Number},
    radius::Number
)::AbstractFloat
    @assert isfinite(mu) && mu > 0
    @assert length(r) == 3 && norm(r) > 0
    @assert length(v) == 3
    @assert isfinite(radius) && radius > 0

    nu = true_anomaly_from_radius(mu, r, v, radius)

    t1 = time_to_next_true_anomaly(mu, r, v,  nu)
    t2 = time_to_next_true_anomaly(mu, r, v, -nu)

    if t1 >= 0 && t2 >= 0
        return min(t1, t2)
    elseif t1 < 0 && t2 < 0
        return max(t1, t2)
    else
        return t1 >= 0 ? t1 : t2
    end
end

function time_to_next_true_anomaly(mu, r, v, nu2)
    sma, ecc, _, _, _, nu1, _ = keplerian_from_state_vectors(mu, r, v)

    return time_to_next_true_anomaly(mu, sma, ecc, nu1, nu2);
end

function time_to_next_true_anomaly(mu, sma, ecc, nu1, nu2)
    mm = mean_motion(mu, sma)
    manom1 = m_from_nu(nu1, ecc)
    manom2 = m_from_nu(nu2, ecc)

    if ecc < 1
        return mod2pi( manom2 - manom1 ) / mm
    end
    return ( manom2 - manom1 ) / mm
end

function mean_motion(mu, sma)
    return sqrt(abs(mu / (sma * sma * sma)))
end

"""
    true_anomaly_from_radius(mu::AbstractFloat, r::AbstractVector{<:AbstractFloat},
                             v::AbstractVector{<:AbstractFloat}, radius::AbstractFloat)
                             -> AbstractFloat

Computes the true anomaly from state vectors and radius.

# Arguments
- `mu`: Gravitational parameter
- `r`: Position vector (3-element vector)
- `v`: Velocity vector (3-element vector)
- `radius`: Desired radius on orbit

# Returns
- True anomaly in radians
"""
function true_anomaly_from_radius(
        mu::AbstractFloat,
        r::AbstractVector{<:AbstractFloat},
        v::AbstractVector{<:AbstractFloat},
        radius::AbstractFloat
    )::AbstractFloat
    @assert isfinite(mu) && mu > 0
    @assert length(r) == 3 && norm(r) > 0
    @assert length(v) == 3
    @assert isfinite(radius) && radius > 0
    sma, ecc = sma_ecc_from_state_vectors(mu, r, v)
    return true_anomaly_from_radius(sma, ecc, radius)
end

"""
    sma_ecc_from_state_vectors(mu::Number, r::AbstractVector{<:Number}, v::AbstractVector{<:Number}) 
    -> Tuple{AbstractFloat, AbstractFloat}

Computes both the semi-major axis and eccentricity from state vectors.

# Arguments
- `mu`: Gravitational parameter
- `r`: Position vector (3-element)
- `v`: Velocity vector (3-element)

# Returns
- Tuple of:
    - Semi-major axis (a)
    - Eccentricity (e)
"""
function sma_ecc_from_state_vectors(
        mu::Number,
        r::AbstractVector{<:Number},
        v::AbstractVector{<:Number}
    )::Tuple{AbstractFloat, AbstractFloat}
    @assert isfinite(mu) && mu > 0
    @assert length(r) == 3 && norm(r) > 0
    @assert length(v) == 3
    h = cross(r, v)
    sma = sma_from_state_vectors(mu, r, v)
    return sma, sqrt(max(1 - dot(h,h) / (sma * mu), 0))
end

"""
    sma_from_state_vectors(mu::AbstractFloat, r::AbstractVector{<:AbstractFloat}, v::AbstractVector{<:AbstractFloat}) 
    -> AbstractFloat

Computes the semi-major axis of an orbit using state vectors.

# Arguments
- `mu`: Gravitational parameter
- `r`: Position vector (3-element)
- `v`: Velocity vector (3-element)

# Returns
- Semi-major axis (a)
"""
function sma_from_state_vectors(
        mu::Number,
        r::AbstractVector{<:Number},
        v::AbstractVector{<:Number}
    )::AbstractFloat
    @assert mu > 0
    @assert length(r) == 3 && norm(r) > 0
    @assert length(v) == 3
    return mu / (2 * mu / norm(r) - dot(v,v))
end

"""
    true_anomaly_from_radius(sma::Number, ecc::Number, radius::Number) -> AbstractFloat

Computes the true anomaly from keplerian elements and desired radius.

# Arguments
- `sma`: Semi-major axis of the orbit
- `ecc`: Eccentricity of the orbit
- `radius`: Desired radius on orbit

# Returns
- True anomaly in radians
"""
function true_anomaly_from_radius(sma::Number, ecc::Number, radius::Number)::AbstractFloat
    @assert isfinite(sma)
    @assert isfinite(ecc) && ecc > 0
    @assert isfinite(radius) && radius > 0
    slr = sma * (1 - ecc * ecc)
    return safe_acos((slr / radius - 1) / ecc)
end

function periapsis_from_state_vectors(mu, r, v)
    h = cross(r, v)
    sma = mu / (2.0 * mu / norm(r) - dot(v,v))
    ecc = sqrt(1 - dot(h,h)/(sma*mu))
    sma * (1 - ecc)
end

function period_from_state_vectors(μ, r, v)
    sma = μ / (2 * μ / norm(r) - dot(v,v))
    return 2*π*sqrt(sma^3/μ)
end
