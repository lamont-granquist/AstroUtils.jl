using Optim
using LinearAlgebra

# Single-impulse transfer from an ellipitical, non-coplanar parking
# orbit to an arbitrary hyperbolic v-infinity target.
#
# mu - gravitational parameter of the body (scalar)
# r0 - reference position on the parking orbit (3x1)
# v0 - reference velocity on the parking orbit (3x1)
# v_inf - target v-infinity (3x1)
#
# vneg - velocity vector on the parking orbit at the impulsive burn (3x1)
# vpos - velocity vector on the hyperbolic orbit at the impuslvie burn (3x1)
# r - radius vector to the point of the impulsive burn (3x1)
# dt - coasting time on the parking orbit to the impulsive burn from the reference
#
# Ocampo, C., & Saudemont, R. R. (2010). Initial Trajectory Model for a Multi-Maneuver Moon-to-Earth Abort Sequence.
# Journal of Guidance, Control, and Dynamics, 33(4), 1184–1194.
#
function hyperb(mu, r0, v0, v_inf)
  # do a 1-dimensional bounded line search to find the rotation that minimizes the ejection burn
  result = optimize(rot -> first(hyperb2(mu, r0, v0, v_inf, rot)), -pi, pi)
  rot = Optim.minimizer(result)
  return hyperb2(mu, r0, v0, v_inf, rot)
end

# Main algorithm, this is from section 3 of the one-impulse case.
#
# rot - angle of rotation around v_inf to select hf of ejection hyperbola
#
# In the case where rot = 0 this recovers the algorithm in section 2 where
# hf is selected to be closest to h0.  For highly inclined cases, it may be sufficiently
# accurate to call this function with rot = 0.
#
# For perfectly coplanar transfers instead of implementing the algorithms in 1 and 2 the
# rotation is instead simply applied to the r1 position of the burn instead of treating the
# circular case specially.
#
function hyperb2(mu, r0, v0, v_inf, rot)
  h0 = cross(r0,v0)

  # semi major axis of parking orbit
  a0 = 1.0 / (2.0 / norm(r0) - norm(v0)^2 / mu)

  # sma of hyperbolic ejection orbit
  af = - mu / norm(v_inf)^2

  # parking orbit angular momentum unit
  h0_hat = h0/norm(h0)

  # eccentricity vector
  ecc = cross(v0,h0)/mu - r0/norm(r0)

  # eccentricity of the parking orbit.
  e0 = norm(ecc)

  # semilatus rectum of parking orbit
  p0 = a0 * ( 1 - e0^2 )

  # parking orbit periapsis position unit vector
  if e0 != 0
    rp0_hat = ecc/e0
  else
    rp0_hat = r0/norm(r0)
  end

  # parking orbit periapsis velocity unit vector
  vp0_hat = cross(h0, rp0_hat)
  vp0_hat = vp0_hat/norm(vp0_hat)

  # direction of hyperbolic v-infinity vector
  v_inf_hat = v_inf/norm(v_inf)

  # 3 cases for finding hf_hat
  dotp = dot(h0_hat, v_inf_hat)
  if abs(dotp) == 1
    # 90 degree plane change case
    hf_hat = cross(rp0_hat, v_inf_hat)
    hf_hat = hf_hat / norm(hf_hat)
  else
    # general case (also works for dotp = 0)
    hf_hat = cross(v_inf_hat, cross(h0_hat, v_inf_hat))
    hf_hat = hf_hat / norm(hf_hat)
  end

  if abs(dot(h0_hat, v_inf_hat)) > eps()
    # if the planes are not coincident, rotate hf_hat by applying rodrigues formula around v_inf_hat
    hf_hat = hf_hat * cos(rot) + cross(v_inf_hat, hf_hat) * sin(rot) + v_inf_hat * dot(v_inf_hat, hf_hat) * (1 - cos(rot))

    # unit vector pointing at the position of the burn on the parking orbit
    r1_hat = sign(dot(h0_hat, v_inf_hat)) * cross(h0_hat, hf_hat)
    r1_hat = r1_hat/norm(r1_hat)
  else
    # unit vector pointing at the position of the burn on the parking orbit
    r1_hat = cross(v_inf_hat, hf_hat)
    r1_hat = r1_hat/norm(r1_hat)

    # if the planes are coincident, rotate r1_hat by applying rodrigues formula around h0_hat
    r1_hat = r1_hat * cos(rot) + cross(h0_hat, r1_hat) * sin(rot) + h0_hat * dot(h0_hat, r1_hat) * (1 - cos(rot))
  end

  # true anomaly of r1 on the parking orbit
  nu_10 = sign(dot(h0_hat, cross(rp0_hat, r1_hat))) * acos(dot(rp0_hat,r1_hat))

  # length of the position vector of the burn on the parking orbit
  r1 = p0 / ( 1 + e0 * cos(nu_10) )

  # position of the burn
  r = r1 * r1_hat

  # constant
  k = - af / r1

  # angle between v_inf and the r1 burn
  delta_nu = acos(min(max(dot(r1_hat, v_inf_hat), -1), 1))

  # eccentricity of the hyperbolic ejection orbit
  sindnu  = sin(delta_nu)
  sin2dnu = sindnu * sindnu
  cosdnu  = cos(delta_nu)
  #ef = max(sqrt(sin2dnu + 2*k^2 + 2*k*(1-cosdnu) + sindnu*sqrt(sin2dnu + 4*k*(1-cosdnu)))/(sqrt(2)*k), 1+eps())
  ef = sqrt(sin2dnu + 2*k^2 + 2*k*(1-cosdnu) + sindnu*sqrt(sin2dnu + 4*k*(1-cosdnu)))/(sqrt(2)*k)

  # semilatus rectum of hyperbolic ejection orbit
  pf = af * ( 1 - ef^2 )

  # true anomaly of the v_inf on the hyperbolic ejection orbit
  nu_inf = acos(-1/ef)

  # true anomaly of the burn on the hyperbolic ejection orbit
  nu_1f = acos(min(max(-1/ef * cosdnu + sqrt(ef*ef-1)/ef * sindnu,-1), 1))

  # turning angle of the hyperbolic orbit
  delta = 2 * asin(1/ef)

  # incoming hyperbolic velocity unit vector
  v_inf_minus_hat = cos(delta) * v_inf_hat + sin(delta) * cross(v_inf_hat, hf_hat)

  # periapsis position and velocity vectors of the hyperbolic ejection orbit
  rpf_hat = v_inf_minus_hat - v_inf_hat
  rpf_hat = rpf_hat / norm(rpf_hat)
  vpf_hat = v_inf_minus_hat + v_inf_hat
  vpf_hat = vpf_hat / norm(vpf_hat)

  # compute the velocity on the hyperbola and the parking orbit
  vpos = sqrt(mu/pf) * ( -sin(nu_1f) * rpf_hat + ( ef + cos(nu_1f) ) * vpf_hat )
  vneg = sqrt(mu/p0) * ( -sin(nu_10) * rp0_hat + ( e0 + cos(nu_10) ) * vp0_hat )

  # compute the cost
  dv = norm( vpos - vneg )

#  if nargout > 3
    # compute nu of the reference position on the parking orbit
    r0_hat = r0/norm(r0)
    nu0 = sign(dot(h0_hat, cross(rp0_hat, r0_hat))) * acos(dot(rp0_hat,r0_hat))

    # mean angular motion of the parking orbit
    n = sqrt( mu / a0^3 )

    # eccentric anomalies of reference position and r1 on the parking orbit
    E0 = atan(sqrt(1 - e0^2) * sin(nu0), e0 + cos(nu0))
    E1 = atan(sqrt(1 - e0^2) * sin(nu_10), e0 + cos(nu_10))

    # mean anomalies of reference position and r1 on the parking orbit
    M0 = E0 - e0 * sin( E0 )
    M1 = E1 - e0 * sin( E1 )

    # coast time on the parking orbit
    dt = ( M1 - M0 ) / n
    if dt < 0
      dt = dt + 2 * pi / n
    end
#  end

  return dv, vneg, vpos, r, dt
end
