module AstroUtils

export solve_and_print_solution, safe_acos
export keplerian_from_state_vectors, state_vectors_from_keplerian
export twobody
export hyperb, hyperb2
export time_to_next_radius, time_to_next_true_anomaly, mean_motion, true_anomaly_from_radius, sma_ecc_from_state_vectors, sma_from_state_vectors
export true_anomaly_from_radius, periapsis_from_state_vectors, period_from_state_vectors
export m_from_nu, m_from_d, m_from_f, m_from_e, d_from_m, f_from_m, e_from_m, nu_from_f, nu_from_e, f_from_nu, e_from_nu, d_from_nu, nu_from_d

include("angles.jl")
include("conversion.jl")
include("hyperb.jl")
include("misc.jl")
include("twobody.jl")
include("utils.jl")

end
