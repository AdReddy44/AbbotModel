function dM_dt = ode_M(t, M, tau_minus)
    dM_dt = -M / tau_minus;
end