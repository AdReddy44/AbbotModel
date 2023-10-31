function dPa_dt = ode_Pa(t, Pa, tau_plus)
    dPa_dt = -Pa / tau_plus;
end