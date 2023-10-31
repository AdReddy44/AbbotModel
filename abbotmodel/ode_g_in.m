%ode for g_in
function dg_in_dt = ode_g_in(t, g_in, tau_in)
    dg_in_dt = -g_in / tau_in;
end