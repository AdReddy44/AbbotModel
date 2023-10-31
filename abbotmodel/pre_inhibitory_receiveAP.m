function g_in = pre_inhibitory_receiveAP(g_in, g_in_bar, tval)
    g_in(tval) = g_in(tval) + g_in_bar;
end
