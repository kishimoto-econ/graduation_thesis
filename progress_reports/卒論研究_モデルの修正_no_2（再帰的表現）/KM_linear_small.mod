var q_hat phi_til k_hat kp_hat b_hat;

varexo eb;

parameters
    alpha m K_bar betap beta a c z
    kss kpss bss qss xss xpss muss phiss
    ;

    alpha = 1/3;
    m = 0.5;
    K_bar = 1;
    betap = 0.99;
    beta = 0.98;
    a = 0.7;
    c = 0.3;
    z = 0.01;

model(linear);

    // (1)
    qss * kss * k_hat + (bss/betap) * b_hat(-1)
    = (qss * kss + (a+c) * kss) * k_hat(-1) + bss * b_hat;

    // (2)
    b_hat = q_hat(+1) + k_hat - eb;

    // (3)
    qss * (1+phiss) * q_hat - (muss*betap-1)*qss*phi_til
    = (beta*(1+phiss)+muss)*qss*q_hat(+1)
    + beta*(a+qss-muss*qss*betap)*phi_til(+1) - muss*qss*eb;

    // (4)
    qss * q_hat
    = betap * alpha * (alpha-1) * (z+kpss) ^ (alpha-2) * kpss * kp_hat
    + betap * qss * q_hat(+1);

    // (5)
    kss * k_hat + m * kpss * kp_hat = 0;

end;

steady_state_model;
    qss = a/(1-betap);
    kpss = (betap*alpha/a)^(1/(1-alpha)) - z;
    kss = K_bar - m*kpss;
    bss = betap*qss*kss;
    xpss = (1/m)*(a*kss + m*(z + kpss)^(alpha));
    phiss = (a*(beta-1) + beta*c)/(a*(1-beta));
    muss = (betap-beta)*beta*c/(a*(1-beta));
    xss = c*kss;
end;

resid;
steady;
check;
model_diagnostics;

shocks;
    var eb = 0.01^2;
end;

stoch_simul(order=1,irf=12,ar=0) k_hat kp_hat b_hat q_hat phi_til;
M_.endo_names(oo_.dr.state_var,:)

