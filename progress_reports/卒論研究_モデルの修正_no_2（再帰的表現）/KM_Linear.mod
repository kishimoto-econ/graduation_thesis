var q_hat phi_til k_hat kp_hat b_hat  x_hat xp_hat mu_til  C_hat Y_hat;

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
    qss * kss * k_hat + (bss/betap) * b_hat(-1) + xss * x_hat
    = (qss * kss + (a+c) * kss) * k_hat(-1) + bss * b_hat;

    // (2)
    b_hat = q_hat(+1) + k_hat - eb;

    // (3)
    x_hat = k_hat(-1);

    // (4)
    phi_til = (beta/betap) * phi_til(+1) + (1/betap) * mu_til;

    // (5)
    qss * (1+phiss) * q_hat + qss * phi_til + beta*c*phi_til(+1)
    = (beta*(1+phiss)+muss)*qss*q_hat(+1)
    + beta*(a+c+qss)*phi_til(+1) + muss*qss*mu_til - muss*qss*eb;

    // (6)
    qss * q_hat
    = betap * alpha * (alpha-1) * (z+kpss) ^ (alpha-2) * kpss * kp_hat
    + betap * qss * q_hat(+1);

    // (7)
    xss * x_hat + m * xpss * xp_hat
    =(a+c) * kss * k_hat(-1) + m * alpha * (z+kpss) ^ (alpha-1)
    * kpss * kp_hat(-1);

    // (8)
    kss * k_hat + m * kpss * kp_hat = 0;

    // (9)
    C_hat = x_hat + m*xp_hat;

    // (10)
    Y_hat = C_hat; 

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
    Css = xss + m*xpss;
    Yss = Css;
end;

resid;
steady;
check;
model_diagnostics;

shocks;
    var eb = 0.01^2;
end;

stoch_simul(order=1,irf=12,ar=0) k_hat kp_hat b_hat q_hat x_hat xp_hat mu_til phi_til Y_hat;
M_.endo_names(oo_.dr.state_var,:)

