var
    x   ${x}$           (long_name='consumption farmer')
    xp  ${x^\prime}$    (long_name='consumption gatherer')
    b   $b$             (long_name='bonds')
    k   ${k}$           (long_name='Land held by farmer')
    kp  ${k^\prime}$    (long_name='Land held by gatherer')
    q   $q$             (long_name='fruit price')
    mu  $\mu$           (long_name='multiplier enforceability constraint')  
    phi $\phi$          (long_name='multiplier consumption constraint farmer')  
    C   $C$             (long_name='aggregate consumption')
    Y   $Y$             (long_name='aggregate output')
    ;

varexo eps $\varepsilon$;

parameters
    alpha       $\alpha$            (long_name='exponent production function gatherers')
    m           $m$                 (long_name='population size gatherers')
    K_bar       $\bar K$            (long_name='land supply')
    beta        $\beta$             (long_name='discount factor farmer')
    betap       ${\beta^\prime}$    (long_name='discount factor farmer')
    a           $a$                 (long_name='tradability share')
    c           $c$                 (long_name='non-tradability share') 
    z           $z$                 (long_name='constant production function gatherers')
    theta       $\theta$            
    ;

alpha = 0.05;  % 1/3 から 0.05 へ変更（gathererの土地価値の感度を最大に）
m = 0.5;
K_bar = 1;
betap = 0.99;
beta = 0.98;
a = 0.3;  
c = 0.1;
z = 0.01;
theta = 0.99;

model;
  q * (k - k(-1)) + (1/betap) * b(-1) + x - c * k(-1)^theta = a * k(-1)^theta + b;
  b = betap * q(+1) * k * (1-eps);
  x = c * k(-1)^theta;
  1 + phi = (beta * (1 + phi(+1)) + mu) / betap;
  q * (1 + phi) + beta * c * phi(+1) * theta * k^(theta-1) = 
      beta * (1 + phi(+1)) * (theta * (a + c) * k(+1)^(theta-1) + q(+1)) + mu * q(+1) * (1-eps);
  q = betap * (alpha * (z + kp)^(alpha-1) + q(+1));
  x + m * xp = (a + c) * k(-1)^theta + m * (z + kp(-1))^alpha;
  k + m * kp = K_bar;
  C = x + m*xp;
  Y = C;
end;

initval;
  phi = (theta * beta * c) / (a * (1 - theta * beta)) - 1;
  mu  = (betap - beta) * (1 + phi);

  % k の「最初の推測値」として 0.5 を与える
  k = 0.5; 

  kp = (K_bar - k) / m;
  q       = (a / (1 - betap)) * k^(theta - 1);
  b       = betap * q * k;
  x       = c * k^theta;
  xp = (a * k^theta) / m + (z + kp)^alpha;
end;

shocks;
  var eps = 0.01;
end;

steady;
check;

stoch_simul(order=1, irf=12, Tex) k kp Y q mu b phi x xp;