function J = mollweide_Jacobian(LAMBDA, PHI, R)

Lambda = LAMBDA * pi / 180;
Phi    = PHI * pi / 180;

t = Phi; dif = 1;
while dif > 1e-13
   t = t - (2*t + sin(2*t) - pi * sin(Phi)) ./ (2 + 2 * cos(2*t));
   dif = max(abs(2*t + sin(2*t) - pi * sin(Phi)));
end

J = R * sqrt(2) / pi * [2 * cos(t), -pi * Lambda .* sin(t) .* cos(Phi) ./ (1 + cos(2 * t));
                        0, pi^2 * cos(t) .* cos(Phi) ./ (2 * (1 + cos(2 * t)))];
