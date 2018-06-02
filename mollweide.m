function [x, y] = mollweide(LAMBDA, PHI, R)
% Name, given name, matrikulation number
% 
% calculates the Mollweide map projection

Lambda = LAMBDA * pi / 180;
Phi    = PHI * pi / 180;

t = Phi; dif = 1;
while dif > 1e-13
    t = t - (2 * t + sin(2 * t) - pi * sin(Phi)) ./ (2 + 2 * cos(2 * t));
    dif = max(abs(2 * t + sin(2 * t) - pi * sin(Phi)));
end

x = sqrt(2) * R / pi * 2 * Lambda .* cos(t);
y = sqrt(2) * R * sin(t);