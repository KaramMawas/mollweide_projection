close all;
clear all;
clc;

% program settings
R = 6378;
scale = 700;
Phi_intervall    = -90:15:90;
Lambda_intervall = -180:15:180;
Phi_intervall_Tissot    = -60:30:60;
Lambda_intervall_Tissot = -180:30:180;

% open new figure with white background
fig = figure('color', [1 1 1]);

% coastlines
load('coast.mat');
[x, y] = mollweide(long, lat, R);
plot(x, y, 'k'); hold on;

% meridians
[Lambda, Phi] = meshgrid(Lambda_intervall, min(Phi_intervall):1:max(Phi_intervall));
[x, y] = mollweide(Lambda, Phi, R);
plot(x, y, 'color', [0.5 0.5 0.5]);

% parallel circles
[Lambda, Phi] = meshgrid(min(Lambda_intervall):1:max(Lambda_intervall), Phi_intervall);
[x, y] = mollweide(Lambda, Phi, R);
plot(x', y', 'color', [0.5 0.5 0.5]);

% Tissot ellipses
for Lambda = Lambda_intervall_Tissot
    for Phi = Phi_intervall_Tissot
        [x, y] = mollweide(Lambda, Phi, R);
        
        % metric matrix of the source
        G = R^2 * [cos(Phi * pi / 180)^2 0; 0 1];
        
        % Jacobian
        J = mollweide_Jacobian(Lambda, Phi, R);
        
        % Cauchy-Green tensor
        C = J' * J;
        
        % solve the general eigenvalue problem
        [F, Lambda_12] = eig(C, G);
        
        % transform the eigenvectors from source to map
        f = J * F;
        
        % length of the semi axes
        lambda1 = sqrt(Lambda_12(1, 1));
        lambda2 = sqrt(Lambda_12(2, 2));
        
        % angle of semi major axis
        ang = atan2(f(2, 1), f(1, 1));
        
        % ATTENTION!
        % consider the following three cases:
        % 1.) no distortion: green circle of radius r = lambda1 = lambda2 == (equal) 1
        % 2.) conformal:     red circle of radius r = lambda1 = lambda2 ~= (not equal) 1
        % 3.) else:          red ellipse with semi major axis a = lambda1 and semi minor axis b = lambda2

        % As example, here: plot Tissot ellipse with semi major/minor axes
        ellipse(scale * lambda1, scale * lambda2, ang, x, y, 'r');
        
        ca = cos(ang);
        sa = sin(ang);
        
        ax(1, 1) = x - scale * ca * lambda1;
        ay(1, 1) = y - scale * sa * lambda1;
        ax(2, 1) = x + scale * ca * lambda1;
        ay(2, 1) = y + scale * sa * lambda1;

        ax(1, 2) = x + scale * sa * lambda2;
        ay(1, 2) = y - scale * ca * lambda2;
        ax(2, 2) = x - scale * sa * lambda2;
        ay(2, 2) = y + scale * ca * lambda2;
       
        plot(ax, ay, 'r');
        
    end
end

title('Mollweide Projection');
axis equal;
axis off;
