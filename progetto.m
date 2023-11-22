W = 1;
D = 1;
e_max = 0.1;


beta = 1;
J0 = 1.7e-4;
Ji =    [1e-5     7e-5    9e-6    2e-6];
phi =   [-0.04   2.9     2.8     -2.6];
J = @(x) J0 + Ji(1) * cos(x + phi(1))+ Ji(2) * cos(2*x + phi(2))+ Ji(3) * cos(3*x + phi(3))+ Ji(4) * cos(4*x + phi(4));
Jdot = @(x) -(Ji(1) * sin(x + phi(1))+ Ji(2) * 2 * sin(2*x + phi(2))+ Ji(3) * 3 * sin(3*x + phi(3))+ Ji(4) * 4 * sin(4*x + phi(4)));
% x1 = theta, x2 = omega, u = Cm
% xdot = [ x2; (u - beta * x2) / J(x1) ] = f(x, u)
% y = x1 = h(x,u)

omega_e = 0;
theta_e = pi / 3;
u_e = 0; % beta * x2 = beta * omega_e = 0
x_e = [theta_e; omega_e];

%In questo caso A_e = [0 1; 0 -1/J(theta_e)]
A_e = [0 1; -(J(theta_e))^-2 * Jdot(theta_e) * (u_e - beta * omega_e)   -beta/J(theta_e)];
B_e = [0; 1/J(theta_e)];
C = [1 0];
D = 0;

modello = ss(A_e, B_e, C, D);
G = tf(modello);

figure(1);
bode(G, {0.5 1e5});

%animation(-1.5*pi:0.1:pi/2, 2);




function [] = animation(thetas, figure_n)
    wheel_center = [0 2];
    wheel_r = 0.6;
    vbarL_bottom = [2.4 0];
    vbarR_bottom = [4.4 0];
    vbar_length = 2;
    bar_length = 2.4;
    bar_c = [2 0];

    % circle: y^2 + (x-2)^2 = 4
    figure(figure_n);

    for i=1:length(thetas)
        clf(figure_n);
        grid on;
        axis([-1 6 0 4]);
        axis equal;
        theta = thetas(i);
        pin1 = wheel_center + wheel_r * [cos(theta), sin(theta)];
        [cx, cy] = findPointOnCircle(pin1, bar_length, bar_c, vbar_length);
        hbar_left = [cx cy];
        vbarL_top = hbar_left + [bar_length-vbar_length 0];
        vbarR_top = hbar_left + [bar_length 0];

        viscircles(wheel_center, wheel_r);
        %viscircles(bar_c, vbar_length);
        line([pin1(1) hbar_left(1)], [pin1(2) hbar_left(2)]);
        line([hbar_left(1) vbarR_top(1)], [hbar_left(2) vbarR_top(2)]);
        line([vbarL_top(1) vbarL_bottom(1)], [vbarL_top(2) vbarL_bottom(2)]);
        line([vbarR_top(1) vbarR_bottom(1)], [vbarR_top(2) vbarR_bottom(2)]);

        pause(0.1);
    end  
end



function [x_c, y_c] = findPointOnCircle(P, D, center, radius)
    [x_out, y_out] = circcirc(P(1), P(2), D, center(1), center(2), radius);
    x_c = x_out(1);
    y_c = y_out(1);
end


