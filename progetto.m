clc; close all; clear all;
set(0,'DefaultFigureWindowStyle','docked');

W = 1;
D = 1;
e_max = 0.1;
S = 0.15;
T_a5 = 0.01;


beta = 1;
J0 = 1.7e-4;
Ji =    [1e-5     7e-5    9e-6    2e-6];
phi =   [-0.04   2.9     2.8     -2.6];
k = [1 2 3 4];
J = @(x) J0 + Ji(1) * cos(x + phi(1))+ Ji(2) * cos(2*x + phi(2))+ Ji(3) * cos(3*x + phi(3))+ Ji(4) * cos(4*x + phi(4));
Jdot = @(x) -(Ji(1) * sin(x + phi(1))+ Ji(2) * 2 * sin(2*x + phi(2))+ Ji(3) * 3 * sin(3*x + phi(3))+ Ji(4) * 4 * sin(4*x + phi(4)));

DD = @(x) sin(0.1*k(1)*x) + sin(0.1*k(2)*x) + sin(0.1*k(3)*x) + sin(0.1*k(4)*x);
NN = @(x) sin((5e4)*k(1)*x) + sin((5e4)*k(2)*x) + sin((5e4)*k(3)*x) + sin((5e4)*k(4)*x);
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

s = tf('s');

modello = ss(A_e, B_e, C, D);
G = tf(modello);

xi = sqrt(log(S)^2/(pi^2+log(S)^2));
Mf = xi*100;

omega_cMin = 3 / (T_a5*xi);

mappatura_specifiche_bode(G, 'G', 5e4, 65, 0.5, 30, omega_cMin, Mf);

phi = deg2rad(-5);
M_star = 0.31; %circa -10dB
omega_star = 5e4;
tau = (cos(phi) - 1/M_star)/(omega_star * sin(phi));
alpha = (M_star - cos(phi))/(omega_star * sin(phi)) / tau;

ritardo = (1+tau * alpha * s)/(1+tau*s);
mu = omega_cMin * 1.5;
R = mu * ritardo;
LL =  R * G;

mappatura_specifiche_bode(LL, 'L', 5e4, 65, 0.5, 30, omega_cMin, Mf);

FF = LL / (1+LL);

figure();
hold on;
T_simulation = 2*T_a5;
[y_step,t_step] = step(W*FF, T_simulation);
plot(t_step,y_step,'b');
grid on, zoom on, hold on;

LV = W; % lim s->0 W*FF(s)

% vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0],[LV*(1+S),LV*(1+S),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento all'5%
patch([T_a5,T_simulation,T_simulation,T_a5],[LV*(1-0.05),LV*(1-0.05),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a5,T_simulation,T_simulation,T_a5],[LV*(1+0.05),LV*(1+0.05),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

ylim([0,LV*2]);

Legend_step = ["Risposta al gradino"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);

%% Check disturbo in uscita

% Funzione di sensitività
SS = 1/(1+LL);
figure();

% Simulazione disturbo in uscita a pulsazione 0.05
tt = 0:1e-2:2e2;
dd = DD(tt);
y_d = lsim(SS,dd,tt);
hold on, grid on, zoom on
plot(tt,dd,'m')
plot(tt,y_d,'b')
grid on
legend('d(t)','y_d(t)')

%% Check disturbo di misura
figure();

% Simulazione disturbo di misura a pulsazione 1000
tt = 0:1e-6:2e-3;
nn = NN(tt);
y_n = lsim(-FF,nn,tt);
hold on, grid on, zoom on
plot(tt,nn,'m')
plot(tt,y_n,'b')
grid on
legend('n(t)','y_n(t')


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

function [omega_plot_min, omega_plot_max, omega_cMax] = mappatura_specifiche_bode(GG, titolo, omega_n, A_n, omega_d, A_d, omega_cMin, Mf)
    figure()
    omega_plot_min = 1e-4;
    omega_plot_max = 5e6;
    omega_cMax = omega_n;
    [Mag, phase, omega] = bode(GG, {omega_plot_min, omega_plot_max});
    
    % Rumore (giallo)
    patch([omega_n, omega_plot_max, omega_plot_max, omega_n], [-A_n, -A_n, 200, 200], 'y', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
    
    % Pulsazione critica (azzurro)
    patch([omega_plot_min, omega_d, omega_d, omega_cMin, omega_cMin, omega_plot_min], [A_d,A_d, 0,0, -200, -200], 'c', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3); 
    grid on, hold on;
    margin(Mag, phase, omega);
    
    % Margine di fase (verde)
    patch([omega_cMin, omega_cMax, omega_cMax, omega_cMin], [-180 + Mf, -180 + Mf, -270, -270], 'g', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3); 
    grid on, hold on;
    title(titolo);
end