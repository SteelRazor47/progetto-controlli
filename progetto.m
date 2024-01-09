clc; close all; clearvars;
set(0,'DefaultFigureWindowStyle','docked');

%% Inizializzazione dei parametri 
W = 1;
D = 1;
e_max = 0.04;
S = 0.10;
T_a5 = 0.008;
A_d = 30;
A_n = 65;
omega_d_max = 0.5;
omega_n_min = 5e4;

beta = 1;
J0 = 1.7;
Ji =    [0.1     0.7    0.09    0.02];
phi =   [-0.04   2.9    2.8     -2.6];
k = 50;
J = @(x) J0 + Ji(1) * cos(x + phi(1))+ Ji(2) * cos(2*x + phi(2))+ Ji(3) * cos(3*x + phi(3))+ Ji(4) * cos(4*x + phi(4));
Jdot = @(x) -(Ji(1) * sin(x + phi(1))+ Ji(2) * 2 * sin(2*x + phi(2))+ Ji(3) * 3 * sin(3*x + phi(3))+ Ji(4) * 4 * sin(4*x + phi(4)));

DD = noise_gen(0.1);
NN = noise_gen(5e4);


%% Valori all'equilibrio
omega_e = 0;
theta_e = pi / 3;
u_e = beta * omega_e + k * theta_e;
x_e = [theta_e; omega_e];

%% Linearizzazione
df2dx1 = (( k*theta_e+beta*omega_e-u_e )*Jdot(theta_e) - k*J(theta_e))/J(theta_e)^2;
A_e = [0 1; df2dx1    -beta/J(theta_e)];
B_e = [0; 1/J(theta_e)];
C_e = [1 0];
D_e = 0;

%% Creazione del modello

s = tf('s');

modello = ss(A_e, B_e, C_e, D_e);
G = tf(modello);

Gdisp = zpk(G);
Gdisp.DisplayFormat = 'Frequency';
display(Gdisp);

xi = sqrt(log(S)^2/(pi^2+log(S)^2));
Mf = max(xi*100, 30);

omega_cMin = 300 / (T_a5*Mf);

mappatura_specifiche_bode(G, 'G', omega_n_min, A_n, omega_d_max, A_d, omega_cMin, Mf);

%% Regolatore statico - proporzionale senza poli nell'origine

% valore minimo prescritto per L(0)
mu_s_error = (D+W)/e_max;
mu_s_dist  = 10^(A_d/20);

% guadagno minimo del regolatore ottenuto come L(0)/G(0)
G_0 = abs(evalfr(G,0));
G_omega_d_max = abs(evalfr(G,1i*omega_d_max));

R_s = max(mu_s_error/G_0,mu_s_dist/G_omega_d_max);

% Sistema esteso
G_e = R_s*G;

mappatura_specifiche_bode(G_e, 'G_e', omega_n_min, A_n, omega_d_max, A_d, omega_cMin, Mf);

%% Regolatore dinamico
% Ci troviamo nel caso B quindi progettiamo ua rete anticipatrice

Mf_star = Mf+5;
omega_c_star = 750;

mag_omega_c_star    = abs(evalfr(G_e,1j*omega_c_star));
arg_omega_c_star    = rad2deg(angle(evalfr(G_e,1j*omega_c_star)));

M_star = 1/mag_omega_c_star;
phi_star = deg2rad(Mf_star - 180 - arg_omega_c_star);

tau = (M_star - cos(phi_star))/(omega_c_star*sin(phi_star));
alpha_tau = (cos(phi_star) - 1/M_star)/(omega_c_star*sin(phi_star));
alpha = alpha_tau / tau;

R_d = (1+tau * s)/(1+alpha*tau*s);
R = R_s * R_d;
LL =  R * G;

mappatura_specifiche_bode(LL, 'L', omega_n_min, A_n, omega_d_max, A_d, omega_cMin, Mf);

%% Check prestazioni in anello chiuso

FF = LL / (1+LL);

figure();
grid on, zoom on, hold on;
T_simulation = 2*T_a5;
[y_step,t_step] = step(W*FF, T_simulation);
plot(t_step,y_step,'b');

LV = evalfr(W*FF,0); % lim s->0 W*FF(s)

% vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0],[LV*(1+S),LV*(1+S),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento al 5%
patch([T_a5,T_simulation,T_simulation,T_a5],[LV*(1-0.05),LV*(1-0.05),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a5,T_simulation,T_simulation,T_a5],[LV*(1+0.05),LV*(1+0.05),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

ylim([0,LV*2]);

Legend_step = ["Risposta al gradino (lin.)"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);
xlabel('t(s)');
ylabel('\theta (rad)');

%% Check disturbo in uscita

% Funzione di sensitività
SS = 1/(1+LL);
figure();

% Simulazione disturbo in uscita a pulsazione 0.1-0.4
tt = 0:1e-2:2e2;
dd = DD(tt);
y_d = lsim(SS,dd,tt);
hold on, grid on, zoom on
plot(tt,dd,'m')
plot(tt,y_d,'b')
legend('d(t)','y_d(t)')
xlabel('t(s)');
ylabel('\theta (rad)');

%% Check disturbo di misura
figure();

% Simulazione disturbo di misura a pulsazione 50k-200k
tt = 0:1e-6:2e-3;
nn = NN(tt);
y_n = lsim(-FF,nn,tt);
hold on, grid on, zoom on
plot(tt,nn,'m')
plot(tt,y_n,'b')
legend('n(t)','y_n(t)')
xlabel('t(s)');
ylabel('\theta (rad)');

%% Check prestazioni linearizzato
x_sim = x_e;
W_sim = W;
T_simulation = 10;
out = sim("SdC_lineare_progetto.slx","StopTime",num2str(T_simulation));
LV = out.y_sim.Data(end);

T_simulation = 3*T_a5;
out = sim("SdC_lineare_progetto.slx","StopTime",num2str(T_simulation));
data = out.y_sim.Data(:);

figure();
hold on; zoom on; grid on;

plot(out.y_sim.Time, data);

% vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0],[LV*(1+S),LV*(1+S),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento al 5%
patch([T_a5,T_simulation,T_simulation,T_a5],[LV*(1-0.05),LV*(1-0.05),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a5,T_simulation,T_simulation,T_a5],[LV*(1+0.05),LV*(1+0.05),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

ylim([0,LV*2]);

Legend_step = ["Risposta al gradino (lin. con dist.)"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);
xlabel('t(s)');
ylabel('\theta (rad)');

%% Check prestazioni non linearizzato
x_sim = x_e;
W_sim = W;
T_simulation = 10;
out = sim("SdC_progetto_fast.slx","StopTime",num2str(T_simulation));
LV = out.y_sim.Data(end) - theta_e;

T_simulation = 3*T_a5;
out = sim("SdC_progetto_fast.slx","StopTime",num2str(T_simulation));
data = out.y_sim.Data(:) - theta_e;

figure();
hold on; zoom on; grid on;

plot(out.y_sim.Time, data);

% vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0],[LV*(1+S),LV*(1+S),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento al 5%
patch([T_a5,T_simulation,T_simulation,T_a5],[LV*(1-0.05),LV*(1-0.05),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a5,T_simulation,T_simulation,T_a5],[LV*(1+0.05),LV*(1+0.05),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

ylim([0,LV*2]);

Legend_step = ["Risposta al gradino (non lin.)"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);
xlabel('t(s)');
ylabel('\theta - \theta_e (rad)');

%% Check prestazioni non linearizzato con disturbi
x_sim = x_e;
W_sim = W;

T_simulation = 10;
out = sim("SdC_progetto.slx","StopTime",num2str(T_simulation));
LV = out.y_sim.Data(end) - theta_e;

T_simulation = 3*T_a5;
out = sim("SdC_progetto.slx","StopTime",num2str(T_simulation));
data = out.y_sim.Data(:) - theta_e;

figure();
hold on; zoom on; grid on;

plot(out.y_sim.Time, data);

% vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0],[LV*(1+S),LV*(1+S),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento al 5%
patch([T_a5,T_simulation,T_simulation,T_a5],[LV*(1-0.05),LV*(1-0.05),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a5,T_simulation,T_simulation,T_a5],[LV*(1+0.05),LV*(1+0.05),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

ylim([0,LV*2]);

Legend_step = ["Risposta al gradino (non lin. senza dist.)"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);
xlabel('t(s)');
ylabel('\theta - \theta_e (rad)');

%% Confronto
T_simulation = 0.01; 
x_sim = x_e;
W_sim = W;
outl = sim("SdC_lineare_progetto.slx","StopTime",num2str(T_simulation));
out = sim("SdC_progetto.slx","StopTime",num2str(T_simulation));
out_nodist = sim("SdC_progetto_fast.slx","StopTime",num2str(T_simulation));
data = out.y_sim.Data(:) - theta_e;
data_nodist = out_nodist.y_sim.Data(:) - theta_e;

figure();
hold on; zoom on; grid on;
plot(outl.y_sim,'b');
plot(out.y_sim.Time, data, 'r');
plot(out_nodist.y_sim.Time, data_nodist, 'm');

Legend_step = ["Risposta al gradino (lin.)"; "Risposta al gradino (non lin. con dist.)"; "Risposta al gradino (non lin. senza dist.)"];
legend(Legend_step);
xlabel('t(s)');
ylabel('\theta - \theta_e (rad)');

%% Confronto regime
T_simulation = 10; 
x_sim = x_e;
W_sim = W;
outl = sim("SdC_lineare_progetto.slx","StopTime",num2str(T_simulation));
out = sim("SdC_progetto.slx","StopTime",num2str(T_simulation));
out_nodist = sim("SdC_progetto_fast.slx","StopTime",num2str(T_simulation));
data = out.y_sim.Data(:) - theta_e;
data_nodist = out_nodist.y_sim.Data(:) - theta_e;

figure();
hold on; zoom on; grid on;
plot(outl.y_sim,'b');
plot(out.y_sim.Time, data, 'r');
plot(out_nodist.y_sim.Time, data_nodist, 'm');

Legend_step = ["Risposta al gradino (lin.)"; "Risposta al gradino (non lin. con dist.)"; "Risposta al gradino (non lin. senza dist.)"];
legend(Legend_step);
xlabel('t(s)');
ylabel('\theta - \theta_e (rad)');
%% Punti opzionali

% Angolo iniziale
figure();
hold on, grid on, zoom on;

W_sim = 0;
theta_range = -pi:pi/3:pi;
for dtheta = theta_range
    x_sim = x_e + [dtheta; 0];
    out = sim("SdC_progetto_fast.slx","StopTime","8");
    data = out.y_sim.Data(:) - theta_e;
    plot(out.y_sim.Time, data);
end
yline([-e_max e_max], 'HandleVisibility', 'off');
legendCell = [  "\Delta\theta = -\pi";   "\Delta\theta = -2\pi/3"; "\Delta\theta = -\pi/3";  "\Delta\theta = 0"; 
                "\Delta\theta = \pi/3"; "\Delta\theta = 2\pi/3"; "\Delta\theta = \pi"];
legend(legendCell);
xlabel('t(s)');
ylabel('\theta - \theta_e (rad)');

% Velocità iniziale
figure();
hold on, grid on, zoom on;

W_sim = 0;
omega_range = 10.^(-1:3);
for vel = omega_range
    x_sim = x_e + [0; vel];
    out = sim("SdC_progetto_fast.slx","StopTime","8");
    data = out.y_sim.Data(:) - theta_e;
    p = plot(out.y_sim.Time, data);

    x_sim = x_e + [0; -vel];
    out = sim("SdC_progetto_fast.slx","StopTime","8");
    data = out.y_sim.Data(:) - theta_e;
    plot(out.y_sim.Time, data, 'Color', p.Color, 'HandleVisibility', 'off');
end
yline([-e_max e_max], 'HandleVisibility', 'off');
legendCell = strcat('\Delta\omega = \pm',string(num2cell(omega_range)));
legend(legendCell);
xlabel('t(s)');
ylabel('\theta - \theta_e (rad)');


%% Gradini senza disturbi - W in [0, 3]

x_sim = x_e;
W_range = 0.01:0.01:3;
LVs = cell(length(W_range), 1);
TAs = cell(length(W_range), 1);
SEs = cell(length(W_range), 1);

T_simulation = 10;
fprintf("%6.2f%%\n", 0);
for i = 1:length(W_range)
    W_sim = W_range(i);
    out = sim("SdC_progetto_fast.slx","StopTime",num2str(T_simulation));
    data = out.y_sim.Data(:) - theta_e;

    LVs{i} = data(end);
    t_ind = find(data >= LVs{i}*1.05 | data <= LVs{i}*0.95, 1, 'last');
    TAs{i} = out.y_sim.Time(t_ind + 1); % +1 perchè a t_ind è fuori dal range
    SEs{i} = max(data)/LVs{i} - 1;

    %progress bar
    fprintf(repmat('\b', 1, 8));
    fprintf("%6.2f%%\n", W_sim/3*100);
end

% Plot  gradini
figure();
data = abs(cell2mat(LVs) - W_range');
plot(W_range, data);
yline(e_max);
xlabel('W(rad)');
ylabel('e_{\infty}(rad)');
grid on, hold on, zoom on;

figure();
plot(W_range, cell2mat(SEs) * 100);
yline(S * 100);
xlabel('W(rad)');
ylabel('Sovraelungazione %');
grid on, hold on, zoom on;

figure();
semilogy(W_range, cell2mat(TAs));
yline(T_a5);
xlabel('W(rad)');
ylabel('Tempo di assestamento(log_{10} s)');
grid on, hold on, zoom on;

%% Range di gradini validi
error = cell2mat(LVs) - W_range';
W_err_range = W_range(error < e_max);
W_se_range = W_range(cell2mat(SEs) < S);
W_ta_range = W_range(cell2mat(TAs)< T_a5);
W_range_valid = intersect(intersect(W_err_range, W_se_range), W_ta_range);
disp(W_range_valid);

%% Animazione
x_sim = x_e;
W_sim = W;
out = sim("SdC_progetto.slx","StopTime","0.01");
animation(out.y_sim);

%% Funzioni
function [out] = noise_gen(omega)
    function [res] = sum(x)
        res = 0;
        for k = 1:4
            res = res + sin(omega * k * x);
        end
    end
    out = @sum;
end

function [] = animation(y_sim)
    frame_length = max(floor(length(y_sim.Data)/60), 1); %circa 15fps per 4s
    thetas = y_sim.Data(1:frame_length:end);
    tt = y_sim.Time(1:frame_length:end);


    wheel_center = [0 2];
    wheel_r = 0.6;
    vbarL_bottom = [2.4 0];
    vbarR_bottom = [4.4 0];
    vbar_length = 2;
    bar_length = 2.4;
    bar_c = [2 0];
    pinbar_length = norm((wheel_center+[0 wheel_r]) - (bar_c + [0 vbar_length]));

    gifFile = 'animazione.gif';
    f = figure();
    for i=1:length(thetas)
        clf(f);
        grid on;
        axis equal;
        axis manual;
        axis([-1 6 -1 4]);

        theta = pi/2-thetas(i);
        pin1 = wheel_center + wheel_r * [cos(theta), sin(theta)];
        [cx, cy] = findPointOnCircle(pin1, pinbar_length, bar_c, vbar_length);
        hbar_left = [cx cy];
        vbarL_top = hbar_left + [bar_length-vbar_length 0];
        vbarR_top = hbar_left + [bar_length 0];

        pos = [wheel_center-(wheel_r*1.25) 2*(wheel_r*1.25) 2*(wheel_r*1.25)];
        rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', [0 0 0 0.1]);

        %pos = [bar_c-vbar_length 2*vbar_length 2*vbar_length];
        %rectangle('Position',pos,'Curvature',[1 1]);
        
        patch([-1,6,6,-1],[0,0,-1,-1],'g','FaceAlpha',0.1,'EdgeAlpha',0.1, 'FaceColor', 'Black');

        vbar_angle = atan((vbarL_top(2) - vbarL_bottom(2))/(vbarL_top(1)-vbarL_bottom(1)));
        dx = 0.2 * cos(vbar_angle);
        dy = 0.2 *sin(vbar_angle);
        line([vbarL_top(1)+dx vbarL_bottom(1)-dx], [vbarL_top(2)+dy vbarL_bottom(2)-dy], 'Color', [1 0 0 0.2], 'LineWidth', 10);
        line([vbarR_top(1)+dx vbarR_bottom(1)-dx], [vbarR_top(2)+dy vbarR_bottom(2)-dy], 'Color', [1 0 0 0.2], 'LineWidth', 10);
        line([hbar_left(1)-0.2 vbarR_top(1)+0.2], [hbar_left(2) vbarR_top(2)], 'Color', [1 0 0 0.2], 'LineWidth', 10);
        hbarL_angle = atan((hbar_left(2) - pin1(2))/(hbar_left(1)-pin1(1)));
        dx = 0.1 * cos(hbarL_angle);
        dy = 0.1 *sin(hbarL_angle);
        line([pin1(1)-dx hbar_left(1)+dx], [pin1(2)-dy hbar_left(2)+dy], 'Color', [1 0 0 0.2], 'LineWidth', 6);

        line([pin1(1) hbar_left(1)], [pin1(2) hbar_left(2)], 'Color', 'Black');
        line([hbar_left(1) vbarR_top(1)], [hbar_left(2) vbarR_top(2)], 'Color', 'Black');
        line([vbarL_top(1) vbarL_bottom(1)], [vbarL_top(2) vbarL_bottom(2)], 'Color', 'Black');
        line([vbarR_top(1) vbarR_bottom(1)], [vbarR_top(2) vbarR_bottom(2)], 'Color', 'Black');


        pos = [pin1-0.05 2*0.05 2*0.05];
        rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', [0 0 0 0.1]);
        %text(pin1(1), pin1(2)+0.3, 'A');

        pos = [hbar_left-0.05 2*0.05 2*0.05];
        rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', [0 0 0 0.1]);
        %text(hbar_left(1), hbar_left(2)+0.3, 'B');

        pos = [vbarL_top-0.09 2*0.09 2*0.09];
        rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', [0 0 0 0.1]);
        %text(vbarL_top(1), vbarL_top(2)+0.3, 'C');
        pos = [vbarR_top-0.09 2*0.09 2*0.09];
        rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', [0 0 0 0.1]);
        %text(vbarR_top(1), vbarR_top(2)+0.3, 'D');

        pos = [vbarL_bottom-0.09 2*0.09 2*0.09];
        rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', [0 0 0 0.1]);
        %text(vbarL_bottom(1), vbarL_bottom(2)+0.3, 'E');
        pos = [vbarR_bottom-0.09 2*0.09 2*0.09];
        rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', [0 0 0 0.1]);
        %text(vbarR_bottom(1), vbarR_bottom(2)+0.3, 'F');

        exportgraphics(gcf, gifFile, Append=true);
        if i == 1
            pause(1);
        else 
            pause((tt(i)-tt(i-1)) * 400);
        end
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