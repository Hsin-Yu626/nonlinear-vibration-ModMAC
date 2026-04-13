%========================================================================
% DESCRIPTION:
% Investigation of the forced and damped dynamics of a
% 5-DOF oscillator with cubic springs (Liao's symmetric structure).
%
% Method: Shooting method with Floquet stability analysis.
%         Both UP-sweep and DOWN-sweep to capture all folded branches.
%
% Based on NLvib (Krack & Gross, Springer 2019).
%========================================================================
clearvars; 
%close all; 
clc;

addpath('../../SRC');
addpath('../../SRC/MechanicalSystems');

%% ========================================================
%  USER SETTINGS
%  ========================================================
F_amp    = 6;         % excitation amplitude [N]
ds       = 0.0005;    % arc-length step size
eps_tol  = 0.05e-7;    % Newton convergence tolerance
stepmax  = 20000;     % maximum continuation steps
label    = '6N';

%% ========================================================
%  SECTION 1: Define the physical system
%  ========================================================
mi  = [150 5 5 5 5];
ki1 = [0 1.5e5 1.5e5 0 0 0];
di  = 0*ki1;
%gamma = 9e10;
%linear
gamma = 0;

nonlinear_elements = cell(1,2);
nonlinear_elements{1} = struct('type','cubicSpring',...
    'stiffness', gamma, 'force_direction', [1;-1;0;0;0]);
nonlinear_elements{2} = struct('type','cubicSpring',...
    'stiffness', gamma, 'force_direction', [1;0;0;-1;0]);

Fex1 = [0; 0; 0; 0; F_amp];
oscillator = ChainOfOscillators(mi, di, ki1, nonlinear_elements, Fex1);

oscillator.K([1 4 5],[1 4 5]) = oscillator.K([1 4 5],[1 4 5]) ...
                               + oscillator.K([1 2 3],[1 2 3]);
oscillator.K(1,1) = oscillator.K(1,1) + 1.5e7;
oscillator.D = 0.6*oscillator.M + 2.5e-7*oscillator.K;

n = oscillator.n;

%% ========================================================
%  SECTION 2: Linear natural frequencies
%  ========================================================
[~, D_modes] = eig(oscillator.K, oscillator.M);
freq_n_all = sort(sqrt(diag(D_modes)) / (2*pi));

fprintf('\n Linear Natural Frequencies (Hz):\n');
for kk = 1:length(freq_n_all)
    fprintf('  Mode %d:  %.4f Hz\n', kk, freq_n_all(kk));
end
fprintf('\n');

%% ========================================================
%  SECTION 3: Shooting parameters & sweep range
%  ========================================================
analysis = 'FRF';
Ntd = 2^8;
Np  = 1;
H   = 7;

Om_s = 42*2*pi;
Om_e = 49*2*pi;

fprintf('Sweep range: %.1f -- %.1f Hz\n', Om_s/2/pi, Om_e/2/pi);
fprintf('ds = %.4f,  eps = %.1e,  stepmax = %d\n\n', ds, eps_tol, stepmax);

%% ========================================================
%  SECTION 4a: UP-SWEEP (42 Hz -> 49 Hz)
%  ========================================================
Q1_up = (-Om_s^2*oscillator.M + 1i*Om_s*oscillator.D ...
         + oscillator.K) \ Fex1;
ys_up = [real(Q1_up); -Om_s*imag(Q1_up)];

Sopt_up = struct(...
    'Dscale',              [1e0*ones(size(ys_up)); Om_s], ...
    'dynamicDscale',       1, ...
    'eps',                 eps_tol, ...
    'stepmax',             stepmax, ...
    'reversaltolerance',   inf, ...
    'errmax',              5, ...
    'dsmin',               ds/50, ...
    'dsmax',               ds*10, ...
    'noconv_stepcor',      'redinc');

fprintf('Starting UP-sweep (%s)...\n', label);
tic;
[X_up, ~, ~] = solve_and_continue(ys_up, ...
    @(X) shooting_residual(X, oscillator, Ntd, Np, analysis), ...
    Om_s, Om_e, ds, Sopt_up);
fprintf('  UP-sweep done in %.1f s.  %d points.\n\n', toc, size(X_up,2));

%% ========================================================
%  SECTION 4b: DOWN-SWEEP (49 Hz -> 42 Hz)
%  ========================================================
Q1_dn = (-Om_e^2*oscillator.M + 1i*Om_e*oscillator.D ...
         + oscillator.K) \ Fex1;
ys_dn = [real(Q1_dn); -Om_e*imag(Q1_dn)];

Sopt_dn = struct(...
    'Dscale',              [1e0*ones(size(ys_dn)); Om_e], ...
    'dynamicDscale',       1, ...
    'eps',                 eps_tol, ...
    'stepmax',             stepmax, ...
    'reversaltolerance',   inf, ...
    'errmax',              5, ...
    'dsmin',               ds/50, ...
    'dsmax',               ds*10, ...
    'noconv_stepcor',      'redinc');

fprintf('Starting DOWN-sweep (%s)...\n', label);
tic;
[X_dn, ~, ~] = solve_and_continue(ys_dn, ...
    @(X) shooting_residual(X, oscillator, Ntd, Np, analysis), ...
    Om_e, Om_s, ds, Sopt_dn);
fprintf('  DOWN-sweep done in %.1f s.  %d points.\n\n', toc, size(X_dn,2));

%% ========================================================
%  SECTION 5: Floquet stability (both sweeps)
%  ========================================================
Om_up_hz = X_up(end,:);
Om_dn_hz = X_dn(end,:);

fprintf('Running Floquet analysis on UP-sweep...\n');
tic;
[a_rms_up, stable_up, ~, BPs_up] = ...
    floquet_analysis(X_up, oscillator, Ntd, Np, analysis, H, n);
fprintf('  UP Floquet done in %.1f s.\n', toc);

fprintf('Running Floquet analysis on DOWN-sweep...\n');
tic;
[a_rms_dn, stable_dn, ~, BPs_dn] = ...
    floquet_analysis(X_dn, oscillator, Ntd, Np, analysis, H, n);
fprintf('  DOWN Floquet done in %.1f s.\n\n', toc);

fprintf('--- Bifurcation Points (UP-sweep) ---\n');
for i = 1:length(BPs_up)
    fprintf('  %s at %.4f Hz,  a_rms = %.3e\n', ...
        BPs_up(i).type, BPs_up(i).Om/2/pi, BPs_up(i).a_rms);
end
fprintf('\n--- Bifurcation Points (DOWN-sweep) ---\n');
for i = 1:length(BPs_dn)
    fprintf('  %s at %.4f Hz,  a_rms = %.3e\n', ...
        BPs_dn(i).type, BPs_dn(i).Om/2/pi, BPs_dn(i).a_rms);
end
fprintf('\n');

%% ========================================================
%  SECTION 6: Plot results (both sweeps combined)
%  ========================================================
figure('Name', ['FRF Shooting — ' label], 'Color', 'w');
hold on; box on; grid on;

% UP stable
hS = plot(Om_up_hz(logical(stable_up))/2/pi, ...
          a_rms_up(logical(stable_up)), ...
          'k-', 'LineWidth', 1.5, 'DisplayName', 'Shooting (stable)');
% DN stable
plot(Om_dn_hz(logical(stable_dn))/2/pi, ...
     a_rms_dn(logical(stable_dn)), ...
     'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% UP unstable
hU = plot(Om_up_hz(~logical(stable_up))/2/pi, ...
          a_rms_up(~logical(stable_up)), ...
          'r--', 'LineWidth', 1.0, 'DisplayName', 'Unstable');
% DN unstable
plot(Om_dn_hz(~logical(stable_dn))/2/pi, ...
     a_rms_dn(~logical(stable_dn)), ...
     'r--', 'LineWidth', 1.0, 'HandleVisibility', 'off');

% Bifurcation markers
BPs_all = [BPs_up, BPs_dn];
plot_bifurcation_points(BPs_all);

hTP = plot(nan, nan, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', ...
    'DisplayName', 'TP');
hNS = plot(nan, nan, 'c^', 'MarkerSize', 6, 'MarkerFaceColor', 'c', ...
    'DisplayName', 'NS');
hPD = plot(nan, nan, 'ms', 'MarkerSize', 6, 'MarkerFaceColor', 'm', ...
    'DisplayName', 'PD');

xlabel('Excitation Frequency (Hz)', 'FontSize', 12);
ylabel('Response Amplitude (RMS)',  'FontSize', 12);
title(['Single point, ' label], 'FontSize', 14, 'FontWeight', 'bold');
legend([hS, hU, hTP, hNS, hPD], 'Location', 'northwest', 'FontSize', 9);
set(gca, 'XLim', sort([Om_s Om_e]/2/pi), 'FontSize', 11);

%% ========================================================
%  SECTION 7: Save data
%  ========================================================
fname = ['data_' label '_' datestr(now,30) '.mat'];
save(fname, 'oscillator', 'ds', 'Om_s', 'Om_e', ...
     'X_up', 'a_rms_up', 'stable_up', 'BPs_up', ...
     'X_dn', 'a_rms_dn', 'stable_dn', 'BPs_dn');
fprintf('Data saved to %s\n', fname);

%% ========================================================
%  LOCAL FUNCTIONS
%  ========================================================

function [a_rms, stable, mucrit, BPs] = ...
    floquet_analysis(X_shoot, oscillator, Ntd, Np, analysis, H, n)

    Npts   = size(X_shoot, 2);
    a_rms  = zeros(1, Npts);
    stable = zeros(1, Npts);
    mucrit = zeros(1, Npts);
    BPs(1:Npts) = struct('Om',[],'a_rms',[],'type','');
    nBP = 0;

    for i = 1:Npts
        [~,~,~, Y, dye_dys] = shooting_residual(X_shoot(:,i), ...
            oscillator, Ntd, Np, analysis);
        Qc        = fft(Y(:,n)) / Ntd;
        a_rms(i)  = sqrt(sum(abs(Qc(1:H+1)).^2)) / sqrt(2) * 2;
        mucrit(i) = eigs(dye_dys, 1, 'lm');
        stable(i) = (abs(mucrit(i)) <= 1);

        if i > 1 && stable(i) ~= stable(i-1)
            nBP = nBP + 1;
            BPs(nBP).Om    = (X_shoot(end,i) + X_shoot(end,i-1)) / 2;
            BPs(nBP).a_rms = (a_rms(i) + a_rms(i-1)) / 2;
            theta = angle(mucrit(i));
            if abs(theta) < 1e-3
                BPs(nBP).type = 'TP';
            elseif abs(abs(theta) - pi) < 1e-3
                BPs(nBP).type = 'PD';
            else
                BPs(nBP).type = 'NS';
            end
        end
    end
    BPs = BPs(1:nBP);
end

function plot_bifurcation_points(BPs)
    iii = find(strcmp({BPs.type}, 'TP'));
    if ~isempty(iii)
        plot([BPs(iii).Om]/2/pi, [BPs(iii).a_rms], ...
            'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    end
    iii = find(strcmp({BPs.type}, 'NS'));
    if ~isempty(iii)
        plot([BPs(iii).Om]/2/pi, [BPs(iii).a_rms], ...
            'c^', 'MarkerSize', 6, 'MarkerFaceColor', 'c');
    end
    iii = find(strcmp({BPs.type}, 'PD'));
    if ~isempty(iii)
        plot([BPs(iii).Om]/2/pi, [BPs(iii).a_rms], ...
            'ms', 'MarkerSize', 6, 'MarkerFaceColor', 'm');
    end
end