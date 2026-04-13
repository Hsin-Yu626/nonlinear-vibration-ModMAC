%========================================================================
% ModMAC Post-Processing for Shooting / HBM Results
%
% ModMAC(j,i) = |Phi_B'*M*Phi_E|^2 / |(Phi_E'*M*Phi_E)*(Phi_B'*M*Phi_B)|
%
% KEY FIX (advisor):  Phi_B must be the COMPLEX ODS, not just the real part.
%   For shooting: Phi_B = q(0) + 1i*qdot(0)   [rows 1:n + 1i*rows n+1:2n]
%   For HBM:     Phi_B = Qcos  + 1i*Qsin       [rows 1:n + 1i*rows n+1:2n]
%   ModMAC is scale-invariant, so the omega factor in qdot doesn't matter.
%
% Figure 1: FRF + ModMAC heatmap  (colour = dominant mode, solid/dashed = stable/unstable)
% Figure 2: ModMAC curves         (solid/dashed = stable/unstable)
%
% Requires in workspace: oscillator, n
%   and at least one of:  X_up (+stable_up), X_dn (+stable_dn), X_shoot
%   optional:             a_rms_up, a_rms_dn, BPs_up, BPs_dn
%========================================================================

%% ========================================================
%  STEP 1: Linear mode shapes
%  ========================================================

[V_modes, D_modes] = eig(full(oscillator.K), full(oscillator.M));
[omega2_sorted, sort_idx] = sort(diag(D_modes));
V = V_modes(:, sort_idx);
for im = 1:n
    V(:,im) = V(:,im) / norm(V(:,im));
end

freq_n = sqrt(omega2_sorted) / (2*pi);
fprintf('\n=== Linear Natural Frequencies ===\n');
for im = 1:n
    fprintf('  Mode %d: %.4f Hz\n', im, freq_n(im));
end

%% ========================================================
%  STEP 2: Build sweep structs -- UP and DOWN kept SEPARATE
%  ========================================================

has_up     = exist('X_up',     'var') && exist('stable_up', 'var');
has_dn     = exist('X_dn',     'var') && exist('stable_dn', 'var');
has_rms_up = exist('a_rms_up', 'var');
has_rms_dn = exist('a_rms_dn', 'var');

if ~has_up && ~exist('X_shoot','var')
    error('No shooting results found. Run LiaoCase_ShootingV2 first.');
end

sweeps = struct('freq',{},'arms',{},'stable',{});

if has_up
    s.freq   = X_up(end,:) / (2*pi);
    s.stable = logical(stable_up(:)');
    if has_rms_up
        s.arms = a_rms_up(:)';
    else
        % Complex amplitude RMS for proxy
        Phi_complex = X_up(1:n,:) + 1i*X_up(n+1:2*n,:);
        s.arms = sqrt(mean(abs(Phi_complex).^2, 1));
        fprintf('  [Warning] a_rms_up missing -- using complex-amplitude RMS as proxy.\n');
    end
    sweeps(end+1) = s;
elseif exist('X_shoot','var')
    s.freq   = X_shoot(end,:) / (2*pi);
    Phi_complex = X_shoot(1:n,:) + 1i*X_shoot(n+1:2*n,:);
    s.arms   = sqrt(mean(abs(Phi_complex).^2, 1));
    s.stable = true(1, size(X_shoot,2));
    sweeps(end+1) = s;
end

if has_dn
    s.freq   = X_dn(end,:) / (2*pi);
    s.stable = logical(stable_dn(:)');
    if has_rms_dn
        s.arms = a_rms_dn(:)';
    else
        Phi_complex = X_dn(1:n,:) + 1i*X_dn(n+1:2*n,:);
        s.arms = sqrt(mean(abs(Phi_complex).^2, 1));
        fprintf('  [Warning] a_rms_dn missing -- using complex-amplitude RMS as proxy.\n');
    end
    sweeps(end+1) = s;
end

fprintf('  Loaded %d sweep(s).\n', numel(sweeps));

for si = 1:numel(sweeps)
    sweeps(si).stable = logical(sweeps(si).stable);
    n_stable   = sum( sweeps(si).stable);
    n_unstable = sum(~sweeps(si).stable);
    fprintf('  Sweep %d: %d stable, %d unstable points\n', si, n_stable, n_unstable);
end

%% ========================================================
%  STEP 3: Compute ModMAC for every point in every sweep
%          ** Using COMPLEX ODS:  Phi_B = X(1:n) + 1i*X(n+1:2n) **
%  ========================================================

M_mat = full(oscillator.M);

for si = 1:numel(sweeps)
    sw = sweeps(si);
    N  = numel(sw.freq);

    % Recover ODS from the matching X matrix
    if si == 1 && has_up
        X_src = X_up;
    elseif si == 1 && ~has_up
        X_src = X_shoot;
    else
        X_src = X_dn;
    end

    MM = zeros(N, n);
    for i_mode = 1:n
        Phi_E = V(:, i_mode);   % real linear mode shape
        for k = 1:N
            %--- KEY FIX: complex ODS (advisor's instruction) ---
            Phi_B = X_src(1:n, k) + 1i * X_src(n+1:2*n, k);

            % Hermitian inner products (MATLAB ' is conjugate transpose)
            num = abs(Phi_B' * M_mat * Phi_E)^2;
            den = abs((Phi_E' * M_mat * Phi_E) * (Phi_B' * M_mat * Phi_B));
            MM(k, i_mode) = num / den;
        end
    end

    sweeps(si).ModMAC = MM;
    [~, sweeps(si).dom] = max(MM, [], 2);  % dominant mode index (Nx1)
end

%% ========================================================
%  STEP 4: Colour palette
%  ========================================================

mode_colors = [0.00 0.45 0.74;    % Mode 1 - blue
               0.85 0.33 0.10;    % Mode 2 - red-orange
               0.93 0.69 0.13;    % Mode 3 - gold
               0.49 0.18 0.56;    % Mode 4 - purple
               0.47 0.67 0.19];   % Mode 5 - green

%% ========================================================
%  STEP 5: Figure 1 -- FRF coloured by dominant ModMAC mode
%          + ModMAC heatmap (tiled layout)
%  ========================================================

set(0,'DefaultFigureWindowStyle','docked');

figure(1); clf;
set(gcf,'Name','FRF + ModMAC heatmap','Color','w');
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl);
hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');

for si = 1:numel(sweeps)
    sw = sweeps(si);
    N  = numel(sw.freq);
    seg_start = 1;
    for k = 2:N+1
        last = (k == N+1);
        if last || (sw.dom(k) ~= sw.dom(k-1)) || (sw.stable(k) ~= sw.stable(k-1))
            seg = seg_start : k-1;
            col = mode_colors(sw.dom(seg_start), :);
            if sw.stable(seg_start), ls = '-';  lw = 1.8;
            else,                    ls = '--'; lw = 1.2;
            end
            plot(ax1, sw.freq(seg), sw.arms(seg), ls, ...
                 'Color', col, 'LineWidth', lw, 'HandleVisibility','off');
            seg_start = k;
        end
    end
end

% Bifurcation markers
BPs_combined = [];
if exist('BPs_up','var'), BPs_combined = [BPs_combined, BPs_up]; end
if exist('BPs_dn','var'), BPs_combined = [BPs_combined, BPs_dn]; end
if ~isempty(BPs_combined)
    plot_bps(ax1, BPs_combined);
end

% Legend
h_leg = gobjects(n+2,1);
for im = 1:n
    h_leg(im) = patch(ax1, nan, nan, mode_colors(im,:), ...
        'EdgeColor', mode_colors(im,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Mode %d', im));
end
h_leg(n+1) = plot(ax1, nan, nan, 'k-',  'LineWidth', 1.8, 'DisplayName', 'Stable');
h_leg(n+2) = plot(ax1, nan, nan, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Unstable');
legend(ax1, h_leg, 'Location','northwest', 'FontSize', 9, 'NumColumns', 2);

xlabel(ax1, 'Excitation Frequency (Hz)', 'FontSize', 11);
ylabel(ax1, 'Response Amplitude (RMS)',  'FontSize', 11);
title(ax1, 'FRF  (colour: dominant ModMAC mode  |  line: stability)', ...
     'FontSize', 11, 'FontWeight', 'bold');
set(ax1, 'FontSize', 10, 'TickDir', 'out');
xlim(ax1, [Om_s Om_e] / (2*pi));

%% ========================================================
%  STEP 6: ModMAC heatmap (first sweep, inside Figure 1)
%  ========================================================

sw1 = sweeps(1);

ax2 = nexttile(tl);

[freq_sorted, sort_idx] = sort(sw1.freq);
ModMAC_sorted = sw1.ModMAC(sort_idx, :);

freq_uniform = linspace(Om_s/(2*pi), Om_e/(2*pi), 1000);
ModMAC_uniform = zeros(1000, n);
for im = 1:n
    ModMAC_uniform(:,im) = interp1(freq_sorted, ModMAC_sorted(:,im), ...
        freq_uniform, 'linear', 'extrap');
end
ModMAC_uniform = max(0, min(1, ModMAC_uniform));

imagesc(ax2, freq_uniform, 1:n, ModMAC_uniform');
set(ax2,'YDir','normal');
colormap(ax2, flip(gray));
cb = colorbar(ax2); cb.Label.String = 'ModMAC';
clim(ax2, [0 1]);
xlim(ax2, [Om_s Om_e] / (2*pi));
xlabel(ax2, 'Excitation frequency (Hz)', 'FontSize', 11);
ylabel(ax2, 'Mode', 'FontSize', 11);
yticks(ax2, 1:n);
set(ax2, 'FontSize', 10, 'TickDir', 'out');

%% ========================================================
%  STEP 7: Figure 2 -- ModMAC curves, solid/dashed = stable/unstable
%  ========================================================

figure(2); clf;
set(gcf,'Name','ModMAC curves','Color','w');
ax3 = gca;
hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');

for si = 1:numel(sweeps)
    sw = sweeps(si);
    for im = 1:n
        col = mode_colors(im,:);
        seg_plot(ax3, sw.freq, sw.ModMAC(:,im)', sw.stable,  col, '-',  1.4);
        seg_plot(ax3, sw.freq, sw.ModMAC(:,im)', ~sw.stable, col, '--', 1.0);
    end
end

h_dummy = gobjects(n+1,1);
for im = 1:n
    h_dummy(im) = plot(ax3, nan, nan, '-', ...
        'Color', mode_colors(im,:), 'LineWidth', 1.4, ...
        'DisplayName', sprintf('Mode %d', im));
end
h_dummy(n+1) = plot(ax3, nan, nan, 'k--', 'LineWidth', 1.0, 'DisplayName', 'Unstable');
legend(ax3, h_dummy, 'Location','best', 'FontSize', 9);

xlabel(ax3, 'Excitation Frequency (Hz)', 'FontSize', 11);
ylabel(ax3, 'ModMAC', 'FontSize', 11);
ylim(ax3, [0 1.05]);
xlim(ax3, [Om_s Om_e] / (2*pi));
set(ax3, 'FontSize', 10, 'TickDir', 'out');

%% ========================================================
%  LOCAL FUNCTIONS
%  ========================================================

function seg_plot(ax, x, y, mask, color, ls, lw)
    idx = find(mask);
    if isempty(idx), return; end
    breaks = [0, find(diff(idx) > 1), numel(idx)];
    for k = 1:numel(breaks)-1
        seg = idx(breaks(k)+1 : breaks(k+1));
        plot(ax, x(seg), y(seg), ls, ...
            'Color', color, 'LineWidth', lw, 'HandleVisibility', 'off');
    end
end

function plot_bps(ax, BPs)
    types  = {'TP',  'NS',  'PD'};
    marker = {'ro',  'c^',  'ms'};
    mfc    = {'r',   'c',   'm' };
    for t = 1:3
        idx = strcmp({BPs.type}, types{t});
        if any(idx)
            plot(ax, [BPs(idx).Om]/2/pi, [BPs(idx).a_rms], ...
                marker{t}, 'MarkerSize', 6, 'MarkerFaceColor', mfc{t}, ...
                'HandleVisibility','off');
        end
    end
end