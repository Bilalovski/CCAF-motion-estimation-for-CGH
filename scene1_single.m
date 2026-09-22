%% scene1_single.m -- CCAF on ONE space-frequency crop of scene_1 (frame 1 -> frame 2)
% Click the centre of an M x M crop on the reconstruction, then (if K < M) the
% centre of a K x K band on that crop's spectrum. The CCAF runs on the crop;
% compare its estimate with the true motion of each object, printed below it.
% Axes: rows = x, columns = y (theta_x / tx act along rows, theta_y / ty along columns).

clear; clc; close all;
this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'lib'));

%% ===================== PARAMETERS =====================
% M : spatial crop size (px, even)
%     - M/2 > the translation: the CCAF finds shifts in [-M/2, M/2)
%       (pretzel |ty| = 100 px -> M >= 256)
%     - M < the object you crop, so the crop holds a single motion
% K : frequency crop size (px, even, K <= M); K = M keeps the whole spectrum
%     - K/2 >= M*pp*sin(theta)/lambda: a tilt theta shifts the crop's spectrum by that many
%       bins, and the band must still overlap between the two frames
%       (M = 256: ~5 bins per 0.1 deg -> the pretzel's 1 deg needs K >= 101)
M = 256;
K = 64;

THETA = -1.2:0.02:1.2;               % deg, search grid for theta_x and theta_y
TZ    = linspace(-5e-3, 5e-3, 21);   % m,   search grid for t_z
%% ======================================================

%% ---- LOADING ----
[H, G, S] = load_scene1_pair(fullfile(this_dir, 'data'), 1, 2);
pp = S.pp;  lambda = S.lambda;  N = S.N;
Hobj = ang_spectrum(H, pp, -S.z_obj, lambda);
Gobj = ang_spectrum(G, pp, -S.z_obj, lambda);

use_gpu = false;
try, use_gpu = gpuDeviceCount > 0; catch, end    % no Parallel Computing Toolbox -> CPU
dev = {'CPU', 'GPU'};
fprintf('CCAF on the %s | M = %d: shifts up to %d px | K = %d: tilts up to %.2f deg\n', ...
        dev{use_gpu+1}, M, M/2, K, asind((K/2)*lambda/(M*pp)));

%% ---- CHOOSE THE CROP ----
a = abs(Hobj);  s = sort(a(:));
figure('Color','w','Position',[60 80 650 650],'Name','crop');
imagesc(a, [0 s(round(0.995*end))]); axis image; colormap gray; hold on;
xlabel('column = y (px)'); ylabel('row = x (px)');
title(sprintf('frame 1: click the centre of your %dx%d crop', M, M));
[c, r] = ginput(1);
r = min(max(round(r), M/2), N-M/2);   c = min(max(round(c), M/2), N-M/2);   % keep the crop inside
rectangle('Position', [c-M/2, r-M/2, M, M], 'EdgeColor', 'g', 'LineWidth', 2);
title(sprintf('frame 1: crop centred at (row %d, col %d)', r, c));  drawnow;
Hs = Hobj(r-M/2+1:r+M/2, c-M/2+1:c+M/2);
Gs = Gobj(r-M/2+1:r+M/2, c-M/2+1:c+M/2);

if K < M    % click the band centre on the crop's spectrum, keep only that band
    Hf = fftshift(fft2(Hs));  Gf = fftshift(fft2(Gs));
    figure('Color','w','Position',[730 80 600 650],'Name','crop spectrum');
    imagesc(log10(abs(Hf))); axis image; colormap(gca, 'hot'); hold on;
    plot(M/2+1, M/2+1, 'c+', 'MarkerSize', 10);
    title(sprintf('crop spectrum (+ = DC): click the centre of the %dx%d band', K, K));
    [fc, fr] = ginput(1);
    fr = min(max(round(fr), K/2), M-K/2);   fc = min(max(round(fc), K/2), M-K/2);
    rectangle('Position', [fc-K/2+0.5, fr-K/2+0.5, K, K], 'EdgeColor', 'c', 'LineWidth', 2);  drawnow;
    band = false(M);  band(fr-K/2+1:fr+K/2, fc-K/2+1:fc+K/2) = true;
    Hc = ifft2(ifftshift(Hf .* band));  Gc = ifft2(ifftshift(Gf .* band));
else        % K = M: the whole spectrum
    fr = M/2+1;  fc = M/2+1;  Hc = Hs;  Gc = Gs;
end

%% ---- CCAF ----
[bp, A] = ccaf_2D_fast_twostage_refined(Hc, Gc, fc, fr, THETA, THETA, TZ, pp, lambda, K, M, use_gpu);

est = [bp.theta_x, bp.theta_y, bp.tx/pp, bp.ty/pp, bp.tz*1e3];
gt  = [S.gt.theta_x_deg; S.gt.theta_y_deg; S.gt.tx_px; S.gt.ty_px; S.gt.tz_mm]';
fprintf('\ntrue motion of each object vs the CCAF estimate for your crop:\n\n');
disp(array2table([gt; est], 'RowNames', [S.names, {'CCAF estimate'}], ...
     'VariableNames', {'theta_x_deg', 'theta_y_deg', 'tx_px', 'ty_px', 'tz_mm'}));

% response landscape (theta_x vs theta_y at the best t_z), every object's true rotation marked
figure('Color','w','Position',[1350 80 560 600],'Name','CCAF landscape');
imagesc(THETA, THETA, A(:,:,bp.i_tz)); axis xy; axis image; colorbar; hold on;
plot(bp.theta_y, bp.theta_x, 'g+', 'MarkerSize', 16, 'LineWidth', 2);
plot(S.gt.theta_y_deg, S.gt.theta_x_deg, 'wo', 'MarkerSize', 10, 'LineWidth', 1.2);
text(S.gt.theta_y_deg + 0.06, S.gt.theta_x_deg, S.names, 'Color', 'w');
xlabel('\theta_y (deg)'); ylabel('\theta_x (deg)');
title(sprintf('CCAF peak (+): (\\theta_x, \\theta_y) = (%+.2f, %+.2f)\\circ   o = true', bp.theta_x, bp.theta_y));
