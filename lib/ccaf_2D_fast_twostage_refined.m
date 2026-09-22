%% ========== CCAF FUNCTION (GPU or CPU) ==========
function [best_params, A_peaks] = ccaf_2D_fast_twostage_refined(H_crop, G_crop, cf_x, cf_y, ...
    theta_x_range, theta_y_range, tz_range, pp, lambda, K, M, use_gpu)
% Two-stage Cross-Chirp Ambiguity Function motion estimator.
%
% Runs on a CUDA GPU when one is available (Parallel Computing Toolbox), and
% otherwise falls back to plain CPU arrays so it still works on a laptop with
% no GPU. The maths is identical either way -- only the array backend changes.
%
%   use_gpu (optional): [] / omitted -> auto-detect; true -> force GPU;
%                       false -> force CPU. Auto-detect never errors if the
%                       Parallel Computing Toolbox is missing.
%
% NOTE on CPU: the search is heavy (n_theta_x * n_theta_y * n_tz FFTs of size
% M and K). On a CPU this can take minutes. If it is too slow, shrink the
% search ranges (fewer points) or M/K -- accuracy of the coarse stage is
% governed by the grid spacing, and stage 2 refines tz around the peak.

    % ---------- pick the array backend (GPU if present, else CPU) ----------
    if nargin < 12 || isempty(use_gpu)
        use_gpu = false;
        try
            use_gpu = gpuDeviceCount > 0;   % errors if no PCT -> caught -> CPU
        catch
            use_gpu = false;
        end
    end
    if use_gpu
        toG    = @(x) gpuArray(x);
        zerosG = @(varargin) zeros(varargin{:}, 'gpuArray');
    else
        toG    = @(x) x;
        zerosG = @(varargin) zeros(varargin{:});
        fprintf(['[ccaf] no GPU -> running on CPU (this can take a while; ' ...
                 'reduce the theta/tz ranges if needed).\n']);
    end

    n_theta_x = length(theta_x_range);
    n_theta_y = length(theta_y_range);
    n_tz = length(tz_range);
    half_K = K / 2;

    alpha_range = -lambda * tz_range / (M * pp^2);
    nu_x_all = -M * pp * deg2rad(-theta_x_range) / lambda;
    nu_y_all = M * pp * deg2rad(-theta_y_range) / lambda;

    crop_rows_shifted = (cf_y - half_K + 1):(cf_y + half_K);
    crop_cols_shifted = (cf_x - half_K + 1):(cf_x + half_K);
    crop_rows_unshifted = mod(crop_rows_shifted - 1 - M/2, M) + 1;
    crop_cols_unshifted = mod(crop_cols_shifted - 1 - M/2, M) + 1;

    [j_grid, i_grid] = meshgrid(1:K, 1:K);
    rows_unsh = crop_rows_unshifted(i_grid);
    cols_unsh = crop_cols_unshifted(j_grid);
    linear_idx = sub2ind([M, M], rows_unsh, cols_unsh);
    linear_idx_gpu = toG(linear_idx(:));

    H_gpu = toG(H_crop);
    G_gpu = toG(G_crop);
    G_fft = fft2(G_gpu);
    G_fft_K = G_fft(linear_idx_gpu);

    wx_full = (-M/2:M/2-1)';
    wy_full = (-M/2:M/2-1);

    QE_all_K = zerosG(K*K, n_tz);
    for i_tz = 1:n_tz
        alpha = alpha_range(i_tz);
        QE_x = exp(1i * pi / M * alpha * wx_full.^2);
        QE_y = exp(1i * pi / M * alpha * wy_full.^2);
        QE_full = fftshift(QE_x * QE_y);
        QE_all_K(:, i_tz) = QE_full(linear_idx(:));
    end

    nx = toG((0:M-1)');
    ny = toG((0:M-1));

    shift_y_all = zerosG(1, M, n_theta_y);
    for i_thy = 1:n_theta_y
        shift_y_all(1, :, i_thy) = exp(2i * pi * nu_y_all(i_thy) * ny / M);
    end

    A_peaks_gpu = zerosG(n_theta_x, n_theta_y, n_tz);

    for i_thx = 1:n_theta_x
        shift_x = exp(2i * pi * nu_x_all(i_thx) * nx / M);
        H_shifted_x = H_gpu .* shift_x;
        H_shifted_all = H_shifted_x .* shift_y_all;

        H_fft_all = fft2(H_shifted_all);

        H_fft_K_all = zerosG(K*K, n_theta_y);
        for i_thy = 1:n_theta_y
            H_fft_slice = H_fft_all(:,:,i_thy);
            H_fft_K_all(:, i_thy) = H_fft_slice(linear_idx_gpu);
        end

        cross_spec_K_all = H_fft_K_all .* conj(G_fft_K);

        for i_tz = 1:n_tz
            QE_K = QE_all_K(:, i_tz);
            spectrum_K_all = cross_spec_K_all .* QE_K;
            spectrum_K_3d = reshape(spectrum_K_all, K, K, n_theta_y);
            spatial_K_all = ifft2(spectrum_K_3d);
            A_K_all = abs(spatial_K_all).^2;
            A_reshaped = reshape(A_K_all, K*K, n_theta_y);
            max_vals = max(A_reshaped, [], 1);
            A_peaks_gpu(i_thx, :, i_tz) = max_vals;
        end
    end

    A_peaks = gather(A_peaks_gpu);
    [~, max_ind] = max(A_peaks(:));
    [i_thx_opt, i_thy_opt, i_tz_opt] = ind2sub(size(A_peaks), max_ind);

    tx_range = pp * (-M/2:M/2-1);

    tz_window = 3;
    n_tz_fine = 31;

    i_tz_min = max(1, i_tz_opt - tz_window);
    i_tz_max = min(n_tz, i_tz_opt + tz_window);
    tz_min = tz_range(i_tz_min);
    tz_max = tz_range(i_tz_max);

    tz_fine = linspace(tz_min, tz_max, n_tz_fine);
    alpha_fine = -lambda * tz_fine / (M * pp^2);

    QE_fine_K = zerosG(K*K, n_tz_fine);
    for i_tz = 1:n_tz_fine
        alpha = alpha_fine(i_tz);
        QE_x = exp(1i * pi / M * alpha * wx_full.^2);
        QE_y = exp(1i * pi / M * alpha * wy_full.^2);
        QE_full = fftshift(QE_x * QE_y);
        QE_fine_K(:, i_tz) = QE_full(linear_idx(:));
    end

    nu_x_opt = nu_x_all(i_thx_opt);
    nu_y_opt = nu_y_all(i_thy_opt);

    shift_x_opt = exp(2i * pi * nu_x_opt * nx / M);
    shift_y_opt = exp(2i * pi * nu_y_opt * ny / M);

    H_shifted_opt = H_gpu .* (shift_x_opt .* shift_y_opt);
    H_fft_opt = fft2(H_shifted_opt);
    cross_spec_opt = H_fft_opt .* conj(G_fft);
    cross_spec_K_opt = cross_spec_opt(linear_idx_gpu);

    best_peak = 0;
    best_tz = tz_range(i_tz_opt);
    best_tx = 0;
    best_ty = 0;

    for i_tz = 1:n_tz_fine
        QE_K = QE_fine_K(:, i_tz);
        spectrum_K = cross_spec_K_opt .* QE_K;

        spectrum_M = zerosG(M, M);
        spectrum_M(linear_idx_gpu) = spectrum_K;

        spatial_M = ifft2(spectrum_M);
        A_M = abs(spatial_M).^2;

        [max_val, max_ind_local] = max(A_M(:));

        if max_val > best_peak
            best_peak = max_val;
            best_tz = tz_fine(i_tz);

            [peak_row, peak_col] = ind2sub([M, M], max_ind_local);
            i_tx = mod(peak_row - 1 + M/2, M) + 1;
            i_ty = mod(peak_col - 1 + M/2, M) + 1;
            best_tx = -tx_range(i_tx);
            best_ty = -tx_range(i_ty);
        end
    end

    best_params.theta_x = theta_x_range(i_thx_opt);
    best_params.theta_y = theta_y_range(i_thy_opt);
    best_params.tz = best_tz;
    best_params.tx = best_tx;
    best_params.ty = best_ty;
    best_params.peak_val = gather(best_peak);
    best_params.i_thx = i_thx_opt;
    best_params.i_thy = i_thy_opt;
    best_params.i_tz = i_tz_opt;

    if use_gpu
        reset(gpuDevice);
    end
end
