function [H, G, S] = load_scene1_pair(data_dir, frame_a, frame_b)
% Load two frames of the scene_1 sequence plus the ground truth for that pair.
%
%   [H, G, S] = load_scene1_pair(data_dir, frame_a, frame_b)
%
%   H, G : frame_a / frame_b, complex field at the hologram plane (N x N, single)
%   S    : pp, lambda, z_obj (m), N, n_frames, frame_a, frame_b
%          names  : 1 x n_obj cellstr (pretzel, pawn, bishop, rook)
%          gt     : ground-truth motion frame_a -> frame_b, one entry per object:
%                   theta_x_deg, theta_y_deg, tx_px, ty_px, tz_mm (1 x n_obj each)
%          label  : object map of frame_a at N/label_ds resolution
%                   (0 = background, j = names{j}); pixel (r,c) of the field is
%                   label(ceil(r/label_ds), ceil(c/label_ds))
%          label_ds
%
% Any frame_a, frame_b in 1..n_frames works, in either order (b < a = reversed motion).

    meta_file = fullfile(data_dir, 'scene1_meta.mat');
    check_file(meta_file);
    M = load(meta_file);
    S.pp = double(M.pp);  S.lambda = double(M.lambda);  S.z_obj = double(M.z_obj);
    S.N = double(M.N);    S.n_frames = double(M.n_frames);
    assert(all(ismember([frame_a frame_b], 1:S.n_frames)), ...
           'frames must be integers in 1..%d (got %g and %g)', S.n_frames, frame_a, frame_b);
    S.frame_a = frame_a;  S.frame_b = frame_b;
    S.names = cellstr(M.gt.names);  S.names = S.names(:)';

    % exact pairwise ground truth, stored as (frame_a, frame_b, object)
    f = {'theta_x_deg','theta_y_deg','tx_px','ty_px','tz_mm'};
    for i = 1:numel(f)
        v = reshape(double(M.gt.(['pair_' f{i}])(frame_a, frame_b, :)), 1, []);
        v(abs(v) < 1e-9) = 0;          % float round-off (e.g. -1e-17) -> prints as 0, not -0
        S.gt.(f{i}) = v;
    end
    S.gt.names = S.names;

    S.label = M.label(:, :, frame_a);
    S.label_ds = double(M.label_ds);

    H = load_frame(data_dir, frame_a);
    G = load_frame(data_dir, frame_b);
end

function U = load_frame(data_dir, k)
    fn = fullfile(data_dir, sprintf('scene1_frame%02d.mat', k));
    check_file(fn);
    F = load(fn, 'H');
    U = F.H;
end

function check_file(fn)
    if ~isfile(fn)
        error('missing %s', fn);
    end
    d = dir(fn);
    if d.bytes < 10e3      % a Git LFS pointer is ~130 bytes of text
        error(['%s is only %d bytes -- probably a Git LFS pointer, not the data. ' ...
               'Run "git lfs pull", or for just some frames e.g. ' ...
               'git lfs pull --include "data/scene1_meta.mat,data/scene1_frame01.mat,data/scene1_frame04.mat"'], ...
              fn, d.bytes);
    end
end
