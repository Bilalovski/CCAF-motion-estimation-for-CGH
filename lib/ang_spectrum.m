function X = ang_spectrum(X,pp,z,wlen)
% applies the angular spectrum propagation method on the wavefield X with
% pixel pitch 'pp', propagation distance 'z' and wavelength 'wlen'.
% maintains pixel pitch size.

qpos = @(N) (wlen/pp/N*(-N/2:N/2-1)).^2; % quadratic position index
X = ifft2(fft2(X) .* ifftshift(exp(2i*pi/wlen*z*sqrt(1 - qpos(size(X,2)) - qpos(size(X,1))'))));

end