function stimulus_fft(type, wl, ratio)
% STIMULUS_FFT uses fast fourier transform to plot the amplitude spectrum of a potential stimulus pattern
%
% Inputs:
%   type   - type of stimulus, e.g. 'bar'/'dog'/'square'/'log' (see below for full list)
%   wl     - wavelength of stimulus in degrees (definition depends on stimulus stype)
%   ratio  - ratio between black and white areas; usually ratio of amplitudes, but definition depends on stimulus type
%
% example usage: stimulus_fft('dog', 10, 1)
%
% Written by Jochen Smolka and John D. Kirwan, Lund University, 2017
%
% This file is published as supplementary to the following article:
%   John D. Kirwan, Michael J. Bok, Jochen Smolka, James J. Foster, José Carlos Hernández, Dan-Eric Nilsson (2017) 
%   The sea urchin Diadema africanum uses low resolution vision to find shelter and deter enemies. Journal of Experimental Biology.

if nargin<3, ratio = 1; end
if nargin<2, wl = 10; end

va  = -180:00.1:179.99;
Fs  = length(va);         % Simulation resolution (points per 360 degrees)
azi = (0:Fs-1)/Fs * 360 - 180;      % azimuth vector
S   = stimulus_create(type, wl, va, ratio);

figure(1);
clf;
plot(azi,S)
title(sprintf('%s, wl=%d', type, wl))
xlabel('azimuth (\circ)')
ylabel('intensity')
axis([-180 180 -1.1 1.1])

Y           = fft(S);
P2          = abs(Y/Fs);
P1          = P2(1:Fs/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f           = Fs*(0:(Fs/2))/Fs/360;

figure(2);
clf;
plot(f, P1)
title(sprintf('%s, wl=%d', type, wl))
xlabel('f (cycles / deg)')
ylabel('|P1(f)|')
xlim([0 1])