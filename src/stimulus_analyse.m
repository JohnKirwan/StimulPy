function varout = stimulus_analyse(pattern_x, pattern_y, pw, filtertype, plotit)
% STIMULUS_ANALYSE calculates statistics for the ionput stimulus pattern.
%   varout = stimulus_analyse(pattern_x, pattern_y, pw, filtertype, plotit)
%
% Inputs:
%   pattern_x  - x (angular) values describing the pattern; in degrees
%   pattern_x  - y (intensity) values describing the pattern; in any units
%   pw         - pattern width (width of white and black bars; this is only used to determine reasonable axis limits for plotting)
%   filtertype - string describing the acceptance angle filter function (e.f. 'g' for gaussian acceptance functions)
%   plotit     - logical, whether to plot results or not; also determines what output is returned
%
% Outputs:
%   varout     - if plotit is true:  axes handle of all subplots
%                if plotit is false: vector of cut-off values of spatial resolution to get 10%/5%/2.5% image contrast
%
% Example usage: varout = stimulus_analyse(va, int, 3*10, 'g', false)
%
% Written by Jochen Smolka and John D. Kirwan, Lund University, 2017
%
% This file is published as supplementary to the following article:
%   John D. Kirwan, Michael J. Bok, Jochen Smolka, James J. Foster, José Carlos Hernández, Dan-Eric Nilsson (2017) 
%   The sea urchin Diadema africanum uses low resolution vision to find shelter and deter enemies. Journal of Experimental Biology.

%% defaults
if nargin<5, plotit = false; end
if nargin<4, filtertype = 'g'; end

%% parameters
praa_acc = 0.1;                     % accuracy with which PRAAs (photoreceptor acceptance angles) are sampled; 
                                    % the sampling here determines the accuracy of cut-off resolutions. 
                                    % Set to 1 or 2 for quick tests, 0.1 for final smooth graphs
praas    = 1:praa_acc:180;          % PRAAs to test (photoreceptor acceptance angles)
acc      = median(diff(pattern_x)); % the accuracy with which the input pattern was sampled in degrees/pixel, same as acc in stimulus_main

%% create filter (not normalised)
switch lower(filtertype)
    case {'g', 'gauss', 'gaussian'}
        % Gaussian-shaped acceptance function: praas are interpreted as FWHMs (full widths at half-maximum) of the Gaussian
        g = @(x, hw) exp( -(x.^2)/(2*(hw/2.35482)^2) );
    case {'s', 'square'}
        % Square-wave-shaped acceptance function: praas are interpreted as the full width of the square acceptance function
        g = @(x, hw) double(abs(x)<hw/2);
    case {'sin', 'cos'}
        % Sinusoidal acceptance function: praas are interpreted as the full width between zero-crossings of the cosine
        g = @(x, hw) cosd(x./hw*180) .* double(abs(x)<hw/2);
    otherwise
        error('Unknown filtertype: %s', filterpara.type)
end

%% apply filter 
int_filt     = nan(length(pattern_y), length(praas)); % preallocate matrix for filtered intensity patterns
int_amps     = nan(1, length(praas));                 % preallocate matrix for amplitudes of filtered patterns (max-min)
int_mod      = nan(1, length(praas));                 % preallocate matrix for modulation of filtered patterns (max-min)/(max+min)
int_distmins = nan(1, length(praas));                 % preallocate matrix for distance between minimums (currently unused)
int_distmaxs = nan(1, length(praas));                 % preallocate matrix for distance between maximums

for i = 1:length(praas)
    fkern   = g(-300:acc:300, praas(i));    % calculate filter kernel
    fkern   = fkern / sum(fkern);           % normalise filter kernel

    % calculate filtered image, and some statistics
    [int_filt(:, i), int_amps(i), int_distmins(i), int_distmaxs(i), int_mod(i)] = stimulus_analyse_one(pattern_x, pattern_y, fkern, round(1/4*length(pattern_x)):round(6/8*length(pattern_x)));
end

pattern_mod = (max(pattern_y) - min(pattern_y)) / (max(pattern_y) + min(pattern_y)); % modulation (Michelson contrast) of the unfiltered pattern
int_relamps = 100 * int_mod / pattern_mod; % calculate remaining modulation (relative to unfiltered modulation)

%% find the cut-off resolutions (acceptance angles where contrast drops below 20/10/5/1%)
i1 = praas(find(int_relamps<20, 1, 'first'));
if isempty(i1), i1 = NaN; end
i2 = praas(find(int_relamps<10, 1, 'first'));
if isempty(i2), i2 = NaN; end
i3 = praas(find(int_relamps<5, 1, 'first'));
if isempty(i3), i3 = NaN; end
i4 = praas(find(int_relamps<1, 1, 'first'));
if isempty(i4), i4 = NaN; end

if ~plotit % if no plotting required, return a vector of cut-off resolutions
    varout = [i1 i2 i3 i4];
    return;
else % otherwise create a figure and return the subplot handles
    h = sub_formatA4();
    varout = h;
end

%% plot subplot 1: intensity profiles for all filtering steps
axes(h{1});
hh = plot(pattern_x, pattern_y, 'k', 'linewidth', 2);   % plot unfiltered pattern
plot(pattern_x, int_filt(:, mod(praas, 1)==0), 'color', [.7 .7 .7]);         % plot every 1 degree in grey
plot(pattern_x, int_filt(:, mod(praas, 5)==0), 'k');    % plot every 5 degrees in black
uistack(hh, 'top'); % bring unfiltered pattern to the top
xlabel('visual angle (\circ)');
ylabel('intensity');
xlim([-2*pw 2*pw]);
title(sprintf('%d\\circ black, 1:5:90\\circ acceptance angles', pw/3));

%% subplot 2: distance between maxima
axes(h{2}); hold on;
plot(praas, int_distmaxs, 'k.');
xlabel('spatial resolution (hw of acc. angle) (\circ)');
ylabel('distance between white peaks (\circ)');
grid on; grid minor;

%% subplot 3: relative amplitude decay
axes(h{3}); hold on;
plot(praas, int_relamps, 'k');
xlabel('opt. resolution (hw of acc. angle) (\circ)');
ylabel('relative contrast (%)');
grid on; grid minor;

%% subplot 4: differential of amplitude decay
axes(h{4}); hold on;
plot(praas(1:end-1), diff(int_relamps)./diff(praas), 'k');
xlabel('opt. resolution (hw of acc. angle) (\circ)');
ylabel('change in modulation (%)');
grid on; grid minor;
% set(gca, 'xdir', 'reverse')

%% figure 2: example filtering steps
figure(2); clf;
subplot(2, 3, 1);
minlvl = min(pattern_y);
maxlvl = max(pattern_y);

reppedimage = repmat(pattern_y, [length(pattern_y) 1]);
scaledimage = 255*(reppedimage - minlvl)/(maxlvl-minlvl);
image(uint8(scaledimage));
axis image off;
title('unfiltered');

for hw = 10:10:50
    subplot(2, 3, hw/10+1);
    reppedimage = repmat(int_filt(:, praas==hw), [1 size(int_filt, 1)]);
    scaledimage = 255*(reppedimage - minlvl)/(maxlvl-minlvl);
    I = uint8(scaledimage)';

    image(I);
    axis image off;
    title(sprintf('%d\\circ opt. resolution', hw));
end

colormap(gray(256))

end % main

%% subfunctions
function [filtpattern, amp, distmins, distmaxs, modulation] = stimulus_analyse_one(pattern_x, pattern_y, gkern, roi)
% STIMULUS_ANALYSE_ONE filters a single pattern with a given optical resolution
% [filtpattern, amp, distmins, distmaxs, modulation] = stimulus_analyse_one(pattern_x, pattern_y, gkern, roi)
% 
% Inputs:
%   pattern_x   - x (angular) values describing the pattern; in degrees
%   pattern_y   - y (intensity) values describing the pattern; in any units
%   gkern       - Gaussian filter kernel
%   roi         - region of interest: a vector of indices into pattern_x/pattern_y. Only values in the region of interest will be considered for the calculation of 
%
% Outputs:
%   filtpattern - the full filtered pattern (same size and units as pattern_y)
%   amp         - amplitude of the pattern (max - min, only considering values inside the ROI, units as pattern_y)
%   distmins    - distance between the minima (in degrees; assuming a symmetric input pattern)
%   distmaxs    - distance between the maxima (in degrees; assuming a symmetric input pattern)
%   modulation  - modulation of the pattern as Michelson contrast: mod = (max - min) / (max + min)

    if nargin < 4, roi = 1:length(pattern_y); end % If no ROI has been provided, use the full pattern
    
    % Filter the pattern
    temp                = conv([pattern_y pattern_y pattern_y], gkern, 'same');     % replicate the pattern 3 times to avoid encroaching black borders
    filtpattern         = temp(length(pattern_y)+1:2*length(pattern_y));
    
    % Extract region of interest
    filtpattern_roi     = filtpattern(roi);                 
    pattern_x_roi       = pattern_x(roi);

    % Calculate further stats
    [filtmax, filtmaxi] = max(filtpattern_roi);
    [filtmin, filtmini] = min(filtpattern_roi);

    amp                 = filtmax - filtmin;
    modulation          = (filtmax - filtmin) / (filtmax + filtmin);
    distmins            = 2 * abs(pattern_x_roi(filtmini(1)));
    distmaxs            = 2 * abs(pattern_x_roi(filtmaxi(1)));
end

function h = sub_formatA4()
% SUB_FORMATA4 creates a figure with 6 equal-sized panels
% h = sub_formatA4()
%
% Outputs:
%   h      - cell array of axes handles, where h{i} if the handle to the i'th panel

    % Constants
    rimo = .1;     % outer border width, in pixels
    rimi = .1;    % width of border between panels, in pixels

    % Set up figure
    fh = figure(1); 
    clf;
    orient portrait;
    set(fh, 'position', [1 53 638 902], 'PaperType', 'A4', 'color', 'w', 'PaperPositionMode', 'auto');

    % Create panels
    pw   = (1-rimo-rimo-rimi)/2;
    ph   = (1-rimo-rimo-2*rimi)/3;
    h{1} = subplot('position', [rimo         rimo+2*ph+2*rimi pw ph]); hold on;
    h{2} = subplot('position', [rimo+pw+rimi rimo+2*ph+2*rimi pw ph]); hold on;
    h{3} = subplot('position', [rimo         rimo+1*ph+1*rimi pw ph]); hold on;
    h{4} = subplot('position', [rimo+pw+rimi rimo+1*ph+1*rimi pw ph]); hold on;
    h{5} = subplot('position', [rimo         rimo+0*ph+0*rimi pw ph]); hold on;
    h{6} = subplot('position', [rimo+pw+rimi rimo+0*ph+0*rimi pw ph]); hold on;
end