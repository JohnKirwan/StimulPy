function stimulus_main
% STIMULUS_MAIN is the main function for creating statistics to compare potential stimuli for phototaxis experiments.
%   1. Enter parameters in the following code block
%   2. Run the function and see resulting pdfs for a comparison of stimuli
%
% Written by Jochen Smolka and John D. Kirwan, Lund University, 2017
%
% This file is published as supplementary to the following article:
%   John D. Kirwan, Michael J. Bok, Jochen Smolka, James J. Foster, José Carlos Hernández, Dan-Eric Nilsson (2017) 
%   The sea urchin Diadema africanum uses low resolution vision to find shelter and deter enemies. Journal of Experimental Biology.
close all;

%% parameters
acc         = 0.1;                % accuracy of spatial sampling, in degrees; use 1 or 2 for testing, 0.1 for final output
acc2        = 0.1;                % accuracy of resolution sampling, in degrees; use 1 or 2 for testing, 0.1 for final output
wls         = 1:acc2:180;       % range of pattern wavelengths to test in plots 5&6
patterns    =   {'square'};        % add several patterns here to calculate pdfs for each one
plotit      = false;
% Pattern type options:
%                'bar'          % a single black bar of width wl/2 on a white background
%                'singlegauss'  % single dark gaussian (with half-width wl/2) on white background
%                'square'       % square wave pattern: black stripe of width wl/2 flanked by two white stripes of half the width
%                '2square'      % square wavelet pattern (Haar wavelet): black stripe and white stripe, each of width wl/2, on grey background
%                'dog'          % difference of Gaussians, central Gaussian with half-width wl/2
%                'dog2'         % difference of Gaussians with different half-width ratio
%                'log'          % Laplacian of Gaussian, v1
%                'mylog'        % Laplacian of Gaussian, v2
%                'cos'          % piece-wise sine with wavelenth wl (as used for Onychophorans); ratio ignored, always 0.5/1/0.5
%                'sin'          % continous sine-wave pattern with wavelength wl
%                'multiple'     % 7 progressively smaller cosines
ratio       = 1;                % amplitude ratio of inner and outer stripes 
                                % for DoG, ratio+1 is the ratio between sigmas
%% constants
va          = -180:acc:180;     % visual angles over which the pattern is observed
examwidth   = 10;               % example pattern width to be plotted in plots 1-4
acc_type    = 'g';              % shape of acceptance angle function ('g' = gaussian, 's' = square/tophat, 'm' = mixture of the two)

%% main loop
% for each pattern type, create two figures, one with diagnostic plots and one with an example pattern filtered at different spatial resolutions
cores = nan(length(patterns), 4); % pre-allocate cut-off resolution; for each pattern width, 
                                  % this will list the animal's spatial resolution, at which pattern contrast drops below 20%/10%/5%

for j = 1:length(patterns)
    close all;
    pattern = patterns{j};
    fprintf('Analysing pattern %s.\n\tPattern width: ', pattern);
    
    % calculate one example, and plot it in detail in plots 1-4
    int = stimulus_create(pattern, examwidth, va, ratio) + 1;       % create the stimulus (int contains the intensity for each visual angle in va)
    ah  = stimulus_analyse(va, int, 3*examwidth, acc_type, plotit);   % calculates diagnostic plots for example pattern, and returns the axes handles to subplots for later use
    
    % calculate cut-off resolutions for different pattern widths
    h = waitbar(0, 'Calculating cut-off resolutions');
    for i = 1:length(wls)
        wl               = wls(i);
        int              = stimulus_create(pattern, wl, va, ratio) + 1;
        physwidth{j}(i)  = 2 * va(find(va>0 & int>=1, 1, 'first')-1); % Distance between the pattern's zero-crossings
        cores(i, 1:4)    = stimulus_analyse(va, int, 3*wl, acc_type);
        waitbar(i/length(wls), h);
    end
    close(h);
    fprintf('done.\n');
    
    % now plot information over all pattern widths in plots 5&6
    
    plot(ah{5}, wls, cores(:, 1), 'g', wls, cores(:, 2), 'k', wls, cores(:, 3), 'b', wls, cores(:, 4), 'r');
    xlabel(ah{5}, 'pattern width (\circ)');
    ylabel(ah{5}, 'spatial res. for varying contr. sens. (\circ)');
    grid(ah{5}, 'on'); 
    grid(ah{5}, 'minor');
    legend(ah{5}, '20%', '10%', '5%', '1%', 'location', 'SE');
    
    x{j} = cores(:, 2);
    y{j} = cores(:, 3)-cores(:, 1);
    plot(ah{6}, x{j}, y{j}, 'k');
    xlabel(ah{6}, 'spatial res. (\circ)');
    ylabel(ah{6}, 'uncertainty (\circ)');
    grid(ah{6}, 'on'); grid(ah{6}, 'minor');   
    save(sprintf('results_%s_%d_%s_%.1f_%.1f.mat', pattern, ratio, acc_type, acc, acc2));
    
    pdfsave(1, sprintf('results_%s_%d_%s_%.1f_%.1f.pdf', pattern, ratio, acc_type, acc, acc2));
    pdfsave(2, sprintf('results_%s_%d_%s_%.1f_%.1f_example.pdf', pattern, ratio, acc_type, acc, acc2));
end

%% plot a comparison between all tested pattern (an overlay of all subplots #6)
figure(3); clf; hold on;
for j = 1:length(x)
    plot(physwidth{j}, x{j});
end
xlabel('pattern width (\circ)');
ylabel('spatial res. (\circ)');
grid on; grid minor;
legend(patterns);
pdfsave(3, sprintf('patterncomp_%d_%s_width.pdf', ratio, acc_type));

figure(4); clf; hold on;
for j = 1:length(x)
    plot(x{j}, y{j});
end
xlabel('spatial res. (\circ)');
ylabel('uncertainty (\circ)');
grid on; grid minor;
legend(patterns);
pdfsave(4, sprintf('patterncomp_%d_%s.pdf', ratio, acc_type));

end % main

%% sub functions
function pdfsave(fignum, filename)
    % Saves figure to pdf
    fh = figure(fignum);

    set(fh, 'Units', 'centimeters');
    pos = get(fh, 'Position');
    set(fh, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)]);

    set(fh, 'renderer', 'painters');
    print(fh, '-dpdf', '-r1200', filename);
    disp(['Saved figure ' num2str(fignum) ' to ' filename]);
end
