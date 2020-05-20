def stimulus_main():
    # STIMULUS_MAIN is the main function for creating statistics to compare potential stimuli for phototaxis experiments.
    # 1. Enter parameters in the following code block
    # 2. Run the function and see resulting pdfs for a comparison of stimuli
    # Written by Jochen Smolka and John D. Kirwan, Lund University, 2017
    # This file is published as supplementary to the following article:
    import matplotlib.pyplot as plt
    import numpy as np
    from progress.bar import Bar
    plt.close("all")
    ## parameters
    acc      = 1; # accuracy of spatial sampling, in degrees; use 1 or 2 for testing, 0.1 for final output
    acc2     = 1; # accuracy of resolution sampling, in degrees; use 1 or 2 for testing, 0.1 for final output
    wls      = np.linspace(1,180,180/acc2+1);# range of pattern wavelengths to test in plots 5&6
    patterns =   {'dog'}; # add several patterns here to calculate pdfs for each one
    # Pattern type options:
    # 'bar'   # a single black bar of width wl/2 on a white background
    # 'singlegauss' # single dark gaussian (with half-width wl/2) on white background
    # 'square' # square wave pattern: black stripe of width wl/2 flanked by two white stripes of half the width
    # '2square'# square wavelet pattern (Haar wavelet): black stripe and white stripe, each of width wl/2, on grey background
    # 'dog'    # difference of Gaussians, central Gaussian with half-width wl/2
    # 'dog2'   # difference of Gaussians with different half-width ratio
    # 'log'    # Laplacian of Gaussian, v1
    # 'mylog'  # Laplacian of Gaussian, v2
    # 'cos'    # piece-wise sine with wavelenth wl (as used for Onychophorans); ratio ignored, always 0.5/1/0.5
    # 'sin'    # continous sine-wave pattern with wavelength wl
    # 'multiple' # 7 progressively smaller cosines
    ratio = 1; # amplitude ratio of inner and outer stripes
               # for DoG, ratio+1 is the ratio between sigmas
    #### constants
    va         = va = np.linspace(-180,180,360/acc+1); # visual angles over which the pattern is observed
    examwidth  = 10;           # example pattern width to be plotted in plots 1-4
    acc_type   = 'g';          # shape of acceptance angle function ('g' = gaussian, 's' = square/tophat, 'm' = mixture of the two)
    #### main loop
    # for each pattern type, create two figures, one with diagnostic plots and one with an example pattern filtered at different spatial resolutions
    cores    = np.empty((np.size(patterns), 4,)); # pre-allocate cut-off resolution; for each pattern width,
    cores[:] = np.nan;
    # this will list the animal's spatial resolution, at which pattern contrast drops below a specified contrast
    for j in patterns: #(np.arange(np.size(patterns))):
        plt.close("all")
        pattern = j;
        #fprintf('Analysing pattern %s.\n\tPattern width: ', pattern);
        print('Analysing pattern \n\tPattern width: ', j);
    
        # calculate one example, and plot it in detail in plots 1-4
        int = stimulus_create(pattern, examwidth, va, ratio) + 1;       # create the stimulus (int contains the intensity for each visual angle in va)
        ah  = stimulus_analyse(va, int, 3*examwidth, acc_type, true);   # calculates diagnostic plots for example pattern, and returns the axes handles to subplots for later use
        # calculate cut-off resolutions for different pattern widths
        h = waitbar(0, 'Calculating cut-off resolutions');
 #      h = Bar()  #python progess bar
        
        for i in (np.arange(np.size(wls))):
            wl   = wls(i);
            int  = stimulus_create(pattern, wl, va, ratio) + 1;
            physwidth[j][i]  = 2 * va(find(va>0 & int>=1, 1, 'first')-1); # Distance between the pattern's zero-crossings
            cores[i, 1:4]    = stimulus_analyse(va, int, 3*wl, acc_type);
            waitbar(i/len(wls), h);
        close(h);
        print('done.\n');
    
        # now plot information over all pattern widths in plots 5&6
        plt.plot(ah{5}, wls, cores(:, 1), 'g', wls, cores(:, 2), 'k', wls, cores(:, 3), 'b', wls, cores(:, 4), 'r');
        plt.xlabel(ah{5}, 'pattern width (\circ)');
        plt.ylabel(ah{5}, 'spatial res. for varying contr. sens. (\circ)');
        plt.grid(ah{5}, 'on');
        plt.grid(ah{5}, 'minor');
        plt.legend(ah{5}, '20#', '10#', '5#', '1#', 'location', 'SE');
    
        x{j} = cores(:, 2);
        y{j} = cores(:, 3)-cores(:, 1);
        plt.plot(ah{6}, x{j}, y{j}, 'k');
        plt.xlabel(ah{6}, 'spatial res. (\circ)');
        plt.ylabel(ah{6}, 'uncertainty (\circ)');
        plt.grid(ah{6}, 'on'); grid(ah{6}, 'minor');
        save(sprintf('results_#s_#d_#s_#.1f_#.1f.mat', pattern, ratio, acc_type, acc, acc2));
    
        pdfsave(1, sprintf('results_#s_#d_#s_#.1f_#.1f.pdf', pattern, ratio, acc_type, acc, acc2));
        pdfsave(2, sprintf('results_#s_#d_#s_#.1f_#.1f_example.pdf', pattern, ratio, acc_type, acc, acc2));
    
    ## plot a comparison between all tested pattern (an overlay of all subplots #6)
    plt.figure(3); plt.clf; hold on;
    for j in (np.arange(np.size(x))):
        plot(physwidth{j}, x{j});
    plt.xlabel('pattern width (\circ)');
    plt.ylabel('spatial res. (\circ)');
    plt.grid on; grid minor;
    plt.legend(patterns);
    pdfsave(3, sprintf('patterncomp_#d_#s_width.pdf', ratio, acc_type));
    f.plt.savefig("foo.pdf", bbox_inches='tight') # save 1 page PDF
                       
    plt.figure(4); pyplot.clf(); hold on;
    for j = np.linspace(1,len(x)):
        plt.plot(x{j}, y{j});
    plt.xlabel('spatial res. (\circ)');
    plt.ylabel('uncertainty (\circ)');
    plt.grid on; plt.grid minor;
    plt.legend(patterns);
    pdfsave(4, print('patterncomp_#d_#s.pdf', ratio, acc_type));
    
    ## sub functions
    function pdfsave(fignum, filename)
        # Saves figure to pdf
        fh = plt.figure(fignum);
    
        plt.set(fh, 'Units', 'centimeters');
        pos = plt.get(fh, 'Position');
        plt.set(fh, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)]);
    
        plt.set(fh, 'renderer', 'painters');
        print(fh, '-dpdf', '-r1200', filename);
        disp(['Saved figure ' num2str(fignum) ' to ' filename]);
        
        plt.show()
