def stimulus_fft(type, wl, ratio):
    # STIMULUS_FFT uses fast fourier transform to plot the amplitude spectrum of a potential stimulus pattern
    # Inputs:
    #   type   - type of stimulus, e.g. 'bar'/'dog'/'square'/'log' (see below for full list)
    #   wl     - wavelength of stimulus in degrees (definition depends on stimulus stype)
    #   ratio  - ratio between black and white areas; usually ratio of amplitudes, but definition depends on stimulus type
    # example usage: stimulus_fft('dog', 10, 1)
    # Written by Jochen Smolka and John D. Kirwan, Lund University, 2017
    # This file is published as supplementary to the following article:
    #   John D. Kirwan, Michael J. Bok, Jochen Smolka, James J. Foster, José Carlos Hernández, Dan-Eric Nilsson (2017)
    #   The sea urchin Diadema africanum uses low resolution vision to find shelter and deter enemies. Journal of Experimental Biology.
    import numpy as np
    import matplotlib.pyplot as plt
     ## set defaults
    if ratio == None:
        ratio=1;
    if wl == None:
        wl = 10;
    if type == None:
        type = 'bar';

    va  = np.arange(-180,179.99,0.1).tolist();
    Fs  = len(va);         # Simulation resolution (points per 360 degrees)
    azi = np.linspace(0,Fs-1)/Fs * 360 - 180;      # azimuth vector
    # check indexing 0 or 1 above
    S   = stimulus_create(type, wl, va, ratio); # *so far, empty array*
    
    plt.figure(1);
    plt.clf();
    plt.plot(azi,S)
    plt.title( print('%s, wl=%d', type, wl))
    plt.xlabel('azimuth (\circ)')
    plt.ylabel('intensity')
    plt.axis([-180, 180, -1.1, 1.1])
    
    Y           = np.fft(S);
    P2          = abs(Y/Fs);
    P1          = P2(np.linspace(1,Fs/2+1,1)); #P2(1:Fs/2+1);   
    P1[1:-1] = 2*P1(np.linspace(2,-1,1));  #Pythonic indexing
    f           = Fs*(np.linspace(0,Fs/2,1))/Fs/360; 
    
    plt.figure(2);
    plt.clf();
    plt.plot(f, P1)
    plt.title( print(wl=%d', type, wl, format(5))) 
    plt.xlabel('f (cycles / deg)')
    plt.ylabel('|P1(f)|')
    plt.xlim([0, 1])

## Matlab is 1 indexed, while Python is zero indexed! Need to change all