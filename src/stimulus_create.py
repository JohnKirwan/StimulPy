def stimulus_create(type, wl, va, ratio):
    # STIMULUS_CREATE generates a 2D stimulus vector for phototaxis experiments of desired type and resolution
    # Inputs:
    #   type   - type of stimulus, e.g. 'bar'/'dog'/'square'/'log' (see below for full list)
    #   wl     - "width" or "wavelength" of pattern in degrees (definition varies with stimulus type)
    #   va     - visual angles over which to create the stimulus (in degrees)
    #   ratio  - ratio between black and white areas; usually ratio of amplitudes, but definition depends on stimulus type
    # Outputs:
    #   outint - stimulus intensity over all angles in va (arbitraty units, normalize later)
    # Pattern type options:
    # 'bar'         # a single black bar of width wl/2 on a white background
    # 'singlegauss' # single dark gaussian (with half-width wl/2) on white background
    # 'square'    # square wave pattern: black stripe of width wl/2 flanked by two white stripes of half the width
    # '2square'   # square wavelet pattern (Haar wavelet): black stripe and white stripe, each of width wl/2, on grey background
    # 'dog'       # difference of Gaussians, central Gaussian with half-width wl/2
    # 'dog2'      # difference of Gaussians with different half-width ratio
    # 'log'       # Laplacian of Gaussian, v1
    # 'mylog'     # Laplacian of Gaussian, v2
    # 'cos'       # piece-wise sine with wavelenth wl (as used for Onychophorans); ratio ignored, always 0.5/1/0.5
    # 'sin'       # continous sine-wave pattern with wavelenth wl
    # 'multiple'  # 7 progressively smaller cosines
    # example usage: outint = stimulus_create('dog', 10, -180:0.1:180, 1)
    # Written by Jochen Smolka and John D. Kirwan, Lund University, 2017

#   import matplotlib.pyplot as plt #for sanity checking
    import numpy as np
    ## set defaults
    if ratio == None:
        ratio=1;
    if va == None:
        va=np.linspace(-180,180,3601);
    if wl == None:
        wl = 10;
    if type == None:
        type = 'bar';
    
    outint = np.zeros(np.size(va)); # preallocate vector of pattern intensity
    # this vector is populated with 
    ## for all piece-wise defined patterns (e.g. square, triangle, sin), define the regions of black and white "stripes"
    sel_whiteleft  = np.logical_and(np.greater_equal(va,-(ratio+1)*wl/2),
                                    np.greater(va,-wl/2) );
    sel_black      = np.logical_and(np.greater_equal(va,-wl/2),
                                    np.greater(va, wl/2));
    sel_whiteright = np.logical_and(np.greater_equal(va,wl/2),
                                    np.greater(va,(ratio+1)*wl/2));
    amps           = [1, 1/ratio]; # relative amplitudes of black and white areas
    
    ## main switch to create patterns
    if type.lower() == 'bar': # single black bar (of width wl/2) on white background
    #switch lower(type)
        ## single bar patterns
        #case 'bar'          
        sel_black           = np.logical_and(np.greater_equal(va,-wl/4,np.less(va, wl/4);
        outint[sel_black]   = -1;
    elif type.lower() == 'singlegauss':  # single dark gaussian (with half-width wl/2) on white background
        #case 'singlegauss'
        sigma  = wl/2 / (2*np.sqrt(2*np.log(2)));   # sigma of Gaussian
        outint = - np.exp(-va**2/(2*sigma**2));   # was -.va^2
    
        ## balanced patterns (integral = 0)
    elif type.lower() == 'square': # square wave pattern: black stripe of width wl/2 flanked by two white stripes of half the width
    #case 'square'                    
        sel_whiteleft       = np.logical_and(np.greater_equal(va,-(ratio+1)*wl/4),np.less(va, -wl/4));
        sel_black           = np.logical_and(np.greater_equal(va,-wl/4),np.less(va, wl/4));
        sel_whiteright      = np.logical_and(np.greater_equal(va, wl/4),np.less(va, (ratio+1)*wl/4));
        outint[sel_black]      = -amps[1];
        outint[sel_whiteleft]  =  amps[2];
        outint[sel_whiteright] =  amps[2];
    elif type.lower() == '2square': # square wavelet pattern (Haar wavelet): black stripe and white stripe, each of width wl/2, on grey background
    #case '2square'                  
        sel_white           = np.logical_and(np.greater_equal(va,-wl/2),np.less(va, 0));
        sel_black           = np.logical_and(np.greater_equal(va,0),np.less(va, wl/2)); # logical of indices of va between 0 and wl/2
        outint[sel_black]   = -amps[1]; # apply this negative ampltiude to the true values of sel_black
        outint[sel_white]   =  amps[2];
    elif type.lower() == 'dog': # difference of Gaussians pattern
        fwhm1   = wl/2;                           # half-width of primary (black) Gaussian; also equals half the distance between white peaks
        fwhm2   = wl/2 * (ratio+1);               # half-width of secondary (white) Gaussian
        sigma1  = fwhm1 / (2*np.sqrt(2*np.log(2)));     # sigma of primary Gaussian
        sigma2  = fwhm2 / (2*np.sqrt(2*np.log(2)));     # sigma of secondary Gaussian
        amp1    = -1;
        amp2    = -amp1 / (fwhm2/fwhm1);
        outint  = amp1 * np.exp(-va**2/(2*sigma1**2)) + amp2 * np.exp(-va**2/(2*sigma2**2));
    elif type.lower() == 'dog2': # difference of Gaussians pattern with different half-width ratio
        fwhm1   = wl;                           # half-width of primary (black) Gaussian; also equals half the distance between white peaks
        fwhm2   = wl * 1.2;                     # half-width of secondary (white) Gaussian
        sigma1  = fwhm1 / (2*np.sqrt(2*np.log(2)));   # sigma of primary Gaussian
        sigma2  = fwhm2 / (2*np.sqrt(2*np.log(2)));   # sigma of secondary Gaussian
        amp1    = -1*6;
        amp2    = -amp1 / (fwhm2/fwhm1);
        outint  = amp1 * np.exp(-va**2/(2*sigma1**2)) + amp2 * np.exp(-va**2/(2*sigma2**2));
    elif type.lower() == 'log': # Laplacian of Gaussian (not working properly?)
        from scipy import ndimage
        sigma   = wl/2 / (2*np.sqrt(2*np.log(2)));      # sigma
        acc     = np.median(np.diff(va));             # degrees/pixel
        #outint  = fspecial('log', np.size(va), sigma/acc);
        outint  = ndimage.gaussian_laplace(outint, sigma/acc) #array and sigma of filter 
    elif type.lower() == 'mylog':  # Laplacian of Gaussian
        sigma   = wl/2 / (2*np.sqrt(2*np.log(2)));
        outint  = -(1-va**2/sigma**2)*np.exp(-va**2/(2*sigma**2)); #was .*exp()
    elif type.lower() == 'cos':  # piece-wise sine with wavelenth wl (as used for Onychophorans)
        sel_whiteleft       = np.logical_and(np.greater_equal(va,-3/4*wl),np.less(va, -1/4*wl);
        sel_black           = np.logical_and(np.greater_equal(va,-1/4*wl),np.less(va,  1/4*wl);
        sel_whiteright      = np.logical_and(np.greater_equal(va, 1/4*wl),np.less(va,  3/4*wl);
        amps                = [1, 1/2]; # relative amplitudes of black and white areas
        outint[sel_black]   = - amps[1] * np.cos(np.deg2rad(np.linspace(-90,90,np.count_nonzero(sel_black))));
        ### CAN'T ASSIGN TO F(X) CALL
        # was linspace and nnz() in matlab
        outint[sel_whiteleft]    =   amps[2] * np.cos(np.deg2rad(np.linspace(-90,90,np.count_nonzero(sel_whiteleft))));
        outint[sel_whiteright]   =   amps[2] * np.cos(np.deg2rad(np.linspace(-90,90,np.count_nonzero(sel_whiteright))));
    elif type.lower() == 'sin': # continous sine-wave pattern with wavelenth wl
        outint =  np.cos(np.deg2rad(va*360/wl));
    elif type.lower() == 'multiple': # piece-wise cosine with 7 pieces (was intended to more closely simulate a continuous sine-wave)
        sel_white2left          = np.logical_and(np.greater_equal(va,-(3*ratio+1)*wl/2,np.less(va,= -(2*ratio+1)*wl/2);
        sel_black2left          = np.logical_and(np.greater_equal(va,-(2*ratio+1)*wl/2,np.less(va,= -(ratio+1)*wl/2);
        sel_black2right         = np.logical_and(np.greater_equal(va,(1*ratio+1)*wl/2,np.less(va,= (2*ratio+1)*wl/2);
        sel_white2right         = np.logical_and(np.greater_equal(va,(2*ratio+1)*wl/2,np.less(va,= (3*ratio+1)*wl/2);
        amps                    = [1, 4/3/ratio, 2/3/ratio, 2/6/ratio];
        outint[sel_black]       = - amps[1] * np.cos(np.deg2rad(np.linspace(-90,90,np.count_nonzero(sel_black))));
        outint[sel_whiteleft]   =   amps[2] * np.cos(np.deg2rad(np.linspace(-90,90,np.count_nonzero(sel_whiteleft))));
        outint[sel_whiteright]  =   amps[2] * np.cos(np.deg2rad(np.linspace(-90,90,np.count_nonzero(sel_whiteright))));
        outint[sel_black2left]  = - amps[3] * np.cos(np.deg2rad(np.linspace(-90,90,np.count_nonzero(sel_black2left))));
        outint[sel_black2right] = - amps[3] * np.cos(np.deg2rad(np.linspace(-90,90,np.count_nonzero(sel_black2right))));
        outint[sel_white2left]  =   amps[4] * np.cos(np.deg2rad(np.linspace(-90,90,np.count_nonzero(sel_white2left))));
        outint[sel_white2right] =   amps[4] * np.cos(np.deg2rad(np.linspace(-90,90,np.count_nonzero(sel_white2right))));
    
    return(outint)
