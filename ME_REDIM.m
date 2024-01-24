function ME_REDIM(fdwi,fbvec,fbval,fmask,TE,outpath,options)
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Joint cumulant moments of relaxation and diffusion coefficients with the REDIM approach. 
    (10.1109/TMI.2019.2933982)     
    %}
    
    arguments
        fdwi                    (1,:)    {mustBeNonzeroLengthText}
        fbvec                   (1,:)    {mustBeNonzeroLengthText}
        fbval                   (1,:)    {mustBeNonzeroLengthText}
        fmask                   string   {mustBeFile}
        TE                      (1,:)    {mustBeNumeric}
        outpath                 string   {mustBeNonzeroLengthText}

        options.approach        (1,1)    {mustBeNonzeroLengthText} = 'wl'
    end

    addpath('third/redim');
    addpath('third/matlab');

    cellfun(@(x)assert(exist(x,'file'),'Input DWI %s does not exist', x),fdwi,'UniformOutput',false);
    cellfun(@(x)assert(exist(x,'file'),'Input Bvec %s does not exist', x),fbvec,'UniformOutput',false);
    cellfun(@(x)assert(exist(x,'file'),'Input Bval %s does not exist', x),fbval,'UniformOutput',false);
    assert(exist(fmask,'file'),'Input mask %s does not exist', fmask);
    assert(exist(options.spectrum,'file'),'Input spectrum %s does not exist', options.spectrum);

    %% load multi-echo dMRI dataset
    redim_pipe(fdwi, fullfile(outpath,options.approach), fmask, fbvec{1}, fbval{1}, TE,options.approach);
