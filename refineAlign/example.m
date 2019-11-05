opts = struct();

% General options
opts.projDir = '/project/directory/'; % Any non-existing directory
opts.iters = 10;

% Executable locations and remote setup
opts.cAligner = '/path/to/cAligner';
opts.EmSART = 'EmSART';
opts.EmSARTRefine = 'EmSARTRefine';
opts.STA = 'SubTomogramAverageMPI';
opts.STA_dir = 'path/to/SubTomogramAverageMPI';
opts.remote = true;
opts.host = 'localhost';

% Reconstruction parameters
opts.projFile = '/path/to/tiltstack.st';
opts.markerFile = '/path/to/markerfile.em';
opts.reconDim = [928 928 350]; %typical 1k size
opts.imDim = [927 959]; % 1k image stack
opts.volumeShifts = [0 0 -100];
opts.maAmount = 1;
opts.maAngle = 0;
opts.voxelSize = 1;

% Averaging parameters
opts.tomoNr = 78;
opts.boxSize = 64;
opts.motl = artia.em.read('/path/to/start/motl.em');
opts.wedge = artia.em.read('/path/to/wedge.em');
opts.mask = artia.em.read('/path/to/mask.em');
opts.reference = artia.em.read('/path/to/reference.em'); 
opts.maskCC = artia.em.read('/path/to/maskCC.em');
opts.angIter = 3;
opts.angIncr = 5;
opts.phiAngIter = 3;
opts.phiAngIncr = 5;
opts.avgLowPass = 12;
opts.avgHighPass = 0;
opts.avgSigma = 3;

% Refinement (projection matching) parameters
opts.groupMode = 'MaxDistance';
opts.maxDistance = 150;
opts.groupSize = 20;
opts.maxShift = 15;
opts.speedUpDist = 60;

% Refinement band pass filter
opts.lowPass = 200;
opts.lowPassSigma = 50;
opts.highPass = 20;
opts.highPassSigma = 10;

% Volume size computation
opts.borderSize = 5*opts.boxSize;

% RUN!
iterative_alignment(opts)