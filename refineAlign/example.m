opts = struct();
opts.cAligner = '/path/to/cAligner';
opts.EmSART = 'EmSART';
opts.EmSARTRefine = 'EmSARTRefine';
opts.STA = 'SubTomogramAverageMPI';
opts.STA_dir = 'path/to/SubTomogramAverageMPI';
opts.remote = true;
opts.host = 'localhost';


opts.projDir = '/project/directory/';
opts.iters = 10;
opts.motl = artia.em.read('/path/to/start/motl.em');
opts.projFile = '/path/to/tiltstack.st';
opts.tomoNr = 1;
opts.markerFile = '/path/to/markerfile.em';
opts.wedge = artia.em.read('/path/to/wedge.em');
opts.mask = artia.em.read('/path/to/mask.em');
opts.reference = artia.em.read('/path/to/reference.em'); 
opts.maskCC = artia.em.read('/path/to/maskCC.em');
opts.reconDim = [928 928 200]; %typical 1k size
opts.imDim = [3708 3838]; % 4k image stack
opts.volumeShifts = [0 0 0];
opts.groupSize = 20;
opts.maxShift = 30;
opts.speedUpDist = 100;
opts.voxelSize = 4;
opts.boxSize = 64;


iterative_alignment(opts)
