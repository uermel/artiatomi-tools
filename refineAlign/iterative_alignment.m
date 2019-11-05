function iterative_alignment(opts)
    
    % Make project dir
    opts.projDir = saneDir(opts.projDir);
    mkdir(opts.projDir);
    
    for i = 1:opts.iters
        % Make directories and names and save files %%%%%%%%%%%%%%%%%%%%%%%
        disp('------------------------------------------')
        disp(['Iter ' num2str(i)])
        disp('Preparing ...')
        iterDir = sprintf([opts.projDir 'iter%d/'], i); 
        motlDir = [iterDir 'motl/'];
        refDir = [iterDir 'ref/'];
        partDir = [iterDir 'parts/'];
        ccmapDir = [iterDir 'ccmap/'];
        shiftDir = [iterDir 'shift/'];
        mkdir(motlDir);
        mkdir(refDir);
        mkdir(partDir);
        mkdir(ccmapDir);
        mkdir(shiftDir);
        
        motlFilePre = [motlDir 'motl_'];
        extractMotlFile = [motlFilePre '1.em'];
        artia.em.write(opts.motl, extractMotlFile);
        
        refFilePre = [refDir 'ref_'];
        artia.em.write(opts.reference, [refFilePre '1.em']);
        
        partFilePre = [partDir 'part_'];
        
        maskFile = [iterDir 'mask.em'];
        maskCCFile = [iterDir 'maskCC.em'];
        wedgeFile = [iterDir 'wedge.em'];
        refFile = [iterDir 'reference.em'];
        currMarkFile = [iterDir 'current_marker.em'];
        artia.em.write(opts.mask, maskFile);
        artia.em.write(opts.maskCC, maskCCFile);
        artia.em.write(opts.wedge, wedgeFile);
        artia.em.write(opts.reference, refFile);
        copyfile(opts.markerFile, currMarkFile);
        mark = artia.marker.read(currMarkFile);
        phi = mark.beamDeclination;
        
        oss10ReconFile = sprintf([iterDir 'tomo_%d_oss10_iter_%d.em'], opts.tomoNr, i);
        oss1ReconFile = sprintf([iterDir 'tomo_%d_oss1_iter_%d.em'], opts.tomoNr, i);
        oss1FilteredReconFile = sprintf([iterDir 'tomo_%d_oss1_filtered_iter_%d.em'], opts.tomoNr, i);
        
        shiftFile = sprintf([shiftDir 'tomo_%d_shifts_iter_%d.em'], opts.tomoNr, i);
        ccmapPre = sprintf([ccmapDir 'tomo_%d_iter_%d_ccmap_'], opts.tomoNr);
        refMotlFile = [motlFilePre '2.em'];
        
        oss10CfgFile = sprintf([iterDir 'tomo_%d_oss10_iter_%d.cfg'], opts.tomoNr, i);
        oss1CfgFile = sprintf([iterDir 'tomo_%d_oss1_iter_%d.cfg'], opts.tomoNr, i);
        oss1FilteredCfgFile = sprintf([iterDir 'tomo_%d_oss1_filtered_iter_%d.cfg'], opts.tomoNr, i);
        refCfgFile = sprintf([iterDir 'tomo_%d_refine_iter_%d.cfg'], opts.tomoNr, i);
        avgCfgFile = sprintf([iterDir 'avg_iter_%d.cfg'], i);
        
        % Make configs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Writing cfgs ...')
        oss10 = recon_blankConfig();
        oss1 = recon_blankConfig();
        oss1_filtered = recon_blankConfig();
        refine = refine_blankConfig();
        avg = avg_blankConfig();
        
        %oss10
        oss10.ProjectionFile = opts.projFile;
        oss10.OutVolumeFile = oss10ReconFile;
        oss10.MarkerFile = currMarkFile;
        oss10.RecDimesions = sprintf('%d %d %d', round(opts.reconDim(1)), round(opts.reconDim(2)), round(opts.reconDim(3)));
        oss10.VolumeShift = sprintf('%d %d %d', opts.volumeShifts(1), opts.volumeShifts(2), opts.volumeShifts(3));
        oss10.PhiAngle = num2str(phi);
        oss10.VoxelSize = num2str(opts.voxelSize);
        oss10.SIRTCount = '10';
        oss10.MagAnisotropy = [num2str(opts.maAmount) ' ' num2str(opts.maAngle)];
        artia.cfg.write(oss10, oss10CfgFile);
        
        %oss1
        oss1.ProjectionFile = opts.projFile;
        oss1.OutVolumeFile = oss1ReconFile;
        oss1.MarkerFile = currMarkFile;
        oss1.RecDimesions = sprintf('%d %d %d', round(opts.reconDim(1)), round(opts.reconDim(2)), round(opts.reconDim(3)));
        oss1.VolumeShift = sprintf('%d %d %d', opts.volumeShifts(1), opts.volumeShifts(2), opts.volumeShifts(3));
        oss1.PhiAngle = num2str(phi);
        oss1.VoxelSize = num2str(opts.voxelSize);
        oss1.SIRTCount = '1';
        oss1.MagAnisotropy = [num2str(opts.maAmount) ' ' num2str(opts.maAngle)];
        artia.cfg.write(oss1, oss1CfgFile);
        
        %oss1_filtered
        oss1_filtered.ProjectionFile = opts.projFile;
        oss1_filtered.OutVolumeFile = oss1FilteredReconFile;
        oss1_filtered.MarkerFile = currMarkFile;
        oss1_filtered.RecDimesions = sprintf('%d %d %d', round(opts.reconDim(1)), round(opts.reconDim(2)), round(opts.reconDim(3)));
        oss1_filtered.VolumeShift = sprintf('%d %d %d', opts.volumeShifts(1), opts.volumeShifts(2), opts.volumeShifts(3));
        oss1_filtered.PhiAngle = num2str(phi);
        oss1_filtered.VoxelSize = num2str(opts.voxelSize);
        oss1_filtered.SIRTCount = '1';
        oss1_filtered.SkipFilter = 'false';
        oss1_filtered.fourFilterLP = num2str(opts.lowPass);
        oss1_filtered.fourFilterLPS = num2str(opts.lowPassSigma);
        oss1_filtered.fourFilterHP = num2str(opts.highPass);
        oss1_filtered.fourFilterHPS = num2str(opts.highPassSigma);
        oss1_filtered.MagAnisotropy = [num2str(opts.maAmount) ' ' num2str(opts.maAngle)];
        artia.cfg.write(oss1_filtered, oss1FilteredCfgFile);
        
        %refine
        refine.CudaDeviceID = '0';
        refine.ProjectionFile = opts.projFile;
        refine.OutVolumeFile = oss1FilteredReconFile;
        refine.MarkerFile = currMarkFile;
        refine.RecDimesions = sprintf('%d %d %d', round(opts.reconDim(1)), round(opts.reconDim(2)), round(opts.reconDim(3)));
        refine.VolumeShift = sprintf('%d %d %d', opts.volumeShifts(1), opts.volumeShifts(2), opts.volumeShifts(3));
        refine.PhiAngle = num2str(phi);
        refine.VoxelSize = num2str(opts.voxelSize);
        refine.SIRTCount = '1';
        refine.SkipFilter = 'false';
        refine.fourFilterLP = num2str(opts.lowPass);
        refine.fourFilterLPS = num2str(opts.lowPassSigma);
        refine.fourFilterHP = num2str(opts.highPass);
        refine.fourFilterHPS = num2str(opts.highPassSigma);
        refine.MagAnisotropy = [num2str(opts.maAmount) ' ' num2str(opts.maAngle)];
        refine.GroupMode = opts.groupMode;
        refine.GroupSize = num2str(opts.groupSize);
        refine.MaxShift = num2str(opts.maxShift);
        refine.MaxDistance = num2str(opts.maxDistance);
        refine.Reference = refFile;
        refine.SizeSubVol = num2str(opts.boxSize);
        refine.SpeedUpDistance = num2str(opts.speedUpDist);
        refine.VoxelSizeSubVol = num2str(opts.voxelSize);
        refine.CCMapFileName = ccmapPre;
        refine.MotiveList = refMotlFile;
        refine.ShiftOutputFile = shiftFile;
        refine.MultiPeakDetection = 'false';
        artia.cfg.write(refine, refCfgFile);
        
        %avg
        avg.MotiveList = motlFilePre;
        avg.WedgeFile = wedgeFile;
        avg.Particles = partFilePre;
        avg.Reference = refFilePre;
        avg.Mask = maskFile;
        avg.MaskCC = maskCCFile;
        avg.NamingConvention = 'TomoParticle';
        avg.StartIteration = '1';
        avg.EndIteration = '2';
        avg.AngIter = num2str(opts.angIter);
        avg.AngIncr = num2str(opts.angIncr);
        avg.PhiAngIter = num2str(opts.phiAngIter);
        avg.PhiAngIncr = num2str(opts.phiAngIncr);
        avg.LowPass = num2str(opts.avgLowPass);
        avg.HighPass = num2str(opts.avgHighPass);
        avg.Sigma = num2str(opts.avgSigma);
        artia.cfg.write(avg, avgCfgFile);
        
        % Run the reconstructions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Running reconstructions ...')
        artia.mpi.run(opts.EmSART, 4, oss10CfgFile, 'runRemote', opts.remote, 'remoteHost', opts.host, 'suppressOutput', false)
        artia.mpi.run(opts.EmSART, 4, oss1CfgFile, 'runRemote', opts.remote, 'remoteHost', opts.host, 'suppressOutput', false)
        artia.mpi.run(opts.EmSART, 4, oss1FilteredCfgFile, 'runRemote', opts.remote, 'remoteHost', opts.host, 'suppressOutput', false)
        
        % Get Parts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tomolist = {};
        tomolist{opts.tomoNr} = oss1ReconFile;
        artia.particle.extract_write(extractMotlFile, 1, tomolist, opts.boxSize/2, 0, 0, 1, partFilePre);        
        
        % Run the averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Running averaging ...')
        artia.mpi.run(opts.STA, 4, avgCfgFile, 'runRemote', opts.remote, 'remoteHost', opts.host, 'suppressOutput', false, 'execDir', opts.STA_dir)
        
        % Do the refinement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Doing the refinement ...')
        artia.mpi.run(opts.EmSARTRefine, 1, refCfgFile, 'runRemote', opts.remote, 'remoteHost', opts.host, 'suppressOutput', false);
        
        % Project the particles with the shift %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Backprojecting particles with shift ...')
        % Get shifts
        shifts = artia.em.read(shiftFile);
        % Get aligned particles
        motive = artia.em.read(refMotlFile);
        motive(8:10, :) = motive(8:10, :) + motive(11:13, :);
        motive(11:13, :) = 0;
        
        % Get marker
        mFile = artia.marker.read(currMarkFile);
        marker = mFile.ali;
        phi = mFile.beamDeclination;
        
        [newMarker, model] = artia.geo.project3D(motive(8:10, :), marker, phi, opts.reconDim, opts.imDim, opts.voxelSize, opts.volumeShifts, opts.maAmount, opts.maAngle, shifts);
        
        mFile.ali = newMarker;
        mFile.model = model;
        mFile.markerCount = size(motive, 2);
        
        newMarkFile = sprintf([iterDir 'tomo_%d_new_marker_unaligned_iter_%d.em'], opts.tomoNr, i);
        alignedMarkFile = sprintf([iterDir 'tomo_%d_new_marker_aligned_iter_%d.em'], opts.tomoNr, i);
        artia.marker.write(mFile, newMarkFile);
        
        % Run the alignment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Running the 3D alignment ...')
        % Figure out the central marker
        positions = opts.motl(8:10, :) + opts.motl(11:13, :);
        centroid = mean(positions, 2);
        dist = sqrt((positions(1, :) - centroid(1)).^2 + (positions(2, :) - centroid(2)).^2 + (positions(3, :) - centroid(3)).^2);
        [~, idx] = min(dist);
        
        % Run first alignment for beam declination
        format1 = '%s -i %s -o %s -a %d -w %d -h %d --maAmount %f --maAngle %f --phi %f --iter 20 --iterSwitch 10 --doPsi --doFixedPsi --doPhi\n';
        com1 = sprintf(format1, opts.cAligner, newMarkFile, alignedMarkFile, idx-1, opts.imDim(1), opts.imDim(2), opts.maAmount, opts.maAngle, mFile.beamDeclination);
        system(com1);
        
        % Load file to get beam declination
        mFile = artia.marker.read(alignedMarkFile);
        
        % Run second alignment for tilt angles + image roation
        format2 = '%s -i %s -o %s -a %d -w %d -h %d --maAmount %f --maAngle %f --phi %f --iter 20 --iterSwitch 10 --doPsi --doTheta\n';
        com2 = sprintf(format2, opts.cAligner, newMarkFile, alignedMarkFile, idx-1, opts.imDim(1), opts.imDim(2), opts.maAmount, opts.maAngle, mFile.beamDeclination);
        system(com2);
        
        % Retrieve the new coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Making new motivelist ...')
        mark = artia.marker.read(alignedMarkFile);
        alignedMarker = mark.ali;
        coords = mark.model;
        phi = mark.beamDeclination;
    
        % Bring in correct orientation
        mx = -coords(2:3:end);
        my = coords(1:3:end);
        mz = coords(3:3:end);
        
        % Compute new volume shifts to center the particles in the volume
        opts.volumeShifts = round([-mean([max(mx), min(mx)]) -mean([max(my), min(my)]) -mean([max(mz), min(mz)])]);
        
        % Correct binning
        x = (mx + opts.volumeShifts(1))./opts.voxelSize;
        y = (my + opts.volumeShifts(2))./opts.voxelSize;
        z = (mz + opts.volumeShifts(3))./opts.voxelSize;
        
        % Volume size computation (max/min coord + opts.borderSize)
        xr = 4*floor((round(max(x) + abs(min(x))) + opts.borderSize)/4);
        yr = 4*floor((round(max(y) + abs(min(y))) + opts.borderSize)/4);
        zr = 4*floor((round(max(z) + abs(min(z))) + opts.borderSize)/4);
        
        % Volume size divisible by 2
        if mod(xr, 2) ~= 0
            xr = xr+1;
        end
        if mod(yr, 2) ~= 0
            yr = yr+1;
        end
        if mod(zr, 2) ~= 0
            zr = zr+1;
        end
        
        % New volume dimensions
        newreconDim = [xr, yr, zr];
        
        % Now compute the final particle coords in the new volume
        % dimensions
        x = x + newreconDim(1)/2 + 0.5;
        y = y + newreconDim(2)/2 + 0.5;
        z = z + newreconDim(3)/2 + 0.5;

        % Assign everything to the particle list and save
        new_motl = zeros(20, size(alignedMarker, 3));
        new_motl(8, :) = round(x);
        new_motl(9, :) = round(y);
        new_motl(10, :) = round(z);
        new_motl(5, :) = opts.tomoNr;
        new_motl(6, :) = 1:size(alignedMarker, 3);
        
        newMotiveFile = sprintf([iterDir 'tomo_%d_new_motl_iter_%d.em'], opts.tomoNr, i);
        artia.em.write(new_motl, newMotiveFile);
        
        % Hand over stuff to the next iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Handing over ...')
        opts.motl = new_motl;
        opts.markerFile = alignedMarkFile;
        opts.reconDim = newreconDim;
        
        disp('Done ...')
    end
end