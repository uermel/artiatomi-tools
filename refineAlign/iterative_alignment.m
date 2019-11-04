function iterative_alignment(opts)
    
    % Make project dir
    opts.projDir = saneDir(opts.projDir);
    mkdir(opts.projDir);
    
    for i = 1:opts.iters
        % Make directories and names and save files %%%%%%%%%%%%%%%%%%%%%
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
        
        % Make configs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        oss1_filtered.fourFilterLP = '200';
        oss1_filtered.fourFilterLPS = '50';
        oss1_filtered.fourFilterHP = '10';
        oss1_filtered.fourFilterHPS = '5';
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
        refine.fourFilterLP = '200';
        refine.fourFilterLPS = '50';
        refine.fourFilterHP = '10';
        refine.fourFilterHPS = '5';
        refine.GroupSize = num2str(opts.groupSize);
        refine.MaxShift = num2str(opts.maxShift);
        refine.Reference = refFile;
        refine.SizeSubVol = num2str(opts.boxSize);
        refine.SpeedUpDistance = num2str(opts.speedUpDist);
        refine.VoxelSizeSubVol = num2str(opts.voxelSize);
        refine.CCMapFileName = ccmapPre;
        refine.MotiveList = refMotlFile;
        refine.ShiftOutputFile = shiftFile;
        refine.MultiPeakDetection = false;
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
        avg.AngIter = '3';
        avg.AngIncr = '5';
        avg.PhiAngIter = '3';
        avg.PhiAngIncr = '5';
        avg.LowPass = '10';
        avg.HighPass = '1';
        avg.Sigma = '2';
        artia.cfg.write(avg, avgCfgFile);
        
        % Run the reconstructions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Running reconstructions ...')
        % Write the run script
        %runScript = [iterDir 'run_script.sh'];
        %fileID = fopen(runScript,'w');
        %fprintf(fileID,'%s\n','cd /home/Group/Software/tomography/kunzFunctions/EmSART/');
        %com = 'mpiexec -n 4 EmSART -u ';
        %fprintf(fileID,'%s%s\n %s%s\n %s%s\n', com, oss10CfgFile, ...
                                               %com, oss1CfgFile, ...
                                               %com, oss1FilteredCfgFile);                                        
        %fclose(fileID);
        
        % Execute
        %[~, ~] = system(['chmod +x '  runScript]);
        %system(['ssh -t borg ' runScript]);
        artia.mpi.run(opts.EmSART, 4, oss10CfgFile, 'runRemote', opts.remote, 'remoteHost', opts.host, 'suppressOutput', false)
        artia.mpi.run(opts.EmSART, 4, oss1CfgFile, 'runRemote', opts.remote, 'remoteHost', opts.host, 'suppressOutput', false)
        artia.mpi.run(opts.EmSART, 4, oss1FilteredCfgFile, 'runRemote', opts.remote, 'remoteHost', opts.host, 'suppressOutput', false)
        
        % Get Parts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tomolist = {};
        tomolist{opts.tomoNr} = oss1ReconFile;
        artia.particle.extract_write(extractMotlFile, 1, tomolist, opts.boxSize/2, 0, 0, 1, partFilePre);
        %writeParts(extractMotlFile, parts, );
        
        
        % Run the averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Running averaging ...')
        % Write the run script
        %runScript = [iterDir 'run_script.sh'];
        %fileID = fopen(runScript,'w');
        %fprintf(fileID,'%s\n','cd /home/kunz/DevelopSicherung/SubTomogramAverageMPI/SubTomogramAverageMPI/SubTomogramAverageMPI');
        %fprintf(fileID,'%s%s\n','mpiexec -n 4 SubTomogramAverageMPI -u ', avgCfgFile);
        %fclose(fileID);

        % Run the command
        %[~, ~] = system(['chmod +x '  runScript]);
        %system(['ssh -t borg ' runScript]);
        artia.mpi.run(opts.STA, 4, avgCfgFile, 'runRemote', opts.remote, 'remoteHost', opts.host, 'suppressOutput', false, 'execDir', opts.STA_dir)
        
        % Do the refinement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Doing the refinement ...')
        % Write the run script
        %runScript = [iterDir 'run_script.sh'];
        %fileID = fopen(runScript,'w');
        %fprintf(fileID,'%s\n','cd /home/Group/Software/tomography/kunzFunctions/EmSART/');
        %com = 'EmSARTRefine -u ';
        %fprintf(fileID,'%s%s\n', com, refCfgFile);
        %fclose(fileID);
        
        % Execute
        %[~, ~] = system(['chmod +x '  runScript]);
        %system(['ssh -t borg ' runScript]);
        artia.mpi.run(opts.EmSARTRefine, 1, refCfgFile, 'runRemote', opts.remote, 'remoteHost', opts.host, 'suppressOutput', false);
        
        % Project the particles with the shift %%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Backprojecting particles with shift ...')
        shifts = artia.em.read(shiftFile);
        motive = artia.em.read(refMotlFile);
        motive(8:10, :) = motive(8:10, :) + motive(11:13, :);
        motive(11:13, :) = 0;
        %[marker, model, phi] = emreadCoords(currMarkFile);
        mFile = artia.marker.read(currMarkFile);
        marker = mFile.ali;
        phi = mFile.beamDeclination;
        %projCount = size(marker, 2);

        %exclude = zeros(projCount, 1);
        %for j = 1:projCount
        %    if marker(2, j, 1) < 0 || marker(3, j, 1) < 0
        %        exclude(j) = 1;
        %    end
        %end
        
        %newMarker = ConvertAmira3DCoordinatesToClicker2DCoordinatesShift(motive, marker, opts.reconDim, [7420 7676], opts.voxelSize, 0, shifts);
        [newMarker, model] = artia.geo.project3D(motive(8:10, :), marker, phi, opts.reconDim, opts.imDim, opts.voxelSize, opts.volumeShifts, 1.016, 42, shifts);
        
        %for j = 1:projCount
        %    if exclude(j) == 1
        %        newMarker(2, j, :) = -1;
        %        newMarker(3, j, :) = -1;
        %        newMarker(4:end, j, :) = 0;
        %    end
        %end
        mFile.ali = newMarker;
        mFile.model = model;
        mFile.markerCount = size(motive, 2);
        
        newMarkFile = sprintf([iterDir 'tomo_%d_new_marker_unaligned_iter_%d.em'], opts.tomoNr, i);
        alignedMarkFile = sprintf([iterDir 'tomo_%d_new_marker_aligned_iter_%d.em'], opts.tomoNr, i);
        artia.marker.write(mFile, newMarkFile);
        %emwrite_silent(newMarker, newMarkFile);
        
        % Run the alignment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        positions = opts.motl(8:10, :) + opts.motl(11:13, :);
        centroid = mean(positions, 2);
        dist = sqrt((positions(1, :) - centroid(1)).^2 + (positions(2, :) - centroid(2)).^2 + (positions(3, :) - centroid(3)).^2);
        %diff = sum((positions - centroid).^2');
        [~, idx] = min(dist);
        
        disp('Running the 3D alignment ...')
        % Write the run script
        %runScript = [iterDir 'run_script.sh'];
        %fileID = fopen(runScript,'w');
        format = '%s -i %s -o %s -a %d -w %d -h %d --maAmount 1.016 --maAngle 42 --phi %f --iter 20 --iterSwitch 10 --doPsi --doTheta\n';
        com = sprintf(format, opts.cAligner, newMarkFile, alignedMarkFile, idx-1, opts.imDim(1), opts.imDim(2), mFile.beamDeclination);
        %fclose(fileID);
        
        % Execute
        %[~, ~] = system(['chmod +x '  runScript]);
        system(com);
        
        % Retrieve the new coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Making new motivelist ...')
        %[alignedMarker, coords, phi] = emreadCoords(alignedMarkFile);
        mark = artia.marker.read(alignedMarkFile);
        alignedMarker = mark.ali;
        coords = mark.model;
        phi = mark.beamDeclination;
    
        % Bring in correct orientation, binning and add half the recon dim
        mx = -coords(2:3:end);
        my = coords(1:3:end);
        mz = coords(3:3:end);
        
        opts.volumeShifts = [-mean([max(mx), min(mx)]) -mean([max(my), min(my)]) -mean([max(mz), min(mz)])];
        
        %x = mx./opts.voxelSize + 0.5 + opts.reconDim(1)/2;
        %y = my./opts.voxelSize + 0.5 + opts.reconDim(2)/2;
        %z = mz./opts.voxelSize;
        
        %newopts.reconDim = [opts.reconDim(1), opts.reconDim(2), round(std(z)) + 50];
        %z = z + newopts.reconDim(3)/2 + 0.5;
        

        %x = -coords(2:3:end)./opts.voxelSize + opts.reconDim(1)/2 + 0.5*opts.voxelSize;
        %y = coords(1:3:end)./opts.voxelSize + opts.reconDim(2)/2 + 0.5*opts.voxelSize;
        %z = coords(3:3:end)./opts.voxelSize;
        %newopts.reconDim = [opts.reconDim(1), opts.reconDim(2), round(std(z)) + 50];
        %z = z + newopts.reconDim(3)/2 + 0.5*opts.voxelSize;
        
        x = (mx + opts.volumeShifts(1))./opts.voxelSize;
        y = (my + opts.volumeShifts(2))./opts.voxelSize;
        z = (mz + opts.volumeShifts(3))./opts.voxelSize;
        
        xr = 4*floor((round(max(x) + abs(min(x))) + 2*opts.boxSize)/4);
        yr = 4*floor((round(max(y) + abs(min(y))) + 2*opts.boxSize)/4);
        zr = 4*floor((round(max(z) + abs(min(z))) + 2*opts.boxSize)/4);
        
        if mod(xr, 2) ~= 0
            xr = xr+1;
        end
        if mod(yr, 2) ~= 0
            yr = yr+1;
        end
        if mod(zr, 2) ~= 0
            zr = zr+1;
        end
        
        newreconDim = [xr, yr, zr];
        
        x = x + newreconDim(1)/2 + 0.5;
        y = y + newreconDim(2)/2 + 0.5;
        z = z + newreconDim(3)/2 + 0.5;

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
        %opts.volumeShifts = [0 0 0];
        
        disp('Done ...')
    end

end