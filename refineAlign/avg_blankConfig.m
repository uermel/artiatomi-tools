function s = avg_blankConfig()
    s = struct();
    s.CudaDeviceID = '0 1 2 3';
    s.MotiveList = '';
    s.WedgeFile = '';
    s.Particles = '';
    s.WedgeIndices = '';
    s.Classes = '';
    s.MultiReference = 'false';
    s.PathWin = '';
    s.PathLinux = '';
    s.Reference = '';
    s.Mask = '';
    s.MaskCC = '';
    s.NamingConvention = 'TomoParticle';
    s.StartIteration = '';
    s.EndIteration = '';
    s.AngIter = '0';
    s.AngIncr = '0';
    s.PhiAngIter = '0';
    s.PhiAngIncr = '0';
    s.LowPass = '0';
    s.HighPass = '0';
    s.Sigma = '0';
    s.ClearAngles = 'false';
    s.BestParticleRatio = '1';
    s.ApplySymmetry = 'transform';
    s.CouplePhiToPsi = 'true';
end