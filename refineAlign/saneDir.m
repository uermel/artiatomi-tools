function dirName = saneDir(dirName)
    if ~(dirName(end) == '/')
        dirName = [dirName '/'];
    end
end