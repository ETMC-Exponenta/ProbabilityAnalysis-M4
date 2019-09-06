function generate()
p = parameters();
fs = struct2table(dir(p.source));
fs = fs(~ismember(fs.name, {'.' '..' '.git' '.gitignore'}), :);
for i = 1 : height(fs)
    fpath = fullfile(fs.folder{i}, fs.name{i});
    if fs.isdir(i)
        copydir(fpath, pwd);
    else
        copyfile(fpath, pwd);
    end
end
dirs = fs(fs.isdir, :);
for i = 1 : height(dirs)
    d = dirs.name{i};
    pcode(d, '-inplace');
    delete(fullfile(d, '*.m'));
end
end


function copydir(dpath, target)
    [~, dname] = fileparts(dpath);
    targetpath = fullfile(target, dname);
    if ~isfolder(targetpath)
        mkdir(targetpath);
    end
    copyfile(dpath, targetpath);
end