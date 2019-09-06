function generate()
% Load project parameters
p = parameters();
% Get source files
fs = struct2table(dir(p.source));
fs = fs(~ismember(fs.name, {'.' '..' '.git' '.gitignore'}), :);
% Copy source files
disp('Copying files...');
for i = 1 : height(fs)
    fpath = fullfile(fs.folder{i}, fs.name{i});
    if fs.isdir(i)
        copydir(fpath, pwd);
    else
        copyfile(fpath, pwd);
    end
end
% Create p-files
disp('p-coding files...');
dirs = fs(fs.isdir, :);
for i = 1 : height(dirs)
    d = dirs.name{i};
    mfiles = struct2table(dir(fullfile(d, '**/*.m')));
    for j = 1 : height(mfiles)
        try
            f = fullfile(mfiles.folder{j}, mfiles.name{j});
            pcode(f, '-inplace');
            delete(f);
        catch err
            warning('Unable to pcode: %s', f);
        end
    end
end
disp('All done!');
end


function copydir(dpath, target)
% Copy directory to target dir
[~, dname] = fileparts(dpath);
targetpath = fullfile(target, dname);
if ~isfolder(targetpath)
    mkdir(targetpath);
end
copyfile(dpath, targetpath);
end