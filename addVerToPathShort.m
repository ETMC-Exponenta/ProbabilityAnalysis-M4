
function addVerToPathShort

rootDir = fileparts(which('addVerToPathShort'));

path2add{1} = [rootDir,'\Ver\AllSV'];
path2add{2} = [rootDir,'\Ver'];
path2add{3} = [rootDir,'\NewLect'];
path2add{4} = [rootDir,'\Geo\Sh3'];
path2add{5} = [rootDir,'\Geo'];
path2add{6} = [rootDir,'\My_Tools'];
path2add{7} = [rootDir,'\WorkProb'];
path2add{8} = [rootDir,'\Method'];


platform = computer;
if (~any(strcmp(platform,{'PCWIN','PCWIN64'})))
    for i = 1:length(path2add)
        path2add{i}(path2add{i} == '\') = '/';
    end
end

for i = 1:length(path2add)
    addpath(path2add{i},'-begin');
end

