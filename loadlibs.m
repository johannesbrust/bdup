function loadlibs
%LOADLIBS Loads the Matlab interface when available
archs = computer('arch');
if strcmp(archs,'glnxa64')
    archl = 'unix_x64';
elseif strcmp(archs,'maci64')
    archl = 'osx_x64';
else
    fprintf(['There aren''t precompiled interface functions for \n',...
               'architecture:',archs,' and Matlab. \n',...
               'Consider recompiling the interface yourself.']);
    return;
end

ssp     = strsplit(pwd,'/');
issp    = find(contains(ssp,'bdup'));
ssj     = join(ssp(1:issp),'/');

addpath(genpath([ssj{:},'/','FORTRAN/','objs/',archl]));

end