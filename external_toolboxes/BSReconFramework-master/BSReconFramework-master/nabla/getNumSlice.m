function numSlice = getNumSlice(pathName)

currentPath = pwd;
cd(pathName);

info = dir([pwd,'/*.DCM']);
list = {info.name};

%list = ls;
scan = zeros(size(list,1),1);
slice = zeros(size(list,1),1);

for iList =1:length(list)
    scan(iList) = str2double( list{iList}(11:12) );
    slice (iList) = str2double( list{iList}(18:20) );
end

numSlice = max(slice);
cd(currentPath);

end