%% for generate the interface defination file 
%% of ClassCalculateA20InOneElement.hpp 
%% and ClassCalculateA20InOneElement.cpp
%% labrary

clear
clc
productPath='./';
hppFileMy='ClassCalculateA2lijInOneElement.hpp';
cppFileMy='ClassCalculateA2lijInOneElement.cpp';

clibgen.generateLibraryDefinition(fullfile(productPath,hppFileMy),...
"SupportingSourceFiles",fullfile(productPath,cppFileMy),...
"IncludePath",productPath,...
"ReturnCArrays",false,...
'Verbose',true);
