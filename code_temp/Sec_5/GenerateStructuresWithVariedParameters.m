function structs = GenerateStructuresWithVariedParameters( params )

% IN:
%   params : cell array of structures with the following fields:
%                   Name    : a string with the name of a field
%                   Values  : cell array with the possible values taken by the field
%

%
% Mauro Maggioni
% @ Copyright, 2006
%

% Create empty structure
for lk = 1:length(params)
    lCurStruct{1}.(params{lk}.Name) = params{lk}.Values{1};
end
lIdx = 1;
% Create structures with all possible arrangements of parameters
for lk = 1:length(params)
    for lj = 1:length(params{lk}.Values)
        if (lk>1) && (lj==1)
            continue;
        end
        for li = 1:length(lCurStruct)
            lNewStruct{lIdx} = lCurStruct{li};
            lNewStruct{lIdx}.(params{lk}.Name) = params{lk}.Values{lj};
            lIdx = lIdx + 1;
        end
    end
    lCurStruct = lNewStruct;
end

structs = lCurStruct;

return