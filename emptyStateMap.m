function mapNullCells = emptyStateMap(size)
if nargin == 0
    size = 3; % default and minimum: for passive, potentiated, and active default states
end
if size < 3
    disp('emptyStateMap:  Setup error. Size argument must be at least 3, for passive, potentiated, and active states')
    return
end
mapNull = zeros(size,2);
for i = 2:size
    mapNull(i,1)=i-1;
end
mapNullCells = {mapNull}; % state maps are always cell arrays, in the event there are more than one per variable
