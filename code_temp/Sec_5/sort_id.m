function new_id = sort_id(old_id)
% Since the index does not really match
% [1,2,2,1] and [2,1,1,2] indicates the same pattern
% So we assign the index of kernels based on the following rules. 
% 1. The kernel is of type 1
% 2. The next kernel which is different from type 1 is assigned type 2
% 3. The next kernel which is different from type {1,2} is assigned type 3
% ...

% a test seq
% old_id = [3,3,2,2,1,2,3,2,1];
% Q = 3;
N = length(old_id);
new_id = zeros(N, 1);
k = 1;
for i = 1:N
    if new_id(i) == 0
        pos = (old_id == old_id(i));
        new_id(pos) = k;
        k = k + 1;
    end
end