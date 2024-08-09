function [phi,phi_cheb] = get_kernel_from_c( c,basis,support )

%
% function [phi,phi_cheb] = get_kernel_from_c( c,basis,support )
%
% Returns a function pointer to the interaction kernel with coefficients c on the dictionary basis{:}.
% It then approximates is by a Chebfun on the specified support (or [0,10] if support not specified). If the package chebfun is
% not available, it issues a warning and lets phi_cheb = phi

% (c) Mauro Maggioni

phi  = @(x) 0*x;                                                                                                                % The true kernel function handle
n = length(basis);
for i = 1:n
    phi = @(x) phi(x) + c(i)*basis{i}(x);
end

% Maybe nargout = 1 implies chebfun is not used and
% ~exist("chebfun") implies chebfun not found 
% (I installed chebfun and it still says not found made me confused for a while) 

if nargout <= 1
    fprintf('chebfun not used')
else
    if exist("chebfun")
        % if nargout>1 && exist("chebfun")
        if nargin<3 || isempty(support),    support = [0,10];   end                                                                 %MM:TBD: should be handled properly
        try
            phi_cheb = chebfun( {phi,0},[support,Inf],'eps',1e-6,'splitting','on' );
        catch
            warning('chebfun failed')
            phi_cheb = phi;
        end
    else
        warning('chebfun not found')
        phi_cheb = phi;
    end
end

end
