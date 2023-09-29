function [vp, rho, eta, lambda] = compute_vp(F, data, sz, varargin)
%% Method using MATLAB lsqr function
% Can set various regularization, with a more flexible framework (building
% the actual A matrix to inverse.
% use compute_vp for a more typical setting
p = inputParser;
addParameter(p,'order_derivative_smoothness', 2);
addParameter(p,'lambda', nan); % reg param same as smoothness in compute_vp_flex
% if nan (default), find the L-curve value
parse(p,varargin{:});

% get non-zero data
[D, id_nnz] = format_data(data(:,1));

% Build the matrix linking a pixel in 3D to its vertical index
U1 = repelem(speye(sz(3)),prod(sz(1:2)),1);

% Build the inverse system
A = F(id_nnz,:)*U1;
b = D;

% Add smoothness regularization
% Lh = LaplacianMatrix1D(sz(3)); Lh = Lh - diag(sum(Lh,2)); old approach
% building the full Lh matrix
L = get_l(sz(3),p.Results.order_derivative_smoothness);

if isnan(p.Results.lambda)
    [U,s,V] = csvd(A);
    reg_corner = l_curve_wo_fig(U,s,b,'mtsv',L,V);
    lambda = reg_corner;
else
    lambda = p.Results.lambda;
end

[U,sm,X,~,~] = cgsvd(A,L);

[vp, rho, eta] = tikhonov(U,sm,X,b,lambda);


end