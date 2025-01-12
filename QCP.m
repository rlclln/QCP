
function X_est = QCP(AnchorPos,r,NL_err,sigma)
% =======================================
% INPUT: 
% AnchorPos: the position of anchors, which size is dim*N, i.e., [dim,N] = size(AnchorPos).
% r: range measurements with measurement noise and possible NLOS error e_i, that is: 
% r_i =norm(AnchorPos_i - x) + e_i + n_i, for i = 1, ...,N. See Equation (3) in our paper.
% sigma£º standard deviation of measurement noise. The size is N*1 or 1*1 which depends on demand,
% that is,N*1 for experiment or 1*1 for i.i.d.
% NL_err: the maximum value of NLOS error, a positive number .
% OUTPUT:
% X_est_nlos, position estimation of target node.


[dim,N] = size(AnchorPos);
% w = (1 - r./sum(r,1)); % weight factors
w = 1;


Phi = [-2.*AnchorPos' ones(N,1)];
square_ans2org = sum(AnchorPos'.^2,2);
sgm = 1*sigma;
g = [(r + sgm).^2 - sum(AnchorPos'.^2,2)];
r_wave = r.^2;
index = find(NL_err >= 2.*r);
if sum(index)
    r_wave(index) = (r(index) - NL_err).^2;
end

cvx_clear
cvx_begin sdp
cvx_solver SeDuMi
cvx_quiet (1)
cvx_precision best
variable x(dim,1) 
variable h(N) nonnegative
variable d(N) nonnegative
variable v(1) nonnegative
variable u(N) nonnegative

    minimize (sum(h))
    subject to
        norm([sqrt(w).*(r_wave-d)],2) <= geo_mean([2.*d,sigma.^2.*h.*ones(N,1)]); 
        Phi*[x;v] <= g;
        d == Phi*[x;v] + square_ans2org;
        {x,1,v} == rotated_lorentz(3); 
cvx_end

if isnan(cvx_optval) || isinf(cvx_optval)
    cvx_clear
    cvx_begin sdp
    cvx_solver SeDuMi
    cvx_quiet (1)
    cvx_precision best
    variable x(dim,1) 
    variable h(N) nonnegative
    variable d(N) nonnegative
    variable v(1) nonnegative
    variable u(N) nonnegative
        minimize (sum(h))
        subject to
            norm([sqrt(w).*(r_wave-d)],2) <= geo_mean([2.*d,sigma.^2.*h.*ones(N,1)]);
            sum(Phi,1)*[x;v] <= sum(g);
            d == Phi*[x;v] + square_ans2org; 
            {x,1,v} == rotated_lorentz(3); 

    cvx_end
end

X_est =x;

end
