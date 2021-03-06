%% microFJ gives the function to zero for finding equilibrium
% sys1 - reference state
% sys2 - reference state for phase shift condition
% M - number of cars (M must divide N)
% N - discretization points
% L - length of road
% tau - momentum constant;
% Returns:
% F - the function to zero for an equilibrium solution
function [F,J, L] = microFJ(sys1,sys2,M,N,L, v0, tau, h)

k = N/M;
mu = sys1(end);
c = sys1(end - 1);     % wavespeed

u = sys1(1:N);
u0 = sys2(1:N);

%% operators & differentiation matrices

D  = fourdif(N,1)*2*pi/M;
D2 = fourdif(N,2)*(2*pi/M)^2;

e0    = ones(N,1);
shift = sparse([[N-k+1:N] [1:(N-k)]],1:N,e0,N,N);                               % shift matrix
e = sparse(1:N,1:N,e0,N,N); 

LN = -c^2*tau*D2+c*D;                              % linear comp. of DE
LD = D;                                                                %  phase condition matrix
ln_consv = sparse(ones(M,1),k*[1:M],ones(M,1),1,N);     % conservation law vector

%% Function
F = [ LN*u + optimalVelocity(h, shift*u, v0) - optimalVelocity(h, u, v0) + mu.*ones(N,1);
    ln_consv*u - L; ...
    (LD*u)'*(u-u0)...
    ];

%% Jacobian computation
if nargout > 1
    J1 = diag(dOptimalVelocity(e*u,v0));
    J2 = diag(dOptimalVelocity(shift*e*u,v0))*shift;
    Jc = (-2*c*tau*D2+D)*u;
    Jmu = ones(N,1);
    
    J  = sparse(...
        [LN-J1+J2, Jc, Jmu;...
        ln_consv, 0, 0;...
        (LD*u)'+ u'*LD - u0'*LD, 0, 0]...
        );
    
    % derivative in v0
    Jv = [ tanh(shift*u-h) - tanh(u-h) ; 0 ; 0];
    J = [J Jv];
    
    % linearization
    if nargout > 2
        L = LN+J1-J2;
    end
end


% optimal velocity derivative
    function v = dOptimalVelocity(headway,v0)
        v = v0*(1-tanh(headway - h).^2);
    end


end