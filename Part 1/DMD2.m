function [Phi,omega,lambda,b,Xdmd,S] = DMD2(X1,X2,r,t)
%Computes the Dynamic Mode Decomposition of X1, X2

%INPUTS:
% X1 = X, data matrix
% X2 = X', shifted data matrix
% Columns of X1 and X2 are state snapshots
% r = target rank of SVD
% dt = time step advancing X1 to X2
% t is input time

% OUTPUTS:
% Phi, the DMD modes
% omega, the continuous-time DMD eigenvalues
% lambda, the discrete-time DMD eigenvalues
% b, a vector of magnitudes of modes Phi
% Xdmd, the data matrix reconstructed by Phi, omega, b

%%%% Step 1: SVD the raw data X, and get the reduced space
dt = t(2) - t(1);

[U,S,V] = svd(X1,'econ');
r = min(r,size(U,2));


%semilogy(diag(S),'ko')

% truncate to rank-r
U_r = U(:,1:r);
S_r = S(1:r,1:r);
V_r = V(:,1:r);

%%%% Step2: Find Atilde, low-rank dynamics
Atilde = U_r'*X2*V_r/S_r;

%%%% Step3: Find eigenvalue, eigenvector of Atilde
[W_r,D] = eig(Atilde);

%%% Step4: Find DMD modes and initial value

%DMD modes
Phi = X2*V_r/S_r*W_r;

% discrete-time eigenvalues
lambda = diag(D);

% continuous-time eigenvalues
omega = log(lambda)/dt;


% Compute DMD mode amplitudes b
x1 = X1(:,1);
b = Phi\x1;

% DMD reconstruction
% 
% mm1 = size(X1, 2);
% time_dynamics = zeros(r,mm1);
% 
% % time vector
% t = (0:mm1-1)*dt;
% 
% for iter = 1:mm1
%     time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
% end
% 
% Xdmd = Phi * time_dynamics;

time_dynamics = zeros(r,length(t));

for iter = 1:length(t)
    time_dynamics(:,iter) = b.*exp(omega*(t(iter)));
end

Xdmd = Phi*time_dynamics;

end

