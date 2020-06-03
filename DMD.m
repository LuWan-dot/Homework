function [Phi,Lambda,b] = DMD(X,Xprime,r)
%%%%%  X is the raw data
%%%%   Xprime is the shifted data
%%%%%  r is the order of reduced space

%%%% Step 1: SVD the raw data X, and get the reduced space
[U,Sigma,V] = svd(X,'econ');
Ur = U(:,1:r);
Sigmar = Sigma(1:r,1:r);
Vr = V(:,1:r);

%%%% Step2: Find Atilde
Atilde = Ur'*Xprime*Vr/Sigmar;

%%%% Step3: Find eigenvalue, eigenvector of Atilde
[W,Lambda] = eig(Atilde);

%%% Step4: Find DMD modes and initial value

%%%% Phi are DMD modes, eigenvectors of A matrix
%%% Lambda are DMD eigenvalues (of A matrix)
%%% b is the mode amplitude

Phi = Xprime*(Vr/Sigmar)*W;

alpha1 = Sigmar*Vr(1,:)';
b = (W*Lambda)\alpha1;

end

