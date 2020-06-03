%%%%%%Develop a DMD model to forecast the future population states
clear all; close all; clc;
%%%%data x1--Snowshoe Hare   x2----Canada Lynx
x1 = [20,20,52,83,64,68,83,12,36,150,110,60,7,10,70,...
    100,92,70,10,11,137,137,18,22,52,83,18,10,9,65];
x2 = [32,50,12,10,13,36,15,12,6,6,65,70,40,9,20,...
    34,45,40,15,15,60,80,26,18,37,50,35,12,12,25];

slices = 30;
t = linspace(0,58,slices);
dt = t(2) - t(1);



% X1 = [x1(1:25);
%     x2(1:25);
%     x1(2:26);
%     x2(2:26);
%     x1(3:27);
%     x2(3:27);
%     x1(4:28);
%     x2(4:28);
%     x1(5:29);
%     x2(5:29);];
% 
% X2 = [x1(2:26);
%     x2(2:26);
%     x1(3:27);
%     x2(3:27);
%     x1(4:28);
%     x2(4:28);
%     x1(5:29);
%     x2(5:29);
%     x1(6:30);
%     x2(6:30)];

X1 = [];
X2 = [];

kk = 10;
for j=1:kk
    
    X1 = [X1;x1(j:29-kk+j);x2(j:29-kk+j)];
    X2 = [X2;x1(j+1:29-kk+j+1);x2(j+1:29-kk+j+1);];
    
end


r = 16;
%% 
%%%%%%%body of DMD

%%%% Step 1: SVD the raw data X, and get the reduced space
[U,Sigma,V] = svd(X1,'econ');
Ur = U(:,1:r);
Sigmar = Sigma(1:r,1:r);
Vr = V(:,1:r);

%%%% Step2: Find Atilde
Atilde = Ur'*X2*Vr/Sigmar;

%%%% Step3: Find eigenvalue, eigenvector of Atilde
[W,Lambda] = eig(Atilde);

%%% Step4: Find DMD modes and initial value

%%%% Phi are DMD modes, eigenvectors of A matrix
%%% Lambda are DMD eigenvalues (of A matrix)
%%% b is the mode amplitude

Phi = X2*(Vr/Sigmar)*W;

 alpha1 = Sigmar*Vr(1,:)';
 b = (W*Lambda)\alpha1;

%b = Phi\X1(:,1);
%%
%%%%% DMD reconstruction for every time point
%%% turn the mode to exp
Lambda_new = diag(Lambda);
omega = log(Lambda_new)/dt;

%y0 = Phi\X;
u_modes = zeros(r,length(t)); %DMD reconstruction for every time point
% x1_modes = zeros(1,length(t));
% x2_modes = zeros(1,length(t));

for iter = 1:length(t)
    
    u_modes(:,iter) = b.*exp(omega*(t(iter)));

end

u_dmd = Phi*u_modes;

f1 = figure();
subplot(2,2,1)
plot(t,x1,'r-',t,abs(u_dmd(1,:)),'b--');
legend('Hare','DMD Hare');

subplot(2,2,2)
plot(t,x2,'r-',t,abs(u_dmd(2,:)),'b--')
legend('Lynx','DMD Lynx');

subplot(2,2,3);
plot(diag(Sigma),'ko');
title('Sigma')

subplot(2,2,4)
plot(real(omega),imag(omega),'ko')
title('Omega')
xlabel('Real')
ylabel('Imagine')
%%

for j=1:length(t)
    error1(j) = norm(u_dmd(1,j)-x1(j));
    error2(j) = norm(u_dmd(2,j)-x2(j));
end
f2 = figure();
plot(t,error1,'r-',t,error2,'b--')
legend('error1','error2')


