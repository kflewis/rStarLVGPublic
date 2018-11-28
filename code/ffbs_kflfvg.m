function [LL,LLt,BSDraw,m,V,FFDraw] = ffbs_kflfvg(y,u,A,B,C,D,F,G,S,m0,V0,bsInd)
%% function [LL,LLt,S,m,V] = ffbs(y,u,A,B,C,D,F,G,S,m0,V0)
%% A function to perform forward-filter-backward-sampling on a state space model of the
%% form:
%%
%%     x_t = Ax_t-1 + B u_t + C w_t
%%     y_t = D x_t + F u_t + G w_t
%%     where w_t ~ i.i.d. N(0,S)
%%     Given that we have an initial estimate of the state m0, and the varaince of that
%%     estimate, V0.  
%%
%% This function will require that all the data (observed y and observed or exogenous u)
%% be given to the function, as well as the properly constructed matrices and initial
%% state estimates and corresponding variances.  The final input variable, bsInd, is an
%% indicator for whether or not we will conduct backward sampling of the state variable
%% or if only the forward filtering will occur (0 means only forward-filtering, 1 means
%% both filtering and backward sampling, the default will be 0).  This may be avoided in order to speed
%% up things like Metropolis-Hastings draws from the parameters where the likelihood
%% value is the most crucial element.
%%
%% We will assume that the y variable will be (M X T), u will be (N X T), xt and m0 will
%% be q X 1 but the matrix of state observations will be (q X T) and Vt and V0 will be
%% (q X q) meaning that the array of state variances will be (q X q X T).
%%
%% Kurt Lewis
%% Created: April 20, 2016
%% Last-updated Time-stamp: <2018-11-09 16:07:57 (m1kfl00)>
%%
    
% Minimal Idiot-proofing the inputs 
if nargin < 11
    help ffbs_kflfvg
    error('You must give the program data and a full set of matrices and initial conditions');
end

if nargin < 12
    % If they don't specify a value for backward sampling, assume you don't do it.
    bsInd = 0 ;
end

if nargin < 13
    % If they don't ask specifically for the check of V convergence, don't
    VconFlag = 0;
end

% Get some indications about the data and inputs
[M,T] = size(y);
[N, T1] = size(u);

% Throw an error if the number of observations for y and u are different
if T ~= T1
    error('Observations of y and u are different quantities, must fix.')
end

% Compare the shapes of m0 and V0, throw error if not compatible
[q,T2] = size(m0);
[T3,T4] = size(V0);

if T3 ~= T4
    error('V0 is not square')
end

if T3 ~= q
    error('m0 and V0 are not of compatible sizes')
end

% Given the size of the observed data and the needs of the filter, we declare some
% variables we will need to populate

% state estimates and covariance
m = NaN(q,T);
V = NaN(q,q,T);
% The log-likelihood at each point in time
LLt = NaN(T,1);
% The draw of the state at each point in time, to be filled in during backward-sampling
BSDraw = NaN(q,T);

% As we construct matrices which are the sum of other parts, we need to account for
% issues of symmetry that crop up when we add a larger number of symmetric matrices
% together that may be slightly asymmetric due to numerical issues.  We need a
% tolerance for this numerical precision, so we set here a value which will correspond
% to the maximum allowable asymmetry
symTol = 1e-8;


%%% FORWARD FILTER 

%% The first step uses m0 and V0

% Define the stuff from the joint Gaussian distribution
mux1 = A*m0 + B*u(:,1);
muy1 = D*A*m0 + (D*B + F)*u(:,1);
SigXX1 = A*V0*A' + C*S*C';
SigXY1 = A*V0*A'*D' + (C*S)*((D*C) + G)';
SigYY1 = D*A*V0*A'*D' + ((D*C) + G)*S*((D*C) + G)';

% SigXX and SigYY must be symmetric, and we will make them that way by setting them
% equal to the sum of the matrix and the traspose divided by 2.
maxDiffX = max(max(abs(SigXX1 - SigXX1')));
maxDiffY = max(max(abs(SigYY1 - SigYY1')));

if maxDiffX < symTol
    SigXX1 = 0.5*(SigXX1 + SigXX1');
else
    error('SigXX1 is too asymmetric, this happens occasionally, start again with a new rng seed')
end

if maxDiffY < symTol
    SigYY1 = 0.5*(SigYY1 + SigYY1');
else
    error('SigYY1 is too asymmetric, this happens occasionally, start again with a new rng seed')
end

% Define the mean and variance from the resulting conditional distribution
m(:,1) = mux1 + SigXY1*(inv(SigYY1))*(y(:,1) - muy1);
V(:,:,1) = SigXX1 - (SigXY1*(inv(SigYY1))*SigXY1');

[VV,DD] = eig(V(:,:,1));

if (min(DD) < -1e-8)
    error('There is an eigenvalue in V1 that is TOO negative, this happens occasionally, start again with a new rng seed');
else
    V(:,:,1) = VV*abs(DD)*VV';
    clear VV DD
end

FFDraw(:,1) = mvnrnd(m(:,1),V(:,:,1));

% Fill in the log-likelihood value for the first step
LLt(1,1) = -(M/2)*log(2*pi) - 0.5*log(det(SigYY1)) ...
    - 0.5*((y(:,1) - muy1)'*(inv(SigYY1))*(y(:,1) - muy1));


%% Steps 2 through T of the forward filter
for i = 2:T
    
    % Define the stuff from the joint Gaussian distribution
    mux = A*m(:,(i-1)) + B*u(:,i);
    muy = D*A*m(:,(i-1)) + (D*B + F)*u(:,i);
    SigXX = A*V(:,:,(i-1))*A' + C*S*C';
    SigXY = A*V(:,:,(i-1))*A'*D' + (C*S)*((D*C) + G)';
    SigYY = D*A*V(:,:,(i-1))*A'*D' + ((D*C) + G)*S*((D*C) + G)';

    % SigXX and SigYY must be symmetric, and we will make them that way by setting them
    % equal to the sum of the matrix and the traspose divided by 2.
    maxDiffX = max(max(abs(SigXX - SigXX')));
    maxDiffY = max(max(abs(SigYY - SigYY')));
     
    if maxDiffX < symTol
        SigXX = 0.5*(SigXX + SigXX');
    else
        error('SigXX is too asymmetric, this happens occasionally, start again with a new rng seed')
    end
     
    if maxDiffY < symTol
        SigYY = 0.5*(SigYY + SigYY');
    else
        error('SigYY is too asymmetric, this happens occasionally, start again with a new rng seed')
    end
    
    % Define the mean and variance from the resulting conditional distribution
    m(:,i) = mux + SigXY*(inv(SigYY))*(y(:,i) - muy);
    V(:,:,i) = SigXX - (SigXY*(inv(SigYY))*SigXY');
    
    [VV,DD] = eig(V(:,:,i));
    
    if (min(DD) < -1e-8)
        error('There is an eigenvalue in Vt that is TOO negative, this happens occasionally, start again with a new rng seed');
    else 
        V(:,:,i) = VV*abs(DD)*VV';
        clear DD VV
    end
        
    FFDraw(:,i) = mvnrnd(m(:,i),V(:,:,i));
    
    % Fill in the log-likelihood value for the step
    LLt(i,1) = -(M/2)*log(2*pi) - 0.5*log(det(SigYY)) ...
    - 0.5*((y(:,i) - muy)'*(inv(SigYY))*(y(:,i) - muy));
    
    clear mux muy SigXX SigXY SigYY maxDiffX maxDiffY
    
end

% Find the total log-likelihood
LL = sum(LLt);


%%% BACKWARD SAMPLER

%% Only proceed if the flag for this part of the code is operative

if bsInd == 1
    
    % First, sample from the T state distribution
    
    % Start by error checking the covariance matrix
    VT = V(:,:,T);
    maxDiffV = max(max(abs(VT - VT')));
    if maxDiffV < symTol
        VT = 0.5*(VT + VT');
    else
        error('VT is TOO asymmetric, this happens occasionally, start again with a new rng seed')
    end
    
    [VV,DD] = eig(VT);
    
    if (min(DD) < -1e-8)
        error('There is an eigenvalue in VT that is TOO negative, this happens occasionally, start again with a new rng seed');
    else
        VT = VV*abs(DD)*VV';
    end

    % Don't resample for the backward sample, we have a T given T draw already
    %BSDraw(:,T) = mvnrnd(m(:,T),VT);
    BSDraw(:,T) = FFDraw(:,T);
    
    % Looping from T-1 to 1, fill in the rest of the draws for the sample of the state
    for j = (T-1):-1:1
        
        [VV,DD] = eig(A*V(:,:,j)*A' + C*S*C') ;
        DDinv = diag(1./diag(DD)) ;
        
        S1 = (V(:,:,j)*A')*VV*DDinv*VV' ;
        mu1 = m(:,j) + S1*(BSDraw(:,(j+1)) - (A*m(:,j) + B*u(:,j+1))) ;
        Sig1 = V(:,:,j) - S1*A*V(:,:,j) ;
        Sig1 = Sig1/2+Sig1'/2 ;
        [VV,DD] = eig(Sig1) ;
        
        if (min(DD) < -1e-8)
            error('There is an eigenvalue in Sig1 that is TOO negative, this happens occasionally, start again with a new rng seed');
        else
            Sig1 = VV*abs(DD)*VV';
            clear VV DD
        end
                
        BSDraw(:,j) = mvnrnd(mu1,Sig1);
        
        clear mu1 Sig1
    end
end

    


    
    



    
    
