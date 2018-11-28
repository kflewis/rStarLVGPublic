function Lpi = priorEval_distribute(theta,estLoc)
% function Lpi = priorEvalTest1(theta,estLoc)
%
% THIS VERSION OF THE PRIOR EVALUATION FUNCTION TAKES ON BOARD CONSTRAINTS ON
% THE POSSIBLE VALUES ALLOWED BY THE IMPLIED THE LAMBDAS THAT RESULT FROM THE
% DISTRIBUTIONS OF THE SIGMAS.  WE IMPOSE CONSTRAINTS THAT LAMBDAS BE WELL
% BEHAVED BY THE MEASURE OF BEING BETWEEN 0 AND AN UPPER BOUND, LIKE 5 OR 10.
%
% This function takes in the theta vector and returns the log-prior value for the model
% so that we can combine it with the likelihood function to perform MH sampling.
%
% Kurt Lewis
% Created: April 27, 2016
% Last-updated Time-stamp: <2018-05-23 00:33:16 (m1kfl00)>
%
    
% First break out the elements of theta
a1 = theta(1);
a2 = theta(2);
ar = theta(3);
b1 = theta(4);
bY = theta(5);
s1 = theta(6);
s2 = theta(7);
s3 = theta(8);
s4 = theta(9);
s5 = theta(10);
aX = theta(11);
rhoG = theta(12);
rhoZ = theta(13);
muG = theta(14);
muZ = theta(15);

% Now that we have seperate entries, we can find the log-prior value for each variable,
% based on normal and inverse gamma distributions
% Most of the parameters priors will be mean zero and variance 1 for now
m = 0;
s = 2;

%% We are setting the priors on the sigmas to be uniform from 0 to 20, inclusive on
%% both ends.
sigLB = 0;
sigHB = 5;

% This will allow sig3HB to be larger than the others
sig3HB = sigHB;

% Put a floor on sig1
sig1LB = sigLB;

% Boundaries on lambdas
lambHB = 5;
lambLB = 0.01;

%%% We are no longer using the inverse gamma priors
%% To give the least informative prior, I will set both the shape and scale of the inverse
%% gamma to 1.
%sh = 1; 
%sc = 1;

%% From here, we will proceed through the theta vector one at a time to build the log
%% prior of the theta draw.

%% Theta(1): a1
%%
%% a1 has no restrictions so it will be evaluated, if estimated against the
%% untruncated normal distribution with mean m and standard deviation s.  If it is not
%% included in the estLoc vector, it will be given the degenerate prior of 1, which in
%% log-prior space is zero.
if (ismember(1,estLoc) == 1)
    l1 = log(normpdf(a1,m,s));
elseif (ismember(1,estLoc) == 0)
    l1 = 0;
end

%% Theta(2): a2
%%
%% a1 has no restrictions so it will be evaluated, if estimated against the
%% untruncated normal distribution with mean m and standard deviation s.  If it is not
%% included in the estLoc vector, it will be given the degenerate prior of 1, which in
%% log-prior space is zero.
if (ismember(2,estLoc) == 1)
    l2 = log(normpdf(a2,m,s));
elseif (ismember(2,estLoc) == 0)
    l2 = 0;
end

%% Theta(3): ar
%%
%% ar has the restriction that we inherit from (H)LW that it must be less than
%% -0.0025.  If it is estimated, it will be given log-prior pdf value associated with
%% the truncated distribution with mean m and standard deviation s.  If it is not
%% included in the estLoc vector, it will be given the degenerate prior of 1, which in
%% log-prior space is zero.
arVal = -0.0025;
pdTar = truncate(makedist('Normal','mu',m,'sigma',s),-Inf,arVal);
if (ismember(3,estLoc) == 1)
    if (ar < arVal)
        l3 = log(pdf(pdTar,ar));
    else
        l3 = -Inf;
    end
elseif (ismember(3,estLoc) == 0)
    l3 = 0;
end

%% Theta(4): b1
%%
%% b1 has to be between 0 and 1 because we are requiring the inflation effects to sum
%% to 1.  Thus, if it is estimated, the prior will be uniform over 0 and 1, which are
%% the standard levels of blow and bhigh that can be respecified if needed.  If it is not
%% included in the estLoc vector, it will be given the degenerate prior of 1, which in
%% log-prior space is zero.

% Need to restrict b1 to be in between blow and bhigh
blow = 0;
bhigh = 1;

if (ismember(4,estLoc) == 1)
    if (b1 <=blow | b1 >= bhigh)
        l4 = -Inf;
    else
        l4 = log(unifpdf(b1,blow,bhigh));
    end
elseif (ismember(4,estLoc) == 0)
    l4 = 0;
end


%% Theta(5): bY
%%
%% bY has the restriction that it must be positive, or in (H)LW parlance, greater than
%% 0.025.  We will keep their requirements instead of just staying positive for consistency.  If bY is
%% estimated, we will use the truncated normal distribution with mean m and standard
%% deviation s.  If it is not included in the estLoc vector, it will be given the degenerate prior of 1, which in
%% log-prior space is zero.
bYval = 0.025;
pdTbY = truncate(makedist('Normal','mu',m,'sigma',s),bYval,Inf);
if (ismember(5,estLoc) == 1)
    if (bY > bYval)
        l5 = log(pdf(pdTbY,bY));
    else
        l5 = -Inf;
    end
elseif (ismember(5,estLoc) == 0)
    l5 = 0;
end

%% Theta(6): s1
%%
%% s1 is distributed inverse gamma, so its restrictions are all inate to that
%% distribution.  If the proposed value is negative, it will be given -Inf prior
%% value.  If it is estimated, we will draw from inverse gamma with shape
%% parameter sh and scale parameter sc.
if (ismember(6,estLoc) == 1)
    if (s1 <= sig1LB)
        l6 = -Inf;
    else 
        l6 = log(unifpdf(s1,sig1LB,sigHB));
    end
elseif (ismember(6,estLoc) == 0)
    l6 = 0;
end


%% Theta(7): s2
%%
%% s2 is distributed inverse gamma, so its restrictions are all inate to that
%% distribution.  If the proposed value is negative, it will be given -Inf prior
%% value.  If it is estimated, we will draw from inverse gamma with shape
%% parameter sh and scale parameter sc.
if (ismember(7,estLoc) == 1)
    if (s2 <= sigLB)
        l7 = -Inf;
    else 
        l7 = log(unifpdf(s2,sigLB,sigHB));
    end
elseif (ismember(7,estLoc) == 0)
    l7 = 0;
end


%% Theta(8): s3
%%
%% s3 is distributed inverse gamma, so its restrictions are all inate to that
%% distribution.  If the proposed value is negative, it will be given -Inf prior
%% value.  If it is estimated, we will draw from inverse gamma with shape
%% parameter sh and scale parameter sc.
if (ismember(8,estLoc) == 1)
    if (s3 <= sigLB)
        l8 = -Inf;
    else 
        l8 = log(unifpdf(s3,sigLB,sig3HB));
    end
elseif (ismember(8,estLoc) == 0)
    l8 = 0;
end


%% Theta(9): s4
%%
%% s4 is distributed inverse gamma, so its restrictions are all inate to that
%% distribution.  If the proposed value is negative, it will be given -Inf prior
%% value.  If it is estimated, we will draw from inverse gamma with shape
%% parameter sh and scale parameter sc.
if (ismember(9,estLoc) == 1)
    if (s4 <= sigLB)
        l9 = -Inf;
    else 
        l9 = log(unifpdf(s4,sigLB,sigHB));
    end
elseif (ismember(9,estLoc) == 0)
    l9 = 0;
end


%% Theta(10): s5
%%
%% s5 is distributed inverse gamma, so its restrictions are all inate to that
%% distribution.  If the proposed value is negative, it will be given -Inf prior
%% value.  If it is estimated, we will draw from inverse gamma with shape
%% parameter sh and scale parameter sc.
if (ismember(10,estLoc) == 1)
    if (s5 <= sigLB)
        l10 = -Inf;
    else 
        l10 = log(unifpdf(s5,sigLB,sigHB));
    end
elseif (ismember(10,estLoc) == 0)
    l10 = 0;
end


%% Theta(11): aX
%%
%% aX has no restrictions so it will be evaluated, if estimated against the
%% untruncated normal distribution, with mean m and standard deviation s.  If it is not
%% included in the estLoc vector, it will be given the degenerate prior of 1, which in
%% log-prior space is zero.
if (ismember(11,estLoc) == 1)
    l11 = log(normpdf(aX,m,s));
elseif (ismember(11,estLoc) == 0)
    l11 = 0;
end


%% Theta(12): rhoG
%%
%% rhoG will be restricted to behave in a reasonable way, which for us means that it is
%% non-negative.  If estimated against the untruncated normal distribution, with mean m
%% and standard deviation s.  If it is not included in the estLoc vector, it will be
%% given the degenerate prior of 1, which in log-prior space is zero.
rhoGval = 0;
pdTrhoG = truncate(makedist('Normal','mu',m,'sigma',s),rhoGval,Inf);
if (ismember(12,estLoc) == 1)
    if (rhoG >= rhoGval)
        l12 = log(pdf(pdTrhoG,rhoG));
    else
        l12 = -Inf;
    end
elseif (ismember(12,estLoc) == 0)
    l12 = 0;
end


%% Theta(13): rhoZ
%%
%% rhoZ has no restrictions so it will be evaluated, if estimated against the
%% untruncated normal distribution, with mean m and standard deviation s.  If it is not
%% included in the estLoc vector, it will be given the degenerate prior of 1, which in
%% log-prior space is zero.
rhoZsig = 1/2;
rhoZval = 0;
pdTrhoZ = truncate(makedist('Normal','mu',1,'sigma',rhoZsig),rhoZval,Inf);
if (ismember(13,estLoc) == 1)
    if (rhoZ >= rhoZval)
        l13 = log(pdf(pdTrhoZ,rhoZ));
    else
        l13 = -Inf;
    end
elseif (ismember(13,estLoc) == 0)
    l13 = 0;
end


%% Theta(14): muG
%%
%% muG has no restrictions so it will be evaluated, if estimated against the
%% untruncated normal distribution, with mean m and standard deviation s.  If it is not
%% included in the estLoc vector, it will be given the degenerate prior of 1, which in
%% log-prior space is zero.
if (ismember(14,estLoc) == 1)
    l14 = log(normpdf(muG,m,s));
elseif (ismember(14,estLoc) == 0)
    l14 = 0;
end

%% Theta(15): muZ
%%
%% muZ has no restrictions so it will be evaluated, if estimated against the
%% untruncated normal distribution, with mean m and standard deviation s.  If it is not
%% included in the estLoc vector, it will be given the degenerate prior of 1, which in
%% log-prior space is zero.
if (ismember(15,estLoc) == 1)
    l15 = log(normpdf(muZ,m,s));
elseif (ismember(15,estLoc) == 0)
    l15 = 0;
end


%% lambda constraints 
lambG = s5/s4;
lambZ = (-ar*s3)/s1;

% Check to see if either lambda violates either bound
if ((lambG < lambLB) || (lambG > lambHB) || ...
    (lambZ < lambLB) || (lambZ > lambHB))
    % If it does violate the constraint, reject the whole draw
    lLamb = -Inf;
else
    % If it doesn't violate the constraint, go on as usual dependent on the
    % other information from theta.
    lLamb = 0;
end
    


% Add up the log priors from each element to have the overall log prior evaluation
Lpi = l1 + l2 + l3 + l4 + l5 + l6 + l7 + l8 + l9 + l10 + l11 ...
      + l12 + l13 + l14 + l15 + lLamb;


