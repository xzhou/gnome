function [Tc Tc_log Tr TD] = Tc(sample_S, sample_R, A)
% the function compute the test statistic based on pairwise counts
% Sampel_R: reference sample
% sample_S: sample of interest
% A: individual haplotype to be tested
nR = size(sample_R,1);
nS = size(sample_S,1);
Len = size(sample_R,2);
pseudocount = 0.5;

C11R = sample_R'*sample_R;
C_1R =repmat(sum(sample_R,1),Len,1);
C01R = C_1R-C11R;
% C01= (1-sample_R')*sample_R;
C10R = C_1R' - C11R;
% C01= sample_R'* (1-sample_R);
C00R = nR*ones(Len,Len)-C11R-C01R-C10R;

C11S = sample_S'*sample_S;
C_1S =repmat(sum(sample_S,1),Len,1);
C01S = C_1S-C11S;
% C01S= (1-sample_S')*sample_S;
C10S = C_1S' - C11S;
% C01S= sample_S'* (1-sample_S);
C00S = nS*ones(Len,Len)-C11S-C01S-C10S;

A11 = A'*A; 
A01 = (1-A')*A;
A10 = A'* (1-A);
A00 = (1-A') * (1-A);
A2 = (2*A'-1)*(2*A-1);
Tc = (C11S/nS-C11R/nR).*(2*A11-1) + (C10S/nS-C10R/nR).*(2*A10-1) + ...
    (C01S/nS-C01R/nR).*(2*A01-1) + (C00S/nS-C00R/nR).*(2*A00-1);
all_r_S = (C11S.*C00S-C10S.*C01S)./sqrt((C11S+C10S).*(C11S+C01S).*(C00S+C10S).*(C00S+C01S));
all_r_R = (C11R.*C00R-C10R.*C01R)./sqrt((C11R+C10R).*(C11R+C01R).*(C00R+C10R).*(C00R+C01R));
D_S = (C11S.*C00S-C10S.*C01S)/nS^2;
D_R = (C11R.*C00R-C10R.*C01R)/nR^2;


C00S = C00S + pseudocount;
C01S = C01S + pseudocount;
C10S = C10S + pseudocount;
C11S = C11S + pseudocount;
C00R = C00R + pseudocount;
C01R = C01R + pseudocount;
C10R = C10R + pseudocount;
C11R = C11R + pseudocount;
Tc_log = log(C11S./C11R).*A11 + log(C10S./C10R).*A10 + ...
    log(C01S./C01R).*A01 + log(C00S./C00R).*A00 ...
    -log(nS/nR);
Tc = sum(sum(Tc));
Tc_log(~isfinite(Tc_log))=0;
Tc_log = sum(sum(Tc_log));
all_r_S(~isfinite(all_r_S))=0;
all_r_R(~isfinite(all_r_R))=0;
Tr = sum(sum((all_r_S-all_r_R).* A2))/2;
TD = sum(sum((D_S-D_R).* A2))/2;

return