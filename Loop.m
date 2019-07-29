

D=1;
uTarget=1;
m=1;
Records=[];

tfs= [5,15,50];
Gammas=[0,-10^-4,-10^-3,-10^-2,-10^-1,-1];

for ggg=1:length(Gammas)
    for tttrt= 1:length(tfs) 
        Gamma=Gammas(ggg);
        tf=tfs(tttrt);
BackForwardSolveQuadWithGamma;

    end
end