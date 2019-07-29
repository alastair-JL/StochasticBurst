
D=1;
uTarget=1;
m=1;
Records=[];

if(~exist('tf','var'))
tf=15
end


UInitialConditions= @(x) 0;

uFunct= @(x,t) 0; %griddedInterpolant({xs,tsQ'},MQ)


farfield= @(t) 0;
Gamma=0;

    xbase = linspace(0,sqrt(4*D*tf)*10,800);
    tbase = linspace(0,tf,350);

    
dx= xbase(2);

StumpyBoundaryXi= @(xl,ul,xr,ur,t) deal(0,1,ur,0);

StumpyBoundaryU= @(xl,ul,xr,ur,t) deal(0,1,farfield(t)-ur,0);


ev= @(m,t,xmesh,umesh) myEventFunction(umesh,uTarget);

    %This part of the code is pretty hacky.
%MQ=[zeros(size(M,1),1),M];
%tsQ=[-1,ts]+1;
    %End hacky.
    
oldTotalCosts=50000;

if(~exist('ElFar','var'))
    ElFar=1/(tf+1/uTarget);
else
    
end
ElFarStart=ElFar;

alphaIn=[];
ResultOut=[];


jumpFromLow=1;

Totalsteps=0;
        integrationAssistVector = 0*xbase; %This is a vector used to weight terms when in order to approximate our integral in polar coordinates.

if(m==1)
deltaR= xbase(2);

    XiInitialConditions= @(x) (x<deltaR).*(deltaR-x).*3/(pi*deltaR^2) ; %Initial condition is a cone of volume 1
        integrationAssistVector = 0*xbase; %This is a vector used to weight terms when in order to approximate our integral in polar coordinates.


Ars=xbase(2:end);

integrationAssistVector(1)= pi*deltaR.^2/3 %%Area of cone. "Height" given by what we multiply by.

integrationAssistVector(2:end)= (pi./deltaR).*...
   ( ( (Ars+deltaR).^3/6)-( (Ars+deltaR).*Ars.*Ars/2-Ars.^3/3) ...
        - ( (Ars-deltaR).^3/6)+( (Ars-deltaR).*Ars.*Ars/2-Ars.^3/3) );


elseif(m==0)
    XiInitialConditions= @(x) (x<dx)./dx; %How the hell do I pick a good ElFar....
integrationAssistVector=integrationAssistVector+2;
integrationAssistVector(1)=1;
integrationAssistVector=integrationAssistVector*xbase(2);
else
    error('bad m value');
end


if(~exist('ElFar','var'))
    ElFar=1/(tf+1/uTarget);
else
    
end
if(isnan(ElFar))
    ElFar=1/(tf+1/uTarget);
end
ElFarStart=ElFar;


for(ppq=1:6)

ExpectedAccuracyAlpha=10^-4;
    backPDE= @(x,t,xi,DuDx) BackwardGoPde(x,t,xi,DuDx,uFunct,D,tf);
    solback = pdepe(m,backPDE,XiInitialConditions,StumpyBoundaryXi,xbase,tbase);
    
    XI = solback(:,:,1);
    
    
    jumpFromLow=0.25*jumpFromLow;
    
       RawIntegral= sum( integrationAssistVector*(XI').^2 )*(tbase(2)-tbase(1));
       ElFar=uTarget./RawIntegral;

       OverData=[nan,nan,nan];
       UnderData=[nan,nan,nan];
       NewData=[nan,nan,nan];
       
       ElFarNow=0;
       ElFarOld=0;
       
       ElFarHigh=inf;
       ElFarLow=0;
       
       wupNow=0;
       wupOld=0;
       
       wupTarget= tf+1/uTarget;
       
       te=0;
       ue=0;

steps=0;

    while ( ((abs(te-tf)/abs(tf))>ExpectedAccuracyAlpha) || ((abs(ue-uTarget)/abs(uTarget))>ExpectedAccuracyAlpha) )
    Totalsteps=Totalsteps+1;
    steps=steps+1;
    XiFunct= griddedInterpolant({xbase,tbase},ElFar*XI');
        
    forwardPDE= @(x,t,u,DuDx) ForwardGoPde(x,t,u,DuDx,XiFunct,D,tf,Gamma);
    [solForward,tsol,sole,te,ie] = pdepe(m,forwardPDE,UInitialConditions,StumpyBoundaryU,xbase,tbase,odeset('Events',ev));
    
    U = solForward(:,:,1);
    
    
    if(isempty(te))
        te=tf;
        ue=(solForward(end,1,1));
        NewData=[ElFar,tf,ue]
        nextElFar= ElFar+(NewData(1)-UnderData(1))*(uTarget-NewData(3))/(NewData(3)-UnderData(3));
        %    nextElFar=sqrt(OverData(1)*NewData(1));
        tooHigh=false;
            if(nextElFar>OverData(1))
                nextElFar=sqrt(OverData(1)*NewData(1));
            end

            if(isnan(nextElFar))
                nextElFar=ElFar*2;
            end
            
            
        UnderData=NewData;
    else
        te=te(1)
        ue=uTarget;
        NewData=[ElFar,te,uTarget]
        %     nextElFar=sqrt(NewData(1)*UnderData(1));
        nextElFar= ElFar+(NewData(1)-OverData(1))*(tf-NewData(2))/(NewData(2)-OverData(2));
        tooHigh=true;
            if(nextElFar<UnderData(1)||rand()<0.2)
                nextElFar=sqrt(NewData(1)*UnderData(1));
            end
    
            if(isnan(nextElFar))
                nextElFar=ElFar*0.5;
            end
            
        OverData=NewData;
    end
     ElFar=nextElFar;
     if(ElFar<0)
        ElFar=-0.1*ElFar;
     end
     
     
        if((OverData(1)-UnderData(1))<ElFar*10^-9)
            warning('The window is crazy small. Lets call this close enough.');
            break;
        end
       
        if(steps>1500)
           error('Things are taking to long. Bailing out.');
        end
        
        NewData
        ElFar
        
    end
    
end

widthSpike=xbase(sum(U(end,:)>(U(end,1)/2)))

if(exist('longRecord','var'))
    longRecord=[ longRecord, [ElFar;m;te;ue;ElFar*ElFar*RawIntegral;widthSpike;Gamma]];   
else
    longRecord=[ElFar;m;te;ue;ElFar*ElFar*RawIntegral;widthSpike;Gamma];
end


xOptim= [-fliplr(xbase),xbase];
uOptim= [fliplr(U(end,:)),U(end,:)]

plot(xOptim+xM,uOptim);

beep;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

function [val,term,updown] = myEventFunction(umesh,uTarget)

val=max(max(max(umesh))-uTarget);
term=1;
updown=1;
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,f,s] = BackwardGoPde(x,t,xi,DuDx,u,D,tf)
c = 1;
f = D*DuDx;
s = xi+2.*u(x,tf-t).*xi;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [c,f,s] = ForwardGoPde(x,t,u,DuDx,xi,D,tf,Gamma)
c = 1;
f = D*DuDx;
s = u.^2 + xi(x,tf-t)+Gamma;
end