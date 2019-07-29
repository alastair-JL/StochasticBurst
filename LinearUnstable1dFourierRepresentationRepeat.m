
reps=300;

eps= reshape(ones(reps,1)*logspace(-10,-30,9),[9*reps,1])
Ts=0*eps;

for(eee=1:length(eps) )

iii=0;

dt=0.1;

numBox=10000;

FourierModes= 80;

Fc= zeros(FourierModes,1); %Cosine values
Fs= zeros(FourierModes,1); %Sin values;

Constant=0; %Constant Concerntration;
D=1;
L=140;

addScaling=eps(eee);
Periodic=true;

noiseRate=sqrt(4*D*addScaling); 


T=0;

if(Periodic)
    SinThing= @() randn(FourierModes,1);
    Kays= (1:FourierModes)'* 2*pi/L; %List of Fourier frequencies.
else
    SinThing= @() 0;
    Kays= (1:FourierModes)'* pi/L; %List of Fourier frequencies.
end


NoiseMultiplier= noiseRate*sqrt( (1-exp(2.*(1-D.*Kays.^2).*dt))./(2.*(-1+D.*Kays.^2)) );


tnext=0.5;
Tfinal=3000;
NextRecord=0.005;

exes= linspace(0,L,numBox+1);
exes=exes(2:end);


while(T<Tfinal)
    T=T+dt;
    
xiS= SinThing().*NoiseMultiplier; %NOTE: How do I decide reasonable noise rate.
xiC= randn(FourierModes,1).*NoiseMultiplier; %NOTE: How do I decide reasonable noise rate.
%xiConstant= NoiseRate*randn(1,1)*sqrt(dt); %NOTE: How do I decide reasonable noise rate.

Fc= Fc.*exp((1-D.*Kays.^2).*dt/2)+xiC;   %Cosine values
Fs= Fs.*exp((1-D.*Kays.^2).*dt/2)+xiS;
Constant=Constant*exp((1-0).*dt/2)+randn(1,1)*noiseRate*sqrt( (1-exp(2.*(1).*dt))./(2.*(-1)) );
%Constant=Constant+xiConstant; %Note, I can set this noise to zero... if needed.

U= (Fc'*cos(Kays * exes) + Fs'*sin(Kays * exes) );

if(max(abs(U))>=1)
    break;
end

end



U= (Fc'*cos(Kays * exes) + Fs'*sin(Kays * exes) );

if(max(U)<max(-U))
    U=-U;
end

[Mx,Idx] = max(U);

fakeU= exp(- (exes-exes(Idx)).^2/(2*T) ) + exp(- (exes-exes(Idx)+L).^2/(2*T) ) +exp(- (exes-exes(Idx)-L).^2/(2*T) );

eee
T
Ts(eee)=T;

end

eps= reshape(eps,[reps,9])
Ts=reshape(Ts,[reps,9]);