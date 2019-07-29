U =zeros(1,2^11);

dt=0.001;
L=pi*100;
K= 2*pi*[0:(length(U)/2-1),(length(U)/2):-1:1]/L
noise=10^-1.1;
exesExperiment=(L*(1:length(U))/length(U));


T=0;
nextT=1;

while(max(U)<1)
Uf=fft(U);
Uf=Uf.*(exp(-dt*K.^2));
U=ifft(Uf,'symmetric');
U=U+U.^2*dt;
U=U+ randn(size(U))*sqrt(dt*exesExperiment(1))*noise;
T=T+dt;
end


Uexperiment=U;
[uM,xM]=max(U);
xM=exesExperiment(xM);
plot([exesExperiment,exesExperiment+L],[Uexperiment,Uexperiment]);

Gamma=0;
tf=T;