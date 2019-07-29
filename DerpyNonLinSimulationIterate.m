reps=150;

eps= reshape(ones(reps,1)*logspace(-1,-5,9),[9*reps,1]);
Ts=0*eps;


dt=0.001;
L=pi*100;
K= 2*pi*[0:(length(U)/2-1),(length(U)/2):-1:1]/L;
exesExperiment=(L*(1:length(U))/length(U));

for(eee=1:length(eps))

U =zeros(1,2^11);

noise=eps(eee);

T=0;
nextT=1;
tic();

while(max(U)<1)
Uf=fft(U);
Uf=Uf.*(exp(-dt*K.^2));
U=ifft(Uf,'symmetric');
U=U+U.^2*dt;
U=U+ randn(size(U))*sqrt(dt*exesExperiment(1))*noise;
T=T+dt;
end

timeDone= toc();
TimeGo= timeDone*(length(eps)-eee)/eee;

Ts(eee)=T;
[eee,T,TimeGo/3600,mod(TimeGo/60,1)*60]
end



eps= reshape(eps,[reps,9]);
Ts=reshape(Ts,[reps,9]);

boxplot(log(Ts));
