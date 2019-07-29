t=linspace(0,110,1000);
x=linspace(-80,80,1001);

delT= t(2)-t(1);


baseU= 0*x;
NormThing = @(x,t) ((4*pi*t)^-0.5)*exp(-x.^2/(4*t));

for(ttt=1:(length(t)-1))
   tnow=t(ttt);
   baseU= baseU + delT*exp(2*(-t(ttt)))*NormThing(x,2*(t(end)-t(ttt)));
   
end

baseU=baseU./baseU(501);

plot(x,baseU)
