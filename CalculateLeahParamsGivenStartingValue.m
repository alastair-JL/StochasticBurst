Atarget=0.1;

gamma= (1+Atarget*Atarget).^2./(2*Atarget*(1-Atarget)^2);
k0= Atarget/(1-Atarget) - gamma*Atarget*Atarget/(1+Atarget*Atarget)

A=linspace(0,2);
y= (1-A).*(k0+(gamma*A.^2)./(1+A.^2))
plot(A,y);
hold on
plot(A,A);

Afinal= A(sum(y>=A))
k0
