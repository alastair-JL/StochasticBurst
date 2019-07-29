%f = .088

% Size of gridcol
width = 256;
% Diffusion rates
da = (10^-4)*width*width;
noise=10^-1;

% 5,000 simulation seconds with 4 steps per simulated second
dt = .05;
stoptime = 120;

t=0;
A=Atarget*ones(width,width);
B = 1-mean(mean(A));


h=figure();

frametime = 0.005;
nextframe = 0;
tic

rangebounder= zeros(width,1);
rangebounder(1,1)=6;

axes('Position',[0 0 1 1])
axis off
hi = image(A);
hi.CDataMapping = 'scaled';

ht = text(3,width-3,'Time = 0');
ht.Color = [.95 .2 .8];

%gamma=2;
%k0=-0.1;
colorbar
trigger=0%.001;

nframes = 1;

filename = 'LeahModelUnder.gif';
nframes = 1;
        hi.CData = [A,rangebounder];
        
%        caxis([0.05,0.15]);
        
        drawnow
       
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        
imwrite(imind,cm,filename,'gif', 'Loopcount',0,'DelayTime',0); 


timeCourseData=zeros(3,ceil(stoptime/frametime));
xyz=0;
Afull=1.00

while t<stoptime
    A = A + (da*my_laplacian(A) + (B.*(k0 + gamma*A.*A./(1+A.*A) ) -A) )*dt+ randn(size(A))*noise*sqrt(dt)/width;
    B = Afull-mean(mean(Atarget));
    k0=k0+trigger;
    %hi.CData = A;   
    t = t+dt;
    ht.String = ['Time = ' num2str(floor(t)) '  k0 = ' num2str(k0)];

    if t > nextframe
        xyz=xyz+1;
        timeCourseData(1,xyz)=t;
        timeCourseData(2,xyz)=min(min(A));
        timeCourseData(3,xyz)=max(max(A));
%     bob=[(B-mean(mean(B)))/sqrt(mean(var(B)))];%[ (A-mean(mean(A)))/sqrt(mean(var(A))) ,(B-mean(mean(B)))/sqrt(mean(var(B)))];%,rangebounder];
        hi.CData = [A,rangebounder];
        drawnow
        nextframe = nextframe + frametime;

        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        imwrite(imind,cm,filename,'gif','WriteMode','append');         
        
        if(isnan(B(1,1))|| min(min(B))<0 )
            break
        end
    
        nframes = nframes+1;    
        
    end

end

delta = toc;
disp([num2str(nframes) ' frames in ' num2str(delta) ' seconds']);

figure()
plot(timeCourseData(1,1:xyz),timeCourseData(2,1:xyz))
hold on
plot(timeCourseData(1,1:xyz),timeCourseData(3,1:xyz))

xlabel('time')
ylabel('u')
axis([0,120,0,6])
