function [f,v,FVdata]=fv_domain_LRT(Data,dt,x,v,fmin,fmax)

fs=1/dt;
vmin=v(1);
vmax=v(end);
Nt=size(Data,1);
Nf=4001;
FXdata= fft(Data,Nf,1);
d=FXdata(1:round(Nf/2),:)';
f=linspace(0,fs/2,round(Nf/2));
idf=find((f>=fmin).*(f<=fmax));
f=f(idf);
d=d(:,idf);
%p=linspace(1/vmax,1/vmin,500);
pstep=5e-5;
p=1/v(end):pstep:1/v(1);


mp=adj_LRT(d,f,x,p);


vp=1./p;
%v=linspace(vmin,vmax,length(p));
m=zeros(length(v),length(f));
for i=1:length(f)
  m(:,i) = interp1(vp,mp(:,i),v);
end

idf=find((f>=fmin).*(f<=fmax));
f=f(idf);
FVdata=m(:,idf);
end