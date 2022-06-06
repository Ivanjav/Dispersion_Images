function [f,v,p,m,FVdata]=fv_domain_HRLRT(Data,dt,x,v,fmin,fmax,S,num_iters)

fs=1/dt;
Nt=size(Data,1);
Nf=4001;
%Nf=Nt;
FXdata= fft(Data,Nf,1);
f=linspace(0,fs/2,round(Nf/2));
idf=find((f>=fmin).*(f<=fmax));
f=f(idf);
d=FXdata(idf,:)';
%d=FXdata';
pstep=5e-5;
p=1/(v(end)+20):pstep:1/v(1);
%p=1/(v(end)):pstep:1/v(1);
m=adj_HRLRT(d,f,x,p,S,num_iters);

vp=1./p;

FVdata=zeros(length(v),length(f));
for i=1:length(f)
  FVdata(:,i) = interp1(vp,m(:,i),v);
end

%idf=find((f>=fmin).*(f<=fmax));
%f=f(idf);
%FVdata=FVdata(:,idf);
end