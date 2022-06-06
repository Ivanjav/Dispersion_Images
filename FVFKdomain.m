function [FKdata,FVdata,frequency,wavenumber,velocity]=FVFKdomain(Data,dx,dt,N,fmin,fmax,vmin,vmax,vstep)

FKdata=fft2(Data,N,N);
wavenumber=linspace(0,1/dx,N);
k_step=wavenumber(2);
FKdata=FKdata(:,end:-1:1);
%FKdata=FKdata(1:round(N/2),round(N/2):end);
fs=1/dt;
fstep=fs/N;
frequency=0:fstep:fmax;
%frequency = fs*(0:1);
velocity=vmin:vstep:vmax;
FVdata=zeros(length(velocity),length(frequency));

%v_step=f_step/k_step;
for i=1:length(frequency)
    indices=(frequency(i)./velocity)/k_step+1;
    %ifreq=floor(frequency(i)*(N/fs))+1
    ifreq=i;
    FVdata(:,i)=(ceil(indices)-indices).*FKdata(ifreq,floor(indices))+...
                (indices-floor(indices)).*FKdata(ifreq,ceil(indices));
    idk=find(abs(indices-round(indices))==0);
    FVdata(idk,i)=FKdata(ifreq,indices(idk));
end

wavenumber=linspace(0,1/(2*dx),(size(FKdata,2)-1)/2+1);
frequency=fmin:fstep:fmax;
FKdata=FKdata(round(frequency*N/fs),1:length(wavenumber));
