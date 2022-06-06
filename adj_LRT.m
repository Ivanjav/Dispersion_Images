function m=adj_LRT(d,f,x,p)
%Adjoint Linear Radon Transform
m=zeros(length(p),length(f));
for i=1:length(f)
   L=exp(1j*2*pi*f(i)*x'*p);
   m(:,i)=L'*d(:,i);
   %m(:,i)=inv(L'*L+100*eye(length(p)))*L'*d(:,i);
end

end