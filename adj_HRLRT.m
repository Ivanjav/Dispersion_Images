function [m,J]=adj_HRLRT(d,f,x,p,S,num_iters)
%Adjoint HRLRT
mi=adj_LRT(d,f,x,p);
J=zeros(num_iters,1);
%J=zeros(length(f),num_iters);
m=mi;

wait = waitbar(0,'Please wait...');
pause(.5)
for k=1:num_iters
waitbar(k/num_iters,wait,'Computing the HRLRT');
Jf=zeros(length(f),1);
for i=1:length(f)
   Wm=diag(1./sqrt(abs(m(:,i))));
   %Wm=diag(sqrt(abs(m(:,i))));
   L=exp(1j*2*pi*f(i)*x'*p);
   r=L*inv(Wm)*Wm*m(:,i)-d(:,i);
   %Wd=diag(1./sqrt(abs(std(r))));
   Wd=diag(1./sqrt(std(r)*ones(length(x),1)));
   
   B=inv(Wm)'*L'*L;
   lambda=S*(sum(B(1,:))^2)/sum((B(1,:)).^2);
   Jf(i)=norm(Wd*(d(:,i)-L*inv(Wm)*Wm*m(:,i)))+lambda*norm(Wm*m(:,i));
   A=(inv(Wm)'*L'*(Wd'*Wd)*L*inv(Wm)+lambda*eye(length(p)))*Wm;
   b=inv(Wm)'*L'*(Wd'*Wd)*d(:,i);
   %m(:,i) = cgs(A,b);
   m(:,i) = cgs(A,b,1e-10,1e3);
   %m(:,i) = pcg(A,b,1e-10,1e3);
end
J(k)=sum(Jf);
end
delete(wait)
end