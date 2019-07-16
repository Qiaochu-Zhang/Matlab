function f=fivohw21(u)
global dx
global N
utem=[u(N-1),u(N),u,u(1)];
Nt=N+1;
u05=zeros(1,Nt);
f05=u05;
for ii=1:Nt
    u05(ii)=-1/6*utem(ii)+5/6*utem(ii+1)+1/3*utem(ii+2);
    f05(ii)=.5*u05(ii)^2;
end
f=zeros(1,N);
for jii=1:N
    f(jii)=-(f05(jii+1)-f05(jii))/dx;
end

