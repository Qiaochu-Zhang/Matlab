clear
format long
NT=[10 20 40 60 80];
er1=zeros(1,5);
er2=er1;
er4=er1;
erg=er1;
for i=1:5
    N=NT(i);
    Nt=N+1;
    Ns=Nt+1;
    h=2*pi/Nt;
    x=zeros(1,Nt);
    u=x;
    du=x;
    du1=x;
    du2=x;
    du4=x;
    D=zeros(Nt,Nt);
for k1=1:Nt
    k=k1-1;
    for j1=1:Nt
        j=j1-1;
        if k1==j1
            D(k1,j1)=0;
        else
            D(k1,j1)=((-1)^(k1-j1))/2/sin((pi/Nt)*(k-j));
        end
    end
end
    for j=1:Nt
        x(j)=(j-1)*h;
        u(j)=exp(sin(x(j)));
        du(j)=u(j)*cos(x(j));
    end
    ut=[u(N),u(Nt),u,u(1),u(2)];
    
    du_g=(D*u')';
    
    for k=1:Nt
        du1(k)=(ut(k+3)-ut(k+2))/h;
        du2(k)=(ut(k+3)-ut(k+1))/2/h;
        du4(k)=(8*(ut(k+3)-ut(k+1))-(ut(k+4)-ut(k)))/12/h;
    end
    figure(i)
    plot(x,du,'k')
    hold on
    plot(x,du1,'r')
    plot(x,du2,'g')
    plot(x,du4,'b')
    hold off
    er1(i)=norm(du1-du)/sqrt(N);
    er2(i)=norm(du2-du)/sqrt(N);
    er4(i)=norm(du4-du)/sqrt(N);
    erg(i)=norm(du_g-du)/sqrt(N);
    figure(10+i)
    plot(x,du,'k')
    hold on
    plot(x,du_g,'r')
    hold off
    
end
figure(7)
plot(NT,er1,NT,er2,NT,er4)
legend('1st-order error','2nd order error','3rd order error')
    
figure(8)
loglog(NT,er1,NT,er2,NT,er4)  
legend('1st-order error','2nd order error','3rd order error')
title('Error vs N, Lagrange approxiamation')
xlabel('N')
ylabel('Error')
figure(100)
semilogy(NT,erg)
title('Error vs N Global approximation, semilogy')
xlabel('N')
ylabel('Error')
figure(101)
loglog(NT,erg)
title('Error vs N Global approximation, loglog')
xlabel('N')
ylabel('Error')
