NT=[50 100 200 400];
global T
T=3;
global dx
global N
global m
app=0;
app1=0;
kkp=100*ones(4,(N+1));
error=ones(1,4);
for i=1:4
% i=4;
    N=NT(i);
    xp=2*pi;
    dx=xp/N;
    dt=0.01;
    cfl=1.5*dt/dx;
%     dt=.8*dx;
    Ntime=round(T/dt);
% cfl=dt/dx
x=ones(1,N);
u_a=x;
u_p=ones(1,N+1);
    for j=1:N
        x(j)=0+(j-1)*dx;
        u_a(j)=(fhw21cos(x(j)+dx)-fhw21cos(x(j)))/dx;
    end
    
    for it=1:Ntime

        uatem=u_a;
%         u_a=fivohw21(uatem)*dt+uatem;

        k1=fivohw21(uatem);
        k2=fivohw21(uatem+.5*k1*dt);
        k3=fivohw21(uatem-k1*dt+k2*dt*2);
        u_a=uatem+1/6*dt*(k1+4*k2+k3);
%     figure(100)
%     plot(u_a)
%     drawnow
    end
%     figure(88)
%     plot(u_a)
    
    u_ab=[u_a(N-1),u_a(N),u_a,u_a(1)];
    x_p=[x,(x(N)+dx)];
    x_ana1=zeros(1,(N+1));
    x_ana2=x_ana1;
    x_ana=x_ana1;
    u_ana=x_ana;
    for m=1:(N+1)
        u_p(m)=-1/6*u_ab(m)+5/6*u_ab(m+1)+1/3*u_ab(m+2);
        ppana1=dx*(m-1)-3/2*T;
        x_ana1(m)=fsolve(@fhw21sin,ppana1);
%         n1=0;
%         n2=0;
%         while fhw21sin(x_ana1(m))>0.001
%             n1=n1+1;
% %             tt11=x_ana1(m);
%             ppana1=ppana1+0.1;
%             x_ana1(m)=fsolve(@fhw21sin,ppana1);
%             if n1>31
%                 app1=500;
%                 break
%             end
%         end
        ppana2=dx*(m-1)-1/2*T;
        x_ana2(m)=fsolve(@fhw21sin,ppana2);
%         while fhw21sin(x_ana2(m))>0.001
% %             tt11=x_ana2(m);
%             n2=n2+1;
%             ppana2=ppana2-0.1;
%             x_ana1(m)=fsolve(@fhw21sin,ppana2);
%             if n2>30
%                 app1=501;
%                 break
%             end
%         end 
        f_ana1=fhw21sin(x_ana1(m));
        f_ana2=fhw21sin(x_ana2(m));
        xxnew=(m-1)*dx;
        if abs(x_ana1(m)-x_ana2(m))<10^(-6)
            x_ana(m)=x_ana1(m);
            kkp(i,m)=0;
        elseif abs(f_ana1)<abs(f_ana2) && abs(f_ana1-f_ana2)>10^(-6)
            x_ana(m)=x_ana1(m);
            kkp(i,m)=1;
        elseif abs(f_ana1)>abs(f_ana2) && abs(f_ana1-f_ana2)>10^(-6)
            x_ana(m)=x_ana2(m);
            kkp(i,m)=2;
        elseif abs(f_ana1)<10^(-6) && abs(f_ana2)<10^(-6)
            xjg=(xxnew-T);
            if  xjg<pi && xjg>0
                x_ana(m)=x_ana1(m);
                kkp(i,m)=3;
            else
                x_ana(m)=x_ana2(m);
                kkp(i,m)=4;
            end
        else
            app=99;
            kkp(i,m)=5;
            break
        end
        u_ana(m)=1+.5*sin(x_ana(m));
    end
%     if app==99
%                   break
%     end
    figure(i)
    plot(x_p,u_p)
    hold on
    plot(x_p,u_ana)
    hold off
    legend('Finite volume','Analytical')
    xlabel('x')
    ylabel('u')
    title(NT(i))
%     figure(i*2)
%     plot(x_p,u_p)
%     drawnow
%         xlabel('x')
%     ylabel('u')
%     figure(i*3)
%     plot(x_p,u_ana)
%     drawnow
%         xlabel('x')
%     ylabel('u')

error(i)=norm(u_ana-u_p)/((N+1)^0.5);

end
  
