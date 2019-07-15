NT=[40 80 160 320];
T=30;
% global dx
for i=1:4
    N=NT(i);
    dx=2/N;
    cfl=0.8;
    dt=cfl*dx;
    x=zeros(1,N);
    u1=x;
    u2=x;
    u1_acu=u1;
    u2_acu=u2;
    u1lft=u1;
    u2lft=u2;
    D=zeros(N,N);
    for j=1:N
        x(j)=-1+(j-1)*dx;
        u1(j)=exp(-x(j)^2/.04);
        if x(j)>=-0.4 && x(j)<=0.4
            u2(j)=1;
        else
            u2(j)=0;
        end
    end
    Nt=round(T/dt);
    Nx=N;
    xgrid=(-1):dx:(-1+(Nx-1)*dx);
    Nxs=Nx+1;
    u1up=u1;
    u2up=u2;
    u1lf0=u1;
    u2lf0=u2;
    for kk1=1:N
        for jj1=1:N
            if kk1==jj1
               D(kk1,jj1)=0;
            else
                D(kk1,jj1)=pi/2*(-1)^(kk1-jj1)/sin(pi/N*(kk1-jj1));
            end
        end
    end   
    
    for kxa=1:Nx
        xt=dx*(kxa-1)-1-T;
%         u1_acu(kxa)=exp(-xt^2/0.04);
        xt_j=abs(xt)+1-2*floor((abs(xt)+1)/2)-1;
        u1_acu(kxa)=exp(-xt_j^2/.04);
        if xt_j>=-0.4 && xt_j<=0.4
            u2_acu(kxa)=1;
        else
            u2_acu(kxa)=0;
        end
    end
    
    for k21=1:Nx
        u1lf0tem=[u1lf0(end),u1lf0,u1lf0(1)];
        u2lf0tem=[u2lf0(end),u2lf0,u2lf0(1)];
        u1lf0(k21)=u1lf0tem(k21+1)-dt/dx/2*(u1lf0tem(k21+2)-u1lf0tem(k21));
        u2lf0(k21)=u2lf0tem(k21+1)-dt/dx/2*(u2lf0tem(k21+2)-u2lf0tem(k21));
    end
    u1lftem1=u1;
    u2lftem1=u2;
    u1lftem2=u1lf0;
    u2lftem2=u2lf0;
    u1lf=u1;
    u2lf=u2;
    u1lw=u1;
    u2lw=u2;
    u1temcn=u1';
    u2temcn=u2';
    u1gl=u1';
    u2gl=u2';
    for kt=1:Nt
        u1temup=[u1up(end),u1up];
        u2temup=[u2up(end),u2up];
        u1lftem=[u1lftem2(end),u1lftem2,u1lftem2(1)];
        u2lftem=[u2lftem2(end),u2lftem2,u2lftem2(1)];
        u1temlf=[u1lf(end),u1lf,u1lf(1)];
        u2temlf=[u2lf(end),u2lf,u2lf(1)];
        u1temlw=[u1lw(end),u1lw,u1lw(1)];
        u2temlw=[u2lw(end),u2lw,u2lw(1)];
        %Crank Nicolson
        cn=dt/4/dx;
        Mcn=diag(ones(1,Nx))+diag(cn*ones(1,(Nx-1)),1)-diag(cn*ones(1,(Nx-1)),-1);
        Mcn(1,Nx)=-cn;
        Mcn(Nx,1)=cn;
        u1cn=Mcn\(Mcn'*u1temcn);
        u2cn=Mcn\(Mcn'*u2temcn);
        u1temcn=u1cn;
        u2temcn=u2cn;
        u1cn0=u1cn';
        u2cn0=u2cn';
        %Global
        u1gltem=u1gl;
        u2gltem=u2gl;
        du1gl=D*u1gltem;
        du2gl=D*u2gltem;
        u1gl=u1gltem-dt*du1gl;
        u2gl=u2gltem-dt*du2gl;
    for kx=1:Nx
        % Upwind
        u1up(kx)=u1temup(kx+1)-dt/dx*(u1temup(kx+1)-u1temup(kx));
        u2up(kx)=u2temup(kx+1)-dt/dx*(u2temup(kx+1)-u2temup(kx));
        %leap frog
        u1lft(kx)=u1lftem1(kx)-dt/dx*(u1lftem(kx+2)-u1lftem(kx));       
        u2lft(kx)=u2lftem1(kx)-dt/dx*(u2lftem(kx+2)-u2lftem(kx));    
        % lax friedrich
        u1lf(kx)=.5*(u1temlf(kx)+u1temlf(kx+2))-dt/dx/2*(-u1temlf(kx)+u1temlf(kx+2));
        u2lf(kx)=.5*(u2temlf(kx)+u2temlf(kx+2))-dt/dx/2*(-u2temlf(kx)+u2temlf(kx+2));
        %lax wendroff
        u1lw(kx)=u1temlw(kx+1)+dt^2/dx^2/2*(u1temlw(kx)+u1temlw(kx+2)-2*u1temlw(kx+1))-dt/dx/2*(-u1temlw(kx)+u1temlw(kx+2));
        u2lw(kx)=u2temlw(kx+1)+dt^2/dx^2/2*(u2temlw(kx)+u2temlw(kx+2)-2*u2temlw(kx+1))-dt/dx/2*(-u2temlw(kx)+u2temlw(kx+2));          
    end
    u1lftem1=u1lftem2;
    u2lftem1=u2lftem2;
    u1lftem2=u1lft;
    u2lftem2=u2lft;
    
    end
    u1lfr=u1lftem1;
    u2lfr=u2lftem1;
    
    figure(10*i)
    plot(xgrid,u1_acu,'k')
    hold on
    plot(xgrid,u1up,'b')
    plot(xgrid,u1lfr,'g')
    plot(xgrid,u1lf,'r')
    plot(xgrid,u1lw,'c')
    plot(xgrid,u1cn0,'m')
%     plot(xgrid,u1gl,'y')
    hold off
    legend('accurate','upwind','leap frog','lax friedrichs','lax wendroff','crank nicolson','global')
    ss=num2str(N);
    title(ss);
    xlabel('x')
    ylabel('u1')
%     figure(100*i)
%    plot(xgrid,u1_acu,'k')
%     figure(100*i+1)
%      plot(xgrid,u2_acu)
    
    figure(10*i+1)
    plot(xgrid,u2_acu,'k')
    hold on
    plot(xgrid,u2up,'b')
    plot(xgrid,u2lfr,'g')
    plot(xgrid,u2lf,'r')
    plot(xgrid,u2lw,'c')
    plot(xgrid,u2cn0,'m')
%     plot(xgrid,u2gl,'y')
    hold off
    legend('accurate','upwind','leap frog','lax friedrichs','lax wendroff','crank nicolson','global')
    title(ss);
    xlabel('x')
    ylabel('u2')
end


