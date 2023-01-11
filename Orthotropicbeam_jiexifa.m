clear;clc;
e1=[16400 1300 16400];
e2=[1300 900 1300];
g=[1180 79 1180];
u=[0.43 0.63 0.43];
h=[20 40 20];l=1000;
 dh=0.3;
p=length(h);
Kt=10;
p=length(h);

for i=1:p
    d1(i)=0;
    for j=1:i
        d1(i)=d1(i)+h(j);
    end
    d0(i)=d1(i)-h(i);
end
for i=1:p
    ya{i}=d0(i):(h(i)/20):d1(i);
    aU{i}=0;
    aV{i}=0;
    aSx{i}=0;
    aSy{i}=0;
    aTxy{i}=0;
end
syms y x
N=13;
for mm=1:2:N
    am=mm*pi/l;
    Y10=0;q1=1;Yph=2*q1/mm/pi*(cos(mm*pi)-1);
   % Y10=0;Yph=-1;
    
    for i=1:p
        %U V Txy Sy
        s11(i)=1/e1(i);
        s12(i)=-u(i)/e1(i);
        s22(i)=1/e2(i);
        s66(i)=1/g(i);

        f1(i)=2*(am^2)*s12(i)+(am^2)*s66(i);
        f2(i)=sqrt(4*(am^4)*(s12(i)^2)+4*(am^4)*s12(i)*s66(i)+(am^4)*(s66(i)^2)-4*(am^4)*s11(i)*s22(i));
        g1(i)=sqrt((f1(i)-f2(i))/2/s11(i));g2(i)=sqrt((f1(i)+f2(i))/2/s11(i));
        g1(i)=double(g1(i));g2(i)=double(g2(i));
        Sz{i}=[(-(g1(i)^2)*s11(i)/am+am*s12(i))*exp(-g1(i)*y),(-(g1(i)^2)*s11(i)/am+am*s12(i))*exp(g1(i)*y),(-(g2(i)^2)*s11(i)/am+am*s12(i))*exp(-g2(i)*y),(-(g2(i)^2)*s11(i)/am+am*s12(i))*exp(g2(i)*y);
            (-g1(i)*s12(i)+s22(i)*(am^2)/g1(i))*exp(-g1(i)*y),(g1(i)*s12(i)-s22(i)*(am^2)/g1(i))*exp(g1(i)*y),(-g2(i)*s12(i)+s22(i)*(am^2)/g2(i))*exp(-g2(i)*y),(g2(i)*s12(i)-s22(i)*(am^2)/g2(i))*exp(g2(i)*y);
            am*g1(i)*exp(-g1(i)*y),-am*g1(i)*exp(g1(i)*y),am*g2(i)*exp(-g2(i)*y),-am*g2(i)*exp(g2(i)*y);
            -am^2*exp(-g1(i)*y),-am^2*exp(g1(i)*y),-am^2*exp(-g2(i)*y),-am^2*exp(g2(i)*y)];
        Dh{i}=subs(Sz{i},y,d1(i));
        D0{i}=subs(Sz{i},y,d0(i));
         Dh{i}=eval(Dh{i});
         D0{i}=eval(D0{i});
    end
    M=eye(4);
   K=[1,0,1/Kt,0;0,1,0,0;0,0,1,0;0,0,0,1];
    for i=p:-1:2
        M=M*Dh{i}*inv(D0{i})*K;
    end
    M=M*Dh{1}*inv(D0{1});
    M11=M(1:2,1:2);
    M12=M(1:2,3:4);
    M21=M(3:4,1:2);
    M22=M(3:4,3:4);
    UV10=inv(M21)*([0;Yph]-M22*[0;Y10]);
    for i=1:p
        UVZYh{i}=eye(4);
        for j=i:-1:2
            UVZYh{i}=UVZYh{i}*Dh{j}*inv(D0{j})*K;
        end
        UVZYh{i}=UVZYh{i}*Dh{1}*inv(D0{1})*[UV10(1);UV10(2);0;0];
        ABCD{i}=inv(Dh{i})*UVZYh{i};
        U{i}=(-(g1(i)^2)*s11(i)/am+am*s12(i))*exp(-g1(i)*y)*ABCD{i}(1)+(-(g1(i)^2)*s11(i)/am+am*s12(i))*exp(g1(i)*y)*ABCD{i}(2)+(-(g2(i)^2)*s11(i)/am+am*s12(i))*exp(-g2(i)*y)*ABCD{i}(3)+(-(g2(i)^2)*s11(i)/am+am*s12(i))*exp(g2(i)*y)*ABCD{i}(4);
        V{i}=(-g1(i)*s12(i)+s22(i)*am^2/g1(i))*exp(-g1(i)*y)*ABCD{i}(1)+(g1(i)*s12(i)-s22(i)*am^2/g1(i))*exp(g1(i)*y)*ABCD{i}(2)+(-g2(i)*s12(i)+s22(i)*am^2/g2(i))*exp(-g2(i)*y)*ABCD{i}(3)+(g2(i)*s12(i)-s22(i)*am^2/g2(i))*exp(g2(i)*y)*ABCD{i}(4);
        Txy{i}=am*g1(i)*exp(-g1(i)*y)*ABCD{i}(1)-am*g1(i)*exp(g1(i)*y)*ABCD{i}(2)+am*g2(i)*exp(-g2(i)*y)*ABCD{i}(3)-am*g2(i)*exp(g2(i)*y)*ABCD{i}(4);
        Sy{i}=-am^2*exp(-g1(i)*y)*ABCD{i}(1)-am^2*exp(g1(i)*y)*ABCD{i}(2)-am^2*exp(-g2(i)*y)*ABCD{i}(3)-am^2*exp(g2(i)*y)*ABCD{i}(4);
        Sx{i}=g1(i)^2*exp(-g1(i)*y)*ABCD{i}(1)+g1(i)^2*exp(g1(i)*y)*ABCD{i}(2)+g2(i)^2*exp(-g2(i)*y)*ABCD{i}(3)+g2(i)^2*exp(g2(i)*y)*ABCD{i}(4);
    end 
    
 
    x1=0*l;x2=0.5*l;
    for i=1:p%´øÈëy×ø±ê
        aSx{i}=aSx{i}+subs(Sx{i}*sin(am*x),{x,y},{x2,ya{i}});
        aSy{i}=aSy{i}+subs(Sy{i}*sin(am*x),{x,y},{x2,ya{i}});
        aTxy{i}=aTxy{i}+subs(Txy{i}*cos(am*x),{x,y},{x1,ya{i}});
        aU{i}=aU{i}+subs(U{i}*cos(am*x),{x,y},{x1,ya{i}});
        aV{i}=aV{i}+subs(V{i}*sin(am*x),{x,y},{x2,ya{i}});
    end
     vUa=[];vVa=[];vSxa=[];vSya=[];vTxya=[];yaa=[];
    for i=1:p
        yaa=[yaa ya{i}];
        vUa=[vUa aU{i}];
        vVa=[vVa aV{i}];
        vSxa=[vSxa aSx{i}];
        vSya=[vSya aSy{i}];
        vTxya=[vTxya aTxy{i}];
    end
end

  Sxf=vSxa(1);
    Txyc=vTxya(0.5*(length(vTxya)+1));
    Vc=vVa(0.5*(length(vVa)+1));
    Ann=[Sxf,Txyc,Vc];
    vpa(Ann,5)
% figure('name','Sx');hold on;plot(vSxa,yaa);%set(gca,'YDir','reverse');
% figure('name','Sy');hold on;plot(vSya,yaa);%set(gca,'YDir','reverse');
% figure('name','Txy');hold on;plot(vTxya,yaa);%set(gca,'YDir','reverse');
% figure('name','U');hold on;plot(vUa,yaa);%set(gca,'YDir','reverse');
% figure('name','V');hold on;plot(vVa,yaa);%set(gca,'YDir','reverse');   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    % figure('name','Sx');plot(ta,vSxb);%set(gca,'YDir','reverse');
%  figure('name','Sy');plot(ta,vSyb);%set(gca,'YDir','reverse');
%  figure('name','Txy');plot(ta,vTxyb);%set(gca,'YDir','reverse');
%  figure('name','U');plot(ta,vUb);%set(gca,'YDir','reverse');
%  figure('name','V');
 %hold on;plot(ta,vVb);set(gca,'xscale','log');%set(gca,'YDir','reverse');















