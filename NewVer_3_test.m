disp('===========================');
disp('√À¿¬¿ 3                    ');
disp('—ËÒÚÂÏ˚ ÒÎÛ˜‡ÈÌ˚ı ‚ÂÎË˜ËÌ  ');
disp('---------------------------');

close all
clear all

n=5;m=7;p=0.8;g=Ber(p,n);P=zeros(n+1);for j=1:n+1 P(:,j)=g(j)*Ber(1/m,j-1, (0:n)'); end, P, pX=sum(P,2)
r=0.8; G=1-(1-r).^(0:n); Wi = dot(G, pX), M = Wi*m
X=CHI2(1:5);figure(1), Show(X); Y=CHI2(10:10:30), figure(2), Show(Y), figure(3), Show(Y,'F:r')
Z = NORM(30,sqrt(60)); figure(2), Show(Z, 'r--'), figure(3), Show(Z, 'F:r--')
k =[0.5;1;2;3;4]; X=CHI2(3); p=Ver(X,k.^2); g = 2*FLaplas(k)-sqrt(2/pi).*k.*exp(-k.^2/2); p=p',g=g'
Y0=Norm_2, Y1 = Y0*[2, 1], Y2 = Y1+[12;3], Y3 = Rot(Y1,-30) + [12; 3]
R=Rect([0;0],[12;3]), [Ver(Y0,R),Ver(Y1,R),Ver(Y2,R),Ver(Y3,R)], ShowAll(R,'Fy',Y0,'k',Y1,'r',Y2,'b',Y3, 'g')
Z=Norm_2([1.01 1.06],0.48), R=RAYL(Z), Ek=MyPar(R,2); C=Circ(Ek); pN=Ver(Z,C), pR=Ver(R,Ek)
x=-10:5:10; y=-1:3; R=Rect(3,2); RR=moveTo(R, [x;y]); S=Area(R); 
X=Norm_2([5 3]); [p, Pk]=Ver(X,RR); Pk, Vk= Pdf(X,[x;y])*S, figure, ShowAll(X,'r',RR,'k')
G = Circ(10) & Rect([-15;-15],[0;0]);  X = Norm_2([5,5]);  N=1000; Pnt = Gen('lpt',G,N); 
pG = sum(dens(X, Pnt))*Area(G)/N,  p = 1 - exp(-2^2/2); p=p/4, ShowAll(Pnt, G,'Fc')
X=Norm_2([5,3]); XgXi = 'Xg=X*sqrt(r);Xi=X*sqrt(1-r);'; r=0.1; eval(XgXi), R=Rect([5,3],[1;2]);
p1=Ver(X,R), n=30; p1_30 = 1-(1-p1)^n, Xg, Xi, Z = Xg + Xi
zalp_1 = ' P=Point;Z=Xi+Gen(Xg); P=AddPoints(P, Gen(Z,n),0);'; zalp = zalp_1; m=1; N=1000;
Stat = 'L=zeros(1,n*m); for i=1:N eval(zalp), u= Impact(P,R); if u>0 L(u)=L(u)+1; end, end'; 
r=0.1;eval(XgXi); eval(Stat); Rn01=sum(L)/N
for i=1:10 ShowAll(X, R, Gen(Xi + Gen(Xg), n),'R'), hold on, end
r=0.9; eval(XgXi);eval(Stat); Rn09=sum(L)/N;
figure, for i=1:10 ShowAll(X, R, Gen(Xi + Gen(Xg),n), 'R'), hold on, end, eval(XgXi); eval(Stat), Rn09
T=[]; t=0:0.1:1; for r=t eval(XgXi), eval(Stat); T(end+1)=sum(L)/N; end, figure, plot(t,T,'k.'), T,hold on
[Fun,Par,rr] = Approx('Par(2)+sqrt(1-x.^2)*(Par(1)-Par(2))', t, T, [0.5,0.5]); Par, rr, Plot(t, Fun(t, Par),'r')
r=0.4; eval(XgXi); n=5; A=zeros(N,n); for i=1:N eval(zalp), A(i,:)=MyCenter(Gen(Z,n),[],1);  end
E=CorrelCoef(A)
X=Norm_2([5,3]); R=Rect([5,3],[1;2]); XgXi = 'Xg=X*sqrt(r);Xi=X* sqrt(r);'; m=1;n=30; r=0.6; eval(XgXi);
R30 = W2mn(Xg,Xi,R,30)
Stat = 'L=zeros(1,n*m); for i=1:N eval(zalp), u= Impact(P,R); if u>0 L(u)=L(u)+1; end, end, p_m=L/N;';
zalp = ' P=Point;Z=Xi+Gen(Xg); P=AddPoints(P, Gen(Z,n),0);'; 
N=1000; eval(Stat); R30s=sum(p_m)
zalp_m = 'P=Point; B=Gen(Xb); for j=1:m Z=Xi+Gen(Xg+B); P=AddPoints(P, Gen(Z, n), j-1); end';
XbXgXi='k=sqrt([r,1-r]); kb=sqrt([rb, 1-rb]);Xg=X*k(1);Xi=X*k(2);Xb=Xg*kb(1);Xg=Xg*kb(2);';
r=0.6; rb=0.5; X=Norm_2([5,3]); eval(XbXgXi); Xb + Xg + Xi 
m=2; n=3; rb=0.5; r=0.6; c = 1;A=zeros(N,n*m); for i=1:N eval(zalp_m), A(i,:)=MyCenter(P,[], c);end
rX=CorrelCoef(A)
Tr= Ballist('V0,c,teta0',800,0.5,45)
T=Set(Tr,'V0,c,teta0',500,0.5,45); [Xc,tetac]=Get(Tr,'Xc,tetac ')
T=Set(Tr, 'Yc!',1.2,'Xc!!',1000,'teta0?',{5,[0.2,88]});   Get(T,'teta0,tetac,Yc') 
t=[];for D=500:250:2000 T=Set(T,'Yc!',1.2,'Xc!!',D,'teta0?',{5,[0.2,88]}); t(end+1)=Get(T,'tetac');end,t
Th=Set(T, 'V0',400,'Yc!',1000,'tc!!',3,'teta0?',{34,[0.2,88]}); 
[X0, Xc,Vc,teta0,tetac,Yc]=Get(Th,'X0,Xc,Vc,teta0,tetac,Yc')
U=100;dH=U*sin(tetac*pi/180);dX=U*cos(tetac*pi/180); 
TdH=Set(Th,'Yc!',1000-dH,'Xc!', -dX, 'tc!!',3,'teta0?',{3,[0.2,88]},'X0',0);
[X0, Xc,Vc,teta0,tetac,Yc]=Get(TdH,' X0, Xc,Vc,teta0,tetac,Yc ')
V0=700; c=0.5; T=Ballist('V0,c', V0 , c); dLV0=Sensit(T,'Xc|V0', 1)
V0=800;c=0.5;T=Ballist('V0,c,teta0',V0 ,c,45);x = [5:2:47];
Lt= Clc(T,'teta0',x,'Xc');dLV0=Sensit(T,'Xc|V0',1,'teta0',x,1);plot(x,dLV0)
dLt0= Sensit(T,'Xc|teta0',0.180/pi,'teta0',x,1); dLc= Sensit(T,'Xc|c',0.01*c,'teta0',x,1);plot(x,dLt0)
rV0=0.227;rteta0=0.6;rc=1.14*c; 
B=[(dLV0*rV0).^2;(dLt0*rteta0).^2;(dLc*rc).^2];E=sum(sqrt(B)); 
v=200:25:800;Lv= Clc(T,'V0',v,'Xc'); dLV0v=Sensit(T,'Xc|V0',10,'V0',v,10)/10;
dLcv=Sensit(T,'Xc|c',0.01*c,'V0',v,1);  dLt0v=Sensit(T,'Xc|teta0',0.180/pi,'V0',v,1);   
B=[(dLV0v*rV0).^2;(dLt0v*rteta0).^2;(dLcv*rc).^2];Ev=sum(sqrt(B)); plot(Lv,Ev,'k')

disp('^^^^^^^^^^^^^^^^^^^^');
disp('  “≈—“ «¿¬≈–ÿ≈Õ 3   ');
disp('vvvvvvvvvvvvvvvvvvvv');
disp('                    ');
