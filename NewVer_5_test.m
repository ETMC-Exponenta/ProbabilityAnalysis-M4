disp('====================');
disp('Семинары             ');
disp('--------------------');

close all, clear all

a = 100; h = 10; d = 20; Om = Rect([a+h]), A = Rect([a-d]),p=Area(A)/Area(Om)
A1=A+[a+h;0];A2=A+[0;a+h]; A3=A2+[a+h;0]; 
B=Rect(a); B1=B+[a+h;0]; B2=B+[0;a+h]; B3=B2+[a+h;0];
Om1 = Om+[a+h;0];Om2=Om+[0;a+h]; Om3=Om2+[a+h;0]; 

y=XY4(Om); Y=XY4(Om3); R = Rect(y(:,1),Y(:,2)), 
Show( R,'Fk',B,'Fw', B1,'Fw', B2, 'Fw',B3,'Fw'), hold on
Show(Om,'y', Om1, 'y', Om2,'y', Om3,'y')
Show(A, 'Fc', A1,'Fr', A2, 'Fg',A3, 'Fy')
q = XY4(A);C = Circ(10,q(:,1)),Show( C, 'r')
V = Gen(A,100); for i = 1:100 g = moveTo(C,V(i)); Show(g,'Fg'),end

K=Rect([-3;0.5],[3;2]),B=Rect([-1;2],[1;3]),C1=Circ(0.5,[-2;0.5]);C2=C1+[4;0]; 
Om=K|B|C1|C2; A={C1,C2,B}; Show(K,'r',A, 'EFc')
[SA,SAi] = Area(A); P = SA/Area(Om), Pi = SAi/Area(Om), Show(Om,'Hk')

N=5000; Pnt = Gen(Om,N); [M, pop] = Impact(Pnt, A); ps =M/Count(Pnt)
ShowAll( Pnt, 'y', Pnt(pop), 'r.', A,'b', Om, 'k')
ShowAll( Pnt, 'y', pop, 'r', A,'b', Om, 'k')

[A,B,C]=Randevent(0.3,0.4,0.5); B=Set(B,C,0.45); D=A+B*C, E=(A+B)*C
p=[0.2 0.04 0.1 0.05 0.3 0.1 0.06 0.15 0.08]; L='1 + 2 + (3+4)*(5+6) + (7+8)*9'; SYS = Randevent(p);
A=Set(SYS,L), B=sum(SYS)

N=100000; R = Rect([-0.5; -0.5], [0.5;0.5]), C=Circ(0.5); Pnt = Gen(R,N); M=Impact(Pnt,C), p=M/N
N=100000; RX=rand(1,N)-0.5; RY=rand(1,N)-0.5; XY=sqrt(RX.^2 +RY.^2); 
Ind = find(XY<0.5); ps = length(Ind)/N

V=[554 594 572 556 546 585 605 604 549 565 573 590 580 557 584 566 583 576 558 579];
Z=[1 0 0 1 1 1 0 0 1 0 1 0 1 1 0 1 0 0 1 0];
[V,Ind]=sort(V); M=Z(Ind); M1=cumsum(1-M);M0=cumsum(M(end:-1:1));M0=M0(end:-1:1);
M=M1; N=M0+M1; P=Ber([], N, M, 0.8); f = M./N; plot(V, f,'ro--', V,P(1,:),'k--',V,P(2,:),'k--')
[Fun,Par,r]=Approx('FGauss(x, Par)', V, f ,[540,10]), hold on, plot(V,f,'ro--',V,Fun(V,Par),'k'), grid
[Fun0,Par0]=Approx('FGauss(x, Par(1:2))*Par(3)', V,P(1,:) ,[540,10,P(2,end)])
[Fun1,Par1]=Approx('FGauss(x, Par(1:2))*(1-Par(3))+Par(3)', V,P(2,:) ,[540,10,P(2,1)])
hold on, v=550:2:610, plot(v,Fun(v,Par),'k', v,Fun0(v,Par0),'g',v,Fun1(v,Par1),'r')
v=580; p=FGauss(v,Par), p0=Fun0(v,Par0), p1=Fun1(v,Par1)
close all, clear all

disp('====================');
disp('Семинар 2             ');
disp('--------------------');

X=CRV('F:2/pi*asin(x/10)',[0,10])
L=10; N=1000; x=rand(1,N)*pi/2;X=sin(x)*L; [Fs,fs,H]=SmartHist(X,[1,10],10); figure, Show(H,fs)
figure, [Fss,fss,H]=SmartHist(X,[1,10],20); Show(H,fss)

A=Rect(20,10);B=Rect(6,4);b=MySize(B);A1=A-b;A0=A+b;D=Rect(30,20);SA=Area(A)
C= XY4(A); Z = moveTo(B,C(:,1)); T=Sect(A, Z); 
Show(A,'2hb', A1,'k-.', A0, 'k-', D, B, 'Fc', Z, 'Fy', T,'EFr')
N=10000;X=Gen(D,N);figure, for i=1:50 Z=moveTo(B,X(i));T=Sect(A,Z);Show(Z,'Hr',T,'Fc',A,'2b'); end
S=zeros(1,N);  for i=1:N Z=moveTo(B,X(i)); T=Sect(A,Z); S(i)=Area(T); end 
S=S/SA; [F,f,H]=SmartHist(S,[],30); figure, Show(H),figure, Show(F, '1.5b')
U=CRV(F), Show(U, 'F:r.'), m = MO(U), p = Ver(U,U>m)

A=Rect(20,10); B=Circ(2); b=MySize(B); A1=A-b; A0=A+b;D=Rect(25,15);SA=Area(A)
C= XY4(A); Z = moveTo(B,C(:,1)); T=Sect(A, Z); 
Show(A,'h2b', A1,'k', A0, 'k', D, B, 'Fc', Z, T, 'EFr')
N=10000;X=Gen(D,N); for i=1:30 Z=moveTo(B,X(i));T=Sect(A,Z);figure(3),Show(Z,'Hr',T,'EFc',A); end
S=zeros(1,N);  for i=1:N Z=moveTo(B,X(i)); T=Sect(A,Z); S(i)=Area(T); end 
S=S/SA; [F,f,H]=SmartHist(S,[],30);  figure, Show(H),figure, Show(F, '1b')
U=CRV(F), Show(U, 'F:9r.'), m = MO(U), p = Ver(U,U>m)


disp('====================');
disp('Семинар 3 - 4             ');
disp('--------------------');

S='1 or 2 or 3 and 4 or 5 and 6 and 7 or 8 or 9 and 10';  Lam = [3, 2, 12, 15, 10, 11, 11, 1, 18, 20]/100; 
E=EXP(Lam);  B=Randevent( Cdf(E, 1)); R1=Set(B, S)
Ft=[];T=0:0.1:20;for t=T; B=Randevent(Cdf(E, t));Ft(end+1)=Value(Set(B,S));end
Elem = WEI(T,Ft)
R_1=Randevent(Cdf(Elem,1))
R_1=Cdf(Elem,1)
m=-1; sigma = 2; a=-3; b=2; z1=0; z2=3; X=NORM(m, sigma); Y=RND([a,b]); Z=X+Y; p = Ver(Z,[z1,z2])

X=Norm_2([1;2],[3,1.5],0.6);R=Rect(2,4); pX = Ver(X,R),ShowAll(R,'Fc',X)
[Y,alpha] = Rot(X, []); pY=Ver(Y,R)

R1=RotAx(R, alpha, [1;2]), p1 = Ver( X, R1), ShowAll(R,'Fc', R1, Y,'r',X, 'k') 
X=NORM(0,2); Y=NORM(0,3); x=Net(X,30,3);y=Net(Y,30,3);
[xx,yy]=meshgrid(x,y); Z=Pdf(X,xx).* Pdf(Y,yy); surfc(xx,yy,Z)

X=NORM(0,2); Y=NORM(0,3); 
x=Net(X,30,3);y=Net(Y,30,3);
[xx,yy]=meshgrid(x,y); 
Z=Pdf(X,xx).* Pdf(Y,yy); surfc(xx,yy,Z)

A=[100:300:1000, 1500:500:3000, 5000;  1, 1.1, 1.5, 1.8, 2.5, 3.5, 5, 6.5, 10];
Z=Norm_2P([3;4], 0.4, A);
p=[];R = Rect(5,4); LL=500:500:5000; for L = LL X=Z(L); p=[p,Ver(X,R)];end
plot(LL, p)



disp('====================');
disp('Семинар 5-6             ');
disp('--------------------');

R=RAYL(3); a=3; r=2; S = 'abs(R-a) : a'; Z=FUN(S,R,a), p0 = Ver( R, r), p1 = Ver(Z, r)
r=2; eps=1; S = 'abs(R-(Sign(R-eps)+1)/2*a) : a,eps'; Z=FUN(S,R,a,eps), p0 = Ver( R, r), p1 = Ver(Z, r)

X=Norm_2([5,3]); XgXi = 'k=sqrt([r,1-r]); Xg=X*k(1);Xi=X*k(2);'; r=0.1; eval(XgXi), Xg, Xi, Xg+Xi
r=0.5; eval(XgXi), Xg, Xi, Xg+Xi
N=10000; R=Rect([5,3],[1;2]); p1=[Ver(X,R), Ver(Xg+Xi,R), Impact(Gen(Xg+Xi,N),R)/N]
zalp = ' P=Point;Z=Xi+Gen(Xg); P=AddPoints(P, Gen(Z,n),0);'; n=30; eval(zalp)
r=0.1; eval(XgXi),eval(zalp),Show(P,'R'),hold on
Stat = 'L=zeros(1,n); for i=1:N eval(zalp), u= Impact(P,R); if u>0 L(u)=L(u)+1; end, end';
r=0.1; eval(XgXi); eval(Stat); R01=sum(L)/N, r=0.9; eval(XgXi); eval(Stat); R09= sum(L)/N
p=Ver(Xg+Xi,R), Rn = 1-(1-p)^n
t=0:0.1:1; T=[]; N=500;for r=t eval(XgXi), eval(Stat); T(end+1)=sum(L)/N; end, figure, plot(t,T,'k.'),
[Fun,Par,rr] = Approx('Par(2)+sqrt(1-x.^2)*(Par(1)-Par(2))', t, T, [0.5,0.5]); 
Par, rr, plot(t,T,'k.', t, Fun(t, Par),'r')
r=0.4; eval(XgXi); n=5; A=zeros(N,n); for i=1:N eval(zalp), A(i,:)=MyCenter(Gen(Z,n),[],1);  end
E=CorrelCoef(A)


close all, clear all


disp('====================');
disp('Семинар 7             ');
disp('--------------------');

Teta = CRV('cos(x)', [0, pi/2]), Fi = RND([0, 2*pi]) 
S = 'b*c*sin(Teta) + a*(c*abs(sin(Fi)) + b*abs(cos(Fi))).*cos(Teta) : a,b,c'; a=9; b=6; c=4;
Pr3 = FUN(S, Teta, Fi, a, b, c), Show(Pr3, 'F:r')
p = Ver(Pr3,Pr3<24*1.5)
A=(a*b*c)^(1/3); C3 = FUN(S,Teta, Fi, A, A, A, {1000}), Show(Pr3, [25, 70], 'Fk', C3, 'Fr')
R=para3(9, 6, 4); cub=para3(6, 6, 6); Sq='ProectS(G,Teta, Fi) : G'; 
Pr_3 = FUN(Sq,Teta, Fi,  R, {1000}), C_3 = FUN(Sq, Teta, Fi, cub, {1000}) , hold on, Show(Pr_3, [25, 70], 'Fr.', C_3, 'Fk.')
a=6;beta=45;K=para3(a);Pr=K*Affinor(131,beta);Show(Pr+a/2*sin(beta*pi/180),'Fp')
Z3 = FUN(Sq,Teta, Fi, Pr, {1000}), Show(Z3, [25, 80], 'Fk')
a=10; beta=pi/4; ksi=45; h=15;  brus=prism(Rect([a a beta]),100);PL=Rot(Plane([0 0 1]),3, beta, 2, ksi);
[T, brus]=Sect(brus, move(PL,[0;0;40])); [Frag, brus]=Sect(brus, move(PL,[0;0;40-h])); 
Show(move(T, [0;0;10]), Frag,'FR', move(brus,[0;0;-10]))

Sy = 'I1=-4:0.01:-a; I2=-a:0.01:a;I3=a:0.01:4;y=[I1,I2,I3]; fy=[dens(X,I1-a),dens(X,I2-a)+dens(X,I2+a),dens(X,I3+a)];';
X=NORM(0,1); a=sqrt(2/pi); eval(Sy); Y=CRV(y,fy), ShowAll(X,Y, 'r')
YY=FUN('X-a*Sign(X) : a', X,a), hold on, Show(YY)

X=NORM(10,15);d=30;h=50;b=(h+d)/2;y=0:d;ff=Pdf(X,y-b)+Pdf(X,b-y);F=Cdf(X,y-b)+1-Cdf(X,b-y);
plot(y, ff, 'k.'), grid, hold on
h=50; d=30; mx =  10; sx = 15; X = NORM(mx,sx); U=FUN('perecrit(X,h,d) : h,d',X,h,d)
[f,x] = Pdf(U); plot(x,f,'b')

Teta = CRV('cos(x)', [0, pi/2]), Fi = RND([0, 2*pi]) 
S = 'b*c*sin(Teta) + a*(c*abs(sin(Fi)) + b*abs(cos(Fi))).*cos(Teta) : a,b,c'; a=9; b=6; c=4;
Pr3 = FUN(S, Teta, Fi, a, b, c)
tic, A=(a*b*c)^(1/3); C3 = FUN(S,Teta, Fi, A, A, A,3000), Show(Pr3, [25, 70], 'Fk', C3, 'Fr'), toc

Sy = 'I1=-4:0.01:-a; I2=-a:0.01:a;I3=a:0.01:4;y=[I1,I2,I3]; fy=[dens(X,I1-a),dens(X,I2-a)+dens(X,I2+a),dens(X,I3+a)];';
X=NORM(0,1); a=sqrt(2/pi); eval(Sy); Y=CRV(y,fy), ShowAll(X,Y, 'r'); T1=[-1,1];[Ver(Y,T1),Ver(X,T1)]
a0=a; a=a0*2; eval(Sy); Y=CRV(y,fy), hold on, Show(Y, 'k'); [Ver(Y,T1),Ver(Y,T1*2)]
a=a0/2; eval(Sy); Y=CRV(y,fy), hold on, Show(Y, 'k'); [Ver(Y,T1),Ver(Y,T1*2)]
R=RAYL(3); a=3; r=2; S = 'abs(R-a) : a'; Z=FUN(S,R,a), p0 = Ver( R, r), p1 = Ver(Z, r)
Sr='I1=0: 0.1:a;I2=a:0.1:4;r=[I1,I2]; fr=[ Pdf(R,a-I1)+Pdf(R,a+I1),Pdf(R,a+I2)];';
s=1; R=RAYL(s); a=1*s, eval(Sr); Yr=CRV(r,fr), ShowAll(R,Yr, 'r'); 
[Ver(R,1),Ver(Yr,Yr<1)]
S='for i=1:n for j=1:n H(j)=Pdf(X,x1(j),y(i)/x1(j));end,g(i)=Trap(1./abs(x1).*H, x1);end';
r=0.5;X=Norm_2([10;15],[2 4],r); n=100; [X1,X2]=X12(X); x1=Net(X1, n); y=linspace(0,500, n); g=[];
eval(S),plot(y, g, 'r'),hold on, Y=X1*X2, Show(Y, 'k.')

r= -0.5; X=Norm_2([10;15],[2 4], r); [X1,X2]=X12(X); Y=X1*X2,eval(S), plot(y, g, 'r'), Show(Y, 'k.')
H=NORM(5,1); V=Fun(Y,H,'H.*Y')



disp('^^^^^^^^^^^^^^^^^^^^');
disp('  ТЕСТ Семинары ЗАВЕРШЕН   ');
disp('vvvvvvvvvvvvvvvvvvvv');
disp('                    ');
