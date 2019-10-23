disp('===========================');
disp('√À¿¬¿ 4                    ');
disp('‘ÛÌÍˆËË ÒÎÛ˜‡ÈÌ˚ı ‚ÂÎË˜ËÌ  ');
disp('---------------------------');

close all
clear all

X=NORM(5,2), Y=RND([4,6]), U=X*2+Y*3, V=X*Y, W=(X-MO(X))*(Y-MO(Y))
y = [0.1, 0.3, 0.9]; X=NORM(1,2); x=Inv(X, y), y = Cdf(X, x)
Sy='x1=linspace(-4,-a,50); x0=linspace(-a,a,50);x2=linspace(a,4,50);x=[x1,x0,x2];';
Sf ='eval(Sy), fy=[Pdf(X, x1-a),Pdf(X, x0-a)+Pdf(X, x0+a),Pdf(X, x2+a)]; Y=CRV(x,fy);'; 
aD=sqrt(2/pi); T=[-1,1]*0.75;
X=NORM(0,1); a=aD; eval(Sf), Dmin =SKO(Y)^2, ShowAll(X,'g', Y,'r'), P0=Ver(X,T), Popt =Ver(Y,T)
a=aD/2; eval(Sf), hold on, Show(Y,'r'), P1 =Ver(Y,T), a=aD*2; eval(Sf), Show(Y,'r'), P2 =Ver(Y,T)
A=aD/2 :0.02:aD*1.5; P=[];for a=A eval(Sf),P(end+1) = Ver(Y,T); end, figure, plot( A,P) 
[pmax,imax] = max(P); aP = A(imax), Pmax = P(imax)
S = 'X-sign(X)*a : a'; Y1 = FUN(S,X, aD); popt1=Ver(Y1, T), Y2 = FUN(S,X, aP); popt2=Ver(Y2, T) 
Sy='I1=-4:0.01:-a+e; I2=-a+e:0.01:-e;I3=-e:0.01:e;I4=e:0.01:a-e;I5=a-e:0.01:4;I=[I1,I2,I3,I4,I5];';
Sf='eval(Sy),fy=[f(X,I1-a),f(X,I2-a)+f(X,I2+a),f(X,I3-a)+f(X,I3)+f(X,I3+a),f(X,I4-a)+f(X,I4+a),f(X,I5+a)];Y=CRV(I,fy);  ';
X=NORM(0,1); a=sqrt(2/pi); e=0.3;T; eval(Sf); p=Ver(Y,T), Show(Y,'r');
S = 'X-Sign(X,eps)*a : a,eps'; Y1=FUN(S, X, aD, eps); p=Ver(Y1,T), hold on,ShowAll(X,'g', Y1,'r')
f=@(x,p1,p2) -Ver(FUN('X-Sign(X,eps)*a : a,eps', X, x, p1), p2);
[xmin,fmin] = fminsearch(@(x) f(x, eps, T),0.7); a_eps= xmin,P_eps =-fmin
f2=@(x,p1)    -Ver(FUN('X-Sign(X,eps)*a : a,eps', X, x(1), x(2)), p1);
[xmin,fmin] = fminsearch(@(x) f2(x, T),[0.7,0.3]); a_eps= xmin(1), eps= xmin(2) , P_eps =-fmin
Sy='I1=0:0.001:eps; I2=eps:0.001:a-eps;I3=a-eps:0.1:4;I=[I1,I2,I3];';
Sf='eval(Sy); fy=[Pdf(R,I1)+Pdf(R,I1+a)+Pdf(R,a-I1),Pdf(R,I2+a)+Pdf(R,a-I2),Pdf(R,a+I3)]; Y=CRV(I, fy);';
R=RAYL(1.5), a=2;eps = a/3; eval(Sf);P=Ver( Y, Y<1), ShowAll(R,Y)
R=RAYL(1.5); a=2;eps=a/3; S= 'abs(R-Sign(R,eps)*a) : a,eps';  Z=FUN(S, R, a, eps); p=Ver(Z,Z<1), Show(Z)
X=Norm_2([1.5,1.5]);Pnt = Gen(X,20000);A=Modul(Pnt);F = SmartHist(A); U=CRV(F),figure,Show(U)
a=2;eps=a/3; S= 'abs(U-Sign(U,eps)*a) : a,eps';  Y=FUN(S, U, a, eps); p=Ver(Y,Y<1), Show(Y)
f2=@(x,p1) -Ver(FUN('abs(Y-Sign(Y,eps)*a) : a,eps', Y, x(1), x(2)), Y<p1);T=1;
[xmin,fmin] = fminsearch(@(x) f2(x, T),[2,0.5]); a_eps= xmin(1), eps= xmin(2), P_eps =-fmin
a=2;eps=a/3;X=Norm_2([1.5,1.5]);Pnt = Gen(X,20000);A=Modul(Pnt);R = abs(A-Sign(A,eps)*a);
F = SmartHist(R); Y=CRV(F), P=Ver(Y,Y<1), Show(Y)
Ind = find(R<1); p=length(Ind)/Count(Pnt), m = mean(R), s = std(R)
Teta = CRV('cos(x)', [0, pi/2]), Fi = RND([0, 2*pi]) 
S = 'b*c*sin(Teta) + a*(c*abs(sin(Fi)) + b*abs(cos(Fi))).*cos(Teta) : a,b,c'; a=9; b=6; c=4;
Pr3 = FUN(S, Teta, Fi, a, b, c)
A=(a*b*c)^(1/3); C3 = Fun(Teta, Fi,S, A, A, A), Show(Pr3, [25, 70], 'Fk', C3, 'Fr')
R=para3(9, 6, 4); cub=para3(6, 6, 6); Sq='ProectS(G,Teta, Fi) : G'; 
Pr_3 = Fun(Teta, Fi, Sq, R), C_3 = Fun(Teta, Fi, Sq, cub) , Show(Pr_3, [25, 70], 'Fr.', C_3, 'Fk.')
a=6;beta=45;K=para3(a);Pr=K*Affinor(131,beta);Show(Pr+a/2*sin(beta*pi/180),'Fp')
Z3 = Fun(Teta, Fi, Sq, Pr), Show(Z3, [25, 80], 'Fk') %>> Z3 = FUN(Sq,Teta, Fi, Pr), Show(Z3, [25, 80], 'Fk')
Si = Facet(Pr), mo = sum(Si)/4

a=10; beta=pi/4; ksi=45; h=15;  brus=prism(Rect([a a beta]),100);PL=Rot(Plane([0 0 1]),3, beta, 2, ksi);
[T, brus]=Sect(brus, move(PL,[0;0;40])); 
[Frag, brus]=Sect(brus, move(PL,[0;0;40-h])); Show(move(T, [0;0;10]), Frag,'FR', move(brus,[0;0;-10]))
Fs = FUN('ProectS(G,Teta, Fi) : G', Teta, Fi, Frag);
Sa='Mina(frag, Teta,Fi) : frag'; Fa = FUN(Sa, Teta, Fi, Frag);
r=     0; X=Norm_2([10;15],[2 4], r); [X1,X2]=X12(X); Y1=X1*X2; Show(Y1, 'b')
r=  0.5; X=Norm_2([10;15],[2 4], r); [X1,X2]=X12(X); Y2=X1*X2; Show(Y2, 'r')
r= -0.5; X=Norm_2([10;15],[2 4], r); [X1,X2]=X12(X); Y3=X1*X2; Show(Y3, 'k')
X1=NORM(10,2); X2=NORM(15,4);Z=X1*X2; Show(Z,'b.')
X=Norm_2([2 1]); [X1,X2]=X12(X), A=Fun(X1,X2,'atan(X2./X1)'), Show(A,'r.'),hold on
X1=EXP(2); X3=EXP(3);  Y1=X1-X1, Y2=X1-X3, Show(Y1,'b',Y2,'r')
X1=POI(1); X2=POI(2); X3=POI(3); X=X1+X2+X3
m=-1; a=-3; b=2; z1=0; z2=3;
P=[]; for sigma = 0.5:0.5:5  P(end+1) = Ver(NORM(m, sigma) + RND([a,b]), [z1, z2]);  end, P
x=[-5:0.1:5]; F=0.5+Laplas(x); f=f_Gauss(x); n=20; 
plot(x,[f; f.*n.*(1-F).^(n-1); f.*n.*F.^(n-1)],'r'), hold on,plot(x,[F;1-(1-F).^n; F.^n],'b')
N=NORM(0,1),X=N(ones(1,20)),Ymin=min(X),Ymax=max(X),Show(Ymin,'F:r.',Ymax,'F:k.')
R1=RAYL([1,1,1,1]);R12=RAYL([1,2]),min1=min(R1),min12=min(R12) ,max1=max(R1)
E123=EXP([1,2,3]); Emin123=min(E123), Emax123=max(E123)

disp('^^^^^^^^^^^^^^^^^^^^');
disp('  “≈—“ «¿¬≈–ÿ≈Õ 4   ');
disp('vvvvvvvvvvvvvvvvvvvv');
disp('                    ');
