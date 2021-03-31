disp('===========================');
disp('ГЛАВА 4                    ');
disp('Функции случайных величин  ');
disp('---------------------------');

close all
clear all

X=NORM(5,2), Y=RND([4,6]), U=X*2+Y*3, V=X*Y
A = NUM( 10,20) , B = NUM( [1;2;3], [4, 9, 16]), C=NUM(X,Y,U,V)
D=A; A_D=A-D, AD=NUM(A,D), ADS = Set(AD,0.6), Sad=sum(AD), Sads=sum(ADS), Pads = prod(AD)
Br=Set(B,1,2,0.5, 3,2,0.4), r12=Correl(Br,1,2)
X0 = NUM(X), Y0 = NUM(Y), V0 = NUM(V), XYV = sum(X0,Y0,V0), xyv = prod(X0,Y0,V0)
BS=sum(B), Bs  =sum(B*[0.1,0.2,0.3]), Br = Set(B,1,2,0.7, 2,3,0.4), SBr=sum(Br)
pB=prod(B), pBr=prod(Br)
clear all

y = [0.1, 0.3, 0.9]; X=NORM(1,2); x=Inv(X, y), y = Cdf(X, x)
X=RND(0, pi/2); L=10; Y = FUN('a*cos(x) : a', X, L), p=Ver(Y,Y<L/2)
Sy='x1=linspace(-4,-a,50); x0=linspace(-a,a,50);x2=linspace(a,4,50);x=[x1,x0,x2];';
Sf ='eval(Sy), fy=[Pdf(X, x1-a),Pdf(X, x0-a)+Pdf(X, x0+a),Pdf(X, x2+a)]; Y=CRV(x,fy);'; 
aD=sqrt(2/pi); T=[-1,1]*0.75;
X=NORM(0,1); a=aD; eval(Sf), Dmin =SKO(Y)^2, ShowAll(X,'g', Y,'r'), P0=Ver(X,T), Popt =Ver(Y,T)
a=aD/2; eval(Sf), hold on, Show(Y,'r'), P1 =Ver(Y,T), a=aD*2; eval(Sf), Show(Y,'r'), P2 =Ver(Y,T)
A=aD/2 :0.02:aD*1.5; P=[];for a=A eval(Sf),P(end+1) = Ver(Y,T); end, figure, plot( A,P) 
[pmax,imax] = max(P); aP = A(imax), Pmax = P(imax)
S = 'X-sign(X)*a : a'; Y1 = FUN(S,X, aD); popt1=Ver(Y1, T), Y2 = FUN(S,X, aP); popt2=Ver(Y2, T) 
Sy = 'I1=-4:0.01:-a+e; I2=-a+e:0.01:-e;I3=-e:0.01:e;I4=e:0.01:a-e;I5=a-e:0.01:4;I=[I1,I2,I3,I4,I5];';
Sf = 'eval(Sy), fy = [f(X,I1-a),f(X,I2-a)+f(X,I2+a),f(X,I3-a)+f(X,I3)+f(X,I3+a), f(X,I4-a)+f(X,I4+a), f(X,I5+a)];';
X = NORM(0,1); a=sqrt(2/pi); e=0.3; T=[-1,1]*0.75; eval(Sf), Y=CRV(I,fy); 
S = 'X-Sign(X,eps)*a : a,eps'; Z = FUN(S,X, a, 0.3); pY=Ver(Y,T), pZ=Ver(Z, T), Show(Y,'r', Z, 'k')
f1=@(x,p0,p1) -Ver(FUN(S, X, x, p0), p1); eps = 0.4; x0 = 1;
[xmin,fmin] = fminsearch(@(x) f1(x, eps, T), x0); a_eps= xmin,P_eps = -fmin

f2=@(x,p1) -Ver(FUN(S, X, x(1), x(2)), T); x0 = [1, 0.5];
[xmin,fmin] = fminsearch(@(x) f2(x, T),x0); a_opt =  xmin(1), e_opt =  xmin(2), P_opt = -fmin    
Sy='I1=0:0.001:eps; I2=eps:0.001:a-eps;I3=a-eps:0.1:4;I=[I1,I2,I3];';
Sf='eval(Sy); fy=[Pdf(R,I1)+Pdf(R,I1+a)+Pdf(R,a-I1),Pdf(R,I2+a)+Pdf(R,a-I2),Pdf(R,a+I3)]; Y=CRV(I, fy);';
R=RAYL(1.5), a=2;eps = a/3; eval(Sf);Y = CRV('f:',I,fy), P=Ver( Y, Y<1), figure, ShowAll(R,Y), hold on
R=RAYL(1.5); a=2;eps=a/3;S='abs(R-Sign(R,eps)*a) : a,eps'; Z=FUN(S,R, a, eps); p=Ver(Z,Z<1), Show(Z,'r')
X=Norm_2([1.5,1.5]); Pnt = Gen(X,20000); A=Modul(Pnt); F = SmartHist(A); U=CRV(F), figure, Show(U)
a=2;eps=a/3; S= 'abs(U-Sign(U,eps)*a) : a,eps';  Y=FUN(S, U, a, eps); p=Ver(Y,Y<1), Show(Y)
f2=@(x,p1)    -Ver(FUN('abs(Y-Sign(Y,eps)*a) : a,eps', Y, x(1), x(2)), Y<p1);T=1;
[xmin,fmin] = fminsearch(@(x) f2(x, T),[2,0.5]); a_eps= xmin(1), eps= xmin(2) , P_eps =-fmin
a=2;eps=a/3;X=Norm_2([1.5,1.5]);Pnt = Gen(X,20000);A=Modul(Pnt);R = abs(A-Sign(A,eps)*a);
F = SmartHist(R); Y=CRV(F), P=Ver(Y,Y<1), Show(Y)
Ind = find(R<1); p=length(Ind)/Count(Pnt), m = mean(R), s = std(R)
clear all
Teta = CRV('cos(x)', [0, pi/2]), Fi = RND([0, 2*pi]) 
S = 'b*c*sin(Teta) + a*(c*abs(sin(Fi)) + b*abs(cos(Fi))).*cos(Teta) : a,b,c'; a=9; b=6; c=4;
Pr3 = FUN(S, Teta, Fi, a, b, c, 1000)
A=(a*b*c)^(1/3); C3 = FUN(S, Teta, Fi, A, A, A,1000), Show(Pr3, [25, 70], 'Fk', C3, 'Fr')
R=para3(9, 6, 4); cub=para3(6, 6, 6); Sq='ProectS(G,Teta, Fi) : G'; 
Pr_3 = FUN(Sq, Teta, Fi, R,500), C_3 = FUN(Sq, Teta, Fi, cub,500) , Show(Pr_3, [25, 70], 'Fr.', C_3, 'Fk.')
beta = 45;Pr=cub*Affinor(131,beta); Show(Pr+a/2*sin(beta*pi/180),'Fp')
Z3 = FUN(Sq,Teta, Fi, Pr,1000), Show(Z3, [25, 80], 'Fk')
Si = Facet(Pr), mo = sum(Si)/4

a=10; beta=pi/4; ksi=45; h=15;  brus=prism(Rect([a a beta]),100);PL=Rot(Plane([0 0 1]),3, beta, 2, ksi);
[T, brus]=Sect(brus, move(PL,[0;0;40])); 
[Frag, brus]=Sect(brus, move(PL,[0;0;40-h])); Show(move(T, [0;0;10]), Frag,'FR', move(brus,[0;0;-10]))
Fs = FUN('ProectS(G,Teta, Fi) : G', Teta, Fi, Frag,500);
Sa='Mina(frag, Teta,Fi) : frag'; Fa = FUN(Sa, Teta, Fi, Frag,500);
for r= [0, 0.5, -0.5] X=Norm_2([10;15],[2 4], r); [X1,X2]=X12(X); Y=X1*X2; Show(Y), hold on, end
X1=NORM(10,2); X2=NORM(15,4);Z=X1*X2; Show(Z,'b.')
clear all

x=[-1:0.01:1]*pi/2; C='krb'; K=[1, 2, 0.5]; 
for i= 1:3 k=K(i); y=k./(cos(x).^2+k^2*sin(x).^2)/pi; c=C(i); plot(x,y,c), hold on, end
X=Norm_2([2 1]); [X1,X2]=X12(X)
A1=FUN('atan(X2./X1)', X1,X2), A2=FUN('atan(X1./X2)', X1,X2), X1=X2; A3=FUN('atan(X2./X1)', X1,X2)
Show(A1,'r.',A2,'b.', A3,'k.')

L1=1; L2=2; X1=EXP(L1); X2=EXP(L2); X201 = EXP(2.01);  Y1 = X2+X2,  Y2 = X2+X201, Y3=X1+X2
y=0:0.1:5;fy=L1*L2/(L2-L1)*(exp(-L1*y)-exp(-L2*y)); Show(Y1,'b',Y2,'r.',Y3,'k'), hold on, plot(y, fy, 'k.')
X1=EXP(2); X3=EXP(3);  Y1=X1-X1,  Y2=X1-X3, Show(Y1,'b',Y2,'r')
X1=POI(1); X2=POI(2); X3=POI(3); X=X1+X2+X3
m=-1; a=-3; b=2; z1=0; z2=3;
P=[]; for sigma = 0.5:0.5:5  P(end+1) = Ver(NORM(m, sigma) + RND([a,b]), [z1, z2]);  end, P
sigma=3;X=NORM(m, sigma); Y=RND([a,b]); ShowAll(X,'r',Y,'k',X+Y),hold on,plot([z1,z2],[0,0])
x=[-5:0.1:5]; F=0.5+Laplas(x); f=f_Gauss(x); n=20; 
plot(x,[f; f.*n.*(1-F).^(n-1); f.*n.*F.^(n-1)],'r'), hold on,plot(x,[F;1-(1-F).^n; F.^n],'b')


x=[-5:0.1:5]; F=0.5+Laplas(x); f=f_Gauss(x); n=20; 
plot(x,[f; f.*n.*(1-F).^(n-1); f.*n.*F.^(n-1)],'r'), hold on,plot(x,[F;1-(1-F).^n; F.^n],'b')
N=NORM(0,1),  X=N(ones(1,20)), Ymin=min(X), Ymax=max(X), Show(Ymin,'F:r.',Ymax,'F:k.')
R1=RAYL([1,1,1,1]);R12=RAYL([1,2]),min1=min(R1),min12=min(R12) ,max1=max(R1)
E123=EXP([1,2,3]); Emin123=min(E123), Emax123=max(E123)

X=rand(3,100000)-0.5;L=[10;15;20];R=MAdd(X,L);g=prod(R); m=mean(g), s=std(g)
A=RND(10+[-0.5,0.5]),B=RND(15+[-0.5,0.5]),C=RND(20+[-0.5,0.5]), V=A*B*C
X=NORM(0, 0.3, [-0.5,0.5]); AN = X + 10, BN=X+15, CN=X+20, VN = AN*BN*CN, W = V + VN
Show( A, 'k', AN, 'r'), figure, Show( V, 'k', VN, 'r'), figure, Show( W, 'r', NORM(6000, 159), 'k.')
aN = NUM(AN),bN = NUM(BN), cN = NUM(CN), vN=aN*bN*cN

disp('^^^^^^^^^^^^^^^^^^^^');
disp('  ТЕСТ ЗАВЕРШЕН 4   ');
disp('vvvvvvvvvvvvvvvvvvvv');
disp('                    ');
