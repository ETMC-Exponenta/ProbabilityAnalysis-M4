disp('====================');
disp('ÃËÀÂÀ 2             ');
disp('Ñëó÷àéíûå âåëè÷èíû  ');
disp('--------------------');

close all
clear all

X2=POI(2); X3=POI(3); X5=X2+X3
X=GEO(0.1,1), X50=min(X,50), X20=min(X,20)
XF=CRV('2/pi*asin(x/10)',[0,10]), Xf=CRV('2/pi./sqrt(100-x.^2)',[0,10])
Show(XF,'F:b',XF,'f:r'),
x=[0:0.5:10]; hold on, F=2/pi*asin(x/10); f=2/pi./sqrt(100-x.^2); plot(x,F,'b.',x,f,'r.')
p1=Ver(XF,[0,6]) , p2=Ver(XF,XF<6), p3=Ver(XF,XF>6)
pA=Ver(Xf, '0.8*Xf/L : L', 10)
L=10; N=5000; x=rand(1,N)*pi/2;X=sin(x)*L; [F,fs,H]=SmartHist(X); Show(F), S=CRV(F)
B=BIN(0.3,10), Y=Gen(B, 5000); [F,f]=SmartHist(Y);ShowAll(F,f,B,'ro'), b = BIN( f )
A=Rect(20,10);B=Rect(6,4);b=MySize(B);A1=A-b;A0=A+b;D=A0*1.1; SA=Area(A);N=5000;X=Gen(D,N);
for i=1:20 Z=moveTo(B,X(i));T=Sect(A,Z); Show( Z,'Hr',T,'Fc'); end
hold on, Show(A, 'k', A1, 'g-', A0, 'b-', D, 'm')
S=zeros(1,N);  for i=1:N Z=moveTo(B,X(i)); T=Sect(A,Z); S(i)=Area(T); end; u=S/SA; 
[F,f,H]=SmartHist(u); Show(H), figure, Show(F), m = mean(u), s = std(u)
P0 = 1-Area(A0)/Area(D), P1 = Area(A1)/Area(D) 
U = CRV( F ), ShowAll(F, U, 'F:r.'), P05 = Ver(U,U>0.5*0.12), p = 1-Cdf(U, [0.5 0.6 0.7]*0.12)
T=EXP(1.5), Show(T,'Fk', T, 'fr'), T3=T+T+T, figure, Show(T,'k',T+T,'r',T3,'b',T3+T,'g')
X=EXP(1.2:0.2:1.8), Y=Prodf(X), W=WEI(Y), Show(X,'Fb', Y, 'r',Y,'Fk',W,'Fr.')
p2 =1- prod(1-exp(-[1.2:0.2:1.8]*2)), pY = Ver(Y,Y>2), pW=Ver(W,W>2)
S='1 or 2 or 3 and 4 or 5 and 6 and 7 or 8 or 9 and 10';  Lam = [3, 2, 12, 15, 10, 11, 11, 1, 18, 20]/100;
E=EXP(Lam); W34=WEI(Prodf(E(3:4))); W57=WEI(Prodf(E(5:7))); W910=WEI(Prodf(E(9:10)),'F');
L=[Lam(1),Lam(2), MyPar(W34,1) , MyPar(W57,1) ,Lam(8), MyPar(W910,1)];
A=[1,1, MyPar(W34,2) , MyPar(W57,2) ,1, MyPar(W910,2)];
T=0:0.1:20; Z=[]; for t=T Z(end+1)=1-exp(-sum((L*t).^A)); end; plot(T,Z), F=CRV(Z,T),
P1=1-exp(-sum((L*1).^A)), B=Randevent(DF(E, 1)); R1=Set(B, S)
Ft=[];T=0:0.1:20;for t=T; B=Randevent(DF(E, t));Ft(end+1)=Value(Set(B,S));end, hold on, plot(T,Ft,'r--')
Elem = WEI(T,Ft), Show(Elem,'F:k.')
R=Reliab(1,{0.03 0.02, {0.12; 0.15}, {0.1; 0.11; 0.11}, 0.01, {0.18; 0.2}}), MyPar(R)
P=Reliab( 0,{0.68 0.72, 0.87 0.62 0.78, 0.65})
N=NORM(10,3), R=RND(N),ShowAll(N,R)
R1=RND(0,1),R2=R1+R1, R3=R2+R1,R4=R3+R1; R6=R4+R2; N=NORM(R6), Show(R1,R2,R3,R4,R6,N ,'.')
m=5;sigma=0.5; X=NORM(m, sigma), Show(X,'k',X,'Fr'), P3sigma = Ver(X, [m-3*sigma,m+3*sigma])
m=5; sigma=0.5; P= FLaplas((5.5-m)/sigma) - FLaplas((4.5-m)/sigma)
E= FLaplas([],0.25)
p=0.9; X09=m+[-1,1]*FLaplas([],p/2)*sigma
X=NORM(2,3); p1=Ver(X,[-1,3]), p2=Ver(X, X<3), p3=Ver(X,[-1,2,3;3,4,5]), p4=Ver(X,[-1,2,3;3,4,5]')

disp('^^^^^^^^^^^^^^^^^^^^');
disp('  ÒÅÑÒ ÇÀÂÅÐØÅÍ 2   ');
disp('vvvvvvvvvvvvvvvvvvvv');
disp('                    ');
