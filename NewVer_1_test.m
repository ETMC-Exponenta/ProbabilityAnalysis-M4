disp('====================');
disp('√À¿¬¿ 1             ');
disp('—ÎÛ˜‡ÈÌ˚Â ÒÓ·˚ÚËˇ   ');
disp('--------------------');

close all
clear all

N=1000; R=100; M=30; P=HyperGeo(N, R, M, 0:4 )
R=HyperGeo(N, [], M, 2 )
P=para3([3,5,9]), Pr=Proect(P,1,30,3,50), Per = Perimeter(Pr), Pnt = Gen(Per,1000, 'lpt')
Show(P,'FP'), figure, Show(Pr,'k'), figure, ShowAll( Per, 'k', Pnt, 'r' )
[S, Si] = Area(Pr), A = Area(Per), ver = Si/S
M = Impact(Pnt,Pr), N=Count(Pnt), R=M/N
[M, pop] = Impact(Pnt,Pr); p2 = Count(pop{2})/Count(Pnt), P2 = M(2)/Count(Pnt)
[A,B,C]= Randevent(0.3, 0.4, 0.5); D=A+B*C, E=(A+B)*C
A=Set(A,C,0.45); B=Set(B,C,0.35); D=A+B*C, E=(A+B)*C
p=[0.2 0.4 0.45 0.25 0.5 0.1 0.6 0.15 0.4]; 
SYS = Randevent(p), A = sum(SYS),B=prod(SYS(1:4))+prod(SYS(5:9))
U = Set(SYS, '1+2+(3+4)*(5+6)+(7+8)*9')
PostH1 = 0.8*0.9/(0.8*0.9 + 0.2*0.2)
[H0,H1]=Randevent(0.2,0.8); A=SetDepend(H0,0.8,H1,0.1); H=Bayes(H1, not(A))
H=Bayes(H1, A)
H=Bayes(H1, A, not(A), not(A))
H=Bayes(H1, A, not(A), not(A), not(A))
H_10=Bayes(H1, A),             H_13=Bayes(H1, A, not(A),3),           H_24=Bayes(H1, A, 2, not(A), 4)
n=50; m=20; Dov=0.9; [p1, p2]=Ber([],n,m,Dov)
V=[554 594 572 556 546 585 605 604 549 565 573 590 580 557 584 566 583 576 558 579];
Z=[1 0 0 1 1 1 0 0 1 0 1 0 1 1 0 1 0 0 1 0];
[V,Ind]=sort(V); M=Z(Ind); M1=cumsum(1-M);M0=cumsum(M(end:-1:1));M0=M0(end:-1:1);
M=M1; N=M0+M1; P=Ber([], N, M, 0.8); f = M./N; plot(V, f,'ro--', V,P(1,:),'k--',V,P(2,:),'k--')
[Fun,Par,r]=Approx('FGauss(x, Par)', V, f ,[540,10]), hold on, plot(V,f,'ro--',V,Fun(V,Par),'k'), grid
[Fun0,Par0]=Approx('FGauss(x, Par(1:2))*Par(3)', V,P(1,:) ,[540,10,P(2,end)])
[Fun1,Par1]=Approx('FGauss(x, Par(1:2))*(1-Par(3))+Par(3)', V,P(2,:) ,[540,10,P(2,1)])
hold on, v=550:2:610, plot(v,Fun(v,Par),'k', v,Fun0(v,Par0),'g',v,Fun1(v,Par1),'r')
v=580; p=FGauss(v,Par), p0=Fun0(v,Par0), p1=Fun1(v,Par1)
p=[0.01, 0.1:0.1:0.9, 0.99]; eps = 0.1; Dov=0.9; 
nB=Ber(p,[],eps,Dov), nG = ceil(p.*(1-p)/eps^2*FLaplas([],Dov/2)^2)
R=Rect(10,10); R1=Rect([2,3],[-3;-3]); R2=Rect([1,2],[2;3]); R3=Rect(6,1); X=Norm_2([5;5],[6,5]); 
n=200; m=20; N=1000; A=zeros(N,3); M=zeros(m,3); 
for i=1:N Pnt=Gen(X,R,n); A(i,:)=Impact(Pnt,{R1, R2, R3}); end
for i=1:m M(i, :)=sum(A == i -1)/N; end
figure(1), ShowAll(R,R1,R2,R3,Pnt)
figure(2), k=[0:m-1]; plot(k, M(k+1,1), 'r--', k, M(k+1,2), 'k--', k, M(k+1,3), 'b--')
L=[Pdf(X,MyCenter(R1))*6, Pdf(X,MyCenter(R2))*2, Pdf(X,MyCenter(R3))*6]*n/Ver(X,R)
hold on, plot(k, Poi(L(1),k), 'r',  k, Poi(L(2),k),'k',  k, Poi(L(3),k))

disp('^^^^^^^^^^^^^^^^^^^^');
disp('  “≈—“ «¿¬≈–ÿ≈Õ 1   ');
disp('vvvvvvvvvvvvvvvvvvvv');
disp('                    ');
