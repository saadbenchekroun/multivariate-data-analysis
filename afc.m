close all; clc; clear all;
%% table de données
X=[13 2 5;20 2 8;10 5 5;7 1 22]
X(:,end+1)=sum(X,2);
X(end+1,:)=sum(X,1)
[n,p]=size(X);
Total=X(n,p);
N=n-1;
P=p-1;
%% Tableau des Frequences nij/n
for i=1:n
for j=1:p
tab(i,j)=X(i,j)/X(n,p);
end
end
tab
tabpct=tab.*100
%% Tableau des frequences lignes nij/ni.=fij/fi.
for i=1:N
for j=1:P
Xl(i,j)=(X(i,j)/X(i,p));
end
end
Xl
Xlpct=Xl.*100
%% Tableau des frequences colonnes nij/n.j=fij/f.j
for i=1:N
for j=1:P
Xc(i,j)=(X(i,j)/X(n,j));
end
end
Xc
Xcpct=Xc.*100
%% table d'effictif théorique
for i=1:N
for j=1:P
Teff_th(i,j)=tab(n,j)*tab(i,p);
end
end
%% matrice de taux de liaison
for i=1:N
for j=1:P
T(i,j)=(tab(i,j)-Teff_th(i,j))/Teff_th(i,j);
end
end
T=T*100
%% distance entre deux lignes
dl2=0;
for j=1:P
dl2=dl2+(Xl(1,j)-Xl(2,j))^2/tab(n,j);
end
dl2
dl=sqrt(dl2)
%% distance ligne et ligne moy
dm2=0;
for j=1:P
dm2=dm2+((Xl(1,j)-tab(n,j))^2)/tab(n,j);
end
dm2
dm=sqrt(dm2)
%% distance entre deux colonnes
dc2=0;
for i=1:N
dc2=dc2+(Xc(i,1)-Xc(i,2))^2/tab(i,p);
end
dc2
dc=sqrt(dc2)
%% matrice de taux de liaison
for i=1:N
for j=1:P
T(i,j)=(tab(i,j)-tab(n,j)*tab(i,p))/(tab(n,j)*tab(i,p));
end
end
T
%% indice de Chi 2 X=PHI2
PHI2=0;
for i=1:N
for j=1:P
PHI2=PHI2+((tab(i,j)-(tab(n,j)*tab(i,p)))^2)/(tab(n,j)*tab(i,p));
end
end
PHI2=total*PHI2
%% ****** Analyse AFC *********
% Matrice transformé
for i=1:N
for j=1:P
XT(i,j)=Xl(i,j)/sqrt(tab(n,j));
YT(i,j)=Xc(i,j)/sqrt(tab(i,p));
end
end
%% centre de gravité de X
Gx=sqrt(tab(n,1:end-1))
Gy=sqrt(tab(1:end-1,p))
%% matrice transformé centré et p
for i=1:N
XTc(i,:)=XT(i,:)-Gx
end
for j=1:P
YTc(:,j)=YT(:,j)-Gy
end
%% matrice de variance
PP=diag(tab(1:end-1,p))
C=sqrt(PP)*XTc
V=C'*C
%% Matrice d’inertie S
for i=1:N
B(i,:)=XT(i,:)*sqrt(tab(i,p));
end
S=B'*B
S1=B*B'
%% vecteurs propres et valeurs propres de S
[vectpX valpX]=eigs(S)
valpX=diag(valpX)
valpX1=valpX(2:end);
vectpX1=vectpX(:,2:end)
%% vecteurs propres et valeurs propres de S*
[vectpY valpY]=eigs(S1)
valpY=diag(valpY)
valpY1=valpY(2:end);
vectpY1=vectpY(:,2:end)
%% inertie de X et Y
inertieX=valpX1./sum(valpX1)*100
inertieCX=inertieX*[1 0]+cumsum(inertieX)*[0 1]
inertieY=valpY1./sum(valpY1)*100
inertieCY=inertieY*[1 0]+cumsum(inertieY)*[0 1]
figure (1)
subplot(1,2,1), bar(inertieX),title('Inertie de X')
subplot(1,2,2), bar(inertieY),title('Inertie de Y')
%% Coordonnées des nuages projetés sur les axes factoriels
CX=XTc*vectpX1
%% relations de transition
vectpY2=B*vectpX/diag(sqrt(valpX))
DY1=zeros(P,2)
for k=1:2
for j=1:P
for i=1:N
DY1(j,k)=DY1(j,k)+(tab(i,j)*CX(i,k))/tab(n,j);
end
DY1(j,k)=(1/sqrt(valpX1(k)))*DY1(j,k)
end
end
DY=DY1
%% Représentation des nuages lignes et colonnes
figure (3)
indiv = {'L', 'Eco', 'EX ', 'Tech'};
for i=1:N
plot(CX(i,1),CX(i,2),'sr')
text(CX(i,1),CX(i,2),indiv{i},'Color','blue','fontsize',12)
hold on
end
grid
hold on;
var={'Univ','Prepa','BTS'};
for i=1:P
plot(DY(i,1),DY(i,2),'*b')
text(DY(i,1),DY(i,2),var{i},'Color','red','fontsize',12)
hold on
end
plot(-0.8:0.01:0.8,0,'.b',0,-0.8:0.01:0.6,'.b')
%% Contribution du profil-ligne i à l’inertie de l’axe k
for i=1:N
for k=1:2
CTRX(i,k)=100*tab(i,p)*(CX(i,k)^2)/valpX1(k);
end
end
%% Contribution du profil-colonne j à l’inertie de l’axe K
for j=1:P
for k=1:2
CTRY(j,k)=100*tab(n,j)*(DY(j,k)^2)/valpY1(k);
end
end
%% reconstitution des données
phy=zeros(N,P)
for i=1:N
for j=1:P
for k=1:2
phy(i,j)=phy(i,j)+tab(i,p)*tab(n,j)*CX(i,k)*DY(j,k)/sqrt(valpX1(k))
end
t(i,j)=phy(i,j)+tab(i,p)*tab(n,j)
end
end
%%erreur de reconstruction
err_abs=abs(t-tab(1:N,1:P))
err_rel=norm(t-tab(1:N,1:P))/norm(t)
%% Qualité de représentation du point-ligne i sur l’axe k
for i=1:N
for k=1:2
qltX(i,k)=(CX(i,k)^2)/(norm(XTc(i,:)))^2;
end
end
%%Qualité de représentation du point-colonne j sur l’axe k
for j=1:P
for k=1:2
qltY(j,k)=(DY(j,k)^2)/(norm(YTc(:,j)))^2;
end
end