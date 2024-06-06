%% Tableau disjonctif complet:
T =[1 0 1 0 1 0 0; 0 1 0 1 0 1 0; 0 1 0 1 0 0 1; 1 0 0 1 1 0 0; 0 1 1 0 0 1 0; 1 0 1 0 0 0 1];
T(:, end+1)=sum (T,2)
T (end+1,:)=sum(T,1)
%%﻿ Tableau des fréquences associé à T:
[n,p]=size (T);
Total T (n,p);
for i=1:n
for j=i:p
F(i,j)=T(i,j)/T (n,p);
end
end
%% Tableau de profils lignes associés à T :
for i=1:N
for j=1:P
1(i,j)=(T(i,j)/T(i,p));
end
end
%% Tableau de profils colonnes associés à T:
for i=1:N
for j=1:P
C(i,j)=((i,j)/T (n,j));
end
end
﻿Cpst =C*100
﻿%% Tableau de burt
t = [1 0 1 0 1 0 0 0 1 0 1 0 1 0:0 1 0 1 0 0 1:1 0 0 1 1 0 0:0 1 1 0 0 1 0:1 0 1 0 0 0 11];
tb=t'*t
%% Caractéristiques des nuages de points N(X) et N(Y)
﻿for i=1:N for j=1:P
XT(i,j)=1(i,j)/sqrt (F(n,j)):
YT(i,j)=C(i,j)/sqrt(F(i,p));
end
end
%% -Centre de gravité des nuages
﻿Gx=sqrt(F (n, 1: end-1))
Gy=sqrt(F(1:end-1,p))
%% -Décomposition factorielle du tableau disjonctif :
%%transforme centre
﻿for i=1:N
XTC (i,:)=XT (i,:)-GX
end
for j=1:P
YTC (:,j)=YT (:,j)-Gy
%%mat de variance
PP=diag (F (1:end-1,p))
C=sqrt (PP) *XTC
V=C'*C
%% -Matrice dinertie S
﻿for i=1:N
B(i,:)=XT (i,:)*sqrt (F(i,p));
end
S=B**B
%% -Vecteurs propres et valeurs propres de S :
﻿[vectpT valpT]=eigs (S) 
valpT=diag (valpT)
valpT1=valpT (2:end);
vectpT1=vectpT (:,2:end)