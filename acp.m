X = xlsread('Exemple 1.xlsx');
[l,c]=size(X)
moyenne=mean(X,1)
ecart=std(X,1)
XC=X-ones(l,1)*moyenne 
XCR=XC./(ones(l,1)*ecart) 
XCOV=cov(X,1) 
Xcoor=cov(XCR,1)
[vectpropre,valpropre]=eigs(Xcoor)
%% Nombre d’axes nécessaires
inertie=diag(valpropre)./sum(diag(valpropre))*100
inertie_cumule=cumsum(inertie)
figure(1)
bar(inertie)
hold on
plot(inertie)
disp('composante principale');
C=XCR*vectpropre
%% représentation des individus dans le plan factoriel
figure(4)
indiv = {' Aix ', 'Bec', 'Cay ', 'Cha ', 'Cri', 'Cyr', 'Evi', 'Fer',...
'Hip' ,'Lau', 'Oge','Ond','Per','Rib', 'Spa',' Tho', ' Ver','Vil','Vit',
'Vol'};
for i=1:20
plot(C(i,1),C(i,2),'*r')
text(C(i,1),C(i,2),indiv{i},'fontsize',16)
hold on
end
grid
hold on
plot(-6:0.1:6,0,'.b',0,-4:0.1:5,'.b')
%% Représentation graphique des variables dans le plan factoriel
corVF=vectpropre*sqrt(valpropre)
figure (1)
var={'HCO3-','SO4-','Cl-','Ca+','Mg+','Na+'}
for i=1:p
plot(corVF(i,1),corVF(i,2),'*b')
text(corVF(i,1),corVF(i,2),var{i},'fontsize',10)
hold on
end
hold on
t=0:0.1:2*pi;
plot (cos(t),sin(t),'b--')
grid
axis equal
%% contribution des individus
Cr_ind=100*((C(:,1:2).^2)/20)*inv(valpropre(1:2,1:2));
%% contribution des variables
Cr_var=100*vectpropre(1:2,1:2).^2;