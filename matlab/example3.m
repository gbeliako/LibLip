% this example illustrates smoothing of  noisy data
[Xp,Yp]=meshgrid(-2:.125:2);
[r,c]=size(Xp);
clear Z;
for i=1:r for j=1:c Z(i,j)=testf(Xp(i,j),Yp(i,j)); end; end;
surfc(Xp,Yp,Z)
view([113,14]);
disp 'This is a test function. We generate a few data and then approximate it'
disp 'Press any key to continue...'
pause
figure

dim=2; Ndata=200;
noise=0.02*rand(Ndata)-0.01;
XD=4*rand(2,200)-2;
clear YD;
for i=1:Ndata YD(i)=testf(XD(1,i),XD(2,i))+noise(i); end;

[Xp,Yp]=meshgrid(-2:.125:2); clear Z
[r,c]=size(Xp);
LC=mlliblip('ComputeLipschitz',dim,Ndata,XD,YD)
for i=1:r for j=1:c Z(i,j)=mlliblip('Value',dim,Ndata,[Xp(i,j),Yp(i,j)],XD,YD,LC); end; end;
surfc(Xp,Yp,Z)
hold on
plot3(XD(1,:),XD(2,:),YD(:),'x')
hold off
view([113,14]);
disp 'Next smoothen the data.'
disp 'Press any key to continue...'
pause

LC=1.2
TD=mlliblip('SmoothLipschitz',dim,Ndata,XD,YD,LC,0,0,0);
for i=1:r for j=1:c Z(i,j)=mlliblip('Value',dim,Ndata,[Xp(i,j),Yp(i,j)],XD,TD,LC); end; end;
surfc(Xp,Yp,Z)
hold on
plot3(XD(1,:),XD(2,:),YD(:),'x')
hold off
view([113,14]);

disp 'Next estimate Lipshitz constant from noisy data using sample splitting'
disp 'Press any key to continue...'
pause

Ratio=0.2;
TD=mlliblip('ComputeLipschitzSplit',dim,Ndata,XD,YD,Ratio,0);
LC=mlliblip('getlipconst')

for i=1:r for j=1:c Z(i,j)=mlliblip('Value',dim,Ndata,[Xp(i,j),Yp(i,j)],XD,TD,LC); end; end;
surfc(Xp,Yp,Z)
hold on
plot3(XD(1,:),XD(2,:),YD(:),'.')
hold off
view([113,14]);


disp 'Finally use locally Lipschitz interpolation of smoothened data.'
disp 'Press any key to continue...'
pause

mlliblip('ComputeLocalLipschitz',dim,Ndata,XD,TD)
for i=1:r for j=1:c Z(i,j)=mlliblip('ValueLocal',dim,Ndata,[Xp(i,j),Yp(i,j)],XD,TD); end; end;
surfc(Xp,Yp,Z)
hold on
plot3(XD(1,:),XD(2,:),YD(:),'x')
hold off
view([113,14]);


