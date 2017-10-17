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
XD=4*rand(2,200)-2;
clear YD;
for i=1:Ndata YD(i)=testf(XD(1,i),XD(2,i)); end;

[Xp,Yp]=meshgrid(-2:.125:2); clear Z
[r,c]=size(Xp);
for i=1:r for j=1:c Z(i,j)=mlliblip('Value',dim,Ndata,[Xp(i,j),Yp(i,j)],XD,YD,2); end; end;
surfc(Xp,Yp,Z)
hold on
plot3(XD(1,:),XD(2,:),YD(:),'x')
hold off
view([113,14]);
disp 'Next calculate the Lipschitz constant from the data.'
disp 'Press any key to continue...'
pause


LC=mlliblip('ComputeLipschitz',dim,Ndata,XD,YD)
%LC=mlliblip('getlipconst')
for i=1:r for j=1:c Z(i,j)=mlliblip('Value',dim,Ndata,[Xp(i,j),Yp(i,j)],XD,YD,LC); end; end;
surfc(Xp,Yp,Z)
hold on
plot3(XD(1,:),XD(2,:),YD(:),'x')
hold off
view([113,14]);
disp 'Next calculate the local Lipschitz constants from the data.'
disp 'Press any key to continue...'
pause

mlliblip('ComputeLocalLipschitz',dim,Ndata,XD,YD)
for i=1:r for j=1:c Z(i,j)=mlliblip('ValueLocal',dim,Ndata,[Xp(i,j),Yp(i,j)],XD,YD); end; end;
surfc(Xp,Yp,Z)
hold on
plot3(XD(1,:),XD(2,:),YD(:),'x')
hold off
view([113,14]);


