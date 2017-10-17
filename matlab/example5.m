% this example illustrates using interpolation and smoothing with bounds
dim=2; Ndata=50;
XD=rand(2,Ndata);
clear YD;
for i=1:Ndata YD(i)=testf1(XD(1,i),XD(2,i)); end;

[Xp,Yp]=meshgrid(0:.05:1); clear Z
[r,c]=size(Xp);
LC=1.5;

for i=1:r for j=1:c Z(i,j)=mlliblip('Value',dim,Ndata,[Xp(i,j),Yp(i,j)],XD,YD,LC); end; end;
surfc(Xp,Yp,Z)
hold on
plot3(XD(1,:),XD(2,:),YD(:),'.')
hold off
%view([113,14]);
disp 'Next interpolate the data using specified bounds.'
disp 'Press any key to continue...'
pause


clear Z;
% bounds can be supplied in either way as below
mlliblip('setbounds','boundup',@boundlow);
for i=1:r for j=1:c Z(i,j)=mlliblip('Value',dim,Ndata,[Xp(i,j),Yp(i,j)],XD,YD,LC); end; end;
surfc(Xp,Yp,Z)
hold on
plot3(XD(1,:),XD(2,:),YD(:),'.')
hold off
%view([113,14]);

disp 'Next Smoothen the data. Press any key to continue...'
pause

LC=1.0;
TD=mlliblip('SmoothLipschitz',dim,Ndata,XD,YD,LC,0,0,0);
for i=1:r for j=1:c Z(i,j)=mlliblip('Value',dim,Ndata,[Xp(i,j),Yp(i,j)],XD,YD,LC); end; end;
surfc(Xp,Yp,Z)
hold on
plot3(XD(1,:),XD(2,:),TD(:),'.')
hold off

mlliblip('clearbounds');
%disp 'Press any key to continue...'
%pause

