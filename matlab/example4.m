% this example illustrates using fast evaluation method with lots of data
clear all
[Xp,Yp]=meshgrid(-2:.125:2);
disp 'Generate plenty of data to use with the fast method.'
dim=2; Ndata=20000;
XD=4*rand(2,Ndata)-2;
clear YD;
for i=1:Ndata YD(i)=testf(XD(1,i),XD(2,i)); end;

[Xp,Yp]=meshgrid(-2:.125:2); clear Z
[r,c]=size(Xp);
LC=8;
% just in case we run it again, clear the memory
mlliblip('stcfreememory');
result=mlliblip('STCbuildlipinterpolant',dim,Ndata,XD,YD,LC)

disp 'The interpolant has been constructed. Now plot the interpolant.'
disp 'Observe the speed. Press any key...'
pause
for i=1:r for j=1:c Z(i,j)=mlliblip('STCValue',[Xp(i,j),Yp(i,j)]); end; end;
surfc(Xp,Yp,Z)
%hold on
%plot3(XD(1,:),XD(2,:),YD(:),'.')
%hold off
view([113,14]);


