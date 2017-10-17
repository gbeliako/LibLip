load aluminiumdata.txt;
Y=aluminiumdata(:,3);
XD=aluminiumdata(:,1:2)'; % we use transpose to convert to C storage convention
plot3(XD(1,:),XD(2,:),Y(:),'x');
disp 'Data loaded. Press any key to continue...';
pause
figure
dim = 2;
mlliblip('ComputeLocalLipschitz',2,60,XD,Y);
[Xp,Yp]=meshgrid(-2.3:.02:0,-0.07:.02:1.13); clear Z
[r,c]=size(Xp);
for i=1:r for j=1:c Z(i,j)=mlliblip('ValueLocal',2,60,[Xp(i,j),Yp(i,j)],XD,Y); end; end;
surfc(Xp,Yp,Z)
hold on
plot3(XD(1,:),XD(2,:),Y(:),'x')
hold off
disp 'Data interpolated. But monotonicity does not hold. '
disp 'The next command enforces monotonicity. Press any key to continue...';
pause
figure
%view([110,40])
cons=[1,1];
for i=1:r for j=1:c Z(i,j)=mlliblip('ValuelocalCons',2,60,cons,[Xp(i,j),Yp(i,j)],XD,Y); end; end;
surfc(Xp,Yp,Z)
hold on
plot3(XD(1,:),XD(2,:),Y(:),'x')
hold off

