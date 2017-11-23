function [Ex,Ey,Ez,Xmid,Ymid,Zmid]=CreatingSEMCADoutputEtot(Snapshot0,Axis0,Axis1,Axis2);
%%

%%Creating SEMCAD outout
%Input: Snapshot0,Axis0, Axis1, Axis2
%%
Ex=reshape(Snapshot0(:,1),max(size(Axis0))-1,max(size(Axis1))-1,max(size(Axis2))-1);
Ey=reshape(Snapshot0(:,2),max(size(Axis0))-1,max(size(Axis1))-1,max(size(Axis2))-1);
Ez=reshape(Snapshot0(:,3),max(size(Axis0))-1,max(size(Axis1))-1,max(size(Axis2))-1);
%%
x_axis=Axis0';
x_axis2=circshift(x_axis,1);
%creating the midpoints of the grid
Xmid=(x_axis+x_axis2)/2;
Xmid(1)=[];
 
y_axis=Axis1';
y_axis2=circshift(y_axis,1);
%creating the midpoints of the grid
Ymid=(y_axis+y_axis2)/2;
Ymid(1)=[];

z_axis=Axis2';
z_axis2=circshift(z_axis,1);
%creating the midpoints of the grid
Zmid=(z_axis+z_axis2)/2;
Zmid(1)=[];

save('SEMCAD.mat','Ex','Ey','Ez','Xmid','Ymid','Zmid')
%%


