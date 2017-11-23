function [Absolute_Modulus_of_B1__0s,Absolute_Modulus_of_B1p_0s,Phase_of_B1__0s,Phase_of_B1p_0s,Xmid,Ymid,Zmid]=CreatingSEMCADB1output(Snapshot0,Axis0,Axis1,Axis2);
%%

%%Creating SEMCAD outout
%Input: Snapshot0,Axis0, Axis1, Axis2
BpField=reshape(Snapshot0(:,1),max(size(Axis0))-1,max(size(Axis1))-1,max(size(Axis2))-1);
%The B1plus field
BmField=reshape(Snapshot0(:,2),max(size(Axis0))-1,max(size(Axis1))-1,max(size(Axis2))-1);
%The B1minus field


Absolute_Modulus_of_B1__0s=abs(BmField);
Absolute_Modulus_of_B1p_0s=abs(BpField);
Phase_of_B1__0s=angle(BmField);
Phase_of_B1p_0s=angle(BpField);
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

save('B1SEMCAD.mat','Absolute_Modulus_of_B1__0s','Absolute_Modulus_of_B1p_0s','Phase_of_B1__0s','Phase_of_B1p_0s','Xmid','Ymid','Zmid')
%%


