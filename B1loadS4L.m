function [B1p,B1m, grid] = B1loadS4L(path)

data = load(path);
B1p = reshape(data.Snapshot0(:,1), [length(data.Axis0)-1 length(data.Axis1)-1 length(data.Axis2)-1]);
B1m = reshape(data.Snapshot0(:,2), [length(data.Axis0)-1 length(data.Axis1)-1 length(data.Axis2)-1]);
grid.xaxis=data.Axis0; 
grid.yaxis=data.Axis1; 
grid.zaxis=data.Axis2; 

end 