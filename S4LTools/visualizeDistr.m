function imageData = visualizeDistr(distribution,x_axis,y_axis,resolution,scale,backGroundColor)
%imageData = visualizeDistr(distribution)
%imageData = visualizeDistr(distribution,x_axis,y_axis,resolution)
%imageData = visualizeDistr(distribution,x_axis,y_axis,resolution,scale)
%imageData = visualizeDistr(distribution,x_axis,y_axis,resolution,scale,backGroundColor)
%
% distribution:     2-D field distribution to image. If complex, absolute value
%                   will be depicted. If no more arguments are given,
%                   distribution is not resampled.
% x_axis:           First axis corresponding to the distribution
% y_axis:           Second axis corresponding to the distribution
% resolution:       Resolution to interpolate to e.g. [0.002 0.002]
% scale:            Scale for imagesc commands. Default uses default of imagesc
% backGroundColor:  Color given to pixels with zero values (presumed
%                   background). Default: white. If empty matrix is passed
%                   [], pixel color with zero values are not changed.
%
%Interpolate and visualize field distribution.

if (nargin == 1)
    x_axis=[];
    y_axis=[];
    resolution=[];
    scale=[];
    backGroundColor = [1.0 1.0 1.0];
elseif (nargin == 4)
    scale=[];
    backGroundColor = [1.0 1.0 1.0];
elseif (nargin == 5)
    backGroundColor = [1.0 1.0 1.0];
elseif (nargin ~= 6)
    error('Wrong number of input arguments');
end

if (~isempty(x_axis))
    
    xmin = min([x_axis]);
    ymin = min([y_axis]);
    xmax = max([x_axis]);
    ymax = max([y_axis]);

    xnew= xmin:resolution(1):xmax;
    ynew= ymin:resolution(2):ymax;

    [Xoldi Yoldi] = meshgrid(y_axis,x_axis);
    [Xnew Ynew] = meshgrid(ynew,xnew);
    imageData = interp2(Xoldi,Yoldi,distribution,Xnew,Ynew);
else
    imageData = distribution;
end

if (isempty(scale))
    imagesc(abs(imageData)); axis equal; set(gca,'XTick',[]); set(gca,'YTick',[])
else
    imagesc(abs(imageData),scale); axis equal; set(gca,'XTick',[]); set(gca,'YTick',[])
end
if (~isempty(backGroundColor))
    mask=double(abs(imageData)>0);
    alpha(mask); set(gca,'color',backGroundColor)
end


