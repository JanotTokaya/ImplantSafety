function visualizeDistr2(distribution,x_axis,y_axis,useColorbar,scale,backGroundColor,mask)
%Visualize field distribution without interpolation (plotting unequal rectangles for each voxel).
%
% distribution:     2-D field distribution to image. If complex, absolute value
%                   will be depicted. 
% x_axis:           First axis corresponding to the distribution
% y_axis:           Second axis corresponding to the distribution
% useColorbar       Draw colorbar along with image. Default = 'true'
% scale:            Scale for imagesc commands. Default uses default of imagesc
% backGroundColor:  Color given to pixels with zero values (presumed
%                   background). Default: white. If empty matrix is passed
%                   [], pixel color with zero values are not changed.
% mask:             2-D distribution of logicals with the same dimensions
%                   as 'distribution'. Indicates mask pattern. Default mask
%                   pattern is to mask pixels with values equal to zero.
%                   To not use a mask in the visualization, pass the value
%                   0 for this parameter.
%
%imageData = visualizeDistr(distribution)
%imageData = visualizeDistr(distribution,x_axis,y_axis,useColorBar)
%imageData = visualizeDistr(distribution,x_axis,y_axis,useColorBar,scale)
%imageData = visualizeDistr(distribution,x_axis,y_axis,useColorBar,scale,backGroundColor)
%imageData = visualizeDistr(distribution,x_axis,y_axis,useColorBar,scale,backGroundColor,mask)
%

useMask=1;

if (nargin == 3)
    scale=[];
    backGroundColor = [1.0 1.0 1.0];
    mask=[];
    useColorbar=1;
elseif (nargin == 4)
    scale=[];
    backGroundColor = [1.0 1.0 1.0];
    mask=[];
elseif (nargin == 5)
    backGroundColor = [1.0 1.0 1.0];
    mask=[];
elseif (nargin == 6)
    mask=[];
elseif (nargin == 7)
    if ((length(mask) == 1) && (mask == 0))
        useMask=0;
    end
elseif (nargin > 7)
    error('Wrong number of input arguments');
end

if (isempty(backGroundColor))
    useMask=0;
    backGroundColor = [1.0 1.0 1.0];
end

if (isempty(scale))
    scale = [0 max(max(distribution))];
end

if (isempty(mask) || ((length(mask) == 1) && (mask == 0)))
    if (useMask==0) 
        mask = ones(size(distribution));
    else
        mask = double(abs(distribution>0));
    end
end
    

map = colormap;
cla

for (k=1:size(distribution,1)-1)
    for (m=1:size(distribution,2)-1)
        if (~useMask || mask(k,m))
            pos = [x_axis(k) y_axis(m) x_axis(k+1)-x_axis(k) y_axis(m+1)-y_axis(m)];
            colIndex = 1+63*max(min((distribution(k,m) - scale(1))/(scale(2)-scale(1)),1),0);
           
            col = squeeze(map(round(colIndex),:));
            rectangle('Position',pos,'FaceColor',col,'EdgeColor','None');
        end
    end
end
set(gca,'YLim',[min(y_axis) max(y_axis)])
set(gca,'XLim',[min(x_axis) max(x_axis)])
set(gca,'Color',backGroundColor)
axis image

if (useColorbar)
    colorbar
    caxis(scale)
else
    colorbar delete
end
    
