image1 = imread('rgb1.png');


%define the scale space
retScaleSpace = scaleSpace(image1,4,3);%4 oct 的图像金字塔
octaveStack = retScaleSpace{1}; %取出oct
accumSigmas = retScaleSpace{2}; %取出accum
octaveDOGStack = calculateDog(octaveStack);%对每层图像金子塔计算DOG
keypoints = calculateKeypoints(octaveDOGStack, image1);%根据DOG计算最大或者最小的点
orientationDef = defineOrientation(keypoints, octaveDOGStack, ...
           octaveStack, image1, accumSigmas);    
descriptor = localDescriptor_v3(orientationDef, keypoints, ...
           accumSigmas, size(image1,1)*2, size(image1,2)*2); 
%to plot the descriptor, uncomment and comment the rest of the code below
plotDescriptor(descriptor, image1, orientationDef, keypoints);