function [keypoints] = sift_extra( image )
%SIFT_EXTRA Summary of this function goes here
%   Detailed explanation goes here
    %define the scale space
    image1 = image;
    retScaleSpace = scaleSpace(image1,4,3);%4 oct ��ͼ�������
    octaveStack = retScaleSpace{1}; %ȡ��oct
    accumSigmas = retScaleSpace{2}; %ȡ��accum
    octaveDOGStack = calculateDog(octaveStack);%��ÿ��ͼ�����������DOG
    keypoints = calculateKeypoints(octaveDOGStack, image1);%����DOG������������С�ĵ�
    orientationDef = defineOrientation(keypoints, octaveDOGStack, ...
               octaveStack, image1, accumSigmas);    
    descriptor = localDescriptor_v3(orientationDef, keypoints, ...
               accumSigmas, size(image1,1)*2, size(image1,2)*2); 
    %to plot the descriptor, uncomment and comment the rest of the code below
    %plotDescriptor(descriptor, image1, orientationDef, keypoints);
    %siftkeypoints =keypoints;
end

