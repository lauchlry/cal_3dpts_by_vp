function aaa = disp_omap(omap, img, OPACITY, keypoints,vp,f,u0,v0)


if ~exist('OPACITY','var')
OPACITY = 0.5;
end

ooo = imresize(double(omap), [size(img,1) size(img,2)], 'nearest');
aaa = uint8(OPACITY*double(img)) + uint8((1-OPACITY)*ooo*256);

    
%add show sift feature 
   keypointDescriptor = keypoints{1};
    
%    orientations = orientationDef{1}; 
    %for each keypoints in octave plots the dots
    xx = [];
    yy = [];
    zz = [];
    xx2 = [];
    yy2 = [];
    zz2 = [];
    xx3 = [];
    yy3 = [];
    zz3 = [];
    for octave = 1:size(keypointDescriptor, 1)
        %for each keypoints layer 
        for kptLayer = 1:size(keypointDescriptor,2)
            
            [rowKpt colKpt] = find(keypointDescriptor{octave,kptLayer} == 1); 
            
             for keypoint = 1:size([rowKpt colKpt],1)
                 %plots the dot 
                 %if this point is in omap[1][2][3]
                 %then:
                    [relRow relCol] = calRelPoint(rowKpt(keypoint), colKpt(keypoint), octave);
                   
                    %%if points on the floor
                    %if ooo(relRow,relCol,1) ~= 0 && ...
                    %    ooo(relRow,relCol, 2) == 0 && ...
                    %    ooo(relRow,relCol, 3) == 0
                    %    sp1 = compute_floor_depth(relRow, relCol, u0, v0, f, vp);
                    %    xx = [xx,sp1(1)]
                    %    yy = [yy,sp1(2)]
                    %    zz = [zz,sp1(3)]
                    %end
                    
                    %%if points on the wall
                    %if ooo(relRow,relCol,1) == 0 && ...
                    %    ooo(relRow,relCol, 2) ~= 0 && ...
                    %    ooo(relRow,relCol, 3) == 0
                        %%step1---line_base:points and vp
                    %    line_base.point1 = [relRow,relCol,1];
                    %    line_base.point2 = [vp{1}(1),vp{1}(2),1];
                       %%step2---sample on line_base
                    %    ls_base = sample_on_line(line_base);
                        %%step3---sample untile first point on the floor occur
                    %    p1 = first_points_onfloor(ls_base,ooo);
                    %    if p1 == [0,0,0]
                    %        break;
                    %    end
                        %%step4---compute point on the wall
                    %    sp2 = compute_wall_depth(relRow, relCol, u0, v0, f, vp, p1);
                    %    xx2 = [xx2,sp2(1)];
                    %    yy2 = [yy2,sp2(2)];
                    %    zz2 = [zz2,sp2(3)];
                    %end
                    %aaa = plotDot(aaa, rowKpt(keypoint), colKpt(keypoint), octave, ooo);
                 %end
                 
                 
             end 
        end 
    end
   
   %%restructe 3d of floor
   [xx, yy, zz] = sampleImage_onfloor(aaa, u0, v0, f, vp, ooo);
   
   m = size(xx,2);
   
   %%restructe 3d of wall
   [xx2, yy2 , zz2] = sampleimage_onwall(aaa, u0, v0, f, vp, ooo);
   
   n = size(xx2,2);
   
   %%
   [xx3, yy3 , zz3] = sampleimage_onwall2(aaa, u0, v0, f, vp, ooo);
   r = size(xx3,2);
   
   INDEX = [ones(m,1); ones(n,1)+1; ones(r,1)+2];
   color_1 = [1 0 0];
   color_2 = [0 0 1];
   color_3 = [0 1 0];
   
   cmap = [color_1; color_2; color_3];
   
   INDEX_color = cmap(INDEX,:);
   xx = [xx, xx2]
   yy = [yy, yy2]
   zz = [zz, zz2]
   xx = [xx, xx3]
   yy = [yy, yy3]
   zz = [zz, zz3]
   
   %color
   
   figure;
   scatter3(xx,yy,zz,50, INDEX_color, 'filled');
   

   %figure;
   hold on;
   %if nargout==0
        %figure; 
        %imshow(aaa)
   %end
    %elseif nargout==1
end

function [xx, yy, zz] = sampleImage_onfloor(image, u0, v0, f, vp, ooo)
    xx = [];
    yy = [];
    zz = [];
   
    [m,n,z] = size(image);
    for i = 2:10:m-1
        for j = 2:13:n-1
             if ooo(i,j, 1) ~= 0 && ...
                ooo(i,j, 2) == 0 && ...
                ooo(i,j, 3) == 0
                    image(i,j,1) = 255; 
                    image(i,j,2) = 255; 
                    image(i,j,3) = 0; 
                    image(i-1,j-1,1) = 255; 
                    image(i-1,j-1,2) = 255; 
                    image(i-1,j-1,3) = 0; 
                    image(i+1,j+1,1) = 255; 
                    image(i+1,j+1,2) = 255; 
                    image(i+1,j+1,3) = 0; 
                    image(i-1,j+1,1) = 255; 
                    image(i-1,j+1,2) = 255; 
                    image(i-1,j+1,3) = 0; 
                    image(i+1,j-1,1) = 255; 
                    image(i+1,j-1,2) = 255; 
                    image(i+1,j-1,3) = 0;       
                    sp1 = compute_floor_depth(j, i, u0, v0, f, vp);
                    xx = [xx,sp1(1)];
                    yy = [yy,sp1(2)];
                    zz = [zz,sp1(3)];
              end

        end
    end
    %figure;
    %scatter3(xx,yy,zz);
    figure; 
    imshow(image);
end

function [xx2, yy2, zz2] = sampleimage_onwall(image, u0, v0, f, vp, ooo)
    xx2 = [];
    yy2 = [];
    zz2 = [];
   
    [m,n,z] = size(image);
    for i = 2:10:m-1
        for j = 2:13:n-1
             if ooo(i,j, 1) == 0 && ...
                ooo(i,j, 2) == 0 && ...
                ooo(i,j, 3) ~= 0
                    image(i,j,1) = 255; 
                    image(i,j,2) = 255; 
                    image(i,j,3) = 0; 
                    image(i-1,j-1,1) = 255; 
                    image(i-1,j-1,2) = 255; 
                    image(i-1,j-1,3) = 0; 
                    image(i+1,j+1,1) = 255; 
                    image(i+1,j+1,2) = 255; 
                    image(i+1,j+1,3) = 0; 
                    image(i-1,j+1,1) = 255; 
                    image(i-1,j+1,2) = 255; 
                    image(i-1,j+1,3) = 0; 
                    image(i+1,j-1,1) = 255; 
                    image(i+1,j-1,2) = 255; 
                    image(i+1,j-1,3) = 0;
                    
                   %%step1---line_base:points and vp
                    line_base.point1 = [j,i,1];
                    line_base.point2 = [vp{1}(1),vp{1}(2),1];
                   %%step2---sample on line_base
                    ls_base = sample_on_line(line_base);
                    %%step3---sample untile first point on the floor occur
                    p1 = first_points_onfloor(ls_base,ooo);
                    if p1 == [0,0,0]
                        continue;
                    end
                    %%step4---compute point on the wall
                    sp2 = compute_wall_depth(j, i, u0, v0, f, vp, p1);
                    xx2 = [xx2,sp2(1)];
                    yy2 = [yy2,sp2(2)];
                    zz2 = [zz2,sp2(3)];
             end
        end
    end
    %figure;
    %scatter3(xx2,yy2,zz2);
    figure; 
    imshow(image);
end

function [xx3, yy3, zz3] = sampleimage_onwall2(image, u0, v0, f, vp, ooo)
    xx3 = [];
    yy3 = [];
    zz3 = [];
   
    [m,n,z] = size(image);
    for i = 2:10:m-1
        for j = 2:13:n-1
             if ooo(i,j, 1) == 0 && ...
                ooo(i,j, 2) ~= 0 && ...
                ooo(i,j, 3) == 0
                    image(i,j,1) = 255; 
                    image(i,j,2) = 255; 
                    image(i,j,3) = 0; 
                    image(i-1,j-1,1) = 255; 
                    image(i-1,j-1,2) = 255; 
                    image(i-1,j-1,3) = 0; 
                    image(i+1,j+1,1) = 255; 
                    image(i+1,j+1,2) = 255; 
                    image(i+1,j+1,3) = 0; 
                    image(i-1,j+1,1) = 255; 
                    image(i-1,j+1,2) = 255; 
                    image(i-1,j+1,3) = 0; 
                    image(i+1,j-1,1) = 255; 
                    image(i+1,j-1,2) = 255; 
                    image(i+1,j-1,3) = 0;
                    
                   %%step1---line_base:points and vp
                    line_base.point1 = [j,i,1];
                    line_base.point2 = [vp{1}(1),vp{1}(2),1];
                   %%step2---sample on line_base
                    ls_base = sample_on_line(line_base);
                    %%step3---sample untile first point on the floor occur
                    p1 = first_points_onfloor(ls_base,ooo);
                    if p1 == [0,0,0]
                        continue;
                    end
                    %%step4---compute point on the wall
                    sp2 = compute_wall_depth(j, i, u0, v0, f, vp, p1);
                    xx3 = [xx3,sp2(1)];
                    yy3 = [yy3,sp2(2)];
                    zz3 = [zz3,sp2(3)];
             end
        end
    end
    %figure;
    %scatter3(xx2,yy2,zz2);
    figure; 
    imshow(image);
end
%liu
function ls_base = sample_on_line(line)
    sample_rate = 5;
    n_sample = ceil(norm(line.point1-line.point2)/sample_rate);
    ls_base.n_sample = n_sample;
    ls_base.sample = [ ...
		linspace(line.point1(1), line.point2(1), n_sample)' ...
		linspace(line.point1(2), line.point2(2), n_sample)' ];
end

%liu
function fpf = first_points_onfloor(line,ooo)
    fpf = [0,0,0];
    for i = 1:size(line.sample)
        x = int32(line.sample(i,1));
        y = int32(line.sample(i,2));
        %fprintf('%d %d ', x, y);
        %fprintf('%d %d\n', size(ooo,2), size(ooo,1));
        if x < size(ooo,2) && ...
           y < size(ooo,1)
        
            if ooo(y,x,1) ~= 0 && ...
              ooo(y,x,2) == 0 && ...
              ooo(y,x,3) == 0
                fpf = [double(x),double(y),double(1)];
                break;
            end
        end
    end
end

%show points one some omap  
function plotDot = plotDot(image, row, col, octave, ooo)
        
        relRow = row; 
        relCol = col; 
        
        if(octave==1)
            relRow = round(row/2); 
            relCol = round(col/2); 
        end

        if(octave>2)                        
            relRow = row * (2^(octave-2)); 
            relCol = col * (2^(octave-2)); 
        end
        if(relRow==1)
            relRow = 2; 
        end
        if(relCol==1)
            relCol = 2; 
        end
        if ooo(relRow,relCol, 1) ~= 0 || ...
           ooo(relRow,relCol, 2) ~= 0 || ...
           ooo(relRow,relCol, 3) ~= 0
            image(relRow,relCol,1) = 255; 
            image(relRow,relCol,2) = 255; 
            image(relRow,relCol,3) = 0; 
            image(relRow-1,relCol-1,1) = 255; 
            image(relRow-1,relCol-1,2) = 255; 
            image(relRow-1,relCol-1,3) = 0; 
            image(relRow+1,relCol+1,1) = 255; 
            image(relRow+1,relCol+1,2) = 255; 
            image(relRow+1,relCol+1,3) = 0; 
            image(relRow-1,relCol+1,1) = 255; 
            image(relRow-1,relCol+1,2) = 255; 
            image(relRow-1,relCol+1,3) = 0; 
            image(relRow+1,relCol-1,1) = 255; 
            image(relRow+1,relCol-1,2) = 255; 
            image(relRow+1,relCol-1,3) = 0;         

        end
        plotDot = image;
        
end

%由于sift对图像金字塔进行操作，返回真实的row和col
function [relRow relCol] = calRelPoint(row, col, octave)
        relRow = row; 
        relCol = col; 
        
        if(octave==1)
            relRow = round(row/2); 
            relCol = round(col/2); 
        end

        if(octave>2)                        
            relRow = row * (2^(octave-2)); 
            relCol = col * (2^(octave-2)); 
        end
        if(relRow==1)
            relRow = 2; 
        end
        if(relCol==1)
            relCol = 2; 
        end
end

%liu
%选一个地面方向的消影点，计算地面上的特征点的空间坐标
function sp = compute_floor_depth(relRow, relCol, u0, v0, f, vp)
    p = [relRow,relCol,1]';
    K = [f,0,u0;0,f,v0;0,0,1];
    inv_K = inv(K);
    vp_pt = [vp{1}(1),vp{1}(2),1]';
    normm = norm((inv_K*vp_pt),2);
    sv = (inv_K*vp_pt)/norm((inv_K*vp_pt),2);
    sp = (inv_K*p)/(sv'*inv_K*p);
end

%liu
%compute the depth of point on the wall using points on the floor
function sp = compute_wall_depth(relRow, relCol, u0, v0, f, vp, point_floor)
    p = [relRow,relCol,1];
    pp = [double(relRow),double(relCol)];
    point_floorr = [point_floor(1),point_floor(2)]; 
    h = norm(pp-point_floorr)/100;
    %h=1;
    K = [f,0,u0;0,f,v0;0,0,1];
    vp_pt = [vp{1}(1),vp{1}(2),1]';
    sv = inv(K)*vp_pt/norm((inv(K)*vp_pt),2);
    sp1 = inv(K)*point_floor'/(sv'*inv(K)*point_floor');
    sp = sp1+h*sv;
end