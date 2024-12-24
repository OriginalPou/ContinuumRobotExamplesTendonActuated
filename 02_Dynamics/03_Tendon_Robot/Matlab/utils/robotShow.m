function robotShow()
    global Y  ;
    centerline=Y(1:3,:);
    tendonlines = zeros(3*num_tendons, size(centerline,2));
    for i = 1 : size(centerline,2)
        p_show = Y(1:3,i);
        R_show = reshape(Y(4:12,i),3,3);
        for j = 1 : num_tendons
            tendonlines(3*j-2 : 3*j, i) = p_show + R_show*r{j};
        end
    end
    disks = zeros(3,4*num_disks);
    for i = 1 : num_disks
        j = round(size(centerline,2) * i / num_disks);
        p_show = Y(1:3,j);
        R_show = reshape(Y(4:12,j),3,3);
        disks(1:3, 4*i-3:4*i-1) = R_show;
        disks(1:3, 4*i) = p_show;
    end
    Visualize(centerline,tendonlines,num_tendons, disks,num_disks);
end