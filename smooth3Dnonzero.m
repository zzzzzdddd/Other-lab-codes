function [smoothed_img1, smoothed_img2] = smooth3Dnonzero(image1, image2, kernel_size, csfi)
    t1max = 1759; 
    t2max = 106;
    t1min = 743;
    t2min = 29;
    CSFthresh = 0.05;
%     kernel = ones(kernel_size) / prod(kernel_size);
    smoothed_img1 = zeros(size(image1));
    smoothed_img2 = zeros(size(image2));
    for x = 1:size(image1, 1)
        for y = 1:size(image1, 2)
            for z = 1:size(image1, 3)
                if image1(x, y, z) ~= 0 & (image1(x, y, z) > t1min) & (image1(x, y, z) < t1max)
                    x_range = max(1, x-floor(kernel_size/2)):min(size(image1,1), x+floor(kernel_size/2));
                    y_range = max(1, y-floor(kernel_size/2)):min(size(image1,2), y+floor(kernel_size/2));
                    z_range = max(1, z-floor(kernel_size/2)):min(size(image1,3), z+floor(kernel_size/2));
                    
                    neighborhood1 = image1(x_range, y_range, z_range);
                    neighborhood2 = image2(x_range, y_range, z_range);
                    csfhood = csfi(x_range, y_range, z_range);

                    non_zero_elements1 = neighborhood1((neighborhood1 < t1max)&(neighborhood1 > t1min)&(csfhood<CSFthresh));
                    non_zero_elements2 = neighborhood2((neighborhood2 < t2max)&(neighborhood2 > t2min)&(csfhood<CSFthresh));
                    
%                     non_zero_elements1 = neighborhood1((neighborhood1 < t1max)&(neighborhood1 > t1min)&(csfhood<CSFthresh));
%                     non_zero_elements2 = neighborhood2((neighborhood2 < t2max)&(neighborhood2 > t2min)&(csfhood<CSFthresh));

                    if ~isempty(non_zero_elements1)
                        smoothed_img1(x, y, z) = mean(non_zero_elements1);
%                     else
%                         smoothed_img1(x, y, z) = image1(x, y, z);
                    end
                    if ~isempty(non_zero_elements2)
                        smoothed_img2(x, y, z) = mean(non_zero_elements2);
%                     else
%                         smoothed_img2(x, y, z) = image2(x, y, z);
                    end                
                end
            end
        end
    end
end