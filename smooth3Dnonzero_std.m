function [smoothed_img1, smoothed_img2] = smooth3Dnonzero_std(image1, image2, sd1, sd2, kernel_size, csfi)
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
                if image1(x, y, z) > 0
                    x_range = max(1, x-floor(kernel_size/2)):min(size(image1,1), x+floor(kernel_size/2));
                    y_range = max(1, y-floor(kernel_size/2)):min(size(image1,2), y+floor(kernel_size/2));
                    z_range = max(1, z-floor(kernel_size/2)):min(size(image1,3), z+floor(kernel_size/2));
                    
                    neighborhood1 = sd1(x_range, y_range, z_range);
                    neighborhood2 = sd2(x_range, y_range, z_range);
                    csfhood = csfi(x_range, y_range, z_range);
%                     non_zero_elements1 = sd1(neighborhood1 > 0);
%                     non_zero_elements2 = sd2(neighborhood2 > 0);
                    
%                     non_zero_elements1 = neighborhood1((neighborhood1 < t1max)&(neighborhood1 > t1min)&(csfhood<CSFthresh));
%                     non_zero_elements2 = neighborhood2((neighborhood2 < t2max)&(neighborhood2 > t2min)&(csfhood<CSFthresh));

                    if ~isempty(neighborhood1((neighborhood1>0)&(csfhood<CSFthresh)))
                        smoothed_img1(x, y, z) =  mean(neighborhood1((neighborhood1>0)&(csfhood<CSFthresh)));
                       
%                     else
%                         smoothed_img1(x, y, z) = image1(x, y, z);
                    end
                    if ~isempty(neighborhood2((neighborhood2>0)&(csfhood<CSFthresh)))
                        smoothed_img2(x, y, z) =  mean(neighborhood2((neighborhood2>0)&(csfhood<CSFthresh)));
%                     else
%                         smoothed_img2(x, y, z) = image2(x, y, z);
                    end                
                end
            end
        end
    end
end