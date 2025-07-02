function [mT1, sdT1, mT2, sdT2, mGMT1, sdGMT1, mWMT1, sdWMT1, mGMT2, sdGMT2, mWMT2, sdWMT2] = MRFclusterwisezscore(region)
    subjID=[XXXX];
         
    MRF_path='/home/dingz/beegfs/all_patients';
    mT1s = [];
    mT2s = [];
    mGMT1s = [];
    mWMT1s = [];
    mGMT2s = [];
    mWMT2s = [];    


    t1max = 1759; 
    t2max = 106;
    t1min = 743;
    t2min = 29;
    for p = subjID
        path = strcat(MRF_path,'/',p,'/MRF_VBM');
        disp(path)
        cd(path);
        f1 = 'MNI_T1.nii';
        f2 = 'MNI_T2.nii';
        f3 = 'MNI_GM_fn.nii';
        f4 = 'MNI_WM_fn.nii';
        a = single(load_untouch_nii(f1).img);
        b = single(load_untouch_nii(f2).img);
        c = single(load_untouch_nii(f3).img);
        d = single(load_untouch_nii(f4).img);
        ac = a;
        bc = b;
        mask = c+d;
        mask(mask>0.95) = 1;
        mask(mask<0.95) = 0;
        c = c.*mask;
        d = d.*mask;
        c(c<0.6) = 0;
        c(c>=0.6) = 1;
        d(d<0.6) = 0;
        d(d>=0.6) = 1;

        good = ((a<t1max)&(a>t1min)&(b<t2max)&(b>t2min));
        a = ac(good).*region(good);
        b = bc(good).*region(good);
        
        e = ac(good).*region(good).*c(good);
        f = ac(good).*region(good).*d(good);
        g = bc(good).*region(good).*c(good);
        h = bc(good).*region(good).*d(good);
        
        mT1s = horzcat(mT1s, mean(a(a>0)));
        mT2s = horzcat(mT2s, mean(b(b>0)));
        mGMT1s = horzcat(mGMT1s, rmmissing(mean(e(e>0))));
        mWMT1s = horzcat(mWMT1s, rmmissing(mean(f(f>0))));
        mGMT2s = horzcat(mGMT2s, rmmissing(mean(g(g>0))));
        mWMT2s = horzcat(mWMT2s, rmmissing(mean(h(h>0))));
    end
    mT1 = mean(mT1s);
    mT2 = mean(mT2s);
    sdT1 = std(mT1s);
    sdT2 = std(mT2s);

    mGMT1 = mean(mGMT1s);
    mGMT2 = mean(mGMT2s);
    sdGMT1 = std(mGMT1s);
    sdGMT2 = std(mGMT2s);

    mWMT1 = mean(mWMT1s);
    mWMT2 = mean(mWMT2s);
    sdWMT1 = std(mWMT1s);
    sdWMT2 = std(mWMT2s);
end