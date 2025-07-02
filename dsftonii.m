mask = tess_mrimask(size(sMRI.Cube), sSurface.tess2mri_interp)

% You need first to right-click on the surface > MRI registration
% > Check MRI/surface registration, so that the field tess2mri_interp is computed for this file.

sMask = sMRI;
sMask.Cube = mask;
out_mri_nii(sMask, 'thicknessmap.nii');