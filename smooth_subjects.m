% List of open inputs
nrun = X; % enter the number of runs here
jobfile = {'T:\Imaging\Multimodal\MRF\Peter\smooth_subjects_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});
