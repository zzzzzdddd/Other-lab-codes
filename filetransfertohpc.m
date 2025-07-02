%%patients


%% volunteers
subjID=["XXX];
cd /Users/irene/Desktop
for p=subjID
    line = strcat("mput ", "/Volumes/eegrvw/Imaging/Multimodal/MRF/Recon_MRF_3T/Normal/", p, "/MRF_VBM/* /mnt/beegfs/dingz/all_patients/", p, "/MRF_VBM");
    writelines([line], "sftp_hpc_commands.rtf", WriteMode="append");
end