Process amsterdam24

# Copy BIDS files 
rsync -av /scratch/mszinte/data/raw_data/prfGaze/dataset_description.json /scratch/mszinte/data/amsterdam24/
rsync -av /scratch/mszinte/data/raw_data/prfGaze/participants.json /scratch/mszinte/data/amsterdam24/
rsync -av /scratch/mszinte/data/raw_data/prfGaze/participants.tsv /scratch/mszinte/data/amsterdam24/

# Task event files
rsync -av /scratch/mszinte/data/raw_data/prfGaze/task-*_events.json /scratch/mszinte/data/amsterdam24/

# convert to nifti per session and subject (compress & make json sidecar)
dcm2niix -z y -b y /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfaniso
dcm2niix -z y -b y /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf
dcm2niix -z y -b y /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfoccl
dcm2niix -z y -b y /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2
dcm2niix -z y -b y /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfstrab
dcm2niix -z y -b y /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3

# rename and restructure in new datarepo (amsterdam24) BIDS 
# ignore ph files 
# sub-XXX/
#    ses-YY/
#        func/
#        fmap/
mkdir  /scratch/mszinte/data/amsterdam24/sub-001 
mkdir  /scratch/mszinte/data/amsterdam24/sub-001/ses-01
mkdir  /scratch/mszinte/data/amsterdam24/sub-001/ses-01/func 
mkdir  /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap
#same for all other sessions ! 

# session mapping sub-001 
# ses-rsprf      = ses-01 (22.07 morning)
# ses-prfoccl    = ses-02 (22.07 midday)
# ses-rs2        = ses-03 (22.07 evening)
# ses-prfaniso   = ses-04 (24.07 morning)
# ses-res3       = ses-05 (25.07 midday)
# ses-prfstrab   = ses-06 (26.07 midday)

# copy func files into func per session 
# sub-001 
# ses-01
# task RSR 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_func-bold_task-RSR_run-01_20240722083343_13.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/func/ses-rsprf_func-bold_task-RSR_run-01_20240722083343_13.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_func-bold_task-RSR_run-01_20240722083343_13.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/func/ses-rsprf_func-bold_task-RSR_run-01_20240722083343_13.json
# task RSL 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_func-bold_task-RSL_run-01_20240722083343_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/func/ses-rsprf_func-bold_task-RSL_run-01_20240722083343_10.nii.gz 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_func-bold_task-RSL_run-01_20240722083343_10.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/func/ses-rsprf_func-bold_task-RSL_run-01_20240722083343_10.json
# task RSC 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_func-bold_task-RSC_run-01_20240722083343_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/func/ses-rsprf_func-bold_task-RSC_run-01_20240722083343_9.nii.gz 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_func-bold_task-RSC_run-01_20240722083343_9.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/func/ses-rsprf_func-bold_task-RSC_run-01_20240722083343_9.json
# task pRF
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_func-bold_task-pRF_run-01_20240722083343_16.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/func/ses-rsprf_func-bold_task-pRF_run-01_20240722083343_16.nii.gz 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_func-bold_task-pRF_run-01_20240722083343_16.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/func/ses-rsprf_func-bold_task-pRF_run-01_20240722083343_16.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_func-bold_task-pRF_run-02_20240722083343_17.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/func/ses-rsprf_func-bold_task-pRF_run-02_20240722083343_17.nii.gz 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_func-bold_task-pRF_run-02_20240722083343_17.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/func/ses-rsprf_func-bold_task-pRF_run-02_20240722083343_17.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_func-bold_task-pRF_run-03_20240722083343_18.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/func/ses-rsprf_func-bold_task-pRF_run-03_20240722083343_18.nii.gz 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_func-bold_task-pRF_run-03_20240722083343_18.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/func/ses-rsprf_func-bold_task-pRF_run-03_20240722083343_18.json

# ses-02 
# task prfoccl
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfoccl/ses-prfoccl_func-bold_task-pRFoccl_run-01_20240722123602_6.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-02/func/ses-prfoccl_func-bold_task-pRFoccl_run-01_20240722123602_6.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfoccl/ses-prfoccl_func-bold_task-pRFoccl_run-01_20240722123602_6.json /scratch/mszinte/data/amsterdam24/sub-001/ses-02/func/ses-prfoccl_func-bold_task-pRFoccl_run-01_20240722123602_6.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfoccl/ses-prfoccl_func-bold_task-pRFoccl_run-02_20240722123602_7.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-02/func/ses-prfoccl_func-bold_task-pRFoccl_run-02_20240722123602_7.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfoccl/ses-prfoccl_func-bold_task-pRFoccl_run-02_20240722123602_7.json /scratch/mszinte/data/amsterdam24/sub-001/ses-02/func/ses-prfoccl_func-bold_task-pRFoccl_run-02_20240722123602_7.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfoccl/ses-prfoccl_func-bold_task-pRFoccl_run-03_20240722123602_8.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-02/func/ses-prfoccl_func-bold_task-pRFoccl_run-03_20240722123602_8.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfoccl/ses-prfoccl_func-bold_task-pRFoccl_run-03_20240722123602_8.json /scratch/mszinte/data/amsterdam24/sub-001/ses-02/func/ses-prfoccl_func-bold_task-pRFoccl_run-03_20240722123602_8.json

# ses-03
# task RSC
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_func-bold_task-RSC_run-01_20240722163325_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-03/func/ses-rs2_func-bold_task-RSC_run-01_20240722163325_9.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_func-bold_task-RSC_run-01_20240722163325_9.json /scratch/mszinte/data/amsterdam24/sub-001/ses-03/func/ses-rs2_func-bold_task-RSC_run-01_20240722163325_9.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_func-bold_task-RSC_run-02_20240722163325_12.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-03/func/ses-rs2_func-bold_task-RSC_run-02_20240722163325_12.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_func-bold_task-RSC_run-02_20240722163325_12.json /scratch/mszinte/data/amsterdam24/sub-001/ses-03/func/ses-rs2_func-bold_task-RSC_run-02_20240722163325_12.json
# task RSL
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_func-bold_task-RSL_run-01_20240722163325_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-03/func/ses-rs2_func-bold_task-RSL_run-01_20240722163325_10.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_func-bold_task-RSL_run-01_20240722163325_10.json /scratch/mszinte/data/amsterdam24/sub-001/ses-03/func/ses-rs2_func-bold_task-RSL_run-01_20240722163325_10.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_func-bold_task-RSL_run-02_20240722163325_13.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-03/func/ses-rs2_func-bold_task-RSL_run-02_20240722163325_13.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_func-bold_task-RSL_run-02_20240722163325_13.json /scratch/mszinte/data/amsterdam24/sub-001/ses-03/func/ses-rs2_func-bold_task-RSL_run-02_20240722163325_13.json
# task RSR 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_func-bold_task-RSR_run-01_20240722163325_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-03/func/ses-rs2_func-bold_task-RSR_run-01_20240722163325_11.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_func-bold_task-RSR_run-01_20240722163325_11.json /scratch/mszinte/data/amsterdam24/sub-001/ses-03/func/ses-rs2_func-bold_task-RSR_run-01_20240722163325_11.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_func-bold_task-RSR_run-02_20240722163325_14.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-03/func/ses-rs2_func-bold_task-RSR_run-02_20240722163325_14.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_func-bold_task-RSR_run-02_20240722163325_14.json /scratch/mszinte/data/amsterdam24/sub-001/ses-03/func/ses-rs2_func-bold_task-RSR_run-02_20240722163325_14.json

# ses-04 
# task prfaniso
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-01_20240724102635_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-04/func/ses-prfaniso_func-bold_task-pRFaniso_run-01_20240724102635_9.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-01_20240724102635_9.json /scratch/mszinte/data/amsterdam24/sub-001/ses-04/func/ses-prfaniso_func-bold_task-pRFaniso_run-01_20240724102635_9.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-02_20240724102635_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-04/func/ses-prfaniso_func-bold_task-pRFaniso_run-02_20240724102635_10.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-02_20240724102635_10.json /scratch/mszinte/data/amsterdam24/sub-001/ses-04/func/sses-prfaniso_func-bold_task-pRFaniso_run-02_20240724102635_10.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-03_20240724102635_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-04/func/ses-prfaniso_func-bold_task-pRFaniso_run-03_20240724102635_11.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-03_20240724102635_11.json /scratch/mszinte/data/amsterdam24/sub-001/ses-04/func/ses-prfaniso_func-bold_task-pRFaniso_run-03_20240724102635_11.json

# ses-05 
# task CAL 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_func-bold_task-CAL_run-01_20240725120124_12.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-05/func/ses-rs3_func-bold_task-CAL_run-01_20240725120124_12.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_func-bold_task-CAL_run-01_20240725120124_12.json /scratch/mszinte/data/amsterdam24/sub-001/ses-05/func/ses-rs3_func-bold_task-CAL_run-01_20240725120124_12.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_func-bold_task-CAL_run-03_20240725120124_13.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-05/func/ses-rs3_func-bold_task-CAL_run-02_20240725120124_13.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_func-bold_task-CAL_run-03_20240725120124_13.json /scratch/mszinte/data/amsterdam24/sub-001/ses-05/func/ses-rs3_func-bold_task-CAL_run-02_20240725120124_13.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_func-bold_task-CAL_run-03_20240725120124_14.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-05/func/ses-rs3_func-bold_task-CAL_run-03_20240725120124_14.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_func-bold_task-CAL_run-03_20240725120124_14.json /scratch/mszinte/data/amsterdam24/sub-001/ses-05/func/ses-rs3_func-bold_task-CAL_run-03_20240725120124_14.json
# task RSC 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_func-bold_task-RSC_run-01_20240725120124_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-05/func/ses-rs3_func-bold_task-RSC_run-01_20240725120124_9.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_func-bold_task-RSC_run-01_20240725120124_9.json /scratch/mszinte/data/amsterdam24/sub-001/ses-05/func/sses-rs3_func-bold_task-RSC_run-01_20240725120124_9.json
# task RSL 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_func-bold_task-RSL_run-01_20240725120124_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-05/func/ses-rs3_func-bold_task-RSL_run-01_20240725120124_10.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_func-bold_task-RSL_run-01_20240725120124_10.json /scratch/mszinte/data/amsterdam24/sub-001/ses-05/func/ses-rs3_func-bold_task-RSL_run-01_20240725120124_10.json
# task RSR 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_func-bold_task-RSR_run-01_20240725120124_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-05/func/ses-rs3_func-bold_task-RSR_run-01_20240725120124_11.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_func-bold_task-RSR_run-01_20240725120124_11.json /scratch/mszinte/data/amsterdam24/sub-001/ses-05/func/ses-rs3_func-bold_task-RSR_run-01_20240725120124_11.json

# ses-06 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-01_20240726132003_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-06/func/ses-prfstrab_func-bold_task-pRFstrab_run-01_20240726132003_9.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-01_20240726132003_9.json /scratch/mszinte/data/amsterdam24/sub-001/ses-06/func/ses-prfstrab_func-bold_task-pRFstrab_run-01_20240726132003_9.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-02_20240726132003_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-06/func/ses-prfstrab_func-bold_task-pRFstrab_run-02_20240726132003_10.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-02_20240726132003_10.json /scratch/mszinte/data/amsterdam24/sub-001/ses-06/func/ses-prfstrab_func-bold_task-pRFstrab_run-02_20240726132003_10.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-03_20240726132003_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-06/func/ses-prfstrab_func-bold_task-pRFstrab_run-03_20240726132003_11.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-03_20240726132003_11.json /scratch/mszinte/data/amsterdam24/sub-001/ses-06/func/ses-prfstrab_func-bold_task-pRFstrab_run-03_20240726132003_11.json

# copy fieldmaps into fmap per session 
# ses-01 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-B0_acq-phdiff_run-1_20240722083343_6_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-B0_acq-phdiff_run-1_20240722083343_6_fieldmaphz.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-B0_acq-phdiff_run-1_20240722083343_6_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-B0_acq-phdiff_run-1_20240722083343_6_fieldmaphz.json 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-B0_acq-phdiff_run-2_20240722083343_15_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-B0_acq-phdiff_run-2_20240722083343_15_fieldmaphz.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-B0_acq-phdiff_run-2_20240722083343_15_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-B0_acq-phdiff_run-2_20240722083343_15_fieldmaphz.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-B0_acq-phdiff_run-3_20240722083343_19_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-B0_acq-phdiff_run-3_20240722083343_19_fieldmaphz.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-B0_acq-phdiff_run-3_20240722083343_19_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-B0_acq-phdiff_run-3_20240722083343_19_fieldmaphz.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-B0_acq-phdiff_run-1_20240722083343_6.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-B0_acq-phdiff_run-1_20240722083343_6.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-B0_acq-phdiff_run-1_20240722083343_6.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-B0_acq-phdiff_run-1_20240722083343_6.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-B0_acq-phdiff_run-2_20240722083343_15.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-B0_acq-phdiff_run-2_20240722083343_15.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-B0_acq-phdiff_run-2_20240722083343_15.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-B0_acq-phdiff_run-2_20240722083343_15.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-B0_acq-phdiff_run-3_20240722083343_19.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-B0_acq-phdiff_run-3_20240722083343_19.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-B0_acq-phdiff_run-3_20240722083343_19.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-B0_acq-phdiff_run-3_20240722083343_19.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-epi_acq-SE_dir-AP_run-01_20240722083343_7.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-epi_acq-SE_dir-AP_run-01_20240722083343_7.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-epi_acq-SE_dir-AP_run-01_20240722083343_7.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-epi_acq-SE_dir-AP_run-01_20240722083343_7.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-epi_acq-SE_dir-AP_run-02_20240722083343_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-epi_acq-SE_dir-AP_run-02_20240722083343_11.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-epi_acq-SE_dir-AP_run-02_20240722083343_11.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-epi_acq-SE_dir-AP_run-02_20240722083343_11.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-epi_acq-SE_dir-AP_run-06_20240722083343_20.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-epi_acq-SE_dir-AP_run-06_20240722083343_20.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-epi_acq-SE_dir-AP_run-06_20240722083343_20.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-epi_acq-SE_dir-AP_run-06_20240722083343_20.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-epi_acq-SE_dir-PA_run-01_20240722083343_8.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-epi_acq-SE_dir-PA_run-01_20240722083343_8.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-epi_acq-SE_dir-PA_run-01_20240722083343_8.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-epi_acq-SE_dir-PA_run-01_20240722083343_8.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-epi_acq-SE_dir-PA_run-02_20240722083343_12.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-epi_acq-SE_dir-PA_run-02_20240722083343_12.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-epi_acq-SE_dir-PA_run-02_20240722083343_12.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-epi_acq-SE_dir-PA_run-02_20240722083343_12.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-epi_acq-SE_dir-PA_run-06_20240722083343_21.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-epi_acq-SE_dir-PA_run-06_20240722083343_21.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rsprf/ses-rsprf_fmap-epi_acq-SE_dir-PA_run-06_20240722083343_21.json /scratch/mszinte/data/amsterdam24/sub-001/ses-01/fmap/ses-rsprf_fmap-epi_acq-SE_dir-PA_run-06_20240722083343_21.json

# ses-02 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfoccl/ses-prfoccl_fmap-B0_acq-phdiff_run-1_20240722123602_9_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-02/fmap/ses-prfoccl_fmap-B0_acq-phdiff_run-1_20240722123602_9_fieldmaphz.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfoccl/ses-prfoccl_fmap-B0_acq-phdiff_run-1_20240722123602_9_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-001/ses-02/fmap/ses-prfoccl_fmap-B0_acq-phdiff_run-1_20240722123602_9_fieldmaphz.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfoccl/ses-prfoccl_fmap-B0_acq-phdiff_run-1_20240722123602_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-02/fmap/ses-prfoccl_fmap-B0_acq-phdiff_run-1_20240722123602_9.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfoccl/ses-prfoccl_fmap-B0_acq-phdiff_run-1_20240722123602_9.json /scratch/mszinte/data/amsterdam24/sub-001/ses-02/fmap/ses-prfoccl_fmap-B0_acq-phdiff_run-1_20240722123602_9.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfoccl/ses-prfoccl_fmap-epi_acq-SE_dir-AP_run-01_20240722123602_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-02/fmap/ses-prfoccl_fmap-epi_acq-SE_dir-AP_run-01_20240722123602_10.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfoccl/ses-prfoccl_fmap-epi_acq-SE_dir-AP_run-01_20240722123602_10.json /scratch/mszinte/data/amsterdam24/sub-001/ses-02/fmap/ses-prfoccl_fmap-epi_acq-SE_dir-AP_run-01_20240722123602_10.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfoccl/ses-prfoccl_fmap-epi_acq-SE_dir-PA_run-01_20240722123602_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-02/fmap/ses-prfoccl_fmap-epi_acq-SE_dir-PA_run-01_20240722123602_11.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfoccl/ses-prfoccl_fmap-epi_acq-SE_dir-PA_run-01_20240722123602_11.json /scratch/mszinte/data/amsterdam24/sub-001/ses-02/fmap/ses-prfoccl_fmap-epi_acq-SE_dir-PA_run-01_20240722123602_11.json

# ses-03 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_fmap-B0_acq-phdiff_run-1_20240722163325_6_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-03/fmap/ses-rs2_fmap-B0_acq-phdiff_run-1_20240722163325_6_fieldmaphz.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_fmap-B0_acq-phdiff_run-1_20240722163325_6_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-001/ses-03/fmap/ses-rs2_fmap-B0_acq-phdiff_run-1_20240722163325_6_fieldmaphz.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_fmap-B0_acq-phdiff_run-1_20240722163325_6.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-03/fmap/ses-rs2_fmap-B0_acq-phdiff_run-1_20240722163325_6.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_fmap-B0_acq-phdiff_run-1_20240722163325_6.json /scratch/mszinte/data/amsterdam24/sub-001/ses-03/fmap/ses-rs2_fmap-B0_acq-phdiff_run-1_20240722163325_6.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_fmap-epi_acq-SE_dir-AP_run-01_20240722163325_7.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-03/fmap/ses-rs2_fmap-epi_acq-SE_dir-AP_run-01_20240722163325_7.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_fmap-epi_acq-SE_dir-AP_run-01_20240722163325_7.json /scratch/mszinte/data/amsterdam24/sub-001/ses-03/fmap/ses-rs2_fmap-epi_acq-SE_dir-AP_run-01_20240722163325_7.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_fmap-epi_acq-SE_dir-AP_run-04_20240722163325_16.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-03/fmap/ses-rs2_fmap-epi_acq-SE_dir-AP_run-04_20240722163325_16.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_fmap-epi_acq-SE_dir-AP_run-04_20240722163325_16.json /scratch/mszinte/data/amsterdam24/sub-001/ses-03/fmap/ses-rs2_fmap-epi_acq-SE_dir-AP_run-04_20240722163325_16.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_fmap-epi_acq-SE_dir-PA_run-01_20240722163325_8.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-03/fmap/ses-rs2_fmap-epi_acq-SE_dir-PA_run-01_20240722163325_8.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_fmap-epi_acq-SE_dir-PA_run-01_20240722163325_8.json /scratch/mszinte/data/amsterdam24/sub-001/ses-03/fmap/ses-rs2_fmap-epi_acq-SE_dir-PA_run-01_20240722163325_8.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_fmap-epi_acq-SE_dir-PA_run-04_20240722163325_17.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-03/fmap/ses-rs2_fmap-epi_acq-SE_dir-PA_run-04_20240722163325_17.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs2/ses-rs2_fmap-epi_acq-SE_dir-PA_run-04_20240722163325_17.json /scratch/mszinte/data/amsterdam24/sub-001/ses-03/fmap/ses-rs2_fmap-epi_acq-SE_dir-PA_run-04_20240722163325_178.json

# ses-04
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfaniso/ses-prfaniso_fmap-B0_acq-phdiff_run-1_20240724102635_6_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-04/fmap/ses-prfaniso_fmap-B0_acq-phdiff_run-1_20240724102635_6_fieldmaphz.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfaniso/ses-prfaniso_fmap-B0_acq-phdiff_run-1_20240724102635_6_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-001/ses-04/fmap/ses-prfaniso_fmap-B0_acq-phdiff_run-1_20240724102635_6_fieldmaphz.json 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfaniso/ses-prfaniso_fmap-B0_acq-phdiff_run-1_20240724102635_6.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-04/fmap/ses-prfaniso_fmap-B0_acq-phdiff_run-1_20240724102635_6.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfaniso/ses-prfaniso_fmap-B0_acq-phdiff_run-1_20240724102635_6.json /scratch/mszinte/data/amsterdam24/sub-001/ses-04/fmap/ses-prfaniso_fmap-B0_acq-phdiff_run-1_20240724102635_6.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfaniso/ses-prfaniso_fmap-epi_acq-SE_dir-AP_run-01_20240724102635_7.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-04/fmap/ses-prfaniso_fmap-epi_acq-SE_dir-AP_run-01_20240724102635_7.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfaniso/ses-prfaniso_fmap-epi_acq-SE_dir-AP_run-01_20240724102635_7.json /scratch/mszinte/data/amsterdam24/sub-001/ses-04/fmap/ses-prfaniso_fmap-epi_acq-SE_dir-AP_run-01_20240724102635_7.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfaniso/ses-prfaniso_fmap-epi_acq-SE_dir-PA_run-01_20240724102635_8.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-04/fmap/ses-prfaniso_fmap-epi_acq-SE_dir-PA_run-01_20240724102635_8.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfaniso/ses-prfaniso_fmap-epi_acq-SE_dir-PA_run-01_20240724102635_8.json /scratch/mszinte/data/amsterdam24/sub-001/ses-04/fmap/ses-prfaniso_fmap-epi_acq-SE_dir-PA_run-01_20240724102635_8.json

# ses-05 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_fmap-B0_acq-phdiff_run-1_20240725120124_6_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-05/fmap/ses-rs3_fmap-B0_acq-phdiff_run-1_20240725120124_6_fieldmaphz.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_fmap-B0_acq-phdiff_run-1_20240725120124_6_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-001/ses-05/fmap/ses-rs3_fmap-B0_acq-phdiff_run-1_20240725120124_6_fieldmaphz.json 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_fmap-B0_acq-phdiff_run-1_20240725120124_6.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-05/fmap/ses-rs3_fmap-B0_acq-phdiff_run-1_20240725120124_6.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_fmap-B0_acq-phdiff_run-1_20240725120124_6.json /scratch/mszinte/data/amsterdam24/sub-001/ses-05/fmap/ses-rs3_fmap-B0_acq-phdiff_run-1_20240725120124_6.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_fmap-epi_acq-SE_dir-AP_run-01_20240725120124_7.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-05/fmap/ses-rs3_fmap-epi_acq-SE_dir-AP_run-01_20240725120124_7.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_fmap-epi_acq-SE_dir-AP_run-01_20240725120124_7.json /scratch/mszinte/data/amsterdam24/sub-001/ses-05/fmap/ses-rs3_fmap-epi_acq-SE_dir-AP_run-01_20240725120124_7.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_fmap-epi_acq-SE_dir-PA_run-01_20240725120124_8.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-05/fmap/ses-rs3_fmap-epi_acq-SE_dir-PA_run-01_20240725120124_8.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-rs3/ses-rs3_fmap-epi_acq-SE_dir-PA_run-01_20240725120124_8.json /scratch/mszinte/data/amsterdam24/sub-001/ses-05/fmap/ses-rs3_fmap-epi_acq-SE_dir-PA_run-01_20240725120124_8.json

# ses-06 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfstrab/ses-prfstrab_fmap-B0_acq-phdiff_run-1_20240726132003_6_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-06/fmap/ses-prfstrab_fmap-B0_acq-phdiff_run-1_20240726132003_6_fieldmaphz.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfstrab/ses-prfstrab_fmap-B0_acq-phdiff_run-1_20240726132003_6_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-001/ses-06/fmap/ses-prfstrab_fmap-B0_acq-phdiff_run-1_20240726132003_6_fieldmaphz.json 
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfstrab/ses-prfstrab_fmap-B0_acq-phdiff_run-1_20240726132003_6.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-06/fmap/ses-prfstrab_fmap-B0_acq-phdiff_run-1_20240726132003_6.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfstrab/ses-prfstrab_fmap-B0_acq-phdiff_run-1_20240726132003_6.json /scratch/mszinte/data/amsterdam24/sub-001/ses-06/fmap/ses-prfstrab_fmap-B0_acq-phdiff_run-1_20240726132003_6.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfstrab/ses-prfstrab_fmap-epi_acq-SE_dir-AP_run-01_20240726132003_7.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-06/fmap/ses-prfstrab_fmap-epi_acq-SE_dir-AP_run-01_20240726132003_7.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfstrab/ses-prfstrab_fmap-epi_acq-SE_dir-AP_run-01_20240726132003_7.json /scratch/mszinte/data/amsterdam24/sub-001/ses-06/fmap/ses-prfstrab_fmap-epi_acq-SE_dir-AP_run-01_20240726132003_7.json
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfstrab/ses-prfstrab_fmap-epi_acq-SE_dir-PA_run-01_20240726132003_8.nii.gz /scratch/mszinte/data/amsterdam24/sub-001/ses-06/fmap/ses-prfstrab_fmap-epi_acq-SE_dir-PA_run-01_20240726132003_8.nii.gz
rsync -av --progress /scratch/mszinte/data/raw_data/prfGaze/sourcedata/sub-001/ses-prfstrab/ses-prfstrab_fmap-epi_acq-SE_dir-PA_run-01_20240726132003_8.json /scratch/mszinte/data/amsterdam24/sub-001/ses-06/fmap/ses-prfstrab_fmap-epi_acq-SE_dir-PA_run-01_20240726132003_8.json


# session mapping sub-001 
#       = ses-01 
#     = ses-02 
#         = ses-03 
#   = ses-04 
#       = ses-05 
#    = ses-06 



#json side cares (use _add_missing_BIDS_metadata_and_save_to_disk?)





































