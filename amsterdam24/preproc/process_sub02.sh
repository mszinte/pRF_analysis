#Process amsterdam24

# sub 02 
# convert to nifti per session and subject (compress & make json sidecar)
dcm2niix -z y -b y /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfaniso
dcm2niix -z y -b y /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf
dcm2niix -z y -b y /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl
dcm2niix -z y -b y /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl2
dcm2niix -z y -b y /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs
dcm2niix -z y -b y /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfstrab
dcm2niix -z y -b y /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs2

# rename and restructure in new datarepo (amsterdam24) BIDS 
# ignore ph files 
# sub-XXX/
#    ses-YY/
#        func/
#        fmap/
mkdir  /scratch/mszinte/data/amsterdam24/sub-02 
mkdir  /scratch/mszinte/data/amsterdam24/sub-02/ses-01
mkdir  /scratch/mszinte/data/amsterdam24/sub-02/ses-01/func 
mkdir  /scratch/mszinte/data/amsterdam24/sub-02/ses-01/fmap
#same for all other sessions ! 

# session mapping sub-002
#  "ses-01": "ses-prfstrab", (22nd morning)
#"ses-02": "ses-rs",  (22nd afternoon)
# "ses-03": "ses-rscalprf", (23rd morning)
#"ses-04": "ses-prfoccl", (23rd afternoon)
#"ses-05": "ses-prfaniso", (25th afternoon)
#"ses-06": "ses-prfoccl2", (26th morning)
# "ses-07": "ses-rs2" (26th afternoon)

# copy func files into func per session 
# ses-01
# task pRFstrab

rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-01_20240722115705_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-01_20240722115705_9.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-01_20240722115705_9.json /scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-01_20240722115705_9.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-02_20240722115705_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-02_20240722115705_10.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-02_20240722115705_10.json /scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-02_20240722115705_10.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-03_20240722115705_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-03_20240722115705_11.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-03_20240722115705_11.json /scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-03_20240722115705_11.json
# ses-02 
# task RSC
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_func-bold_task-RSC_run-01_20240722153846_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSC_run-01_20240722153846_9.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_func-bold_task-RSC_run-01_20240722153846_9.json /scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSC_run-01_20240722153846_9.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_func-bold_task-RSC_run-02_20240722153846_12.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSC_run-02_20240722153846_12.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_func-bold_task-RSC_run-02_20240722153846_12.json /scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSC_run-02_20240722153846_12.json
# task RSL
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_func-bold_task-RSL_run-01_20240722153846_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSL_run-01_20240722153846_10.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_func-bold_task-RSL_run-01_20240722153846_10.json /scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSL_run-01_20240722153846_10.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_func-bold_task-RSL_run-02_20240722153846_13.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSL_run-02_20240722153846_13.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_func-bold_task-RSL_run-02_20240722153846_13.json /scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSL_run-02_20240722153846_13.json
# task RSR
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_func-bold_task-RSR_run-01_20240722153846_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSR_run-01_20240722153846_11.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_func-bold_task-RSR_run-01_20240722153846_11.json /scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSR_run-01_20240722153846_11.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_func-bold_task-RSR_run-02_20240722153846_14.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSR_run-02_20240722153846_14.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_func-bold_task-RSR_run-02_20240722153846_14.json /scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSR_run-02_20240722153846_14.json
# ses-03
#ses-rscalprf
#task RSC
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-RSC_run-01_20240723090248_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-RSC_run-01_20240723090248_9.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-RSC_run-01_20240723090248_9.json /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-RSC_run-01_20240723090248_9.json
# task RSL
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-RSL_run-01_20240723090248_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-RSL_run-01_20240723090248_10.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-RSL_run-01_20240723090248_10.json /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-RSL_run-01_20240723090248_10.json
# task RSR
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-RSR_run-01_20240723090248_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-RSR_run-01_20240723090248_11.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-RSR_run-01_20240723090248_11.json /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-RSR_run-01_20240723090248_11.json
# task pRF 
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-pRF_run-01_20240723090248_19.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-pRF_run-01_20240723090248_19.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-pRF_run-01_20240723090248_19.json /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-pRF_run-01_20240723090248_19.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-pRF_run-02_20240723090248_20.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-pRF_run-02_20240723090248_20.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-pRF_run-02_20240723090248_20.json /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-pRF_run-02_20240723090248_20.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-pRF_run-03_20240723090248_21.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-pRF_run-03_20240723090248_21.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-pRF_run-03_20240723090248_21.json /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-pRF_run-03_20240723090248_21.json
# task CAL 
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-CAL_run-01_20240723090248_14.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-CAL_run-01_20240723090248_14.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-CAL_run-01_20240723090248_14.json /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-CAL_run-01_20240723090248_14.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-CAL_run-02_20240723090248_15.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sses-rscalprf_func-bold_task-CAL_run-02_20240723090248_15.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-CAL_run-02_20240723090248_15.json /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-CAL_run-02_20240723090248_15.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-CAL_run-03_20240723090248_16.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-CAL_run-03_20240723090248_16.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_func-bold_task-CAL_run-03_20240723090248_16.json /scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-CAL_run-03_20240723090248_16.json

#ses 04
# task prfoccl
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl/ses-prfoccl_func-bold_task-pRFoccl_run-01_20240723165122_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-04/func/ses-prfoccl_func-bold_task-pRFoccl_run-01_20240723165122_9.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl/ses-prfoccl_func-bold_task-pRFoccl_run-01_20240723165122_9.json /scratch/mszinte/data/amsterdam24/sub-02/ses-04/func/ses-prfoccl_func-bold_task-pRFoccl_run-01_20240723165122_9.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl/ses-prfoccl_func-bold_task-pRFoccl_run-02_20240723165122_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-04/func/ses-prfoccl_func-bold_task-pRFoccl_run-02_20240723165122_10.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl/ses-prfoccl_func-bold_task-pRFoccl_run-02_20240723165122_10.json /scratch/mszinte/data/amsterdam24/sub-02/ses-04/func/ses-prfoccl_func-bold_task-pRFoccl_run-02_20240723165122_10.json

# ses-05 
# task ses-prfaniso
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-01_20240725125934_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/ses-prfaniso_func-bold_task-pRFaniso_run-01_20240725125934_9.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-01_20240725125934_9.json /scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/ses-prfaniso_func-bold_task-pRFaniso_run-01_20240725125934_9.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-02_20240725125934_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/ses-prfaniso_func-bold_task-pRFaniso_run-02_20240725125934_10.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-02_20240725125934_10.json /scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/ses-prfaniso_func-bold_task-pRFaniso_run-02_20240725125934_10.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-03_20240725125934_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/ses-prfaniso_func-bold_task-pRFaniso_run-03_20240725125934_11.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-03_20240725125934_11.json /scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/ses-prfaniso_func-bold_task-pRFaniso_run-03_20240725125934_11.json
# ses-06
#task prfoccl2
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl2/ses-prfoccl2_func-bold_task-pRFoccl_run-01_20240726121029_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/ses-prfoccl2_func-bold_task-pRFoccl_run-01_20240726121029_9.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl2/ses-prfoccl2_func-bold_task-pRFoccl_run-01_20240726121029_9.json /scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/ses-prfoccl2_func-bold_task-pRFoccl_run-01_20240726121029_9.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl2/ses-prfoccl2_func-bold_task-pRFoccl_run-02_20240726121029_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/ses-prfoccl2_func-bold_task-pRFoccl_run-02_20240726121029_10.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl2/ses-prfoccl2_func-bold_task-pRFoccl_run-02_20240726121029_10.json /scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/ses-prfoccl2_func-bold_task-pRFoccl_run-02_20240726121029_10.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl2/ses-prfoccl2_func-bold_task-pRFoccl_run-03_20240726121029_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/ses-prfoccl2_func-bold_task-pRFoccl_run-03_20240726121029_11.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl2/ses-prfoccl2_func-bold_task-pRFoccl_run-03_20240726121029_11.json /scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/ses-prfoccl2_func-bold_task-pRFoccl_run-03_20240726121029_11.json

# ses 07 
# task RSC
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs2/ses-rs2_func-bold_task-RSC_run-01_20240726140128_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/ses-rs2_func-bold_task-RSC_run-01_20240726140128_9.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs2/ses-rs2_func-bold_task-RSC_run-01_20240726140128_9.json /scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/ses-rs2_func-bold_task-RSC_run-01_20240726140128_9.json
#task RSL 
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs2/ses-rs2_func-bold_task-RSL_run-01_20240726140128_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/ses-rs2_func-bold_task-RSL_run-01_20240726140128_10.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs2/ses-rs2_func-bold_task-RSL_run-01_20240726140128_10.json /scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/ses-rs2_func-bold_task-RSL_run-01_20240726140128_10.json
#task RSR
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs2/ses-rs2_func-bold_task-RSR_run-01_20240726140128_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/ses-rs2_func-bold_task-RSR_run-01_20240726140128_11.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs2/ses-rs2_func-bold_task-RSR_run-01_20240726140128_11.json /scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/ses-rs2_func-bold_task-RSR_run-01_20240726140128_11.json


# RENAME FUNC FILES 
# ses-01 task pRFstrab
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-01_20240722115705_9.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/sub-02_ses-01_task-pRFstrab_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-01_20240722115705_9.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/sub-02_ses-01_task-pRFstrab_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-02_20240722115705_10.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/sub-02_ses-01_task-pRFstrab_run-02_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-02_20240722115705_10.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/sub-02_ses-01_task-pRFstrab_run-02_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-03_20240722115705_11.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/sub-02_ses-01_task-pRFstrab_run-03_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-03_20240722115705_11.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-01/func/sub-02_ses-01_task-pRFstrab_run-03_bold.json"


# ses-02 task RSC, RSL, RSR

mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSC_run-01_20240722153846_9.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/sub-02_ses-02_task-RSC_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSC_run-01_20240722153846_9.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/sub-02_ses-02_task-RSC_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSC_run-02_20240722153846_12.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/sub-02_ses-02_task-RSC_run-02_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSC_run-02_20240722153846_12.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/sub-02_ses-02_task-RSC_run-02_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSL_run-01_20240722153846_10.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/sub-02_ses-02_task-RSL_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSL_run-01_20240722153846_10.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/sub-02_ses-02_task-RSL_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSL_run-02_20240722153846_13.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/sub-02_ses-02_task-RSL_run-02_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSL_run-02_20240722153846_13.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/sub-02_ses-02_task-RSL_run-02_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSR_run-01_20240722153846_11.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/sub-02_ses-02_task-RSR_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSR_run-01_20240722153846_11.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/sub-02_ses-02_task-RSR_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSR_run-02_20240722153846_14.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/sub-02_ses-02_task-RSR_run-02_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/ses-rs_func-bold_task-RSR_run-02_20240722153846_14.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-02/func/sub-02_ses-02_task-RSR_run-02_bold.json"


# ses-03 ses-rscalprf

mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-RSC_run-01_20240723090248_9.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-RSC_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-RSC_run-01_20240723090248_9.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-RSC_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-RSL_run-01_20240723090248_10.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-RSL_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-RSL_run-01_20240723090248_10.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-RSL_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-RSR_run-01_20240723090248_11.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-RSR_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-RSR_run-01_20240723090248_11.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-RSR_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-pRF_run-01_20240723090248_19.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-pRF_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-pRF_run-01_20240723090248_19.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-pRF_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-pRF_run-02_20240723090248_20.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-pRF_run-02_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-pRF_run-02_20240723090248_20.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-pRF_run-02_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-pRF_run-03_20240723090248_21.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-pRF_run-03_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-pRF_run-03_20240723090248_21.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-pRF_run-03_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-CAL_run-01_20240723090248_14.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-CAL_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-CAL_run-01_20240723090248_14.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-CAL_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-CAL_run-02_20240723090248_15.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-CAL_run-02_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-CAL_run-02_20240723090248_15.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-CAL_run-02_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-CAL_run-03_20240723090248_16.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-CAL_run-03_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/ses-rscalprf_func-bold_task-CAL_run-03_20240723090248_16.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-03/func/sub-02_ses-03_task-CAL_run-03_bold.json"


# ses-04  task prfoccl
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-04/func/ses-prfoccl_func-bold_task-pRFoccl_run-01_20240723165122_9.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-04/func/sub-02_ses-04_task-pRFoccl_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-04/func/ses-prfoccl_func-bold_task-pRFoccl_run-01_20240723165122_9.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-04/func/sub-02_ses-04_task-pRFoccl_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-04/func/ses-prfoccl_func-bold_task-pRFoccl_run-02_20240723165122_10.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-04/func/sub-02_ses-04_task-pRFoccl_run-02_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-04/func/ses-prfoccl_func-bold_task-pRFoccl_run-02_20240723165122_10.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-04/func/sub-02_ses-04_task-pRFoccl_run-02_bold.json"


# ses-05 task ses-prfaniso
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/ses-prfaniso_func-bold_task-pRFaniso_run-01_20240725125934_9.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/sub-02_ses-05_task-pRFaniso_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/ses-prfaniso_func-bold_task-pRFaniso_run-01_20240725125934_9.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/sub-02_ses-05_task-pRFaniso_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/ses-prfaniso_func-bold_task-pRFaniso_run-02_20240725125934_10.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/sub-02_ses-05_task-pRFaniso_run-02_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/ses-prfaniso_func-bold_task-pRFaniso_run-02_20240725125934_10.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/sub-02_ses-05_task-pRFaniso_run-02_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/ses-prfaniso_func-bold_task-pRFaniso_run-03_20240725125934_11.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/sub-02_ses-05_task-pRFaniso_run-03_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/ses-prfaniso_func-bold_task-pRFaniso_run-03_20240725125934_11.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-05/func/sub-02_ses-05_task-pRFaniso_run-03_bold.json"


# ses-06 task prfoccl2
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/ses-prfoccl2_func-bold_task-pRFoccl_run-01_20240726121029_9.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/sub-02_ses-06_task-pRFoccl_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/ses-prfoccl2_func-bold_task-pRFoccl_run-01_20240726121029_9.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/sub-02_ses-06_task-pRFoccl_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/ses-prfoccl2_func-bold_task-pRFoccl_run-02_20240726121029_10.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/sub-02_ses-06_task-pRFoccl_run-02_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/ses-prfoccl2_func-bold_task-pRFoccl_run-02_20240726121029_10.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/sub-02_ses-06_task-pRFoccl_run-02_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/ses-prfoccl2_func-bold_task-pRFoccl_run-03_20240726121029_11.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/sub-02_ses-06_task-pRFoccl_run-03_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/ses-prfoccl2_func-bold_task-pRFoccl_run-03_20240726121029_11.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-06/func/sub-02_ses-06_task-pRFoccl_run-03_bold.json"


# ses-07 RSC, RSL, RSR
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/ses-rs2_func-bold_task-RSC_run-01_20240726140128_9.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/sub-02_ses-07_task-RSC_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/ses-rs2_func-bold_task-RSC_run-01_20240726140128_9.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/sub-02_ses-07_task-RSC_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/ses-rs2_func-bold_task-RSL_run-01_20240726140128_10.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/sub-02_ses-07_task-RSL_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/ses-rs2_func-bold_task-RSL_run-01_20240726140128_10.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/sub-02_ses-07_task-RSL_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/ses-rs2_func-bold_task-RSR_run-01_20240726140128_11.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/sub-02_ses-07_task-RSR_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/ses-rs2_func-bold_task-RSR_run-01_20240726140128_11.json" "/scratch/mszinte/data/amsterdam24/sub-02/ses-07/func/sub-02_ses-07_task-RSR_run-01_bold.json"



# copy fieldmaps AND RENAME into fmap per session (use only LAST ones for each type if multiple runs)
# ses-01 ses-prfstrab
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfstrab/ses-prfstrab_fmap-B0_acq-phdiff_run-1_20240722115705_6_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-01/fmap/sub-02_ses-01_fieldmap.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfstrab/ses-prfstrab_fmap-B0_acq-phdiff_run-1_20240722115705_6_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-02/ses-01/fmap/sub-02_ses-01_fieldmap.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfstrab/ses-prfstrab_fmap-B0_acq-phdiff_run-1_20240722115705_6.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-01/fmap/sub-02_ses-01_magnitude.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfstrab/ses-prfstrab_fmap-epi_acq-SE_dir-AP_run-01_20240722115705_7.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-01/fmap/sub-02_ses-01_dir-AP_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfstrab/ses-prfstrab_fmap-epi_acq-SE_dir-AP_run-01_20240722115705_7.json /scratch/mszinte/data/amsterdam24/sub-02/ses-01/fmap/sub-02_ses-01_dir-AP_epi.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfstrab/ses-prfstrab_fmap-epi_acq-SE_dir-PA_run-01_20240722115705_8.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-01/fmap/sub-02_ses-01_dir-PA_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfstrab/ses-prfstrab_fmap-epi_acq-SE_dir-PA_run-01_20240722115705_8.json /scratch/mszinte/data/amsterdam24/sub-02/ses-01/fmap/sub-02_ses-01_dir-PA_epi.json


# ses-02 ses-rs
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_fmap-B0_acq-phdiff_run-2_20240722153846_15_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-02/fmap/sub-02_ses-02_fieldmap.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_fmap-B0_acq-phdiff_run-2_20240722153846_15_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-02/ses-02/fmap/sub-02_ses-02_fieldmap.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_fmap-B0_acq-phdiff_run-2_20240722153846_15.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-02/fmap/sub-02_ses-02_magnitude.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_fmap-epi_acq-SE_dir-AP_run-04_20240722153846_16.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-02/fmap/sub-02_ses-02_dir-AP_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_fmap-epi_acq-SE_dir-AP_run-04_20240722153846_16.json /scratch/mszinte/data/amsterdam24/sub-02/ses-02/fmap/sub-02_ses-02_dir-AP_epi.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_fmap-epi_acq-SE_dir-PA_run-04_20240722153846_17.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-02/fmap/sub-02_ses-02_dir-PA_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs/ses-rs_fmap-epi_acq-SE_dir-PA_run-04_20240722153846_17.json /scratch/mszinte/data/amsterdam24/sub-02/ses-02/fmap/sub-02_ses-02_dir-PA_epi.json

# ses-03 ses-rscalprf
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_fmap-B0_acq-phdiff_run-3_20240723090248_22_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-03/fmap/sub-02_ses-03_fieldmap.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_fmap-B0_acq-phdiff_run-3_20240723090248_22_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-02/ses-03/fmap/sub-02_ses-03_fieldmap.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_fmap-B0_acq-phdiff_run-3_20240723090248_22.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-03/fmap/sub-02_ses-03_magnitude.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_fmap-epi_acq-SE_dir-AP_run-04_20240723090248_23.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-03/fmap/sub-02_ses-03_dir-AP_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_fmap-epi_acq-SE_dir-AP_run-04_20240723090248_23.json /scratch/mszinte/data/amsterdam24/sub-02/ses-03/fmap/sub-02_ses-03_dir-AP_epi.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_fmap-epi_acq-SE_dir-PA_run-04_20240723090248_24.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-03/fmap/sub-02_ses-03_dir-PA_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rscalprf/ses-rscalprf_fmap-epi_acq-SE_dir-PA_run-04_20240723090248_24.json /scratch/mszinte/data/amsterdam24/sub-02/ses-03/fmap/sub-02_ses-03_dir-PA_epi.json

# ses-04 ses-prfoccl
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl/ses-prfoccl_fmap-B0_acq-phdiff_run-1_20240723165122_6_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-04/fmap/sub-02_ses-04_fieldmap.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl/ses-prfoccl_fmap-B0_acq-phdiff_run-1_20240723165122_6_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-02/ses-04/fmap/sub-02_ses-04_fieldmap.json 
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl/ses-prfoccl_fmap-B0_acq-phdiff_run-1_20240723165122_6.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-04/fmap/sub-02_ses-04_magnitude.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl/ses-prfoccl_fmap-epi_acq-SE_dir-AP_run-01_20240723165122_7.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-04/fmap/sub-02_ses-04_dir-AP_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl/ses-prfoccl_fmap-epi_acq-SE_dir-AP_run-01_20240723165122_7.json /scratch/mszinte/data/amsterdam24/sub-02/ses-04/fmap/sub-02_ses-04_dir-AP_epi.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl/ses-prfoccl_fmap-epi_acq-SE_dir-PA_run-01_20240723165122_8.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-04/fmap/sub-02_ses-04_dir-PA_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl/ses-prfoccl_fmap-epi_acq-SE_dir-PA_run-01_20240723165122_8.json /scratch/mszinte/data/amsterdam24/sub-02/ses-04/fmap/sub-02_ses-04_dir-PA_epi.json

# ses-05 ses-prfaniso
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfaniso/ses-prfaniso_fmap-B0_acq-phdiff_run-1_20240725125934_6_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-05/fmap/sub-02_ses-05_fieldmap.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfaniso/ses-prfaniso_fmap-B0_acq-phdiff_run-1_20240725125934_6_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-02/ses-05/fmap/sub-02_ses-05_fieldmap.json 
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfaniso/ses-prfaniso_fmap-B0_acq-phdiff_run-1_20240725125934_6.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-05/fmap/sub-02_ses-05_magnitude.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfaniso/ses-prfaniso_fmap-epi_acq-SE_dir-AP_run-01_20240725125934_7.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-05/fmap/sub-02_ses-05_dir-AP_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfaniso/ses-prfaniso_fmap-epi_acq-SE_dir-AP_run-01_20240725125934_7.json /scratch/mszinte/data/amsterdam24/sub-02/ses-05/fmap/sub-02_ses-05_dir-AP_epi.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfaniso/ses-prfaniso_fmap-epi_acq-SE_dir-PA_run-01_20240725125934_8.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-05/fmap/sub-02_ses-05_dir-PA_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfaniso/ses-prfaniso_fmap-epi_acq-SE_dir-PA_run-01_20240725125934_8.json /scratch/mszinte/data/amsterdam24/sub-02/ses-05/fmap/sub-02_ses-05_dir-PA_epi.json

# ses-06 ses-prfoccl2
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl2/ses-prfoccl2_fmap-B0_acq-phdiff_run-1_20240726121029_6_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-06/fmap/sub-02_ses-06_fieldmap.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl2/ses-prfoccl2_fmap-B0_acq-phdiff_run-1_20240726121029_6_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-02/ses-06/fmap/sub-02_ses-06_fieldmap.json 
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl2/ses-prfoccl2_fmap-B0_acq-phdiff_run-1_20240726121029_6.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-06/fmap/sub-02_ses-06_magnitude.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl2/ses-prfoccl2_fmap-epi_acq-SE_dir-AP_run-01_20240726121029_7.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-06/fmap/sub-02_ses-06_dir-AP_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl2/ses-prfoccl2_fmap-epi_acq-SE_dir-AP_run-01_20240726121029_7.json /scratch/mszinte/data/amsterdam24/sub-02/ses-06/fmap/sub-02_ses-06_dir-AP_epi.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl2/ses-prfoccl2_fmap-epi_acq-SE_dir-PA_run-01_20240726121029_8.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-06/fmap/sub-02_ses-06_dir-PA_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-prfoccl2/ses-prfoccl2_fmap-epi_acq-SE_dir-PA_run-01_20240726121029_8.json /scratch/mszinte/data/amsterdam24/sub-02/ses-06/fmap/sub-02_ses-06_dir-PA_epi.json


# ses-07 ses-rs2
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs2/ses-rs2_fmap-B0_acq-phdiff_run-1_20240726140128_6_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-07/fmap/sub-02_ses-07_fieldmap.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs2/ses-rs2_fmap-B0_acq-phdiff_run-1_20240726140128_6_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-02/ses-07/fmap/sub-02_ses-07_fieldmap.json 
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs2/ses-rs2_fmap-B0_acq-phdiff_run-1_20240726140128_6.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-07/fmap/sub-02_ses-07_magnitude.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs2/ses-rs2_fmap-epi_acq-SE_dir-AP_run-01_20240726140128_7.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-07/fmap/sub-02_ses-07_dir-AP_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs2/ses-rs2_fmap-epi_acq-SE_dir-AP_run-01_20240726140128_7.json /scratch/mszinte/data/amsterdam24/sub-02/ses-07/fmap/sub-02_ses-07_dir-AP_epi.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs2/ses-rs2_fmap-epi_acq-SE_dir-PA_run-01_20240726140128_8.nii.gz /scratch/mszinte/data/amsterdam24/sub-02/ses-07/fmap/sub-02_ses-07_dir-PA_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-rs2/ses-rs2_fmap-epi_acq-SE_dir-PA_run-01_20240726140128_8.json /scratch/mszinte/data/amsterdam24/sub-02/ses-07/fmap/sub-02_ses-07_dir-PA_epi.json


#json side cares 
#see update_metadata.py 


#copy event files (and rename)
# CAL: https://github.com/mszinte/deepmreyecalib/blob/amsterdam-24/experiment_code/data/
# Rename sub-02_ses-01_task-DeepMReyeCalibAms_run-01_events.tsv to sub-02_ses-03_task-CAL_run-01_events.tsv

#prf experiments: https://github.com/mszinte/pRFexp/tree/pRFexp_altVision/data/


#copy anatomy 
#sub -02: from deepmreye
rsync -av --progress /scratch/mszinte/data/deepmreye/sub-02/ses-02/anat /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-anat/
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-002/ses-anat/anat /scratch/mszinte/data/amsterdam24/sub-02/ses-08/
#and rename


python fmriprep_sbatch.py /scratch/mszinte/data amsterdam24 sub-02 30 anat_only_y aroma_n fmapfree_n skip_bids_val_y cifti_output_170k_y fsaverage_y 12 sina.kling@etu.univ-amu.fr 327 b327 fmriprep-25.2.0.simg































