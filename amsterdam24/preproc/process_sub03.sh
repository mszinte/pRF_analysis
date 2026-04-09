#Process amsterdam24

# sub 03
# convert to nifti per session and subject (compress & make json sidecar)
dcm2niix -z y -b y /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfaniso
dcm2niix -z y -b y /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfoccl
dcm2niix -z y -b y /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfstrab
dcm2niix -z y -b y /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs
dcm2niix -z y -b y /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2
dcm2niix -z y -b y /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs3


# rename and restructure in new datarepo (amsterdam24) BIDS 
# ignore ph files 
# sub-XXX/
#    ses-YY/
#        func/
#        fmap/
mkdir  /scratch/mszinte/data/amsterdam24/sub-03 
mkdir  /scratch/mszinte/data/amsterdam24/sub-03/ses-01
mkdir  /scratch/mszinte/data/amsterdam24/sub-03/ses-01/func 
mkdir  /scratch/mszinte/data/amsterdam24/sub-03/ses-01/fmap
#same for all other sessions ! 

# session mapping sub-003
#"sub-003": {
#"ses-01": "ses-prf_strab", (22nd afternoon)
#"ses-02": "ses-rs", (23 afternoon)
#"ses-03": "ses-prfaniso", (23 evening)
#"ses-04": "ses-rs2", (24 morning)
#"ses-05": "ses_occl", (25 midday)
#"ses-06": "ses-rs3", (25 afternoon)

# copy func files into func per session 
# ses-01
# task pRFstrab
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-01_20240722145747_12.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-01_20240722145747_12.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-01_20240722145747_12.json /scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-01_20240722145747_12.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-02_20240722145747_13.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-02_20240722145747_13.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-02_20240722145747_13.json /scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-02_20240722145747_13.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-03_20240722145747_14.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-03_20240722145747_14.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfstrab/ses-prfstrab_func-bold_task-pRFstrab_run-03_20240722145747_14.json /scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-03_20240722145747_14.json
# ses-02 
# task RSC
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_func-bold_task-RSC_run-01_20240723160853_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSC_run-01_20240723160853_9.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_func-bold_task-RSC_run-01_20240723160853_9.json /scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSC_run-01_20240723160853_9.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_func-bold_task-RSC_run-02_20240723160853_12.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSC_run-02_20240723160853_12.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_func-bold_task-RSC_run-02_20240723160853_12.json /scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSC_run-02_20240723160853_12.json
# task RSL
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_func-bold_task-RSL_run-01_20240723160853_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSL_run-01_20240723160853_10.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_func-bold_task-RSL_run-01_20240723160853_10.json /scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSL_run-01_20240723160853_10.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_func-bold_task-RSL_run-02_20240723160853_13.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSL_run-02_20240723160853_13.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_func-bold_task-RSL_run-02_20240723160853_13.json /scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSL_run-02_20240723160853_13.json
# task RSR
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_func-bold_task-RSR_run-01_20240723160853_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSR_run-01_20240723160853_11.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_func-bold_task-RSR_run-01_20240723160853_11.json /scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSR_run-01_20240723160853_11.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_func-bold_task-RSR_run-02_20240723160853_14.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSR_run-02_20240723160853_14.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_func-bold_task-RSR_run-02_20240723160853_14.json /scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSR_run-02_20240723160853_14.json
# ses-03
# task ses-prfaniso
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-01_20240723172444_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/ses-prfaniso_func-bold_task-pRFaniso_run-01_20240723172444_9.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-01_20240723172444_9.json /scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/ses-prfaniso_func-bold_task-pRFaniso_run-01_20240723172444_9.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-02_20240723172444_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/ses-prfaniso_func-bold_task-pRFaniso_run-02_20240723172444_10.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-02_20240723172444_10.json /scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/ses-prfaniso_func-bold_task-pRFaniso_run-02_20240723172444_10.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-03_20240723172444_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/ses-prfaniso_func-bold_task-pRFaniso_run-03_20240723172444_11.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfaniso/ses-prfaniso_func-bold_task-pRFaniso_run-03_20240723172444_11.json /scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/ses-prfaniso_func-bold_task-pRFaniso_run-03_20240723172444_11.json

#ses 04
#ses-rs2
#task RSC
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-RSC_run-01_20240724090329_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-RSC_run-01_20240724090329_9.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-RSC_run-01_20240724090329_9.json /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-RSC_run-01_20240724090329_9.json
# task RSL
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-RSL_run-01_20240724090329_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-RSL_run-01_20240724090329_10.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-RSL_run-01_20240724090329_10.json /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-RSL_run-01_20240724090329_10.json
# task RSR
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-RSR_run-01_20240724090329_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-RSR_run-01_20240724090329_11.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-RSR_run-01_20240724090329_11.json /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-RSR_run-01_20240724090329_11.json
# task pRF 
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-pRF_run-01_20240724090329_15.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-pRF_run-01_20240724090329_15.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-pRF_run-01_20240724090329_15.json /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-pRF_run-01_20240724090329_15.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-pRF_run-02_20240724090329_16.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-pRF_run-02_20240724090329_16.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-pRF_run-02_20240724090329_16.json /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-pRF_run-02_20240724090329_16.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-pRF_run-03_20240724090329_17.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-pRF_run-03_20240724090329_17.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-pRF_run-03_20240724090329_17.json /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-pRF_run-03_20240724090329_17.json
# task CAL 
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-CAL_run-01_20240724090329_12.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-CAL_run-01_20240724090329_12.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-CAL_run-01_20240724090329_12.json /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-CAL_run-01_20240724090329_12.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-CAL_run-03_20240724090329_13.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-CAL_run-02_20240724090329_13.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-CAL_run-03_20240724090329_13.json /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-CAL_run-02_20240724090329_13.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-CAL_run-03_20240724090329_14.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-CAL_run-03_20240724090329_14.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_func-bold_task-CAL_run-03_20240724090329_14.json /scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-CAL_run-03_20240724090329_14.json

# ses-05 
# task prfoccl
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfoccl/ses-prfoccl_func-bold_task-pRFoccl_run-01_20240726124505_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/ses-prfoccl_func-bold_task-pRFoccl_run-01_20240726124505_9.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfoccl/ses-prfoccl_func-bold_task-pRFoccl_run-01_20240726124505_9.json /scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/ses-prfoccl_func-bold_task-pRFoccl_run-01_20240726124505_9.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfoccl/ses-prfoccl_func-bold_task-pRFoccl_run-02_20240726124505_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/ses-prfoccl_func-bold_task-pRFoccl_run-02_20240726124505_10.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfoccl/ses-prfoccl_func-bold_task-pRFoccl_run-02_20240726124505_10.json /scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/ses-prfoccl_func-bold_task-pRFoccl_run-02_20240726124505_10.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfoccl/ses-prfoccl_func-bold_task-pRFoccl_run-03_20240726124505_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/ses-prfoccl_func-bold_task-pRFoccl_run-03_20240726124505_11.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfoccl/ses-prfoccl_func-bold_task-pRFoccl_run-03_20240726124505_11.json /scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/ses-prfoccl_func-bold_task-pRFoccl_run-03_20240726124505_11.json


# ses-06 ses-rs3
# task RSC
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs3/ses-rs3_func-bold_task-RSC_run-01_20240726143047_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/ses-rs3_func-bold_task-RSC_run-01_20240726143047_9.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs3/ses-rs3_func-bold_task-RSC_run-01_20240726143047_9.json /scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/ses-rs3_func-bold_task-RSC_run-01_20240726143047_9.json
#task RSL 
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs3/ses-rs3_func-bold_task-RSL_run-01_20240726143047_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/ses-rs3_func-bold_task-RSL_run-01_20240726143047_10.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs3/ses-rs3_func-bold_task-RSL_run-01_20240726143047_10.json /scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/ses-rs3_func-bold_task-RSL_run-01_20240726143047_10.json
#task RSR
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs3/ses-rs3_func-bold_task-RSR_run-01_20240726143047_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/ses-rs3_func-bold_task-RSR_run-01_20240726143047_11.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs3/ses-rs3_func-bold_task-RSR_run-01_20240726143047_11.json /scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/ses-rs3_func-bold_task-RSR_run-01_20240726143047_11.json


# RENAME FUNC FILES 
# ses-01 task pRFstrab
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-01_20240722145747_12.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/sub-03_ses-01_task-pRFstrab_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-01_20240722145747_12.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/sub-03_ses-01_task-pRFstrab_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-02_20240722145747_13.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/sub-03_ses-01_task-pRFstrab_run-02_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-02_20240722145747_13.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/sub-03_ses-01_task-pRFstrab_run-02_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-03_20240722145747_14.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/sub-03_ses-01_task-pRFstrab_run-03_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/ses-prfstrab_func-bold_task-pRFstrab_run-03_20240722145747_14.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-01/func/sub-03_ses-01_task-pRFstrab_run-03_bold.json"


# ses-02 task RSC, RSL, RSR
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSC_run-01_20240723160853_9.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/sub-03_ses-02_task-RSC_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSC_run-01_20240723160853_9.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/sub-03_ses-02_task-RSC_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSC_run-02_20240723160853_12.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/sub-03_ses-02_task-RSC_run-02_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSC_run-02_20240723160853_12.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/sub-03_ses-02_task-RSC_run-02_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSL_run-01_20240723160853_10.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/sub-03_ses-02_task-RSL_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSL_run-01_20240723160853_10.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/sub-03_ses-02_task-RSL_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSL_run-02_20240722153846_13.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/sub-03_ses-02_task-RSL_run-02_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSL_run-02_20240722153846_13.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/sub-03_ses-02_task-RSL_run-02_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSR_run-01_20240723160853_11.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/sub-03_ses-02_task-RSR_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSR_run-01_20240723160853_11.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/sub-03_ses-02_task-RSR_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSR_run-02_20240723160853_14.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/sub-03_ses-02_task-RSR_run-02_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/ses-rs_func-bold_task-RSR_run-02_20240723160853_14.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-02/func/sub-03_ses-02_task-RSR_run-02_bold.json"

# ses-03 task ses-prfaniso
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/ses-prfaniso_func-bold_task-pRFaniso_run-01_20240723172444_9.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/sub-03_ses-03_task-pRFaniso_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/ses-prfaniso_func-bold_task-pRFaniso_run-01_20240723172444_9.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/sub-03_ses-03_task-pRFaniso_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/ses-prfaniso_func-bold_task-pRFaniso_run-02_20240723172444_10.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/sub-03_ses-03_task-pRFaniso_run-02_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/ses-prfaniso_func-bold_task-pRFaniso_run-02_20240723172444_10.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/sub-03_ses-03_task-pRFaniso_run-02_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/ses-prfaniso_func-bold_task-pRFaniso_run-03_20240723172444_11.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/sub-03_ses-03_task-pRFaniso_run-03_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/ses-prfaniso_func-bold_task-pRFaniso_run-03_20240723172444_11.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-03/func/sub-03_ses-03_task-pRFaniso_run-03_bold.json"

# ses-04 ses-rscalprf
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-RSC_run-01_20240724090329_9.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-RSC_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-RSC_run-01_20240724090329_9.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-RSC_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-RSL_run-01_20240724090329_10.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-RSL_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-RSL_run-01_20240724090329_10.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-RSL_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-RSR_run-01_20240724090329_11.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-RSR_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-RSR_run-01_20240724090329_11.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-RSR_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-pRF_run-01_20240724090329_15.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-pRF_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-pRF_run-01_20240724090329_15.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-pRF_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-pRF_run-02_20240724090329_16.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-pRF_run-02_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-pRF_run-02_20240724090329_16.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-pRF_run-02_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-pRF_run-03_20240724090329_17.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-pRF_run-03_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-pRF_run-03_20240724090329_17.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-pRF_run-03_bold.json"

mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-CAL_run-01_20240724090329_12.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-CAL_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-CAL_run-01_20240724090329_12.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-CAL_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-CAL_run-02_20240724090329_13.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-CAL_run-02_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-CAL_run-02_20240724090329_13.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-CAL_run-02_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-CAL_run-03_20240724090329_14.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-CAL_run-03_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/ses-rs2_func-bold_task-CAL_run-03_20240724090329_14.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-04/func/sub-03_ses-04_task-CAL_run-03_bold.json"


# ses-05  task prfoccl
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/ses-prfoccl_func-bold_task-pRFoccl_run-01_20240726124505_9.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/sub-03_ses-05_task-pRFoccl_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/ses-prfoccl_func-bold_task-pRFoccl_run-01_20240726124505_9.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/sub-03_ses-05_task-pRFoccl_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/ses-prfoccl_func-bold_task-pRFoccl_run-02_20240726124505_10.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/sub-03_ses-05_task-pRFoccl_run-02_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/ses-prfoccl_func-bold_task-pRFoccl_run-02_20240726124505_10.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/sub-03_ses-05_task-pRFoccl_run-02_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/ses-prfoccl_func-bold_task-pRFoccl_run-03_20240726124505_11.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/sub-03_ses-05_task-pRFoccl_run-03_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/ses-prfoccl_func-bold_task-pRFoccl_run-03_20240726124505_11.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-05/func/sub-03_ses-05_task-pRFoccl_run-03_bold.json"


# ses-06 RSC, RSL, RSR
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/ses-rs3_func-bold_task-RSC_run-01_20240726143047_9.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/sub-03_ses-06_task-RSC_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/ses-rs3_func-bold_task-RSC_run-01_20240726143047_9.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/sub-03_ses-06_task-RSC_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/ses-rs3_func-bold_task-RSL_run-01_20240726143047_10.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/sub-03_ses-06_task-RSL_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/ses-rs3_func-bold_task-RSL_run-01_20240726143047_10.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/sub-03_ses-06_task-RSL_run-01_bold.json"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/ses-rs3_func-bold_task-RSR_run-01_20240726143047_11.nii.gz" "/scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/sub-03_ses-06_task-RSR_run-01_bold.nii.gz"
mv "/scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/ses-rs3_func-bold_task-RSR_run-01_20240726143047_11.json" "/scratch/mszinte/data/amsterdam24/sub-03/ses-06/func/sub-03_ses-06_task-RSR_run-01_bold.json"



# copy fieldmaps AND RENAME into fmap per session (use only LAST ones for each type if multiple runs)
# ses-01 ses-prfstrab
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfstrab/ses-prfstrab_fmap-B0_acq-phdiff_run-1_20240722145747_9_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-01/fmap/sub-03_ses-01_fieldmap.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfstrab/ses-prfstrab_fmap-B0_acq-phdiff_run-1_20240722145747_9_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-03/ses-01/fmap/sub-03_ses-01_fieldmap.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfstrab/ses-prfstrab_fmap-B0_acq-phdiff_run-1_20240722145747_9.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-01/fmap/sub-03_ses-01_magnitude.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfstrab/ses-prfstrab_fmap-epi_acq-SE_dir-AP_run-01_20240722145747_10.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-01/fmap/sub-03_ses-01_dir-AP_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfstrab/ses-prfstrab_fmap-epi_acq-SE_dir-AP_run-01_20240722145747_10.json /scratch/mszinte/data/amsterdam24/sub-03/ses-01/fmap/sub-03_ses-01_dir-AP_epi.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfstrab/ses-prfstrab_fmap-epi_acq-SE_dir-PA_run-01_20240722145747_11.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-01/fmap/sub-03_ses-01_dir-PA_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfstrab/ses-prfstrab_fmap-epi_acq-SE_dir-PA_run-01_20240722145747_11.json /scratch/mszinte/data/amsterdam24/sub-03/ses-01/fmap/sub-03_ses-01_dir-PA_epi.json


# ses-02 ses-rs
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_fmap-B0_acq-phdiff_run-2_20240723160853_15_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-02/fmap/sub-03_ses-02_fieldmap.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_fmap-B0_acq-phdiff_run-2_20240723160853_15_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-03/ses-02/fmap/sub-03_ses-02_fieldmap.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_fmap-B0_acq-phdiff_run-2_20240723160853_15.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-02/fmap/sub-03_ses-02_magnitude.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_fmap-epi_acq-SE_dir-AP_run-04_20240723160853_16.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-02/fmap/sub-03_ses-02_dir-AP_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_fmap-epi_acq-SE_dir-AP_run-04_20240723160853_16.json /scratch/mszinte/data/amsterdam24/sub-03/ses-02/fmap/sub-03_ses-02_dir-AP_epi.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_fmap-epi_acq-SE_dir-PA_run-04_20240723160853_17.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-02/fmap/sub-03_ses-02_dir-PA_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs/ses-rs_fmap-epi_acq-SE_dir-PA_run-04_20240723160853_17.json /scratch/mszinte/data/amsterdam24/sub-03/ses-02/fmap/sub-03_ses-02_dir-PA_epi.json

# ses-03 ses-prfaniso
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfaniso/ses-prfaniso_fmap-B0_acq-phdiff_run-1_20240723172444_6_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-03/fmap/sub-03_ses-03_fieldmap.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfaniso/ses-prfaniso_fmap-B0_acq-phdiff_run-1_20240723172444_6_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-03/ses-03/fmap/sub-03_ses-03_fieldmap.json 
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfaniso/ses-prfaniso_fmap-B0_acq-phdiff_run-1_20240723172444_6.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-03/fmap/sub-03_ses-03_magnitude.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfaniso/ses-prfaniso_fmap-epi_acq-SE_dir-AP_run-01_20240723172444_7.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-03/fmap/sub-03_ses-03_dir-AP_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfaniso/ses-prfaniso_fmap-epi_acq-SE_dir-AP_run-01_20240723172444_7.json /scratch/mszinte/data/amsterdam24/sub-03/ses-03/fmap/sub-03_ses-03_dir-AP_epi.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfaniso/ses-prfaniso_fmap-epi_acq-SE_dir-PA_run-01_20240723172444_8.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-03/fmap/sub-03_ses-03_dir-PA_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfaniso/ses-prfaniso_fmap-epi_acq-SE_dir-PA_run-01_20240723172444_8.json /scratch/mszinte/data/amsterdam24/sub-03/ses-03/fmap/sub-03_ses-03_dir-PA_epi.json


# ses-04 ses-rs2
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_fmap-B0_acq-phdiff_run-1_20240724090329_6_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-04/fmap/sub-03_ses-04_fieldmap.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_fmap-B0_acq-phdiff_run-1_20240724090329_6_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-03/ses-04/fmap/sub-03_ses-04_fieldmap.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_fmap-B0_acq-phdiff_run-1_20240724090329_6.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-04/fmap/sub-03_ses-04_magnitude.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_fmap-epi_acq-SE_dir-AP_run-01_20240724090329_7.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-04/fmap/sub-03_ses-04_dir-AP_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_fmap-epi_acq-SE_dir-AP_run-01_20240724090329_7.json /scratch/mszinte/data/amsterdam24/sub-03/ses-04/fmap/sub-03_ses-04_dir-AP_epi.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_fmap-epi_acq-SE_dir-PA_run-01_20240724090329_8.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-04/fmap/sub-03_ses-04_dir-PA_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs2/ses-rs2_fmap-epi_acq-SE_dir-PA_run-01_20240724090329_8.json /scratch/mszinte/data/amsterdam24/sub-03/ses-04/fmap/sub-03_ses-04_dir-PA_epi.json

# ses-05 ses-prfoccl
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfoccl/ses-prfoccl_fmap-B0_acq-phdiff_run-1_20240726124505_6_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-05/fmap/sub-03_ses-05_fieldmap.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfoccl/ses-prfoccl_fmap-B0_acq-phdiff_run-1_20240726124505_6_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-03/ses-05/fmap/sub-03_ses-05_fieldmap.json 
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfoccl/ses-prfoccl_fmap-B0_acq-phdiff_run-1_20240726124505_6.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-05/fmap/sub-03_ses-05_magnitude.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfoccl/ses-prfoccl_fmap-epi_acq-SE_dir-AP_run-01_20240726124505_7.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-05/fmap/sub-03_ses-05_dir-AP_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfoccl/ses-prfoccl_fmap-epi_acq-SE_dir-AP_run-01_20240726124505_7.json /scratch/mszinte/data/amsterdam24/sub-03/ses-05/fmap/sub-03_ses-05_dir-AP_epi.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfoccl/ses-prfoccl_fmap-epi_acq-SE_dir-PA_run-01_20240726124505_8.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-05/fmap/sub-03_ses-05_dir-PA_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-prfoccl/ses-prfoccl_fmap-epi_acq-SE_dir-PA_run-01_20240726124505_8.json /scratch/mszinte/data/amsterdam24/sub-03/ses-05/fmap/sub-03_ses-05_dir-PA_epi.json


# ses-06 ses-rs3
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs3/ses-rs3_fmap-B0_acq-phdiff_run-1_20240726143047_6_fieldmaphz.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-06/fmap/sub-03_ses-06_fieldmap.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs3/ses-rs3_fmap-B0_acq-phdiff_run-1_20240726143047_6_fieldmaphz.json /scratch/mszinte/data/amsterdam24/sub-03/ses-06/fmap/sub-03_ses-06_fieldmap.json 
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs3/ses-rs3_fmap-B0_acq-phdiff_run-1_20240726143047_6.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-06/fmap/sub-03_ses-06_magnitude.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs3/ses-rs3_fmap-epi_acq-SE_dir-AP_run-01_20240726143047_7.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-06/fmap/sub-03_ses-06_dir-AP_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs3/ses-rs3_fmap-epi_acq-SE_dir-AP_run-01_20240726143047_7.json /scratch/mszinte/data/amsterdam24/sub-03/ses-06/fmap/sub-03_ses-06_dir-AP_epi.json
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs3/ses-rs3_fmap-epi_acq-SE_dir-PA_run-01_20240726143047_8.nii.gz /scratch/mszinte/data/amsterdam24/sub-03/ses-06/fmap/sub-03_ses-06_dir-PA_epi.nii.gz
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-rs3/ses-rs3_fmap-epi_acq-SE_dir-PA_run-01_20240726143047_8.json /scratch/mszinte/data/amsterdam24/sub-03/ses-06/fmap/sub-03_ses-06_dir-PA_epi.json


#json side cares 
#see update_metadata.py 


#copy event files (and rename)
# CAL: https://github.com/mszinte/deepmreyecalib/blob/amsterdam-24/experiment_code/data/
# Rename sub-03_ses-01_task-DeepMReyeCalibAms_run-01_events.tsv to sub-03_ses-03_task-CAL_run-01_events.tsv

#prf experiments: https://github.com/mszinte/pRFexp/tree/pRFexp_altVision/data/


#copy anatomy 
#sub -03: from deepmreye
rsync -av --progress /scratch/mszinte/data/deepmreye/sub-03/ses-02/anat /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-anat/
rsync -av --progress /scratch/mszinte/data/amsterdam24/sourcedata/sub-003/ses-anat/anat /scratch/mszinte/data/amsterdam24/sub-03/ses-07/
#and rename


python fmriprep_sbatch.py /scratch/mszinte/data amsterdam24 sub-03 30 anat_only_y aroma_n fmapfree_n skip_bids_val_y cifti_output_170k_y fsaverage_y 12 sina.kling@etu.univ-amu.fr 327 b327 fmriprep-25.2.0.simg































