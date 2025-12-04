

### Model Validation metrics summary 



gam_cv <- readRDS("models/objects/Tracks_GAMM_CV_Models_mp.rds")

gam_cv
gam_cv$metrics   
gam_cv$summary   


rstatix::get_summary_stats(gam_cv$metrics)
