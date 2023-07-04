#' Preprocess sequencing data
#' 
#'
#' @param tumour Path to a single or multiple BAM files
#' @param normal Path to a normal matched BAM file
#' @param chromosomes Select chromosomes to analyze. Default NULL
#' @param env_hatchet Hatchet conda enviroment
#' @param config Hatchet configuration for each step
#' @export

run_meddic2=function(
    cn=NULL,
    env_medicc2=build_default_python_enviroment_list()$env_medicc2,
    ...
){

    run_main=function(
        .env
){
        .this.env=environment()
        append_env(to=.this.env,from=.env)

        set_main(.env=.this.env)

        .main$out_files=list(
            branch_length=paste0(
                out_file_dir,"/",input_id,
                ".MEDICC2_branch_lengths.tsv"
            ),
            cn_profiles_heatmap=paste0(
                out_file_dir,"/",input_id,
                ".MEDICC2_cn_profiles_heatmap.pdf"
            ),
            cn_profiles=paste0(
                out_file_dir,"/",input_id,
                ".MEDICC2_cn_profiles.pdf"
            ),
            cn_profiles_tsv=paste0(
                out_file_dir,"/",input_id,
                ".MEDICC2_final_cn_profiles.tsv"
            ),
            final_tree_new=paste0(
                out_file_dir,"/",input_id,
                ".MEDICC2_final_tree.new"
            ),
            final_tree_new=paste0(
                out_file_dir,"/",input_id,
                ".MEDICC2_final_tree.new"
            ),
            final_tree_png=paste0(
                out_file_dir,"/",input_id,
                ".MEDICC2_final_tree.png"
            ),
            final_tree_xml=paste0(
                out_file_dir,"/",input_id,
                ".MEDICC2_final_tree.xml"
            ),
            distance_tsv=paste0(
                out_file_dir,"/",input_id,
                ".MEDICC2_pairwise_distances.tsv"
            ),
            summry_tsv=paste0(
                out_file_dir,"/",input_id,
                ".MEDICC2_summary.tsv"
            ),
            cn_events=paste0(
                out_file_dir,"/",input_id,
                ".MEDICC2_copynumber_events_df.tsv"
            )

        )

        .main$exec_code=paste(
            set_conda_envir(env=env_medicc2),
            " python -m medicc2 ", cn, out_file_dir, " --plot \"both\" --events "
        )

        run_job(.env=.this.env)
        .main.step=.main$steps[[fn_id]]

        .env$.main <- .main   
    }
         
    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      var
    )

    launch(.env=.base.env)

}
