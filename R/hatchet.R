#' Preprocess sequencing data
#' 
#'
#' @param tumour Path to a single or multiple BAM files
#' @param normal Path to a normal matched BAM file
#' @param chromosomes Select chromosomes to analyze. Default NULL
#' @param env_hatchet Hatchet conda enviroment
#' @param config Hatchet configuration for each step
#' @export

run_hatchet=function(
    tumour=NULL,
    normal=NULL,
    chromosomes=NULL,
    env_hatchet=build_default_python_enviroment_list()$env_hatchet,
    config=build_default_hatchet_config(),
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    ...
){

    run_main=function(
        .env
){
        .this.env=environment()
        append_env(to=.this.env,from=.env)

        set_main(.env=.this.env)

        .main$out_files$hatchet_config=paste0(
            out_file_dir,"/",input_id,
            ".hatchet.config"
        )

        build_ini_hatchet=function(
            config_file=NULL
        ){
            connection<-file(config_file, "w")
            writeLines("[run]", connection)

            out=lapply(1:length(config),FUN=function(x){
                writeLines(paste0(names(config[x])," = ",stringr::str_to_title(config[[x]]$run)), connection)
            })
    
            writeLines("", connection)

            if(!is.null(threads)){
                writeLines(paste0("processes = ",threads), connection)
            }

            writeLines(paste0("chromosomes = ",chromosomes), connection)
            writeLines(paste0("normal = ",normal), connection)
            writeLines(paste0("bams = ",paste0(tumour,collapse="\t")), connection)
            writeLines(paste0("samples = ",paste0(Vectorize(get_file_name)(tumour),collapse="\t")), connection)
            writeLines(paste0("output = ",out_file_dir), connection)
            writeLines(paste0("reference = ",ref_genome), connection)

            writeLines("", connection)

            out=lapply(1:length(config),FUN=function(x){
                step=config[[x]]
                step_name=names(config[x])
                step_state=config[[x]]$run
                if(step$run){
                    if(length(step$config)>0){
                        writeLines(paste0("[",step_name,"]"),connection)
                        out=lapply(1:length(step$config),function(y){
                                config_name=names(step$config[y])
                                config_value=step$config[[y]]
                                writeLines(
                                    paste0(config_name," = ",
                                    config_value), connection
                                )
                                
                        })
                        writeLines("\n", connection)
                    }
                    
                }
                
            })

            close(connection)

        }

        build_ini_hatchet(
            config_file=.main$out_files$hatchet_config
            )
       
        .main$exec_code=paste(
            set_conda_envir(env=env_hatchet),
            " python -m hatchet run ",.main$out_files$hatchet_config
        )

        run_job(.env=.this.env)
        .main.step=.main$steps[[fn_id]]

        .env$.main <- .main   
    }
         
    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars=NULL
    )

    launch(.env=.base.env)

}

