#' Preprocess sequencing data
#' 
#'
#' @param bam Path to a single or multiple BAM files
#' @param chromosomes Select chromosomes to analyze. Default NULL
#' @param env_hatchet Hatchet conda enviroment
#' @param config Hatchet configuration for each step
#' @export

run_hatchet=function(
    bam=NULL,
    chromosomes=NULL,
    env_hatchet=build_default_python_enviroment_list()$env_hatchet,
    config=build_default_hatchet_config(),
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
            config_file=NULL,
            threads=NULL,
            chromosomes=NULL,
            config=NULL
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
            writeLines(paste0("bam = ",paste0(bam,collapse="\t")), connection)
            writeLines(paste0("sample = ",paste0(Vectorize(get_file_name)(bam),collapse="\t")), connection)

            writeLines("", connection)

            out=lapply(1:length(config),FUN=function(x){
                step=config[[x]]
                step_name=names(config[x])
                step_state=config[[x]]$run
                if(step$run){
                    if(length(step$config)>0){
                        writeLines(paste0("[",names(step_name),"]"),connection)
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
            config_file=.main$out_files$hatchet_config,
            threads=threads,
            chromosomes=chromosomes,
            config=config
            )
       
        .main$exec_code=paste(
            set_conda_envir(env=env_hatchet),
            " python -m hatchet ",.main$out_files$hatchet_config
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

