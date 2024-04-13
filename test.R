
devtools::document()
tmp=cp_data(
    origin=c("/lustre/home/regmova/pip-22.1.2.tar.gz","/lustre/home/regmova/GR_SIGNAL.txt","/lustre/home/regmova/output.txt","/lustre/home/regmova/PAA_stdout.log"),
    output_dir="/mnt/gpfs/live/rd01__/ritd-ag-project-rd015z-ammha58/mnt",
    threads=1,verbose=TRUE,overwrite=TRUE,remote="output_dir"
)
