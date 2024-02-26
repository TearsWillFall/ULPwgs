
#' Build dataframe with matched instrument_id patterns and instrument names
#' 
#'
#' @param instruments List with flowcell_id, flowcell_type and instrument information
#' @export

build_instrument_id = function(instruments=list(
  instrument=c("MiSeq",
  "Genome Analyzer IIx",
  "MiSeq",
  "HiSeq 1500",
  "HiSeq 1500",
  "HiSeq 2500",
  "HiSeq 2500",
  "HiSeq 3000",
  "HiSeq 3000",
  "HiSeq 4000",
  "HiSeq X",
  "NextSeq",
  "NextSeq",
  "MiniSeq",
  "NovaSeq 6000"
),
  pattern=c(
    "HWI-M[0-9]{4}$",
    "HWUSI",
    "M[0-9]{5}$",
    "HWI-C[0-9]{5}$",
    "C[0-9]{5}$",
    "HWI-D[0-9]{5}$",
    "D[0-9]{5}$",
    "J[0-9]{5}$",
    "K[0-9]{5}$",
    "K[0-9]{5}$",
    "E[0-9]{5}$",
    "NB[0-9]{6}$",
    "NS[0-9]{6}$",
    "MN[0-9]{5}$",
    "A[0-9]{5}$"
  )
)){
    data.frame(instruments,stringsAsFactors = FALSE)
  }





#' Build dataframe with matched flowcell_id patterns, flowcell_type and instrument
#' 
#'
#' @param flowcells List with flowcell_id, flowcell_type and instrument information
#' @export

build_flowcell_id=function(flowcells=list(
    instrument=c(
      "HiSeq 1500",
      "HiSeq 2000",
      "HiSeq 2500",
      "HiSeq 1000",
      "HiSeq 1500",
      "HiSeq 2000",
      "HiSeq 2500",
      "HiSeq 1500",
      "HiSeq 2500",
      "HiSeq 1500",
      "HiSeq 2500",
      "HiSeq 1500",
      "HiSeq 2500",
      "HiSeq 4000",
      "HiSeq 4000",
      "HiSeq X",
      "HiSeq X",
      "HiSeq X",
      "NextSeq",
      "NextSeq",
      "NextSeq",
      "NextSeq",
      "MiSeq",
      "MiSeq",
      "MiSeq",
      "MiSeq",
      "NovaSeq 6000",
      "NovaSeq 6000",
      "NovaSeq 6000",
      "NovaSeq 6000",
      "NovaSeq 6000",
      "NovaSeq 6000",
      "NovaSeq 6000",
      "NovaSeq 6000"

    ),flowcell=c(
      "MiSeq flow cell",
      "High Output (8-lane) v4 flow cell",
      "High Output (8-lane) v4 flow cell",
      "High Output (8-lane) v4 flow cell",
      "High Output (8-lane) v3 flow cell",
      "High Output (8-lane) v3 flow cell",
      "High Output (8-lane) v3 flow cell",
      "High Output (8-lane) v3 flow cell",
      "Rapid Run (2-lane) v1 flow cell",
      "Rapid Run (2-lane) v1 flow cell",
      "Rapid Run (2-lane) v2 flow cell",
      "Rapid Run (2-lane) v2 flow cell",
      "Rapid Run (2-lane) v2 flow cell",
      "Rapid Run (2-lane) v2 flow cell",
      "(8-lane) v1 flow cell",
      "(8-lane) v1 flow cell",
      "(8-lane) flow cell",
      "(8-lane) flow cell",
      "(8-lane) flow cell",
        "High output flow cell",
        "High output flow cell",
        "High output flow cell",
        "Mid output flow cell",
        "MiSeq flow cell",
        "MiSeq flow cell",
        "MiSeq nano flow cell",
        "MiSeq micro flow cell",
        "S2 flow cell",
        "S2 flow cell",
        "S4 flow cell",
        "S4 flow cell",
        "S4 flow cell",
        "SP flow cell",
        "SP flow cell",
        "SP flow cell"

        


    ),pattern=c(
      "000000000.[A-Z][A-Z][A-Z][A-Z]$",
      "C[A-Z,0-9]{4}ANXX$",
      "C[A-Z,0-9]{4}ANXX$",
      "C[A-Z,0-9]{4}ANXX$",
      "C[A-Z,0-9]{4}ACXX$",
      "C[A-Z,0-9]{4}ACXX$",
      "C[A-Z,0-9]{4}ACXX$",
      "C[A-Z,0-9]{4}ACXX$",
      "H[A-Z,0-9]{4}ADXX$",
      "H[A-Z,0-9]{4}ADXX$",
      "H[A-Z,0-9]{4}BCXX$",
      "H[A-Z,0-9]{4}BCXX$",
      "H[A-Z,0-9]{4}BCXY$",
      "H[A-Z,0-9]{4}BCXY$",
      "H[A-Z,0-9]{4}BBXX$",
      "H[A-Z,0-9]{4}BBXY$",
      "H[A-Z,0-9]{4}CCXX$",
      "H[A-Z,0-9]{4}CCXY$",
      "H[A-Z,0-9]{4}ALXX$",
      "H[A-Z,0-9]{4}BGXX$",
      "H[A-Z,0-9]{4}BGXY$",
      "H[A-Z,0-9]{4}BGX2$",
      "H[A-Z,0-9]{4}AFXX$",
      "A[A-Z,0-9]{4}[A-Z][A-Z][A-Z][A-Z]$",
      "B[A-Z,0-9]{4}[A-Z][A-Z][A-Z][A-Z]$",
      "D[A-Z,0-9]{4}[A-Z][A-Z][A-Z][A-Z]$",
      "G[A-Z,0-9]{4}[A-Z][A-Z][A-Z][A-Z]$",
      "H[A-Z,0-9]{4}DMXX$",
      "H[A-Z,0-9]{4}DMXY$",
      "H[A-Z,0-9]{4}DSXX$",
      "H[A-Z,0-9]{4}DSX[1-9]$",
      "H[A-Z,0-9]{4}DSXY$",
      "H[A-Z,0-9]{4}DRXX$",
      "H[A-Z,0-9]{4}DRX2$",
      "H[A-Z,0-9]{4}DRXY$"
    )
  )
  ) {
      data.frame(flowcells,stringsAsFactors = FALSE)
}
