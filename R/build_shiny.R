#' @export

build_ui <- function() {
  dashboardPage(
    skin = "black",
    build_ui_header(),
    build_ui_sidebar(),
    build_ui_body()
  )
}

#' @export

build_ui_header <- function() {
  dashboardHeader(
    titleWidth=0
    
  )
}
#' @export


build_ui_body <- function() {
  shinydashboard::dashboardBody(
    uiOutput("UI")
  )
}

#' @export

build_ui_sidebar <- function() {
  dashboardSidebar(
     tags$script(JS("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")),
     width = "0px")
 }


#' @export



