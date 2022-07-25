#' @export

build_ui <- function() {
  shinydashboard::dashboardPage(
    skin = "black",
    build_ui_header(),
    build_ui_sidebar(),
    build_ui_body()
  )
}

#' @export

build_ui_header <- function() {
  shinydashboard::dashboardHeader(
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
  shinydashboard::dashboardSidebar(
     shiny::tags$script(htmlwidgets::JS("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")),
     width = "0px")
 }


#' @export



