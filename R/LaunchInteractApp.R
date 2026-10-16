#' Launch INTERACT Viewer Shiny Application
#'
#' @description This function launches the INTERACT Viewer Shiny dashboard.
#'
#' @param results A list containing the results of the INTERACT pipeline (which includes the precomputed data used by the Shiny app).
#' @param launch.browser Logical. If TRUE, launches the app in your default web browser (recommended for Positron). Default is TRUE.
#' @param port Integer. Optional port to listen on.
#' @param host String. The IP address to listen on (default: "127.0.0.1").
#'
#' @return This function does not return a value. It launches a Shiny application.
#'
#' @importFrom shiny shinyOptions runApp
#' @export
LaunchInteractApp = function(results, launch.browser = TRUE, port = NULL, host = "127.0.0.1") {
  if (missing(results) || is.null(results)) {
    stop("Please provide an INTERACT 'results' object to LaunchInteractApp().")
  }
  if (!is.list(results) || is.null(results$precompute)) {
    stop("The provided 'results' object does not contain required '$precompute' data.")
  }

  app.pkgs = c("shinydashboard", "shinyWidgets", "shinycssloaders", "shinyjs", "DT", "plotly", "survival", "survminer", "circlize")
  missing.pkgs = app.pkgs[!vapply(app.pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing.pkgs) > 0) {
    stop("The following packages are required to run the INTERACT Shiny app: ",
         paste(missing.pkgs, collapse = ", "),
         "\nPlease install them with: install.packages(c(",
         paste(sprintf("'%s'", missing.pkgs), collapse = ", "), "))")
  }

  appDir = system.file("app", package = "INTERACT")
  if (appDir == "") {
    appDir = "inst/app"
  }

  message("Starting INTERACT Shiny Viewer...")
  if (!is.null(results$interactions)) {
    message(sprintf("Loaded %d samples and %s interactions.",
                    length(results$precompute$sample.ids.vec),
                    format(nrow(results$interactions), big.mark = ",")))
  }

  shiny::shinyOptions(results = results)
  shiny::runApp(appDir, launch.browser = launch.browser, port = port, host = host, display.mode = "normal")
}

