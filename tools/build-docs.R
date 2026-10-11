# Run from the repository root after installing the current package and pkgdown.
root <- normalizePath(".", winslash = "/", mustWork = TRUE)
stage <- tempfile("emhawkes-docs-")
dir.create(stage)
files <- c("DESCRIPTION", "NAMESPACE", "README.md", "NEWS.md", "_pkgdown.yml",
           "R", "src", "man", "vignettes")
if (!all(file.copy(file.path(root, files), stage, recursive = TRUE))) {
  stop("Could not stage the documentation sources.")
}
articles <- file.path(stage, "vignettes", "articles")
dir.create(articles)
if (!file.copy(file.path(root, "vignettes_ext", "extended_example.Rmd"), articles)) {
  stop("Could not stage the extended article.")
}
override <- list(destination = file.path(root, "docs"))
if (identical(Sys.getenv("GITHUB_ACTIONS"), "true")) {
  pkgdown::build_site_github_pages(pkg = stage, dest_dir = override$destination,
                                  new_process = FALSE, install = FALSE)
} else {
  pkgdown::build_site(pkg = stage, override = override, preview = FALSE,
                      new_process = FALSE, install = FALSE)
}
