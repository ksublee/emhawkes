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
article_sources <- list.files(file.path(root, "vignettes_ext"),
                              pattern = "[.]Rmd$", full.names = TRUE)
if (!length(article_sources) || !all(file.copy(article_sources, articles))) {
  stop("Could not stage the website articles.")
}
override <- list(destination = file.path(root, "docs"))
if (identical(Sys.getenv("GITHUB_ACTIONS"), "true")) {
  pkgdown::build_site_github_pages(pkg = stage, dest_dir = override$destination,
                                  new_process = FALSE, install = FALSE)
} else {
  pkgdown::build_site(pkg = stage, override = override, preview = FALSE,
                      new_process = FALSE, install = FALSE)
}

# Website articles are staged under vignettes/articles for pkgdown. Point their
# source links back to the actual repository files in vignettes_ext.
for (article in basename(article_sources)) {
  page <- file.path(override$destination, "articles", sub("[.]Rmd$", ".html", article))
  html <- readLines(page, warn = FALSE, encoding = "UTF-8")
  html <- gsub("/vignettes/articles/", "/vignettes_ext/", html, fixed = TRUE)
  writeLines(html, page, useBytes = TRUE)
}
