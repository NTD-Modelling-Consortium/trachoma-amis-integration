version <- paste(version[c("major", "minor")], collapse = ".")
underscored_version <- chartr(".", "_", version)
file_name <- paste0("renv_", underscored_version, ".lock")
if (file.exists(file_name)) {
    if (file.exists("renv.lock")) {
        file.remove("renv.lock")
    }
    file.copy(file_name, "renv.lock")
    renv::restore()
} else {
    print(paste("No renv lock file for R", version))
}
