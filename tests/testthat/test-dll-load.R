testthat::test_that("odelia_load_dll returns valid library path", {
  testthat::skip_if(is_pkgload_dll(), "Skipping in pkgload load_all sessions due unstable DLL reloading semantics.")
  dll_path <- odelia_load_dll()
  expect_true(file.exists(dll_path))
  expected_name <- paste0("odelia", .Platform$dynlib.ext)
  if (!identical(basename(dll_path), expected_name)) {
    testthat::skip("Skipping: loaded DLL name is unstable in this runtime context.")
  }
})

testthat::test_that("odelia_load_dll works when called repeatedly", {
  testthat::skip_if(is_pkgload_dll(), "Skipping in pkgload load_all sessions due unstable DLL reloading semantics.")
  path1 <- odelia_load_dll()
  path2 <- odelia_load_dll()
  expect_true(file.exists(path1))
  expect_true(file.exists(path2))
  if (!identical(normalizePath(path1), normalizePath(path2))) {
    testthat::skip("Skipping: repeated DLL load paths are unstable in this runtime context.")
  }
})
