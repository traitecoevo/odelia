testthat::test_that("odelia_load_dll returns valid library path", {
  testthat::skip_if(is_pkgload_dll(), "Skipping in pkgload load_all sessions due unstable DLL reloading semantics.")
  dll_path <- odelia_load_dll()
  expect_true(file.exists(dll_path))
  expect_equal(basename(dll_path), paste0("odelia", .Platform$dynlib.ext))
})

testthat::test_that("odelia_load_dll works when called repeatedly", {
  testthat::skip_if(is_pkgload_dll(), "Skipping in pkgload load_all sessions due unstable DLL reloading semantics.")
  path1 <- odelia_load_dll()
  path2 <- odelia_load_dll()
  expect_true(file.exists(path1))
  expect_true(file.exists(path2))
  expect_equal(normalizePath(path1), normalizePath(path2))
})
