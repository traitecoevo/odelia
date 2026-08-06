# Optional operator splitting on the drainage demonstrator. The fast layers drain
# by a stiff power law with a closed-form recession; the same model is run with the
# adaptive sub-cycle (unsplit) and with the exact-flow + ROS34PW2 split, so the
# cost of resolving the drainage explicitly can be read off directly.

testthat::test_that("split matches unsplit at the accuracy the macro step delivers", {
  L <- 5; M <- 10; tt <- seq(0, 20, by = 1.0); tol <- 1e-5
  us <- odelia:::drainage_mri(1000, L, M, FALSE, "erk33a", tt, tol)
  sp <- odelia:::drainage_mri(1000, L, M, TRUE,  "erk33a", tt, tol)
  ref <- odelia:::drainage_reference(1000, L, M, tt, 1e-10)$states
  # both resolve the same macro scheme: they agree within the macro error
  expect_lt(max(abs(sp$states - us$states)), 1e-3)
  # and split is no less accurate than unsplit against the tight reference
  expect_lt(max(abs(sp$states - ref)), 1.2 * max(abs(us$states - ref)))
})

testthat::test_that("splitting cuts fast-step cost and the win is robust to stiffness", {
  L <- 5; M <- 10; tt <- seq(0, 20, by = 1.0); tol <- 1e-5
  speedup <- function(c) {
    us <- odelia:::drainage_mri(c, L, M, FALSE, "erk33a", tt, tol)$n_fast
    sp <- odelia:::drainage_mri(c, L, M, TRUE,  "erk33a", tt, tol)$n_fast
    us / sp
  }
  s_lo <- speedup(1); s_hi <- speedup(1000)
  expect_gt(s_lo, 3)                          # exact drainage flow is much cheaper
  expect_gt(s_hi, 3)
  expect_lt(abs(s_hi - s_lo) / s_lo, 0.2)     # and the win barely moves with stiffness
})

testthat::test_that("the exact drainage recession preserves positivity", {
  L <- 5; M <- 8; tt <- seq(0, 20, by = 1.0)
  sp <- odelia:::drainage_mri(1000, L, M, TRUE, "erk33a", tt, 1e-5)
  expect_true(all(sp$states[, (M + 1):(M + L)] > 0))   # soil (fast) block is the tail
})

testthat::test_that("reverse mode through the split inner matches finite difference", {
  # record->replay through Strang(exact flow, ROS34PW2, exact flow): the ROS
  # Jacobian is passive, the exact drainage recession and stages tape. Checked
  # against a frozen-schedule central difference over [c, initial state].
  # The match is exact here because the residual r(g-u) is linear in u, so the
  # passive Jacobian is the true one; a residual nonlinear in u would differ from
  # this finite difference by the omitted d(Jacobian) terms (a W-method approximation).
  L <- 5; M <- 8; tt <- seq(0, 20, by = 1.0)
  for (c in c(10, 100, 1000)) {
    g <- odelia:::drainage_gradient_split(c, L, M, "erk33a", tt, 1e-6, 1e-6)
    expect_length(g$grad_adjoint, 1 + L + M)          # c + initial [u; x]
    expect_lt(max(abs(g$grad_adjoint - g$grad_fd)), 1e-6)
    expect_gt(abs(g$grad_adjoint[1]), 1e-6)           # a non-trivial dJ/dc through the flow
  }
})
