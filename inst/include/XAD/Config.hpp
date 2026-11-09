/*******************************************************************************

   Configuration settings for the build

   This file is part of XAD, a comprehensive C++ library for
   automatic differentiation.

   Copyright (C) 2010-2025 Xcelerit Computing Ltd.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Affero General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Affero General Public License for more details.

   You should have received a copy of the GNU Affero General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

******************************************************************************/

#pragma once

/* These options can be changed in client code, after XAD has already been compiled */

// Use strong inlining for higher performance - but compiles significantly slower
#ifndef XAD_USE_STRONG_INLINE
// #define XAD_USE_STRONG_INLINE
#endif

// Allow conversion operator from active type to integers, potentially missing some
// AAD variable dependency tracking
#ifndef XAD_ALLOW_INT_CONVERSION
// #define XAD_ALLOW_INT_CONVERSION
#endif


/******* The following options should not be touched after compilation of XAD */

// keep track of freed-up slots in the tape and re-use them
// #define XAD_TAPE_REUSE_SLOTS

// Disable thread-local tape usage
// #define XAD_NO_THREADLOCAL

// Reduce memory usage in the tape, at a slight performance cost
// odelia local patch: on, so every package reading these headers stores a
// recorded operation the same way -- its multiplier and its slot in two arrays,
// 12 bytes, where the paired layout pads each to 16. A sweep is bound by streaming
// that array, so it is faster this way as well as smaller. Part of the layout
// odelia's compiled tape and every consumer must agree on; see ARCHITECTURE.md.
#define XAD_REDUCED_MEMORY

// Enable JIT compilation support (record-once, compile-once, evaluate-many)
// #define XAD_ENABLE_JIT
