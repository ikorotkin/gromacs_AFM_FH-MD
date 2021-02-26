/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014,2015,2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \libinternal \file
 * \brief
 * Include file for configuration macros from the build system.
 *
 * This header is not installed, so headers must not reference macros defined
 * here.
 *
 * \inlibraryapi
 */
#ifndef GMX_CONFIG_H
#define GMX_CONFIG_H

/* TODO: For now, disable Doxygen warnings from here */
/*! \cond */

/* Work around broken calloc() */
#define GMX_BROKEN_CALLOC 0

/* Do not optimize FFTW setups (not needed with SSE FFT kernels) */
#define GMX_DISABLE_FFTW_MEASURE 0

/* Use FFTW3 FFT library */
#define GMX_FFT_FFTW3 1

/* Use MKL FFT library */
#define GMX_FFT_MKL 0

/* Use built in fftpack FFT library */
#define GMX_FFT_FFTPACK 0

/* Target platform is x86 or x86_64 */
#define GMX_TARGET_X86 1

/* Target platform is BlueGene/Q */
#define GMX_TARGET_BGQ 0

/** Define if we are building natively on Windows */
#define GMX_NATIVE_WINDOWS 0

/** Define if we are building for Cygwin */
#define GMX_CYGWIN 0

/* GCC bug in AVX maskload/maskstore arguments - worked around internally */
#define GMX_SIMD_X86_AVX_GCC_MASKLOAD_BUG 0

/* SSE2 was selected for SIMD instruction set level */
#define GMX_SIMD_X86_SSE2 0

/* SSE4.1 was selected as SIMD instructions */
#define GMX_SIMD_X86_SSE4_1 0

/* AVX 128-bit FMA was selected as SIMD instructions */
#define GMX_SIMD_X86_AVX_128_FMA 0

/* AVX 256-bit was selected as SIMD instructions */
#define GMX_SIMD_X86_AVX_256 0

/* AVX2 256-bit SIMD instruction set level was selected */
#define GMX_SIMD_X86_AVX2_256 0

/* MIC (Xeon Phi) SIMD instruction set level was selected */
#define GMX_SIMD_X86_MIC 0

/* AVX-512F foundation level instruction SIMD */
#define GMX_SIMD_X86_AVX_512 0

/* AVX-512ER foundation level instruction SIMD */
#define GMX_SIMD_X86_AVX_512_KNL 0

/* 32-bit ARM NEON SIMD instruction set level was selected */
#define GMX_SIMD_ARM_NEON 0

/* ARM (AArch64) NEON Advanced SIMD instruction set level was selected */
#define GMX_SIMD_ARM_NEON_ASIMD 0

/* IBM QPX was selected as SIMD instructions (e.g. BlueGene/Q) */
#define GMX_SIMD_IBM_QPX 0

/* IBM VMX was selected as SIMD instructions (Power 6 and later) */
#define GMX_SIMD_IBM_VMX 0

/* IBM VSX was selected as SIMD instructions (Power 7 and later) */
#define GMX_SIMD_IBM_VSX 0

/* Fujitsu Sparc64 HPC-ACE SIMD acceleration */
#define GMX_SIMD_SPARC64_HPC_ACE 0

/* Reference SIMD implementation for testing */
#define GMX_SIMD_REFERENCE 0

/* String for SIMD instruction choice (for writing to log files and stdout) */
#define GMX_SIMD_STRING "NONE"

/* Calling convention string (if any) for routines with SIMD variable args */
#define gmx_simdcall  

/* Target mantissa accuracy for SIMD single precision math */
#define GMX_SIMD_ACCURACY_BITS_SINGLE 22

/* Target mantissa accuracy for SIMD double precision math */
#define GMX_SIMD_ACCURACY_BITS_DOUBLE 44

/* Whether a double-precision configuration may target accuracy equivalent to single precision */
#define GMX_RELAXED_DOUBLE_PRECISION 0

/* Integer byte order is big endian. */
#define GMX_INTEGER_BIG_ENDIAN 0

/* Use our own instead of system XDR libraries */
#define GMX_INTERNAL_XDR 0

/* Compile to use TNG library */
#define GMX_USE_TNG 1

/* Add support for tracing using Extrae */
#define HAVE_EXTRAE 0

/* Use MPI (with mpicc) for parallelization */
#define GMX_LIB_MPI 1

/* Use threads_mpi for parallelization */
#define GMX_THREAD_MPI 0

/* Make a parallel version of GROMACS using message passing
   (MPI or thread_mpi) */
#define GMX_MPI (GMX_LIB_MPI || GMX_THREAD_MPI)

/* MPI_IN_PLACE exists for collective operations */
#define MPI_IN_PLACE_EXISTS 1

/* Use OpenMP multithreading */
#define GMX_OPENMP 1

/* Use the Portable Hardware Locality package (hwloc) */
#define GMX_HWLOC 0

/* Can and should use nice(3) to set priority */
#define GMX_USE_NICE 1

/* Maximum number of OpenMP threads supported */
#define GMX_OPENMP_MAX_THREADS 32

/* Use if we cannot rename checkpoints */
#define GMX_NO_RENAME 0

/* Use (modified) Gamess-UK for QM-MM calculations */
#define GMX_QMMM_GAMESS 0

/* Use (modified) Gaussian0x for QM-MM calculations */
#define GMX_QMMM_GAUSSIAN 0

/* Use (modified) Mopac 7 for QM-MM calculations */
#define GMX_QMMM_MOPAC 0

/* Use ORCA for QM-MM calculations */
#define GMX_QMMM_ORCA 0

/* Use cycle counters */
#define GMX_CYCLECOUNTERS 1

/* Use sub-counters */
#define GMX_CYCLE_SUBCOUNTERS 0

/* Compile with plugin support */
#define GMX_USE_PLUGINS 1

/* Fallback path for VMD plug-ins */
#define GMX_VMD_PLUGIN_PATH "/usr/local/lib/vmd/plugins/*/molfile"

/* Define when pthreads are used */
#define THREAD_PTHREADS

/* Define when Windows threads are used */
/* #undef THREAD_WINDOWS */

/* Define native atomic operations are found */
#define TMPI_ATOMICS

/* Define for busy wait option  */
/* See gmxpre-config.h.cmakein for explanation for the #ifdef */
#ifndef TMPI_WAIT_FOR_NO_ONE
/* #undef TMPI_WAIT_FOR_NO_ONE */
#endif

/* Define for copy buffer option */
/* #undef TMPI_COPY_BUFFER */

/* Define for tmpi warnings option */
/* #undef TMPI_WARNINGS */

/* Define for profiling option */
/* #undef TMPI_PROFILE */

/* Define for Linux pthread_setaffinity_np */
#define HAVE_PTHREAD_SETAFFINITY

/* Define for X-Windows */
#define GMX_X11 0

/* Enable x86 gcc inline assembly */
#define GMX_X86_GCC_INLINE_ASM 1

/* Define constants useful for handling GPU support */
#define GMX_GPU_NONE   0
#define GMX_GPU_CUDA   1
#define GMX_GPU_OPENCL 2
/* Which kind of GPU support is configured */
#define GMX_GPU GMX_GPU_NONE

/* CUDA runtime API version (identical to CUDART_VERSION from cuda_runtime_api.h) */
/* #undef GMX_CUDA_VERSION */

/* Use a single compilation unit when compiling the CUDA (non-bonded) kernels.  */
#define GMX_CUDA_NB_SINGLE_COMPILATION_UNIT 0

/* Use NVML */
#define HAVE_NVML 0

/* Define relative path to OpenCL kernels */
#define OCL_INSTALL_DIR "share/gromacs/opencl"

/* Define to 1 if fseeko (and presumably ftello) exists and is declared. */
#define HAVE_FSEEKO 1

/* Define to 1 if _fseeki64 (and presumably _fseeki64) exists and is declared. */
#define HAVE__FSEEKI64 0

/* Have io.h (windows)*/
#define HAVE_IO_H 0

/* Define to 1 if you have the posix_memalign() function. */
#define HAVE_POSIX_MEMALIGN 1

/* Define to 1 if you have the memalign() function. */
#define HAVE_MEMALIGN 0

/* Define to 1 if you have the MSVC _aligned_malloc() function. */
#define HAVE__ALIGNED_MALLOC 0

/* Define to 1 if you have the clock_gettime() function. */
#define HAVE_CLOCK_GETTIME 1

/* Define to 1 if you have the gettimeofday() function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the rdtscp instruction. */
/* #undef HAVE_RDTSCP */

/* Define to 1 if you have the fsync() function. */
#define HAVE_FSYNC 1

/* Define to 1 if you have the Windows _commit() function. */
#define HAVE__COMMIT 0

/* Define to 1 if you have the fileno() function. */
#define HAVE_FILENO 1

/* Define to 1 if you have the _fileno() function. */
#define HAVE__FILENO 0

/* Define to 1 if you have the sigaction() function. */
#define HAVE_SIGACTION 1

/* Define for the GNU __builtin_clz() function. */
#define HAVE_BUILTIN_CLZ 1

/* Define for the GNU __builtin_clzll() function. */
#define HAVE_BUILTIN_CLZLL 1

/* Define for the MSVC _BitScanReverse() function. */
#define HAVE_BITSCANREVERSE 0

/* Define for the MSVC _BitScanReverse64() function. */
#define HAVE_BITSCANREVERSE64 0

/* Define for the IBM xlc++ __cntlz4() function. */
#define HAVE_CNTLZ4 0

/* Define for the IBM xlc++ __cntlz8() function. */
#define HAVE_CNTLZ8 0

/* Define to 1 if yo have the <unistd.h> header file. */
#define HAVE_UNISTD_H

/* Define to 1 if yo have the <pwd.h> header file. */
#define HAVE_PWD_H 1

/* Define to 1 if yo have the <dirent.h> header file. */
#define HAVE_DIRENT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H

/* Define to 1 if you have the <sched.h> header */
#define HAVE_SCHED_H

/* Define to 1 if mm_malloc.h is present, otherwise 0 */
#define HAVE_MM_MALLOC_H 1

/* Define to 1 if malloc.h is present, otherwise 0 */
#define HAVE_MALLOC_H 1

/* Define to 1 if xmmintrin.h is present, otherwise 0 */
#define HAVE_XMMINTRIN_H 1

/* Define to 1 if you have the POSIX <regex.h> header file. */
#define HAVE_POSIX_REGEX 1

/* Define to 1 if you have the C++11 <regex> header file. */
#define HAVE_CXX11_REGEX 1

/* Define to 1 if you have the sysconf() function */
#define HAVE_SYSCONF

/* Define to 1 if you have the all the affinity functions in sched.h */
#define HAVE_SCHED_AFFINITY 1

/* Define to 1 if _mm_malloc() is present in either mm_malloc.h,
 * malloc.h or xmmintrin.h, and 0 otherwise. Note that you need to
 * conditionally include the three headers too before using _mm_malloc().
 */
#define HAVE__MM_MALLOC 1

/* Define if SIGUSR1 is present */
#define HAVE_SIGUSR1 1

/* Enable gromacs quotes */
#define GMX_COOL_QUOTES 1

/* default name mangling maybe wrong on exotic plattforms */
#define F77_FUNC(name,NAME) name ## _

/* Define if we have pipes */
#define HAVE_PIPES 1

/* Define if we have feenableexcept */
#define HAVE_FEENABLEEXCEPT 1

/* Define if we have zlib */
#define HAVE_ZLIB 0

/*! \endcond */

#endif
