cmake_minimum_required(VERSION 3.20)
project(FastGAJointIndex VERSION 0.1.0 LANGUAGES C CXX)

if(NOT UNIX)
  message(FATAL_ERROR "This implementation requires Linux or macOS (POSIX I/O and pthreads).")
endif()
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_EXTENSIONS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_compile_options(-UNDEBUG)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
set(OPENMP_DEFAULT OFF)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(OPENMP_DEFAULT ON)
endif()
option(SJ_ENABLE_OPENMP "Use OpenMP for the libsais suffix-array sorters of the PFP tools" ${OPENMP_DEFAULT})
option(SJ_SANITIZE_ADAPTER "Instrument native frontends, adapter/driver, provenance writer and codec test with ASAN/UBSAN" OFF)
if(SJ_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS C CXX)
endif()

# Never modify vendored upstream sources. Verify and patch build-local copies.
include(cmake/PrepareVendor.cmake)
set(FG "${CMAKE_BINARY_DIR}/upstream/fastga")
set(BWT "${CMAKE_BINARY_DIR}/upstream/bigbwt")
set(SRC "${PROJECT_SOURCE_DIR}/src")

function(sj_native target)
  target_include_directories(${target} PRIVATE "${FG}" "${SRC}")
  # Keep upstream assertions enabled, including in Release builds.
  target_compile_options(${target} PRIVATE -fno-strict-aliasing -UNDEBUG)
endfunction()
# libdeflate writes and reads the zlib-format block streams (pinned in vendor/).
set(DEFLATE "${PROJECT_SOURCE_DIR}/vendor/libdeflate")
add_library(sj_deflate STATIC "${DEFLATE}/lib/utils.c" "${DEFLATE}/lib/x86/cpu_features.c" "${DEFLATE}/lib/arm/cpu_features.c"
  "${DEFLATE}/lib/deflate_compress.c" "${DEFLATE}/lib/deflate_decompress.c" "${DEFLATE}/lib/zlib_compress.c"
  "${DEFLATE}/lib/zlib_decompress.c" "${DEFLATE}/lib/adler32.c" "${DEFLATE}/lib/crc32.c")
target_include_directories(sj_deflate PUBLIC "${DEFLATE}")
set_target_properties(sj_deflate PROPERTIES POSITION_INDEPENDENT_CODE ON)
function(sj_program target)
  sj_native(${target})
  target_link_libraries(${target} PRIVATE sj_deflate Threads::Threads ZLIB::ZLIB m)
endfunction()
add_library(gene OBJECT "${FG}/gene_core.c" "${FG}/ONElib.c" "${FG}/GDB.c")
add_library(alignment OBJECT "${FG}/align.c" "${FG}/alncode.c")
add_library(annotations OBJECT "${FG}/ANO.c")
add_library(rsd OBJECT "${FG}/RSDsort.c")
add_library(msd OBJECT "${FG}/MSDsort.c")
target_compile_definitions(msd PRIVATE LCPs)
add_library(fastk_native OBJECT "${FG}/libfastk.c")
add_library(fastk_shared OBJECT "${FG}/libfastk.c")
target_compile_definitions(fastk_shared PRIVATE SXGC_SHARED_GIX)
foreach(target gene alignment annotations rsd msd fastk_native fastk_shared)
  sj_native(${target})
endforeach()
add_executable(FastGA-baseline "${FG}/FastGA.c" $<TARGET_OBJECTS:fastk_native>
  $<TARGET_OBJECTS:gene> $<TARGET_OBJECTS:alignment> $<TARGET_OBJECTS:rsd>)
add_executable(FastGA-shared-core "${FG}/FastGA.c" "${SRC}/shared_gix.cpp" $<TARGET_OBJECTS:fastk_shared>
  $<TARGET_OBJECTS:gene> $<TARGET_OBJECTS:alignment> $<TARGET_OBJECTS:rsd>)
target_compile_definitions(FastGA-shared-core PRIVATE SXGC_SHARED_GIX)
add_executable(FastGA-shared "${SRC}/fastga_frontend.cpp")
add_executable(GIXmake-shared "${SRC}/fasta_index.cpp")
add_executable(build-manifest "${PROJECT_SOURCE_DIR}/tools/build_manifest.cpp")
target_include_directories(build-manifest PRIVATE "${SRC}")
if(NOT APPLE)
  find_package(OpenSSL 1.1 REQUIRED COMPONENTS Crypto)
endif()
foreach(target FastGA-shared GIXmake-shared build-manifest)
  target_include_directories(${target} PRIVATE "${PROJECT_SOURCE_DIR}/vendor/nlohmann-json/include")
  target_link_libraries(${target} PRIVATE Threads::Threads)
  if(NOT APPLE)
    target_link_libraries(${target} PRIVATE OpenSSL::Crypto)
  endif()
endforeach()
# Remove the obsolete generated Python runtime from earlier builds.
file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/fastga-runtime")
add_executable(GIXmake "${FG}/GIXmake.c" $<TARGET_OBJECTS:msd> $<TARGET_OBJECTS:annotations>
  $<TARGET_OBJECTS:gene> $<TARGET_OBJECTS:fastk_native>)
target_compile_definitions(GIXmake PRIVATE LCPs)
add_executable(GIXmake-direct "${FG}/GIXmake.c" "${SRC}/direct_stream.c" "${SRC}/joint_filter.cpp"
  $<TARGET_OBJECTS:msd> $<TARGET_OBJECTS:annotations> $<TARGET_OBJECTS:gene> $<TARGET_OBJECTS:fastk_native>)
target_compile_definitions(GIXmake-direct PRIVATE LCPs SXGC_DIRECT_STREAM)
add_executable(GIXmake-marks "${FG}/GIXmake.c" "${SRC}/pfp_marks.cpp" $<TARGET_OBJECTS:msd>
  $<TARGET_OBJECTS:annotations> $<TARGET_OBJECTS:gene> $<TARGET_OBJECTS:fastk_native>)
target_compile_definitions(GIXmake-marks PRIVATE LCPs SXGC_PFP_MARKS)
add_executable(FAtoGDB "${FG}/FAtoGDB.c" $<TARGET_OBJECTS:gene> $<TARGET_OBJECTS:annotations>)
add_executable(ALNtoPAF "${FG}/ALNtoPAF.c" $<TARGET_OBJECTS:gene> $<TARGET_OBJECTS:alignment>)
add_executable(ALNtoPSL "${FG}/ALNtoPSL.c" $<TARGET_OBJECTS:gene> $<TARGET_OBJECTS:alignment>)
add_executable(gixshare "${SRC}/shared_gix.cpp" $<TARGET_OBJECTS:fastk_native> $<TARGET_OBJECTS:gene>)
target_compile_definitions(gixshare PRIVATE SXGC_PACK_MAIN)
add_executable(pfp-emit "${SRC}/pfp_emit.cpp" $<TARGET_OBJECTS:fastk_native> $<TARGET_OBJECTS:gene>)
foreach(target FastGA-baseline FastGA-shared-core GIXmake GIXmake-direct GIXmake-marks FAtoGDB ALNtoPAF ALNtoPSL gixshare pfp-emit)
  sj_program(${target})
endforeach()
add_library(joint_sais OBJECT "${PROJECT_SOURCE_DIR}/vendor/libsais/src/libsais.c"
  "${PROJECT_SOURCE_DIR}/vendor/libsais/src/libsais64.c")
target_include_directories(joint_sais PRIVATE "${PROJECT_SOURCE_DIR}/vendor/libsais/include")
if(SJ_ENABLE_OPENMP)
  target_compile_definitions(joint_sais PRIVATE LIBSAIS_OPENMP)
  target_link_libraries(joint_sais PRIVATE OpenMP::OpenMP_C)
endif()
# Anchor-based copy finder: no suffix array. libsais now serves only the Big-BWT parse/dictionary sorters (pfp-tools).
add_executable(joint-copy "${SRC}/joint_copy.cpp")
sj_program(joint-copy)

# The PFP tools remain separate processes, as in the validated implementation.
add_library(bwt_utils OBJECT "${BWT}/utils.c")
add_library(bwt_gsa OBJECT "${BWT}/gsa/gsacak.c")
add_library(bwt_xerrors OBJECT "${BWT}/xerrors.c")
foreach(target bwt_utils bwt_gsa bwt_xerrors)
  target_compile_definitions(${target} PRIVATE M64)
  target_include_directories(${target} PRIVATE "${BWT}")
endforeach()
add_executable(newscan "${BWT}/newscan.cpp" $<TARGET_OBJECTS:bwt_utils>)
target_compile_definitions(newscan PRIVATE NOTHREADS)
# pscan is the multithreaded parser used by the builder; newscan stays as the serial reference.
add_executable(pscan "${BWT}/pscan.cpp" $<TARGET_OBJECTS:bwt_utils> $<TARGET_OBJECTS:bwt_xerrors>)
add_executable(bwtparse "${BWT}/bwtparse.c" $<TARGET_OBJECTS:bwt_utils> $<TARGET_OBJECTS:bwt_gsa> $<TARGET_OBJECTS:joint_sais>)
target_compile_definitions(bwtparse PRIVATE M64 SJ_LIBSAIS)
target_include_directories(bwtparse PRIVATE "${PROJECT_SOURCE_DIR}/vendor/libsais/include")
if(SJ_ENABLE_OPENMP)
  # bwtparse: libsais_int_omp suffix sorting plus OpenMP chunk/scatter loops
  target_compile_definitions(bwtparse PRIVATE SJ_LIBSAIS_OMP)
  target_link_libraries(bwtparse PRIVATE OpenMP::OpenMP_C)
endif()
# pfbwt builds the dictionary SA/LCP with libsais (parallel generalized SA);
# gSACA-K stays linked for -G (legacy serial path) and -C (equality check).
add_executable(pfbwt "${BWT}/pfbwt.cpp" $<TARGET_OBJECTS:bwt_utils> $<TARGET_OBJECTS:bwt_gsa> $<TARGET_OBJECTS:bwt_xerrors>
  $<TARGET_OBJECTS:joint_sais>)
target_compile_definitions(pfbwt PRIVATE M64)
target_include_directories(pfbwt PRIVATE "${PROJECT_SOURCE_DIR}/vendor/libsais/include" "${SRC}/pfbwt")
# dict_sais.hpp widens 32-bit libsais results in place inside the 64-bit arrays.
target_compile_options(pfbwt PRIVATE -fno-strict-aliasing)
if(SJ_ENABLE_OPENMP)
  target_compile_definitions(pfbwt PRIVATE LIBSAIS_OPENMP)
  target_link_libraries(pfbwt PRIVATE OpenMP::OpenMP_CXX)
endif()
foreach(target newscan pscan bwtparse pfbwt)
  target_include_directories(${target} PRIVATE "${BWT}")
  target_link_libraries(${target} PRIVATE Threads::Threads m)
  set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pfp-tools")
endforeach()

include(CTest)
option(SJ_PYTHON_TESTS "Run complete Python-driven correctness suites (requires Python 3.11+)" OFF)
if(SJ_PYTHON_TESTS)
  if(NOT BUILD_TESTING)
    message(FATAL_ERROR "SJ_PYTHON_TESTS requires BUILD_TESTING=ON")
  endif()
  find_package(Python3 3.11 REQUIRED COMPONENTS Interpreter)
endif()
if(BUILD_TESTING)
  add_executable(joint-codec-test "${PROJECT_SOURCE_DIR}/tests/test_codec.cpp"
    $<TARGET_OBJECTS:fastk_native> $<TARGET_OBJECTS:gene>)
  sj_program(joint-codec-test)
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/test-results/tmp")
  add_test(NAME codec COMMAND joint-codec-test "${CMAKE_BINARY_DIR}/test-results/weight-frame.bin")
  # libsais dictionary SA/LCP of pfbwt against gSACA-K, on synthetic dictionaries.
  add_executable(pfbwt-dict-test "${PROJECT_SOURCE_DIR}/tests/test_pfbwt_dict.cpp"
    $<TARGET_OBJECTS:bwt_gsa> $<TARGET_OBJECTS:joint_sais>)
  target_compile_definitions(pfbwt-dict-test PRIVATE M64)
  target_compile_options(pfbwt-dict-test PRIVATE -fno-strict-aliasing)
  target_include_directories(pfbwt-dict-test PRIVATE "${BWT}" "${SRC}/pfbwt" "${PROJECT_SOURCE_DIR}/vendor/libsais/include")
  if(SJ_ENABLE_OPENMP)
    target_compile_definitions(pfbwt-dict-test PRIVATE LIBSAIS_OPENMP)
    target_link_libraries(pfbwt-dict-test PRIVATE OpenMP::OpenMP_CXX)
  endif()
  add_test(NAME pfbwt-dict COMMAND pfbwt-dict-test)
  set_tests_properties(codec PROPERTIES ENVIRONMENT "SXGC_VERIFY_TMPDIR=${CMAKE_BINARY_DIR}/test-results/tmp;ASAN_OPTIONS=detect_leaks=0;UBSAN_OPTIONS=halt_on_error=1")
endif()
if(SJ_PYTHON_TESTS)
  add_test(NAME build-provenance COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/test_build_provenance.py"
    --bin "${CMAKE_BINARY_DIR}" --report "${CMAKE_BINARY_DIR}/test-results/build-provenance.json")
  add_test(NAME pfp-parse-equivalence COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/test_pfp_parse.py"
    --pfp-tools "${CMAKE_BINARY_DIR}/pfp-tools" --report "${CMAKE_BINARY_DIR}/test-results/pfp-parse.json")
  add_test(NAME pfp-sort-equivalence COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/test_pfp_sort.py"
    --pfp-tools "${CMAKE_BINARY_DIR}/pfp-tools" --report "${CMAKE_BINARY_DIR}/test-results/pfp-sort.json")
  add_test(NAME thread-budget COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/test_thread_budget.py"
    --bin "${CMAKE_BINARY_DIR}" --report "${CMAKE_BINARY_DIR}/test-results/thread-budget.json")
  add_test(NAME joint-fused-equivalence COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/test_joint_fused.py"
    --bin "${CMAKE_BINARY_DIR}" --pfp-tools "${CMAKE_BINARY_DIR}/pfp-tools" --report "${CMAKE_BINARY_DIR}/test-results/joint-fused.json")
  add_test(NAME stream-ranges COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/test_stream_ranges.py"
    --bin "${CMAKE_BINARY_DIR}" --pfp-tools "${CMAKE_BINARY_DIR}/pfp-tools" --report "${CMAKE_BINARY_DIR}/test-results/stream-ranges.json")
  add_test(NAME merge-fasta COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/test_merge_fasta.py")
  add_test(NAME joint-fixtures COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/test_joint.py"
    --bin "${CMAKE_BINARY_DIR}" --pfp-tools "${CMAKE_BINARY_DIR}/pfp-tools" --threads 4
    --report "${CMAKE_BINARY_DIR}/test-results/joint.json")
  add_test(NAME joint-copy-finder COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/test_joint_copy.py"
    --bin "${CMAKE_BINARY_DIR}" --report "${CMAKE_BINARY_DIR}/test-results/joint-copy.json")
  add_test(NAME fasta-entry-point COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/test_fasta_cli.py"
    --bin "${CMAKE_BINARY_DIR}" --report "${CMAKE_BINARY_DIR}/test-results/fasta-cli.json")
  set_tests_properties(fasta-entry-point PROPERTIES
    ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0;UBSAN_OPTIONS=halt_on_error=1")
  add_test(NAME native-index-cli COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/test_index_cli.py"
    --bin "${CMAKE_BINARY_DIR}" --report "${CMAKE_BINARY_DIR}/test-results/index-cli.json")
  set_tests_properties(native-index-cli PROPERTIES
    ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0;UBSAN_OPTIONS=halt_on_error=1")
  add_test(NAME construction-comparison COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/test_construction.py"
    --bin "${CMAKE_BINARY_DIR}" --report "${CMAKE_BINARY_DIR}/test-results/construction.json")
endif()
if(SJ_SANITIZE_ADAPTER)
  foreach(target FastGA-shared GIXmake-shared FastGA-shared-core gixshare joint-codec-test build-manifest)
    if(TARGET ${target})
      target_compile_options(${target} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
      target_link_options(${target} PRIVATE -fsanitize=address,undefined)
    endif()
  endforeach()
  if(SJ_PYTHON_TESTS)
    add_test(NAME sanitized-runtime COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/test_runtime.py"
      --bin "${CMAKE_BINARY_DIR}" --instrumented "$<TARGET_FILE:FastGA-shared-core>"
      --report "${CMAKE_BINARY_DIR}/test-results/runtime.json")
  endif()
endif()
add_custom_target(provenance ALL
  COMMAND build-manifest "${PROJECT_SOURCE_DIR}" "${CMAKE_BINARY_DIR}"
    "${CMAKE_C_COMPILER}" "${CMAKE_CXX_COMPILER}" "$<BOOL:${SJ_ENABLE_OPENMP}>"
  DEPENDS build-manifest FastGA-baseline FastGA-shared GIXmake-shared FastGA-shared-core GIXmake GIXmake-direct GIXmake-marks
    FAtoGDB ALNtoPAF ALNtoPSL gixshare pfp-emit joint-copy newscan pscan bwtparse pfbwt
  VERBATIM)
