


#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# New batch/sweep version
#
# Usage:
#   ./run_sweep_new.sh
#
# It loops over:
#   FORCING_ROOT/scale_*/*.nc
#
# and runs each experiment with:
#   beta s_tau damp1_days u_thres flux_extra alpha const1 ppow aexp tag forcing_file out_dir
# =========================================================

# =========================================================
# 0. paths / params
# =========================================================

SRC="${REPO_ROOT:-$(dirname "$0")/../..}/src/fortran/1d_model.f90"
EXE="./d.out"

FORCING_ROOT="${FORCING_ROOT:-./data/forcing_sweep}"
OUTPUT_ROOT="${OUTPUT_ROOT:-./output/forcing_sweep}"

# -------- new-version parameters --------
beta="0.233"
s_tau="0.5"
damp1_days="5"
u_thres="7"
flux_extra="16"

# -------- newly added parameters --------
alpha="0.5"
const1="2.2e-3"
ppow="4.0"
aexp="0.3"

# =========================================================
# 1. load modules
# =========================================================

echo "[INFO] Trying GCC + generic NetCDF first..."
module purge

# Load HDF5 first
(module load hdf5/1.12.2 || module load hdf5/1.10.9) >/dev/null 2>&1 || true

NETCDF_LOADED=0

for m in \
  netcdf-fortran \
  netcdf-c \
  netcdf \
  netcdf/4.9.0 \
  netcdf/4.9.0+intel-2020
do
  if module load "$m" >/dev/null 2>&1; then
    echo "[INFO] loaded module: $m"
    NETCDF_LOADED=1
    break
  fi
done

if [[ $NETCDF_LOADED -eq 0 ]]; then
  while IFS= read -r m; do
    [[ -n "$m" ]] || continue
    if module load "$m" >/dev/null 2>&1; then
      echo "[INFO] loaded module: $m"
      NETCDF_LOADED=1
      break
    fi
  done < <(module -t avail 2>&1 | grep -E '^(netcdf(-c|-fortran)?)(/|$)' | head -n 10)
fi

[[ $NETCDF_LOADED -eq 1 ]] || {
  echo "[ERR] failed to load any NetCDF module"
  exit 2
}

echo "[INFO] which gfortran : $(command -v gfortran || true)"
echo "[INFO] which ifort    : $(command -v ifort || true)"
echo "[INFO] which nc-config: $(command -v nc-config || true)"
echo "[INFO] which nf-config: $(command -v nf-config || true)"
echo "[INFO] nc-config --prefix: $(nc-config --prefix 2>/dev/null || true)"

# =========================================================
# 2. decide whether we need to build
# =========================================================

need_build=0

if [[ ! -x "$EXE" ]]; then
  need_build=1
elif [[ "$SRC" -nt "$EXE" ]]; then
  need_build=1
fi

# =========================================================
# 3. compile only when needed
# =========================================================

if [[ $need_build -eq 1 ]]; then
  echo "[INFO] building $EXE"

  NETCDF_PREFIX="$(nc-config --prefix 2>/dev/null || true)"
  USE_GFORTRAN=1

  if [[ -n "${NETCDF_PREFIX}" && "${NETCDF_PREFIX}" == *"+intel"* ]]; then
    echo "[INFO] Detected Intel-built NetCDF (${NETCDF_PREFIX}); will use Intel Fortran compiler."
    USE_GFORTRAN=0
  fi

  if command -v gfortran >/dev/null 2>&1 && [[ $USE_GFORTRAN -eq 1 ]]; then
    FC=gfortran
    FFLAGS="-O3 -fPIC -ffree-line-length-none"

    echo "[INFO] Using $FC"

    if command -v nf-config >/dev/null 2>&1; then
      INCFLAGS="$(nf-config --fflags)"
      LIBFLAGS="$(nf-config --flibs)"
      LIBDIR="$(nf-config --prefix 2>/dev/null || true)/lib"
    else
      command -v nc-config >/dev/null 2>&1 || {
        echo "[ERR] nc-config missing"
        exit 2
      }

      INC_C="$(nc-config --includedir)"
      LIBS_C="$(nc-config --libs)"
      PREFIX="$(nc-config --prefix)"

      INC_F="$(dirname "$(find "$PREFIX" -name netcdf.mod 2>/dev/null | head -n1 || true)")"

      [[ -n "$INC_F" && -n "$(find "$PREFIX" -name 'libnetcdff.*' 2>/dev/null | head -n1)" ]] || {
        echo "[ERR] NetCDF-Fortran headers/libs not found under $PREFIX"
        exit 3
      }

      INCFLAGS="-I${INC_C} -I${INC_F}"
      LIBFLAGS="-lnetcdff ${LIBS_C}"
      LIBDIR="$(dirname "$(find "$PREFIX" -name 'libnetcdff.*' 2>/dev/null | head -n1)")"
    fi

    echo "[INFO] INCFLAGS=$INCFLAGS"
    echo "[INFO] LIBFLAGS=$LIBFLAGS"

    "$FC" $FFLAGS $INCFLAGS "$SRC" $LIBFLAGS -Wl,-rpath,"$LIBDIR" -o "$EXE" 2>&1 | tee compile.log

  else
    echo "[INFO] Using Intel stack."

    module purge
    module load netcdf

    FC=ifort
    FFLAGS="-O3 -fp-model precise -mcmodel=medium -shared-intel"

    command -v nc-config >/dev/null 2>&1 || {
      echo "[ERR] nc-config missing"
      exit 2
    }

    INC_C="$(nc-config --includedir)"
    LIBS_C="$(nc-config --libs)"
    PREFIX="$(nc-config --prefix)"
    INC_F="$(dirname "$(find "$PREFIX" -name netcdf.mod 2>/dev/null | head -n1 || true)")"
    LIBDIR="$(dirname "$(find "$PREFIX" -name 'libnetcdff.*' 2>/dev/null | head -n1)")"

    [[ -n "$INC_F" && -n "$LIBDIR" ]] || {
      echo "[ERR] Intel NetCDF-Fortran bits missing under $PREFIX"
      exit 3
    }

    INCFLAGS="-I${INC_C} -I${INC_F}"
    LIBFLAGS="-lnetcdff ${LIBS_C}"

    echo "[INFO] INCFLAGS=$INCFLAGS"
    echo "[INFO] LIBFLAGS=$LIBFLAGS"
    echo "[INFO] Compiling $SRC with $FC"

    "$FC" $FFLAGS $INCFLAGS "$SRC" $LIBFLAGS -Wl,-rpath,"$LIBDIR" -o "$EXE" 2>&1 | tee compile.log
  fi
else
  echo "[INFO] existing executable is up to date: $EXE"
fi

[[ -x "$EXE" ]] || {
  echo "[ERR] compile failed: $EXE not found"
  exit 4
}

# =========================================================
# 4. runtime library path
# =========================================================

if command -v nf-config >/dev/null 2>&1; then
  NF_PREFIX="$(nf-config --prefix 2>/dev/null || true)"
  if [[ -n "$NF_PREFIX" && -d "$NF_PREFIX/lib" ]]; then
    export LD_LIBRARY_PATH="$NF_PREFIX/lib:${LD_LIBRARY_PATH:-}"
  fi
fi

if command -v nc-config >/dev/null 2>&1; then
  NC_PREFIX="$(nc-config --prefix 2>/dev/null || true)"
  if [[ -n "$NC_PREFIX" && -d "$NC_PREFIX/lib" ]]; then
    export LD_LIBRARY_PATH="$NC_PREFIX/lib:${LD_LIBRARY_PATH:-}"
  fi
fi

# =========================================================
# 5. batch run over forcing_sweep
# =========================================================

mkdir -p "$OUTPUT_ROOT"

echo "[INFO] Starting sweep..."
echo "[INFO] FORCING_ROOT=$FORCING_ROOT"
echo "[INFO] OUTPUT_ROOT=$OUTPUT_ROOT"

echo "[INFO] Parameters:"
echo "       beta        = $beta"
echo "       s_tau       = $s_tau"
echo "       damp1_days  = $damp1_days"
echo "       u_thres     = $u_thres"
echo "       flux_extra  = $flux_extra"
echo "       alpha       = $alpha"
echo "       const1      = $const1"
echo "       ppow        = $ppow"
echo "       aexp        = $aexp"

shopt -s nullglob

for scale_dir in "$FORCING_ROOT"/scale_*; do
  [[ -d "$scale_dir" ]] || continue

  scale_name="$(basename "$scale_dir")"

  echo "======================================================"
  echo "[INFO] Processing scale group: $scale_name"

  forcing_files=("$scale_dir"/*.nc)

  if [[ ${#forcing_files[@]} -eq 0 ]]; then
    echo "[WARN] No .nc files found in $scale_dir"
    continue
  fi

  for forcing_file in "${forcing_files[@]}"; do
    [[ -f "$forcing_file" ]] || continue

    base="$(basename "$forcing_file" .nc)"

    # Example:
    #   forcing_total_1000_ens03.nc -> ens03
    ens_name="${base##*_}"

    out_dir="$OUTPUT_ROOT/$scale_name/$ens_name"
    tag="${scale_name}_${ens_name}"

    echo "------------------------------------------------------"
    echo "[INFO] Running experiment:"
    echo "       forcing = $forcing_file"
    echo "       output  = $out_dir"
    echo "       tag     = $tag"

    # Remove old output for this one experiment only
    rm -rf "$out_dir"
    mkdir -p "$out_dir"

    "$EXE" \
      "$beta" \
      "$s_tau" \
      "$damp1_days" \
      "$u_thres" \
      "$flux_extra" \
      "$alpha" \
      "$const1" \
      "$ppow" \
      "$aexp" \
      "$tag" \
      "$forcing_file" \
      "$out_dir" \
      > "$out_dir/log.txt" 2>&1

    rc=$?

    if [[ $rc -ne 0 ]]; then
      echo "[ERR] Run failed for $forcing_file"
      echo "[ERR] See log: $out_dir/log.txt"
      exit $rc
    fi

    echo "[INFO] Finished $scale_name / $ens_name"
  done
done

echo "[DONE] All runs finished."