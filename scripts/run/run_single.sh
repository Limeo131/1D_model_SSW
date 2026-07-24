#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# Usage:
#   ./run.sh [beta s_tau damp1_days u_thres flux_extra alpha const1 ppow aexp [tag]]
# =========================================================

if [[ $# -gt 10 ]]; then
  echo "Usage: $0 [beta s_tau damp1_days u_thres flux_extra alpha const1 ppow aexp [tag]]" >&2
  exit 2
fi

beta="${1:-}"
s_tau="${2:-}"
damp1_days="${3:-}"
u_thres="${4:-}"
flux_extra="${5:-}"
alpha="${6:-}"
const1="${7:-}"
ppow="${8:-}"
aexp="${9:-}"
tag="${10:-}"

# =========================================================
# 0. paths
# =========================================================
SRC="${REPO_ROOT:-$(dirname "$0")/../..}/src/fortran/1d_model.f90"
EXE="./d.out"

# =========================================================
# 1. load modules for whirl
# =========================================================
echo "[INFO] Trying GCC + generic NetCDF first..."
module purge

# Load HDF5 first, as in the whirl-working version
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
# 5. remove old output files only if tag is given
# =========================================================

echo "[INFO] removing old output nc files in current directory..."

rm -f aout*.nc
rm -f uout*.nc
rm -f urout*.nc
rm -f fzout*.nc
rm -f fznout*.nc
rm -f abudget*.nc
rm -f ubudget*.nc

echo "[INFO] old output files removed."

# =========================================================
# 6. run
# =========================================================
# echo "[INFO] Running:"
# echo "       $EXE $beta $s_tau $damp1_days $u_thres $flux_extra $alpha $const1 $ppow $aexp $tag"

LOGFILE="run_${tag:-notag}.txt"

echo "[INFO] Running:"
echo "       $EXE $beta $s_tau $damp1_days $u_thres $flux_extra $alpha $const1 $ppow $aexp $tag"
echo "[INFO] writing runtime output to $LOGFILE"

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
  > "$LOGFILE" 2>&1