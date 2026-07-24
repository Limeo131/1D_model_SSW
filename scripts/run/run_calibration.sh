# #!/usr/bin/env bash
# set -euo pipefail

# if [[ $# -gt 10 ]]; then
#   echo "Usage: $0 [beta s_tau damp1_days u_thres flux_extra alpha const1 ppow aexp [tag]]" >&2
#   exit 2
# fi

# beta="${1:-}"
# s_tau="${2:-}"
# damp1_days="${3:-}"
# u_thres="${4:-}"
# flux_extra="${5:-}"
# alpha="${6:-}"
# const1="${7:-}"
# ppow="${8:-}"
# aexp="${9:-}"
# tag="${10:-}"

# SRC="/mnt/winds/home/smliu01/1.reproduce_nfl20_10.2020-1.2021/1d_refine/1d_para/graybox_mixing_2026/1d_model_urad_3.2026.f90"
# EXE="./d.out"

# NETCDF_ROOT="/software/netcdf-4.9.0"
# INC="-I${NETCDF_ROOT}/include"
# LIBS="-L${NETCDF_ROOT}/lib -lnetcdff -lnetcdf"
# RPATH="-Wl,-rpath,${NETCDF_ROOT}/lib"

# if command -v nf-config >/dev/null 2>&1; then
#   INC="$(nf-config --fflags) -I\"$(nf-config --includedir)\""
#   LIBS="$(nf-config --flibs)"
#   RPATH="-Wl,-rpath,\"$(nf-config --libdir)\""
# fi

# need_build=0
# if [[ ! -x "$EXE" ]]; then
#   need_build=1
# elif [[ "$SRC" -nt "$EXE" ]]; then
#   need_build=1
# fi

# if [[ $need_build -eq 1 ]]; then
#   echo "[INFO] building d.out with gfortran + ${NETCDF_ROOT}"
#   gfortran -O3 -ffree-line-length-none $INC "$SRC" $LIBS $RPATH -o "$EXE"
# fi

# export LD_LIBRARY_PATH="${NETCDF_ROOT}/lib:${LD_LIBRARY_PATH:-}"

# # If tag is not given, do not delete tagged files
# if [[ -n "$tag" ]]; then
#   UOUT="uout_${tag}.nc"
#   AOUT="aout_${tag}.nc"
#   UROUT="urout_${tag}.nc"
#   FZOUT="fzout_${tag}.nc"
#   FZNOUT="fznout_${tag}.nc"
#   ABUDGET="abudget_${tag}.nc"
#   UBUDGET="ubudget_${tag}.nc"

#   for f in "$UOUT" "$AOUT" "$UROUT" "$FZOUT" "$FZNOUT" "$ABUDGET" "$UBUDGET"; do
#     if [[ -e "$f" ]]; then
#       echo "[INFO] removing existing file: $f"
#       rm -f "$f"
#     fi
#   done
# fi

# echo "[INFO] Running: $EXE $beta $s_tau $damp1_days $u_thres $flux_extra $alpha $const1 $ppow $aexp $tag"
# "$EXE" "$beta" "$s_tau" "$damp1_days" "$u_thres" "$flux_extra" "$alpha" "$const1" "$ppow" "$aexp" "$tag"


#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# Single-run version for new server
#
# Usage:
#   ./run_one_new.sh [beta s_tau damp1_days u_thres flux_extra alpha const1 ppow aexp [tag]]
# =========================================================

if [[ $# -gt 10 ]]; then
  echo "Usage: $0 [beta s_tau damp1_days u_thres flux_extra alpha const1 ppow aexp [tag]]" >&2
  exit 2
fi

# =========================================================
# 0. parameters
# =========================================================

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
# 1. paths
# =========================================================
# Fortran source path
SRC="/nas/winds-home/smliu01/1.reproduce_nfl20_10.2020-1.2021/1d_refine/1d_para/graybox_mixing_2026/1d_model_urad_3.2026.f90"

EXE="./d.out"

# =========================================================
# 2. load modules
# =========================================================

echo "[INFO] Trying GCC + generic NetCDF first..."
module purge

# Load HDF5 first if available
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
# 3. decide whether rebuild is needed
# =========================================================

need_build=0

if [[ ! -x "$EXE" ]]; then
  need_build=1
elif [[ "$SRC" -nt "$EXE" ]]; then
  need_build=1
fi

# =========================================================
# 4. compile
# =========================================================

if [[ $need_build -eq 1 ]]; then
  echo "[INFO] building $EXE"

  NETCDF_PREFIX="$(nc-config --prefix 2>/dev/null || true)"
  USE_GFORTRAN=1

  # If NetCDF was built with Intel, use ifort
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
      LIBDIR="$(dirname "$(find "$PREFIX" -name 'libnetcdff.*' 2>/dev/null | head -n1)")"

      [[ -n "$INC_F" && -n "$LIBDIR" ]] || {
        echo "[ERR] NetCDF-Fortran headers/libs not found under $PREFIX"
        exit 3
      }

      INCFLAGS="-I${INC_C} -I${INC_F}"
      LIBFLAGS="-lnetcdff ${LIBS_C}"
    fi

    echo "[INFO] INCFLAGS=$INCFLAGS"
    echo "[INFO] LIBFLAGS=$LIBFLAGS"
    echo "[INFO] LIBDIR=$LIBDIR"

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
    echo "[INFO] LIBDIR=$LIBDIR"
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
# 5. runtime library path
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
# 6. remove old output files only when tag is given
# =========================================================

if [[ -n "$tag" ]]; then
  UOUT="uout_${tag}.nc"
  AOUT="aout_${tag}.nc"
  UROUT="urout_${tag}.nc"
  FZOUT="fzout_${tag}.nc"
  FZNOUT="fznout_${tag}.nc"
  ABUDGET="abudget_${tag}.nc"
  UBUDGET="ubudget_${tag}.nc"

  for f in "$UOUT" "$AOUT" "$UROUT" "$FZOUT" "$FZNOUT" "$ABUDGET" "$UBUDGET"; do
    if [[ -e "$f" ]]; then
      echo "[INFO] removing existing file: $f"
      rm -f "$f"
    fi
  done
fi

# =========================================================
# 7. run
# =========================================================

echo "[INFO] Running:"
echo "       $EXE $beta $s_tau $damp1_days $u_thres $flux_extra $alpha $const1 $ppow $aexp $tag"

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
  "$tag"

echo "[DONE] run finished."