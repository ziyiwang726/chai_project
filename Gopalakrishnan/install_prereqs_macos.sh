#!/usr/bin/env bash
set -euo pipefail

TOOLS_URL="https://mac.r-project.org/tools/"

log() {
  printf '[install-prereqs] %s\n' "$*"
}

fail() {
  printf '[install-prereqs] ERROR: %s\n' "$*" >&2
  exit 1
}

if [[ "$(uname -s)" != "Darwin" ]]; then
  fail "This script is for macOS only."
fi

if ! command -v curl >/dev/null 2>&1; then
  fail "curl is required but not found."
fi

log "Checking Xcode Command Line Tools..."
if xcode-select -p >/dev/null 2>&1; then
  log "Xcode Command Line Tools already installed."
else
  log "Installing Xcode Command Line Tools (GUI prompt may appear)..."
  xcode-select --install || true
  log "If installation prompt appeared, complete it, then rerun this script."
fi

if [[ -x /opt/gfortran/bin/gfortran ]]; then
  log "GNU Fortran already installed at /opt/gfortran/bin/gfortran"
else
  log "Fetching latest gfortran package name from ${TOOLS_URL}"
  PAGE="$(curl -fsSL "${TOOLS_URL}")" || fail "Failed to download ${TOOLS_URL}"

  PKG_LIST="$(
    printf '%s\n' "${PAGE}" \
      | perl -nE 'while(/href=["'"'"']([^"'"'"']*gfortran-[0-9][^"'"'"']*\.pkg)["'"'"']/g){$u=$1; $u =~ s#.*/##; say $u}' \
      | awk '!seen[$0]++'
  )"
  [[ -n "${PKG_LIST}" ]] || fail "Could not parse a gfortran pkg from ${TOOLS_URL}"

  SORTED_CANDIDATES="$(
    printf '%s\n' "${PKG_LIST}" \
      | awk '
          {
            v=$0
            sub(/^gfortran-/, "", v)
            sub(/\.pkg$/, "", v)
            n=split(v, a, /[^0-9]+/)
            key=""
            for (i=1; i<=8; i++) {
              num=(i<=n && a[i] != "") ? a[i] : 0
              key=key sprintf("%05d", num)
            }
            print key "|" $0
          }
        ' \
      | sort -r \
      | cut -d'|' -f2
  )"

  PKG_NAME=""
  while IFS= read -r candidate; do
    [[ -z "${candidate}" ]] && continue
    if curl -fsI "${TOOLS_URL}${candidate}" >/dev/null 2>&1; then
      PKG_NAME="${candidate}"
      break
    fi
  done <<< "${SORTED_CANDIDATES}"
  [[ -n "${PKG_NAME}" ]] || fail "Found package names, but none resolved (HTTP 200) at ${TOOLS_URL}"

  PKG_URL="${TOOLS_URL}${PKG_NAME}"
  PKG_PATH="/tmp/${PKG_NAME}"

  log "Downloading ${PKG_URL}"
  curl -fL "${PKG_URL}" -o "${PKG_PATH}" || fail "Download failed."

  log "Installing ${PKG_NAME} (sudo required)..."
  sudo installer -pkg "${PKG_PATH}" -target / || fail "Installer failed."
fi

if [[ -x /opt/gfortran/bin/gfortran ]]; then
  log "GNU Fortran installed:"
  /opt/gfortran/bin/gfortran --version | head -n 1
else
  fail "gfortran is still missing at /opt/gfortran/bin/gfortran"
fi

if command -v R >/dev/null 2>&1; then
  log "R FLIBS:"
  R CMD config FLIBS

  log "Checking R FLIBS -L directories exist..."
  flibs="$(R CMD config FLIBS)"
  bad_dirs=0
  while IFS= read -r dir; do
    [[ -z "${dir}" ]] && continue
    if [[ ! -d "${dir}" ]]; then
      printf '[install-prereqs] MISSING: %s\n' "${dir}" >&2
      bad_dirs=1
    fi
  done < <(printf '%s\n' "${flibs}" | grep -Eo -- '-L[^[:space:]]+' | sed 's/^-L//')

  if [[ "${bad_dirs}" -eq 0 ]]; then
    log "All FLIBS library directories exist."
  else
    fail "Some FLIBS directories are missing. Reinstall gfortran for your R architecture."
  fi
else
  log "R not found in PATH; skipping R FLIBS checks."
fi

log "Prerequisite setup complete."
log "Now run: Rscript Gopalakrishnan/install_packages.R"
