#!/usr/bin/env bash
# install.sh  —  LSDTopoTools3 beta installer for Ubuntu / WSL2
#
# Usage
# -----
#   Automatic (downloads the latest release from GitHub):
#     bash install.sh
#
#   From a local tarball (useful when the repo is private and you have
#   been emailed the .tar.gz directly):
#     bash install.sh /path/to/lsdtt3-ubuntu-22.04.tar.gz
#
# What it does
# ------------
#   1. Installs the required shared libraries via apt.
#   2. Unpacks the pre-built executables to ~/.local/lsdtt3/bin.
#   3. Adds that directory to PATH in ~/.bashrc (once).
#
# Requirements
# ------------
#   Ubuntu 22.04 or later (including inside WSL2).
#   sudo access to install apt packages.

set -euo pipefail

# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------
INSTALL_DIR="$HOME/.local/lsdtt3"
GITHUB_REPO="LSDtopotools/lsdtt3"
ASSET_NAME="lsdtt3-ubuntu-22.04.tar.gz"

# Runtime libraries that the lsdtt3 executables link against.
# These are the same packages installed by the GitHub Actions build, so
# the binary ABI is guaranteed to match.
RUNTIME_PACKAGES=(
    libgdal30            # GDAL raster/vector I/O
    libfftw3-double3     # FFTW3 spectral analysis
    libopencv-core4.5d   # OpenCV image processing (core)
    libopencv-imgproc4.5d
    libgomp1             # OpenMP runtime
)

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------
info()  { echo "[lsdtt3] $*"; }
error() { echo "[lsdtt3] ERROR: $*" >&2; exit 1; }

check_ubuntu() {
    if ! grep -qi ubuntu /etc/os-release 2>/dev/null; then
        error "This installer is for Ubuntu. For other distributions, build from source."
    fi
}

install_apt_packages() {
    info "Installing required runtime libraries..."
    sudo apt-get update -qq

    # Try to install the exact versioned packages; if they are not found
    # (e.g. Ubuntu 24.04 ships newer versions), fall back to the -dev meta-
    # packages which always pull in the correct runtime libraries.
    if sudo apt-get install -y "${RUNTIME_PACKAGES[@]}" 2>/dev/null; then
        info "Runtime libraries installed."
    else
        info "Versioned packages not found, falling back to -dev packages..."
        sudo apt-get install -y libgdal-dev libfftw3-dev libopencv-dev libgomp1
    fi
}

download_latest_release() {
    local tmpdir
    tmpdir=$(mktemp -d)

    info "Fetching latest release from GitHub..."

    # GitHub API: latest release download URL for our asset
    local url
    url=$(curl -fsSL \
        "https://api.github.com/repos/${GITHUB_REPO}/releases/latest" \
        | grep "browser_download_url" \
        | grep "$ASSET_NAME" \
        | cut -d '"' -f 4)

    if [[ -z "$url" ]]; then
        error "Could not find ${ASSET_NAME} in the latest release. \
The repository may be private. Ask the developer to share the .tar.gz directly \
and re-run: bash install.sh /path/to/${ASSET_NAME}"
    fi

    info "Downloading: $url"
    curl -fsSL -o "$tmpdir/$ASSET_NAME" "$url"
    echo "$tmpdir/$ASSET_NAME"
}

unpack_and_install() {
    local tarball="$1"

    info "Installing to ${INSTALL_DIR}/bin ..."
    mkdir -p "$INSTALL_DIR/bin"

    # The tarball unpacks to lsdtt3-ubuntu-22.04/bin/lsdtt-*
    tar xz --strip-components=2 \
        -C "$INSTALL_DIR/bin" \
        -f "$tarball" \
        --wildcards "*/bin/lsdtt-*"

    chmod +x "$INSTALL_DIR"/bin/*
    info "Installed $(ls "$INSTALL_DIR/bin" | wc -l) executables."
}

add_to_path() {
    local pathline='export PATH="$HOME/.local/lsdtt3/bin:$PATH"'

    # Write to ~/.bashrc (non-login interactive shells — most WSL2 terminals).
    if ! grep -qF 'lsdtt3/bin' ~/.bashrc 2>/dev/null; then
        echo "" >> ~/.bashrc
        echo "# LSDTopoTools3" >> ~/.bashrc
        echo "$pathline" >> ~/.bashrc
        info "Added lsdtt3 to PATH in ~/.bashrc"
    else
        info "PATH entry already present in ~/.bashrc"
    fi

    # Also write to ~/.bash_profile (login shells — some WSL2 / Windows Terminal
    # configurations source this instead of ~/.bashrc).
    if ! grep -qF 'lsdtt3/bin' ~/.bash_profile 2>/dev/null; then
        echo "" >> ~/.bash_profile
        echo "# LSDTopoTools3" >> ~/.bash_profile
        echo "$pathline" >> ~/.bash_profile
    fi
}

# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------
check_ubuntu
install_apt_packages

if [[ $# -ge 1 && -f "$1" ]]; then
    # Local tarball supplied as argument
    TARBALL="$1"
    info "Using local tarball: $TARBALL"
else
    TARBALL=$(download_latest_release)
fi

unpack_and_install "$TARBALL"
add_to_path

echo ""
echo "================================================================"
echo " LSDTopoTools3 installed successfully."
echo ""
echo " Restart your terminal, or run:"
echo "   source ~/.bashrc"
echo ""
echo " Then verify the installation:"
echo "   lsdtt-channel-extraction --help"
echo "================================================================"
