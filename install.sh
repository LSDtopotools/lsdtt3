#!/usr/bin/env bash
# install.sh - LSDTopoTools3 beta installer for Ubuntu / WSL2
#
# Usage
# -----
#   Automatic, downloads the latest release from GitHub:
#     bash install.sh
#
#   From a local tarball:
#     bash install.sh /path/to/lsdtt3-ubuntu-22.04.tar.gz
#
# What it does
# ------------
#   1. Installs required shared libraries via apt.
#   2. Unpacks pre-built executables to ~/.local/lsdtt3/bin.
#   3. Safely adds that directory to PATH in ~/.bashrc.
#
# Safety notes
# ------------
#   - This script does NOT overwrite ~/.bashrc.
#   - This script does NOT create or modify ~/.bash_profile.
#   - Before editing ~/.bashrc, it creates a timestamped backup.
#   - It adds a clearly marked block to ~/.bashrc only if one is not present.

set -euo pipefail

# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------
INSTALL_DIR="$HOME/.local/lsdtt3"
GITHUB_REPO="LSDtopotools/lsdtt3"
ASSET_NAME="lsdtt3-ubuntu-22.04.tar.gz"

RUNTIME_PACKAGES=(
    libgdal30
    libfftw3-double3
    libopencv-core4.5d
    libopencv-imgproc4.5d
    libgomp1
)

FALLBACK_PACKAGES=(
    libgdal-dev
    libfftw3-dev
    libopencv-dev
    libgomp1
)

BASHRC="$HOME/.bashrc"

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------
info() {
    echo "[lsdtt3] $*"
}

warn() {
    echo "[lsdtt3] WARNING: $*" >&2
}

error() {
    echo "[lsdtt3] ERROR: $*" >&2
    exit 1
}

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# -------------------------------------------------------------------------
# Checks
# -------------------------------------------------------------------------
check_ubuntu() {
    if [[ ! -r /etc/os-release ]]; then
        error "Cannot read /etc/os-release. This installer expects Ubuntu."
    fi

    if ! grep -qi '^ID=ubuntu\|^ID_LIKE=.*ubuntu' /etc/os-release; then
        error "This installer is for Ubuntu. For other distributions, build from source."
    fi
}

check_requirements() {
    command_exists sudo || error "sudo is required but was not found."
    command_exists apt-get || error "apt-get is required but was not found."
    command_exists tar || error "tar is required but was not found."
    command_exists grep || error "grep is required but was not found."
    command_exists cut || error "cut is required but was not found."

    if [[ $# -eq 0 ]]; then
        command_exists curl || error "curl is required for downloading releases. Install it with: sudo apt-get install curl"
    fi
}

# -------------------------------------------------------------------------
# Install apt dependencies
# -------------------------------------------------------------------------
install_apt_packages() {
    info "Installing required runtime libraries..."
    sudo apt-get update -qq

    if sudo apt-get install -y "${RUNTIME_PACKAGES[@]}" 2>/dev/null; then
        info "Runtime libraries installed."
    else
        info "Versioned runtime packages not found."
        info "Falling back to development packages..."
        sudo apt-get install -y "${FALLBACK_PACKAGES[@]}"
    fi
}

# -------------------------------------------------------------------------
# Download or use local tarball
# -------------------------------------------------------------------------
download_latest_release() {
    local tmpdir
    local url

    tmpdir="$(mktemp -d)"

    info "Fetching latest release from GitHub..."

    url="$(
        curl -fsSL "https://api.github.com/repos/${GITHUB_REPO}/releases/latest" |
            grep "browser_download_url" |
            grep "$ASSET_NAME" |
            cut -d '"' -f 4
    )"

    if [[ -z "${url:-}" ]]; then
        error "Could not find ${ASSET_NAME} in the latest release.
The repository may be private. Ask the developer to share the .tar.gz directly
and re-run:

  bash install.sh /path/to/${ASSET_NAME}"
    fi

    info "Downloading: $url"
    curl -fsSL -o "$tmpdir/$ASSET_NAME" "$url"

    echo "$tmpdir/$ASSET_NAME"
}

get_tarball() {
    if [[ $# -ge 1 ]]; then
        if [[ -f "$1" ]]; then
            info "Using local tarball: $1"
            echo "$1"
        else
            error "Specified tarball does not exist: $1"
        fi
    else
        download_latest_release
    fi
}

# -------------------------------------------------------------------------
# Install binaries
# -------------------------------------------------------------------------
unpack_and_install() {
    local tarball="$1"

    info "Installing to ${INSTALL_DIR}/bin ..."
    mkdir -p "$INSTALL_DIR/bin"

    tar xz \
        --strip-components=2 \
        -C "$INSTALL_DIR/bin" \
        -f "$tarball" \
        --wildcards "*/bin/lsdtt-*"

    if ! compgen -G "$INSTALL_DIR/bin/lsdtt-*" >/dev/null; then
        error "No lsdtt-* executables were found after unpacking. The tarball may have an unexpected layout."
    fi

    chmod +x "$INSTALL_DIR"/bin/lsdtt-*
    info "Installed $(find "$INSTALL_DIR/bin" -maxdepth 1 -type f -name 'lsdtt-*' | wc -l) executables."
}

# -------------------------------------------------------------------------
# Safely modify ~/.bashrc
# -------------------------------------------------------------------------
ensure_bashrc_exists() {
    if [[ -f "$BASHRC" ]]; then
        return
    fi

    warn "$BASHRC does not exist."

    if [[ -f /etc/skel/.bashrc ]]; then
        info "Creating $BASHRC from /etc/skel/.bashrc"
        cp /etc/skel/.bashrc "$BASHRC"
    else
        info "Creating minimal $BASHRC"
        cat > "$BASHRC" <<'EOF'
# ~/.bashrc
# If not running interactively, do not do anything.
case $- in
    *i*) ;;
      *) return;;
esac
EOF
    fi
}

backup_bashrc() {
    local backup

    backup="$HOME/.bashrc.lsdtt3-backup.$(date +%Y%m%d-%H%M%S)"
    cp "$BASHRC" "$backup"
    info "Backed up existing ~/.bashrc to: $backup"
}

add_to_path() {
    local begin_marker="# >>> LSDTopoTools3 >>>"
    local end_marker="# <<< LSDTopoTools3 <<<"

    ensure_bashrc_exists

    if grep -qF "$begin_marker" "$BASHRC"; then
        info "LSDTopoTools3 PATH block already present in ~/.bashrc"
        return
    fi

    if grep -qF "$HOME/.local/lsdtt3/bin" "$BASHRC" || grep -qF '.local/lsdtt3/bin' "$BASHRC"; then
        info "An lsdtt3 PATH entry already appears to exist in ~/.bashrc"
        return
    fi

    backup_bashrc

    cat >> "$BASHRC" <<'EOF'

# >>> LSDTopoTools3 >>>
# Add LSDTopoTools3 executables to PATH.
if [ -d "$HOME/.local/lsdtt3/bin" ]; then
    case ":$PATH:" in
        *":$HOME/.local/lsdtt3/bin:"*) ;;
        *) export PATH="$HOME/.local/lsdtt3/bin:$PATH" ;;
    esac
fi
# <<< LSDTopoTools3 <<<
EOF

    info "Added LSDTopoTools3 PATH block to ~/.bashrc"
}

# -------------------------------------------------------------------------
# Detect potentially problematic ~/.bash_profile
# -------------------------------------------------------------------------
check_bash_profile_notice() {
    local bash_profile="$HOME/.bash_profile"

    if [[ ! -f "$bash_profile" ]]; then
        return
    fi

    if grep -qF ". ~/.bashrc" "$bash_profile" || grep -qF "source ~/.bashrc" "$bash_profile" || grep -qF ". \$HOME/.bashrc" "$bash_profile" || grep -qF "source \$HOME/.bashrc" "$bash_profile"; then
        return
    fi

    warn "You have a ~/.bash_profile, but it does not appear to source ~/.bashrc."
    warn "Some login shells may not pick up the LSDTopoTools3 PATH setting."
    warn "The installer will not modify ~/.bash_profile automatically, to avoid changing your shell startup behaviour."
    warn "If needed, add this to ~/.bash_profile manually:"
    echo ""
    echo '    if [ -f ~/.bashrc ]; then'
    echo '        . ~/.bashrc'
    echo '    fi'
    echo ""
}

# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------
main() {
    check_ubuntu
    check_requirements "$@"
    install_apt_packages

    local tarball
    tarball="$(get_tarball "$@")"

    unpack_and_install "$tarball"
    add_to_path
    check_bash_profile_notice

    echo ""
    echo "================================================================"
    echo " LSDTopoTools3 installed successfully."
    echo ""
    echo " To use it in this terminal now, run:"
    echo ""
    echo "   source ~/.bashrc"
    echo ""
    echo " Or close and reopen your terminal."
    echo ""
    echo " Then verify the installation:"
    echo ""
    echo "   lsdtt-channel-extraction --help"
    echo ""
    echo " If the command is not found, check that ~/.bashrc is being loaded"
    echo " by your terminal/shell configuration."
    echo "================================================================"
}

main "$@"
