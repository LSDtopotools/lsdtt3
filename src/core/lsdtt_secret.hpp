// lsdtt_secret.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Credential resolution and masking.
//
// Testing and benchmarking notes: docs/testing_and_benchmarking.md
//
// This is the ONLY place in lsdtt3 that handles secrets. It is deliberately
// small so that it is the only thing anyone has to audit.
//
// *** THE RULES ***
//
// 1. Never compile a key into the binary. A string literal in an executable is
//    recoverable with `strings`.
//
// 2. Never accept a key as a command-line argument. argv is visible to every
//    user on the machine via `ps`, and lands in shell history. Note this
//    module offers no way to pass a key as an argument, by design.
//
// 3. Never put a key in a parameter file. lsdtt3 CLIs are parameter-file
//    driven, users routinely commit their parameter files alongside analyses,
//    and print_parameters() echoes every parameter at the end of a successful
//    run. A key held as a parameter would be committed to a repository and
//    printed to a log.
//
// 4. Mask a key everywhere it might be logged. GDAL error messages embed the
//    failing URL, so a key in a query string leaks to stderr on any failure
//    unless it is masked first. Use mask_secret() before printing anything
//    that could contain a URL.
//
// Resolution order (first found wins), matching python_helpers/download_dem_tiles.py
// and ../lsdviztools3/lsdviztools3/download/apikey.py so that a user who has
// already configured the Python helpers needs no extra setup:
//
//   1. Environment variable          e.g. OT_API_KEY
//   2. Key file named by an env var  e.g. OT_API_KEY_FILE
//   3. Interactive prompt            only when stdin is a terminal
//
// The plan originally proposed a --api-key-file command-line flag. That was
// dropped: lsdtt_param_parser::parse_args_deferred() interprets extra argv
// entries as <path> <parameter_file>, so a new flag would require changing the
// shared parser used by every CLI. An environment variable naming the file
// achieves the same thing with no blast radius.
// -----------------------------------------------------------------------------

#pragma once

#include <string>

namespace lsdtt {

// -----------------------------------------------------------------------------
// lsdtt_secret_request
// -----------------------------------------------------------------------------
struct lsdtt_secret_request
{
    // Environment variable holding the secret itself, e.g. "OT_API_KEY".
    std::string env_var;

    // Environment variable holding a PATH to a file containing the secret,
    // e.g. "OT_API_KEY_FILE". The first non-blank, non-comment line is used.
    std::string file_env_var;

    // Shown when prompting. Should name the provider so the user knows what is
    // being asked for.
    std::string prompt;

    // Prompting is skipped when stdin is not a terminal regardless of this
    // flag; set false to suppress it even on a terminal (e.g. in tests).
    bool allow_prompt = true;
};


// -----------------------------------------------------------------------------
// lsdtt_secret_result
// -----------------------------------------------------------------------------
struct lsdtt_secret_result
{
    bool found = false;

    // The secret. Treat as sensitive: never log, never echo, never write to a
    // parameter file. Pass through mask_secret() before printing anything that
    // may contain it.
    std::string value;

    // Which mechanism supplied it: "environment", "key file", "prompt", or
    // "" when not found. Safe to log -- names the source, never the value.
    std::string source;

    // Non-fatal advisories, e.g. a key file readable by other users. Safe to
    // log; contains no secret material.
    std::string warning;

    // Why resolution failed, with actionable guidance. Safe to log.
    std::string message;
};


// -----------------------------------------------------------------------------
// resolve_secret
//
// Never throws, never prints the secret, and never reads argv.
// -----------------------------------------------------------------------------
lsdtt_secret_result resolve_secret(const lsdtt_secret_request& request);


// -----------------------------------------------------------------------------
// mask_secret
//
// Replace every occurrence of `secret` in `text` with "***". Returns `text`
// unchanged when `secret` is empty, so it is safe to call unconditionally.
//
// Call this on ANY string that might contain a credential before it reaches
// stdout, stderr, a log file or an exception message -- particularly URLs,
// which carry API keys in query strings.
// -----------------------------------------------------------------------------
std::string mask_secret(const std::string& text, const std::string& secret);


// -----------------------------------------------------------------------------
// mask_url_query
//
// Replace a URL's query string with "?***". Use for URLs whose credential is
// the whole query rather than a known substring -- notably Azure SAS tokens,
// where the signature is generated per request so there is no fixed secret to
// pass to mask_secret(). A SAS token grants time-limited read access to the
// blob, so it must not be written to logs even though obtaining one needs no
// account.
// -----------------------------------------------------------------------------
std::string mask_url_query(const std::string& url);


// -----------------------------------------------------------------------------
// read_secret_from_file
//
// First non-blank, non-comment ('#') line of a file, trimmed. Empty return
// means nothing usable was found. Exposed for testing.
// -----------------------------------------------------------------------------
std::string read_secret_from_file(const std::string& path, std::string& warning);


// True if standard input is an interactive terminal. Prompting when this is
// false would hang a CI run instead of failing cleanly.
bool stdin_is_terminal();

} // namespace lsdtt
