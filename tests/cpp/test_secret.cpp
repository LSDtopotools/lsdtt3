// test_secret.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Phase 3 of docs/fetch_raster_cli_plan.md: tests for credential handling.
//
// lsdtt_secret is the only place in lsdtt3 that touches secrets, so it gets its
// own test driver rather than being folded into the fetch tests. Everything
// here runs offline and deterministically.
//
// The properties under test are security properties, not conveniences:
//   - a secret never appears in any field intended for logging
//   - masking actually removes every occurrence
//   - resolution order is exactly env var, then key file, then prompt
//   - prompting never happens when stdin is not a terminal (it would hang CI)
//   - a world-readable key file is reported
//
// Usage: lsdtt-test-secret [work_dir]
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#ifndef _WIN32
#  include <sys/stat.h>
#endif

#include "lsdtt_secret.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

void set_env(const char* name, const char* value)
{
#ifdef _WIN32
    _putenv_s(name, value ? value : "");
#else
    if (value) setenv(name, value, 1);
    else       unsetenv(name);
#endif
}

void write_file(const std::string& path, const std::string& content)
{
    std::ofstream out(path);
    out << content;
}


// =============================================================================
// mask_secret
// =============================================================================
void test_masking()
{
    const std::string key = "abc123SECRET";
    const std::string url =
        "https://portal.opentopography.org/API/globaldem?demtype=COP30&API_Key=" + key;

    const std::string masked = lsdtt::mask_secret(url, key);

    check("mask: the key is gone",
          masked.find(key) == std::string::npos, masked);
    check("mask: replaced with ***",
          masked.find("***") != std::string::npos, masked);
    check("mask: the rest of the URL survives for diagnosis",
          masked.find("demtype=COP30") != std::string::npos, masked);

    // A key repeated in one message (e.g. URL echoed twice in an error) must be
    // removed everywhere, not just the first time.
    const std::string twice = url + " retrying " + url;
    const std::string masked_twice = lsdtt::mask_secret(twice, key);
    check("mask: every occurrence is replaced",
          masked_twice.find(key) == std::string::npos, masked_twice);

    // An empty secret must not match at every position and shred the text.
    const std::string untouched = lsdtt::mask_secret(url, "");
    check("mask: empty secret leaves text unchanged",
          untouched == url);

    check("mask: absent secret leaves text unchanged",
          lsdtt::mask_secret("no key here", key) == "no key here");
}


// =============================================================================
// Key files
// =============================================================================
void test_key_file_reading(const std::string& dir)
{
    std::string warning;

    {
        const std::string path = dir + "/plain_key.txt";
        write_file(path, "mykey123\n");
        check("key file: plain key read",
              lsdtt::read_secret_from_file(path, warning) == "mykey123");
    }
    {
        // Users annotate key files; comments and blank lines must be skipped so
        // a commented file does not silently yield "#" as the key.
        const std::string path = dir + "/commented_key.txt";
        write_file(path, "# OpenTopography key for project X\n\n   \nrealkey456\n");
        check("key file: comments and blank lines skipped",
              lsdtt::read_secret_from_file(path, warning) == "realkey456");
    }
    {
        const std::string path = dir + "/whitespace_key.txt";
        write_file(path, "   spacedkey   \n");
        check("key file: surrounding whitespace trimmed",
              lsdtt::read_secret_from_file(path, warning) == "spacedkey");
    }
    {
        const std::string path = dir + "/empty_key.txt";
        write_file(path, "# only a comment\n\n");
        check("key file: no usable key returns empty",
              lsdtt::read_secret_from_file(path, warning).empty());
    }
    {
        check("key file: missing file returns empty",
              lsdtt::read_secret_from_file(dir + "/does_not_exist.txt", warning).empty());
    }

#ifndef _WIN32
    {
        // A key readable by other users on a shared machine is a real exposure.
        // Warn rather than refuse: the user may be mid-analysis and unable to
        // act immediately, and refusing would just push them to an env var.
        const std::string path = dir + "/loose_key.txt";
        write_file(path, "loosekey789\n");
        chmod(path.c_str(), 0644);

        std::string loose_warning;
        const std::string value = lsdtt::read_secret_from_file(path, loose_warning);
        check("key file: world-readable key still read", value == "loosekey789");
        check("key file: world-readable key warns",
              !loose_warning.empty() &&
              loose_warning.find("chmod 600") != std::string::npos,
              loose_warning);
        check("key file: the warning does not contain the key",
              loose_warning.find("loosekey789") == std::string::npos,
              loose_warning);

        const std::string tight = dir + "/tight_key.txt";
        write_file(tight, "tightkey\n");
        chmod(tight.c_str(), 0600);

        // chmod is a silent no-op on some mounts -- notably NTFS under WSL,
        // where the build tree often lives. Verify the mode actually took
        // before asserting on behaviour that depends on it, otherwise this
        // reports a filesystem limitation as a code defect.
        struct stat tight_stat{};
        const bool mode_applied =
            (stat(tight.c_str(), &tight_stat) == 0) &&
            ((tight_stat.st_mode & (S_IRGRP | S_IROTH)) == 0);

        if (mode_applied)
        {
            std::string tight_warning;
            lsdtt::read_secret_from_file(tight, tight_warning);
            check("key file: 0600 key file produces no warning",
                  tight_warning.empty(), tight_warning);
        }
        else
        {
            std::cout << "SKIP  key file: 0600 key file produces no warning"
                         "   [chmod does not take effect on this filesystem]\n";
        }
    }
#endif
}


// =============================================================================
// Resolution order
// =============================================================================
void test_resolution_order(const std::string& dir)
{
    const std::string key_path = dir + "/order_key.txt";
    write_file(key_path, "from_the_file\n");

    lsdtt::lsdtt_secret_request request;
    request.env_var      = "LSDTT_TEST_SECRET";
    request.file_env_var = "LSDTT_TEST_SECRET_FILE";
    request.allow_prompt = false;   // never prompt in a test

    // --- environment variable wins ---
    set_env("LSDTT_TEST_SECRET", "from_the_env");
    set_env("LSDTT_TEST_SECRET_FILE", key_path.c_str());
    {
        const auto result = lsdtt::resolve_secret(request);
        check("order: environment variable takes precedence over key file",
              result.found && result.value == "from_the_env" &&
              result.source == "environment",
              result.source + " / " + (result.found ? "found" : "not found"));
    }

    // --- key file used when the env var is unset ---
    set_env("LSDTT_TEST_SECRET", nullptr);
    {
        const auto result = lsdtt::resolve_secret(request);
        check("order: key file used when the environment variable is unset",
              result.found && result.value == "from_the_file" &&
              result.source == "key file",
              result.source);
    }

    // --- a named but unreadable key file is reported, not silently skipped ---
    set_env("LSDTT_TEST_SECRET_FILE", (dir + "/nope.txt").c_str());
    {
        const auto result = lsdtt::resolve_secret(request);
        check("order: unreadable key file is a reported failure",
              !result.found && !result.message.empty() &&
              result.message.find("nope.txt") != std::string::npos,
              result.message);
    }

    // --- nothing configured ---
    set_env("LSDTT_TEST_SECRET_FILE", nullptr);
    {
        const auto result = lsdtt::resolve_secret(request);
        check("order: nothing configured yields a clear failure",
              !result.found && !result.message.empty(), result.message);
        check("order: failure message names both environment variables",
              result.message.find("LSDTT_TEST_SECRET") != std::string::npos &&
              result.message.find("LSDTT_TEST_SECRET_FILE") != std::string::npos,
              result.message);
        check("order: failure message states keys are not read from argv or params",
              result.message.find("command-line argument") != std::string::npos &&
              result.message.find("parameter file") != std::string::npos,
              result.message);
    }

    // --- prompting must never happen off a terminal ---
    // Under CTest stdin is not a TTY, so this exercises the real guard: without
    // it, a CI run would block forever on a prompt instead of failing.
    {
        lsdtt::lsdtt_secret_request prompting = request;
        prompting.allow_prompt = true;
        prompting.prompt = "should never be shown: ";

        const auto result = lsdtt::resolve_secret(prompting);
        const bool safe = lsdtt::stdin_is_terminal() ? true : !result.found;
        check("order: no prompt when stdin is not a terminal",
              safe,
              "stdin_is_terminal=" + std::to_string(lsdtt::stdin_is_terminal()));
    }
}


// =============================================================================
// Loggable fields never carry the secret
// =============================================================================
void test_no_secret_in_loggable_fields(const std::string& dir)
{
    const std::string secret_value = "TOPSECRETKEYVALUE";
    const std::string key_path = dir + "/leak_key.txt";
    write_file(key_path, secret_value + "\n");

    lsdtt::lsdtt_secret_request request;
    request.env_var      = "LSDTT_TEST_LEAK";
    request.file_env_var = "LSDTT_TEST_LEAK_FILE";
    request.allow_prompt = false;

    set_env("LSDTT_TEST_LEAK", secret_value.c_str());
    {
        const auto result = lsdtt::resolve_secret(request);
        const std::string loggable = result.source + result.warning + result.message;
        check("leakage: env-resolved secret absent from source/warning/message",
              loggable.find(secret_value) == std::string::npos, loggable);
    }

    set_env("LSDTT_TEST_LEAK", nullptr);
    set_env("LSDTT_TEST_LEAK_FILE", key_path.c_str());
    {
        const auto result = lsdtt::resolve_secret(request);
        const std::string loggable = result.source + result.warning + result.message;
        check("leakage: file-resolved secret absent from source/warning/message",
              loggable.find(secret_value) == std::string::npos, loggable);
    }

    set_env("LSDTT_TEST_LEAK_FILE", nullptr);
}

} // namespace


int main(int argc, char** argv)
{
    const std::string work_dir = (argc > 1) ? argv[1] : ".";

    std::cout << "lsdtt3 credential handling tests "
                 "(docs/fetch_raster_cli_plan.md Phase 3)\n\n";

    test_masking();
    std::cout << "\n";
    test_key_file_reading(work_dir);
    std::cout << "\n";
    test_resolution_order(work_dir);
    std::cout << "\n";
    test_no_secret_in_loggable_fields(work_dir);

    std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks
              << " checks passed\n";
    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return (g_n_failed == 0) ? 0 : 1;
}
