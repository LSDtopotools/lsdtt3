// test_secret.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Part of the flow-routing benchmark (docs/testing_and_benchmarking.md): tests for credential handling.
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
//   - prompting never happens when LSDTT_NONINTERACTIVE is set, EVEN WHEN stdin
//     is a terminal -- which it is inside an interactive `ctest` run, because
//     CTest passes the invoking terminal through. Tested on a real
//     pseudo-terminal, so the guard is exercised rather than assumed.
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
#  include <csignal>
#  include <fcntl.h>
#  include <sys/stat.h>
#  include <sys/wait.h>
#  include <unistd.h>
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
// =============================================================================
// The prompt on a REAL terminal
//
// An interactive `ctest` run gives every test the invoking terminal as stdin,
// so "stdin is not a terminal" cannot be what stops a prompt there. These
// checks put a pseudo-terminal on stdin and call resolve_secret() in a FORKED
// CHILD with an alarm, so that a regression fails the check instead of hanging
// the driver -- the very failure mode under test. (Pre-typed input cannot be
// relied on: the prompt sets the terminal with TCSAFLUSH, which discards it.)
//
//   1. With LSDTT_NONINTERACTIVE=1 and nothing typed, the child must return
//      not-found at once. If the switch were ignored it would block at the
//      prompt until the alarm kills it.
//   2. Control: with the switch unset, the parent types a key once the child
//      is waiting, and the child must read it from the prompt. That proves the
//      pseudo-terminal genuinely reaches the prompt, so check 1 is not passing
//      vacuously.
// =============================================================================
#ifndef _WIN32
enum class pty_outcome { refused, prompted, blocked, not_a_terminal, other };

std::string describe(pty_outcome o)
{
    switch (o)
    {
        case pty_outcome::refused:        return "refused without prompting";
        case pty_outcome::prompted:       return "prompted and read the typed key";
        case pty_outcome::blocked:        return "BLOCKED at the prompt (killed by alarm)";
        case pty_outcome::not_a_terminal: return "stdin was not a terminal in the child";
        default:                          return "unexpected outcome";
    }
}

pty_outcome resolve_on_pty(const lsdtt::lsdtt_secret_request& request,
                           int master, int slave,
                           bool noninteractive, bool type_key)
{
    std::cout.flush();
    const pid_t pid = fork();
    if (pid < 0) return pty_outcome::other;

    if (pid == 0)
    {
        dup2(slave, STDIN_FILENO);
        set_env(lsdtt::noninteractive_env_var, noninteractive ? "1" : nullptr);
        if (!lsdtt::stdin_is_terminal()) _exit(13);
        alarm(5);
        const auto r = lsdtt::resolve_secret(request);
        std::cout.flush();
        if (!r.found) _exit(10);
        if (r.source == "prompt" && r.value == "typed_key") _exit(11);
        _exit(12);
    }

    int status = 0;
    bool done = false;
    for (int i = 0; i < 100 && !done; ++i)          // up to ~10 s
    {
        done = (waitpid(pid, &status, WNOHANG) == pid);
        if (done) break;
        if (type_key)
        {
            const char line[] = "typed_key\n";
            if (write(master, line, sizeof(line) - 1) < 0) break;
        }
        usleep(100000);
    }
    if (!done)
    {
        kill(pid, SIGKILL);
        waitpid(pid, &status, 0);
        return pty_outcome::blocked;
    }
    if (WIFSIGNALED(status)) return pty_outcome::blocked;
    switch (WIFEXITED(status) ? WEXITSTATUS(status) : -1)
    {
        case 10: return pty_outcome::refused;
        case 11: return pty_outcome::prompted;
        case 13: return pty_outcome::not_a_terminal;
        default: return pty_outcome::other;
    }
}
#endif

void test_prompt_on_pseudo_terminal(const lsdtt::lsdtt_secret_request& prompting)
{
#ifdef _WIN32
    (void)prompting;
    std::cout << "NOTE  pseudo-terminal prompt test not available on Windows\n";
#else
    const int master = posix_openpt(O_RDWR | O_NOCTTY);
    const bool opened = master >= 0 && grantpt(master) == 0 &&
                        unlockpt(master) == 0 && ptsname(master) != nullptr;
    const int slave = opened ? open(ptsname(master), O_RDWR | O_NOCTTY) : -1;
    if (slave < 0)
    {
        if (master >= 0) close(master);
        std::cout << "NOTE  no pseudo-terminal available; prompt test skipped\n";
        return;
    }

    const pty_outcome refused  = resolve_on_pty(prompting, master, slave,
                                                /*noninteractive=*/true,
                                                /*type_key=*/false);
    const pty_outcome prompted = resolve_on_pty(prompting, master, slave,
                                                /*noninteractive=*/false,
                                                /*type_key=*/true);
    close(slave);
    close(master);

    check("pty: with LSDTT_NONINTERACTIVE=1 a terminal is NOT prompted",
          refused == pty_outcome::refused, describe(refused));
    check("pty: control -- without it, the same terminal IS prompted",
          prompted == pty_outcome::prompted, describe(prompted));
#endif
}


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

    // --- LSDTT_NONINTERACTIVE forbids the prompt ---
    // Deterministic whatever stdin is: with the variable set, a prompting
    // request must come back not-found without touching stdin. CTest sets the
    // variable on every test; it is set explicitly here so the check does not
    // depend on how the driver was launched.
    lsdtt::lsdtt_secret_request prompting = request;
    prompting.allow_prompt = true;
    prompting.prompt = "lsdtt-test-secret prompt: ";

    set_env(lsdtt::noninteractive_env_var, "1");
    {
        check("noninteractive: the switch reads as set for \"1\"",
              lsdtt::noninteractive_requested());
        const auto result = lsdtt::resolve_secret(prompting);
        check("noninteractive: a prompting request is refused, not prompted",
              !result.found && result.message.find("no credential found") == 0,
              result.message);
    }
    set_env(lsdtt::noninteractive_env_var, "0");
    check("noninteractive: \"0\" does not set the switch",
          !lsdtt::noninteractive_requested());
    set_env(lsdtt::noninteractive_env_var, nullptr);
    check("noninteractive: unset does not set the switch",
          !lsdtt::noninteractive_requested());

    // --- prompting must never happen off a terminal ---
    // Only meaningful when stdin is not a terminal. With the switch unset and
    // stdin a terminal, this request WOULD prompt and block, so it is not run.
    if (!lsdtt::stdin_is_terminal())
    {
        const auto result = lsdtt::resolve_secret(prompting);
        check("order: no prompt when stdin is not a terminal", !result.found,
              result.source);
    }
    else
    {
        std::cout << "NOTE  stdin is a terminal; the not-a-terminal check is "
                     "covered by the pseudo-terminal test below instead\n";
    }

    test_prompt_on_pseudo_terminal(prompting);
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
                 " (docs/testing_and_benchmarking.md)\n\n";

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
