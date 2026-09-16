// lsdtt_secret.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// See lsdtt_secret.hpp for the rules this file exists to enforce.

#include "lsdtt_secret.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#ifdef _WIN32
#  include <io.h>
#  include <windows.h>
#else
#  include <sys/stat.h>
#  include <termios.h>
#  include <unistd.h>
#endif

namespace lsdtt {

namespace {

std::string trim(const std::string& text)
{
    const std::string whitespace = " \t\r\n";
    const std::size_t begin = text.find_first_not_of(whitespace);
    if (begin == std::string::npos) return "";
    const std::size_t end = text.find_last_not_of(whitespace);
    return text.substr(begin, end - begin + 1);
}


std::string env_value(const std::string& name)
{
    if (name.empty()) return "";
    // MSVC deprecates getenv in favour of _dupenv_s (C4996). std::getenv is
    // correct and portable here: the value is copied and trimmed immediately,
    // and nothing writes to the returned buffer. Suppressed locally rather than
    // defining _CRT_SECURE_NO_WARNINGS for the whole build.
#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable : 4996)
#endif
    const char* raw = std::getenv(name.c_str());
#ifdef _MSC_VER
#  pragma warning(pop)
#endif
    if (raw == nullptr) return "";
    return trim(raw);
}


// -----------------------------------------------------------------------------
// prompt_for_secret
//
// Reads a line from the terminal with echo disabled, so the key does not appear
// on screen or in a screen-shared session. Returns empty on any failure --
// callers treat that as "not found" rather than proceeding with a partial key.
// -----------------------------------------------------------------------------
std::string prompt_for_secret(const std::string& prompt_text)
{
    std::cout << prompt_text << std::flush;

    std::string value;

#ifdef _WIN32
    HANDLE handle = GetStdHandle(STD_INPUT_HANDLE);
    DWORD  original_mode = 0;
    const bool have_mode = (handle != INVALID_HANDLE_VALUE) &&
                           GetConsoleMode(handle, &original_mode);
    if (have_mode)
    {
        SetConsoleMode(handle, original_mode & ~ENABLE_ECHO_INPUT);
    }

    std::getline(std::cin, value);

    if (have_mode)
    {
        SetConsoleMode(handle, original_mode);
    }
#else
    termios original{};
    const bool have_termios = (tcgetattr(STDIN_FILENO, &original) == 0);
    if (have_termios)
    {
        termios quiet = original;
        quiet.c_lflag &= static_cast<tcflag_t>(~ECHO);
        tcsetattr(STDIN_FILENO, TCSAFLUSH, &quiet);
    }

    std::getline(std::cin, value);

    if (have_termios)
    {
        // Restore echo before anything else can write to the terminal, or the
        // user is left with a terminal that silently swallows their typing.
        tcsetattr(STDIN_FILENO, TCSAFLUSH, &original);
    }
#endif

    // The newline the user typed was not echoed, so emit one to keep subsequent
    // output tidy.
    std::cout << "\n";

    return trim(value);
}

} // anonymous namespace


// =============================================================================
// Public API
// =============================================================================

bool stdin_is_terminal()
{
#ifdef _WIN32
    return _isatty(_fileno(stdin)) != 0;
#else
    return isatty(STDIN_FILENO) != 0;
#endif
}


std::string read_secret_from_file(const std::string& path, std::string& warning)
{
    warning.clear();
    if (path.empty()) return "";

    std::ifstream stream(path);
    if (!stream)
    {
        return "";
    }

#ifndef _WIN32
    // A key file readable by other users is a real exposure, but refusing to
    // read it would be worse than telling the user: they may be mid-analysis on
    // a shared machine and unable to act immediately.
    struct stat file_stat{};
    if (stat(path.c_str(), &file_stat) == 0)
    {
        if (file_stat.st_mode & (S_IRGRP | S_IROTH))
        {
            warning = "key file '" + path + "' is readable by other users; "
                      "consider: chmod 600 " + path;
        }
    }
#endif

    std::string line;
    while (std::getline(stream, line))
    {
        const std::string trimmed = trim(line);
        if (trimmed.empty())        continue;
        if (trimmed[0] == '#')      continue;
        return trimmed;
    }

    return "";
}


lsdtt_secret_result resolve_secret(const lsdtt_secret_request& request)
{
    lsdtt_secret_result result;

    // ---- 1. Environment variable holding the secret ----
    const std::string from_env = env_value(request.env_var);
    if (!from_env.empty())
    {
        result.found  = true;
        result.value  = from_env;
        result.source = "environment";
        return result;
    }

    // ---- 2. Environment variable naming a key file ----
    const std::string key_file = env_value(request.file_env_var);
    if (!key_file.empty())
    {
        std::string warning;
        const std::string from_file = read_secret_from_file(key_file, warning);
        if (!from_file.empty())
        {
            result.found   = true;
            result.value   = from_file;
            result.source  = "key file";
            result.warning = warning;
            return result;
        }

        // Naming a file that cannot be read is a configuration mistake worth
        // reporting precisely, rather than silently falling through to a
        // prompt that may not be possible.
        result.message = "key file '" + key_file + "' (from " +
                         request.file_env_var +
                         ") could not be read, or contains no usable key. "
                         "The first non-blank line that does not start with "
                         "'#' is used as the key.";
        return result;
    }

    // ---- 3. Interactive prompt ----
    // Only on a terminal: prompting in CI would hang the run rather than
    // failing cleanly, which is far worse than a clear error.
    if (request.allow_prompt && stdin_is_terminal())
    {
        const std::string prompt_text =
            request.prompt.empty() ? std::string("Enter API key: ")
                                   : request.prompt;
        const std::string typed = prompt_for_secret(prompt_text);
        if (!typed.empty())
        {
            result.found  = true;
            result.value  = typed;
            result.source = "prompt";
            return result;
        }
        result.message = "no key entered at the prompt.";
        return result;
    }

    result.message =
        "no credential found. Set the " + request.env_var +
        " environment variable, or set " + request.file_env_var +
        " to the path of a file whose first line is the key "
        "(create it with restrictive permissions: chmod 600). "
        "The key is deliberately never read from a command-line argument or a "
        "parameter file, because both are routinely logged or committed.";
    return result;
}


std::string mask_url_query(const std::string& url)
{
    const std::size_t question = url.find('?');
    if (question == std::string::npos) return url;
    return url.substr(0, question) + "?***";
}


std::string mask_secret(const std::string& text, const std::string& secret)
{
    // An empty secret would otherwise match at every position and corrupt the
    // text, so the guard is load-bearing, not defensive noise.
    if (secret.empty()) return text;

    std::string masked = text;
    std::size_t position = masked.find(secret);
    while (position != std::string::npos)
    {
        masked.replace(position, secret.size(), "***");
        position = masked.find(secret, position + 3);
    }
    return masked;
}

} // namespace lsdtt
