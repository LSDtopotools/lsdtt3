// lsdtt_exe_path.hpp
//
// Cross-platform helpers for locating the running executable and the data files
// bundled with lsdtt3. Used so the cosmo tools can find the shipped
// crn_parameters.nc without the user supplying a path.
//
// Layouts searched (relative to the executable's directory), covering both the
// install tree (bin/ + share/) and the build tree (build/ + build/share/):
//   <exe_dir>/../share/lsdtt3/crn_parameters.nc   (installed / packaged)
//   <exe_dir>/share/lsdtt3/crn_parameters.nc      (in-build)
//
#ifndef LSDTT_EXE_PATH_HPP
#define LSDTT_EXE_PATH_HPP

#include <filesystem>
#include <string>
#include <system_error>

#if defined(_WIN32)
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#include <climits>
#include <cstdlib>
#else
#include <unistd.h>
#include <climits>
#endif

namespace lsdtt
{

// Absolute path of the running executable, or an empty path on failure.
inline std::filesystem::path executable_path()
{
#if defined(_WIN32)
    std::wstring buf(MAX_PATH, L'\0');
    DWORD n = GetModuleFileNameW(nullptr, buf.data(),
                                 static_cast<DWORD>(buf.size()));
    if (n == 0)
        return {};
    // Grow if truncated (long paths).
    while (n == buf.size())
    {
        buf.resize(buf.size() * 2);
        n = GetModuleFileNameW(nullptr, buf.data(),
                               static_cast<DWORD>(buf.size()));
        if (n == 0)
            return {};
    }
    buf.resize(n);
    return std::filesystem::path(buf);
#elif defined(__APPLE__)
    char buf[PATH_MAX];
    uint32_t size = sizeof(buf);
    if (_NSGetExecutablePath(buf, &size) != 0)
        return {};  // buffer too small; PATH_MAX is enough in practice
    char real[PATH_MAX];
    if (::realpath(buf, real) != nullptr)
        return std::filesystem::path(real);
    return std::filesystem::path(buf);
#else
    char buf[PATH_MAX];
    ssize_t n = ::readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (n <= 0)
        return {};
    buf[n] = '\0';
    return std::filesystem::path(buf);
#endif
}

// Directory containing the running executable, or an empty path on failure.
inline std::filesystem::path executable_dir()
{
    std::filesystem::path exe = executable_path();
    return exe.empty() ? std::filesystem::path{} : exe.parent_path();
}

// Locate a bundled data file (given its path relative to share/lsdtt3), searching
// the install and build layouts. Returns "" if not found.
inline std::string find_bundled_data(const std::string& rel_to_share_lsdtt3)
{
    namespace fs = std::filesystem;
    const fs::path exe_dir = executable_dir();
    if (exe_dir.empty())
        return "";

    const fs::path candidates[] = {
        exe_dir / ".." / "share" / "lsdtt3" / rel_to_share_lsdtt3,
        exe_dir / "share" / "lsdtt3" / rel_to_share_lsdtt3,
    };
    for (const fs::path& c : candidates)
    {
        std::error_code ec;
        if (fs::exists(c, ec) && !ec)
        {
            fs::path norm = fs::weakly_canonical(c, ec);
            return (ec ? c : norm).string();
        }
    }
    return "";
}

// Path to the bundled crn_parameters.nc, or "" if it cannot be located.
inline std::string default_crn_nc_path()
{
    return find_bundled_data("crn_parameters.nc");
}

} // namespace lsdtt

#endif // LSDTT_EXE_PATH_HPP
