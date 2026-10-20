// lsdtt_http_fetch.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Small HTTP requests through GDAL's CPLHTTPFetch (the only HTTP client lsdtt3
// uses), shared by the raster and lithology fetchers. Errors report the URL with
// its query masked, because some services put credentials there.
// GDAL-free header.
// -----------------------------------------------------------------------------

#pragma once

#include <string>

namespace lsdtt
{

struct lsdtt_http_request_options
{
    int timeout_seconds = 60;          // whole request
    int connect_timeout_seconds = 0;   // 0 = GDAL's default
    int max_retry = 3;                 // transient failures, 2 s apart
    std::string user_agent;            // empty = GDAL's default
};

// True if GDAL was built with curl, so HTTP requests can be made at all.
bool http_fetch_available();

// GET (post_body empty) or POST (JSON body) `url`; the response body goes to
// `body`. Returns false with `error` set on a transport failure or an empty
// response. Never throws.
bool http_request_body(const std::string& url, const std::string& post_body,
                       const lsdtt_http_request_options& opts, std::string& body,
                       std::string& error);

} // namespace lsdtt
