// lsdtt_http_fetch.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_http_fetch.hpp"

#include <cpl_conv.h>
#include <cpl_error.h>
#include <cpl_http.h>
#include <cpl_string.h>

#include "lsdtt_gdal_init.hpp"
#include "lsdtt_secret.hpp"

namespace lsdtt
{

bool http_fetch_available()
{
    gdal_init();
    return CPLHTTPEnabled() != 0;
}

bool http_request_body(const std::string& url, const std::string& post_body,
                       const lsdtt_http_request_options& o, std::string& body,
                       std::string& error)
{
    gdal_init();
    char** opts = nullptr;
    opts = CSLAddString(opts, ("TIMEOUT=" + std::to_string(o.timeout_seconds)).c_str());
    opts = CSLAddString(opts, ("MAX_RETRY=" + std::to_string(o.max_retry)).c_str());
    opts = CSLAddString(opts, "RETRY_DELAY=2");
    if (o.connect_timeout_seconds > 0)
        opts = CSLAddString(opts,
                            ("CONNECTTIMEOUT=" + std::to_string(o.connect_timeout_seconds)).c_str());
    if (!o.user_agent.empty())
        opts = CSLAddString(opts, ("USERAGENT=" + o.user_agent).c_str());
    if (!post_body.empty())
    {
        opts = CSLAddString(opts, ("POSTFIELDS=" + post_body).c_str());
        opts = CSLAddString(opts, "HEADERS=Content-Type: application/json");
    }

    CPLPushErrorHandler(CPLQuietErrorHandler);
    CPLHTTPResult* response = CPLHTTPFetch(url.c_str(), opts);
    CPLPopErrorHandler();
    CSLDestroy(opts);

    if (response == nullptr)
    {
        error = "no response from " + mask_url_query(url);
        return false;
    }
    if (response->nStatus != 0)
    {
        error = "request to " + mask_url_query(url) + " failed: " +
                (response->pszErrBuf ? response->pszErrBuf : "transport error");
        CPLHTTPDestroyResult(response);
        return false;
    }
    if (response->nDataLen <= 0 || response->pabyData == nullptr)
    {
        error = "empty response from " + mask_url_query(url);
        CPLHTTPDestroyResult(response);
        return false;
    }
    body.assign(reinterpret_cast<const char*>(response->pabyData),
                static_cast<std::size_t>(response->nDataLen));
    CPLHTTPDestroyResult(response);
    return true;
}

} // namespace lsdtt
