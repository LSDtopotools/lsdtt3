// lsdtt_outlet_clusters.cpp

#include "lsdtt_outlet_clusters.hpp"

#include <algorithm>
#include <cctype>
#include <map>
#include <stdexcept>

namespace lsdtt
{

namespace
{

std::string lowercase_copy(std::string text)
{
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c)
                   {
                       return static_cast<char>(std::tolower(c));
                   });
    return text;
}

} // namespace

lsdtt_outlet_cluster_mode parse_outlet_cluster_mode(const std::string& mode_text)
{
    const std::string mode = lowercase_copy(mode_text);
    if (mode == "auto") return lsdtt_outlet_cluster_mode::auto_mode;
    if (mode == "integer") return lsdtt_outlet_cluster_mode::integer_mode;
    if (mode == "string") return lsdtt_outlet_cluster_mode::string_mode;

    throw std::runtime_error(
        "outlet_cluster_mode must be one of: auto, integer, string.");
}

bool outlet_cluster_value_is_integer(const std::string& value)
{
    if (value.empty()) return false;

    std::size_t i = 0;
    if (value[0] == '-' || value[0] == '+')
    {
        i = 1;
    }
    if (i >= value.size()) return false;

    for (; i < value.size(); ++i)
    {
        if (!std::isdigit(static_cast<unsigned char>(value[i]))) return false;
    }
    return true;
}

std::string make_unit_pair_id(int unit_a, int unit_b)
{
    const int ua = std::min(unit_a, unit_b);
    const int ub = std::max(unit_a, unit_b);
    return std::to_string(ua) + "|" + std::to_string(ub);
}

std::vector<int> assign_unit_ids_from_keys(const std::vector<std::string>& unit_keys)
{
    std::map<std::string, int> key_to_unit;
    std::vector<int> unit_ids;
    unit_ids.reserve(unit_keys.size());

    int next_id = 1;
    for (const std::string& key : unit_keys)
    {
        auto it = key_to_unit.find(key);
        if (it == key_to_unit.end())
        {
            key_to_unit[key] = next_id;
            unit_ids.push_back(next_id);
            ++next_id;
        }
        else
        {
            unit_ids.push_back(it->second);
        }
    }
    return unit_ids;
}

} // namespace lsdtt
