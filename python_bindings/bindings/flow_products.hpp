#pragma once

#include "binding_common.hpp"

namespace lsdtt3_python {
void bind_flow_products(nb::class_<FlowInfo>& flow);
void bind_flow_network(nb::module_& module, nb::class_<FlowInfo>& flow);
}
