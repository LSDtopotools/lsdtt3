#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <cpu/lsdtt_landslide_crn_cpu.hpp>
#include <cpu/lsdtt_landslide_event_cpu.hpp>
#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {
void same(const Raster& a, const Raster& b, const char* name) {
    const auto& x=a.native().get_meta(); const auto& y=b.native().get_meta();
    if(x.n_rows!=y.n_rows||x.n_cols!=y.n_cols||x.cell_size!=y.cell_size||x.x_min!=y.x_min||x.y_min!=y.y_min||x.crs_wkt!=y.crs_wkt) throw std::invalid_argument(std::string(name)+" must share the DEM grid");
}
class LandslideResult {
public:
 LandslideResult(lsdtt::lsdtt_raster_cpu z,lsdtt::lsdtt_raster_cpu t,lsdtt::lsdtt_raster_cpu d,const lsdtt::LandslideEventResult& r):z_(std::move(z)),t_(std::move(t)),d_(std::move(d)),n_(r.n_landslides),failed_(r.n_failed),requested_(r.requested_volume_m3),achieved_(r.achieved_volume_m3),summary_(r.summary()){}
 Raster elevation()const{return Raster(z_);} Raster failure_time()const{return Raster(t_);} Raster excavation_depth()const{return Raster(d_);} int n_landslides()const{return n_;} int n_failed()const{return failed_;} double requested_volume_m3()const{return requested_;} double achieved_volume_m3()const{return achieved_;} const std::string& summary()const{return summary_;}
private: lsdtt::lsdtt_raster_cpu z_,t_,d_; int n_,failed_; double requested_,achieved_; std::string summary_;
};
LandslideResult trigger_landslides(const Raster& dem,const FlowInfo& flow,double volume,double event_time,int threshold,const std::string& model,const std::string& preset,int min_pixels,double max_area,double anneal,std::uint64_t seed,const Raster* active_mask) {
    const auto& fm=flow.native().meta(); const auto& dm=dem.native().get_meta();
    if(fm.n_rows!=dm.n_rows||fm.n_cols!=dm.n_cols||fm.cell_size!=dm.cell_size||fm.x_min!=dm.x_min||fm.y_min!=dm.y_min) throw std::invalid_argument("flow must share the DEM grid");
    if(volume<=0||threshold<1||min_pixels<1) throw std::invalid_argument("target_volume_m3, channel_threshold_pixels, and min_pixels must be positive");
    if(active_mask) same(dem,*active_mask,"active_mask");
    lsdtt::LandslideSizeParams p; if(!lsdtt::parse_landslide_size_model(model,p.model)) throw std::invalid_argument("size_model must be 'inverse_gamma' or 'truncated_pareto'"); if(!lsdtt::AreaVolumeScaling::from_name(preset,p.av)) throw std::invalid_argument("unknown area_volume_preset");
    p.n_min_pixels=min_pixels; double cap=threshold*dm.cell_size*dm.cell_size; p.area_max_m2=max_area>0?max_area:cap; if(p.area_max_m2>cap) throw std::invalid_argument("max_area_m2 cannot exceed the channel-threshold area");
    lsdtt::LandslideSizeSampler sampler(p,dm.cell_size); lsdtt::LandslidePlacementConfig place; place.channel_threshold_pixels=threshold; lsdtt::LandslidePlacer placer(flow.native(),dm,place,active_mask?&active_mask->native():nullptr);
    std::vector<std::uint8_t> mask=dem.native().get_mask(); std::vector<float> times(mask.size(),static_cast<float>(dm.no_data_value)),depth(mask.size(),0.f); lsdtt::lsdtt_raster_cpu tfail(std::move(times),std::vector<std::uint8_t>(mask),dm),exc(std::move(depth),std::move(mask),dm),z=dem.native();
    lsdtt::LandslideEventConfig cfg; cfg.target_volume_m3=volume; cfg.event_time_yr=event_time; cfg.anneal_time_yr=anneal; cfg.seed=seed; auto result=lsdtt::trigger_landslide_event(z,tfail,exc,placer,sampler,cfg); return LandslideResult(std::move(z),std::move(tfail),std::move(exc),result);
}
lsdtt::LandslideCrnConfig config(double erosion,double production,double shield,const std::string& nuclide,const std::string& scheme,double now) { if(erosion<=0||production<=0||shield<=0) throw std::invalid_argument("erosion, production, and shielding must be positive"); lsdtt::LandslideCrnConfig c; c.background_eff_erate_g_cm2_yr=erosion;c.production_scaling=production;c.total_shielding=shield;c.nuclide=nuclide;c.muon_scheme=scheme;c.t_now_yr=now;return c; }
Raster landslide_crn_concentration(const Raster& depth,const Raster& failure,double erosion,double production,double shield,const std::string& nuclide,const std::string& scheme,double now){same(depth,failure,"failure_time");return Raster(lsdtt::compute_landslide_surface_concentration(depth.native(),failure.native(),config(erosion,production,shield,nuclide,scheme,now)));}
nb::dict landslide_crn_concentrations(const Raster& depth,const Raster& failure,double erosion,double production,double shield,const std::string& scheme,double now){same(depth,failure,"failure_time");auto r=lsdtt::compute_landslide_surface_concentration_multi(depth.native(),failure.native(),config(erosion,production,shield,"Be10",scheme,now));nb::dict out;out["Be10"]=Raster(std::move(r.Be10));out["Al26"]=Raster(std::move(r.Al26));out["C14"]=Raster(std::move(r.C14));return out;}
}
void bind_landslides(nb::module_& m){
 nb::class_<LandslideResult>(m,"LandslideResult",R"doc(In-memory stochastic landslide-event outputs.

Authors
-------
B. Gailleton
)doc").def_prop_ro("elevation",&LandslideResult::elevation).def_prop_ro("failure_time",&LandslideResult::failure_time).def_prop_ro("excavation_depth",&LandslideResult::excavation_depth).def_prop_ro("n_landslides",&LandslideResult::n_landslides).def_prop_ro("n_failed",&LandslideResult::n_failed).def_prop_ro("requested_volume_m3",&LandslideResult::requested_volume_m3).def_prop_ro("achieved_volume_m3",&LandslideResult::achieved_volume_m3).def_prop_ro("summary",&LandslideResult::summary);
 m.def("trigger_landslides",&trigger_landslides,nb::arg("dem"),nb::arg("flow"),nb::arg("target_volume_m3"),nb::arg("event_time_yr")=0.,nb::arg("channel_threshold_pixels")=1000,nb::arg("size_model")="inverse_gamma",nb::arg("area_volume_preset")="larsen_soil_scar",nb::arg("min_pixels")=1,nb::arg("max_area_m2")=0.,nb::arg("anneal_time_yr")=0.,nb::arg("seed")=1u,nb::arg("active_mask")=nullptr,R"doc(Trigger a volume-budgeted stochastic landslide event.

Parameters
----------
dem : Raster
    Depression-conditioned elevation raster.
flow : FlowInfo
    Same-grid D8 routing topology.
target_volume_m3 : float
    Requested total scar volume.
event_time_yr, channel_threshold_pixels, size_model, area_volume_preset
    Event timing and native placement/size controls.
min_pixels, max_area_m2, anneal_time_yr, seed, active_mask
    Scar bounds, recovery, reproducibility, and optional host restriction.

Returns
-------
LandslideResult
    Excavated surface, failure-time/depth rasters, and budget diagnostics.

Authors
-------
B. Gailleton
)doc");
 m.def("landslide_crn_concentration",&landslide_crn_concentration,nb::arg("excavation_depth"),nb::arg("failure_time"),nb::arg("background_erosion_rate_g_cm2_yr")=.001,nb::arg("production_scaling")=1.,nb::arg("total_shielding")=1.,nb::arg("nuclide")="Be10",nb::arg("muon_scheme")="BraucherBorchers",nb::arg("report_time_yr")=0.,R"doc(Compute one-nuclide surface concentration after landslide excavation.

Parameters
----------
excavation_depth, failure_time : Raster
    Same-grid event rasters from ``trigger_landslides``.
background_erosion_rate_g_cm2_yr, production_scaling, total_shielding : float
    Background CRN conditions.
nuclide, muon_scheme, report_time_yr
    Nuclide, production scheme, and absolute reporting time.

Returns
-------
Raster
    Surface concentration in atoms g-1.

Authors
-------
B. Gailleton
)doc");
 m.def("landslide_crn_concentrations",&landslide_crn_concentrations,nb::arg("excavation_depth"),nb::arg("failure_time"),nb::arg("background_erosion_rate_g_cm2_yr")=.001,nb::arg("production_scaling")=1.,nb::arg("total_shielding")=1.,nb::arg("muon_scheme")="BraucherBorchers",nb::arg("report_time_yr")=0.,R"doc(Compute Be10, Al26, and C14 concentrations after landslide excavation.

Parameters
----------
excavation_depth, failure_time : Raster
    Same-grid event rasters.
background_erosion_rate_g_cm2_yr, production_scaling, total_shielding, muon_scheme, report_time_yr
    Background CRN conditions.

Returns
-------
dict
    Concentration rasters keyed by nuclide.

Authors
-------
B. Gailleton
)doc");
}
} // namespace lsdtt3_python
