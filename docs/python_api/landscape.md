# Landscape And Surface Workflows

## Landscape evolution

::: lsdtt3._api_docs.LandscapeModel

`incision_threshold_specific_area` is measured as `A/dx` in metres. The former
`incision_threshold_area` spelling is deprecated but retains the same value
during its compatibility window.

## Flat-surface detection

::: lsdtt3._api_docs.find_flat_surfaces

`relief_window_side_pixels` is the full odd side length of the local window, not
a radius. `relief_window_pixels` remains a deprecated compatibility keyword.

## Valley centrelines

::: lsdtt3._api_docs.extract_valley_centrelines

Valley identification belongs upstream. Supply `valley_mask` explicitly; the
centreline operation does not derive a HAND break. A classification result's
`floodplain_mask` can be passed directly when that classification is the chosen
in-memory valley definition.
