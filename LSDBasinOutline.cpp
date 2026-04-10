//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDBasinOutline
// Land Surface Dynamics Basin Outline
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for extracting various metrics related to basin outline.
//  This object differs from LSDBasin in that (a) drainage divides are defined
//  in-between DEM pixels and (b) outline data storage is structured in a way
//  to enable analysis of ridgelines.
//
// Developed by:
//  Rutger Vervoordeldonk
//  Simon M. Mudd
//
// Copyright (C) 2025 Rutger Vervoordeldonk 2025
//
// Lead developer can be contacted by simon.m.mudd _at_ ed.ac.uk
//
//    Simon Mudd
//    University of Edinburgh
//    School of GeoSciences
//    Drummond Street
//    Edinburgh, EH8 9XP
//    Scotland
//    United Kingdom
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the
// GNU General Public License along with this program;
// if not, write to:
// Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor,
// Boston, MA 02110-1301
// USA
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDBasinOutline.cpp
// cpp file for the LSDBasinOutline object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Rutger Vervoordeldonk, University of Edinburgh
// NOTE 27/6/25: this object is to replace LSDBasinOutlines, which holds
// outlines for many basins and is therefore confusing. This object only holds
// information for one basin at a time.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <algorithm>
#include <set>
#include <cmath>
#include <ctime>
#include <numbers>
#include <limits>
#include "TNT/tnt.h"
#include "LSDBasinOutline.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDBasin.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDRasterInfo.hpp"
#include "LSDRasterSpectral.hpp"
#include "LSDGeometry.hpp"
#include "LSDChiTools.hpp"
#include "LSDStatsTools.hpp"
#include "LSDRidgeline.hpp"
#include "fftw-3.3.4/api/fftw3.h"
using namespace std;
using namespace TNT;

LSDBasinOutline& LSDBasinOutline::operator=(const LSDBasinOutline& rhs){
    if (&rhs != this){
        NRows = rhs.NRows;
        NCols = rhs.NCols;
        NDataNodes = rhs.NDataNodes;
        XMinimum = rhs.XMinimum;
        YMinimum = rhs.YMinimum;

        NRowsDiv = rhs.NRowsDiv;
        NColsDiv = rhs.NColsDiv;
        NDataNodesDiv = rhs.NDataNodesDiv;
        XMinimumDiv = rhs.XMinimumDiv;
        YMinimumDiv = rhs.YMinimumDiv;

        boundary_conditions = rhs.boundary_conditions;
        is_periodic_EW = rhs.is_periodic_EW;
        is_periodic_NS = rhs.is_periodic_NS;

        GeoReferencingStrings = rhs.GeoReferencingStrings;

        DataResolution = rhs.DataResolution;
        NoDataValue = rhs.NoDataValue;

        basin_junction_number = rhs.basin_junction_number;
        basin_outlet_node = rhs.basin_outlet_node;

        outline_nodes = rhs.outline_nodes;
        outline_cols = rhs.outline_cols;
        outline_rows = rhs.outline_rows;

        outline_distance_cw = rhs.outline_distance_cw;
        outline_distance_ccw = rhs.outline_distance_ccw;
        outline_elevation = rhs.outline_elevation;
    }
    return *this;
}


void LSDBasinOutline::create(void){
    /*
    vector<vector<int>> emptyvec_two_levels;
    vector<int> emptyvec_one_level;
    vector<int> emptyvec_float;

    basin_outlines_nodes = emptyvec_two_levels;
    basin_outlines_rows = emptyvec_two_levels;
    basin_outlines_cols = emptyvec_two_levels;

    basin
    */
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDBasinOutline::create(LSDFlowInfo& FI, int& junction_number, 
    int& outlet_node){
    // Store some info
    basin_junction_number = junction_number;
    basin_outlet_node = outlet_node;

    // Store grid information
    NRows = FI.get_NRows();
    NCols = FI.get_NCols(); 
    NDataNodes = FI.get_NDataNodes();
    XMinimum = FI.get_XMinimum();
    YMinimum = FI.get_YMinimum();
    DataResolution = FI.get_DataResolution();
    GeoReferencingStrings = FI.get_GeoReferencingStrings();
    NoDataValue = FI.get_NoDataValue();
    boundary_conditions = FI.get_BoundaryConditions();

    NRowsDiv = NRows+1;
    NColsDiv = NCols+1;
    NDataNodesDiv = NRowsDiv*NColsDiv;
    XMinimumDiv = XMinimum - DataResolution/2;
    YMinimumDiv = YMinimum - DataResolution/2;

    // Check for periodic boundary conditions
    is_periodic_NS = false;
    is_periodic_EW = false;

    if(boundary_conditions.size() == 4){
        if (boundary_conditions[0] == "p" && boundary_conditions[2] == "p"){
            is_periodic_NS = true;
        }
        if (boundary_conditions[1] == "p" && boundary_conditions[3] == "p"){
            is_periodic_EW = true;
        }
    }

    if (is_periodic_NS && is_periodic_EW){
        cout << "LSDBasinOutline::create: both NS and EW boundaries periodic!" << endl;
        cout << "Can't deal with that. Turning off periodicity NS boundaries." << endl;
        is_periodic_NS = false;
    }

    // These functions will deal with edge tracing and result in filled 
    // outline_rows and outline_cols vectors that describe the basin outline.
    if (is_periodic_EW || is_periodic_NS){
        edge_tracing_periodic_boundaries(FI);
    }
    else {
        edge_tracing_normal_boundaries(FI);
    }

    // At this point, the vectors just start at the top left due to edge tracing
    // but we want it to start at the basin outlet.

    // Get the linear row col indices of the basin outlet.
    int outlet_row, outlet_col;
    FI.retrieve_current_row_and_col(outlet_node, outlet_row, outlet_col);

    // Get the pixel immediately downstram of the basin outlet.
    // This pixel always exists provided we are dealing with a complete basin.
    int ds_outlet_node, ds_outlet_row, ds_outlet_col;
    FI.retrieve_receiver_information(outlet_node, ds_outlet_node, 
        ds_outlet_row, ds_outlet_col);

    // Get the corners of the outlet pixel that mark the start and the end of
    // the edge-based basin outline.
    int startrowDiv, startcolDiv;
    int endrowDiv, endcolDiv;
    int rowchange = ds_outlet_row - outlet_row;
    int colchange = ds_outlet_col - outlet_col;

    // There are three options:
    // * we are at a baselevel (on the DEM edge, if no holes in DEM). 
    //      -> start and end on the outlet pixel corners that are on the
    //         DEM edge. 
    // * we have diagonal flow going out of the basin. 
    //      -> the start and end pixels corners are the same.
    // * we have horizontal or vertical flow out of the basin
    //      -> start on the two adjacent pixel corners next to where the
    //         flow is going.   

    // 1) Baselevel:
    if (rowchange == 0 && colchange == 0) {
        // Let's find the pixel corners on the DEM edge
        // We do this by iterating through the neighbours and identifying
        // when we change from a nodatavalue pixel (off the DEM) to a 
        // datavalue pixel (on the DEM). Nodatavalue either when there 
        // is no data (duh) or when we are out-of-bounds of the DEM.
        vector<pair<int,int>> neighbor_offsets = {
            {-1,-1},{-1, 0},{-1, 1},
            { 0, 1},{ 1, 1},{ 1, 0},
            { 1,-1},{ 0,-1}
        };
        for (size_t i = 0; i < neighbor_offsets.size(); ++i) {
            // neighbor 1
            int dr1 = neighbor_offsets[i].first;
            int dc1 = neighbor_offsets[i].second;
            int n1r_raw = outlet_row + dr1;
            int n1c_raw = outlet_col + dc1;
            int n1r = wrap_index(n1r_raw, NRows, is_periodic_NS);
            int n1c = wrap_index(n1c_raw, NCols, is_periodic_EW);
            bool nodata1 = (n1r == NoDataValue || n1c == NoDataValue)
                        || FI.get_NodeIndex_from_row_col(n1r, n1c) == NoDataValue;


            // neighbor 2 (next in CW)
            size_t ni = (i+1) % neighbor_offsets.size();
            int dr2 = neighbor_offsets[ni].first;
            int dc2 = neighbor_offsets[ni].second;
            int n2r_raw = outlet_row + dr2;
            int n2c_raw = outlet_col + dc2;
            int n2r = wrap_index(n2r_raw, NRows, is_periodic_NS);
            int n2c = wrap_index(n2c_raw, NCols, is_periodic_EW);
            bool nodata2 = (n2r == NoDataValue || n2c == NoDataValue)
                        || FI.get_NodeIndex_from_row_col(n2r, n2c) == NoDataValue;

            // entry point onto DEM
            if (nodata1 && !nodata2) {
                int sum_dr = dr1 + dr2;
                int sum_dc = dc1 + dc2;
                int raw_row = outlet_row + (sum_dr >= 0 ? 1 : 0);
                int raw_col = outlet_col + (sum_dc >= 0 ? 1 : 0);
                startrowDiv = wrap_index_divgrid(raw_row, NRowsDiv, is_periodic_NS);
                startcolDiv = wrap_index_divgrid(raw_col, NColsDiv, is_periodic_EW);
            }

            // exit point off DEM
            if (!nodata1 && nodata2) {
                int sum_dr = dr1 + dr2;
                int sum_dc = dc1 + dc2;
                int raw_row = outlet_row + (sum_dr >= 0 ? 1 : 0);
                int raw_col = outlet_col + (sum_dc >= 0 ? 1 : 0);
                endrowDiv = wrap_index_divgrid(raw_row, NRowsDiv, is_periodic_NS);
                endcolDiv = wrap_index_divgrid(raw_col, NColsDiv, is_periodic_EW);
            }

        }
    }

    // 2) Diagonal flow
    else if (rowchange != 0 && colchange != 0) {
        int dr = (rowchange > 0 ? 1 : 0);
        int dc = (colchange > 0 ? 1 : 0);
        startrowDiv = endrowDiv = wrap_index_divgrid(outlet_row + dr, NRowsDiv, is_periodic_NS);
        startcolDiv = endcolDiv = wrap_index_divgrid(outlet_col + dc, NColsDiv, is_periodic_EW);
    }
    // 3) Horizontal flow (two adjacent corners on E or W side)
    else if (rowchange == 0 && colchange != 0) {
        if (colchange > 0) {
            // flow → E : corners SE then NE
            startrowDiv = wrap_index_divgrid(outlet_row + 1, NRowsDiv, is_periodic_NS);
            startcolDiv = wrap_index_divgrid(outlet_col + 1, NColsDiv, is_periodic_EW);
            endrowDiv   = wrap_index_divgrid(outlet_row, NRowsDiv, is_periodic_NS);
            endcolDiv   = wrap_index_divgrid(outlet_col + 1, NColsDiv, is_periodic_EW);
        } else {
            // flow → W : corners NW then SW
            startrowDiv = wrap_index_divgrid(outlet_row, NRowsDiv, is_periodic_NS);
            startcolDiv = wrap_index_divgrid(outlet_col, NColsDiv, is_periodic_EW);
            endrowDiv   = wrap_index_divgrid(outlet_row + 1, NRowsDiv, is_periodic_NS);
            endcolDiv   = wrap_index_divgrid(outlet_col, NColsDiv, is_periodic_EW);
        }
    }
    // 4) Vertical flow (two adjacent corners on N or S side)
    else if (rowchange != 0 && colchange == 0) {
        if (rowchange > 0) {
            // flow → S : corners SW then SE
            startrowDiv = wrap_index_divgrid(outlet_row + 1, NRowsDiv, is_periodic_NS);
            startcolDiv = wrap_index_divgrid(outlet_col, NColsDiv, is_periodic_EW);
            endrowDiv   = wrap_index_divgrid(outlet_row + 1, NRowsDiv, is_periodic_NS);
            endcolDiv   = wrap_index_divgrid(outlet_col + 1, NColsDiv, is_periodic_EW);
        } else {
            // flow → N : corners NE then NW
            startrowDiv = wrap_index_divgrid(outlet_row, NRowsDiv, is_periodic_NS);
            startcolDiv = wrap_index_divgrid(outlet_col + 1, NColsDiv, is_periodic_EW);
            endrowDiv   = wrap_index_divgrid(outlet_row, NRowsDiv, is_periodic_NS);
            endcolDiv   = wrap_index_divgrid(outlet_col, NColsDiv, is_periodic_EW);
        }
    }

    // Locate start and end in the old outline
    int start_idx = -1;
    int end_idx   = -1;
    int M = outline_rows.size();
    for (int i = 0; i < M; ++i) {
        if (outline_rows[i] == startrowDiv && outline_cols[i] == startcolDiv) {
            start_idx = i;
        }
        if (outline_rows[i] == endrowDiv   && outline_cols[i] == endcolDiv) {
            end_idx = i;
        }
    }

    // Sanity check
    if (start_idx < 0 || end_idx < 0) {
        cout << endl << "ERROR LSDBasinOutline::create: could not find start or end corner in outline" << endl;
        cout << "Basin junction number is: " << basin_junction_number << endl;
        
        cout << "The start_idx is: " << start_idx << endl;
        cout << "The end_idx is: " << end_idx << endl;
        cout << "The start row/col (divide grid) is: " << startrowDiv << "/" << startcolDiv << endl;
        cout << "The end row/col (divide grid) is: " << endrowDiv << "/" << endcolDiv<< endl;
        cout << "The outline consists of row/col coordinates:" << endl;
        
        for(size_t i=0; i<outline_rows.size(); i++){
            cout << outline_rows[i] << "/" << outline_cols[i] << endl;
        }
        
        
    }

    // Build a new outline that walks from start → end (inclusive)
    vector<int> new_rows, new_cols;
    new_rows.reserve(M);   // at most M points
    new_cols.reserve(M);

    // Diagonal case: start==end → do a full cycle then repeat the start corner
    if (start_idx == end_idx) {
        for (int k = 0; k < M; ++k) {
            int idx = (start_idx + k) % M;
            new_rows.push_back(outline_rows[idx]);
            new_cols.push_back(outline_cols[idx]);
        }
        // close the loop by repeating the start corner
        new_rows.push_back(outline_rows[start_idx]);
        new_cols.push_back(outline_cols[start_idx]);
    }
    // General case: walk from start → end (inclusive)
    else {
        // General case: walk from start → end (inclusive), wrapping if needed
        int i = start_idx; 
        do { 
            new_rows.push_back(outline_rows[i]); 
            new_cols.push_back(outline_cols[i]); 
            if (i == end_idx) break; 
            i = (i + 1) % M; 
        } while (true);
    }

    // Swap old and new
    outline_rows.swap(new_rows);
    outline_cols.swap(new_cols);

    // Store number of nodes in outline
    NNodes_outline = outline_rows.size();

    // Get the outline in terms of node indices
    for(int i=0; i<NNodes_outline; i++){
        outline_nodes.push_back(
            retrieve_node_from_row_and_column_dividegrid(
                outline_rows[i], outline_cols[i]
            )
        );
    }

    // Get the distance along the outline (clockwise and counter-clockwise).
    float starting_distance=0;
    outline_distance_cw.resize(NNodes_outline);
    outline_distance_ccw.resize(NNodes_outline);
    for(int i=0; i<NNodes_outline; i++){
        // The distance is then relatively trivial as all outline nodes
        // are 4-connected along px edges the size of DataResolution.
        outline_distance_cw[i] = i*DataResolution + starting_distance;
        outline_distance_ccw[NNodes_outline-1-i] = 
            i*DataResolution + starting_distance;
    }   
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Function called by create to do edge tracing on basin in non-periodic 
// grid.
//
// RV 15/09/2025
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDBasinOutline::edge_tracing_normal_boundaries(LSDFlowInfo& FI){
    // Get all the nodes in this basin
    vector<int> BasinNodes = FI.get_upslope_nodes(basin_outlet_node);

    // Let's get the outline of this basin
    // The following lines are identical to LSDBasin's get_edge_boundary 
    // function, but copying the lines of code in here avoids making
    // an LSDBasin, which requires a JN.
    
    // Get the row and colum positions of all basin nodes, store in r and c 
    // vectors and find out what the max/min row/col indices are
    vector<int> r,c;
    r.reserve(BasinNodes.size());
    c.reserve(BasinNodes.size());

    int min_r = NRows, max_r = 0;
    int min_c = NCols, max_c = 0;
    for (int node : BasinNodes) {
        int row, col;
        FI.retrieve_current_row_and_col(node, row, col);
        r.push_back(row);
        c.push_back(col);
        if (row < min_r) min_r = row;
        if (row > max_r) max_r = row;
        if (col < min_c) min_c = col;
        if (col > max_c) max_c = col;
    }

    // Expand max/min row/col indices by 1 pixel each direction
    min_r = min_r-1;
    max_r = max_r+1;
    min_c = min_c-1;
    max_c = max_c+1;
  
    // Build double resolution binary grid. The double resolution grid will
    // have 9 pixels associated with each DEM pixel. 4 of these lie on the DEM
    // pixel corners, and that is what we want for edge-based boundaries.
    // This array is 'cropped' around the basin (based on the max/min row/col
    // indices)
    

    int NRowsDouble = (max_r - min_r +1) * 2 + 1;
    int NColsDouble = (max_c - min_c +1) * 2 + 1;
    Array2D<int> binary_array(NRowsDouble, NColsDouble, 0);

    const int offsets[9][2] = {
        {-1, -1}, {-1, 0}, {-1, 1},
        { 0, -1}, { 0, 0}, { 0, 1},
        { 1, -1}, { 1, 0}, { 1, 1}
    };

    for (size_t k = 0; k < r.size(); k++){
        int base_r = (r[k] - min_r) * 2;
        int base_c = (c[k] - min_c) * 2;
        for (int i = 0; i < 9; i++){
        int dr = base_r + offsets[i][0];
        int dc = base_c + offsets[i][1];


        binary_array[dr][dc] = 1; 
        }
    }

    // Now run the edge tracing algorithm.
    // This edge tracing was first coded up in the LSDIndexRaster, so need to 
    // make a dummy indexraster to access the tracing function.
    LSDIndexRaster dummy_raster;
    pair<vector<int>,vector<int>> boundarypixels =  
        dummy_raster.trace_boundary_pixels(binary_array, 4);
    vector<int> rr_boundarypixels = boundarypixels.first;
    vector<int> cc_boundarypixels = boundarypixels.second;

    // Discard of all the coordinates that are NOT on the corners of the 
    // original pixels. Coordinate pairs on the corners will have both row and
    // column coordinate being odd owing to the 0.5*2 offset from original 
    // pixel in both row and col directions.
    outline_rows.reserve(rr_boundarypixels.size());
    outline_cols.reserve(cc_boundarypixels.size());
    for(size_t k = 0; k < rr_boundarypixels.size(); k++){
        int rr = rr_boundarypixels[k];
        int cc = cc_boundarypixels[k];
        if ((rr & 1) && (cc & 1)){
        outline_rows.push_back(((rr+1) >> 1)+min_r);
        outline_cols.push_back(((cc+1) >> 1)+min_c);
        }
    }
    
    // Get rid of the final node since that is the same as the first node.
    outline_rows.pop_back();
    outline_cols.pop_back();
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Function called by create to do edge tracing on basin in non-periodic 
// grid.
//
// RV 15/09/2025
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDBasinOutline::edge_tracing_periodic_boundaries(LSDFlowInfo& FI){
    // Overall strategy of tracing periodic boundaries is by first making a 
    // binary array that contains the "unwrapped" basin.
    // This is done by iterating up the basin stack, and every time we cross 
    // a periodic boundary in one direction we move to an "adjacent" grid, and
    // if we cross the periodic boundary in the other direction, we move back 
    // to the original grid. 
    // For each node we know its position by storing its row index, column index
    // and the grid_index. The grid_index of the outlet is 0.

    // Get the nodes upslope of the outlet node
    // This is in the stacking order
    vector<int> BasinNodes = FI.get_upslope_nodes(basin_outlet_node);
    
    // Get the row, column and grid_index positions of all basin nodes
    vector<int> r,c;
    r.reserve(BasinNodes.size());
    c.reserve(BasinNodes.size());
    vector<int> grid_index(BasinNodes.size(), -1);
    grid_index[0] = 0; // grid index of outlet is 0.
    int first_SVector_index = FI.get_index_SVector(basin_outlet_node);

    for(int i=0; i<BasinNodes.size(); i++){
        // Get the current row and col and append
        int thisNode = BasinNodes[i];
        int thisRow, thisCol;
        FI.retrieve_current_row_and_col(thisNode, thisRow, thisCol);
        r.push_back(thisRow);
        c.push_back(thisCol);

        // The index of the receiver node in the BasinNodes vector
        // (note the BasinNodes vector is essentially a slide of the stack or
        // SVector in LSDFlowInfo)
        int receiverNode;
        FI.retrieve_receiver_information(thisNode, receiverNode);
        int receiverIndex = FI.get_index_SVector(receiverNode)-first_SVector_index;

        // Look up the grid_index of the receiver 
        int thisGridIndex = grid_index[receiverIndex];

        // If we cross a periodic boundary, move to adjacent grid by amending
        // thisGridIndex
        int flow_dir = FI.get_LocalFlowDirection(thisRow, thisCol);
        if (is_periodic_EW &&  thisCol == NCols-1){
            if (flow_dir == 1 || flow_dir == 2 || flow_dir == 3){
                // on periodic East boundary and flow to East, so
                // we are on -1 grid index compared to its receiver.
                thisGridIndex -= 1;
            }
        }
        if (is_periodic_EW && thisCol == 0){
            if (flow_dir == 5 || flow_dir == 6 || flow_dir == 7){
                // on periodic West boundary and flow to West, so
                // we are on +1 grid index compared to its receiver.
                thisGridIndex += 1;
            }
        }
        if (is_periodic_NS && thisRow == NRows-1){
            if (flow_dir == 3 || flow_dir == 4 || flow_dir == 5){
                // on periodic South boundary and flow to South, so
                // we are on +1 grid index compared to its receiver.
                thisGridIndex += 1;
            }
        }
        if (is_periodic_NS && thisRow == 0){
            if (flow_dir == 7 || flow_dir == 0 || flow_dir == 1){
                // on periodic North boundary and flow to North, so
                // we are on -1 grid index compared to its receiver.
                thisGridIndex -= 1;
            }
        }
        grid_index[i] = thisGridIndex; 
    }

    // Now we know for each node on "which grid" it is located
    // and what it row and column on that grid is.
    // Get the min and max grid indices.
    int min_grid = *min_element(grid_index.begin(), grid_index.end());
    int max_grid = *max_element(grid_index.begin(), grid_index.end());

    // Determine the unwrapped extents 
    // Coordinate mapping depends on boundary conditions.
    // For EW periodic:
    //   Each col in the copy sits at c[i] + (grid_index[i]-min_grid)*NCols
    //   Each row sits at original row
    // For NS periodic:
    //   Each row in the copy sits at r[i] + (grid_index[i]-min_grid)*NRows
    //   Each col sits at original row
    // Compute unwrapped min/max row and col indices
    int unwrapped_min_row = numeric_limits<int>::max();
    int unwrapped_max_row = numeric_limits<int>::min();
    int unwrapped_min_col = numeric_limits<int>::max();
    int unwrapped_max_col = numeric_limits<int>::min();

    for (size_t i = 0; i < BasinNodes.size(); ++i) {
        int row = r[i];
        int col = c[i];
        int gidx = grid_index[i];

        int urow = is_periodic_NS ? row + (gidx - min_grid) * NRows : row;
        int ucol = is_periodic_EW ? col + (gidx - min_grid) * NCols : col;

        unwrapped_min_row = min(unwrapped_min_row, urow);
        unwrapped_max_row = max(unwrapped_max_row, urow);
        unwrapped_min_col = min(unwrapped_min_col, ucol);
        unwrapped_max_col = max(unwrapped_max_col, ucol);
    }

    // Expand by 1 DEM pixel in every direction
    int r0, r1, c0, c1;
    r0 = unwrapped_min_row - 1;
    r1 = unwrapped_max_row + 1;
    c0 = unwrapped_min_col - 1;
    c1 = unwrapped_max_col + 1;


    // Build the double-resolution grid
    int nR = (r1 - r0 + 1) * 2 + 1;
    int nC = (c1 - c0 + 1) * 2 + 1;
    Array2D<int> binary2(nR, nC, 0);

    static const int offs[9][2] = {
        {-1, -1}, {-1,  0}, {-1,  1},
        { 0, -1}, { 0,  0}, { 0,  1},
        { 1, -1}, { 1,  0}, { 1,  1}
    };

    for (size_t i = 0; i < BasinNodes.size(); ++i) {
        int base_r, base_c;
        if (is_periodic_EW) {
            int uR = r[i];
            int uC = c[i] + (grid_index[i] - min_grid) * NCols;
            base_r = (uR - r0) * 2;
            base_c = (uC - c0) * 2;
        }
        else { // NS periodic
            int uR = r[i] + (grid_index[i] - min_grid) * NRows;
            int uC = c[i];
            base_r = (uR - r0) * 2;
            base_c = (uC - c0) * 2;
        }

        for (int k = 0; k < 9; ++k) {
            int dr = base_r + offs[k][0];
            int dc = base_c + offs[k][1];
            if (dr >= 0 && dr < nR && dc >= 0 && dc < nC) {
            binary2[dr][dc] = 1;
            }
        }
    }

    // Run the edge-tracing algorithm
    LSDIndexRaster dummy_raster;
    pair<vector<int>,vector<int>> boundarypixels =  
        dummy_raster.trace_boundary_pixels(binary2, 4);
    vector<int> rr_boundarypixels = boundarypixels.first;
    vector<int> cc_boundarypixels = boundarypixels.second;

    // Set divide grid dimensions
    int NRowsDiv = NRows + (is_periodic_NS  ? 0 : 1);
    int NColsDiv = NCols + (is_periodic_EW  ? 0 : 1);

    // Prepare divide-grid based outline rows and cols.
    outline_rows.reserve(rr_boundarypixels.size());
    outline_cols.reserve(cc_boundarypixels.size());

    // Helper for wrapping back to a periodic divide grid
    auto wrap = [&](int idx, int M) {
        int m = idx % M;
        return (m < 0 ? m + M : m);
    };

    for (size_t k = 0; k < rr_boundarypixels.size(); ++k) {
        int rr = rr_boundarypixels[k];
        int cc = cc_boundarypixels[k];

        // keep only the true pixel corners  
        if ((rr & 1) && (cc & 1)) {
            // unwrapped cell‐corner coordinates  
            int corner_r = ((rr + 1) >> 1) + r0;
            int corner_c = ((cc + 1) >> 1) + c0;

            // fold back onto original periodic domain in each direction  
            if (is_periodic_NS) corner_r = wrap(corner_r, NRowsDiv);
            if (is_periodic_EW) corner_c = wrap(corner_c, NColsDiv);

            outline_rows.push_back(corner_r);
            outline_cols.push_back(corner_c);
        }
    }
    
    // Drop the last point (duplicate of first)
    outline_rows.pop_back();
    outline_cols.pop_back();

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Function to retrieve the node from row and column of dividegrid. Returns int
//
// RV 19/02/2025
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDBasinOutline::retrieve_node_from_row_and_column_dividegrid(int row, 
    int col){
    int node = row * NColsDiv + col; 
    return node;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Function to retrieve the row and column from a node in the divide grid.
// Passes on row/col on to empty ints that you need to provide.
//
// RV 19/02/2025
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDBasinOutline::retrieve_row_and_column_from_node_dividegrid(int node, 
int& row, int& col){
    int cr, cc;
    cr = node / NColsDiv;
    cc = node % NColsDiv;
    row = cr;
    col = cc; 
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// A simple helper function to wrap a DEM index / flow grid index in case 
// boundaries are periodic.
//
// RV 19/02/2025
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDBasinOutline::wrap_index(int idx, int max, bool is_periodic) {
    if (is_periodic) {
        if (idx < 0) return max - 1;
        if (idx >= max) return 0;
    } else {
        if (idx < 0 || idx >= max) return NoDataValue;
    }
    return idx;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// A simple helper function to wrap a divide grid index in case 
// boundaries are periodic.
//
// RV 19/02/2025
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDBasinOutline::wrap_index_divgrid(int idx, int max, bool is_periodic) {
    if (is_periodic) {
        if (idx < 0) return max - 2;
        if (idx >= max-1) return 0;
    } else {
        if (idx < 0 || idx >= max) return NoDataValue;
    }
    return idx;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Function to convert a node position in the divide grid with a row and column 
// to  lat/long coordinate. Taken from SMM's 
// LSDJunctionNetwork::get_lat_and_long_locations.
//
// RV 21/02/2025
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDBasinOutline::get_lat_and_long_locations_dividegrid(int dividerow, 
    int dividecol, double& latitude, double& longitude, 
    LSDCoordinateConverterLLandUTM Converter, LSDJunctionNetwork JN){
    // get the x and y locations of the node
    double x_loc, y_loc;
    x_loc = XMinimumDiv + float(dividecol)*DataResolution + 0.5*DataResolution;
    y_loc = YMinimumDiv + float(NRowsDiv-dividerow)*DataResolution-0.5*DataResolution;

    // get the UTM zone of the node
    int UTM_zone;
    bool is_north;
    JN.get_UTM_information(UTM_zone, is_north);

    if(UTM_zone == NoDataValue){
        latitude = NoDataValue;
        longitude = NoDataValue;
    } else {
        // set the default ellipsoid to WGS84
        int eId = 22;

        double xld = double(x_loc);
        double yld = double(y_loc);

        // use the converter to convert to lat and long
        double Lat,Long;
        Converter.UTMtoLL(eId, yld, xld, UTM_zone, is_north, Lat, Long);

        latitude = Lat;
        longitude = Long;
    }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Get lats and longs for basin outline nodes.
//
// RV 08/12/2025
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDBasinOutline::get_latitudes_longitudes_basin_outline(vector<double>& latitudes, vector<double>& longitudes, LSDJunctionNetwork& JN){
    // Get the converter
    LSDCoordinateConverterLLandUTM Converter;

    // Clear the supplied lists
    latitudes.clear();
    longitudes.clear();

    // Loop through each outline node and push into lats and longs
    for(int i=0; i<NNodes_outline; i++){
        // Get the node
        int divcol = outline_cols[i];
        int divrow = outline_rows[i];

        // Call latlong function
        double this_latitude, this_longitude;
        get_lat_and_long_locations_dividegrid(divrow, divcol, this_latitude, this_longitude, Converter, JN);

        // Push back lat long
        latitudes.push_back(this_latitude);
        longitudes.push_back(this_longitude);
    }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Function to interpolate the four pixel values (of a supplied LSDRaster in
// the reference of the flow grid) that surround each node in a basin outline.
//
// RV 27/06/2025 based on similar function LSDDivides by RV 19/02/2025
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<float> LSDBasinOutline::interpolate_pixel_values_on_outline(
    LSDRaster& data_raster){
    
    vector<float> interpolated_values(NNodes_outline, NoDataValue);
    // Loop through outline nodes
    for(int i=0; i<NNodes_outline; i++){
        // Get the node
        int divcol = outline_cols[i];
        int divrow = outline_rows[i];
        //cout << divrow << "/" << divcol << endl;

        // Get the x and y for this node in the reference of the flow grid
        float x_interp = divcol-0.5;
        float y_interp = divrow-0.5;

        // Get the x's and y's indices of the 4 pixels of the flow grid
        int x1 = static_cast<int>(floor(x_interp));
        int y1 = static_cast<int>(floor(y_interp));
        int x2 = x1 + 1;
        int y2 = y1 + 1;

        // Deal with periodic boundaries and grid edges
        int NColsDataRaster = data_raster.get_NCols();
        if(is_periodic_EW){
            // Periodic east to west and west to east!
            // If the lower column index is -1, it should point to the final col
            if(x1 == -1){x1 = NColsDataRaster;}
            // If the higher column index is final col+1, it should point to col0
            if(x2 == NColsDataRaster){x2=0;}
        } else {
            // Closed boundaries. Set non-existing coordinate to NoDataValue.
            if(x1 < 0){x1 = NoDataValue;}
            if(x2 >= NColsDataRaster){x2 = NoDataValue;}
        }
        int NRowsDataRaster = data_raster.get_NRows();
        if(is_periodic_NS){
            // Periodic north to south and south to north!
            // If the lower row index is -1, it should point to the final row
            if(y1 == -1){y1 = NRowsDataRaster;}
            // If the higher row index is final row+1, it should point to row0
            if(y2 == NRowsDataRaster){y2=0;}
        } else {
            // Closed boundaries. Set non-existing coordinate to NoDataValue.
            if(y1 < 0){y1 = NoDataValue;}
            if(y2 >= NRowsDataRaster){y2 = NoDataValue;}
        }

        // Initialize corner values
        float Q11 = NoDataValue, Q12 = NoDataValue, Q21 = NoDataValue, Q22 = NoDataValue;

        // Flags to track validity
        bool valid_Q11 = (x1 != NoDataValue && y1 != NoDataValue);
        bool valid_Q12 = (x1 != NoDataValue && y2 != NoDataValue);
        bool valid_Q21 = (x2 != NoDataValue && y1 != NoDataValue);
        bool valid_Q22 = (x2 != NoDataValue && y2 != NoDataValue);

        // Safely access data only if coordinates are valid
        if(valid_Q11) Q11 = data_raster.get_data_element(y1, x1);
        if(valid_Q12) Q12 = data_raster.get_data_element(y2, x1);
        if(valid_Q21) Q21 = data_raster.get_data_element(y1, x2);
        if(valid_Q22) Q22 = data_raster.get_data_element(y2, x2);

        // Re-check validity based on actual data values
        valid_Q11 = valid_Q11 && Q11 != NoDataValue;
        valid_Q12 = valid_Q12 && Q12 != NoDataValue;
        valid_Q21 = valid_Q21 && Q21 != NoDataValue;
        valid_Q22 = valid_Q22 && Q22 != NoDataValue;

        int valid_count = valid_Q11 + valid_Q12 + valid_Q21 + valid_Q22;

        if(valid_count == 4){
            // Bilinear interpolation
            float fxy1 = (Q11*(x2 - x_interp) + Q21*(x_interp - x1)) / (x2 - x1);
            float fxy2 = (Q12*(x2 - x_interp) + Q22*(x_interp - x1)) / (x2 - x1);
            float fxy = (fxy1*(y2 - y_interp) + fxy2*(y_interp - y1)) / (y2 - y1);
            interpolated_values[i] = fxy;
        }
        else if(valid_count == 3){
            // Average the three valid values
            float sum = 0.0;
            int count = 0;
            if(valid_Q11){ sum += Q11; count++; }
            if(valid_Q12){ sum += Q12; count++; }
            if(valid_Q21){ sum += Q21; count++; }
            if(valid_Q22){ sum += Q22; count++; }
            interpolated_values[i] = sum / count;
        }
        else if(valid_count == 2){
            // Linear interpolation fallback
            if(valid_Q11 && valid_Q21){
                interpolated_values[i] = (Q11*(x2 - x_interp) + Q21*(x_interp - x1)) / (x2 - x1);
            }
            else if(valid_Q12 && valid_Q22){
                interpolated_values[i] = (Q12*(x2 - x_interp) + Q22*(x_interp - x1)) / (x2 - x1);
            }
            else if(valid_Q11 && valid_Q12){
                interpolated_values[i] = (Q11*(y2 - y_interp) + Q12*(y_interp - y1)) / (y2 - y1);
            }
            else if(valid_Q21 && valid_Q22){
                interpolated_values[i] = (Q21*(y2 - y_interp) + Q22*(y_interp - y1)) / (y2 - y1);
            }
            else{
                cout << "Basin outline of basin with junction number " << basin_junction_number << ": Edge fallback failed at outline position " << i << endl;
                continue;
            }
        }
        else if(valid_count == 1){
            // Use the only valid value
            if(valid_Q11) interpolated_values[i] = Q11;
            else if(valid_Q12) interpolated_values[i] = Q12;
            else if(valid_Q21) interpolated_values[i] = Q21;
            else if(valid_Q22) interpolated_values[i] = Q22;
        }
        else{
            cout << "Basin outline of junction number " << basin_junction_number << " has missing data for pixel interpolation! Leaving outline position " << i << " as NoDataValue." << endl;
            continue;
        }

    }
    return interpolated_values;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Function to map values of a raster onto the basin outline node 
//
// RV 27/06/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDBasinOutline::map_pixel_values_onto_outline(LSDRaster& data_raster, 
    vector<float>& sampled_inbasin, vector<float>& sampled_outbasin,
    LSDFlowInfo& FI){

    sampled_inbasin.clear();
    sampled_inbasin.resize(NNodes_outline);
    sampled_outbasin.clear();
    sampled_outbasin.resize(NNodes_outline);
    // Loop through outline nodes
    for(int i=0; i<NNodes_outline; i++){
        // Get indices of pixels in and out of basin
        vector<int> px_in_basin, px_out_basin;
        get_pixels_around_basin_outline_node(i, px_in_basin, px_out_basin, FI);

        // First process the in basin
        // Get the first value
        int this_row, this_col;
        FI.retrieve_current_row_and_col(px_in_basin[0], this_row, this_col);
        float this_value = data_raster.get_data_element(this_row, this_col);
        float current_value = this_value;

        // If there are subsequent values only keep them if lower than currrent
        for(int j=1; j<px_in_basin.size(); j++){
            FI.retrieve_current_row_and_col(px_in_basin[j], this_row,this_col);
            this_value = data_raster.get_data_element(this_row, this_col);
            if(this_value>current_value){current_value=this_value;}
            // If the current value is a NoDataValue it might not store 
            // this value, since NoDataValue is typically highly negative.
            // So add extra logic to store thsivalue if current value is NDV
            else {if(current_value==NoDataValue){current_value=this_value;}}
        }

        // Add current value to mapped values
        sampled_inbasin[i] = current_value;

        // Now process out of basin
        // Get the first value
        FI.retrieve_current_row_and_col(px_out_basin[0], this_row, this_col);
        this_value = data_raster.get_data_element(this_row, this_col);
        current_value = this_value;

        // If there are subsequent values only keep them if lower than currrent
        for(int j=1; j<px_out_basin.size(); j++){
            FI.retrieve_current_row_and_col(px_out_basin[j], this_row,this_col);
            this_value = data_raster.get_data_element(this_row, this_col);
            if(this_value>current_value){current_value=this_value;}
            else {if(current_value==NoDataValue){current_value=this_value;}}
        }

        // Add current value to mapped values
        sampled_outbasin[i] = current_value;
    }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Function to get the 4 pixels surrounding a basin outline node, 
// distinguishing between pixels in the basin and pixels outside the basin.
//
// RV 27/06/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDBasinOutline::get_pixels_around_basin_outline_node(
    int outline_node_index, vector<int>& in_basin, vector<int>& out_basin,
    LSDFlowInfo& FI){
    
    int i = outline_node_index;
    in_basin.clear();
    out_basin.clear();

    // Get the nodes of the previous and next outline nodes too
    int next_i, prev_i;
    if(i==NNodes_outline-1){
        // Next node is at the start of outline
        next_i = 0; 
        // But if diagonal outlet flow, this is actually the same node,
        // so advance one more.
        if(outline_nodes[next_i] == outline_nodes[i]){
            next_i += 1;
        }
    } else {next_i = i+1;}
    if (i==0){
        // Previous node is at the back of outline
        prev_i = NNodes_outline-1;
        // But if diagonal outlet flow, this is actually the same node,
        // so go one further back.
        if(outline_nodes[i] == outline_nodes[prev_i]){
            prev_i -= 1;
        }
    } else {prev_i= i-1;}

    // Get the row and col of prev_i, i, next_i
    int prev_i_row = outline_rows[prev_i];
    int prev_i_col = outline_cols[prev_i];
    int i_row = outline_rows[i];
    int i_col = outline_cols[i];
    int next_i_row = outline_rows[next_i];
    int next_i_col = outline_cols[next_i];

    // Now check if the outline edges coming in and out of i have direction
    // north (0), east (1), south (2) or west (3).
    // Local lambdas for wrapped differences in case of periodic boundaries
    auto wrapped_row_diff = [&](int from, int to) -> int {
        if (!is_periodic_NS) return to - from;
        int diff = to - from;
        if (diff > NRows / 2) diff -= NRows;
        if (diff < -NRows / 2) diff += NRows;
        return diff;
    };

    auto wrapped_col_diff = [&](int from, int to) -> int {
        if (!is_periodic_EW) return to - from;
        int diff = to - from;
        if (diff > NCols / 2) diff -= NCols;
        if (diff < -NCols / 2) diff += NCols;
        return diff;
    };

    int edge_direction_in, edge_direction_out;

    int d_row_out = wrapped_row_diff(i_row, next_i_row);
    int d_col_out = wrapped_col_diff(i_col, next_i_col);

    if (d_row_out < 0) {
        edge_direction_out = 0; // North
    } else if (d_row_out > 0) {
        edge_direction_out = 2; // South
    } else if (d_col_out > 0) {
        edge_direction_out = 1; // East
    } else if (d_col_out < 0) {
        edge_direction_out = 3; // West
    } else {
        cout << "Warning: edge_direction_out not found!" << endl;
    }

    int d_row_in = wrapped_row_diff(prev_i_row, i_row);
    int d_col_in = wrapped_col_diff(prev_i_col, i_col);

    if (d_row_in < 0) {
        edge_direction_in = 0; // North
    } else if (d_row_in > 0) {
        edge_direction_in = 2; // South
    } else if (d_col_in > 0) {
        edge_direction_in = 1; // East
    } else if (d_col_in < 0) {
        edge_direction_in = 3; // West
    } else {
        cout << "Warning: edge_direction_in not found!" << endl;
    }


    // Get the x and y for this node in the reference of the flow grid
    float x_interp = i_col-0.5;
    float y_interp = i_row-0.5;

    // Get the x's and y's indices of the 4 pixels of the flow grid
    // these coordinates look as follows.
    //    x1  x2
    // y1  O  O
    // y2  O  O
    int raw_x1 = static_cast<int>(floor(x_interp));
    int raw_y1 = static_cast<int>(floor(y_interp));
    int raw_x2 = raw_x1 + 1;
    int raw_y2 = raw_y1 + 1;

    // Account for periodicity
    int x1 = wrap_index(raw_x1, NCols, is_periodic_EW);
    int x2 = wrap_index(raw_x2, NCols, is_periodic_EW);
    int y1 = wrap_index(raw_y1, NRows, is_periodic_NS);
    int y2 = wrap_index(raw_y2, NRows, is_periodic_NS);

    // Figure out which of the four nodes surrounding the edge are in the 
    // basin, using the fact that outlines are traced CLOCKWISE and that
    // we now the direction of the outline edges that come in and go out
    // of the outline node i.
    // Clockwise means that in basin pixel is on RHS when facing 
    // direction of edge.
    set<int> dem_pixels_in_basin, dem_pixels_out_basin;
    int pix_inb, pix_outb;

    // Lets start with outgoing node
    if(edge_direction_out == 0){
        // North: (y1,x2) in basin, (y1, x1) out of basin
        if(y1 != NoDataValue && x2 != NoDataValue){
            pix_inb = FI.retrieve_node_from_row_and_column(y1,x2);
        } else {
            pix_inb = NoDataValue;
        }
        if(y1 != NoDataValue && x1 != NoDataValue){
            pix_outb = FI.retrieve_node_from_row_and_column(y1,x1);
        } else {
            pix_outb = NoDataValue;
        }
    } else if (edge_direction_out == 1){
        // East: (y2,x2) in basin, (y1,x2) out of basin
        if(y2 != NoDataValue && x2 != NoDataValue){
            pix_inb = FI.retrieve_node_from_row_and_column(y2,x2);
        } else {
            pix_inb = NoDataValue;
        }
        if(y1 != NoDataValue && x2 != NoDataValue){
            pix_outb = FI.retrieve_node_from_row_and_column(y1,x2);
        } else {
            pix_outb = NoDataValue;
        }
    } else if (edge_direction_out == 2){
        // South: (y2,x1) in basin, (y2,x2) out of basin
        if(y2 != NoDataValue && x1 != NoDataValue){
            pix_inb = FI.retrieve_node_from_row_and_column(y2,x1);
        } else {
            pix_inb = NoDataValue;
        }
        if(y2 != NoDataValue && x2 != NoDataValue){
            pix_outb = FI.retrieve_node_from_row_and_column(y2,x2);
        } else {
            pix_outb = NoDataValue;
        }
    } else if (edge_direction_out == 3){
        // West: (y1,x1) in basin, (y2,x1) out of basin
        if(y1 != NoDataValue && x1 != NoDataValue){
            pix_inb = FI.retrieve_node_from_row_and_column(y1,x1);
        } else {
            pix_inb = NoDataValue;
        }
        if(y2 != NoDataValue && x1 != NoDataValue){
            pix_outb = FI.retrieve_node_from_row_and_column(y2,x1);
        } else {
            pix_outb = NoDataValue;
        }
    }
    if (pix_inb != NoDataValue){ dem_pixels_in_basin.insert(pix_inb); }
    if (pix_outb != NoDataValue){ dem_pixels_out_basin.insert(pix_outb); }

    // Now the ingoing node
    if(edge_direction_in == 0){
        // North: (y2,x2) in basin, (y2, x1) out of basin
        if(y2 != NoDataValue && x2 != NoDataValue){
            pix_inb = FI.retrieve_node_from_row_and_column(y2,x2);
        } else {
            pix_inb = NoDataValue;
        }
        if(y2 != NoDataValue && x1 != NoDataValue){
            pix_outb = FI.retrieve_node_from_row_and_column(y2,x1);
        } else {
            pix_outb = NoDataValue;
        }
    } else if (edge_direction_in == 1){
        // East: (y2,x1) in basin, (y1,x1) out of basin
        if(y2 != NoDataValue && x1 != NoDataValue){
            pix_inb = FI.retrieve_node_from_row_and_column(y2,x1);
        } else {
            pix_inb = NoDataValue;
        }
        if(y1 != NoDataValue && x1 != NoDataValue){
            pix_outb = FI.retrieve_node_from_row_and_column(y1,x1);
        } else {
            pix_outb = NoDataValue;
        }
    } else if (edge_direction_in == 2){
        // South: (y1,x1) in basin, (y1,x2) out of basin
        if(y1 != NoDataValue && x1 != NoDataValue){
            pix_inb = FI.retrieve_node_from_row_and_column(y1,x1);
        } else {
            pix_inb = NoDataValue;
        }
        if(y1 != NoDataValue && x2 != NoDataValue){
            pix_outb = FI.retrieve_node_from_row_and_column(y1,x2);
        } else {
            pix_outb = NoDataValue;
        }
    } else if (edge_direction_in == 3){
        // West: (y1,x2) in basin, (y2,x2) out of basin
        if(y1 != NoDataValue && x2 != NoDataValue){
            pix_inb = FI.retrieve_node_from_row_and_column(y1,x2);
        } else {
            pix_inb = NoDataValue;
        }
        if(y2 != NoDataValue && x2 != NoDataValue){
            pix_outb = FI.retrieve_node_from_row_and_column(y2,x2);
        } else {
            pix_outb = NoDataValue;
        }
    }
    if (pix_inb != NoDataValue){ dem_pixels_in_basin.insert(pix_inb); }
    if (pix_outb != NoDataValue){ dem_pixels_out_basin.insert(pix_outb); }

    // Corner logic with guards
    int corner_pix;
    if(edge_direction_in != edge_direction_out){
        auto safe_insert = [&](int y, int x, bool in_basin){
            if(y != NoDataValue && x != NoDataValue){
                corner_pix = FI.retrieve_node_from_row_and_column(y,x);
                if(corner_pix != NoDataValue){
                    if(in_basin){
                        dem_pixels_in_basin.insert(corner_pix);
                    } else {
                        dem_pixels_out_basin.insert(corner_pix);
                    }
                }
            }
        };
        if(edge_direction_in==0 && edge_direction_out==1){
            safe_insert(y1, x1, false);
        }
        if(edge_direction_in==0 && edge_direction_out==3){
            safe_insert(y1, x2, true);
        }
        if(edge_direction_in==1 && edge_direction_out==0){
            safe_insert(y2, x2, true);
        }
        if(edge_direction_in==1 && edge_direction_out==2){
            safe_insert(y1, x2, false);
        }
        if(edge_direction_in==2 && edge_direction_out==1){
            safe_insert(y2, x1, true);
        }
        if(edge_direction_in==2 && edge_direction_out==3){
            safe_insert(y2, x2, false);
        }
        if(edge_direction_in==3 && edge_direction_out==2){
            safe_insert(y1, x1, true);
        }
        if(edge_direction_in==3 && edge_direction_out==0){
            safe_insert(y2, x1, false);
        }
    }
  
    // If the outlet node was accidentally put in the out_basin set,
    // remove it and ensure it is in the in_basin set.
    // This error can arise from poorly defined incoming edge direction.
    // As outline node on the other side of the outlet may be diagonally
    // on the other side of the outlet pixel.
    int outlet_node = basin_outlet_node;   
    if (dem_pixels_out_basin.find(outlet_node) != dem_pixels_out_basin.end()) {
        dem_pixels_out_basin.erase(outlet_node);
        dem_pixels_in_basin.insert(outlet_node);
    }

    // Build the vector
    out_basin.assign(dem_pixels_out_basin.begin(), dem_pixels_out_basin.end());
    in_basin.assign(dem_pixels_in_basin.begin(), dem_pixels_in_basin.end());
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This gets all the linear indices of the nodes on the flowgrid that are 
// surrounding the outline on the outside of the basin.
//
// RV 27/06/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDBasinOutline::get_all_pixels_surrounding_outline_nodes_outside_basin(LSDFlowInfo& FI){
    // Initiate set to output (later converted to vector)
    set<int> nodes_to_return;
    
    // Walk through outline
    for(int i=0; i<NNodes_outline; i++){
        vector<int> px_in_basin, px_out_basin;
        get_pixels_around_basin_outline_node(i, px_in_basin, px_out_basin, FI);
        // insert the px_out_basin to the nodes to return
        for(size_t n=0; n<px_out_basin.size(); n++){
            nodes_to_return.insert(px_out_basin[n]);
        }
    }

    // Return
    vector<int> return_vec(nodes_to_return.begin(), nodes_to_return.end());
    return return_vec;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Function to get a vector of where the triple divides are located on the 
// basin outline.
//
// RV 27/06/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDBasinOutline::get_triple_divides_vector(LSDFlowInfo& FI, 
    LSDJunctionNetwork& JN){
    // Initiate vector
    vector<int> triple_divides_vector(NNodes_outline, 0);

    // Loop through outline nodes
    for(int i=0; i<NNodes_outline; i++){
        // Get indices of pixels in and out of basin
        vector<int> px_in_basin, px_out_basin;
        get_pixels_around_basin_outline_node(i, px_in_basin, px_out_basin, FI);
        vector<int> px_all;
        px_all.insert(px_all.end(), px_in_basin.begin(), px_in_basin.end());
        px_all.insert(px_all.end(), px_out_basin.begin(), px_out_basin.end());

        // We need to consider both the pixels in and out of the basin
        // and count the number of distinct basins these pixels belong to.
        // Note that, even though the basin outline outlines one basin,
        // if there are multiple px_in_basin, these may still belong to 
        // different basins since we are considering nested basins too. 
        set<int> all_basin_juncs;
        for(int k=0; k<px_all.size(); k++){
            // Get the node of the nearest channel (down-d8)
            int cnode=JN.get_nodeindex_of_nearest_channel(px_all[k], 1, FI);
            
            // Get the first upstream junction from that point on the channel
            // This junction uniquely defines basins (and therefore divides).
            int junc;
            if(cnode == NoDataValue){
                // If there is no channel, cnode will be NoDataValue.
                // Assign a junction number of -1 that stands for unchanneled
                // hillslope near DEM edge.
                // This ensures we can pick up triple-divides on say two 
                // parallel basin outlines where the space in between is filled
                // by a hillslope. 
                // Thus this introduces the assumption into the analysis that
                // DEM pixels not flowing into a channel are hillslopes.
                junc = -1;
            } else {
                junc=JN.find_upstream_junction_from_channel_nodeindex(cnode, FI);
            }
            // Insert in the set
            all_basin_juncs.insert(junc);
        }
        
        // If there are more than 2 different neighbouring basins, we have 
        // a triple divide junction ...
        if (all_basin_juncs.size()>2){
            // ... except if there is a diagonal flow node running through
            // this divide node (i.e. this means one basin is diagonally 
            // connected, two others lie to the side. The two divides 
            // between the diagonnally connected basin and each of the two
            // basins to the side touch each outher rather than junction.)

            // Check for each DEM pixel surrounding this node whether 
            // its downstream node is diagonally across through the divide node
            bool found_diagonal_flow = false;
            for(int k=0; k<px_all.size(); k++){
                int receiver_node;
                FI.retrieve_receiver_information(px_all[k], receiver_node);
                auto it = find(px_all.begin(), px_all.end(), receiver_node);
                if(it!=px_all.end()){
                    // Receiver is also neighbouring the divide node
                    int row, col;
                    FI.retrieve_current_row_and_col(px_all[k], row, col);
                    int fd = FI.get_LocalFlowDirection(row, col);
                    if (fd==7 || fd==1 || fd==5 || fd==3){
                        // Have diagonal flow!
                        found_diagonal_flow = true;
                    }
                }
            } 
            
            // If there's no diagonal flow, then we have a real triple divide
            if(!found_diagonal_flow){
                triple_divides_vector[i] = 1;
            }          
        }
    }
    return triple_divides_vector;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Finds the distance to nearest channel of each node in the basin outline.
// Returns a vector of NNodes_outline length.
//
// RV 30/06/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDBasinOutline::find_distance_of_outline_to_nearest_channel(
    LSDJunctionNetwork& JN, LSDFlowInfo& FI, LSDRaster& FlowDist,
    vector<float>& outline_to_dist_chan_inbasin_vector,
    vector<float>& outline_to_dist_chan_outbasin_vector){
    // Make empty vector which is to be returned
    outline_to_dist_chan_inbasin_vector.clear();
    outline_to_dist_chan_inbasin_vector.resize(NNodes_outline);
    outline_to_dist_chan_outbasin_vector.clear();
    outline_to_dist_chan_outbasin_vector.resize(NNodes_outline);

    // Loop over all outline nodes
    for(int i=0;i<NNodes_outline;i++){
        // Get the nodes that are inside and out of the basin
        vector<int> px_in_basin, px_out_basin;
        get_pixels_around_basin_outline_node(i, px_in_basin, px_out_basin, FI);

        // dem_pixels_in_basin and dem_pixels_out_basin should now (combined)
        // have four nodes! These are the surrounding nodes.
        // Lets get the distances to nearest channel in_basin and out_basin.
        
        // Get the first value
        int this_row, this_col;
        float this_dist = JN.get_distance_to_nearest_channel(px_in_basin[0], 1, 
            FI, FlowDist);
        float current_dist_inb = this_dist;

        // If there are subsequent values only keep them if lower than currrent
        for(int j=1; j<px_in_basin.size(); j++){
            this_dist = JN.get_distance_to_nearest_channel(px_in_basin[j], 1, 
                FI, FlowDist);
            if(this_dist<current_dist_inb){current_dist_inb=current_dist_inb;}
        }

        // Add current value to mapped values
        outline_to_dist_chan_inbasin_vector[i] = current_dist_inb;

        // Get the first value
        this_dist = JN.get_distance_to_nearest_channel(px_out_basin[0], 1, 
            FI, FlowDist);
        float current_dist_outb = this_dist;

        // If there are subsequent values only keep them if lower than currrent
        for(int j=1; j<px_out_basin.size(); j++){
            this_dist = JN.get_distance_to_nearest_channel(px_out_basin[j], 1, 
                FI, FlowDist);
            if(this_dist<current_dist_outb){current_dist_outb=current_dist_outb;}
        }

        // Add current value to mapped values
        outline_to_dist_chan_outbasin_vector[i] = current_dist_outb;
    }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Function to fill the basin outline elevation vector. IMPORTANT: should be 
// called before any function that requires outline elevations.
//
// RV 27/06/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDBasinOutline::fill_outline_elevation_vector(LSDRaster& elevation_raster){
    outline_elevation = interpolate_pixel_values_on_outline(elevation_raster);
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Function to write a specified part of the basin outline to a NEW 
// LSDRidgeline object. Returns the object.
// 
// Can also write in reverse (i.e. counter-clockwise), in this case the indices
// should still be givin in forward reference.
//
// Includes the start index, excludes the stop index.
//
// RV 06/10/2025
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDRidgeline LSDBasinOutline::write_outline_portion_to_LSDRidgeline(
    int start_index, int stop_index, bool write_in_reverse
){
    vector<int> ridgeline_nodes;
    if(!write_in_reverse){
        vector<int> temp(outline_nodes.begin()+start_index, 
                                        outline_nodes.begin()+stop_index);
        ridgeline_nodes = temp;
    } else {
        vector<int> temp(outline_nodes.begin()+stop_index+1,
                                 outline_nodes.begin()+start_index+1);
        reverse(temp.begin(), temp.end());
        ridgeline_nodes= temp;
    }
    

    LSDRidgeline ridgeline(
        NRows, NCols, NDataNodes, XMinimum, YMinimum, DataResolution, 
        GeoReferencingStrings, NoDataValue, boundary_conditions,
        ridgeline_nodes
    );

    return ridgeline;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Function returning two LSDRidgelines (returned in vector) by splitting
// the basin outline by max elevation.
//
// RV 26/08/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<LSDRidgeline> LSDBasinOutline::split_outline_into_ridgelines_by_max_elevation(){
    if(outline_elevation.empty()){
        cout << "Trying to split your basin outline into two ridgelines by max";
        cout << " elevation, but you need to fill_outline_elevation_vector";
        cout << " first!" << endl;
    }

    // Getting the maximum elevation index.
    auto max_elev_iter = max_element(outline_elevation.begin(), 
        outline_elevation.end());
    int current_index = distance(outline_elevation.begin(), max_elev_iter);

    vector<LSDRidgeline> ridgelines;
    ridgelines.push_back(write_outline_portion_to_LSDRidgeline(0, current_index,false));
    ridgelines.push_back(write_outline_portion_to_LSDRidgeline(NNodes_outline-1, current_index-1, true));
    return ridgelines;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Function returning two LSDRidgelines (returned in vector) by splitting
// the basin outline by minimum distance to head of longest channel.
//
// RV 26/08/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<LSDRidgeline> LSDBasinOutline::split_outline_into_ridgelines_by_longest_channel(
    LSDFlowInfo FI, LSDJunctionNetwork JN){
    // Determine what the channel head node is of the longest channel upstream
    // of the junction node that defines the basin. (May also be the junction
    // node itself, if the junction is a channel head.)
    LSDRaster distance_from_outlet = FI.distance_from_outlet();
    vector<int> longest_channel = JN.generate_longest_nodeindex_channel_in_basin(
        basin_junction_number, FI,distance_from_outlet);
    int ch_node = longest_channel[0];

    // Now we have to loop through all nodes, calculating the Euclidean
    // distance from the node to the channel head node.
    float current_euclidean_distance;
    int current_index;
    for(int i=0; i<NNodes_outline; i++){
        // Retrieve the row and col of the outline node and channel head
        int outline_row = outline_rows[i];
        int outline_col = outline_cols[i];
        int ch_row, ch_col;
        FI.retrieve_current_row_and_col(ch_node, ch_row, ch_col);

        // Now these rows and cols are not in the same reference frame
        // becasue divide grid is shifted half grid cell size with respect to 
        // flow grid. We can map both references onto a double-resolution grid.
        outline_row *= 2;
        outline_col *= 2;
        ch_row = ch_row*2+1; // +1 indicates half grid cell size shift
        ch_col = ch_col*2+1;

        float euclidean_distance = sqrt(pow((outline_row-ch_row), 2) 
            + pow((outline_col - ch_col), 2));

        if (i==0){
            current_euclidean_distance = euclidean_distance;
            current_index = i;
        } else if (euclidean_distance < current_euclidean_distance){
            current_euclidean_distance = euclidean_distance;
            current_index = i;
        }      
    } 

    vector<LSDRidgeline> ridgelines;
    ridgelines.push_back(write_outline_portion_to_LSDRidgeline(0, current_index, false));
    ridgelines.push_back(write_outline_portion_to_LSDRidgeline(NNodes_outline-1, current_index-1, true));
    return ridgelines;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Function returning two LSDRidgelines (returned in vector) by splitting
// the basin outline by equal length.
//
// RV 26/08/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<LSDRidgeline> LSDBasinOutline::split_outline_into_ridgelines_by_equal_length(){
    // Get the perimeter of the basin outline, and half the length
    float perimeter = outline_distance_cw.back();
    float current_distance = 0;
    int current_index = 0;
    while (current_distance < perimeter/2){
        current_index++;
        current_distance = outline_distance_cw[current_index];
    } 
    vector<LSDRidgeline> ridgelines;
    ridgelines.push_back(write_outline_portion_to_LSDRidgeline(0, current_index, false));
    ridgelines.push_back(write_outline_portion_to_LSDRidgeline(NNodes_outline-1, current_index-1, false));
    return ridgelines;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// If the basin is range-spanning and has a transverse divide on its left
// side (i.e., seen clockwise along outline before hitting the main drainage
// divide), it returns an LSDRidgeline of this transverse divide.
//
// If it does not have a transverse divide, or if this outline is not 
// range-spanning, it returns an empty LSDRidgeline.
//
// RV 06/10/2025
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDRidgeline LSDBasinOutline::get_left_transverse_divide_as_ridgeline(
    vector<int> range_spanning_basins_N, vector<int> range_spanning_basins_S,
    LSDFlowInfo FI, LSDJunctionNetwork JN
) {
    // Retrieve the start end nodes of the transverse divide
    // find_left_transverse_divide_nodes
    vector<int> start_end_nodes = find_left_transverse_divide_nodes(
        range_spanning_basins_N,
        range_spanning_basins_S,
        FI, JN
    );

    // If no transverse divide found then return empty ridgeline
    // Note that find_left_transverse_divide_nodes automatically returns
    // an empty vector if this basin is not range spanning to start with.
    if (start_end_nodes.size()!=2){
        LSDRidgeline RL;
        return RL;
    }
      
    // Otherwise use write_outline_portion_to_LSDRidgeline to make an
    // LSDRidgeline of the right transverse divide
    // Note the reverse ridgeline making so to start at the TDJ.
    LSDRidgeline RL = write_outline_portion_to_LSDRidgeline(
        start_end_nodes[0], 
        start_end_nodes[1],
        false);
    return RL;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// If the basin is range-spanning and has a transverse divide on its right
// side (i.e., seen clockwise along outline before hitting the main drainage
// divide), it returns an LSDRidgeline of this transverse divide.
//
// If it does not have a transverse divide, or if this outline is not 
// range-spanning, it returns an empty LSDRidgeline.
//
// RV 26/08/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDRidgeline LSDBasinOutline::get_right_transverse_divide_as_ridgeline(
    vector<int> range_spanning_basins_N, vector<int> range_spanning_basins_S,
    LSDFlowInfo FI, LSDJunctionNetwork JN
) {
    // Retrieve the start end nodes of the transverse divide
    // find_right_transverse_divide_nodes
    vector<int> start_end_nodes = find_right_transverse_divide_nodes(
        range_spanning_basins_N,
        range_spanning_basins_S,
        FI, JN
    );

    // If no transverse divide found then return empty ridgeline
    // Note that find_right_transverse_divide_nodes automatically returns
    // an empty vector if this basin is not range spanning to start with.
    if (start_end_nodes.size()!=2){
        LSDRidgeline RL;
        return RL;
    }
      
    // Otherwise use write_outline_portion_to_LSDRidgeline to make an
    // LSDRidgeline of the right transverse divide
    // Note the reverse ridgeline making so to start at the TDJ.
    LSDRidgeline RL = write_outline_portion_to_LSDRidgeline(
        start_end_nodes[0], 
        start_end_nodes[1],
        true);
    return RL;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This gets the portion of the basin outline that is on the main drainage 
// divide, inclusive of both triple divide junctions where the basin enters and
// leaves the main drainage divide.
//
// The MDD is always exported in the direction of the basin outline, that is 
// clockwise.
//
// If the basin does not touch the main drainage divide, an empty LSDRidgeline
// is returned.
//
// RV 26/08/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDRidgeline LSDBasinOutline::get_outline_nodes_on_main_drainage_divide_as_ridgeline(
    vector<int> range_spanning_basins_N, vector<int> range_spanning_basins_S,
    LSDFlowInfo FI, LSDJunctionNetwork JN
){
    // Get the outline indices of the nodes where the MDD starts/ends
    vector<int> start_end_nodes = find_main_drainage_divide_nodes(
        range_spanning_basins_N,
        range_spanning_basins_S,
        FI, JN
    );

    // If not on MDD, start_end_nodes will be empty. Return empty LSDRidgeline.
    if (start_end_nodes.size()!=2){
        LSDRidgeline RL;
        return RL;
    }

    // Write LSDRidgeline that contains the MDD portion of this outline
    LSDRidgeline RL = write_outline_portion_to_LSDRidgeline(
        start_end_nodes[0], 
        start_end_nodes[1],
        false);
    return RL;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This gets the portion of the basin outline that is shared with another basin
// inclusive of both triple divide junctions where the basin joins and
// leaves the portion that is shared with another basin.
//
// The interfluvial divide is always exported in the direction of the basin 
// outline, that is clockwise.
//
// If the basin does not touch the other basin, an empty LSDRidgeline
// is returned.
//
// RV 30/01/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDRidgeline LSDBasinOutline::get_interfluvial_divide_as_ridgeline(
    int otherBasin_junction_number, LSDFlowInfo FI, LSDJunctionNetwork JN
){
    // Get the outline indices of the nodes where the interfluvial divide starts/ends
    vector<int> start_end_nodes = find_interfluvial_divide_nodes(
        otherBasin_junction_number,
        FI, JN
    );

    // If no interfluvial divide between the two basins, start_end_nodes will be empty. Return empty LSDRidgeline.
    if (start_end_nodes.size()!=2){
        LSDRidgeline RL;
        return RL;
    }

    // Write LSDRidgeline that contains the MDD portion of this outline
    LSDRidgeline RL = write_outline_portion_to_LSDRidgeline(
        start_end_nodes[0], 
        start_end_nodes[1],
        false);
    return RL;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This scans the outline for triple divide junctions where the outline 
// (potentially) joins and leaves the left transverse divide.
// By iterating through the outline clockwise (forward).
//
// Returns a vector of outline indices.
// A transverse divide is always at least 1 pixel edge long, so will have a
// distinct joining and leaving node.
// If the basin is not range-spanning or if it has no transverse divide on this
// side, an emtpy vector is returned.
//
// RV 26/08/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDBasinOutline::find_left_transverse_divide_nodes(
    vector<int> range_spanning_basins_N, vector<int> range_spanning_basins_S,
    LSDFlowInfo FI, LSDJunctionNetwork JN
) {
    // Find out which side the current basin is in by comparing basin_junction_
    // number to range_spanning_basins_N and range_spanning_basins_S
    bool basin_is_north = false; 
    bool basin_is_south = false;
    if(find(range_spanning_basins_N.begin(), range_spanning_basins_N.end(), 
        basin_junction_number) != range_spanning_basins_N.end()){
        basin_is_north = true;
    } else if (find(range_spanning_basins_S.begin(), range_spanning_basins_S.end(), 
        basin_junction_number) != range_spanning_basins_S.end()){
        basin_is_south = true;
    } else {
        return {};
    }

    vector<int> rspb_thisside, rspb_otherside;
    if(basin_is_north){
        rspb_thisside=range_spanning_basins_N;
        rspb_otherside=range_spanning_basins_S;
    }
    if(basin_is_south){
        rspb_thisside=range_spanning_basins_S;
        rspb_otherside=range_spanning_basins_N;
    }
    
    vector<int> nodes_to_return;
    // Iterate over the nodes in the basin outline, forwards
    for(int i=0; i<NNodes_outline; i++){
        // Get the DEM pixels outside of the basin surrounding outline node
        vector<int> px_in_basin, px_out_basin;
        get_pixels_around_basin_outline_node(i,px_in_basin, px_out_basin, FI);

        /*
        cout << "Iterating, basin outline node " << i << endl;
        cout << "There are " << px_out_basin.size() << " pixels outside the basin: " << endl;
        for(int node : px_out_basin){
            int row, col;
            FI.retrieve_current_row_and_col(node, row, col);
            cout << node << ": " << row << "/" << col << endl;
            }
        cout << "There are " << px_in_basin.size() << " pixels insdie the basin " << endl;
        for(int node : px_in_basin){
            int row, col;
            FI.retrieve_current_row_and_col(node, row, col);
            cout << node << ": " << row << "/" << col << endl;
            }
        */
        
        // Loop through the pixels outside of the basin
        int N_pixels_in_adjacent_rsb=0;
        vector<int> outside_nodes_not_in_adjacent_rsb;
        for(int j=0; j<px_out_basin.size(); j++){
            // Get the baselevel junction of this neighbouring pixel
            int thisNode = px_out_basin[j];
            int thisNodesBL = FI.get_baselevel_of_node(thisNode);
            int thisBLJunc = JN.get_Junction_of_Node(thisNodesBL, FI);
            //cout << "This node is " << thisNode << endl;
            //cout << "This nodes baselevel is " << thisNodesBL << endl;
            //cout << "This nodes baselevel junction is " << thisBLJunc << endl;
          
            // Check if out basin pixel is in range-spanning basin on this side
            if (find(rspb_thisside.begin(), rspb_thisside.end(), thisBLJunc)
                != rspb_thisside.end()){
                N_pixels_in_adjacent_rsb++;
            } else {
                outside_nodes_not_in_adjacent_rsb.push_back(thisNode);
            }
        }
        //cout << "I got " << N_pixels_in_adjacent_rsb << " pixels in adjacent rsb." << endl;

        // If the current nodes to return is empty and
        // if the node has a surrounding pixel outside the 
        // current basin that is in another range spanning basin on the 
        // same side, we have found the start of the TVD
        if(nodes_to_return.empty() && N_pixels_in_adjacent_rsb>0){
            bool join = true;
            // Ok, need to deal with an edge case:
            // This outline node has a surrounding pixel in a neighbouring range 
            // spanning basin, but it is separated by a diagonal flow edge that 
            // is in neither range-spanning basin. 
            // We can identify such cases as:
            // (1) 3 pixels outside the basin
            // (2) one of which in a different range-spanning basin and
            // (3) the 2 remaining pixels are connected by a 
            //     diagonal flow edge.
            if(px_out_basin.size() == 3 &&
               N_pixels_in_adjacent_rsb == 1 &&
               outside_nodes_not_in_adjacent_rsb.size() == 2
            ){
                int node1 = outside_nodes_not_in_adjacent_rsb[0];
                int node1r, node1c, node1receiver;
                FI.retrieve_current_row_and_col(node1, node1r, node1c);
                FI.retrieve_receiver_information(node1, node1receiver);
                
                int node2 = outside_nodes_not_in_adjacent_rsb[1];
                int node2r, node2c, node2receiver;
                FI.retrieve_current_row_and_col(node2, node2r, node2c);
                FI.retrieve_receiver_information(node2, node2receiver);

                if(node1receiver == node2 || node2receiver == node1){
                    // Ok, the two nodes outside basin draining same edge
                    // are each others receiver.
                    // Now check if they are diagonally connected.
                    if(node1r != node2r && node1c != node2c){
                        // Diagonally connected. Set join to false.
                        join = false;
                    }
                }
            } 
            
            if(join){
                nodes_to_return.push_back(i);
                //cout << "Recording joining node." << endl;
            }            
        }
            
        // Stop if the node has a surrounding pixel outside the current
        // basin that is NOT touching a range-spanning basin
        // In this case, we have just left the Transverse divide
        if(nodes_to_return.size() == 1 && N_pixels_in_adjacent_rsb==0){
            nodes_to_return.push_back(i);
            //cout << "Recording leaving node." << endl;
            return nodes_to_return;
            
        }
    }

    // If not returned previously, a transverse divide was not found
    // on this side.
    // This can happen due to it being the first basin along a non-periodic
    // boundary, or due to interstitial unchanneled space in-between range-spanning
    // basins (that stretches all the way from mdd to boundary).
    //cout << "Basin " << basin_junction_number << " I found no left tvd!" << endl;
    return nodes_to_return;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This scans the outline for triple divide junctions where the outline 
// (potentially) joins and leaves the right transverse divide.
// By iterating through the outline counter-clockwise.
//
// Returns a vector of outline indices.
// A transverse divide is always at least 1 pixel edge long, so will have a
// distinct joining and leaving node.
// If the basin is not range-spanning or if it has no transverse divide on this
// side, an emtpy vector is returned.
//
// RV 26/08/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDBasinOutline::find_right_transverse_divide_nodes(
    vector<int> range_spanning_basins_N, vector<int> range_spanning_basins_S,
    LSDFlowInfo FI, LSDJunctionNetwork JN
) {
    // Find out which side the current basin is in by comparing basin_junction_
    // number to range_spanning_basins_N and range_spanning_basins_S
    bool basin_is_north = false; 
    bool basin_is_south = false;
    if(find(range_spanning_basins_N.begin(), range_spanning_basins_N.end(), 
        basin_junction_number) != range_spanning_basins_N.end()){
        basin_is_north = true;
    } else if (find(range_spanning_basins_S.begin(), range_spanning_basins_S.end(), 
        basin_junction_number) != range_spanning_basins_S.end()){
        basin_is_south = true;
    } else {
        return {};
    }

    vector<int> rspb_thisside, rspb_otherside;
    if(basin_is_north){
        rspb_thisside=range_spanning_basins_N;
        rspb_otherside=range_spanning_basins_S;
    }
    if(basin_is_south){
        rspb_thisside=range_spanning_basins_S;
        rspb_otherside=range_spanning_basins_N;
    }
    
    vector<int> nodes_to_return;
    // Iterate over the nodes in the basin outline, backwards
    for(int i=NNodes_outline-1; i>=0; i--){
        // Get the DEM pixels outside of the basin surrounding outline node
        vector<int> px_in_basin, px_out_basin;
        get_pixels_around_basin_outline_node(i,px_in_basin, px_out_basin, FI);

        // Loop through the pixels outside of the basin
        int N_pixels_in_adjacent_rsb=0;
        vector<int> outside_nodes_not_in_adjacent_rsb;
        for(int j=0; j<px_out_basin.size(); j++){
            // Get the baselevel junction of this neighbouring pixel
            int thisNode = px_out_basin[j];
            int thisNodesBL = FI.get_baselevel_of_node(thisNode);
            int thisBLJunc = JN.get_Junction_of_Node(thisNodesBL, FI);

            // Check if out basin pixel is in range-spanning basin on this side
            if (find(rspb_thisside.begin(), rspb_thisside.end(), thisBLJunc)
                != rspb_thisside.end()){
                N_pixels_in_adjacent_rsb++;
            } else {
                outside_nodes_not_in_adjacent_rsb.push_back(thisNode);
            }
        }

        // If the current nodes to return is empty and
        // if the node has a surrounding pixel outside the 
        // current basin that is in another range spanning basin on the 
        // same side, we have found the start of the TDJ
        if(nodes_to_return.empty() && N_pixels_in_adjacent_rsb>0){
            bool join = true;
            // Ok, need to deal with an edge case:
            // This outline node has a surrounding pixel in a neighbouring range 
            // spanning basin, but it is separated by a diagonal flow edge that 
            // is in neither range-spanning basin. 
            // We can identify such cases as:
            // (1) 3 pixels outside the basin
            // (2) one of which in a different range-spanning basin and
            // (3) the 2 remaining pixels are connected by a 
            //     diagonal flow edge.
            if(px_out_basin.size() == 3 &&
               N_pixels_in_adjacent_rsb == 1 &&
               outside_nodes_not_in_adjacent_rsb.size() == 2
            ){
                int node1 = outside_nodes_not_in_adjacent_rsb[0];
                int node1r, node1c, node1receiver;
                FI.retrieve_current_row_and_col(node1, node1r, node1c);
                FI.retrieve_receiver_information(node1, node1receiver);
                
                int node2 = outside_nodes_not_in_adjacent_rsb[1];
                int node2r, node2c, node2receiver;
                FI.retrieve_current_row_and_col(node2, node2r, node2c);
                FI.retrieve_receiver_information(node2, node2receiver);

                if(node1receiver == node2 || node2receiver == node1){
                    // Ok, the two nodes outside basin draining same edge
                    // are each others receiver.
                    // Now check if they are diagonally connected.
                    if(node1r != node2r && node1c != node2c){
                        // Diagonally connected. Set join to false.
                        join = false;
                    }
                }
            } 
            
            if(join){
                nodes_to_return.push_back(i);
                //cout << "Recording joining node." << endl;
            }            
        }
            
        // Stop if the node has a surrounding pixel outside the current
        // basin that is NOT touching a range-spanning basin
        // In this case, we have just left the Transverse divide
        if(nodes_to_return.size() == 1 && N_pixels_in_adjacent_rsb==0){
            nodes_to_return.push_back(i);
            return nodes_to_return;
        }
    }

    // If not returned previously, a transverse divide was not found
    // on this side.
    return nodes_to_return;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This scans the outline for the points where the outline 
// (potentially) joins and leaves the main drainage divide.
//
// The start/end nodes are returned in a vector. If the returned vector if empty,
// the basin outline does not touch or stretch along the MDD.
//
// RV 26/08/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDBasinOutline::find_main_drainage_divide_nodes(
    vector<int> range_spanning_basins_N, vector<int> range_spanning_basins_S,
    LSDFlowInfo FI, LSDJunctionNetwork JN
) {
    // Find out which side the current basin is in by taking the current basin
    // junction number and tracing it down to the baselevel junction.
    // Then compare this baselevel junction to range_spanning_basins_N and 
    // range_spanning_basins_S.
    int basin_bl_junction;
    basin_bl_junction = JN.find_base_level_node_of_junction(basin_junction_number);

    bool basin_is_north = false; 
    bool basin_is_south = false;
    if(find(range_spanning_basins_N.begin(), range_spanning_basins_N.end(), 
        basin_bl_junction) != range_spanning_basins_N.end()){
        basin_is_north = true;
    } else if (find(range_spanning_basins_S.begin(), range_spanning_basins_S.end(), 
        basin_bl_junction) != range_spanning_basins_S.end()){
        basin_is_south = true;
    } else {
        cout << "I am trying to find the triple divide junctions where this ";
        cout << "basin outline joins and leaves the main drainage divide, but";
        cout << "it is looking like this basin does not actually touch on the"; 
        cout << "main drainage divide. Returning an empty vector."<<endl;
        vector<int> empty_vec;
        return empty_vec;
    }

    vector<int> rspb_thisside, rspb_otherside;
    if(basin_is_north){
        rspb_thisside=range_spanning_basins_N;
        rspb_otherside=range_spanning_basins_S;
    }
    if(basin_is_south){
        rspb_thisside=range_spanning_basins_S;
        rspb_otherside=range_spanning_basins_N;
    }

    // Make a vector for joining indices and leaving indices
    vector<int> join_indices;
    vector<int> leave_indices;
    bool on_mdd = false;

    // Iterate over the nodes in the basin outline
    for(int i=0; i<NNodes_outline; i++){
        // Get the DEM pixels outside of the basin surrounding outline node
        vector<int> px_in_basin, px_out_basin;
        get_pixels_around_basin_outline_node(i,px_in_basin, px_out_basin, FI);

        // Loop through the pixels outside of the basin
        // And get the number of different basins on the same side
        int N_pixels_draining_opposite_side=0;
        vector<int> nodes_outside_basin_draining_same_side;
        for(int j=0; j<px_out_basin.size(); j++){
            // Get the baselevel junction of this neighbouring pixel
            int thisNode = px_out_basin[j];
            int thisNodesBL = FI.get_baselevel_of_node(thisNode);
            int thisBLJunc = JN.get_Junction_of_Node(thisNodesBL, FI);

            // See if any flow to opposite side
            if (find(rspb_otherside.begin(), rspb_otherside.end(), thisBLJunc)
               != rspb_otherside.end()){
               N_pixels_draining_opposite_side++;
            }

            // Identify those drainaing to the same side
            if (find(rspb_thisside.begin(), rspb_thisside.end(), thisBLJunc)
               != rspb_thisside.end()){
               nodes_outside_basin_draining_same_side.push_back(thisNode);
            }
        }

        int N_pixels_outside_basin = px_out_basin.size();

        // If we are on the mdd but we find that the 
        // current node does NOT have any surounding
        // pixels in a range-spanning basin on the opposite side, we know
        // that we have left the MDD and found the leaving node (the 
        // "leaving node" is the first node not on the MDD anymore).
        if(on_mdd && N_pixels_draining_opposite_side==0){
            int leaving_index = i;
            // Here we need to treat an edge case. If the node before
            // the leaving node is on the same location as the joining node
            // then we need to mark the leaving node one earlier.
            // This becasue it has essentially "lapped back" to the node
            // where it joins, where a diagonal flow edge connects the basin
            // pixel. But when lapping back (i.e. visiting the second time)
            // you actually have that flow edge separating the divide and 
            // the basin on the other side.
            if(outline_nodes[leaving_index-1]==outline_nodes[join_indices.back()]
              && leaving_index-1 != join_indices.back()){
                leaving_index--;
            }

            // Push back the leaving index
            leave_indices.push_back(leaving_index);
            on_mdd = false;
        }

        // If we are not on the mdd and
        // the pixels surrounding the node outside the 
        // current basin is in a range-spanning basin on the opposite side.
        // Then, we have found the first TDJ (we are joining the MDD)
        // Push back this index to the return vector
        if(!on_mdd && N_pixels_draining_opposite_side>0){
            int joining_index = i;
            bool join = true;
            
            // Ok, need to deal with an edge case:
            // The case where the node has a diagonal flow edge from another
            // basin on the same side separating the node from the basin on the
            // opposite side. We can detect these if there are 
            // (1) 3 pixels outside the basin
            // (2) one of which is on the other side of the MDD and
            // (3) the 2 remaining pixels on the same side are connected by a 
            //     diagonal flow edge.
            if(N_pixels_outside_basin == 3 &&
               N_pixels_draining_opposite_side == 1 &&
               nodes_outside_basin_draining_same_side.size()==2
            ){
                int node1 = nodes_outside_basin_draining_same_side[0];
                int node1r, node1c, node1receiver;
                FI.retrieve_current_row_and_col(node1, node1r, node1c);
                FI.retrieve_receiver_information(node1, node1receiver);
                
                int node2 = nodes_outside_basin_draining_same_side[1];
                int node2r, node2c, node2receiver;
                FI.retrieve_current_row_and_col(node2, node2r, node2c);
                FI.retrieve_receiver_information(node2, node2receiver);

                if(node1receiver == node2 || node2receiver == node1){
                    // Ok, the two nodes outside basin draining same edge
                    // are each others receiver.
                    // Now check if they are diagonally connected.
                    if(node1r != node2r && node1c != node2c){
                        // Diagonally connected. Set join to false.
                        join = false;
                    }
                }
            }
            
            if(join){
                join_indices.push_back(joining_index);
                on_mdd = true;
            }
        }
    }    

    if(join_indices.empty()){
        cout << "LSDBasinOutline::find_main_drainage_divide_nodes: " ;
        cout << "Check if this is correct! A range-spanning basin that never ";
        cout << "joins the mdd..." << endl;
        return {};
    }

    if(leave_indices.empty()){
        // add the final index + 1 as the final node
        cout << "LSDBasinOutline::find_main_drainage_divide_nodes: " ;
        cout << "You have reached a previously unencountered edge case where ";
        cout << "the outline only joins but never leaves the MDD. Be careful ";
        cout << "and check if your results look OK!" << endl;
        leave_indices.push_back(NNodes_outline);
    }

    return {join_indices.front(), leave_indices.back()};
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This scans the outline for the points where the outline 
// (potentially) joins and leaves the outline of another, specified basin.
//
// The start/end nodes are returned in a vector. If the returned vector if empty,
// the basin outline does not touch or stretch along the outline of another basin.
//
// RV 26/08/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDBasinOutline::find_interfluvial_divide_nodes(
    int otherBasin_junction_number, LSDFlowInfo FI, LSDJunctionNetwork JN
) {
    
    // Make a vector for joining indices and leaving indices
    vector<int> join_indices;
    vector<int> leave_indices;
    bool on_interfluvial_divide = false;

    //cout << "chekcing for nesting!" << endl;
    // Check if the otherBasin_junction_number is NOT nested in the current basin junction number, or
    // vice versa
    if (
        JN.is_junction_upstream(basin_junction_number, otherBasin_junction_number) ||
        JN.is_junction_upstream(otherBasin_junction_number, basin_junction_number)
     
    ){
        //cout << "nesting identified" << endl;
        return {};
    }

    //cout << "Iterating over outline nodes" << endl;
    // Iterate over the nodes in the basin outline
    for(int i=0; i<NNodes_outline; i++){
        //cout << "i="<<i<<"/"<<NNodes_outline<<endl;
        // Get the DEM pixels surrounding the outline node on the outside of the basin
        vector<int> px_in_basin, px_out_basin;
        get_pixels_around_basin_outline_node(i,px_in_basin, px_out_basin, FI);  

        // Now check if any of these pixels are in the other basin
        // and record the nodes that are in neither the current and the other basin.
        int N_pixels_draining_other_basin=0;
        vector<int> nodes_in_neither_basin;

        for(int j=0; j<px_out_basin.size(); j++){
            // Get the first channelised node down the d8 flow path from thisNode
            int thisNode = px_out_basin[j];
            int thisNodesNearestChannelNode = JN.get_nodeindex_of_nearest_channel(thisNode, 1, FI);

            // Check for the case that there is no channelised node
            if(thisNodesNearestChannelNode == NoDataValue){
                nodes_in_neither_basin.push_back(thisNode);
                continue;
            }

            // Now get the first upstream junction of that channelised node
            int thisNodesJunction = JN.find_upstream_junction_from_channel_nodeindex(thisNodesNearestChannelNode, FI);

            // And check if this junction is upstream of or equal to the junction number of the other basin            
            if(JN.is_junction_upstream(otherBasin_junction_number, thisNodesJunction)){
                // Yes, this node is in the other basin
                N_pixels_draining_other_basin++;
            }       
            else {
                // No, the node is in neither the current nor the other basin
                nodes_in_neither_basin.push_back(thisNode);
            }     
        }
        int N_pixels_outside_basin = px_out_basin.size();

        //cout << "N_px_draining_other_basin:" << N_pixels_draining_other_basin << endl;
        //cout << "N_pixels_outside_basin:" << N_pixels_outside_basin << endl;
        //cout << "N_pixels_in_neither_basin:"<< nodes_in_neither_basin.size() << endl;

        // If we are on the interfluvial divide between this and the other basin
        // and the current node deos not have any outside pixels in the other 
        // basin, then we know that we have left the interfluvial divide and 
        // record the leaving node. The leaving node is the first node not on the
        // interfluvial divide anymore.
        if(on_interfluvial_divide && N_pixels_draining_other_basin==0){
            int leaving_index = i;
            // Here we need to treat an edge case. If the node before
            // the leaving node is on the same location as the joining node
            // then we need to mark the leaving node one earlier.
            // This becasue it has essentially "lapped back" to the node
            // where it joins, where a diagonal flow edge connects the basin
            // pixel. But when lapping back (i.e. visiting the second time)
            // you actually have that flow edge separating the divide and 
            // the basin on the other side.
            if(outline_nodes[leaving_index-1]==outline_nodes[join_indices.back()]
              && leaving_index-1 != join_indices.back()){
                leaving_index--;
            }

            // Push back the leaving index
            leave_indices.push_back(leaving_index);
            //cout << "Leaving the interfluvial divide" << endl;
            on_interfluvial_divide = false;
        }

        // If we are not on the interfluvial divide yet and one or more pixels
        // surrunding the current divide node on the outside is in the other basin, 
        // then we record this as the leaving node by pushing back this index to 
        // the return vector.
        if(!on_interfluvial_divide && N_pixels_draining_other_basin>0){
            int joining_index = i;
            bool join = true;
            
            // Ok, need to deal with an edge case:
            // The case where the node has a diagonal flow edge from another
            // basin separating the outline node from the other basin. 
            // We can detect these if there are 
            // (1) 3 pixels outside the basin
            // (2) one of which is in the other basin
            // (3) the 2 remaining pixels on the same side are connected by a 
            //     diagonal flow edge.
            if(N_pixels_outside_basin == 3 &&
               N_pixels_draining_other_basin == 1 &&
               nodes_in_neither_basin.size()==2
            ){
                int node1 = nodes_in_neither_basin[0];
                int node1r, node1c, node1receiver;
                FI.retrieve_current_row_and_col(node1, node1r, node1c);
                FI.retrieve_receiver_information(node1, node1receiver);
                
                int node2 = nodes_in_neither_basin[1];
                int node2r, node2c, node2receiver;
                FI.retrieve_current_row_and_col(node2, node2r, node2c);
                FI.retrieve_receiver_information(node2, node2receiver);

                if(node1receiver == node2 || node2receiver == node1){
                    // Ok, the two nodes outside basin draining same edge
                    // are each others receiver.
                    // Now check if they are diagonally connected.
                    if(node1r != node2r && node1c != node2c){
                        // Diagonally connected. Set join to false.
                        join = false;
                    }
                }
            }
            
            if(join){
                //cout << "joining the interfluvial divide" << endl;
                join_indices.push_back(joining_index);
                on_interfluvial_divide = true;
            }
        }
    }    

    if(join_indices.size() == 1 && leave_indices.empty()){
        // add the final index + 1 as the final node
        cout << "LSDBasinOutline::find_interfluvial_divide_nodes: " ;
        cout << "You have reached a previously unencountered edge case where ";
        cout << "the outline only joins but never leaves the interfluvial divide! Be careful ";
        cout << "and check if your results look OK!" << endl;
        leave_indices.push_back(NNodes_outline);
    }

    if(join_indices.empty() && leave_indices.empty()){
        return {};
    }

    return {join_indices.front(), leave_indices.back()};
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Prints basin outlines to csv with elevations 
//
// RV 26/08/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDBasinOutline::print_basin_outline_to_csv(string csv_name,
    LSDJunctionNetwork JN, vector<vector<float>> additional_mapped_data ){
    
    // Get the coordinate converter
    LSDCoordinateConverterLLandUTM Converter;

    // Open the csv
    ofstream csv_out;
    csv_out.open(csv_name.c_str(), ios::app); // ios:app prevents overwriting
                                   // the csv, instead appends to existing csv.
    csv_out.precision(9);

    // Run through all nodes
    for(int i=0; i<NNodes_outline; i++){
        // Junction number 
        csv_out << basin_junction_number << ",";

        // Latitude and longitude
        double lat, lon;
        get_lat_and_long_locations_dividegrid(outline_rows[i],
            outline_cols[i], lat, lon, Converter, JN);
        csv_out << lat << "," << lon <<  ",";

        // The clockwise distance
        csv_out << outline_distance_cw[i] << ",";

        // The counterclockwise distance
        csv_out << outline_distance_ccw[i] << ",";

        // The elevations
        csv_out << outline_elevation[i];

        // The additional data if supplied
        if (!additional_mapped_data.empty()){
            for(int j=0; j<additional_mapped_data.size(); j++){
                csv_out << "," << additional_mapped_data[j][i];
            }
        }

        // Finish line
        csv_out << endl;
    }
    csv_out.close();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Prints basin outlines to csv with elevations. This version prints divide
// row and cols instead of latitude/longitude.
//
// RV 26/08/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDBasinOutline::print_basin_outline_to_csv_rowscols(string csv_name,
    vector<vector<float>> additional_mapped_data ){

    // Open the csv
    ofstream csv_out;
    csv_out.open(csv_name.c_str(), ios::app); // ios:app prevents overwriting
                                   // the csv, instead appends to existing csv.
    csv_out.precision(9);

    // Run through all nodes
    for(int i=0; i<NNodes_outline; i++){
        // Junction number 
        csv_out << basin_junction_number << ",";

        // Divide row and col
        csv_out << outline_rows[i] << "," << outline_cols[i] <<  ",";

        // The clockwise distance
        csv_out << outline_distance_cw[i] << ",";

        // The counterclockwise distance
        csv_out << outline_distance_ccw[i] << ",";

        // The elevations
        csv_out << outline_elevation[i] << ",";

        // The additional data if supplied
        if (!additional_mapped_data.empty()){
            for(int j=0; j<additional_mapped_data.size(); j++){
                csv_out << "," << additional_mapped_data[j][i];
            }
        }

        // Finish line
        csv_out << endl;
    }
    csv_out.close();
}

