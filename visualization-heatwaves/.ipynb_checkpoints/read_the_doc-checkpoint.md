# heatwaves.py (v0.1)

A quick and dirty reference manual

(In process...........)

### List of functions

#### Heatwave calculation

This step operates on the E-OBS maximum temperature dataset at 0.5deg of resolution. For 
each date within a heatwave period, it will calculate a `context`. A `context` is a period
of time surrounding each day belonging to a heatwave period. In this case, we are using
as context the two dates around a given date as in the example:

**Context for 1-8-2003**

`context = [30-7-2003, 31-7-2003, 1-8-2003, 2-8-2003, 3-8-2003]`

This means that a date is not checked individually, but within a small group. Hence, the
bigger the `context` the more robust the verification of a heatwave day. This is an 
important concept in this step, and the functions below are structured around it. 


- `gendates`: For a given date, returns all dates in the past X years (i.e. climatology)


- `context_around_date`: For each date belonging to a `context` (i.e. -2, -1, 0, 1, 2), 
    creates a dictionary holding all dates in a single structure.


- `get_context_dates_dict`: Wrapper function. Applies the past two to each date in the
    selected heatwave period. 


- `get_context_arrays_from_climatology`: Extracts from the big E-OBS volume the slices
    selected in the previous functions, filtered for a given extent.


- `calculate_statistics`: Applies the definition of a heatwave to these extracted data. In
    here we are using **the 90th percentile of the climatology** to decide whether a point
    on a given date is a heatwave or not. Values not being heatwaves will be masked. This
    operation is done for a whole `context`, so that in the next function there is some
    more information to decide if it is a heatwave.


- `is_heatwave`: Operates between dates in the `context` to identify actual heatwaves. 
    Wherever these operations add up `2` means that part is candidate to be labeled
    as heatwave.


- `heatwave_calc`: Wrapper function. Contains the pipeline and returns is heatwave? 
    calculation.

#### Event Synchronization

This step implements the schema shown in Fig. 3 of (Malik et al., 2011). This pipeline 
receives at the input the binarized E-OBS `tx` dataset and will apply some intermediate
operations to it to create a directed and an undirected adjacency matrices. All these
intermediate operations **have been vectorized** with respect to the original version. First,
it will compute `J_forward` and `J_backward` (i.e. J_ij and J_ji in the paper), which
check if `today` and `tomorrow` (for the whole raster) are both part of a heatwave. Then 
will calculate `big_q` and `small_q`. The former provides an indication of the 'strength of
the heatwave synchronization' between each pair of points, and the latter the delay
behaviour (i.e. whether an event at `i` always precedes an event at `j`). Finally, these
are turned into the `AQ` and `Aq` adjacency matrices that are understood by the 
complext network.

- `forward_calc`: Checks if `today` contains a heatwave event that
    precedes one happening `tomorrow`. 


- `backward_calc`: Checks if `tomorrow` contains a heatwave event
    that was preceded by `today`


- `big_q_calc`: Calculates the *strength of the event synchronization*.
    This means there is an **UNdirected** link between each pair of points. 


- `small_q_calc`: Calculates the *delay behaviour*. This means there is a **directed**
    link between each pair of points. The direction comes from the **sign** of this 
    operation. 


- `adjacency_big_q_calc`: Calculates the UNdirected adjacency matrix (symmetric). Here (I think)
    a `1` means there is a link between `i` and `j` but is not (necessarily) direct. 


- `adjacency_big_q_calc`: Calculates the directed adjacency matrix (non-symmetric). A `1`
    means there is a direct link from `i` to `j` and not viceversa.


- `eventsync_calc`: Wrapper function. Contains the pipeline and returns the adjacency
    matrices.

#### Modelling with complex networks

- `complexnetwork_modelling`: Wrapper function. For each date in a heatwave period, fetches
    its associated adjacency matrices and creates a directed graph and an undirected graph.
    Five network coefficients (i.e. betweenness centrality, degree of centrality, clustering
    coefficient, input degree, and output degree) are calculated. The function packs the
    results and returns them for later visualization.


- `unpack_and_stack_rasters`: This function reshapes the results of the modelling with complex
    networks and stacks them in numpy arrays to ease the subsequent operations.


- `complexnetwork_summary_coefs`: Calculates the mean of the stacks created in the previous
    function. Currently, these averages are not used anywhere else (too flat), but this function could
    be modified to obtain more metrics if needed.


#### Visualization and utilities

- `get_bounds_for_extent_name`: For a given extent name (i.e. Iberia, Mediterranean, Europe) 
    returns the coordinates of the associated bounding box.


- `generate_panel_titles_from_dates`: Formats panel titles


- `prepare_coefficients_for_plotting`: Groups the network coefficients in four parts and averages, so that
    them, so that the evolution of a heatwave can be explored in a visualization.


- `save_as_raster_files`: Takes the stacks from the function `unpack_and_stack_rasters` for each
    coefficient and saves them in NetCDF format in another folder to make visualizations later on.


- `plot_coefficients_per_period`: Plots the coefficients grouped in parts in the function
    `prepare_coefficients_for_plotting`. The figures are saved in a folder. 


- `maximize`: Just a function to save the visualizations maximized, thus avoiding cluttering.

#### References

- Malik, N., Bookhagen, B., Marwan, N. et al. Analysis of spatial and temporal extreme monsoonal rainfall over South Asia using complex networks. Clim Dyn 39, 971–987 (2012). https://doi.org/10.1007/s00382-011-1156-4