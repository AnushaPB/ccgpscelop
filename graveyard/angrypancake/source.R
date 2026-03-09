
  # Get sources
  fut_targets = future_map_bin
  cur_targets = cur_ref_targets
  find_source <- function(fut_targets, cur_targets) {
    # Identify cells of type
    type <- unique(values(cur_targets, na.rm = TRUE))[1]
    cur_targets[cur_targets == fut_targets | cur_targets == 0.5 | fut_targets == 0.5] <- NA
    fut_targets[fut_targets == cur_targets | fut_targets == 0.5 | cur_targets == 0.5] <- NA
    cells_current <- which(values(cur_targets) == 1)
    cells_future <- which(values(fut_targets) == type)

    # Get coordinates of those cells
    coords_current <- xyFromCell(cur_targets, cells_current)
    coords_future <- xyFromCell(fut_targets, cells_future)

    # Convert to SpatVector
    v_current <- vect(coords_current, type = "points", crs = crs(cur_targets))
    v_future  <- vect(coords_future,  type = "points", crs = crs(future_map_bin))

    # Find nearest type cell in future for each type cell in current
    # (using fast nearest-neighbor search)
    nn <- terra::nearest(v_future, v_current)

    # Tabulate how often each CURRENT cell was chosen as the nearest "source"
    nn_df  <- as.data.frame(nn) %>% filter(distance < 100000)
    to_idx <- nn_df$to_id                       # indices into v_current / cells_current
    counts <- tabulate(to_idx, nbins = length(cells_current))
    
    source_count <- rast(cur_targets)
    values(source_count) <- NA                 # default 0 everywhere
    source_count[cells_current] <- counts     # counts only at current type-1 cells
    source_count[source_count == 0] <- NA              # set 0 counts to NA
    tp(source_count)
  }
