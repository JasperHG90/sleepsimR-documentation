### Nested loop plot
#' @keywords internal
nlp <- function(data, methodvar, by, stats, target, top, col = NULL, linesize=1) {
  ### Compute internal data
  opts <- lapply(X = by, FUN = function(x) levels(data[[x]]))
  names(opts) <- by
  dgms <- do.call(expand.grid, opts)
  dgms[[".scenario"]] <- seq(nrow(dgms))
  data <- merge(x = data, y = dgms)
  data <- data[order(data[[".scenario"]]), ]
  
  ### Compute limits and placement of nested loop plot labels
  limits <- range(data[["est"]], na.rm = TRUE)
  delta <- diff(range(data[["est"]])) / 10
  placement <- vector(mode = "list", length = length(by))
  for (i in seq_along(placement)) {
    if (i == 1) {
      if (top) {
        placement[[i]] <- c(round(limits[2], digits = 2) + delta, round(limits[2], digits = 2) + 2 * delta)
      } else {
        placement[[i]] <- c(round(limits[1], digits = 2) - 2 * delta, round(limits[1], digits = 2) - delta)
      }
    } else {
      if (top) {
        placement[[i]] <- c(placement[[i - 1]][2] + delta, placement[[i - 1]][2] + 2 * delta)
      } else {
        placement[[i]] <- c(placement[[i - 1]][1] - 2 * delta, placement[[i - 1]][1] - delta)
      }
    }
  }
  
  ### Rescale variables included in the nested loop plot
  for (i in seq_along(by)) {
    data[[paste0(".", by[i])]] <- scales::rescale(x = as.numeric(data[[by[i]]]), to = placement[[i]])
  }
  
  ### Build basic plot
  if (!is.null(methodvar)) {
    methodvar <- rlang::sym(methodvar)
    gg <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = .scenario, y = est, group = !!methodvar)) +
      ggplot2::geom_hline(yintercept = target, linetype = "dotted") +
      ggplot2::geom_step(mapping = ggplot2::aes(color = !!methodvar), size=linesize)
  } else {
    gg <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = .scenario, y = est)) +
      ggplot2::geom_hline(yintercept = target, linetype = "dotted") +
      ggplot2::geom_step(color = col, size=linesize)
  }
  gg <- gg +
    ggplot2::labs(x = paste0(paste(vapply(X = by, FUN = function(x) length(levels(data[[x]])), FUN.VALUE = numeric(1)), collapse = " x "), " = ", 
                             max(data[[".scenario"]]), " ordered scenarios"), y = stats)
  
  ### Build and add legends of nested loop plot
  for (i in seq_along(by)) {
    .tmp <- rlang::sym(paste0(".", by[i]))
    if(by[i] == "zeta") {
      lab <- latex2exp::TeX(paste0("$\\", by[i], "$: ", paste(levels(data[[by[i]]]), collapse = ", ")))
    } else {
      lab <- latex2exp::TeX(paste0("$", toupper(by[i]), "$: ", paste(levels(data[[by[i]]]), collapse = ", ")))
    }
    
    print(lab)
    gg <- gg +
      ggplot2::geom_step(mapping = ggplot2::aes(y = !!.tmp), size=linesize * (2/3)) +
      ggplot2::annotate(geom = "text", x = 1, y = placement[[i]][2] + delta / 2, 
                        label = lab, 
                        hjust = 0, vjust = 0.5)
  }
  
  ### Return plot
  return(gg)
}
