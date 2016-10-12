requireNamespace('Cairo')
requireNamespace('dplyr')
requireNamespace('ggplot2')
requireNamespace('grid')
requireNamespace('gtable')

#' @include heatmap.R scales.R theme_nature.R
#' @importFrom grid unit
NULL

enrichment_by_cluster <- function(enrichments, cluster_density, scale_name) {
    (heatmap(ggplot(enrichments, aes(x=cluster, y=pheno, fill=score))) +
     scale_heatmap(name=scale_name) +
     labs(x='Enhancer module', y='Phenotype') +
     theme_nature +
     theme(axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           legend.position='right',
           plot.margin=unit(c(2, rep(0, 3)), 'mm')))
}

parse_enhancer_enrichments <- function(filename) {
    enrichments <- read.delim(filename, header=FALSE, sep=' ')
    enrichments$pheno <- factor(toupper(sub('^[^-]*-', '', enrichments$V1)))
    enrichments$pheno <- factor(enrichments$pheno, levels=rev(phenotype_ordering))
    subset(enrichments, V6 > 0 & V7 > 0)
}

plot_enhancer_enrichments_by_eid <- function(filename) {
    data(roadmap_sample_info)
    enrichments <- parse_enhancer_enrichments(filename)
    enrichments$eid <- factor(enrichments$V3, levels=eid_ordering)
    enrichments <- transform(enrichments, val=(V5 - V6)/ sqrt(V7))

    my_ggplot <- (heatmap(ggplot(enrichments, aes(x=eid, y=pheno, fill=val))) +
                  scale_heatmap(name='z-score') +
                  scale_y_discrete(limits=rev(levels(enrichments$pheno))) +
                  labs(x='Reference epigenome', y='Phenotype') +
                  theme(axis.text.x=element_blank(),
                        legend.position='bottom'))
    my_gtable <- ggplotGrob(my_ggplot)

    ## Add tissue colors
    tissue_grob <- ggplotGrob(epigenome_by_tissue())$grobs[[4]]
    my_gtable <- gtable::gtable_add_rows(my_gtable, unit(2, 'mm'), 0)
    my_gtable <- gtable::gtable_add_grob(my_gtable, tissue_grob, t=1, l=4, r=4)

    ## Add tissue legend
    my_legend <- roadmap_tissue_legend(list(legend.position='bottom'), list(direction='horizontal', nrow=2))
    my_gtable <- gtable::gtable_add_rows(my_gtable, unit(10, 'mm'))
    my_gtable <- gtable::gtable_add_grob(my_gtable, my_legend, t=dim(my_gtable)[1], l=1, r=dim(my_gtable)[2])

    Cairo(type='pdf', file=sub('.in$', '.pdf', filename), width=210, height=50, units='mm')
    grid::grid.draw(my_gtable)
    dev.off()
}

plot_enhancer_enrichments <- function(filename, cluster_density, plot_log_fold=FALSE, plot_fold=FALSE, flip=FALSE) {
    enrichments <- parse_enhancer_enrichments(filename)
    enrichments$cluster <- factor(enrichments$V3, levels=row.names(cluster_density))
    if (plot_fold) {
        enrichments <- transform(enrichments, score=V5 / V6)
        label <- 'Fold enrichment'
    }
    else if (plot_log_fold) {
        enrichments <- transform(enrichments, score=log10(V5) - log10(V6))
        label <- 'Log fold enrichment'
    }
    else {
        enrichments <- transform(enrichments, score=(V5 - V6) / sqrt(V7))
        label <- 'z-score'
    }

    my_density <- (density_by_cluster(cluster_density, keep=unique(enrichments$cluster)) +
                   theme(legend.position='right',
                         axis.text.y=element_text(margin=margin(2))))
    my_gtable <- gtable:::rbind.gtable(ggplotGrob(enrichment_by_cluster(enrichments, cluster_density, label)),
                                       ggplotGrob(my_density),
                                       size='last')

    ## Add tissue colors
    cluster_by_tissue_grob <- ggplotGrob(cluster_by_tissue(cluster_density, keep=unique(enrichments$cluster)))
    my_gtable <- gtable::gtable_add_rows(my_gtable, unit(4, 'mm'), 7)
    my_gtable <- gtable::gtable_add_grob(my_gtable, cluster_by_tissue_grob, t=8, l=4)

    ## Add tissue legend
    my_gtable <- gtable::gtable_add_rows(my_gtable, unit(20, 'mm'))
    my_legend <- roadmap_tissue_legend(list(legend.position='bottom'), list(direction='horizontal', nrow=2))
    my_gtable <- gtable::gtable_add_grob(my_gtable, my_legend, t=-1, l=1, r=-1)

    Cairo(type='pdf', file=sub('.in$', '.pdf', filename), width=190, height=80, units='mm')
    grid::grid.draw(my_gtable)
    dev.off()
}

plot_diagnostic <- function(filename, cluster_density) {
    enrichments <- parse_enhancer_enrichments(filename)
    enrichments$sig <- enrichments$V8 < 0.005
    enrichments$cluster <- factor(enrichments$V3, levels=row.names(cluster_density))
    p <- (ggplot(enrichments, aes(x=cluster, color=cluster, alpha=sig)) +
          geom_point(aes(y=V5), size=0.1) +
          geom_linerange(aes(ymin=V6-sqrt(V7), ymax=V6+sqrt(V7)), size=0.25) +
          scale_x_discrete(drop=FALSE) +
          scale_alpha_manual(values=c('TRUE'=1, 'FALSE'=0.25)) +
          scale_color_manual(values=color_by_cluster_top(cluster_density)) +
          coord_flip() +
          labs(x='Overlaps', y='Cluster') +
          facet_wrap(~ pheno, nrow=1) +
          theme_nature +
          theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
                panel.grid.major=element_line(size=0.1, color='gray80')))
    Cairo(type='pdf', file=sub('.in$', '-diagnostic.pdf', filename), width=190, height=270, units='mm')
    print(p)
}

plot_enrichment_and_pve <- function(enrichment_file, pve_file, cluster_density) {
    pve <- (read.table(pve_file, sep=' ') %>%
            dplyr::mutate(cluster=factor(V1, levels=row.names(cluster_density))) %>%
            dplyr::select(cluster, pve=V3, se=V4) %>%
            dplyr::filter(cluster %in% enrichments$cluster))
    pve_plot <- (ggplot(pve, aes(x=cluster, y=pve, ymin = pve - se, ymax = pve + se, color=cluster)) +
                 geom_point(size=.25) +
                 geom_linerange(size=.25) +
                 scale_color_manual(values=color_by_cluster_top(cluster_density)) +
                 labs(y='T1D PVE') +
                 theme_nature +
                 theme(axis.text.x=element_blank(),
                       axis.title.x=element_blank(),
                       plot.margin=unit(c(2, 0, 0, 0), 'mm')))

    enrichments <- (parse_enhancer_enrichments(enrichment_file) %>%
                    dplyr::mutate(cluster=factor(V3, levels=row.names(cluster_density)),
                                  sig=V8 < 0.002) %>%
                    dplyr::group_by(cluster) %>%
                    dplyr::filter(cluster %in% pve$cluster & TRUE %in% sig & pheno == "T1D") %>%
                    dplyr::select(pheno, cluster, obs=V5, null_mean=V6, null_var=V7, sig))
    diagnostic_plot <- (ggplot(enrichments, aes(x=factor(cluster), y=obs,
                                                ymin=null_mean - sqrt(null_var),
                                                ymax=null_mean + sqrt(null_var), color=cluster,
                                                alpha=sig, shape=sig)) +
                        geom_point(size=0.25) +
                        geom_linerange(size=0.25) +
                        scale_x_discrete(drop=FALSE) +
                        scale_alpha_manual(values=c('TRUE'=1, 'FALSE'=0.25)) +
                        scale_shape_manual(values=c('TRUE'=19, 'FALSE'=1)) +
                        scale_color_manual(values=color_by_cluster_top(cluster_density)) +
                        labs(y='T1D overlaps') +
                        theme_nature +
                        theme(axis.text.x=element_blank(),
                              axis.title.x=element_blank()))


    my_density <- (density_by_cluster(cluster_density, keep=unique(enrichments$cluster)) +
                   theme(legend.position='right',
                         axis.text.y=element_text(margin=margin(2))))
    my_gtable <- gtable:::rbind.gtable(gtable::gtable_add_cols(ggplotGrob(pve_plot), unit(1, "lines")),
                                       gtable::gtable_add_cols(ggplotGrob(diagnostic_plot), unit(1, "lines")),
                                       ggplotGrob(my_density),
                                       size='last')

    # Fix up the widths/heights since these get clobbered by bind
    my_gtable$widths <- grid:::unit.list(my_gtable$widths)
    my_gtable$heights <- grid:::unit.list(my_gtable$heights)
    widths <- list(unit(length(unique(pve$cluster)), 'null'))
    heights <- lapply(list(6, 6, 30), function(x) unit(x, 'null'))
    my_gtable$widths[my_gtable$layout$l[grepl("panel", my_gtable$layout$name)]] <- widths
    my_gtable$heights[my_gtable$layout$t[grepl("panel", my_gtable$layout$name)]] <- heights

    ## Add tissue colors
    cluster_by_tissue_grob <- ggplotGrob(cluster_by_tissue(cluster_density, keep=unique(enrichments$cluster)))
    my_gtable <- gtable::gtable_add_rows(my_gtable, unit(4, 'mm'), 7)
    my_gtable <- gtable::gtable_add_rows(my_gtable, unit(7, 'mm'), 13)
    my_gtable <- gtable::gtable_add_grob(my_gtable, cluster_by_tissue_grob, t=14, l=4)

    ## Add tissue legend
    my_gtable <- gtable::gtable_add_rows(my_gtable, unit(20, 'mm'))
    my_legend <- roadmap_tissue_legend(list(legend.position='bottom'), list(direction='horizontal', nrow=4))
    my_gtable <- gtable::gtable_add_grob(my_gtable, my_legend, t=-1, l=1, r=-1)

    Cairo(type='pdf', file=sub('.in$', '-pve.pdf', sig_enrichment_file), width=100, height=120, units='mm')
    grid::grid.draw(my_gtable)
    dev.off()
}
