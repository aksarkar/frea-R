requireNamespace('Cairo')
requireNamespace('ggplot2')
requireNamespace('dplyr')

#' @include heatmap.R scales.R

verify_impg_preprocess <- function(orig_file, impg_file) {
    pgc_scz <- (read.table(gzfile(orig_file), header=T, sep='\t',
                           comment.char='', nrow=10e6) %>%
                mutate(snp=snpid, pos=bp, ref=a1, alt=a2,
                       z=sqrt(qchisq(p, 1, lower.tail=F)) * sign(or - 1)) %>%
                select(snp, pos, ref, alt, z))
    my_pre <- read.table(gzfile(impg_file), header=F, sep=' ',
                         comment.char='', nrow=9e6,
                         col.names=c('snp', 'pos', 'ref', 'alt', 'z'))
    return(list(pgc_scz, my_pre))
}

parse_impg_corr <- function(ind_file, impg_file, nrow=10e6) {
    ind_imputed <- read.table(ind_file, header=F, sep=' ', comment.char='',
                                  colClasses=c(rep('character', 4), 'numeric'),
                                  col.names=c('snp', 'pos', 'ref', 'alt', 'z'),
                                  nrow=nrow)
    impg_imputed <- read.table(impg_file, header=F,
                                   colClasses=c(rep('character', 4), rep('numeric', 2)),
                                   col.names=c('snp', 'pos', 'ref', 'alt', 'z', 'r2'),
                                   nrow=nrow, comment.char='')
    impg_corr <- (impg_imputed %>% dplyr::filter(r2 > 0.6) %>%
                  dplyr::inner_join(ind_imputed, by=c('snp', 'pos', 'ref', 'alt')))
    return(impg_corr)
}

plot_impg_corr <- function(impg_corr) {
    my_plot <- (ggplot(data=impg_corr, aes(y=z.x, x=z.y)) +
                labs(y='Imputed z-score', x='True z-score') +
                geom_hex(aes(fill=log10(..value..)), binwidth=c(.25,.25)) +
                coord_equal() +
                scale_heatmap(name='Log count') +
                theme_nature)
    Cairo(file='impg-holdout-corr.pdf', type='pdf', height=89, width=89, units='mm')
    print(my_plot)
    dev.off()
}
