## file to establish plotting conventions

## plotting utility functions


## remove x and y axis information
remove_x_info <- function(plot, text = FALSE) {

    plot <- plot + theme(axis.title.x = element_blank())
    if (text) plot <- plot + theme(axis.text.x = element_blank())
    return(plot)

}

remove_y_info <- function(plot, text = FALSE) {

    plot <- plot + theme(axis.title.y = element_blank())
    if (text) plot <- plot + theme(axis.text.y = element_blank())
    return(plot)

}

remove_legend <- function(plot) {

    plot <- plot + theme(legend.position = "none")
    return(plot)

}

second_axis <- function(plot) {

    plot <- plot + guides(x.sec = "axis", y.sec = "axis")
    return(plot)

}


## create custom theme

theme_jabo <- function(){

    list(

        ## adjust theme elements
        theme(panel.grid.major = element_blank(),

              panel.grid.minor = element_blank(),

              panel.background = element_blank(),

              axis.title = element_text(size = 12),

              text = element_text(family = "sans"),

              plot.title = element_text(size = 15, face = "bold"),

              plot.subtitle = element_text(size = 14),

              legend.title = element_text(face = "bold"),

              axis.line = element_line(size = 0.3),

              axis.text.x.top = element_blank(),

              axis.text.y.right = element_blank(),

              axis.ticks.x.top = element_blank(),

              axis.ticks.y.right = element_blank()

              ##legend.key = element_rect(fill = "transparent")

          )

    )
    
}
