
theme_awesome <- function (base_family = "Avenir") {
  ggthemes::theme_fivethirtyeight(base_family = base_family) %+replace% 
    theme(
      panel.grid.major.x = element_blank() ,
      panel.grid.major.y = element_blank() ,
      plot.title = element_text(hjust=0.5, size=28,face="bold"),
      plot.subtitle = element_text(hjust=0.5, size=20, face="italic"),
      legend.position = "bottom",
      legend.text=element_text(size=16, face="bold"),
      axis.text.x = element_text(size=12,face="bold"),
      axis.text.y = element_text(size=12,face="bold"),
      axis.title = element_text(size=16,face="bold"), 
      axis.title.x = element_blank(),
      plot.background=element_rect(fill="white"),
      panel.background=element_rect(fill="white"),
      legend.background=element_rect(fill="white")
    )   
}

theme_fancy <- function() {
  theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(hjust=0.5, size=28,face="bold"),
          plot.subtitle = element_text(hjust=0.5, size=20, face="italic"),
          legend.position = "bottom",
          plot.background=element_rect(fill="white"),
          panel.background=element_rect(fill="white"),
          legend.background=element_rect(fill="white"),
          panel.border = element_blank()
          )
}