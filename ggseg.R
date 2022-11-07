
tiff('surface_area.png',units="cm",width=16,height=16,res=600)

sa <- read_excel("ggseg_surface_area.xlsx")
sa$Group <- as.character(sa$groups) 
sa$Group <- factor(sa$Group, levels=unique(sa$Group))

p1 <- sa %>%
  group_by(Group) %>%
 ggseg(atlas=dk, 
             colour="white", 
             hemi="left",
             mapping=aes(fill=Beta,)) +
  facet_wrap(~Group, ncol=1) +
  theme(legend.position = "right") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x=element_blank(), 
       axis.ticks.y=element_blank()) +
  scale_fill_gradientn(colours = c("royalblue","firebrick","goldenrod"),na.value="grey")


p1
dev.off()

thickness <- read_excel("ggseg_thickness.xlsx")
thickness$Group <- as.character(thickness$groups) 
thickness$Group <- factor(thickness$Group, levels=unique(thickness$Group))


tiff('thickness.png',units="cm",width=16,height=16,res=600)
p2 <- thickness %>%
  group_by(Group) %>%
  ggseg(atlas=dk, 
        colour="white", 
        hemi="left",
        mapping=aes(fill=Beta,)) +
  facet_wrap(~Group, ncol=1) +
  theme(legend.position = "right") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x=element_blank(), 
        axis.ticks.y=element_blank()) +
  scale_fill_gradientn(colours = c("royalblue","firebrick","goldenrod"),na.value="grey")
p2
dev.off()
