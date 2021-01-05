map_mtb <- function(Y, idx_mtb, MTB, pal, pal_rev, scale_range) {

  oo_fia = merge(MTB, Y, by = "NAME")

  theMap <- leaflet(options = leafletOptions(zoomControl = FALSE))  %>%
    #addTiles() %>% 
    #addTiles()  %>% 
    #addProviderTiles("Esri.WorldTerrain") %>%
    #addProviderTiles("Esri.WorldPhysical") %>%
    addProviderTiles("Esri.WorldShadedRelief") %>%
    addPolygons(data = oo_fia, 
                fillColor = ~pal(oo_fia@data$score), 
                stroke = F, 
                fillOpacity = .7  ) %>% #, 
                #color = "black", 
                #weight = 0.1,
                #dashArray = "" #, 
                #highlight = highlightOptions(weight = 5, 
                #                            color = ~colorNumeric(palette="Greys", domain=oo_fia@data$score)(oo_fia@data$score), 
                #                            dashArray = "", 
                #                            fillOpacity = 0.8, 
                #                            bringToFront = TRUE),
                #labelOptions = labelOptions( style = list("font-weight" = "normal", padding = "3px 8px"), 
                 #                            textsize = "13px", 
                #                             direction = "auto")
    addPolygons(data = germany_0,
                fillOpacity = 0, 
                color = "black",
                weight = 1)  %>%
    addLegend(pal = pal_rev,
              values = scale_range,
              position = "bottomright",
              title = FALSE, 
              labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
  
  theMap
}
