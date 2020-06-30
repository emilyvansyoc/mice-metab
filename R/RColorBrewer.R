### QUICK CODE FOR RCOLORBREWER

require(RColorBrewer)

# show all palettes
display.brewer.all(colorblindFriendly = TRUE)

# pick one palette
display.brewer.pal(n = 12, name = "Paired")

# return the hex code
hex <- brewer.pal(n = 12, name = "Paired")

# save the desired colors to your variables
tissuecols <- c("Plasma" = hex[1], "Tumor" = hex[2])
treatmentcols <- c("EX" = hex[3], "SED" = hex[4], "AL" = hex[5], "ER" = hex[6])
treatmentIntcols <- c("EX-AL" = hex[7], "EX-ER" = hex[8], "SED-AL" = hex[9], "SED-ER" = hex[10])


# accent colors
accentcols <- c(hex[8:12])

## let's get more colors than RColorBrewer allows
getPalette <- colorRampPalette(brewer.pal(n = 9, name = "Paired"))


