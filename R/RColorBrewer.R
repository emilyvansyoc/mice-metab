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
#treatmentcols <- c("EX" = hex[3], "SED" = hex[4], "AL" = hex[5], "ER" = hex[6])
treatmentIntcols <- c("PA-AL" = hex[7], "PA-ER" = hex[8], "SED-AL" = hex[9], "SED-ER" = hex[10])
timecols <- c("Day 7" = hex[3], "Day 21" = hex[4], "Day 35" = hex[5])


# accent colors
accentcols <- c(hex[8:12])

## let's get more colors than RColorBrewer allows
getPalette <- colorRampPalette(brewer.pal(n = 9, name = "Paired"))

#### ---- GREYSCALE ----

grey <- brewer.pal(n = 4, name = "Greys")

# save colors to variables
treatmentGreys <- c("SED-AL" = hex[1], "PA-AL" = grey[2], "SED-ER" = hex[3], "PA-ER" = hex[4])

