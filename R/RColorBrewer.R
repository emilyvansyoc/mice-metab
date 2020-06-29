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
treatmentIntcols <- c("EX_AL" = hex[7], "EX_ER" = hex[8], "SED_AL" = hex[9], "SED_ER" = hex[10])


# accent colors
accentcols <- c(hex[10:12])
