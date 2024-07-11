require(gWidgets2)
library(gWidgets2tcltk)

w <- gwindow("Hello...", visible=FALSE)       ## a parent container
g <- ggroup (cont = w)                        ## A box container
b <- gbutton("Click me for a message", cont=g, expand=TRUE)  ## some control

addHandlerClicked(b, function(...) {          ## adding a callback to an event
  gmessage("Hello world!", parent=w)          ## a dialog
})

visible(w) <- TRUE                            ## a method call

### GUI for CLD  ver 0.2

if (.Platform$OS.type == "unix")
{
  setwd("~/Dropbox/AQUABC/scripts")
  inputs_path = "~/AQUABC/INPUTS"
  outputs_path = "~/AQUABC/OUTPUTS"
  exe_path = "~/AQUABC"
} else
{
  setwd(file.path(Sys.getenv("USERPROFILE"), "Dropbox/AQUABC/scripts"))
  inputs_path = "C:/AQUABC/INPUTS"
  outputs_path = "C:/AQUABC/OUTPUTS"
  exe_path = "C:/AQUABC"
}


## cairoDevice allows plotting in GTK ggraphics widgets
library(cairoDevice)

###  set working directory to the scripts to load include files
setwd(file.path(Sys.getenv("USERPROFILE"), "Dropbox/AQUABC/scripts"))
