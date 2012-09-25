## setMethod("widget", "Protocol",
##           function(object, input, output = NULL, gg = ggobi(), apply_cb = NULL) 
## {
##   # attempt to translate slots to a GUI
  
##   # get types, defaults, and current values of parameters
##   types <- getClass(class(object))@slots
##   types <- types[!(names(types) %in% slotNames("Protocol"))]
##   proto <- attributes(getClass(class(object))@prototype)
##   params <- parameters(object)
##   slots <- names(types)
  
##   # build radio groups
##   choices <- sapply(slots, function(slot)
##     extends("character", types[[slot]]) && length(eval(types[[slot]])) > 1)
##   radios <- sapply(slots[choices], function(slot) {
##     group <- gtkVBox(TRUE, 5)
##     radios <- list()
##     options <- eval(proto[[slot]])
##     for (option in options) {
##       button <- gtkRadioButton(radios, option)
##       radios <- c(radios, button)
##       group$packStart(group, button, FALSE, FALSE)
##     }
##     radios[[match(eval(params[[slot]]), options)[1]]]$setActive(TRUE)
##     frm <- gtkFrame(slot)
##     frm$add(group)
##     frm
##   })
  
##   # get checkboxes
  
##   toggles <- sapply(slots, function(slot) extends("logical", types[[slot]]))
##   checks <- sapply(slots[toggles], function(slot) {
##     input <- gtkCheckButton(slot)
##     input$setActive(eval(params[[slot]]))
##     input
##   })
  
##   # get label/entry pairs
  
##   sizegroup <- gtkSizeGroup()
##   sizegroup$setMode("horizontal")
##   entries <- !choices & !toggles
##   boxes <- sapply(slots[entries], function(slot) {
##     lab <- gtkLabel(slot)
##     lab["xalign"] <- 0
##     input <- gtkEntry()
##     input$setText(as.character(params[[slot]]))
##     box <- gtkHBox(FALSE, 5)
##     box$packStart(lab, FALSE, FALSE)
##     box$packStart(input, TRUE, TRUE)
##     sizegroup$addWidget(input)
##     box
##   })
  
##   applybutton <- gtkButton(stock = "gtk-apply")
##   gSignalConnect(applybutton, "clicked", function(wid) {
##     # put stuff back into object
##     # then perform
##     output <- perform(object, input)
##     explorer <- explore(object, input, output)
##     if (!is.null(apply_cb))
##       apply_cb(object, output)
##   })
##   helpbutton <- gtkButton(stock = "gtk-help")
##   gSignalConnect(helpbutton, "clicked", function(wid) {
##     help(paste(class(object), "class", sep="-"))
##   })
  
##   box <- gtkVBox(FALSE, 5)
##   widgets <- c(radios, checks, boxes)
##   names(widgets) <- c(slots[choices], slots[toggles], slots[entries]) 
##   sapply(widgets[slots], box$packStart, FALSE, FALSE)
##   buttonbox <- gtkButtonBox()
##   buttonbox$setStyle("edge")
##   buttonbox$add(helpbutton)
##   buttonbox$add(applybutton)
##   box$packStart(buttonbox, FALSE, FALSE)
  
##   explorer <- explore(object, input, output, gg)
  
##   main <- gtkHBox(FALSE, 5)
##   main$packStart(box, FALSE, FALSE)
##   main$packStart(explorer, TRUE, TRUE)
  
##   main
## })

## progress_bar <- function(bar) {
##   cursor <- gdkCursorNew("watch")
##   window <- bar$getToplevel()$window
##   prev_frac <- 0
##   function(fraction, task) {
##     if (!missing(task)) {
##       bar$setText(task)
##       update_gui()
##     }
##     if (!missing(count)) {
##       if (fraction == 0) # begin
##         window$setCursor(cursor)
##       else if (fraction == 1) { # finished
##         window$setCursor(NULL)
##         fraction <- 0
##       } else if (fraction - prev_frac < 0.01) # no need to update
##         return()
##       bar$setFraction(fraction)
##       prev_frac <<- fraction
##       update_gui()
##     }
##   }
## }

## update_gui <- function() while(gtkEventsPending()) gtkMainIteration()
