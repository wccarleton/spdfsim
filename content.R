pagetitle <- h1(strong("SPDF Simulator"),style="color: white")

shorthead <- h4(em("A simple app for simulating a Summed Probability Density Function.",style="color: white"))

plotspacetitle <- h2(strong("Plot"))

paramstitle <- h2(strong("Parameters"))

explainertext1 <- p("This R Shiny app is intended for exploring summed probability density functions (SPDFs). An SPDF is simply a sum of radiocarbon-date densities. Archaeologists and palaeoclimatologists have been using these functions as proxies for past processes like human population level change and climate dynamics. They can be very easy to misinterpret, however. This app is intended to allow interested users to see what sorts of SPDFs can be generated given a set of simple parameters. For more information about SPDFs see",a(href="https://www.sciencedirect.com/science/article/pii/S0305440311002482?via%3Dihub",target="_blank",strong("Williams (2012)"),style="color: white"),style="color: white")

explainertext2 <- p("To explore SPDFs, simply select the start date, span of time, number of radiocarbon-dated events, and functional form of a theoretical process (along with any relevant parameters for that process) by entering the required information in the relevant fields. The app should responsively produce a plot displaying the true process (bottom plot panel), a random sample of that process with no chronological uncertainty (middle panel), and an SPDF based on simulated calibrated radiocarbon dates corresponding to the sample of event times (top panel). Note that the procedure can take a little while depending on the size of the sample (N. Dates parameter).",style="color: white")

contact <- p("Concerns or questions? Email wcarleton[at]ice.mpg.de",style="color: white")
