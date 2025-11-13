# Phase_Optimisation
Optimises the geometry of an acoustic metadiffuser using the minimisation algorithm `fmincon`.

The objective function seeks to minimise the function $e = \arg(R_{meta})-\arg(R_{QWR})$ where $e$ is the error $R_{meta}$ is the reflection coefficient of the metadiffuser and $R_{QWR}$ is the reflection coefficient of the chosen Schroeder diffuser.