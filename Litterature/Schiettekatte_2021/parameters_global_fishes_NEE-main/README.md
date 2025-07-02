
# Code and data of supplementary methods of “Biological trade-offs underpin coral reef ecosystem functioning”

The project uses a drake workflow, which binds all elements of the
analyses together with a plan. To reproduce output, make sure all
packages are installed (listed in R/packages.R) and run
`drake::r_make()`. CAUTION: the code uses 50 cores at times and takes a
couple of days to run. If you try to run this project on a computer with
less than 50 cores, it will crash. You can run the code on a computer
with less cores, but you will need to decrease the number of cores used
in the code. Be aware that decreasing the number of cores will
substantially increase computation time.
