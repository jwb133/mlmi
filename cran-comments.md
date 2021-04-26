## Test environments
* local Windows install, R 4.0.5
* Travis-CI on R version 4.0.2 (2020-06-22)
* win-builder on devel 4.1.0 alpha (2021-04-25 r80223)
* win-builder on release version 4.0.5

## R CMD check results
On CRAN the previous version has a note: 'LazyData' is specified without a 'data' directory but I'm not sure how to fix this. CRAN check on r-devel-windows-x86_64-gcc10-UCRT showing only showing some packages not available, but I think they are still available.
