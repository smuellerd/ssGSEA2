## Internal assignment to avoid data.table note about `.`
`.` = list

`%dopar%` = foreach::`%dopar%`

## Prevent note "no visible binding for global variable"
utils::globalVariables(c("gs.i"))
