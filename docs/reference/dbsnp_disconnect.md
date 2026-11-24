# Close the connection behind a [`dbsnp_tbl()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_tbl.md) result

Close the connection behind a
[`dbsnp_tbl()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_tbl.md)
result

## Usage

``` r
dbsnp_disconnect(tbl)
```

## Arguments

- tbl:

  A lazy table returned by
  [`dbsnp_tbl()`](https://github.com/vjcitn/YESCDS/reference/dbsnp_tbl.md)
  (or any `tbl_sql`).

## Value

`NULL`, invisibly.
