# River-line fixtures

`os_open_rivers_extract/` is a 24-link extract of OS Open Rivers
(British National Grid E 368500-374500, N 350500-356000), laid out as the inside
of the national download zip (`Data/oprvrs_gb.gpkg`, layer `watercourse_link`,
and `Doc/licence.txt`). It is stored unzipped because `*.zip` is gitignored;
the tests zip it into a cache directory at run time, so
`lsdtt-fetch-river-lines` reads it through `river_cache_dir` with no network.

Contains OS data (c) Crown copyright and database rights 2026. Licensed under
the Open Government Licence v3.0 (http://os.uk/opendata/licence).
