# Stack segment summaries from a list of true (ground-truth) paths

Stack segment summaries from a list of true (ground-truth) paths

## Usage

``` r
summarize_segments(path_list, group_label)
```

## Arguments

- path_list:

  a list of path objects, each with a \`\$segments\` tibble.

- group_label:

  label to attach to every row (e.g. experimental group).

## Value

a tibble with columns \`durations\`, \`states\`, \`speeds\`, \`label\`,
\`path_id\`.
